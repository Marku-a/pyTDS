"""
D1 Fantasia: cardio-respiratory coupling (RESP -> HR), 20 young vs 20 elderly.

Implements section "D1 Fantasia" of SPEC.md. System = [resp, hr] at 4 Hz, MAX_WINDOW_S = 300.
Tests per variant: within-subject pair_test(resp_i, hr_i) and cross-subject null
pair_test(resp_i, hr_j), j = next eligible record in sorted order (wrap around).

Deviations / choices (see report):
- With --n-subjects < 40 (and --smoke) the subset is balanced young/elderly (alternating
  sorted young and sorted elderly records) instead of the first n sorted records, which
  would all be elderly. With all subjects the order is plain sorted order.
- RR filter (deviation from the spec's "previous accepted RR"): each RR is compared with the
  median of the last 5 accepted RRs (20 % rule and 0.3-2.0 s range kept); after 5 consecutive
  rejections the reference resets to the median of those 5. Records with a longest interpolated
  HR gap > 60 s or a longest RESP NaN gap > 60 s are skipped (reported in meta). RESP NaNs are
  linearly interpolated before filtering. Sampling rate is read from the header (f2y02 is 333 Hz);
  resampling to 4 Hz uses the exact rational ratio.
- Beat symbols = wfdb annotation symbols whose label code is flagged by wfdb's is_qrs table.
- Records are downloaded whole (.hea/.dat/.ecg, ~7 MB each).
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from common import (VARIANTS, auc, group_effect, pair_test, params_row, stable_seed,  # noqa: E402
                    variant_params)

DATASET = "fantasia"
FS = 4.0
MAX_WINDOW_S = 300
MAX_SAMPLES = 110 * 60 * 4  # 26,400
MIN_SAMPLES = 60 * 60 * 4
BASE = "https://physionet-open.s3.amazonaws.com"
VERSION = "fantasia/1.0.0"
RECORDS = [f"f{g}{a}{i:02d}" for g in (1, 2) for a in ("y", "o") for i in range(1, 11)]


# ---------------------------------------------------------------- data
def list_keys(prefix: str) -> dict[str, int]:
    keys, token = {}, None
    while True:
        url = f"{BASE}/?list-type=2&prefix={prefix}" + (f"&continuation-token={token}" if token else "")
        xml = subprocess.run(["curl", "-sSf", "--retry", "3", url], capture_output=True, text=True, check=True).stdout
        for k, s in re.findall(r"<Key>([^<]*)</Key>.*?<Size>(\d+)</Size>", xml):
            keys[k] = int(s)
        if "<IsTruncated>true</IsTruncated>" in xml:
            token = re.search(r"<NextContinuationToken>([^<]*)<", xml).group(1)
        else:
            return keys


def download(data_dir: str, records: list[str]) -> dict[str, int]:
    os.makedirs(data_dir, exist_ok=True)
    keys = list_keys(VERSION + "/")
    prov = {}
    for rec in records:
        for ext in ("hea", "dat", "ecg"):
            fn = f"{rec}.{ext}"
            dst = os.path.join(data_dir, fn)
            if not os.path.exists(dst):
                subprocess.run(["curl", "-sSf", "--retry", "3", "-o", dst, f"{BASE}/{VERSION}/{fn}"], check=True)
            prov[fn] = os.path.getsize(dst)
    return prov


def beat_symbols() -> set[str]:
    from wfdb.io import annotation as a
    t = a.ann_label_table
    return {t.loc[c, "symbol"] for c in t.index if a.is_qrs[c]}


def _longest_run(mask: np.ndarray) -> int:
    best = cur = 0
    for m in mask:
        cur = cur + 1 if m else 0
        best = max(best, cur)
    return best


def preprocess(data_dir: str, rec: str, beats: set[str]):
    """Return ((resp, hr), reason, info) at 4 Hz; data is None if the record is skipped."""
    import wfdb
    from fractions import Fraction
    from scipy.signal import butter, filtfilt, resample_poly

    path = os.path.join(data_dir, rec)
    r = wfdb.rdrecord(path)
    ann = wfdb.rdann(path, "ecg")
    fs = float(r.fs)
    afs = float(getattr(ann, "fs", None) or fs)
    info = {"fs": fs}
    resp = r.p_signal[:, r.sig_name.index("RESP")].astype(float)
    nan = ~np.isfinite(resp)
    info["resp_nan_count"] = int(nan.sum())
    info["resp_nan_longest_s"] = _longest_run(nan) / fs
    if info["resp_nan_longest_s"] > 60:
        return None, f"RESP NaN gap {info['resp_nan_longest_s']:.0f} s > 60 s", info
    if nan.any():
        if nan.all():
            return None, "RESP all NaN", info
        idx = np.arange(len(resp))
        resp[nan] = np.interp(idx[nan], idx[~nan], resp[~nan])
    sel = np.array([s in beats for s in ann.symbol], bool)
    bt = ann.sample[sel] / afs
    rr = np.diff(bt)
    tt, hh, acc, rej = [], [], [], []
    n_ok = 0
    for k, v in enumerate(rr):
        ok = 0.3 <= v <= 2.0
        if ok and acc:
            ref = float(np.median(acc[-5:]))
            ok = abs(v - ref) / ref <= 0.20
        if ok:
            tt.append(bt[k + 1])
            hh.append(60.0 / v)
            acc.append(v)
            rej = []
            n_ok += 1
        else:
            rej.append(v)
            if len(rej) >= 5:  # stuck reference: reset to median of the 5 rejected RRs
                acc = [float(np.median(rej[-5:]))]
                rej = []
    info["rr_accepted_frac"] = n_ok / max(1, len(rr))
    tt, hh = np.array(tt), np.array(hh)
    if len(tt) < 2:
        return None, "too few accepted beats", info
    grid_hr = np.arange(np.ceil(tt[0] * FS), np.floor(tt[-1] * FS) + 1) / FS  # no extrapolation
    hr = np.interp(grid_hr, tt, hh)
    b, a = butter(4, 1.0, btype="low", fs=fs)
    fr = Fraction(FS).limit_denominator(10000) / Fraction(fs).limit_denominator(10000)
    rs = resample_poly(filtfilt(b, a, resp), fr.numerator, fr.denominator)
    grid_rs = np.arange(len(rs)) / FS
    lo, hi = max(grid_hr[0], grid_rs[0]), min(grid_hr[-1], grid_rs[-1])
    tg = np.arange(int(np.ceil(lo * FS - 1e-9)), int(np.floor(hi * FS + 1e-9)) + 1) / FS
    resp_c = np.interp(tg, grid_rs, rs)[:MAX_SAMPLES]
    hr_c = np.interp(tg, grid_hr, hr)[:MAX_SAMPLES]
    if len(resp_c) < MIN_SAMPLES:
        return None, f"{len(resp_c) / FS / 60:.1f} min < 60 min", info
    t_end = tg[len(resp_c) - 1]
    pts = np.concatenate([[tg[0]], tt[(tt > tg[0]) & (tt < t_end)], [t_end]])
    info["hr_gap_longest_s"] = float(np.max(np.diff(pts)))
    if info["hr_gap_longest_s"] > 60:
        return None, f"interpolated HR gap {info['hr_gap_longest_s']:.0f} s > 60 s", info
    assert np.isfinite(resp_c).all() and np.isfinite(hr_c).all()
    return (resp_c, hr_c), "", info


def age_sex(data_dir: str, rec: str) -> tuple[str, str]:
    txt = open(os.path.join(data_dir, rec + ".hea")).read()
    m = re.search(r"Age:\s*(\S+)\s+Sex:\s*(\S+)", txt)
    return (m.group(1), m.group(2)) if m else ("", "")


# ---------------------------------------------------------------- worker
def run_subject(args):
    rec, resp, hr, hr_next, next_rec, n_null = args
    system = [resp, hr]
    vp = variant_params(system, FS, MAX_WINDOW_S, seed=stable_seed(DATASET, rec))
    prm_rows, rows = [], []
    n = min(len(resp), len(hr_next))
    for v in VARIANTS:
        p, info = vp[v]
        prm_rows.append({"subject": rec, "variant": v, **params_row(p, FS),
                         "calibration_failed": bool(info.get("calibration_failed", False)),
                         "window_clipped": bool(info.get("window_clipped", False))})
        tests = {"within": (resp, hr, rec), "cross": (resp[:n], hr_next[:n], f"{rec}>{next_rec}")}
        for test, (x, y, tag) in tests.items():
            res = pair_test(x, y, p, n_null=n_null, seed=stable_seed(DATASET, rec, test))
            rows.append({"subject": rec, "variant": v, "test": test, "pair": tag, "score": res["score"],
                         "chance95": res["chance95"], "excess": res["score"] - res["chance95"], "p": res["p"],
                         "lag_s": -res["lag_samples"] / FS, "n_windows": res["n_windows"]})
    return prm_rows, rows


# ---------------------------------------------------------------- summary
def summarize(rows: pd.DataFrame, prm: pd.DataFrame) -> pd.DataFrame:
    out = []
    for v in VARIANTS:
        r = rows[rows.variant == v]
        w = r[r.test == "within"].set_index("subject")
        c = r[r.test == "cross"].set_index("subject")
        pv = prm[prm.variant == v]
        det = w[w.p < 0.05]
        young = w[w.index.str[2] == "y"]
        old = w[w.index.str[2] == "o"]
        ge_x = group_effect(young.excess, old.excess)
        ge_s = group_effect(young.score, old.score)
        lag = det.lag_s.dropna()
        out.append({
            "variant": v, "n_subjects": len(w), "detect": float((w.p < 0.05).mean()),
            "xsubj_fpr": float((c.p < 0.05).mean()),
            "auc_score": auc(w.score, c.score), "auc_excess": auc(w.excess, c.excess),
            "yo_excess_d": ge_x["d"], "yo_excess_auc": ge_x["auc"],
            "yo_score_d": ge_s["d"], "yo_score_auc": ge_s["auc"],
            "young_excess_mean": ge_x.get("mean_a", np.nan), "old_excess_mean": ge_x.get("mean_b", np.nan),
            "n_young": ge_x["n_a"], "n_old": ge_x["n_b"],
            "lag_s_median": float(lag.median()) if len(lag) else np.nan,
            "lag_s_q25": float(lag.quantile(.25)) if len(lag) else np.nan,
            "lag_s_q75": float(lag.quantile(.75)) if len(lag) else np.nan,
            "n_lag": len(lag),
            "window_s_median": float(pv.window_s.median()), "max_lag_s_median": float(pv.max_lag_s.median()),
            "tolerance_s_median": float(pv.tolerance_s.median()),
            "calibration_failed": int(pv.calibration_failed.sum()),
        })
    return pd.DataFrame(out)


def select(records: list[str], n: int | None) -> list[str]:
    records = sorted(records)
    if n is None or n >= len(records):
        return records
    y = [r for r in records if r[2] == "y"]
    o = [r for r in records if r[2] == "o"]
    pick = []
    for k in range(n):
        src = y if k % 2 == 0 else o
        pick.append(src[k // 2])
    return sorted(pick)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/home/user/data_gt/fantasia")
    ap.add_argument("--out", default=os.path.join(HERE, "results"))
    ap.add_argument("--n-subjects", type=int, default=None)
    ap.add_argument("--n-null", type=int, default=100)
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--smoke", action="store_true")
    a = ap.parse_args()
    if a.smoke:
        a.n_subjects, a.n_null = 3, 20
        a.out = os.path.join(HERE, "results", "smoke")
    os.makedirs(a.out, exist_ok=True)
    t_start = time.time()

    recs = select(RECORDS, a.n_subjects)
    prov = download(a.data_dir, recs)
    beats = beat_symbols()
    data, skipped, pre = {}, {}, {}
    for rec in recs:
        d, why, pinfo = preprocess(a.data_dir, rec, beats)
        pre[rec] = pinfo
        if d is None:
            skipped[rec] = why
            print(f"skip {rec}: {why}")
        else:
            data[rec] = d
    names = sorted(data)
    jobs = []
    for k, rec in enumerate(names):
        nxt = names[(k + 1) % len(names)]
        jobs.append((rec, data[rec][0], data[rec][1], data[nxt][1], nxt, a.n_null))
    prm_rows, rows = [], []
    with ProcessPoolExecutor(max_workers=a.jobs) as ex:
        for pr, rr in ex.map(run_subject, jobs):
            prm_rows += pr
            rows += rr
            print(f"done {pr[0]['subject']}  ({time.time() - t_start:.0f}s)", flush=True)
    prm, rows = pd.DataFrame(prm_rows), pd.DataFrame(rows)
    summ = summarize(rows, prm)
    prm.to_csv(os.path.join(a.out, "fantasia_params.csv"), index=False)
    rows.to_csv(os.path.join(a.out, "fantasia_rows.csv"), index=False)
    summ.to_csv(os.path.join(a.out, "fantasia_summary.csv"), index=False)
    wall = time.time() - t_start
    meta = {"dataset": DATASET, "fs": FS, "max_window_s": MAX_WINDOW_S, "n_null": a.n_null,
            "subjects": names, "skipped": skipped, "preprocessing": pre, "wall_time_s": wall,
            "age_sex": {r: age_sex(a.data_dir, r) for r in names},
            "files": prov, "total_bytes": int(sum(prov.values())), "source": f"{BASE}/{VERSION}/"}
    with open(os.path.join(a.out, "fantasia_meta.json"), "w") as f:
        json.dump(meta, f, indent=1)
    with pd.option_context("display.width", 250, "display.max_columns", 50, "display.float_format", "{:.3f}".format):
        print(summ.to_string(index=False))
    print(f"wall {wall:.0f}s, downloaded/available {sum(prov.values()) / 1e6:.1f} MB, subjects {len(names)}")


if __name__ == "__main__":
    main()
