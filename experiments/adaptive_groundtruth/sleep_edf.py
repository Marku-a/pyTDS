"""
D2 Sleep-EDF Expanded: brain-rhythm network (10 EEG band-power signals) vs sleep stage.

Implements SPEC.md section D2. Night 1 of the first 20 subjects (sorted) of sleep-cassette. Per night,
10 log band-power signals at 1 Hz (delta/theta/alpha/sigma/beta x Fpz-Cz/Pz-Oz); 45 pairs; per-stage
real and circular-shift chance TDS scores from per-window stable labels.

Notes / deviations (see report):
- Hypnogram labels are mapped from the actual annotation descriptions: "Sleep stage W" -> W,
  "Sleep stage 1"/"Sleep stage 2" -> LS, "Sleep stage 3"/"Sleep stage 4" -> DS, "Sleep stage R" -> REM,
  anything else ("Sleep stage ?", "Movement time") -> NaN. Annotation onsets are taken as seconds from the
  start of the PSG recording (the pair shares its start time in this dataset).
- "non-W epoch" for the crop = any second labelled LS, DS or REM.
- Band edges are half-open [lo, hi) on the 0.5 Hz FFT bins of the 2-s window.
- Chance: window membership of the real windows is reused for the shifted-y windows.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (VARIANTS, paired_effect, params_row, stable_seed,  # noqa: E402
                    variant_params, window_labels)

DATASET = "sleep-edfx"
BASE = "https://physionet-open.s3.amazonaws.com"
PREFIX = "sleep-edfx/1.0.0/sleep-cassette/"
FS = 1.0
MAX_WINDOW_S = 300
N_SUBJECTS_DEFAULT = 20
CHANNELS = ["EEG Fpz-Cz", "EEG Pz-Oz"]
BANDS = [("delta", 0.5, 4), ("theta", 4, 8), ("alpha", 8, 12), ("sigma", 12, 16), ("beta", 16, 30)]
STAGES = ("W", "LS", "DS", "REM")
STAGE_MAP = {"Sleep stage W": "W", "Sleep stage 1": "LS", "Sleep stage 2": "LS",
             "Sleep stage 3": "DS", "Sleep stage 4": "DS", "Sleep stage R": "REM"}
HERE = os.path.dirname(os.path.abspath(__file__))


def list_keys() -> list[tuple[str, int]]:
    import urllib.parse
    keys, token = [], None
    while True:
        url = f"{BASE}/?list-type=2&prefix={PREFIX}"
        if token:
            url += "&continuation-token=" + urllib.parse.quote(token, safe="")
        xml = subprocess.run(["curl", "-sSf", "--retry", "3", url], capture_output=True, text=True, check=True).stdout
        keys += [(k, int(s)) for k, s in re.findall(r"<Key>([^<]*)</Key>.*?<Size>(\d+)</Size>", xml, flags=re.S)]
        m = re.search(r"<NextContinuationToken>([^<]*)</NextContinuationToken>", xml)
        if "<IsTruncated>true</IsTruncated>" in xml and m:
            token = m.group(1)
        else:
            return keys


def select_subjects(n: int | None) -> list[dict]:
    keys = dict(list_keys())
    psg = sorted(k.split("/")[-1] for k in keys if re.search(r"SC4\d\d1E0-PSG\.edf$", k))
    out = []
    for f in psg:
        tag = f[:6]
        hyp = sorted(k.split("/")[-1] for k in keys if k.split("/")[-1].startswith(tag) and k.endswith("-Hypnogram.edf"))
        if hyp:
            out.append({"subject": tag, "psg": f, "hyp": hyp[0], "psg_bytes": keys[PREFIX + f],
                        "hyp_bytes": keys[PREFIX + hyp[0]]})
    return out[:N_SUBJECTS_DEFAULT if n is None else n]


def download(name: str, dst: str) -> None:
    if os.path.exists(dst) and os.path.getsize(dst) > 0:
        return
    tmp = dst + ".part"
    subprocess.run(["curl", "-sSf", "--retry", "3", "-o", tmp, f"{BASE}/{PREFIX}{name}"], check=True)
    os.replace(tmp, dst)


def band_powers(x: np.ndarray, fs_in: float) -> np.ndarray:
    """(n_bands, n_sec) log10 band power; 2-s Hann windows stepped 1 s; column k = window starting at k s."""
    w = int(round(2 * fs_in))
    step = int(round(fs_in))
    win = np.hanning(w)
    freqs = np.fft.rfftfreq(w, 1 / fs_in)
    n = (len(x) - w) // step + 1
    out = np.empty((len(BANDS), n))
    for s in range(0, n, 4000):
        e = min(n, s + 4000)
        seg = np.lib.stride_tricks.sliding_window_view(x, w)[s * step:(e - 1) * step + 1:step]
        p = np.abs(np.fft.rfft(seg * win, axis=1)) ** 2
        for b, (_, lo, hi) in enumerate(BANDS):
            out[b, s:e] = np.log10(p[:, (freqs >= lo) & (freqs < hi)].sum(axis=1) + 1e-12)
    return out


def load_night(info: dict, data_dir: str) -> dict:
    """Cached band powers + per-second stage codes (0 W, 1 LS, 2 DS, 3 REM, -1 other) + provenance."""
    cache = os.path.join(data_dir, info["subject"] + ".npz")
    if not os.path.exists(cache):
        import mne
        mne.set_log_level("ERROR")
        psg, hyp = os.path.join(data_dir, info["psg"]), os.path.join(data_dir, info["hyp"])
        download(info["hyp"], hyp)
        download(info["psg"], psg)
        raw = mne.io.read_raw_edf(psg, include=CHANNELS, preload=True, verbose="ERROR")
        assert raw.ch_names == CHANNELS or sorted(raw.ch_names) == sorted(CHANNELS), raw.ch_names
        sig = []
        for ch in CHANNELS:
            x = raw.get_data(picks=[ch])[0] * 1e6
            sig.append(band_powers(x, raw.info["sfreq"]))
        bp = np.vstack(sig)  # (10, n_sec): fpz bands then pz bands
        n_sec = bp.shape[1]
        ann = mne.read_annotations(hyp)
        st = np.full(int(np.ceil(raw.times[-1])) + 1, -1, dtype=np.int8)
        descs = {}
        for on, du, de in zip(ann.onset, ann.duration, ann.description):
            descs[de] = descs.get(de, 0) + 1
            code = {"W": 0, "LS": 1, "DS": 2, "REM": 3}.get(STAGE_MAP.get(de), -1)
            st[int(round(on)):int(round(on + du))] = code
        np.savez(cache, bp=bp, stage=st[:n_sec], psg_bytes=os.path.getsize(psg), descs=json.dumps(descs))
        os.remove(psg)
    z = np.load(cache)
    return {"bp": z["bp"], "stage": z["stage"], "psg_bytes": int(z["psg_bytes"]), "descs": json.loads(str(z["descs"]))}


def build_system(night: dict) -> tuple[list[np.ndarray], np.ndarray, list[str]]:
    bp, st = night["bp"], night["stage"]
    n = min(bp.shape[1], len(st))
    bp, st = bp[:, :n], st[:n]
    nonw = np.flatnonzero(st >= 1)
    if len(nonw) == 0:
        raise ValueError("no sleep epochs")
    a, b = max(0, nonw[0] - 1800), min(n, nonw[-1] + 1 + 1800)
    names = [f"{ch.split()[-1].replace('-', '')}_{band}" for ch in CHANNELS for band, _, _ in BANDS]
    return [bp[i, a:b].astype(float) for i in range(bp.shape[0])], st[a:b], names


def window_stage(centres: np.ndarray, st: np.ndarray, W: int) -> np.ndarray:
    """Stage index (0..3 per STAGES) of each window if >= 80 % of its samples carry it, else -1."""
    out = np.full(len(centres), -1)
    for i, c in enumerate(centres):
        s = int(c) - W // 2
        seg = st[max(0, s):s + W]
        if len(seg) < W:
            continue
        for k in range(4):
            if np.mean(seg == k) >= 0.8:
                out[i] = k
                break
    return out


def stage_scores(lbl: np.ndarray, ws: np.ndarray) -> list[float]:
    return [100.0 * float(np.mean(lbl[ws == k])) if np.sum(ws == k) >= 5 else float("nan") for k in range(4)]


def run_subject(args) -> dict:
    info, data_dir, n_null_unused = args
    sub = info["subject"]
    night = load_night(info, data_dir)
    sigs, st, names = build_system(night)
    T = len(sigs[0])
    assert all(np.isfinite(s).all() for s in sigs)
    vp = variant_params(sigs, FS, MAX_WINDOW_S, seed=stable_seed(DATASET, sub))
    prow, rows = [], []
    for v in VARIANTS:
        p, inf = vp[v]
        prow.append({"subject": sub, "variant": v, **params_row(p, FS),
                     "calibration_failed": bool(inf.get("calibration_failed", False)),
                     "window_clipped": bool(inf.get("window_clipped", False))})
        W = p.window
        lo = max(2 * W, 5 * p.max_lag)
        if 2 * lo >= T:
            lo = T // 3
        for (i, j) in combinations(range(len(sigs)), 2):
            pair = f"{names[i]}|{names[j]}"
            c, lbl, _ = window_labels(sigs[i], sigs[j], p)
            ws = window_stage(c, st, W)
            real = stage_scores(lbl, ws)
            rng = np.random.default_rng(stable_seed(DATASET, sub, pair, "chance"))
            ch = []
            for _ in range(5):
                sh = int(rng.integers(lo, T - lo + 1))
                _, l0, _ = window_labels(sigs[i], np.roll(sigs[j], sh), p)
                ch.append(stage_scores(l0, ws))
            ch = np.nanmean(np.array(ch), axis=0) if np.isfinite(ch).any() else np.full(4, np.nan)
            for k, s in enumerate(STAGES):
                rows.append({"subject": sub, "variant": v, "pair": pair, "stage": s,
                             "n_windows": int(np.sum(ws == k)), "real": real[k], "chance": ch[k],
                             "excess": real[k] - ch[k]})
    return {"subject": sub, "params": prow, "rows": rows, "T": T, "psg_bytes": night["psg_bytes"],
            "descs": night["descs"]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/home/user/data_gt/sleep-edfx")
    ap.add_argument("--out", default=os.path.join(HERE, "results"))
    ap.add_argument("--n-subjects", type=int, default=None)
    ap.add_argument("--n-null", type=int, default=100)  # unused: D2 chance uses 5 shifts per spec
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--smoke", action="store_true")
    a = ap.parse_args()
    if a.smoke:
        a.n_subjects, a.n_null, a.out = 3, 20, os.path.join(HERE, "results", "smoke")
    os.makedirs(a.data_dir, exist_ok=True)
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    subs = select_subjects(a.n_subjects)
    with ProcessPoolExecutor(a.jobs) as ex:
        res = list(ex.map(run_subject, [(s, a.data_dir, a.n_null) for s in subs]))
    params = pd.DataFrame([r for x in res for r in x["params"]])
    rows = pd.DataFrame([r for x in res for r in x["rows"]])
    params.to_csv(os.path.join(a.out, "sleep_edf_params.csv"), index=False)
    rows.to_csv(os.path.join(a.out, "sleep_edf_rows.csv"), index=False)

    srows = []
    for v in VARIANTS:
        rv = rows[rows.variant == v]
        night = rv.groupby(["subject", "stage"])[["real", "chance", "excess", "n_windows"]].agg(
            {"real": "mean", "chance": "mean", "excess": "mean", "n_windows": "mean"}).reset_index()
        piv = {m: night.pivot(index="subject", columns="stage", values=m) for m in ("real", "excess", "n_windows")}
        d = {"variant": v, "n_subjects": len(subs)}
        for m in ("real", "excess"):
            for s1 in ("LS", "W", "REM"):
                e = paired_effect(piv[m][s1].values, piv[m]["DS"].values)
                d.update({f"{m}_{s1}_vs_DS_{k}": val for k, val in e.items()})
        for s in STAGES:
            d[f"mean_windows_{s}"] = float(piv["n_windows"][s].mean())
            d[f"mean_real_{s}"] = float(piv["real"][s].mean())
            d[f"mean_excess_{s}"] = float(piv["excess"][s].mean())
        pv = params[params.variant == v]
        for c in ("window_s", "max_lag_s", "tolerance_s"):
            d[f"median_{c}"] = float(pv[c].median())
        d["calibration_failed"] = int(pv.calibration_failed.sum())
        srows.append(d)
    summ = pd.DataFrame(srows)
    summ.to_csv(os.path.join(a.out, "sleep_edf_summary.csv"), index=False)
    wall = time.time() - t0
    descs: dict = {}
    for x in res:
        for k, c in x["descs"].items():
            descs[k] = descs.get(k, 0) + c
    meta = {"dataset": DATASET, "wall_time_s": wall, "jobs": a.jobs, "n_null": a.n_null,
            "max_window_s": MAX_WINDOW_S, "annotation_labels_seen": descs,
            "files": [{"psg": s["psg"], "psg_bytes": s["psg_bytes"], "hyp": s["hyp"], "hyp_bytes": s["hyp_bytes"]}
                      for s in subs],
            "samples_per_night": {x["subject"]: x["T"] for x in res}}
    with open(os.path.join(a.out, "sleep_edf_meta.json"), "w") as f:
        json.dump(meta, f, indent=1)

    pd.set_option("display.width", 250, "display.max_columns", 50)
    cols = ["variant", "median_window_s", "median_max_lag_s", "median_tolerance_s"]
    for m in ("real", "excess"):
        for s1 in ("LS", "W", "REM"):
            cols += [f"{m}_{s1}_vs_DS_mean_diff"]
    print(summ[cols].round(2).to_string(index=False))
    print(summ[["variant"] + [f"{m}_{s1}_vs_DS_{k}" for m in ("excess",) for s1 in ("LS", "W", "REM")
                              for k in ("dz", "frac_pos", "sign_p")]].round(3).to_string(index=False))
    print(summ[["variant"] + [f"mean_windows_{s}" for s in STAGES]].round(1).to_string(index=False))
    print("annotation labels:", descs)
    print(f"wall time {wall:.1f}s")


if __name__ == "__main__":
    main()
