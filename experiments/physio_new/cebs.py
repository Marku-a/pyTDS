"""CEBS (PhysioNet cebsdb 1.0.0): ECG -> seismocardiogram (SCG) with aTDS.

    python experiments/physio_new/cebs.py <raw_dir> [out_dir]

<raw_dir> holds b001..b020 and p001..p020 (.hea/.dat) from https://physionet.org/files/cebsdb/1.0.0/
See SPEC_COMMON.md / SPEC_cebs.md. TDS via this repo's pyTDS only (through experiments/atds_real/atds21.py).
Pair: x = ECG II (0.5-30 Hz), y = SCG (1-30 Hz), both at 250 Hz, first 270 s. Reported delay = -tau0 (positive = ECG leads).
Units = 40 records (20 basal + 20 post-music); controls pair x of record i with y of record i+13 (mod 40; never the same subject).
Reference: median over beats of time from R peak (detected on ECG II) to max |SCG| in 20-200 ms after R.
"""
import csv
import json
import os
import sys
import time
from dataclasses import replace
from multiprocessing import Pool

import numpy as np
from scipy.signal import butter, find_peaks, sosfiltfilt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "atds_real"))
import atds21 as A  # noqa: E402
from bidmc import auc, boot_auc, wrap, _first_below_e  # noqa: E402  (helpers only; unchanged)

FS = 250
FS_OUT = 250
FS_RAW = 5000
SECONDS = 270
IDS = [f"{g}{i:03d}" for g in "bp" for i in range(1, 21)]
WORKERS = 1
A.WORKERS = 1
DIR = None
SIG = {}
REF_TOL_S = 0.1       # "within tolerance" of the reference delay (assumption: the spec states none)
TOL2_S = 0.05


# ------------------------------------------------------------------ WFDB reader (formats 16 and 212, single segment)
def read_wfdb(d, rec, n_frames):
    lines = [l.split() for l in open(os.path.join(d, rec + ".hea")).read().splitlines() if l and not l.startswith("#")]
    nsig = int(lines[0][1]); sigs = lines[1:1 + nsig]
    fmt = sigs[0][1].split("x")[0].split(":")[0]
    assert all(s[1].split("x")[0].split(":")[0] == fmt for s in sigs), "mixed formats"
    names = [s[-1] for s in sigs]
    path = os.path.join(d, sigs[0][0])
    if fmt == "16":
        raw = np.fromfile(path, "<i2", count=n_frames * nsig).reshape(-1, nsig).astype(float)
    elif fmt == "212":
        nb = int(np.ceil(n_frames * nsig * 1.5)) + 3
        b = np.fromfile(path, np.uint8, count=nb).astype(np.int32)
        b = b[: len(b) // 3 * 3].reshape(-1, 3)
        s0 = b[:, 0] | ((b[:, 1] & 0x0F) << 8); s1 = b[:, 2] | ((b[:, 1] & 0xF0) << 4)
        v = np.stack([s0, s1], 1).reshape(-1)
        v = np.where(v >= 2048, v - 4096, v)[: n_frames * nsig]
        raw = v.reshape(-1, nsig).astype(float)
    else:
        raise ValueError(fmt)
    return {n: raw[:, k] for k, n in enumerate(names)}


# ------------------------------------------------------------------ preprocessing
def zs(x):
    return (x - x.mean()) / x.std()


def preprocess(d, rec):
    n = FS_RAW * SECONDS
    s = read_wfdb(d, rec, n)
    ecg, scg = s["II"], s["SCG"]
    assert len(ecg) == n, (rec, len(ecg))
    e = sosfiltfilt(butter(2, [0.5, 30], "band", fs=FS_RAW, output="sos"), ecg - ecg.mean())
    c = sosfiltfilt(butter(2, [1, 30], "band", fs=FS_RAW, output="sos"), scg - scg.mean())
    q = FS_RAW // FS
    e250, c250 = e[::q], c[::q]
    # R peaks on ECG II (polarity independent): 5-15 Hz band-pass, squared derivative energy, peaks >= 0.3 s apart (250 Hz)
    f = sosfiltfilt(butter(2, [5, 15], "band", fs=FS, output="sos"), e250)
    env = np.convolve(np.diff(f, prepend=f[0]) ** 2, np.ones(7) / 7, "same")
    pk, h = find_peaks(env, distance=int(0.3 * FS), height=np.percentile(env, 99) * 0.02)
    pk = pk[h["peak_heights"] > 0.3 * np.median(h["peak_heights"])]
    # refine to the local extremum of |ECG| (R wave) within +-60 ms
    w = int(0.06 * FS)
    R = np.array([max(0, k - w) + int(np.argmax(np.abs(e250[max(0, k - w):k + w + 1]))) for k in pk])
    R = np.unique(R)
    return zs(e250), zs(c250), R, {"n_beats": int(len(R))}


def reference(ecg, scg, R):
    """Median over beats of the time (s) from R to max |SCG| in 20-200 ms after R; median RR (s)."""
    lo, hi = int(round(0.02 * FS)), int(round(0.2 * FS))
    d = [(int(np.argmax(np.abs(scg[k + lo:k + hi + 1]))) + lo) / FS for k in R if k + hi + 1 <= len(scg)]
    rr = np.diff(R) / FS
    return {"rr_s": float(np.median(rr)), "ref_delay_s": float(np.median(d)), "ref_delay_iqr_s": float(np.subtract(*np.percentile(d, [75, 25]))),
            "n_beats": int(len(d))}


# ------------------------------------------------------------------ scoring
def run(p, x, y):
    tau = A.time_delay_interaction(x, y, p)[0]
    lab = A.stable_label(tau, p)
    dl = A.delay_of(tau, lab)
    return float(A.tds_score(lab)), None if dl is None else -dl / FS_OUT


def xcorr(x, y, lag):
    n = min(len(x), len(y)); x, y = zs(x[:n]), zs(y[:n])
    ks = np.arange(-lag, lag + 1)
    cc = np.array([np.mean(x[max(0, -k):n - max(0, k)] * y[max(0, k):n - max(0, -k)]) for k in ks])
    k = int(ks[int(np.argmax(np.abs(cc)))])
    return float(np.max(np.abs(cc))), k / FS_OUT      # y[t+k] ~ x[t], k>0: first signal leads


def _load(rec):
    x, y, R, info = preprocess(DIR, rec)
    return rec, x, y, R, info


def _subject_params(rec):
    x, y = SIG[rec]
    p = A.window_params([x, y])
    per = [A.dominant_period(s) for s in (x, y)]
    per = [q for q in per if q is not None]
    tc = max(_first_below_e(x), _first_below_e(y))
    return {"id": rec, "window": p.window, "max_lag": p.max_lag, "cap_applied": bool(p.max_lag < (p.window - 1) // 2),
            "min_period_samples": min(per) if per else None, "tc": tc, "r_hat": max(1, int(round(tc / 2)))}


def _score(args):
    (i, j), settings, lag = args
    x, y = SIG[IDS[i]][0], SIG[IDS[j]][1]
    n = min(len(x), len(y)); x, y = x[:n], y[:n]
    o = {"x": IDS[i], "y": IDS[j], "true": i == j}
    for name, p in settings.items():
        o[f"{name}_score"], o[f"{name}_delay_s"] = run(p, x, y)
    o["xcorr_score"], o["xcorr_delay_s"] = xcorr(x, y, lag)
    return o


def pd_(p):
    return {"window": p.window, "step": p.overlap, "max_lag": p.max_lag, "tolerance": p.tolerance,
            "window_s": p.window / FS_OUT, "max_lag_s": p.max_lag / FS_OUT, "tolerance_s": p.tolerance / FS_OUT}


def metrics(rows, m, ref, cond=None):
    sel = (lambda r: True) if cond is None else (lambda r: r["x"][0] == cond)
    c = [r for r in rows if r["true"] and sel(r)]; u = [r for r in rows if not r["true"] and sel(r)]
    a = [r[f"{m}_score"] for r in c]; b = [r[f"{m}_score"] for r in u]; thr = float(np.percentile(b, 95))
    S = {"true_mean": float(np.mean(a)), "true_median": float(np.median(a)), "control_mean": float(np.mean(b)),
         "control_median": float(np.median(b)), "control_95pct": thr, "auc": auc(a, b), "auc_ci95": boot_auc(a, b),
         "n_true_above_control95": int(sum(v > thr for v in a)), "n_true": len(a)}
    e = []
    for r in c:
        dl = r[f"{m}_delay_s"]
        if dl is None:
            continue
        e.append(float(wrap(dl - ref[r["x"]]["ref_delay_s"], ref[r["x"]]["rr_s"])))
    S["n_with_delay"] = len(e)
    S["err_vs_ref_s_median_abs"] = float(np.median(np.abs(e))) if e else None
    S["err_vs_ref_s_median_signed"] = float(np.median(e)) if e else None
    S[f"n_within_{REF_TOL_S}s_of_ref"] = int(sum(abs(v) <= REF_TOL_S for v in e))
    S[f"n_within_{TOL2_S}s_of_ref"] = int(sum(abs(v) <= TOL2_S for v in e))
    S["median_delay_s"] = float(np.median([r[f"{m}_delay_s"] for r in c if r[f"{m}_delay_s"] is not None])) if e else None
    return S


def main(d, out_dir):
    global DIR, SIG
    DIR = d
    os.makedirs(out_dir, exist_ok=True)
    t0 = time.time()
    with Pool(WORKERS) as pool:
        loaded = pool.map(_load, IDS)
    SIG = {r: (x, y) for r, x, y, _, _ in loaded}
    info = {r: i for r, _, _, _, i in loaded}
    ref = {r: dict(reference(SIG[r][0], SIG[r][1], R), **info[r]) for r, _, _, R, _ in loaded}
    n = len(IDS)
    with Pool(WORKERS) as pool:
        sp = pool.map(_subject_params, IDS)
    r_hat = int(np.median([s["r_hat"] for s in sp]))
    L = int(np.median([s["window"] for s in sp])); lag = int(np.median([s["max_lag"] for s in sp]))
    p0 = replace(A.PUBLISHED, window=L, overlap=L // 2, max_lag=min(lag, (L - 1) // 2), tolerance=0)
    atds, fail = A.calib(p0, [SIG[r] for r in IDS])
    settings = {"fixed": A.PUBLISHED,
                "rescaled": replace(A.PUBLISHED, window=60 * r_hat, overlap=30 * r_hat, max_lag=30 * r_hat, tolerance=r_hat),
                "atds": atds}
    k = 13
    assert np.gcd(k, n) == 1
    jobs = [(i, i) for i in range(n)] + [(i, (i + k) % n) for i in range(n)]
    with Pool(WORKERS) as pool:
        rows = pool.map(_score, [(ij, settings, atds.max_lag) for ij in jobs])
    methods = ["fixed", "rescaled", "atds", "xcorr"]
    rr = float(np.median([v["rr_s"] for v in ref.values()]))
    S = {"n_records": n, "control_shift_k": k, "r_hat_median": r_hat,
         "r_hat_range": [min(s["r_hat"] for s in sp), max(s["r_hat"] for s in sp)],
         "params": {m: pd_(p) for m, p in settings.items()}, "atds_calibration_failed": bool(fail),
         "xcorr_lag_samples": atds.max_lag,
         "reference": {"median_rr_s": rr, "median_ref_delay_s": float(np.median([v["ref_delay_s"] for v in ref.values()])),
                       "median_ref_delay_iqr_within_record_s": float(np.median([v["ref_delay_iqr_s"] for v in ref.values()]))},
         "reference_by_condition": {c: {"rr_s_median": float(np.median([ref[r]["rr_s"] for r in IDS if r[0] == c])),
                                        "ref_delay_s_median": float(np.median([ref[r]["ref_delay_s"] for r in IDS if r[0] == c]))}
                                    for c in "bp"},
         "atds_cap_applied_n_records": int(sum(s["cap_applied"] for s in sp)),
         "atds_lag_over_half_rr_median": float(np.median([s["max_lag"] / FS_OUT / (ref[s["id"]]["rr_s"] / 2) for s in sp])),
         "common_lag_over_median_rr": atds.max_lag / FS_OUT / rr,
         "rescaled_lag_over_median_rr": settings["rescaled"].max_lag / FS_OUT / rr,
         "fixed_lag_over_median_rr": A.PUBLISHED.max_lag / FS_OUT / rr,
         "metrics": {m: metrics(rows, m, ref) for m in methods},
         "metrics_basal": {m: metrics(rows, m, ref, "b") for m in methods},
         "metrics_post": {m: metrics(rows, m, ref, "p") for m in methods}}
    S["basal_vs_post"] = {}
    for m in methods:
        a = [r[f"{m}_score"] for r in rows if r["true"] and r["x"][0] == "b"]
        b = [r[f"{m}_score"] for r in rows if r["true"] and r["x"][0] == "p"]
        S["basal_vs_post"][m] = {"basal_mean": float(np.mean(a)), "post_mean": float(np.mean(b)),
                                 "auc_basal_gt_post": auc(a, b), "auc_ci95": boot_auc(a, b)}
    persub = {s["id"]: s for s in sp}
    for r in rows:
        if r["true"]:
            rf = ref[r["x"]]
            r.update(rr_s=rf["rr_s"], ref_delay_s=rf["ref_delay_s"], n_beats=rf["n_beats"],
                     persub_window=persub[r["x"]]["window"], persub_max_lag=persub[r["x"]]["max_lag"],
                     persub_cap_applied=persub[r["x"]]["cap_applied"], persub_tc=persub[r["x"]]["tc"])
    S["wall_seconds"] = time.time() - t0
    json.dump({"summary": S, "records": rows}, open(os.path.join(out_dir, "cebs.json"), "w"), indent=1)
    keys = list(dict.fromkeys(k_ for r in rows for k_ in r))
    with open(os.path.join(out_dir, "cebs.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)
    for title, key in (("ALL 40", "metrics"), ("BASAL 20", "metrics_basal"), ("POST 20", "metrics_post")):
        print(title)
        print(f"{'setting':9s} {'win/lag/tol (samples)':>22s} {'AUC [95% CI]':>22s} {'true mean':>9s} {'ctl95':>6s} {'>ctl95':>6s} {'delay':>5s} {'medAE s':>8s} {'<=0.1s':>6s} {'<=0.05':>6s}")
        for m in methods:
            q = S[key][m]
            ps = "%d/%d/%d" % (settings[m].window, settings[m].max_lag, settings[m].tolerance) if m in settings else f"lag {atds.max_lag}"
            nn = q["n_true"]
            print(f"{m:9s} {ps:>22s} {q['auc']:.3f} [{q['auc_ci95'][0]:.2f},{q['auc_ci95'][1]:.2f}] {q['true_mean']:9.1f} "
                  f"{q['control_95pct']:6.1f} {q['n_true_above_control95']:>3d}/{nn} {q['n_with_delay']:>5d} "
                  f"{q['err_vs_ref_s_median_abs'] if q['err_vs_ref_s_median_abs'] is not None else float('nan'):8.3f} "
                  f"{q[f'n_within_{REF_TOL_S}s_of_ref']:>6d} {q[f'n_within_{TOL2_S}s_of_ref']:>6d}")
    print(json.dumps({k_: S[k_] for k_ in ("params", "reference", "reference_by_condition", "basal_vs_post", "atds_cap_applied_n_records",
                                          "atds_lag_over_half_rr_median", "common_lag_over_median_rr", "r_hat_median", "r_hat_range",
                                          "atds_calibration_failed")}, indent=1))
    print("wall seconds", round(S["wall_seconds"], 1))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "results"))
