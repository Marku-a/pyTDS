"""New physiological test: ECG -> finger pulse (PPG) coupling in the BIDMC PPG and Respiration dataset.

PhysioNet `bidmc` 1.0.0 (Pimentel et al., IEEE TBME 2017; Goldberger et al., Circulation 2000):
53 ICU recordings of 8 min at 125 Hz with ECG lead II and the finger photoplethysmogram (PLETH).
The pulse at the finger follows each heartbeat by the pulse arrival time (PAT, typically 0.2-0.6 s),
so the two signals are coupled with a known physiological delay. TDS has not been applied to these data before.

Design (fixed before scoring): ECG band-passed 0.5-30 Hz, PPG 0.5-8 Hz (2nd-order Butterworth, zero phase).
Coupled pairs: ECG and PPG of the same patient (53). Uncoupled controls: ECG of patient i with PPG of patient
i+17 (mod 53), so every PPG is used once. Settings: published (60/30/+-30/+-1 samples) and aTDS v2.1 set on the
pair itself. Reference delay: median beat-by-beat PAT from the R peak to the steepest PPG upstroke, and the
full-record cross-correlation peak within half a beat (the waveform alignment). Because the heartbeat repeats,
a delay is only defined modulo one beat period; delays are compared modulo the median RR interval.

    python experiments/atds_real/bidmc.py <bidmc_csv_dir> [out_dir] [--common]

<bidmc_csv_dir> holds bidmc_XX_Signals.csv from https://physionet.org/files/bidmc/1.0.0/bidmc_csv/
TDS via this repo's pyTDS only (through atds21).
"""
import csv
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.signal import butter, filtfilt, find_peaks

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import atds21 as A  # noqa: E402

FS = 125
IDS = list(range(1, 54))
SHIFT = 17
DIR = None


def bp(x, lo, hi):
    b, a = butter(2, [lo / (FS / 2), hi / (FS / 2)], "band")
    return filtfilt(b, a, x)


def load(d, i):
    df = pd.read_csv(os.path.join(d, f"bidmc_{i:02d}_Signals.csv"))
    df.columns = [c.strip() for c in df.columns]
    return A.interp_nan(df["II"].values), A.interp_nan(df["PLETH"].values)


def reference(ecg, ppg):
    """Median beat-by-beat PAT (samples) and median RR interval (samples)."""
    e = bp(ecg, 5, 30)
    if np.abs(e.min()) > e.max():
        e = -e
    r, _ = find_peaks(e, distance=int(0.3 * FS), height=np.percentile(e, 95) * 0.5)
    d = np.gradient(bp(ppg, 0.5, 8))
    pat = [int(np.argmax(d[a:b])) for a, b in zip(r[:-1], r[1:]) if b - a > 10]
    return float(np.median(pat)), float(np.median(np.diff(r))), int(len(r))


def xcorr_delay(x, y, period):
    """Full-record cross-correlation peak within half a beat, in pyTDS sign convention (negative = y lags)."""
    x = (x - x.mean()) / x.std(); y = (y - y.mean()) / y.std(); h = int(period // 2)
    cc = [np.mean(x[max(0, -k):len(x) - max(0, k)] * y[max(0, k):len(y) - max(0, -k)]) for k in range(-h, h + 1)]
    return float(-np.arange(-h, h + 1)[int(np.argmax(np.abs(cc)))])


def wrap(d, period):
    """Difference d wrapped into [-period/2, period/2)."""
    return (d + period / 2) % period - period / 2


def pair(d, i, j=None):
    """Filtered (ECG_i, PPG_j) as a dict usable by the plotting code."""
    j = i if j is None else j
    ex, _ = load(d, i); _, py = load(d, j)
    x, y = bp(ex, 0.5, 30), bp(py, 0.5, 8); n = min(len(x), len(y))
    out = dict(s1=x[:n], s2=y[:n], dt=1.0 / FS, names=("ECG", "finger pulse"), expected=None)
    if i == j:
        pat, rr, _ = reference(ex, py)
        e = wrap(-pat, rr)
        out["expected"] = ("delay", -e)
    return out


def run(p, x, y):
    tau = A.time_delay_interaction(x, y, p)[0]
    lab = A.stable_label(tau, p)
    return float(A.tds_score(lab)), A.delay_of(tau, lab), int(np.sum(lab)), int(len(tau))


def job(ij):
    i, j = ij
    ex, _ = load(DIR, i); _, py = load(DIR, j)
    x, y = bp(ex, 0.5, 30), bp(py, 0.5, 8); n = min(len(x), len(y)); x, y = x[:n], y[:n]
    pa = A.window_params([x, y]); pa, fail = A.calib(pa, [(x, y)])
    o = {"ecg": i, "ppg": j, "coupled": i == j, "n": n, "atds_window": pa.window, "atds_step": pa.overlap,
         "atds_max_lag": pa.max_lag, "atds_tolerance": pa.tolerance, "calibration_failed": bool(fail)}
    for name, p in (("published", A.PUBLISHED), ("atds", pa)):
        s, dl, ns, nw = run(p, x, y)
        o.update({f"{name}_score": s, f"{name}_delay": dl, f"{name}_n_stable": ns, f"{name}_n_windows": nw})
    if i == j:
        pat, rr, nb = reference(ex, py)
        o.update(pat_samples=pat, rr_samples=rr, n_beats=nb, xcorr_delay=xcorr_delay(x, y, rr))
        for name in ("published", "atds"):
            dl = o[f"{name}_delay"]
            o[f"{name}_err_pat_s"] = None if dl is None else float(wrap(dl - (-pat), rr) / FS)
            o[f"{name}_err_xcorr_s"] = None if dl is None else float(wrap(dl - o["xcorr_delay"], rr) / FS)
    return o


def auc(a, b):
    a, b = np.asarray(a)[:, None], np.asarray(b)[None, :]
    return float((np.sum(a > b) + 0.5 * np.sum(a == b)) / (a.size * b.size))


def boot_auc(a, b, n=5000, seed=0):
    rng = np.random.default_rng(seed); a, b = np.asarray(a), np.asarray(b)
    v = [auc(rng.choice(a, a.size), rng.choice(b, b.size)) for _ in range(n)]
    return [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]


def summarize(rows):
    c = [r for r in rows if r["coupled"]]; u = [r for r in rows if not r["coupled"]]
    S = {"n_coupled": len(c), "n_controls": len(u),
         "median_pat_s": float(np.median([r["pat_samples"] for r in c]) / FS),
         "median_rr_s": float(np.median([r["rr_samples"] for r in c]) / FS),
         "atds_window_s_median": float(np.median([r["atds_window"] for r in c]) / FS),
         "atds_window_s_range": [min(r["atds_window"] for r in c) / FS, max(r["atds_window"] for r in c) / FS],
         "atds_max_lag_s_median": float(np.median([r["atds_max_lag"] for r in c]) / FS),
         "atds_tolerance_samples_median": float(np.median([r["atds_tolerance"] for r in c]))}
    for m in ("published", "atds"):
        a = [r[f"{m}_score"] for r in c]; b = [r[f"{m}_score"] for r in u]; thr = float(np.percentile(b, 95))
        ep = [r[f"{m}_err_pat_s"] for r in c if r[f"{m}_err_pat_s"] is not None]
        ex = [r[f"{m}_err_xcorr_s"] for r in c if r[f"{m}_err_xcorr_s"] is not None]
        S[m] = {"coupled_mean": float(np.mean(a)), "coupled_median": float(np.median(a)),
                "control_mean": float(np.mean(b)), "control_95pct": thr,
                "auc": auc(a, b), "auc_ci95": boot_auc(a, b),
                "n_coupled_above_control95": int(sum(x > thr for x in a)),
                "n_with_delay": len(ep),
                "err_vs_pat_s_median": float(np.median(ep)) if ep else None,
                "err_vs_pat_s_iqr": [float(np.percentile(ep, 25)), float(np.percentile(ep, 75))] if ep else None,
                "n_within_0p1s_of_pat": int(sum(abs(e) <= 0.1 for e in ep)),
                "err_vs_xcorr_s_median": float(np.median(ex)) if ex else None,
                "n_within_0p04s_of_xcorr": int(sum(abs(e) <= 0.04 for e in ex))}
    return S


def main(d, out_dir):
    global DIR
    DIR = d
    os.makedirs(out_dir, exist_ok=True)
    t0 = time.time()
    jobs = [(i, i) for i in IDS] + [(i, IDS[(k + SHIFT) % len(IDS)]) for k, i in enumerate(IDS)]
    with Pool(A.WORKERS) as pool:
        rows = pool.map(job, jobs)
    S = summarize(rows); S["wall_seconds"] = time.time() - t0
    json.dump({"summary": S, "records": rows}, open(os.path.join(out_dir, "bidmc_ecg_ppg.json"), "w"), indent=1)
    keys = sorted({k for r in rows for k in r}, key=lambda k: list(rows[0]).index(k) if k in rows[0] else 99)
    with open(os.path.join(out_dir, "bidmc_ecg_ppg.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)
    print(json.dumps(S, indent=1))




# ---------------------------------------------------------------- one parameter set for the whole cohort
def _pair_signals(ij):
    i, j = ij
    ex, _ = load(DIR, i); _, py = load(DIR, j)
    x, y = bp(ex, 0.5, 30), bp(py, 0.5, 8); n = min(len(x), len(y))
    return x[:n], y[:n]


def _pair_window(i):
    x, y = _pair_signals((i, i))
    p = A.window_params([x, y])
    return p.window, p.max_lag


def _score_common(args):
    ij, p = args
    x, y = _pair_signals(ij)
    s, dl, ns, nw = run(p, x, y)
    return {"ecg": ij[0], "ppg": ij[1], "coupled": ij[0] == ij[1], "score": s, "delay": dl, "n_stable": ns}


def main_common(d, out_dir):
    """aTDS with ONE parameter set for all 106 pairs: window and lag = median over patients of the per-patient
    aTDS values (rhythm cap included); tolerance calibrated on rotations pooled over the true pairs (calib uses
    at most 50 pairs, seeded). Fixed TDS for reference. Delays compared with the per-patient reference of main()."""
    global DIR
    DIR = d
    t0 = time.time()
    with Pool(A.WORKERS) as pool:
        wl = pool.map(_pair_window, IDS)
    L = int(np.median([w for w, _ in wl])); lag = int(np.median([g for _, g in wl]))
    from dataclasses import replace
    p0 = replace(A.PUBLISHED, window=L, overlap=L // 2, max_lag=min(lag, (L - 1) // 2), tolerance=0)
    pairs = [_pair_signals((i, i)) for i in IDS]
    p, fail = A.calib(p0, pairs)
    del pairs
    jobs = [(i, i) for i in IDS] + [(i, IDS[(k + SHIFT) % len(IDS)]) for k, i in enumerate(IDS)]
    with Pool(A.WORKERS) as pool:
        rows = pool.map(_score_common, [(ij, p) for ij in jobs])
    ref = {r["ecg"]: r for r in json.load(open(os.path.join(out_dir, "bidmc_ecg_ppg.json")))["records"] if r["coupled"]}
    c = [r for r in rows if r["coupled"]]; u = [r for r in rows if not r["coupled"]]
    a = [r["score"] for r in c]; b = [r["score"] for r in u]; thr = float(np.percentile(b, 95))
    ex = [wrap(r["delay"] - ref[r["ecg"]]["xcorr_delay"], ref[r["ecg"]]["rr_samples"]) / FS for r in c if r["delay"] is not None]
    S = {"params": {"window": p.window, "step": p.overlap, "max_lag": p.max_lag, "tolerance": p.tolerance,
                    "calibration_failed": bool(fail)},
         "coupled_mean": float(np.mean(a)), "coupled_median": float(np.median(a)), "control_mean": float(np.mean(b)),
         "control_95pct": thr, "auc": auc(a, b), "auc_ci95": boot_auc(a, b),
         "n_coupled_above_control95": int(sum(x > thr for x in a)), "n_with_delay": len(ex),
         "n_within_0p04s_of_xcorr": int(sum(abs(e) <= 0.04 for e in ex)),
         "err_vs_xcorr_s_median": float(np.median(ex)) if ex else None, "wall_seconds": time.time() - t0}
    json.dump({"summary": S, "records": rows}, open(os.path.join(out_dir, "bidmc_common_params.json"), "w"), indent=1)
    print(json.dumps(S, indent=1))


if __name__ == "__main__":
    out = sys.argv[2] if len(sys.argv) > 2 and not sys.argv[2].startswith("--") else os.path.join(HERE, "results")
    if "--common" in sys.argv:
        main_common(sys.argv[1], out)
    else:
        main(sys.argv[1], out)
