"""Fantasia (PhysioNet fantasia 1.0.0): respiration -> heart rate (respiratory sinus arrhythmia) with aTDS.

    python experiments/physio_new/fantasia.py <raw_dir> [out_dir]

<raw_dir> holds the 40 records f{1,2}{y,o}{01..10}.hea/.dat from https://physionet.org/files/fantasia/1.0.0/
See SPEC_COMMON.md / SPEC_fantasia.md. TDS via this repo's pyTDS only (through experiments/atds_real/atds21.py).
Pair: x = respiration, y = heart rate (4 Hz, first 60 min). Reported delay = -tau0 (positive = respiration leads).
"""
import csv
import json
import os
import sys
import time
from dataclasses import replace
from multiprocessing import Pool

import numpy as np
from scipy.ndimage import uniform_filter1d
from scipy.signal import butter, find_peaks, sosfiltfilt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "atds_real"))
import atds21 as A  # noqa: E402
from bidmc import auc, boot_auc, wrap, _first_below_e  # noqa: E402  (helpers only; unchanged)

FS = 250
FS_OUT = 4
MINUTES = 60
IDS = [f"f{f}{g}{i:02d}" for g in "yo" for f in (1, 2) for i in range(1, 11)]
WORKERS = 2
DIR = None
SIG = {}
REF_TOL_S = 1.0       # "within tolerance" of the reference delay (assumption: the spec states none)
TOL2_S = 0.5


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
DETREND = float(os.environ.get("FANTASIA_DETREND_S", "0"))
TAG = f"_detrend{int(DETREND)}s" if DETREND > 0 else ""


def zs(x):
    return (x - x.mean()) / x.std()


def preprocess(d, rec):
    n = FS * 60 * MINUTES
    s = read_wfdb(d, rec, n)
    resp, ecg = s["RESP"], s["ECG"]
    assert len(resp) == n
    t_out = np.arange(0, MINUTES * 60, 1.0 / FS_OUT)
    t = np.arange(n) / FS
    # respiration: 0.05-1 Hz band-pass (zero phase), resampled to 4 Hz
    r = sosfiltfilt(butter(2, [0.05, 1.0], "band", fs=FS, output="sos"), resp - resp.mean())
    resp4 = np.interp(t_out, t, r)
    # R peaks: 5-15 Hz band-pass, squared derivative energy (polarity independent), peaks >= 0.4 s apart
    e = sosfiltfilt(butter(2, [5, 15], "band", fs=FS, output="sos"), ecg - ecg.mean())
    env = np.convolve(np.diff(e, prepend=e[0]) ** 2, np.ones(25) / 25, "same")
    pk, h = find_peaks(env, distance=int(0.4 * FS), height=np.percentile(env, 99) * 0.02)
    pk = pk[h["peak_heights"] > 0.3 * np.median(h["peak_heights"])]
    bt = pk / FS
    rr = np.diff(bt); tb = bt[1:]
    med = np.array([np.median(rr[max(0, k - 5):k + 6]) for k in range(len(rr))])
    ok = (rr > 0.3) & (rr < 2.5) & (np.abs(rr - med) < 0.3 * med)       # drop missed/extra beats
    hr = np.interp(t_out, tb[ok], 60.0 / rr[ok])
    if DETREND > 0:   # optional: remove trends slower than DETREND s from both signals (same for all settings)
        k = int(DETREND * FS_OUT)
        resp4 = resp4 - uniform_filter1d(resp4, k, mode="nearest")
        hr = hr - uniform_filter1d(hr, k, mode="nearest")
    return zs(resp4), zs(hr), {"n_beats": int(len(pk)), "frac_rr_rejected": float(1 - ok.mean())}


def reference(resp, hr):
    """Breath-locked HR response: latency (s) of the HR maximum after inspiration onset and the median breath period (s)."""
    tr, _ = find_peaks(-resp, distance=int(1.0 * FS_OUT), prominence=0.5)
    # parabolic refinement of trough position
    tt = []
    for k in tr:
        if 0 < k < len(resp) - 1:
            a, b, c = resp[k - 1], resp[k], resp[k + 1]; den = a - 2 * b + c
            tt.append((k + (0.5 * (a - c) / den if den != 0 else 0.0)) / FS_OUT)
    tt = np.array(tt); Tb = float(np.median(np.diff(tt)))
    off = np.arange(0, 1.5 * Tb, 0.05); t = np.arange(len(hr)) / FS_OUT
    seg = [np.interp(t0 + off, t, hr) for t0 in tt if t0 + 1.5 * Tb < t[-1]]
    seg = np.array([s - s.mean() for s in seg]); m = seg.mean(0)
    return {"breath_period_s": Tb, "n_breaths": int(len(seg)), "hr_max_latency_s": float(off[int(np.argmax(m))]),
            "rsa_amp_z": float(m.max() - m.min())}


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
    x, y, info = preprocess(DIR, rec)
    return rec, x, y, info


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


def metrics(rows, m, ref):
    c = [r for r in rows if r["true"]]; u = [r for r in rows if not r["true"]]
    a = [r[f"{m}_score"] for r in c]; b = [r[f"{m}_score"] for r in u]; thr = float(np.percentile(b, 95))
    S = {"true_mean": float(np.mean(a)), "true_median": float(np.median(a)), "control_mean": float(np.mean(b)),
         "control_median": float(np.median(b)), "control_95pct": thr, "auc": auc(a, b), "auc_ci95": boot_auc(a, b),
         "n_true_above_control95": int(sum(v > thr for v in a)), "n_true": len(a)}
    err = {"hr_latency": [], "xcorr": []}
    for r in c:
        dl = r[f"{m}_delay_s"]
        if dl is None:
            continue
        rf = ref[r["x"]]; Tb = rf["breath_period_s"]
        err["hr_latency"].append(float(wrap(dl - rf["hr_max_latency_s"], Tb)))
        err["xcorr"].append(float(wrap(dl - r["xcorr_delay_s"], Tb)))
    S["n_with_delay"] = len(err["hr_latency"])
    if m == "xcorr":
        S["n_with_delay"] = len(c)
        err = {"hr_latency": [float(wrap(r["xcorr_delay_s"] - ref[r["x"]]["hr_max_latency_s"], ref[r["x"]]["breath_period_s"])) for r in c],
               "xcorr": []}
    for k, e in err.items():
        S[f"err_vs_{k}_s_median_abs"] = float(np.median(np.abs(e))) if e else None
        S[f"n_within_{REF_TOL_S}s_of_{k}"] = int(sum(abs(v) <= REF_TOL_S for v in e))
        S[f"n_within_{TOL2_S}s_of_{k}"] = int(sum(abs(v) <= TOL2_S for v in e))
    return S


def main(d, out_dir):
    global DIR, SIG
    DIR = d
    os.makedirs(out_dir, exist_ok=True)
    t0 = time.time()
    with Pool(WORKERS) as pool:
        loaded = pool.map(_load, IDS)
    SIG = {r: (x, y) for r, x, y, _ in loaded}
    info = {r: i for r, _, _, i in loaded}
    ref = {r: reference(x, y) for r, (x, y) in SIG.items()}
    n = len(IDS)
    with Pool(WORKERS) as pool:
        sp = pool.map(_subject_params, IDS)
    # settings
    r_hat = int(np.median([s["r_hat"] for s in sp]))
    L = int(np.median([s["window"] for s in sp])); lag = int(np.median([s["max_lag"] for s in sp]))
    p0 = replace(A.PUBLISHED, window=L, overlap=L // 2, max_lag=min(lag, (L - 1) // 2), tolerance=0)
    atds, fail = A.calib(p0, [SIG[r] for r in IDS])
    settings = {"fixed": A.PUBLISHED,
                "rescaled": replace(A.PUBLISHED, window=60 * r_hat, overlap=30 * r_hat, max_lag=30 * r_hat, tolerance=r_hat),
                "atds": atds}
    k = n // 3
    assert np.gcd(k, n) == 1
    jobs = [(i, i) for i in range(n)] + [(i, (i + k) % n) for i in range(n)]
    with Pool(WORKERS) as pool:
        rows = pool.map(_score, [(ij, settings, atds.max_lag) for ij in jobs])
    methods = ["fixed", "rescaled", "atds", "xcorr"]
    S = {"n_subjects": n, "control_shift_k": k, "r_hat_median": r_hat,
         "r_hat_range": [min(s["r_hat"] for s in sp), max(s["r_hat"] for s in sp)],
         "params": {m: pd_(p) for m, p in settings.items()}, "atds_calibration_failed": bool(fail),
         "xcorr_lag_samples": atds.max_lag,
         "reference": {"median_breath_period_s": float(np.median([v["breath_period_s"] for v in ref.values()])),
                       "median_hr_max_latency_s": float(np.median([v["hr_max_latency_s"] for v in ref.values()]))},
         "atds_cap_applied_n_subjects": int(sum(s["cap_applied"] for s in sp)),
         "atds_lag_over_half_breath_period_median": float(np.median([s["max_lag"] / FS_OUT / (ref[s["id"]]["breath_period_s"] / 2) for s in sp])),
         "common_lag_s_over_median_breath_period": atds.max_lag / FS_OUT / float(np.median([v["breath_period_s"] for v in ref.values()])),
         "rescaled_lag_s_over_median_breath_period": settings["rescaled"].max_lag / FS_OUT / float(np.median([v["breath_period_s"] for v in ref.values()])),
         "fixed_lag_s_over_median_breath_period": A.PUBLISHED.max_lag / FS_OUT / float(np.median([v["breath_period_s"] for v in ref.values()])),
         "metrics": {m: metrics(rows, m, ref) for m in methods}}
    # young vs elderly (true pairs and controls pooled by x-subject group)
    grp = lambda r: "young" if "y" in r["x"][2] else "elderly"
    S["young_vs_elderly"] = {}
    for m in methods:
        thr = S["metrics"][m]["control_95pct"]
        S["young_vs_elderly"][m] = {g: {"n": sum(1 for r in rows if r["true"] and grp(r) == g),
                                        "true_mean": float(np.mean([r[f"{m}_score"] for r in rows if r["true"] and grp(r) == g])),
                                        "true_median": float(np.median([r[f"{m}_score"] for r in rows if r["true"] and grp(r) == g])),
                                        "n_above_control95": int(sum(r[f"{m}_score"] > thr for r in rows if r["true"] and grp(r) == g))}
                                    for g in ("young", "elderly")}
        a = [r[f"{m}_score"] for r in rows if r["true"] and grp(r) == "young"]
        b = [r[f"{m}_score"] for r in rows if r["true"] and grp(r) == "elderly"]
        S["young_vs_elderly"][m]["auc_young_gt_elderly"] = auc(a, b)
        S["young_vs_elderly"][m]["auc_young_gt_elderly_ci95"] = boot_auc(a, b)
    S["reference_by_group"] = {g: {"hr_latency_median_s": float(np.median([ref[r]["hr_max_latency_s"] for r in IDS if grp({"x": r}) == g])),
                                   "breath_period_median_s": float(np.median([ref[r]["breath_period_s"] for r in IDS if grp({"x": r}) == g])),
                                   "rsa_amp_z_median": float(np.median([ref[r]["rsa_amp_z"] for r in IDS if grp({"x": r}) == g]))}
                               for g in ("young", "elderly")}
    # delay found compared modulo one breath period (done in metrics via wrap); also per-subject rows
    for r in rows:
        if r["true"]:
            rf = ref[r["x"]]
            r.update(breath_period_s=rf["breath_period_s"], ref_hr_max_latency_s=rf["hr_max_latency_s"], n_breaths=rf["n_breaths"],
                     rsa_amp_z=rf["rsa_amp_z"], **{f"info_{k_}": v for k_, v in info[r["x"]].items()})
    persub = {s["id"]: s for s in sp}
    for r in rows:
        if r["true"]:
            r.update(persub_window=persub[r["x"]]["window"], persub_max_lag=persub[r["x"]]["max_lag"],
                     persub_cap_applied=persub[r["x"]]["cap_applied"], persub_tc=persub[r["x"]]["tc"])
    S["wall_seconds"] = time.time() - t0
    json.dump({"summary": S, "records": rows}, open(os.path.join(out_dir, f"fantasia{TAG}.json"), "w"), indent=1)
    keys = list(dict.fromkeys(k_ for r in rows for k_ in r))
    with open(os.path.join(out_dir, f"fantasia{TAG}.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)
    # compact table
    print(f"{'setting':9s} {'win/lag/tol (samples)':>22s} {'AUC [95% CI]':>22s} {'true mean':>9s} {'ctl95':>6s} {'>ctl95':>6s} {'delay':>5s} {'medAE s':>8s} {'<=1s':>5s}")
    for m in methods:
        q = S["metrics"][m]
        ps = "%d/%d/%d" % (settings[m].window, settings[m].max_lag, settings[m].tolerance) if m in settings else f"lag {atds.max_lag}"
        print(f"{m:9s} {ps:>22s} {q['auc']:.3f} [{q['auc_ci95'][0]:.2f},{q['auc_ci95'][1]:.2f}] {q['true_mean']:9.1f} "
              f"{q['control_95pct']:6.1f} {q['n_true_above_control95']:>3d}/{n} {q['n_with_delay']:>5d} "
              f"{q['err_vs_hr_latency_s_median_abs'] if q['err_vs_hr_latency_s_median_abs'] is not None else float('nan'):8.2f} "
              f"{q[f'n_within_{REF_TOL_S}s_of_hr_latency']:>5d}")
    print(json.dumps({k_: S[k_] for k_ in ("reference", "reference_by_group", "young_vs_elderly", "atds_cap_applied_n_subjects",
                                          "atds_lag_over_half_breath_period_median", "common_lag_s_over_median_breath_period",
                                          "rescaled_lag_s_over_median_breath_period", "fixed_lag_s_over_median_breath_period",
                                          "r_hat_median", "r_hat_range", "atds_calibration_failed")}, indent=1))
    print("wall seconds", round(S["wall_seconds"], 1))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "results"))
