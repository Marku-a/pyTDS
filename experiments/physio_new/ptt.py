"""Pulse transit time during exercise (PhysioNet pulse-transit-time-ppg 1.1.0): ECG -> finger PPG with aTDS.

    python experiments/physio_new/ptt.py <raw_dir> [out_dir]

<raw_dir> holds s{1..22}_{sit,walk,run}.hea/.dat/.atr from https://physionet.org/files/pulse-transit-time-ppg/1.1.0/
See SPEC_COMMON.md / SPEC_ptt.md. TDS via this repo's pyTDS only (through experiments/atds_real/atds21.py).
Pair: x = ECG (0.5-30 Hz), y = distal-phalanx infrared PPG pleth_2 (0.5-8 Hz), both 125 Hz.
Reported delay = -tau0 (positive = ECG leads). One worker process (no Pool).
"""
import csv
import json
import os
import sys
import time
from dataclasses import replace

import numpy as np
from scipy.signal import butter, sosfiltfilt
from scipy.stats import skew

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "atds_real"))
import atds21 as A  # noqa: E402
from bidmc import auc, boot_auc, wrap, _first_below_e  # noqa: E402  (helpers only; unchanged)

A.WORKERS = 1
FS_RAW = 500
FS = 125
NSUB = 22
ACTS = ["sit", "walk", "run"]
SUBS = list(range(1, NSUB + 1))
RECS = [(s, a) for a in ACTS for s in SUBS]       # ordered by activity, then subject
PPG_CH = "pleth_2"                                # infrared, distal phalanx (fingertip side)
SHIFT = 7                                         # gcd(7, 22) = 1: every PPG used once as control
REF_TOL_S = 0.1       # "within tolerance" of the reference delay (assumption: the spec states none)
TOL2_S = 0.05


def name(s, a):
    return f"s{s}_{a}"


# ------------------------------------------------------------------ readers
def read_wfdb(d, rec, want):
    """Format-212 single-segment reader returning only the channels in `want`."""
    lines = [l.split() for l in open(os.path.join(d, rec + ".hea")).read().splitlines() if l and not l.startswith("#")]
    nsig, nfr = int(lines[0][1]), int(lines[0][3])
    sigs = lines[1:1 + nsig]
    assert all(s[1].split("x")[0].split(":")[0] == "212" for s in sigs), "unexpected format"
    names = [s[-1] for s in sigs]
    b = np.fromfile(os.path.join(d, sigs[0][0]), np.uint8, count=int(np.ceil(nfr * nsig * 1.5)) + 3).astype(np.int32)
    b = b[: len(b) // 3 * 3].reshape(-1, 3)
    v = np.stack([b[:, 0] | ((b[:, 1] & 0x0F) << 8), b[:, 2] | ((b[:, 1] & 0xF0) << 4)], 1).reshape(-1)
    v = np.where(v >= 2048, v - 4096, v)[: nfr * nsig].reshape(-1, nsig)
    out = {}
    for w in want:
        k = names.index(w); s = sigs[k]
        gain = float(s[2].split("/")[0].split("(")[0]); base = int(s[2].split("(")[1].split(")")[0]) if "(" in s[2] else int(s[4])
        out[w] = (v[:, k] - base) / gain
    return out


def read_atr(d, rec):
    w = np.fromfile(os.path.join(d, rec + ".atr"), "<u2"); i = 0; t = 0; out = []
    while i < len(w):
        c, dt = int(w[i]) >> 10, int(w[i]) & 1023; i += 1
        if c == 0 and dt == 0:
            break
        if c == 59:
            t += (int(w[i]) << 16) + int(w[i + 1]); i += 2; continue
        if c in (60, 61, 62):
            continue
        if c == 63:
            i += (dt + 1) // 2; continue
        t += dt
        if c in (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 34, 35, 37, 38, 41):   # beat annotations
            out.append(t)
    return np.array(out) / FS_RAW


# ------------------------------------------------------------------ preprocessing / reference
def zs(x):
    return (x - x.mean()) / x.std()


def preprocess(d, rec):
    s = read_wfdb(d, rec, ["ecg", PPG_CH])
    ecg = sosfiltfilt(butter(2, [0.5, 30], "band", fs=FS_RAW, output="sos"), s["ecg"] - s["ecg"].mean())
    ppg = sosfiltfilt(butter(2, [0.5, 8], "band", fs=FS_RAW, output="sos"), s[PPG_CH] - s[PPG_CH].mean())
    q = FS_RAW // FS
    return zs(ecg[::q]), zs(ppg[::q]), ppg


def reference(d, rec, x, y, ppg500):
    """Median beat-by-beat PAT (R peak -> steepest PPG upstroke before the next R peak), median RR, and the
    full-record cross-correlation delay within half a beat."""
    r = read_atr(d, rec); rr = np.diff(r); rr_med = float(np.median(rr))
    dp = np.gradient(ppg500)
    sgn = 1.0 if skew(dp) > 0 else -1.0       # true PPG polarity: fast rise, slow fall -> positively skewed derivative
    dp = sgn * dp
    pat = []
    for k in range(len(r) - 1):
        if not (0.3 < rr[k] < 2.0):
            continue
        a, b = int((r[k] + 0.05) * FS_RAW), int((r[k] + rr[k]) * FS_RAW)
        if b >= len(dp):
            continue
        pat.append(a / FS_RAW - r[k] + np.argmax(dp[a:b]) / FS_RAW)
    lag = max(1, int(rr_med / 2 * FS))
    xs, xd = xcorr(x, y, lag)
    return {"rr_median_s": rr_med, "hr_bpm": 60.0 / rr_med, "n_beats": int(len(r)), "n_pat": len(pat),
            "pat_median_s": float(np.median(pat)), "pat_iqr_s": float(np.subtract(*np.percentile(pat, [75, 25]))),
            "ppg_polarity": sgn, "xcorr_ref_delay_s": xd, "xcorr_ref_score": xs}


# ------------------------------------------------------------------ scoring
def run(p, x, y):
    tau = A.time_delay_interaction(x, y, p)[0]
    lab = A.stable_label(tau, p)
    dl = A.delay_of(tau, lab)
    return float(A.tds_score(lab)), None if dl is None else -dl / FS


def xcorr(x, y, lag):
    """Normalised cross-correlation (mean of products over the overlap) via FFT, argmax |C| within +-lag samples."""
    n = min(len(x), len(y)); x, y = zs(x[:n]), zs(y[:n]); lag = min(lag, n // 2)
    m = 1 << int(np.ceil(np.log2(2 * n)))
    c = np.fft.irfft(np.conj(np.fft.rfft(x, m)) * np.fft.rfft(y, m), m)   # c[k] = sum_t x[t] y[t+k]
    ks = np.arange(-lag, lag + 1)
    cc = c[ks % m] / (n - np.abs(ks))
    i = int(np.argmax(np.abs(cc)))
    return float(np.abs(cc[i])), int(ks[i]) / FS


def pd_(p):
    return {"window": p.window, "step": p.overlap, "max_lag": p.max_lag, "tolerance": p.tolerance,
            "window_s": p.window / FS, "max_lag_s": p.max_lag / FS, "tolerance_s": p.tolerance / FS}


def common_set(sig, recs):
    """aTDS common set over `recs`: median window and lag of the per-record window_params, calibrated tolerance."""
    ps = [A.window_params(list(sig[r])) for r in recs]
    L = int(np.median([p.window for p in ps])); lag = int(np.median([p.max_lag for p in ps]))
    p0 = replace(A.PUBLISHED, window=L, overlap=L // 2, max_lag=min(lag, (L - 1) // 2), tolerance=0)
    return A.calib(p0, [sig[r] for r in recs])


def metrics(rows, m, ref):
    c = [r for r in rows if r["true"]]; u = [r for r in rows if not r["true"]]
    a = [r[f"{m}_score"] for r in c]; b = [r[f"{m}_score"] for r in u]; thr = float(np.percentile(b, 95))
    S = {"true_mean": float(np.mean(a)), "true_median": float(np.median(a)), "control_mean": float(np.mean(b)),
         "control_median": float(np.median(b)), "control_95pct": thr, "auc": auc(a, b), "auc_ci95": boot_auc(a, b),
         "n_true_above_control95": int(sum(v > thr for v in a)), "n_true": len(a)}
    def errs(rs):
        e = {"pat": [], "xcorr_ref": []}
        for r in rs:
            dl = r[f"{m}_delay_s"]
            if dl is None:
                continue
            rf = ref[(r["subject"], r["activity"])]
            e["pat"].append(float(wrap(dl - rf["pat_median_s"], rf["rr_median_s"])))
            e["xcorr_ref"].append(float(wrap(dl - rf["xcorr_ref_delay_s"], rf["rr_median_s"])))
        return e
    e = errs(c)
    S["n_with_delay"] = len(e["pat"])
    for k, v in e.items():
        S[f"err_vs_{k}_s_median_abs"] = float(np.median(np.abs(v))) if v else None
        S[f"n_within_{REF_TOL_S}s_of_{k}"] = int(sum(abs(q) <= REF_TOL_S for q in v))
        S[f"n_within_{TOL2_S}s_of_{k}"] = int(sum(abs(q) <= TOL2_S for q in v))
    S["per_activity"] = {}
    for act in ACTS:
        ca = [r for r in c if r["activity"] == act]; ua = [r for r in u if r["activity"] == act]
        ea = errs(ca); dls = [r[f"{m}_delay_s"] for r in ca if r[f"{m}_delay_s"] is not None]
        aa = [r[f"{m}_score"] for r in ca]; bb = [r[f"{m}_score"] for r in ua]
        S["per_activity"][act] = {
            "true_mean": float(np.mean(aa)), "control_mean": float(np.mean(bb)), "control_95pct": float(np.percentile(bb, 95)),
            "auc": auc(aa, bb), "n_true_above_control95": int(sum(v > np.percentile(bb, 95) for v in aa)), "n_with_delay": len(dls),
            "found_delay_median_s": float(np.median(dls)) if dls else None,
            "ref_pat_median_s": float(np.median([ref[(r["subject"], act)]["pat_median_s"] for r in ca])),
            "ref_xcorr_median_s": float(np.median([ref[(r["subject"], act)]["xcorr_ref_delay_s"] for r in ca])),
            "ref_rr_median_s": float(np.median([ref[(r["subject"], act)]["rr_median_s"] for r in ca])),
            "err_vs_pat_s_median_abs": float(np.median(np.abs(ea["pat"]))) if ea["pat"] else None,
            "n_within_tol_of_pat": int(sum(abs(q) <= REF_TOL_S for q in ea["pat"]))}
    return S


def main(d, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    t0 = time.time()
    SIG, ref, ppg_ = {}, {}, {}
    for s, a in RECS:
        rec = name(s, a)
        x, y, p5 = preprocess(d, rec)
        SIG[(s, a)] = (x, y)
        ref[(s, a)] = reference(d, rec, x, y, p5)
    # rescaling factor and per-record parameters
    rh, tcs, pr = {}, {}, {}
    for k, (x, y) in SIG.items():
        tcs[k] = max(_first_below_e(x), _first_below_e(y)); rh[k] = max(1, int(round(tcs[k] / 2)))
        pr[k] = A.window_params([x, y])
    r_hat = int(np.median(list(rh.values())))
    settings = {"fixed": {a: A.PUBLISHED for a in ACTS},
                "rescaled": {a: replace(A.PUBLISHED, window=60 * r_hat, overlap=30 * r_hat, max_lag=30 * r_hat, tolerance=r_hat) for a in ACTS}}
    L = int(np.median([p.window for p in pr.values()])); lag = int(np.median([p.max_lag for p in pr.values()]))
    p0 = replace(A.PUBLISHED, window=L, overlap=L // 2, max_lag=min(lag, (L - 1) // 2), tolerance=0)
    p_all, fail_all = A.calib(p0, [SIG[k] for k in RECS])
    settings["atds"] = {a: p_all for a in ACTS}
    fails = {"atds": bool(fail_all)}
    settings["atds_act"] = {}
    for a in ACTS:
        pa, fa = common_set(SIG, [(s, a) for s in SUBS])
        settings["atds_act"][a] = pa; fails.setdefault("atds_act", {})[a] = bool(fa)
    print("settings", {m: {a: pd_(p) for a, p in v.items()} for m, v in settings.items()}, flush=True)
    xlag = p_all.max_lag
    rows = []
    for a in ACTS:
        for i, s in enumerate(SUBS):
            for ctl in (False, True):
                sy = s if not ctl else SUBS[(i + SHIFT) % NSUB]
                x, y = SIG[(s, a)][0], SIG[(sy, a)][1]
                n = min(len(x), len(y)); x, y = x[:n], y[:n]
                o = {"subject": s, "activity": a, "y_subject": sy, "true": not ctl}
                for m, ps in settings.items():
                    o[f"{m}_score"], o[f"{m}_delay_s"] = run(ps[a], x, y)
                o["xcorr_score"], o["xcorr_delay_s"] = xcorr(x, y, xlag)
                rows.append(o)
        print("scored", a, round(time.time() - t0), "s", flush=True)
    for r in rows:
        if r["true"]:
            k = (r["subject"], r["activity"]); rf = ref[k]
            r.update({f"ref_{q}": v for q, v in rf.items()}, persub_window=pr[k].window, persub_max_lag=pr[k].max_lag,
                     persub_tc=tcs[k], persub_r_hat=rh[k])
    methods = ["fixed", "rescaled", "atds", "atds_act", "xcorr"]
    S = {"n_records": len(RECS), "control_shift_k": SHIFT, "ppg_channel": PPG_CH, "r_hat_median": r_hat,
         "r_hat_range": [min(rh.values()), max(rh.values())],
         "r_hat_median_by_activity": {a: int(np.median([rh[(s, a)] for s in SUBS])) for a in ACTS},
         "params": {m: {a: pd_(p) for a, p in v.items()} for m, v in settings.items()},
         "calibration_failed": fails, "xcorr_lag_samples": xlag,
         "reference_by_activity": {a: {q: float(np.median([ref[(s, a)][q] for s in SUBS]))
                                       for q in ("rr_median_s", "pat_median_s", "xcorr_ref_delay_s")} for a in ACTS},
         "n_ppg_inverted": int(sum(v["ppg_polarity"] < 0 for v in ref.values())),
         "metrics": {m: metrics(rows, m, ref) for m in methods}}
    # sit -> run: does the found delay shorten as the reference does?
    sr = {"ref_pat": [], "ref_xcorr": []}
    for s in SUBS:
        sr["ref_pat"].append(ref[(s, "run")]["pat_median_s"] - ref[(s, "sit")]["pat_median_s"])
        sr["ref_xcorr"].append(ref[(s, "run")]["xcorr_ref_delay_s"] - ref[(s, "sit")]["xcorr_ref_delay_s"])
    S["sit_to_run"] = {"ref_pat_change_median_s": float(np.median(sr["ref_pat"])),
                       "ref_pat_n_shorter": int(sum(v < 0 for v in sr["ref_pat"])),
                       "ref_xcorr_change_median_s": float(np.median(sr["ref_xcorr"])),
                       "ref_xcorr_n_shorter": int(sum(v < 0 for v in sr["ref_xcorr"]))}
    for m in methods:
        ch = []
        for s in SUBS:
            ds = {r["activity"]: r[f"{m}_delay_s"] for r in rows if r["true"] and r["subject"] == s and r["activity"] in ("sit", "run")}
            if ds.get("sit") is not None and ds.get("run") is not None:
                ch.append(ds["run"] - ds["sit"])
        S["sit_to_run"][m] = {"n_subjects_both": len(ch), "change_median_s": float(np.median(ch)) if ch else None,
                              "n_shorter": int(sum(v < 0 for v in ch))}
    S["wall_seconds"] = time.time() - t0
    json.dump({"summary": S, "records": rows}, open(os.path.join(out_dir, "ptt.json"), "w"), indent=1)
    keys = list(dict.fromkeys(k_ for r in rows for k_ in r))
    with open(os.path.join(out_dir, "ptt.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)
    n = len(RECS)
    print(f"{'setting':9s} {'win/lag/tol (samples)':>22s} {'AUC [95% CI]':>22s} {'true mean':>9s} {'ctl95':>6s} {'>ctl95':>7s} {'delay':>5s} {'medAE s':>8s} {'<=tol':>5s}")
    for m in methods:
        q = S["metrics"][m]
        if m == "xcorr":
            ps = f"lag {xlag}"
        elif m == "atds_act":
            ps = "/".join("%d,%d,%d" % (settings[m][a].window, settings[m][a].max_lag, settings[m][a].tolerance) for a in ACTS)
        else:
            p = settings[m]["sit"]; ps = "%d/%d/%d" % (p.window, p.max_lag, p.tolerance)
        e = q["err_vs_pat_s_median_abs"]
        print(f"{m:9s} {ps:>22s} {q['auc']:.3f} [{q['auc_ci95'][0]:.2f},{q['auc_ci95'][1]:.2f}] {q['true_mean']:9.1f} "
              f"{q['control_95pct']:6.1f} {q['n_true_above_control95']:>4d}/{n} {q['n_with_delay']:>5d} "
              f"{e if e is not None else float('nan'):8.3f} {q[f'n_within_{REF_TOL_S}s_of_pat']:>5d}")
    print(json.dumps({k_: S[k_] for k_ in ("reference_by_activity", "sit_to_run", "r_hat_median", "r_hat_median_by_activity",
                                          "calibration_failed", "n_ppg_inverted")}, indent=1))
    for m in methods:
        print(m, {a: {k_: (round(v, 3) if isinstance(v, float) else v) for k_, v in q.items() if k_ in
                      ("auc", "true_mean", "n_with_delay", "found_delay_median_s", "err_vs_pat_s_median_abs", "n_within_tol_of_pat")}
                  for a, q in S["metrics"][m]["per_activity"].items()})
    print("wall seconds", round(S["wall_seconds"], 1))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "results"))
