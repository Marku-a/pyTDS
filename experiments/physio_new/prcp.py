"""Posture change and baroreflex (PhysioNet prcp 1.0.0): SBP -> RR interval with aTDS.

    python experiments/physio_new/prcp.py <raw_dir> [out_dir]

<raw_dir> holds the 10 records NNNNN.{hea,dat,wabp,wqrs,anI} from https://physionet.org/files/prcp/1.0.0/
See SPEC_COMMON.md / SPEC_prcp.md. TDS via this repo's pyTDS only (through experiments/atds_real/atds21.py).
Pair: x = systolic BP, y = RR interval (4 Hz, whole record, 60-s moving average removed).
Reported delay = -tau0 (positive = SBP leads RR).
"""
import csv
import json
import os
import sys
import time
from dataclasses import replace
from multiprocessing import Pool

import numpy as np
from scipy.ndimage import median_filter, uniform_filter1d

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "atds_real"))
import atds21 as A  # noqa: E402
from bidmc import auc, boot_auc, _first_below_e  # noqa: E402  (helpers only; unchanged)

FS = 250
FS_OUT = 4
IDS = "12726 12734 12744 12754 12755 12814 12815 12819 12821 13960".split()
A.WORKERS = 1
WORKERS = 1
DIR = None
SIG = {}
EV = {}
REF_TOL_S = 1.0       # "within tolerance" of the reference delay (assumption: the spec states none)
TOL2_S = 0.5
PRE_S = 240.0         # "before standing" = this many seconds preceding a stand-up annotation


# ------------------------------------------------------------------ WFDB readers
def read_header(d, rec):
    lines = [l.split() for l in open(os.path.join(d, rec + ".hea")).read().splitlines() if l and not l.startswith("#")]
    nsig = int(lines[0][1]); sigs = lines[1:1 + nsig]
    names = [s[-1] for s in sigs]
    gains = []
    for s in sigs:
        g = s[2].split("/")[0]; base = 0.0
        if "(" in g:
            g, b = g.split("("); base = float(b.rstrip(")"))
        gains.append((float(g), base))
    assert all(s[1] == "16" for s in sigs)
    return nsig, names, gains, sigs[0][0]


def read_signal(d, rec, name):
    nsig, names, gains, fn = read_header(d, rec)
    raw = np.fromfile(os.path.join(d, fn), "<i2").reshape(-1, nsig).astype(float)
    k = names.index(name); g, b = gains[k]
    return (raw[:, k] - b) / g


def read_ann(path):
    """WFDB annotation file -> list of [sample, code, aux text] (sample units)."""
    w = np.fromfile(path, "<u2"); t = 0; out = []; i = 0
    while i < len(w):
        v = int(w[i]); i += 1; code = v >> 10; dt = v & 1023
        if v == 0:
            break
        if code == 59:                                   # SKIP: signed 32-bit interval, high word first
            dd = (int(w[i]) << 16) | int(w[i + 1]); dd = dd - (1 << 32) if dd >= 1 << 31 else dd
            i += 2; t += dd; continue
        if code in (60, 61, 62):
            continue
        if code == 63:                                   # AUX text of the previous annotation
            if out:
                out[-1][2] = bytes(w[i:i + (dt + 1) // 2].view(np.uint8)[:dt]).decode(errors="replace")
            i += (dt + 1) // 2; continue
        t += dt; out.append([t, code, ""])
    return out


def zs(x):
    return (x - x.mean()) / x.std()


# ------------------------------------------------------------------ preprocessing
def preprocess(d, rec):
    abp = read_signal(d, rec, "ABP")
    T = len(abp) / FS
    q = np.array([a[0] for a in read_ann(os.path.join(d, rec + ".wqrs")) if a[1] in (1, 5)]) / FS
    on = np.array([a[0] for a in read_ann(os.path.join(d, rec + ".wabp")) if a[1] == 1])
    # systolic BP of each pulse: max ABP between this onset and the next
    nxt = np.append(on[1:], len(abp))
    sbp = np.array([abp[a:b].max() if b > a else np.nan for a, b in zip(on, nxt)])
    ton = on / FS
    sbp_med = median_filter(np.nan_to_num(sbp, nan=np.nanmedian(sbp)), 11, mode="nearest")
    okS = np.isfinite(sbp) & (sbp > 50) & (sbp < 260) & (np.abs(sbp - sbp_med) < 0.25 * sbp_med)
    # RR_k = q[k] - q[k-1], stamped at q[k]
    rr = np.diff(q); tb = q[1:]
    med = median_filter(rr, 11, mode="nearest")
    okR = (rr > 0.3) & (rr < 2.5) & (np.abs(rr - med) < 0.3 * med)
    t_out = np.arange(0, T, 1.0 / FS_OUT)
    s4 = np.interp(t_out, ton[okS], sbp[okS]); r4 = np.interp(t_out, tb[okR], rr[okR])
    k = 60 * FS_OUT
    s4 = s4 - uniform_filter1d(s4, k, mode="nearest"); r4 = r4 - uniform_filter1d(r4, k, mode="nearest")
    # beat-domain series for the reference: SBP_n = pulse following R_n (onset within 0.6 s), RR_n = q[n]-q[n-1]
    m = np.searchsorted(ton, q, side="left")
    has = (m < len(ton)) & (np.where(m < len(ton), ton[np.minimum(m, len(ton) - 1)], 1e9) - q < 0.6)
    mm = np.minimum(m, len(ton) - 1)
    sb = np.where(has & okS[mm], sbp[mm], np.nan)
    rb = np.full(len(q), np.nan); rb[1:] = np.where(okR, rr, np.nan)
    ev = [(a / FS, c, s) for a, c, s in read_ann(os.path.join(d, rec + ".anI")) if c == 22]
    info = {"duration_s": T, "n_qrs": int(len(q)), "n_pulses": int(len(on)), "frac_sbp_rejected": float(1 - okS.mean()),
            "frac_rr_rejected": float(1 - okR.mean()), "median_rr_s": float(np.median(rr[okR]))}
    return zs(s4), zs(r4), (sb, rb), ev, info


def reference(beats, rr_med, jmax=5, det=60):
    """Beat-domain xcorr SBP_n vs RR_{n+j}, j=0..jmax (series detrended by a 60-beat moving average)."""
    sb, rb = beats
    def dt(v):
        w = np.isfinite(v); vi = np.interp(np.arange(len(v)), np.nonzero(w)[0], v[w])
        return vi - uniform_filter1d(vi, det, mode="nearest"), w
    s, ws = dt(sb); r, wr = dt(rb)
    cc = []
    for j in range(jmax + 1):
        n = len(s) - j; ok = ws[:n] & wr[j:j + n]
        cc.append(float(np.corrcoef(s[:n][ok], r[j:j + n][ok])[0, 1]))
    j = int(np.argmax(cc))
    return {"ref_lag_beats": j, "ref_delay_s": j * rr_med, "ref_r_max": cc[j], "ref_r_by_lag": cc}


def stand_episodes(ev):
    """(stand_start, back_to_supine) in s from the annotation texts."""
    eps, cur = [], None
    for t, _, s in ev:
        sl = s.lower()
        if "stand" in sl and cur is None:
            cur = t
        elif cur is not None and "supine" in sl:
            eps.append((cur, t)); cur = None
    return eps


# ------------------------------------------------------------------ scoring
def run(p, x, y, detail=False):
    tau, tv, _ = A.time_delay_interaction(x, y, p)
    lab = np.asarray(A.stable_label(tau, p)).astype(bool)
    dl = A.delay_of(tau, lab)
    out = (float(A.tds_score(lab)), None if dl is None else -dl / FS_OUT)
    return out + ((np.asarray(tv) / FS_OUT, lab),) if detail else out


def xcorr(x, y, lag):
    n = min(len(x), len(y)); x, y = zs(x[:n]), zs(y[:n])
    ks = np.arange(-lag, lag + 1)
    cc = np.array([np.mean(x[max(0, -k):n - max(0, k)] * y[max(0, k):n - max(0, -k)]) for k in ks])
    k = int(ks[int(np.argmax(np.abs(cc)))])
    return float(np.max(np.abs(cc))), k / FS_OUT      # y[t+k] ~ x[t], k>0: first signal leads


def _load(rec):
    x, y, beats, ev, info = preprocess(DIR, rec)
    return rec, x, y, beats, ev, info


def _subject_params(rec):
    x, y = SIG[rec]
    p = A.window_params([x, y])
    tc = max(_first_below_e(x), _first_below_e(y))
    return {"id": rec, "window": p.window, "max_lag": p.max_lag, "cap_applied": bool(p.max_lag < (p.window - 1) // 2),
            "tc": tc, "r_hat": max(1, int(round(tc / 2)))}


def _score(args):
    (i, j), settings, lag = args
    x, y = SIG[IDS[i]][0], SIG[IDS[j]][1]
    n = min(len(x), len(y)); x, y = x[:n], y[:n]
    o = {"x": IDS[i], "y": IDS[j], "true": i == j}
    eps = stand_episodes(EV[IDS[i]]) if i == j else []
    for name, p in settings.items():
        if i == j:
            o[f"{name}_score"], o[f"{name}_delay_s"], (tv, lab) = run(p, x, y, True)
            pre = np.zeros(len(tv), bool); st = np.zeros(len(tv), bool)
            for a, b in eps:
                pre |= (tv >= a - PRE_S) & (tv < a); st |= (tv >= a) & (tv < b)
            o[f"{name}_frac_before"] = float(100 * lab[pre].mean()) if pre.any() else None
            o[f"{name}_frac_standing"] = float(100 * lab[st].mean()) if st.any() else None
        else:
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
    err = [r[f"{m}_delay_s"] - ref[r["x"]]["ref_delay_s"] for r in c if r[f"{m}_delay_s"] is not None]
    S["n_with_delay"] = len(err)
    S["err_vs_ref_s_median_abs"] = float(np.median(np.abs(err))) if err else None
    S[f"n_within_{REF_TOL_S}s_of_ref"] = int(sum(abs(v) <= REF_TOL_S for v in err))
    S[f"n_within_{TOL2_S}s_of_ref"] = int(sum(abs(v) <= TOL2_S for v in err))
    fb = [r[f"{m}_frac_before"] for r in c if r.get(f"{m}_frac_before") is not None and r.get(f"{m}_frac_standing") is not None]
    fs = [r[f"{m}_frac_standing"] for r in c if r.get(f"{m}_frac_before") is not None and r.get(f"{m}_frac_standing") is not None]
    S["stable_pct_before_standing_mean"] = float(np.mean(fb)) if fb else None
    S["stable_pct_during_standing_mean"] = float(np.mean(fs)) if fs else None
    S["n_subjects_stable_lower_when_standing"] = int(sum(s < b_ for b_, s in zip(fb, fs)))
    S["n_subjects_with_standing"] = len(fb)
    return S


def main(d, out_dir):
    global DIR, SIG, EV
    DIR = d
    os.makedirs(out_dir, exist_ok=True)
    t0 = time.time()
    loaded = [_load(r) for r in IDS]
    SIG = {r: (x, y) for r, x, y, *_ in loaded}
    EV = {r: ev for r, _, _, _, ev, _ in loaded}
    info = {r: i for r, *_, i in loaded}
    ref = {r: reference(b, info[r]["median_rr_s"]) for r, _, _, b, _, _ in loaded}
    n = len(IDS)
    sp = [_subject_params(r) for r in IDS]
    r_hat = int(max(1, round(np.median([s["tc"] for s in sp]) / 2)))
    L = int(np.median([s["window"] for s in sp])); lag = int(np.median([s["max_lag"] for s in sp]))
    p0 = replace(A.PUBLISHED, window=L, overlap=L // 2, max_lag=min(lag, (L - 1) // 2), tolerance=0)
    atds, fail = A.calib(p0, [SIG[r] for r in IDS])
    settings = {"fixed": A.PUBLISHED,
                "rescaled": replace(A.PUBLISHED, window=60 * r_hat, overlap=30 * r_hat, max_lag=30 * r_hat, tolerance=r_hat),
                "atds": atds}
    k = n // 3
    assert np.gcd(k, n) == 1
    jobs = [(i, i) for i in range(n)] + [(i, (i + k) % n) for i in range(n)]
    rows = [_score((ij, settings, atds.max_lag)) for ij in jobs]
    methods = ["fixed", "rescaled", "atds", "xcorr"]
    S = {"n_subjects": n, "control_shift_k": k, "r_hat": r_hat, "tc_hat_median": float(np.median([s["tc"] for s in sp])),
         "params": {m: pd_(p) for m, p in settings.items()}, "atds_calibration_failed": bool(fail),
         "xcorr_lag_samples": atds.max_lag,
         "reference": {"median_ref_delay_s": float(np.median([v["ref_delay_s"] for v in ref.values()])),
                       "ref_lag_beats": {r: v["ref_lag_beats"] for r, v in ref.items()},
                       "median_rr_s": float(np.median([v["median_rr_s"] for v in info.values()]))},
         "atds_cap_applied_n_subjects": int(sum(s["cap_applied"] for s in sp)),
         "metrics": {m: metrics(rows, m, ref) for m in methods},
         "note": "n=10 subjects, 10 controls (k=3): AUC bootstrap CIs are wide"}
    persub = {s["id"]: s for s in sp}
    for r in rows:
        if r["true"]:
            r.update({k_: v for k_, v in ref[r["x"]].items() if k_ != "ref_r_by_lag"},
                     **{f"info_{k_}": v for k_, v in info[r["x"]].items()},
                     persub_window=persub[r["x"]]["window"], persub_max_lag=persub[r["x"]]["max_lag"],
                     persub_cap_applied=persub[r["x"]]["cap_applied"], persub_tc=persub[r["x"]]["tc"],
                     n_stand_episodes=len(stand_episodes(EV[r["x"]])))
    S["wall_seconds"] = time.time() - t0
    json.dump({"summary": S, "reference_by_lag": {r: v["ref_r_by_lag"] for r, v in ref.items()}, "records": rows},
              open(os.path.join(out_dir, "prcp.json"), "w"), indent=1)
    keys = list(dict.fromkeys(k_ for r in rows for k_ in r))
    with open(os.path.join(out_dir, "prcp.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)
    print(f"{'setting':9s} {'win/lag/tol (samples)':>22s} {'AUC [95% CI]':>18s} {'true mean':>9s} {'ctl95':>6s} {'>ctl95':>6s} {'delay':>5s} {'medAE s':>8s} {'<=1s':>5s} {'stab%pre':>8s} {'stab%stand':>10s}")
    f_ = lambda v: float("nan") if v is None else v
    for m in methods:
        q = S["metrics"][m]
        ps = "%d/%d/%d" % (settings[m].window, settings[m].max_lag, settings[m].tolerance) if m in settings else f"lag {atds.max_lag}"
        print(f"{m:9s} {ps:>22s} {q['auc']:.3f} [{q['auc_ci95'][0]:.2f},{q['auc_ci95'][1]:.2f}] {q['true_mean']:9.1f} "
              f"{q['control_95pct']:6.1f} {q['n_true_above_control95']:>3d}/{n} {q['n_with_delay']:>5d} "
              f"{f_(q['err_vs_ref_s_median_abs']):8.2f} {q[f'n_within_{REF_TOL_S}s_of_ref']:>5d} "
              f"{f_(q['stable_pct_before_standing_mean']):8.1f} {f_(q['stable_pct_during_standing_mean']):10.1f}")
    print(json.dumps({k_: S[k_] for k_ in ("reference", "r_hat", "tc_hat_median", "atds_cap_applied_n_subjects", "atds_calibration_failed")}))
    print("wall seconds", round(S["wall_seconds"], 1))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "results"))
