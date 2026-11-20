"""Cerebral vasoregulation (PhysioNet cves 1.0.0, transcranial-doppler/*-tcd-0, Day 2 supine protocol) with aTDS.

    python experiments/physio_new/cves.py <raw_dir> [out_dir]

<raw_dir> holds sNNNN-tcd-0.hea/.dat/.marker from https://physionet.org/files/cves/1.0.0/data/transcranial-doppler/
See SPEC_COMMON.md / SPEC_cves.md. TDS via this repo's pyTDS only (through experiments/atds_real/atds21.py).
Pairs: (co2_cbfv) x = end-tidal CO2, y = CBFV; (abp_cbfv) x = mean ABP, y = CBFV. 4 Hz. Delay = -tau0 (positive = x leads).
Segment: 300 s before hyperventilation start to the end of CO2 rebreathing (hyperventilation + rebreathing manoeuvres).
"""
import csv
import glob
import json
import os
import sys
import time
from dataclasses import replace

import numpy as np
from scipy.ndimage import uniform_filter1d
from scipy.signal import butter, find_peaks, sosfiltfilt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "atds_real"))
import atds21 as A  # noqa: E402
from bidmc import auc, boot_auc, _first_below_e  # noqa: E402  (helpers only; unchanged)

A.WORKERS = 1
FS_IN = 50
FS = 4
MAX_SUBJECTS = 40
PRE_S = 300.0
MAN_S = 180.0
TOL1_S, TOL2_S = 3.0, 5.0     # "within tolerance" of the reference delay (assumption: the spec states none)
XLAG_DIFF_S = 10.0            # lag range of the first-difference cross-correlation reference (ABP pair)
CH_ALIAS = {"cerebral_blood_velocity_right": "R", "cerebral_blood_velocity_left": "L", "abp": "abp", "co2": "co2"}


# ------------------------------------------------------------------ readers
def read_hea(path):
    L = [l.split() for l in open(path).read().splitlines() if l and not l.startswith("#")]
    nsig, n = int(L[0][1]), int(L[0][3])
    ch = []
    for s in L[1:1 + nsig]:
        g, rest = s[2].split("(")
        ch.append({"name": s[-1], "gain": float(g), "base": float(rest.split(")")[0])})
    return n, ch


def read_markers(path):
    """WFDB annotation stream (16-bit words: 6-bit type, 10-bit time increment; SKIP=59, AUX=63). Returns [(sample, label)]."""
    d = open(path, "rb").read(); i = 0; t = 0; ev = []
    while i + 1 < len(d):
        w = d[i] | (d[i + 1] << 8); i += 2; typ = w >> 10; n = w & 1023
        if typ == 59:
            t += ((d[i] | (d[i + 1] << 8)) << 16) | (d[i + 2] | (d[i + 3] << 8)); i += 4
        elif typ == 63:
            s = d[i:i + n].decode("latin1").rstrip("\0"); i += n + (n & 1)
            if ev:
                ev[-1] = (ev[-1][0], s)
        elif w == 0:
            break
        else:
            t += n; ev.append((t, None))
    return [(t, s) for t, s in ev if s]


def zs(x):
    return (x - x.mean()) / x.std()


def lp(x, fc=0.5):
    return sosfiltfilt(butter(2, fc, "low", fs=FS_IN, output="sos"), x)


def etco2(co2):
    """End-tidal CO2 (mmHg) at 50 Hz: capnogram end-expiratory peaks, linearly interpolated."""
    amp = np.percentile(co2, 99) - np.percentile(co2, 1)
    pk, _ = find_peaks(co2, distance=int(0.7 * FS_IN), prominence=max(3.0, 0.15 * amp))
    if len(pk) < 20:
        return None
    return np.interp(np.arange(len(co2)), pk, co2[pk])


def load(d, rec):
    n, ch = read_hea(os.path.join(d, rec + ".hea"))
    raw = np.fromfile(os.path.join(d, rec + ".dat"), "<i2", count=n * len(ch)).reshape(-1, len(ch)).astype(float)
    sig, bad = {}, {}
    for k, c in enumerate(ch):
        a = CH_ALIAS.get(c["name"])
        if a:
            sig[a] = (raw[:, k] - c["base"]) / c["gain"]
            bad[a] = float(np.mean(np.abs(raw[:, k]) >= 32766))
    return sig, bad


def man_times(mk):
    """Hyperventilation start (first HYPERVENT, else estimated 240 s before rebreathing) and rebreathing start (first *HYPOVENT*
    after it), plus the earliest VALSALVA marker after rebreathing (end bound)."""
    hv = [t / FS_IN for t, s in mk if s == "HYPERVENT"]
    hy = [t / FS_IN for t, s in mk if s.endswith("HYPOVENT")]
    if not hy:
        return None
    t_hy = hy[0]
    est = not hv
    t_hv = (hv[0] if hv and hv[0] < t_hy else t_hy - 240.0)
    est = est or not (hv and hv[0] < t_hy)
    va = [t / FS_IN for t, s in mk if s == "VALSALVA" and t / FS_IN > t_hy + 60]
    t_end = t_hy + MAN_S + 20
    if va:
        t_end = min(t_end, va[0])
    return t_hv, t_hy, t_end, est


def preprocess(d, rec):
    names = {c["name"] for c in read_hea(os.path.join(d, rec + ".hea"))[1]}
    if "co2" not in names or "abp" not in names:
        return None, "no co2/abp channel in header"
    sig, bad = load(d, rec)
    if not all(k in sig for k in ("abp", "co2")) or not ({"R", "L"} & set(sig)):
        return None, "missing channel"
    mt = man_times(read_markers(os.path.join(d, rec + ".marker")))
    if mt is None:
        return None, "no rebreathing marker"
    t_hv, t_hy, t_end, est = mt
    good = [k for k in ("R", "L") if k in sig and bad[k] < 0.01 and np.std(sig[k]) > 1]
    if not good or bad["abp"] > 0.01 or bad["co2"] > 0.01:
        return None, "saturated/flat channel"
    et = etco2(sig["co2"])
    if et is None:
        return None, "no capnogram peaks"
    n = len(sig["abp"]); t_out = np.arange(0, n / FS_IN, 1.0 / FS); t = np.arange(n) / FS_IN
    cb = np.mean([lp(sig[k]) for k in good], axis=0)
    raw4 = {"co2": np.interp(t_out, t, lp(et)), "abp": np.interp(t_out, t, lp(sig["abp"])), "cbfv": np.interp(t_out, t, cb)}
    ma = {k: v - uniform_filter1d(v, 120 * FS, mode="nearest") for k, v in raw4.items()}     # 120-s moving average removed
    i0, i1 = int(max(0, t_hv - PRE_S) * FS), int(min(t_end, n / FS_IN) * FS)
    seg = {k: v[i0:i1] for k, v in ma.items()}
    if len(seg["co2"]) < 400 * FS or not all(np.std(v) > 0 for v in seg.values()):
        return None, "segment too short/flat"
    info = {"t_hv_s": t_hv, "t_hy_s": t_hy, "hv_start_estimated": bool(est), "seg_start_s": i0 / FS, "seg_len_s": (i1 - i0) / FS,
            "cbfv_sides": "".join(good), "n_man_end_s": t_end}
    return {"x_co2": zs(seg["co2"]), "x_abp": zs(seg["abp"]), "y": zs(seg["cbfv"]), "raw": raw4, "seg0": i0}, info


# ------------------------------------------------------------------ reference delays
def _cross(v, target, rising, t0):
    """First time >= t0 (s) where v crosses target in the stated direction."""
    for k in range(int(t0 * FS), len(v)):
        if (v[k] >= target) if rising else (v[k] <= target):
            return k / FS
    return None


def step_ref(raw, t_hv, t_hy, hv_est):
    """Step-response delays (s, positive = CBFV lags CO2): per manoeuvre (hyperventilation, rebreathing) smoothed (5 s) ETCO2 and CBFV.
    onset63 = CBFV 63 % time - CO2 change onset (first 10 % of its change); t63 = CBFV 63 % time - CO2 63 % time."""
    co2, cb = (uniform_filter1d(raw[k], 5 * FS, mode="nearest") for k in ("co2", "cbfv"))
    out = []
    for name, t0 in (("hv", t_hv), ("rebreath", t_hy)):
        if name == "hv" and hv_est:
            continue
        b = slice(int((t0 - 30) * FS), int(t0 * FS)); e = slice(int((t0 + MAN_S - 30) * FS), int((t0 + MAN_S) * FS))
        if e.stop > len(co2) or b.start < 0:
            continue
        dc, db = co2[e].mean() - co2[b].mean(), cb[e].mean() - cb[b].mean()
        if abs(dc) < 3 or abs(db) < 3 or np.sign(dc) != np.sign(db):    # need a real, same-direction CO2 and CBFV change
            continue
        c10 = _cross(co2, co2[b].mean() + 0.1 * dc, dc > 0, t0)
        c63 = _cross(co2, co2[b].mean() + 0.63 * dc, dc > 0, t0)
        f63 = _cross(cb, cb[b].mean() + 0.63 * db, db > 0, t0)
        if None in (c10, c63, f63):
            continue
        out.append({"man": name, "onset63": f63 - c10, "t63": f63 - c63, "dco2": float(dc), "dcbfv": float(db)})
    return out


def diff_xcorr_lag(x, y, lag):
    dx, dy = np.diff(x), np.diff(y); n = len(dx); dx, dy = zs(dx), zs(dy)
    ks = np.arange(-lag, lag + 1)
    cc = np.array([np.mean(dx[max(0, -k):n - max(0, k)] * dy[max(0, k):n - max(0, -k)]) for k in ks])
    return float(ks[int(np.argmax(np.abs(cc)))]) / FS


# ------------------------------------------------------------------ scoring
def run(p, x, y):
    tau = A.time_delay_interaction(x, y, p)[0]
    lab = A.stable_label(tau, p)
    dl = A.delay_of(tau, lab)
    return float(A.tds_score(lab)), None if dl is None else -dl / FS


def xcorr(x, y, lag):
    n = min(len(x), len(y)); x, y = zs(x[:n]), zs(y[:n])
    ks = np.arange(-lag, lag + 1)
    cc = np.array([np.mean(x[max(0, -k):n - max(0, k)] * y[max(0, k):n - max(0, -k)]) for k in ks])
    k = int(ks[int(np.argmax(np.abs(cc)))])
    return float(np.max(np.abs(cc))), k / FS         # y[t+k] ~ x[t], k>0: first signal leads


def pd_(p):
    return {"window": p.window, "step": p.overlap, "max_lag": p.max_lag, "tolerance": p.tolerance,
            "window_s": p.window / FS, "max_lag_s": p.max_lag / FS, "tolerance_s": p.tolerance / FS}


def metrics(rows, m, refs):
    c = [r for r in rows if r["true"]]; u = [r for r in rows if not r["true"]]
    a = [r[f"{m}_score"] for r in c]; b = [r[f"{m}_score"] for r in u]; thr = float(np.percentile(b, 95))
    S = {"true_mean": float(np.mean(a)), "true_median": float(np.median(a)), "control_mean": float(np.mean(b)),
         "control_median": float(np.median(b)), "control_95pct": thr, "auc": auc(a, b), "auc_ci95": boot_auc(a, b),
         "n_true_above_control95": int(sum(v > thr for v in a)), "n_true": len(a)}
    S["n_with_delay"] = int(sum(r[f"{m}_delay_s"] is not None for r in c))
    for rn in refs:
        e = [r[f"{m}_delay_s"] - r[f"ref_{rn}_s"] for r in c if r[f"{m}_delay_s"] is not None and r[f"ref_{rn}_s"] is not None]
        S[f"n_with_delay_and_ref_{rn}"] = len(e)
        S[f"err_vs_{rn}_s_median_abs"] = float(np.median(np.abs(e))) if e else None
        S[f"err_vs_{rn}_s_median_signed"] = float(np.median(e)) if e else None
        S[f"n_within_{TOL1_S:g}s_of_{rn}"] = int(sum(abs(v) <= TOL1_S for v in e))
        S[f"n_within_{TOL2_S:g}s_of_{rn}"] = int(sum(abs(v) <= TOL2_S for v in e))
        d_ = [r[f"{m}_delay_s"] for r in c if r[f"{m}_delay_s"] is not None]
        S["delay_found_median_s"] = float(np.median(d_)) if d_ else None
    return S


def run_pair(pair, ids, sig, refs, extra_rows):
    n = len(ids)
    sp = []
    for r in ids:
        x, y = sig[r]
        p = A.window_params([x, y])
        tc = max(_first_below_e(x), _first_below_e(y))
        sp.append({"id": r, "window": p.window, "max_lag": p.max_lag, "cap_applied": bool(p.max_lag < (p.window - 1) // 2),
                   "tc": tc, "r_hat": max(1, int(round(tc / 2)))})
    r_hat = int(np.median([s["r_hat"] for s in sp]))
    L = int(np.median([s["window"] for s in sp])); lag = int(np.median([s["max_lag"] for s in sp]))
    p0 = replace(A.PUBLISHED, window=L, overlap=L // 2, max_lag=min(lag, (L - 1) // 2), tolerance=0)
    atds, fail = A.calib(p0, [sig[r] for r in ids])
    settings = {"fixed": A.PUBLISHED,
                "rescaled": replace(A.PUBLISHED, window=60 * r_hat, overlap=30 * r_hat, max_lag=30 * r_hat, tolerance=r_hat),
                "atds": atds}
    k = n // 3
    while np.gcd(k, n) != 1:
        k += 1
    jobs = [(i, i) for i in range(n)] + [(i, (i + k) % n) for i in range(n)]
    rows = []
    for i, j in jobs:
        x, y = sig[ids[i]][0], sig[ids[j]][1]
        m_ = min(len(x), len(y)); x, y = x[:m_], y[:m_]
        o = {"pair": pair, "x": ids[i], "y": ids[j], "true": i == j}
        for name, p in settings.items():
            o[f"{name}_score"], o[f"{name}_delay_s"] = run(p, x, y)
        o["xcorr_score"], o["xcorr_delay_s"] = xcorr(x, y, atds.max_lag)
        if i == j:
            for rn in refs:
                o[f"ref_{rn}_s"] = refs[rn][ids[i]]
            o.update(extra_rows[ids[i]])
            ps = next(s for s in sp if s["id"] == ids[i])
            o.update(persub_window=ps["window"], persub_max_lag=ps["max_lag"], persub_cap_applied=ps["cap_applied"], persub_tc=ps["tc"])
        rows.append(o)
    methods = ["fixed", "rescaled", "atds", "xcorr"]
    S = {"n_subjects": n, "control_shift_k": k, "r_hat_median": r_hat, "r_hat_range": [min(s["r_hat"] for s in sp), max(s["r_hat"] for s in sp)],
         "params": {m: pd_(p) for m, p in settings.items()}, "atds_calibration_failed": bool(fail),
         "xcorr_lag_samples": atds.max_lag, "atds_cap_applied_n_subjects": int(sum(s["cap_applied"] for s in sp)),
         "reference_median_s": {rn: (float(np.median([v for v in refs[rn].values() if v is not None]))
                                     if any(v is not None for v in refs[rn].values()) else None) for rn in refs},
         "reference_n_available": {rn: int(sum(v is not None for v in refs[rn].values())) for rn in refs},
         "metrics": {m: metrics(rows, m, list(refs)) for m in methods}}
    return S, rows, settings, atds


def main(d, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    t0 = time.time()
    recs = sorted(os.path.basename(f)[:-4] for f in glob.glob(os.path.join(d, "*-tcd-0.hea")))
    ids, data, info, excluded = [], {}, {}, {}
    for r in recs:
        if len(ids) >= MAX_SUBJECTS:
            break
        out, inf = preprocess(d, r)
        if out is None:
            excluded[r] = inf; continue
        ids.append(r); data[r] = out; info[r] = inf
    print("subjects used", len(ids), "excluded", excluded)
    refs_co2 = {"step_t63": {}, "step_onset": {}}
    extra = {}
    for r in ids:
        sr = step_ref(data[r]["raw"], info[r]["t_hv_s"], info[r]["t_hy_s"], info[r]["hv_start_estimated"])
        refs_co2["step_t63"][r] = float(np.median([s["t63"] for s in sr])) if sr else None
        refs_co2["step_onset"][r] = float(np.median([s["onset63"] for s in sr])) if sr else None
        extra[r] = {**info[r], "ref_n_manoeuvres": len(sr), "ref_manoeuvres": ",".join(s["man"] for s in sr)}
    refs_abp = {"diff_xcorr": {r: diff_xcorr_lag(data[r]["x_abp"], data[r]["y"], int(XLAG_DIFF_S * FS)) for r in ids}}
    out = {"summary": {"n_subjects_used": len(ids), "excluded": excluded, "subjects": ids, "pairs": {}}, "records": []}
    for pair, key, refs in (("co2_cbfv", "x_co2", refs_co2), ("abp_cbfv", "x_abp", refs_abp)):
        sig = {r: (data[r][key], data[r]["y"]) for r in ids}
        S, rows, settings, atds = run_pair(pair, ids, sig, refs, extra)
        out["summary"]["pairs"][pair] = S; out["records"] += rows
        print(f"\n=== {pair}  n={len(ids)}  r_hat={S['r_hat_median']}  calib_failed={S['atds_calibration_failed']}  ref medians {S['reference_median_s']} (n {S['reference_n_available']})")
        for rn in refs:
            print(f" reference {rn}")
            print(f"{'setting':9s} {'win/lag/tol':>12s} {'AUC [95% CI]':>18s} {'trueMean':>8s} {'ctl95':>6s} {'>ctl95':>6s} {'delay':>5s} {'medAE s':>8s} {f'<={TOL1_S:g}s':>5s} {f'<={TOL2_S:g}s':>5s} {'med delay':>9s}")
            for m in ("fixed", "rescaled", "atds", "xcorr"):
                q = S["metrics"][m]
                ps = "%d/%d/%d" % (settings[m].window, settings[m].max_lag, settings[m].tolerance) if m in settings else f"lag {atds.max_lag}"
                ae = q[f"err_vs_{rn}_s_median_abs"]
                md = q["delay_found_median_s"]
                print(f"{m:9s} {ps:>12s} {q['auc']:.3f} [{q['auc_ci95'][0]:.2f},{q['auc_ci95'][1]:.2f}] {q['true_mean']:8.1f} {q['control_95pct']:6.1f} "
                      f"{q['n_true_above_control95']:>3d}/{len(ids)} {q['n_with_delay']:>5d} {ae if ae is not None else float('nan'):8.2f} "
                      f"{q[f'n_within_{TOL1_S:g}s_of_{rn}']:>5d} {q[f'n_within_{TOL2_S:g}s_of_{rn}']:>5d} {md if md is not None else float('nan'):9.2f}")
    out["summary"]["wall_seconds"] = time.time() - t0
    json.dump(out, open(os.path.join(out_dir, "cves.json"), "w"), indent=1)
    keys = list(dict.fromkeys(k_ for r in out["records"] for k_ in r))
    with open(os.path.join(out_dir, "cves.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(out["records"])
    print("wall seconds", round(out["summary"]["wall_seconds"], 1))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "results"))
