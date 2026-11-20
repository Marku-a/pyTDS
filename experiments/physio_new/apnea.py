"""Apnea: breathing envelope -> SpO2 coupling in PhysioNet ucddb 1.0.0 (see SPEC_COMMON.md / SPEC_apnea.md).

    python experiments/physio_new/apnea.py <raw_dir> [out_dir]

<raw_dir> holds ucddbNNN.rec and ucddbNNN_respevt.txt from https://physionet.org/files/ucddb/1.0.0/
TDS via this repo's pyTDS only (through experiments/atds_real/atds21.py). Delays: seconds, positive = envelope leads SpO2.
"""
import csv
import json
import os
import re
import sys
import time
from dataclasses import replace
from multiprocessing import Pool

import numpy as np
from scipy.ndimage import uniform_filter1d
from scipy.signal import butter, filtfilt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "atds_real"))
import atds21 as A  # noqa: E402

FS = 8
WORKERS = 2
A.WORKERS = WORKERS
PAD = 30            # event period = event +- 30 s
REF_LAG = 60        # reference (b) lag range, s
TOL_S = 10.0        # "within tolerance" of the delay reference (not stated in the spec; chosen here)
SIG = {}            # subject -> dict(x, y, mask, events)  (inherited by forked workers)


# ---------------------------------------------------------------- data
def read_edf(path):
    with open(path, "rb") as f:
        h = f.read(256)
        ns = int(h[252:256]); nrec = int(h[236:244]); dur = float(h[244:252])
        hh = f.read(256 * ns)
        start = h[176:184].decode().strip()
        def fld(off, w): return [hh[off * ns + i * w: off * ns + (i + 1) * w].decode().strip() for i in range(ns)]
        label = fld(0, 16); pmin = fld(16 + 80 + 8, 8); pmax = fld(16 + 80 + 16, 8)
        dmin = fld(16 + 80 + 24, 8); dmax = fld(16 + 80 + 32, 8); nsamp = [int(v) for v in fld(16 + 80 + 40 + 80, 8)]
        raw = np.fromfile(f, dtype="<i2")
    per = sum(nsamp); nrec = min(nrec, len(raw) // per); raw = raw[:nrec * per].reshape(nrec, per)
    out = {}
    for name in ("Flow", "SpO2"):
        i = label.index(name); o = sum(nsamp[:i])
        assert nsamp[i] == int(FS * dur), (name, nsamp[i])
        d = raw[:, o:o + nsamp[i]].reshape(-1).astype(float)
        g = (float(pmax[i]) - float(pmin[i])) / (float(dmax[i]) - float(dmin[i]))
        out[name] = (d - float(dmin[i])) * g + float(pmin[i])
    hms = [int(v) for v in start.split(".")]
    return out["Flow"], out["SpO2"], hms[0] * 3600 + hms[1] * 60 + hms[2]


def read_events(path, t0):
    ev = []
    for line in open(path, errors="ignore"):
        m = re.match(r"^\s*(\d\d):(\d\d):(\d\d)\s+((?:APNEA|HYP)-\S+)\s+(?:[A-Za-z]\S*\s+)?(\d+(?:\.\d+)?)", line)
        if not m:
            continue
        t = int(m[1]) * 3600 + int(m[2]) * 60 + int(m[3])
        ev.append(((t - t0) % 86400, float(m[5]), m[4]))
    return ev


def movmean(x, n):
    return uniform_filter1d(x, size=n, mode="nearest")


def preprocess(flow, spo2):
    b, a = butter(2, [0.05 / (FS / 2), 2.0 / (FS / 2)], "band")
    f = filtfilt(b, a, A.interp_nan(flow))
    env = np.sqrt(np.maximum(movmean(f ** 2, 4 * FS), 0.0))
    s = np.asarray(spo2, float).copy(); s[(s < 50) | (s > 100)] = np.nan; s = A.interp_nan(s)
    k = 300 * FS
    return env - movmean(env, k), s - movmean(s, k)


def load(sid):
    flow, spo2, t0 = read_edf(os.path.join(DIR, f"{sid}.rec"))
    n = min(len(flow), len(spo2))
    x, y = preprocess(flow[:n], spo2[:n])
    ev = [e for e in read_events(os.path.join(DIR, f"{sid}_respevt.txt"), t0) if e[0] * FS < n]
    mask = np.zeros(n, bool)
    for o, d, _ in ev:
        mask[max(0, int((o - PAD) * FS)):min(n, int((o + d + PAD) * FS))] = True
    return dict(x=x, y=y, mask=mask, events=ev, t0=t0)


# ---------------------------------------------------------------- references
def ref_a(x, y, events):
    """Per event: SpO2 minimum in [onset, end+60 s] minus envelope minimum in the event, in s; median over events."""
    d = []
    for o, dur, _ in events:
        a = int(o * FS); e = int((o + dur) * FS); n = len(x)
        if e - a < 2 or a >= n: continue
        ty = a + int(np.argmin(y[a:min(n, e + 60 * FS)])); tx = a + int(np.argmin(x[a:min(n, e)]))
        d.append((ty - tx) / FS)
    return (float(np.median(d)) if d else None), len(d)


def _xc(x, y, lags):
    """sum_t x[t] y[t+k] for k in lags via FFT (linear)."""
    n = len(x); nf = 1 << int(np.ceil(np.log2(2 * n)))
    c = np.fft.irfft(np.conj(np.fft.rfft(x, nf)) * np.fft.rfft(y, nf), nf)
    return np.array([c[k] if k >= 0 else c[nf + k] for k in lags])


def ref_b(x, y, mask, lag_s=REF_LAG):
    """Event-period cross-correlation: both signals z-scored on event periods and zeroed elsewhere; lag of max |C|, s."""
    if mask.sum() < 10 * FS: return None, None
    lags = np.arange(-lag_s * FS, lag_s * FS + 1)
    out = []
    for both in (True, False):
        xm = np.where(mask, (x - x[mask].mean()) / x[mask].std(), 0.0)
        ym = (y - y[mask].mean()) / y[mask].std(); ym = np.where(mask, ym, 0.0) if both else ym
        C = _xc(xm, ym, lags) / np.sqrt(np.sum(xm ** 2) * np.sum(ym ** 2))
        out.append(float(lags[int(np.argmax(np.abs(C)))] / FS))
    return out[0], out[1]


def xcorr_full(x, y, lag):
    """Full-record normalised xcorr of the z-scored pair; delay (s, + = x leads) = argmax |C| within +-lag; score = max |C|."""
    xz = (x - x.mean()) / x.std(); yz = (y - y.mean()) / y.std(); lags = np.arange(-lag, lag + 1)
    C = _xc(xz, yz, lags) / len(x)
    k = int(np.argmax(np.abs(C)))
    return float(abs(C[k])), float(lags[k] / FS)


# ---------------------------------------------------------------- scoring
def pair_xy(ij):
    i, j = ij
    n = min(len(SIG[i]["x"]), len(SIG[j]["y"]))
    return SIG[i]["x"][:n], SIG[j]["y"][:n], SIG[i]["mask"][:n]


def score_pair(args):
    ij, p = args
    x, y, mask = pair_xy(ij)
    tau, t, _ = A.time_delay_interaction(x, y, p)
    lab = np.asarray(A.stable_label(tau, p)).astype(bool)
    inm = mask[np.clip(np.round(t).astype(int), 0, len(mask) - 1)]
    dl = A.delay_of(tau, lab)
    din = A.delay_of(tau[inm], lab[inm]) if inm.any() else None
    return {"score": float(A.tds_score(lab)), "delay_s": None if dl is None else -dl / FS,
            "delay_in_events_s": None if din is None else -din / FS, "n_windows": int(len(tau)), "n_stable": int(lab.sum()),
            "pct_stable_in_events": float(100 * lab[inm].mean()) if inm.any() else None,
            "pct_stable_outside": float(100 * lab[~inm].mean()) if (~inm).any() else None}


def xcorr_pair(args):
    ij, lag = args
    x, y, _ = pair_xy(ij)
    s, d = xcorr_full(x, y, lag)
    return {"score": s, "delay_s": d}


def sid_window(sid):
    p = A.window_params([SIG[sid]["x"], SIG[sid]["y"]])
    return p.window, p.max_lag


def sid_rhat(sid):
    def fb(v):
        r = A.acf(v, len(v) // 4); k = np.nonzero(r < np.exp(-1.0))[0]
        return float(k[0]) if len(k) else float(len(r) - 1)
    tc = max(fb(SIG[sid]["x"]), fb(SIG[sid]["y"]))
    return max(1, int(round(tc / 2)))


def load_job(sid):
    return sid, load(sid)


def auc(a, b):
    a, b = np.asarray(a)[:, None], np.asarray(b)[None, :]
    return float((np.sum(a > b) + 0.5 * np.sum(a == b)) / (a.size * b.size))


def boot_auc(a, b, n=5000, seed=0):
    rng = np.random.default_rng(seed); a, b = np.asarray(a), np.asarray(b)
    v = [auc(rng.choice(a, a.size), rng.choice(b, b.size)) for _ in range(n)]
    return [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]


def metrics(c, u, refs):
    a = [r["score"] for r in c]; b = [r["score"] for r in u]; thr = float(np.percentile(b, 95))
    S = {"coupled_mean": float(np.mean(a)), "coupled_median": float(np.median(a)), "control_mean": float(np.mean(b)),
         "control_median": float(np.median(b)), "control_95pct": thr, "auc": auc(a, b), "auc_ci95": boot_auc(a, b),
         "n_above_control95": int(sum(v > thr for v in a)), "n_coupled": len(a), "n_with_delay": sum(r["delay_s"] is not None for r in c)}
    for name, ref in refs.items():
        e = [abs(r["delay_s"] - ref[r["sid"]]) for r in c if r["delay_s"] is not None and ref[r["sid"]] is not None]
        S[f"abs_err_vs_{name}_median_s"] = float(np.median(e)) if e else None
        S[f"n_within_{TOL_S:g}s_of_{name}"] = int(sum(v <= TOL_S for v in e))
        S[f"median_delay_minus_{name}_s"] = float(np.median([r["delay_s"] - ref[r["sid"]] for r in c if r["delay_s"] is not None
                                                             and ref[r["sid"]] is not None])) if e else None
    for tag, rows in (("coupled", c), ("control", u)):
        i = [r["pct_stable_in_events"] for r in rows if r.get("pct_stable_in_events") is not None]
        o = [r["pct_stable_outside"] for r in rows if r.get("pct_stable_outside") is not None]
        if i: S[f"{tag}_pct_stable_in_events_mean"] = float(np.mean(i)); S[f"{tag}_pct_stable_outside_mean"] = float(np.mean(o))
    return S


def main(d, out_dir):
    global DIR
    DIR = d
    os.makedirs(out_dir, exist_ok=True)
    t_start = time.time()
    ids = [l.strip() for l in open(os.path.join(d, "IDS"))] if os.path.exists(os.path.join(d, "IDS")) else \
        sorted(f[:-4] for f in os.listdir(d) if re.fullmatch(r"ucddb\d+\.rec", f))
    with Pool(WORKERS) as pool:
        for sid, s in pool.imap(load_job, ids):
            SIG[sid] = s
    n = len(ids); k = max(1, round(n / 3))
    print(f"loaded {n} subjects; control shift {k}", flush=True)

    refa, refb, refbx, nev = {}, {}, {}, {}
    for sid in ids:
        s = SIG[sid]
        refa[sid], nev[sid] = ref_a(s["x"], s["y"], s["events"])
        refb[sid], refbx[sid] = ref_b(s["x"], s["y"], s["mask"])
    refs = {"ref_a": refa, "ref_b": refb}

    # parameters
    with Pool(WORKERS) as pool:
        wl = pool.map(sid_window, ids)
        rh = pool.map(sid_rhat, ids)
    L = int(np.median([w for w, _ in wl])); lag = min(int(np.median([g for _, g in wl])), (L - 1) // 2)
    p0 = replace(A.PUBLISHED, window=L, overlap=L // 2, max_lag=lag, tolerance=0)
    pa, fail = A.calib(p0, [pair_xy((s, s))[:2] for s in ids])
    r = max(1, int(np.median(rh)))
    prs = replace(A.PUBLISHED, window=60 * r, overlap=30 * r, max_lag=30 * r, tolerance=r)
    params = {"fixed": A.PUBLISHED, "rescaled": prs, "atds": pa}
    pinfo = {k_: {"window": p.window, "step": p.overlap, "max_lag": p.max_lag, "tolerance": p.tolerance,
                  "window_s": p.window / FS, "step_s": p.overlap / FS, "max_lag_s": p.max_lag / FS, "tolerance_s": p.tolerance / FS}
             for k_, p in params.items()}
    pinfo["rescaled"]["r_hat"] = r; pinfo["rescaled"]["r_hat_range"] = [int(min(rh)), int(max(rh))]
    pinfo["atds"]["calibration_failed"] = bool(fail)
    pinfo["atds"]["per_subject_window_s_range"] = [min(w for w, _ in wl) / FS, max(w for w, _ in wl) / FS]
    pinfo["atds"]["per_subject_lag_s_median"] = float(np.median([g for _, g in wl]) / FS)
    pinfo["xcorr"] = {"max_lag": pa.max_lag, "max_lag_s": pa.max_lag / FS}
    print(json.dumps(pinfo), flush=True)

    jobs = [(s, s) for s in ids] + [(ids[i], ids[(i + k) % n]) for i in range(n)]
    rows = [{"x_subject": a, "y_subject": b, "coupled": a == b} for a, b in jobs]
    for name, p in params.items():
        with Pool(WORKERS) as pool:
            res = pool.map(score_pair, [(ij, p) for ij in jobs])
        for row, rs in zip(rows, res):
            row.update({f"{name}_{kk}": v for kk, v in rs.items()})
        print(name, "done", round(time.time() - t_start), "s", flush=True)
    with Pool(WORKERS) as pool:
        res = pool.map(xcorr_pair, [(ij, pa.max_lag) for ij in jobs])
    for row, rs in zip(rows, res):
        row.update({f"xcorr_{kk}": v for kk, v in rs.items()})
    for row in rows:
        s = row["x_subject"]; row.update(ref_a_s=refa[s], ref_b_s=refb[s], ref_b_xmask_s=refbx[s], n_events=nev[s])

    summary = {"n_subjects": n, "control_shift": k, "params": pinfo, "tolerance_for_within_s": TOL_S,
               "ref_a_median_s": float(np.median([v for v in refa.values() if v is not None])),
               "ref_b_median_s": float(np.median([v for v in refb.values() if v is not None])),
               "corr_ref_a_ref_b": None}
    va = [(refa[s], refb[s]) for s in ids if refa[s] is not None and refb[s] is not None]
    if len(va) > 2: summary["corr_ref_a_ref_b"] = float(np.corrcoef(*zip(*va))[0, 1])
    summary["ref_a_minus_ref_b_abs_median_s"] = float(np.median([abs(a - b) for a, b in va])) if va else None
    for name in ("fixed", "rescaled", "atds", "xcorr"):
        c = [{"sid": r_["x_subject"], **{kk[len(name) + 1:]: v for kk, v in r_.items() if kk.startswith(name + "_")}}
             for r_ in rows if r_["coupled"]]
        u = [{"sid": r_["x_subject"], **{kk[len(name) + 1:]: v for kk, v in r_.items() if kk.startswith(name + "_")}}
             for r_ in rows if not r_["coupled"]]
        summary[name] = metrics(c, u, refs)
    # probe: subject 002 (all settings use the common sets; the per-subject aTDS set is reported separately)
    probe = {}
    if "ucddb002" in SIG:
        s = "ucddb002"; row = next(r_ for r_ in rows if r_["coupled"] and r_["x_subject"] == s)
        pp = A.window_params([SIG[s]["x"], SIG[s]["y"]])
        pp, pf = A.calib(pp, [pair_xy((s, s))[:2]])
        own = score_pair(((s, s), pp))
        x, y, mask = pair_xy((s, s))
        tau, t, _ = A.time_delay_interaction(x, y, pp); lab = np.asarray(A.stable_label(tau, pp)).astype(bool)
        from collections import Counter
        top = Counter((-tau[lab] / FS).round(1).tolist()).most_common(8)
        probe = {"own_params": {"window_s": pp.window / FS, "max_lag_s": pp.max_lag / FS, "tolerance_s": pp.tolerance / FS,
                                "calibration_failed": bool(pf)}, "own": own, "own_top_stable_delays_s": top,
                 "own_median_stable_delay_s": float(np.median(-tau[lab] / FS)),
                 "own_stable_delay_iqr_s": [float(np.percentile(-tau[lab] / FS, q)) for q in (25, 75)],
                 "ref_a_s": refa[s], "ref_b_s": refb[s], "ref_b_xmask_s": refbx[s], "n_events": nev[s],
                 "common": {nm: {kk: row[f"{nm}_{kk}"] for kk in ("delay_s", "delay_in_events_s", "pct_stable_in_events", "pct_stable_outside")}
                            for nm in ("fixed", "rescaled", "atds")}, "xcorr_delay_s": row["xcorr_delay_s"]}
        # reference (b) restricted to |lag|<=lag of the aTDS setting and event-by-event spread of (a)
        evs = [e for e in SIG[s]["events"]]
        ds = []
        for o, dur, _ in evs:
            a = int(o * FS); e = int((o + dur) * FS)
            if e - a >= 2:
                ds.append((a + int(np.argmin(y[a:min(len(y), e + 60 * FS)])) - a - int(np.argmin(x[a:e]))) / FS)
        probe["ref_a_event_delays_percentiles_s"] = [float(np.percentile(ds, q)) for q in (10, 25, 50, 75, 90)] if ds else None
    summary["probe_ucddb002"] = probe
    summary["wall_seconds"] = time.time() - t_start

    json.dump({"summary": summary, "records": rows}, open(os.path.join(out_dir, "apnea.json"), "w"), indent=1)
    keys = list(rows[0]);
    with open(os.path.join(out_dir, "apnea.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)

    print("\nsetting   | coupled mean/med | ctrl mean/med/p95 | AUC [CI] | >p95 | n delay | medAbsErr a/b (s) | within %gs a/b | stable%% in/out (coupled) | ctrl in/out" % TOL_S)
    for name in ("fixed", "rescaled", "atds", "xcorr"):
        m = summary[name]
        print(f"{name:9s} | {m['coupled_mean']:.3f}/{m['coupled_median']:.3f} | {m['control_mean']:.3f}/{m['control_median']:.3f}/{m['control_95pct']:.3f} | "
              f"{m['auc']:.2f} [{m['auc_ci95'][0]:.2f},{m['auc_ci95'][1]:.2f}] | {m['n_above_control95']}/{m['n_coupled']} | {m['n_with_delay']} | "
              f"{m['abs_err_vs_ref_a_median_s']} / {m['abs_err_vs_ref_b_median_s']} | {m[f'n_within_{TOL_S:g}s_of_ref_a']}/{m[f'n_within_{TOL_S:g}s_of_ref_b']} | "
              f"{m.get('coupled_pct_stable_in_events_mean', float('nan')):.1f}/{m.get('coupled_pct_stable_outside_mean', float('nan')):.1f} | "
              f"{m.get('control_pct_stable_in_events_mean', float('nan')):.1f}/{m.get('control_pct_stable_outside_mean', float('nan')):.1f}")
    print("refs: median a", summary["ref_a_median_s"], "median b", summary["ref_b_median_s"], "corr(a,b)", summary["corr_ref_a_ref_b"])
    print("probe 002:", json.dumps(probe))
    print("wall seconds", round(summary["wall_seconds"]))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "results"))
