"""Treadmill exercise (PhysioNet treadmill-exercise-cardioresp 1.0.1): VO2->VCO2 and speed->HR, fixed vs rescaled vs aTDS vs xcorr.

    python experiments/physio_new/treadmill.py <raw_dir>

See SPEC_COMMON.md / SPEC_treadmill.md. TDS via pyTDS only (through experiments/atds_real/atds21.py).
Preprocessing: per test, mean over duplicated time stamps, linear interpolation to a 1-Hz grid, minus a centred 60-s moving
average (variant: no detrending, aTDS and fixed only). Delay convention: positive = first signal leads (= -tau0), seconds = samples.
"""
import csv
import json
import os
import sys
import time
from dataclasses import replace
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.ndimage import uniform_filter1d

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "atds_real"))
import atds21 as A  # noqa: E402

A.WORKERS = 2
NSAMPLE, SEED, REFLAG, TOL_S = 200, 0, 60, 5.0
PAIRS = {"VO2_VCO2": ("VO2", "VCO2"), "Speed_HR": ("Speed", "HR")}
DATA: dict = {}


# ------------------------------------------------------------------ data
def load(d):
    df = pd.read_csv(os.path.join(d, "test_measure.csv"))
    sel = []
    for tid, g in df.groupby("ID_test", sort=True):
        if g[["HR", "VO2", "VCO2"]].isna().any().any() or g.time.max() < 480:
            continue
        sel.append(tid)
    n_valid = len(sel)
    rng = np.random.default_rng(SEED)
    sel = [sel[k] for k in sorted(rng.choice(len(sel), NSAMPLE, replace=False))]
    sel = [sel[k] for k in rng.permutation(len(sel))]
    out = {}
    for tid, g in df[df.ID_test.isin(sel)].groupby("ID_test"):
        g = g.groupby("time").mean(numeric_only=True)
        t = np.arange(int(g.index.max()) + 1)
        out[tid] = {c: np.interp(t, g.index.values, g[c].values) for c in ("Speed", "HR", "VO2", "VCO2")}
        out[tid]["subject"] = int(g.ID.iloc[0])
    return sel, out, n_valid


def prep(x, detrend):
    return x - uniform_filter1d(x, 60, mode="nearest") if detrend else x.copy()


def signals(tid, pair, detrend):
    a, b = PAIRS[pair]
    return prep(DATA[tid][a], detrend), prep(DATA[tid][b], detrend)


def pair_xy(i, j, pair, detrend):
    x, _ = signals(IDS[i], pair, detrend); _, y = signals(IDS[j], pair, detrend)
    n = min(len(x), len(y)); return x[:n], y[:n]


# ------------------------------------------------------------------ references (cross-correlation of first differences)
def xc_lag(x, y, h, diff=False):
    """Lag k in [-h,h] maximising |corr(x[t], y[t+k])|; returns (delay=k, max|C|). Positive = x leads y."""
    if diff:
        x, y = np.diff(x), np.diff(y)
    if x.std() == 0 or y.std() == 0:
        return None, 0.0
    x = (x - x.mean()) / x.std(); y = (y - y.mean()) / y.std(); n = len(x)
    ks = np.arange(-h, h + 1)
    cc = np.array([np.mean(x[max(0, -k):n - max(0, k)] * y[max(0, k):n - max(0, -k)]) for k in ks])
    m = int(np.argmax(np.abs(cc)))
    return float(ks[m]), float(abs(cc[m]))


def step_delay(speed, hr):
    """Median time from each speed step (>0.3 km/h within 5 s) until smoothed HR reaches 63 % of its change; None if no steps."""
    ds = speed[5:] - speed[:-5]; idx = np.nonzero(np.abs(ds) > 0.3)[0]
    if len(idx) == 0:
        return None, 0
    ev = [idx[0]]
    for k in idx[1:]:
        if k - ev[-1] > 30: ev.append(k)
    ev = [e + 2 for e in ev]; h = uniform_filter1d(hr, 5, mode="nearest"); res = []
    for q, e in enumerate(ev):
        end = ev[q + 1] if q + 1 < len(ev) else len(hr)
        if e < 10 or end - e < 40: continue
        base = h[e - 10:e].mean(); lev = h[max(e + 20, end - 30):end].mean(); ch = lev - base
        if abs(ch) < 3: continue
        hit = np.nonzero((h[e:end] - base) / ch >= 0.63)[0]
        if len(hit): res.append(float(hit[0]))
    return (float(np.median(res)) if res else None), len(res)


# ------------------------------------------------------------------ scoring jobs
def run(p, x, y):
    tau = A.time_delay_interaction(x, y, p)[0]; lab = A.stable_label(tau, p)
    dl = A.delay_of(tau, lab)
    return float(A.tds_score(lab)), (None if dl is None else -dl), int(np.sum(lab)), int(len(tau))


def job(a):
    kind, p, i, j, pair, detrend = a
    x, y = pair_xy(i, j, pair, detrend)
    if kind == "xcorr":
        d, c = xc_lag(x, y, p, False)
        return {"score": 100.0 * c, "delay": d, "n_stable": None}
    s, dl, ns, nw = run(p, x, y)
    return {"score": s, "delay": dl, "n_stable": ns}


def jwin(args):
    i, pair, detrend = args
    x, y = pair_xy(i, i, pair, detrend); p = A.window_params([x, y]); return p.window, p.max_lag


def jrhat(args):
    i, pair, detrend = args
    x, y = pair_xy(i, i, pair, detrend)
    def fb(v):
        r = A.acf(v, len(v) // 4); k = np.nonzero(r < np.exp(-1.0))[0]
        return float(k[0]) if len(k) else float(len(r) - 1)
    return max(1, int(round(max(fb(x), fb(y)) / 2)))


# ------------------------------------------------------------------ stats
def auc(a, b):
    a, b = np.asarray(a)[:, None], np.asarray(b)[None, :]
    return float((np.sum(a > b) + 0.5 * np.sum(a == b)) / (a.size * b.size))


def boot_auc(a, b, n=2000, seed=0):
    rng = np.random.default_rng(seed); a, b = np.asarray(a), np.asarray(b)
    v = [auc(rng.choice(a, a.size), rng.choice(b, b.size)) for _ in range(n)]
    return [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]


def metrics(true, ctrl, ref, p, fail, lagmax):
    a = [r["score"] for r in true]; b = [r["score"] for r in ctrl]; thr = float(np.percentile(b, 95))
    err = [abs(r["delay"] - ref[k]) for k, r in enumerate(true) if r["delay"] is not None and ref[k] is not None]
    dl = [r["delay"] for r in true if r["delay"] is not None]
    return {"true_mean": float(np.mean(a)), "true_median": float(np.median(a)), "ctrl_mean": float(np.mean(b)),
            "ctrl_median": float(np.median(b)), "ctrl_95pct": thr, "auc": auc(a, b), "auc_ci95": boot_auc(a, b),
            "n_true_above_ctrl95": int(sum(x > thr for x in a)), "n_with_delay": len(dl),
            "delay_median_s": float(np.median(dl)) if dl else None,
            "delay_iqr_s": [float(np.percentile(dl, 25)), float(np.percentile(dl, 75))] if dl else None,
            "n_delay_at_lag_edge": int(sum(abs(d) >= lagmax for d in dl)),
            "n_delay_vs_ref": len(err), "median_abs_err_vs_ref_s": float(np.median(err)) if err else None,
            "n_within_5s_of_ref": int(sum(e <= TOL_S for e in err)),
            "params_samples": p, "params_seconds": p, "calibration_failed": bool(fail)}


# ------------------------------------------------------------------ main
def settings(pool, pair, detrend, which):
    n = len(IDS); ids = range(n)
    wl = pool.map(jwin, [(i, pair, detrend) for i in ids])
    L = int(np.median([w for w, _ in wl])); lag = int(np.median([g for _, g in wl]))
    p0 = replace(A.PUBLISHED, window=L, overlap=L // 2, max_lag=min(lag, (L - 1) // 2), tolerance=0)
    pa, fail = A.calib(p0, [pair_xy(i, i, pair, detrend) for i in ids])
    out = {"fixed": (A.PUBLISHED, False), "atds": (pa, fail)}
    if "rescaled" in which:
        r = int(np.median(pool.map(jrhat, [(i, pair, detrend) for i in ids])))
        out["rescaled"] = (replace(A.PUBLISHED, window=60 * r, overlap=30 * r, max_lag=30 * r, tolerance=r), False)
    if "xcorr" in which:
        out["xcorr"] = (pa.max_lag, False)
    return out


def pdict(kind, p):
    if kind == "xcorr":
        return {"max_lag": int(p)}
    return {"window": p.window, "step": p.overlap, "max_lag": p.max_lag, "tolerance": p.tolerance}


def main(raw):
    global IDS
    t0 = time.time()
    sel, DATA_, n_valid = load(raw); DATA.update(DATA_); IDS = sel; n = len(IDS); k = n // 3
    ctrl_j = [(i + k) % n for i in range(n)]
    same_subj = int(sum(DATA[IDS[i]]["subject"] == DATA[IDS[ctrl_j[i]]]["subject"] for i in range(n)))
    S = {"n_valid_tests": n_valid, "n_selected": n, "seed": SEED, "selected": IDS, "control_shift": k,
         "controls_same_subject": same_subj, "reference": "cross-correlation of first differences within +-60 s (not an annotation); "
         "tolerance for 'within' = 5 s", "preprocessing": "1 Hz interpolation, minus 60-s centred moving average (variant: none)"}
    records, csvrows, rows_out = [], [], {}
    with Pool(2) as pool:
        for pair in PAIRS:
            ref = [xc_lag(*pair_xy(i, i, pair, True), REFLAG, diff=True) for i in range(n)]
            refd = [r[0] for r in ref]
            sd = [step_delay(DATA[IDS[i]]["Speed"], DATA[IDS[i]]["HR"]) for i in range(n)] if pair == "Speed_HR" else None
            S[f"{pair}_ref_delay_median_s"] = float(np.median([d for d in refd if d is not None]))
            S[f"{pair}_ref_delay_iqr_s"] = [float(np.percentile([d for d in refd if d is not None], q)) for q in (25, 75)]
            if sd:
                v = [s[0] for s in sd if s[0] is not None]
                S["Speed_HR_step_response"] = {"n_tests_with_step_delay": len(v), "median_s": float(np.median(v)) if v else "n/a",
                                               "iqr_s": [float(np.percentile(v, 25)), float(np.percentile(v, 75))] if v else "n/a"}
            for detrend in (True, False):
                which = ("rescaled", "xcorr") if detrend else ()
                st = settings(pool, pair, detrend, which)
                for name, (p, fail) in st.items():
                    jobs = [(name if name == "xcorr" else "tds", p, i, i, pair, detrend) for i in range(n)] + \
                           [(name if name == "xcorr" else "tds", p, i, ctrl_j[i], pair, detrend) for i in range(n)]
                    res = pool.map(job, jobs, chunksize=10)
                    true, ctrl = res[:n], res[n:]
                    key = f"{pair}|{'detrended' if detrend else 'raw'}|{name}"
                    pd_ = pdict(name, p)
                    m = metrics(true, ctrl, refd, pd_, fail, pd_["max_lag"])
                    m["params_seconds"] = {a: float(b) for a, b in pd_.items()}
                    rows_out[key] = m
                    for i in range(n):
                        for typ, r, jj in (("true", true[i], i), ("control", ctrl[i], ctrl_j[i])):
                            csvrows.append({"pair": pair, "preproc": "detrended" if detrend else "raw", "setting": name, "type": typ,
                                            "x_test": IDS[i], "y_test": IDS[jj], "score": r["score"], "delay_s": r["delay"],
                                            "n_stable": r["n_stable"], "ref_delay_s": refd[i] if typ == "true" else None,
                                            "ref_xcorr_c": ref[i][1] if typ == "true" else None,
                                            "ref_step_delay_s": (sd[i][0] if sd and typ == "true" else None)})
                    print(key, pd_, "AUC %.3f" % m["auc"], flush=True)
    S["wall_seconds"] = time.time() - t0
    rd = os.path.join(HERE, "results"); os.makedirs(rd, exist_ok=True)
    json.dump({"summary": S, "settings": rows_out, "records": csvrows}, open(os.path.join(rd, "treadmill.json"), "w"), indent=1)
    with open(os.path.join(rd, "treadmill.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(csvrows[0])); w.writeheader(); w.writerows(csvrows)
    print(json.dumps({k: v for k, v in S.items() if k != "selected"}, indent=1))
    print("%-34s %6s %6s %6s %13s %5s %5s %7s %6s %5s" % ("setting", "tmean", "cmean", "c95", "AUC", "n>c95", "ndel", "medErr", "n<=5s", "edge"))
    for k_, m in rows_out.items():
        print("%-34s %6.1f %6.1f %6.1f %5.3f[%.2f,%.2f] %5d %5d %7s %6d %5d" % (
            k_, m["true_mean"], m["ctrl_mean"], m["ctrl_95pct"], m["auc"], *m["auc_ci95"], m["n_true_above_ctrl95"], m["n_with_delay"],
            "%.1f" % m["median_abs_err_vs_ref_s"] if m["median_abs_err_vs_ref_s"] is not None else "n/a", m["n_within_5s_of_ref"],
            m["n_delay_at_lag_edge"]))


if __name__ == "__main__":
    main(sys.argv[1])
