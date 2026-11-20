"""Gait in Parkinson's disease (PhysioNet gaitpdb 1.0.0): L/R total force and within-foot heel->toe, fixed vs rescaled vs aTDS vs xcorr.

    python experiments/physio_new/gait.py <raw_dir>

See SPEC_COMMON.md / SPEC_gait.md. TDS via pyTDS only (through experiments/atds_real/atds21.py).
Preprocessing: z-scoring only (100 Hz). Delay convention: positive = first signal leads (= -tau0); delays in samples/100 = s.
Heel/toe sensor of the left foot: the sensor layout is not in the files, so the heel sensor = the column whose rise (20 % of its own
max) comes earliest after the left heel strike and the toe sensor = the latest, taken as the median over all subjects (one choice).
"""
import csv
import json
import os
import sys
import time
from dataclasses import replace

import numpy as np
from scipy.stats import mannwhitneyu

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "atds_real"))
import atds21 as A  # noqa: E402

A.WORKERS = 1
FS, SEED, TOL_S = 100.0, 0, 0.05
PAIRS = {"LR_total": 0, "heel_toe": 1}
DATA: dict = {}
IDS: list = []


# ------------------------------------------------------------------ data
def load(d):
    names = sorted(f for f in os.listdir(d) if f.endswith("_01.txt") and f[:2] in ("Ga", "Ju", "Si"))
    out = {}
    for f in names:
        a = np.loadtxt(os.path.join(d, f))
        if a.ndim != 2 or a.shape[1] < 19:
            continue
        out[f[:-4]] = {"L": a[:, 1:9].astype(np.float32), "totL": a[:, 17], "totR": a[:, 18],
                       "group": "control" if "Co" in f[:6] else "parkinson"}
    return sorted(out), out


def upcross(x, thr):
    return np.nonzero((x[:-1] < thr) & (x[1:] >= thr))[0] + 1


def downcross(x, thr):
    return np.nonzero((x[:-1] >= thr) & (x[1:] < thr))[0] + 1


def stances(tot):
    """(heel strike, toe off) index pairs of one foot from 20 % of max total force."""
    thr = 0.2 * tot.max(); up = upcross(tot, thr); dn = downcross(tot, thr); res = []
    for u in up:
        k = dn[dn > u]
        if len(k): res.append((u, k[0]))
    return up, res


def sensor_rises(s):
    """Per stance, rise index (20 % of sensor max) of each of the 8 left sensors relative to heel strike (nan if none)."""
    up, st = stances(s["totL"]); thr = 0.2 * s["L"].max(0); rows = []
    for u, e in st:
        seg = s["L"][u:e]; r = []
        for c in range(8):
            h = np.nonzero(seg[:, c] >= thr[c])[0]; r.append(h[0] if len(h) else np.nan)
        rows.append(r)
    return np.array(rows, float).reshape(-1, 8)


def choose_sensors():
    med = np.array([np.nanmedian(r, 0) for r in (sensor_rises(DATA[k]) for k in IDS) if len(r) > 3])
    m = np.nanmedian(med, 0)
    return int(np.nanargmin(m)), int(np.nanargmax(m)), m


def wrap(d, S):
    return (d + S / 2.0) % S - S / 2.0


def reference(s, heel, toe):
    """Exact reference (samples): L->R delay, stride time, heel->toe delay (all medians over the record)."""
    upL, st = stances(s["totL"]); upR, _ = stances(s["totR"])
    stride = float(np.median(np.diff(upL))) if len(upL) > 2 else None
    lr = [upR[upR > u][0] - u for u in upL if (upR > u).any()]
    ht = []
    for u, e in st:
        seg = s["L"][u:e]; thr = 0.2 * s["L"].max(0)
        h = np.nonzero(seg[:, heel] >= thr[heel])[0]; t = np.nonzero(seg[:, toe] >= thr[toe])[0]
        if len(h) and len(t): ht.append(t[0] - h[0])
    return {"stride": stride, "lr": float(np.median(lr)) if lr else None, "ht": float(np.median(ht)) if ht else None,
            "n_strides": int(len(upL))}


def zs(v):
    v = np.asarray(v, float); sd = v.std()
    return (v - v.mean()) / (sd if sd > 0 else 1.0)


def signal(key, pair, heel, toe, which):
    s = DATA[key]
    if pair == "LR_total":
        return s["totL"] if which == 0 else s["totR"]
    return s["L"][:, heel if which == 0 else toe].astype(float)


def pair_xy(i, j, pair):
    x = signal(IDS[i], pair, HEEL, TOE, 0); y = signal(IDS[j], pair, HEEL, TOE, 1)
    n = min(len(x), len(y)); return zs(x[:n]), zs(y[:n])


# ------------------------------------------------------------------ scoring
def xc_lag(x, y, h):
    """Lag k in [-h,h] maximising |C(k)|, C(k)=mean(x[t] y[t+k]) of the z-scored pair; positive = x leads y."""
    n = len(x); f = 1 << int(np.ceil(np.log2(2 * n)))
    c = np.fft.irfft(np.conj(np.fft.rfft(x, f)) * np.fft.rfft(y, f), f)
    cc = np.concatenate([c[-h:], c[:h + 1]]) / n; ks = np.arange(-h, h + 1); m = int(np.argmax(np.abs(cc)))
    return float(ks[m]), float(abs(cc[m]))


def job(kind, p, i, j, pair):
    x, y = pair_xy(i, j, pair)
    if kind == "xcorr":
        d, c = xc_lag(x, y, p); return {"score": 100.0 * c, "delay": d, "n_stable": None}
    tau = A.time_delay_interaction(x, y, p)[0]; lab = A.stable_label(tau, p); dl = A.delay_of(tau, lab)
    return {"score": float(A.tds_score(lab)), "delay": None if dl is None else -dl, "n_stable": int(np.sum(lab))}


def fb(v):
    r = A.acf(v, len(v) // 4); k = np.nonzero(r < np.exp(-1.0))[0]
    return float(k[0]) if len(k) else float(len(r) - 1)


def auc(a, b):
    a, b = np.asarray(a)[:, None], np.asarray(b)[None, :]
    return float((np.sum(a > b) + 0.5 * np.sum(a == b)) / (a.size * b.size))


def boot_auc(a, b, n=2000, seed=0):
    rng = np.random.default_rng(seed); a, b = np.asarray(a), np.asarray(b)
    v = [auc(rng.choice(a, a.size), rng.choice(b, b.size)) for _ in range(n)]
    return [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]


def metrics(true, ctrl, ref, strides, p, fail, lagmax):
    a = [r["score"] for r in true]; b = [r["score"] for r in ctrl]; thr = float(np.percentile(b, 95))
    err = [abs(wrap(r["delay"] - ref[k], strides[k])) / FS for k, r in enumerate(true)
           if r["delay"] is not None and ref[k] is not None and strides[k]]
    dl = [r["delay"] / FS for r in true if r["delay"] is not None]
    return {"true_mean": float(np.mean(a)), "true_median": float(np.median(a)), "ctrl_mean": float(np.mean(b)),
            "ctrl_median": float(np.median(b)), "ctrl_95pct": thr, "auc": auc(a, b), "auc_ci95": boot_auc(a, b),
            "n_true_above_ctrl95": int(sum(x > thr for x in a)), "n_with_delay": len(dl),
            "delay_median_s": float(np.median(dl)) if dl else None,
            "n_delay_at_lag_edge": int(sum(abs(d) * FS >= lagmax for d in dl)),
            "n_delay_vs_ref": len(err), "median_abs_err_vs_ref_s": float(np.median(err)) if err else None,
            "n_within_tol_of_ref": int(sum(e <= TOL_S for e in err)), "tol_s": TOL_S,
            "params_samples": p, "params_seconds": {k: float(v) / FS for k, v in p.items()}, "calibration_failed": bool(fail)}


def pdict(kind, p):
    if kind == "xcorr":
        return {"max_lag": int(p)}
    return {"window": p.window, "step": p.overlap, "max_lag": p.max_lag, "tolerance": p.tolerance}


def settings(pair):
    n = len(IDS)
    wl = []
    for i in range(n):
        p = A.window_params(list(pair_xy(i, i, pair))); wl.append((p.window, p.max_lag))
    L = int(np.median([w for w, _ in wl])); lag = int(np.median([g for _, g in wl]))
    p0 = replace(A.PUBLISHED, window=L, overlap=L // 2, max_lag=min(lag, (L - 1) // 2), tolerance=0)
    pa, fail = A.calib(p0, [pair_xy(i, i, pair) for i in range(n)])
    rh = [max(1, int(round(max(fb(x), fb(y)) / 2))) for x, y in (pair_xy(i, i, pair) for i in range(n))]
    r = int(np.median(rh))
    return {"fixed": (A.PUBLISHED, False),
            "rescaled": (replace(A.PUBLISHED, window=60 * r, overlap=30 * r, max_lag=30 * r, tolerance=r), False),
            "atds": (pa, fail), "xcorr": (pa.max_lag, False)}, {"r_hat": r, "per_subject_window_median": L,
                                                               "per_subject_lag_median": lag}


def main(raw):
    global IDS, HEEL, TOE
    t0 = time.time()
    IDS, D = load(raw); DATA.update(D); n = len(IDS); k = n // 3
    HEEL, TOE, rise_med = choose_sensors()
    ctrl_j = [(i + k) % n for i in range(n)]
    refs = [reference(DATA[key], HEEL, TOE) for key in IDS]
    strides = [r["stride"] for r in refs]
    S = {"n_subjects": n, "n_control": sum(DATA[x]["group"] == "control" for x in IDS), "control_shift": k, "seed": SEED,
         "heel_sensor_L": HEEL + 1, "toe_sensor_L": TOE + 1, "median_rise_offset_per_sensor_samples": [float(v) for v in rise_med],
         "median_record_samples": float(np.median([len(DATA[x]["totL"]) for x in IDS])),
         "median_stride_s": float(np.median([s for s in strides if s])) / FS, "tol_s": TOL_S,
         "controls_same_subject": 0, "reference": "20 % of max force crossings; delays compared modulo stride time"}
    records, csvrows, rows_out, extra = [], [], {}, {}
    grp = np.array([DATA[x]["group"] for x in IDS])
    for pair in PAIRS:
        rk = "lr" if pair == "LR_total" else "ht"
        ref = [r[rk] for r in refs]
        okr = [v for v in ref if v is not None]
        S[f"{pair}_ref_delay_median_s"] = float(np.median(okr)) / FS
        S[f"{pair}_ref_delay_iqr_s"] = [float(np.percentile(okr, q)) / FS for q in (25, 75)]
        st, info = settings(pair); S[f"{pair}_scale_info"] = info
        for name, (p, fail) in st.items():
            kind = "xcorr" if name == "xcorr" else "tds"
            true = [job(kind, p, i, i, pair) for i in range(n)]
            ctrl = [job(kind, p, i, ctrl_j[i], pair) for i in range(n)]
            pd_ = pdict(name, p); m = metrics(true, ctrl, ref, strides, pd_, fail, pd_["max_lag"])
            if name == "atds":
                m["lag_over_median_stride"] = pd_["max_lag"] / np.median([s for s in strides if s])
                m["rhythm_cap_applies"] = bool(pd_["max_lag"] <= 0.6 * np.median([s for s in strides if s]))
            # controls (healthy) vs Parkinson's comparison of delays found and scores
            g = {}
            for gn in ("control", "parkinson"):
                sel = [r for r, gg in zip(true, grp) if gg == gn]
                dd = [r["delay"] / FS for r in sel if r["delay"] is not None]
                g[gn] = {"n": len(sel), "score_median": float(np.median([r["score"] for r in sel])),
                         "n_delay": len(dd), "delay_median_s": float(np.median(dd)) if dd else None,
                         "delay_iqr_s": [float(np.percentile(dd, q)) for q in (25, 75)] if dd else None}
            sc = ([r["score"] for r, gg in zip(true, grp) if gg == "control"], [r["score"] for r, gg in zip(true, grp) if gg == "parkinson"])
            g["score_mannwhitney_p"] = float(mannwhitneyu(*sc).pvalue)
            d1 = [r["delay"] / FS for r, gg in zip(true, grp) if gg == "control" and r["delay"] is not None]
            d2 = [r["delay"] / FS for r, gg in zip(true, grp) if gg == "parkinson" and r["delay"] is not None]
            g["delay_mannwhitney_p"] = float(mannwhitneyu(d1, d2).pvalue) if len(d1) > 2 and len(d2) > 2 else None
            m["group_comparison"] = g
            key = f"{pair}|{name}"; rows_out[key] = m
            for i in range(n):
                for typ, r, jj in (("true", true[i], i), ("control", ctrl[i], ctrl_j[i])):
                    csvrows.append({"pair": pair, "setting": name, "type": typ, "x_subject": IDS[i], "y_subject": IDS[jj],
                                    "group": DATA[IDS[i]]["group"], "score": r["score"],
                                    "delay_s": None if r["delay"] is None else r["delay"] / FS, "n_stable": r["n_stable"],
                                    "ref_delay_s": ref[i] / FS if typ == "true" and ref[i] is not None else None,
                                    "stride_s": strides[i] / FS if strides[i] else None})
            print(key, pd_, "AUC %.3f" % m["auc"], flush=True)
    S["wall_seconds"] = time.time() - t0
    rd = os.path.join(HERE, "results"); os.makedirs(rd, exist_ok=True)
    json.dump({"summary": S, "settings": rows_out, "reference_rows": {IDS[i]: refs[i] for i in range(n)}, "records": csvrows},
              open(os.path.join(rd, "gait.json"), "w"), indent=1)
    with open(os.path.join(rd, "gait.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(csvrows[0])); w.writeheader(); w.writerows(csvrows)
    print(json.dumps(S, indent=1))
    print("%-18s %6s %6s %6s %13s %5s %5s %7s %6s %5s" % ("setting", "tmean", "cmean", "c95", "AUC", "n>c95", "ndel", "medErr", "n<=tol", "edge"))
    for k_, m in rows_out.items():
        print("%-18s %6.1f %6.1f %6.1f %5.3f[%.2f,%.2f] %5d %5d %7s %6d %5d" % (
            k_, m["true_mean"], m["ctrl_mean"], m["ctrl_95pct"], m["auc"], *m["auc_ci95"], m["n_true_above_ctrl95"], m["n_with_delay"],
            "%.3f" % m["median_abs_err_vs_ref_s"] if m["median_abs_err_vs_ref_s"] is not None else "n/a", m["n_within_tol_of_ref"],
            m["n_delay_at_lag_edge"]))
    for k_, m in rows_out.items():
        print(k_, json.dumps(m["group_comparison"]))


if __name__ == "__main__":
    main(sys.argv[1])
