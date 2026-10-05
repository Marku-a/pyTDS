"""R1.4: pyTDS core (circular pbc_xcorr) vs the thesis port at the thesis setting TDSParams(window=60, overlap=30,
max_lag=30, tolerance=1), stability_anchor "any" and "end".

Real links for the 35 D1 subjects (35 x 140, flatten order F, link order of thesis_port) and the classic
benchmark (data/data-delta.txt, data-sigma.txt; thesis_port = MATLAB Time_Delay_Interaction/StableLabel).

Window-time mapping: pyTDS (window_anchor="center") returns t_vec = start + L//2 = start+30 as a 0-based index; the
thesis time is the 1-based sample (k-1)*30+30 = start+30 (1-based) = 0-based index start+29. We therefore use
stage = stages[t_vec - 1] (identical to thesis_port's stages[t-1]). pyTDS and the port produce the same number of
windows, floor(N/30)-1.

Usage: python r1_pytds_vs_port.py DATA.npz   -> ../results/r1_pytds_vs_port.json
"""
import json, os, sys
from concurrent.futures import ProcessPoolExecutor
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", ".."))
from thesis_port import _order, _good, links_strength, time_delay_interaction, stable_label, STAGES
from pyTDS.core import tds
from pyTDS.params import TDSParams

DATA = sys.argv[1]
OUT = os.path.join(HERE, "..", "results")
ANCHORS = ("any", "end")


def one(i):
    z = np.load(DATA); n = z["names"][i]
    X, S = z["X_" + n], z["S_" + n]
    res = {"port": links_strength(X, S).flatten(order="F")}
    for a in ANCHORS:
        P = TDSParams(window=60, overlap=30, max_lag=30, tolerance=1, stability_anchor=a)
        out = np.zeros((35, 4))
        for li, (p, q) in enumerate(_order(_good())):
            r = tds(X[:, p], X[:, q], P)
            st = S[r["t_vec"].astype(int) - 1]
            for k, s in enumerate(STAGES):
                m = st == s
                out[li, k] = 100 * r["stbl_lbl"][m].mean() if m.any() else np.nan
        res[a] = out.flatten(order="F")
    return res


def cmp(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    d = a[m] - b[m]
    return {"n": int(m.sum()), "mean_diff_pytds_minus_port": float(d.mean()), "mean_abs_diff": float(np.abs(d).mean()),
            "pearson_r": float(np.corrcoef(a[m], b[m])[0, 1]), "max_abs_diff": float(np.abs(d).max()),
            "nan_pytds": int(np.isnan(a).sum()), "nan_port": int(np.isnan(b).sum())}


if __name__ == "__main__":
    with ProcessPoolExecutor(4) as ex:
        res = list(ex.map(one, range(35)))
    port = np.concatenate([r["port"] for r in res])
    out = {"setting": "window=60, overlap=30, max_lag=30, tolerance=1", "real_links": {}}
    for a in ANCHORS:
        out["real_links"][a] = cmp(np.concatenate([r[a] for r in res]), port)
        # per-stage breakdown (columns 35*k..35*k+34 per subject)
        pa = np.array([r[a] for r in res]); pp = np.array([r["port"] for r in res])
        out["real_links"][a]["per_stage_mean_diff"] = {s: float(np.nanmean(pa[:, 35 * k:35 * k + 35] - pp[:, 35 * k:35 * k + 35]))
                                                       for k, s in enumerate(STAGES)}
    data = os.path.join(HERE, "..", "..", "data")
    delta = np.loadtxt(os.path.join(data, "data-delta.txt"))[:, 1]
    sigma = np.loadtxt(os.path.join(data, "data-sigma.txt"))[:, 1]
    n = min(len(delta), len(sigma))
    sigma, delta = sigma[:n], delta[:n]       # sigma has one extra row; truncated to common length, tds(sigma, delta)
    tau, _ = time_delay_interaction(sigma, delta)
    b = {"n_samples": int(n), "port_score": float(100 * stable_label(tau).mean())}
    for a in ANCHORS:
        b["pytds_" + a] = float(tds(sigma, delta, TDSParams(window=60, overlap=30, max_lag=30, tolerance=1, stability_anchor=a))["score"])
    b["pytds_default_params"] = float(tds(sigma, delta)["score"])
    b["original_matlab"] = 43.6
    out["benchmark"] = b
    json.dump(out, open(os.path.join(OUT, "r1_pytds_vs_port.json"), "w"), indent=1)
    print(json.dumps(out, indent=1))
