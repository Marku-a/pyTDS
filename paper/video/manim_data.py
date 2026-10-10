"""Data for the 3Blue1Brown-style aTDS video, computed with the real pyTDS / aTDS code.

Synthetic signals are generated here (seeded). Real-data panels use committed aggregates
(paper/figures/figdata.json, reader.json) and the R5 worked examples (examples.json, kept out of git).

usage: python manim_data.py --examples examples.json --out video_data.json
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path[:0] = [ROOT, os.path.join(ROOT, "experiments", "adaptive_groundtruth"), os.path.join(ROOT, "paper", "scripts")]
from adaptive import acf, adaptive_params, bartlett_factor  # noqa: E402
from figdata import r  # noqa: E402
from pyTDS.core import tds  # noqa: E402
from pyTDS.params import TDSParams  # noqa: E402
from pyTDS.utils import pbc_xcorr, zscore  # noqa: E402

PUB = TDSParams(window=60, overlap=30, max_lag=30, tolerance=1)


def ar(n, tc, rng):
    a = np.exp(-1 / tc)
    x = np.zeros(n)
    e = rng.standard_normal(n)
    for i in range(1, n):
        x[i] = a * x[i - 1] + e[i]
    return x / x.std()


def coupled(n, tc, delay, c, seed):
    rng = np.random.default_rng(seed)
    x = ar(n, tc, rng)
    y = c * np.roll(x, delay) + np.sqrt(1 - c * c) * ar(n, tc, rng)
    return x, y


def lag_view(res):
    """pyTDS tau < 0 means s1 leads; return 'how much signal 2 lags signal 1' = -tau."""
    return {"centre": np.asarray(res["t_vec"]), "lag": -np.asarray(res["tau"]), "stable": np.asarray(res["stbl_lbl"]),
            "score": float(res["score"])}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--examples", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = {}

    # Ch2-3: a fast coupled pair (memory ~4 s, signal 2 follows by 6 s), 900 s
    x, y = coupled(900, 4, 6, 0.8, 11)
    s = 300
    lags, C = pbc_xcorr(zscore(x[s:s + 60]), zscore(y[s:s + 60]), 30)
    out["fast"] = {"x": x, "y": y, "win_start": s, "xc_lag": -np.asarray(lags), "xc_C": C, "true_delay": 6,
                   "tds": lag_view(tds(x, y, PUB))}

    # Ch4, Ch9: a slow coupled pair (memory ~30 s, delay 10 s, coupling 0.5), 4 h
    xs, ys = coupled(14400, 30, 10, 0.5, 2)
    p, info = adaptive_params([xs, ys], window_rule="bartlett", tolerance="calibrated", base=PUB, max_window=1200,
                              n_shifts=10, max_pairs=10, seed=0)
    nulls = {k: [tds(xs, np.roll(ys, 3000 + i * 1500), q)["score"] for i in range(5)] for k, q in (("tds", PUB), ("atds", p))}
    out["slow"] = {"x": xs[:900], "y": ys[:900], "true_delay": 10, "B": bartlett_factor(xs),
                   "tds": lag_view(tds(xs, ys, PUB)), "atds": lag_view(tds(xs, ys, p)),
                   "atds_params": {"window": p.window, "step": p.overlap, "max_lag": p.max_lag, "tolerance": p.tolerance},
                   "null_mean": {k: float(np.mean(v)) for k, v in nulls.items()},
                   "null_by_tol": {str(k): v for k, v in info["null_score_by_tolerance"].items()},
                   "raw_window": float(30 * bartlett_factor(xs))}

    # Ch5: memory of a fast and a slow signal
    xf = ar(3000, 2, np.random.default_rng(5))
    out["memory"] = {"fast": {"x": xf[:300], "acf": acf(xf, 120), "B": bartlett_factor(xf)},
                     "slow": {"x": xs[:300], "acf": acf(xs, 120), "B": bartlett_factor(xs)}}

    # Ch9-10: committed aggregates and the worked examples
    reader = json.load(open(os.path.join(ROOT, "paper", "figures", "reader.json")))
    fig = json.load(open(os.path.join(ROOT, "paper", "figures", "figdata.json")))
    ex = json.load(open(a.examples))
    out["e1"] = reader["synthetic"]["E1_weak_coupling_vs_slowness"]
    conn = {}
    for v in ("published", "v2_calibrated"):
        conn[v] = {st: float(np.nanmean([q for q in fig["connectivity"][v]["strict"][st] if q is not None])) for st in ("LS", "awake", "REM", "DS")}
    out["sleep_stage_means"] = conn
    g = ex["garmin"]
    out["garmin"] = {k: {"centre": m["centre"], "lag": [-q for q in m["lag_hr_leads_speed_s"]], "stable": m["stable"],
                         "params": m["params"], "detected": m["detected"], "p": m["p"]} for k, m in g["methods"].items()}
    out["garmin"]["length_min"] = g["length_min"]
    json.dump(r(out, 5), open(a.out, "w"), separators=(",", ":"))
    sl = out["slow"]
    print("fast TDS score", round(out["fast"]["tds"]["score"], 1), "| slow B", round(sl["B"], 1), "TDS", round(sl["tds"]["score"], 1),
          "aTDS", round(sl["atds"]["score"], 1), sl["atds_params"], "null", sl["null_mean"], "| memory B fast/slow",
          round(out["memory"]["fast"]["B"], 1), round(out["memory"]["slow"]["B"], 1))


if __name__ == "__main__":
    main()
