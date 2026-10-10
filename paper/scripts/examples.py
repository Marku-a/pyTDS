"""R5 worked examples on actual signals: how TDS and aTDS process the same real data.

The examples are chosen by the rule fixed in ANALYSIS_PLAN.md (R5) before they were computed.
Output contains short raw traces, so it is written OUTSIDE git (the artifact publishes it privately).

usage: PYTHONPATH=<pyTDS> python examples.py --ncom NCOM35.npz --garmin-repo <GarminTDS> \
           --garmin-npz runs.npz --garmin-rows rows_real.csv --out examples.json
"""
import argparse
import json
import os
import sys
from dataclasses import replace

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import r2_run  # noqa: E402  (also puts experiments/adaptive_groundtruth and the repo root on sys.path)
from adaptive import acf, adaptive_params, bartlett_factor  # noqa: E402
from common import published_params, stable_seed, variant_params  # noqa: E402
from figdata import r  # noqa: E402
from pyTDS.core import tds  # noqa: E402
from pyTDS.params import TDSParams  # noqa: E402
from pyTDS.utils import pbc_xcorr, zscore  # noqa: E402

RES = os.path.join(HERE, "..", "results")
SPAN, PRE = 2400, 600


def xcorr_at(x, y, centre, p):
    s = int(centre) - p.window // 2
    lags, C = pbc_xcorr(zscore(x[s:s + p.window]), zscore(y[s:s + p.window]), p.max_lag)
    return {"centre": int(centre), "lags": lags, "C": C}


def sleep_example(ncom):
    sc = pd.read_csv(os.path.join(RES, "r2_scores.csv"))
    t = sc[(sc.rule == "strict") & (sc.thesis_link == 1) & sc.variant.isin(["published", "v2_calibrated"])]
    piv = t.pivot_table(index=["subject", "pair"], columns=["variant", "stage"], values="excess")
    d_ls_ds = {v: piv[(v, "LS")] - piv[(v, "DS")] for v in ("published", "v2_calibrated")}
    link_gain = (d_ls_ds["v2_calibrated"] - d_ls_ds["published"]).rename("gain").reset_index()
    conn = {v: d_ls_ds[v].groupby("subject").mean() for v in d_ls_ds}       # mean over links (NaN-skipping)
    subj_gain = (conn["v2_calibrated"] - conn["published"]).dropna()
    med = float(subj_gain.median())
    s_idx = int((subj_gain - med).abs().sort_index().idxmin())               # tie -> lower index
    lg = link_gain[link_gain.subject == s_idx].dropna(subset=["gain"]).set_index("pair").gain
    lmed = float(lg.median())
    pair = (lg - lmed).abs().idxmin()
    i, j = map(int, pair.split("-"))

    z = np.load(ncom)
    name = str(z["names"][s_idx])
    X, S = np.asarray(z["X_" + name], float), np.asarray(z["S_" + name]).astype(str)
    T = min(len(X), len(S))
    X, S = X[:T], S[:T]
    code = {s: k for k, s in enumerate(r2_run.STAGES)}
    st = np.array([code.get(s, -1) for s in S])
    trans = [k for k in range(PRE, T) if S[k - 1] == "LS" and S[k] == "DS"]
    t0 = trans[0]
    start = (t0 - PRE) // 150 * 150
    start = min(start, T - SPAN)
    end = start + SPAN

    sig = [X[:, k] for k in range(10)]
    seed = stable_seed("ncom35", s_idx)
    vp = variant_params(sig, fs=1.0, max_window_s=300, seed=seed)
    prm = pd.read_csv(os.path.join(RES, "r2_params.csv"))
    methods = {}
    for v, label in (("published", "TDS"), ("v2_calibrated", "aTDS")):
        p = vp[v][0]
        row = prm[(prm.subject == s_idx) & (prm.variant == v)].iloc[0]
        assert (p.window, p.overlap, p.max_lag, p.tolerance) == (row.window, row.step, row.max_lag, row.tolerance)
        res = tds(X[:, i], X[:, j], p)
        c = np.asarray(res["t_vec"])
        ws = r2_run.window_stage(c, st, p.window)
        strict, _ = r2_run.stage_masks(ws)
        lab = np.full(len(c), -1)
        for k in range(4):
            lab[strict[k]] = k
        keep = (c >= start) & (c < end)
        picks = [start + 300, t0, t0 + 600]
        xc = [xcorr_at(X[:, i], X[:, j], c[np.argmin(np.abs(c - q))], p) for q in picks]
        stage_scores = sc[(sc.subject == s_idx) & (sc.pair == pair) & (sc.variant == v) & (sc.rule == "strict")]
        methods[label] = {"variant": v, "params": {"window": p.window, "step": p.overlap, "max_lag": p.max_lag,
                                                   "tolerance": p.tolerance},
                          "centre": c[keep], "tau": np.asarray(res["tau"])[keep],
                          "stable": np.asarray(res["stbl_lbl"])[keep], "strict_stage": lab[keep],
                          "xcorr": xc,
                          "stage_scores": stage_scores[["stage", "n_windows", "real", "chance", "excess"]].to_dict("records"),
                          "whole_night_score": float(res["score"]),
                          "night": {"centre": c, "tau": np.asarray(res["tau"]), "stable": np.asarray(res["stbl_lbl"]),
                                    "strict_stage": lab}}

    _, info = adaptive_params(sig, window_rule="bartlett", tolerance="calibrated", base=published_params(1.0),
                              max_window=300, n_shifts=10, max_pairs=10, seed=seed)
    deriv = {"cap": 300, "raw_window_by_signal": dict(zip(r2_run.SIGNALS, [round(30 * b) for b in info["bartlett"]])),
             "null_score_by_tolerance": {str(k): v for k, v in info["null_score_by_tolerance"].items()},
             "chosen_tolerance": info["tolerance"], "alpha_pct": 5.0,
             "acf": {r2_run.SIGNALS[k]: acf(X[:, k], 600) for k in (i, j)},
             "bartlett_factor": {r2_run.SIGNALS[k]: bartlett_factor(X[:, k]) for k in (i, j)}}
    return {"rule": {"subject_gain_median": med, "subject_gain": float(subj_gain[s_idx]),
                     "n_subjects_eligible": int(len(subj_gain)), "link_gain_median": lmed,
                     "link_gain": float(lg[pair])},
            "subject": s_idx, "pair": [r2_run.SIGNALS[i], r2_run.SIGNALS[j]],
            "span": [int(start), int(end)], "transition": int(t0),
            "stage": S[start:end].tolist(),
            "hypnogram_30s": S[::30].tolist(), "x": X[start:end, i], "y": X[start:end, j],
            "methods": methods, "atds_derivation": deriv}


def garmin_example(repo, npz, rows_csv):
    sys.path.insert(0, repo)
    from experiments.adaptive_compare import compare as C   # noqa: E402
    from garmintds import lag                               # noqa: E402
    from garmintds.network import span_segments             # noqa: E402

    rows = pd.read_csv(rows_csv)
    x = rows[(rows.kind != "xrun") & (rows.a == "hr") & (rows.b == "speed") &
             rows.variant.isin(["frozen", "v2_calibrated"])].copy()
    x["det"] = (x.p < 0.05) & (x.lag_s < 0) & (~x.pinned.astype(bool))
    piv = x.pivot(index="run", columns="variant", values="det")
    dis = piv[piv.frozen != piv.v2_calibrated].index.astype(str).tolist()
    runs = {ru.id: ru for ru in C.load_real(type(C.Path())(npz))}
    lens = sorted((runs[d].n, d) for d in dis)
    n_len, rid = lens[(len(lens) - 1) // 2]                   # median; even count -> shorter
    run = runs[rid]
    out = {"n_disagreeing_runs": len(dis), "length_min": round(run.n / 60, 1), "methods": {}}
    hr, sp = run.grid["hr"], run.grid["speed"]
    for v, label in (("frozen", "TDS (app setting)"), ("v2_calibrated", "aTDS")):
        row = x[(x.run.astype(str) == rid) & (x.variant == v)].iloc[0]
        p = replace(lag.PARAMS, window=int(row.window), overlap=int(row.step), max_lag=int(row.max_lag),
                    tolerance=int(row.tolerance))
        segs = span_segments(run, "hr", "speed", (0, run.n), p.window + 4 * p.overlap)
        rr = lag.lead_lag(hr, sp, segs, p)
        assert abs(rr["score"] - row.score) < 1e-9, (rr["score"], row.score)
        null = lag.chance_scores(hr, sp, segs, p, n=100, seed=C.stable_seed(rid, "hr", "speed"))
        wl = lag.window_labels(hr, sp, segs, p)
        out["methods"][label] = {"variant": v, "detected": bool(row.det),
                                 "params": {"window": p.window, "step": p.overlap, "max_lag": p.max_lag,
                                            "tolerance": p.tolerance},
                                 "score": rr["score"], "lag_s": rr["lag_s"], "p": float(row.p),
                                 "chance95": float(row.chance95), "null_scores": null,
                                 "centre": [w[0] for w in wl], "stable": [w[1] for w in wl],
                                 "lag_hr_leads_speed_s": [w[2] for w in wl]}
        out["segments"] = [list(map(int, s)) for s in segs]
    out["hr"], out["speed"] = hr, sp
    return out


def main():
    ap = argparse.ArgumentParser()
    for k in ("--ncom", "--garmin-repo", "--garmin-npz", "--garmin-rows", "--out"):
        ap.add_argument(k, required=True)
    a = ap.parse_args()
    out = {"sleep": sleep_example(a.ncom), "garmin": garmin_example(a.garmin_repo, a.garmin_npz, a.garmin_rows)}
    json.dump(r(out, 5), open(a.out, "w"), separators=(",", ":"))
    s, g = out["sleep"], out["garmin"]
    print(a.out, os.path.getsize(a.out) // 1024, "KB")
    print("sleep: subject", s["subject"], "pair", s["pair"], "rule", r(s["rule"]), "span", s["span"],
          {k: m["params"] for k, m in s["methods"].items()})
    print("garmin:", g["length_min"], "min;", {k: (m["detected"], round(m["score"], 1), m["p"])
                                                for k, m in g["methods"].items()})


if __name__ == "__main__":
    main()
