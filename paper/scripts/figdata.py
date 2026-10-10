"""Figure data for the aTDS paper.

Two outputs:
- --agg-out (committed, public): aggregates only, from paper/results and the stored thesis/reproduction outputs.
  No time series, subject names, Garmin activity ids or dates.
- --ex-out  (kept out of git): the worked examples on actual signals (R5), written by examples.py.

usage: python figdata.py --msc <msc-project> --agg-out ../figures/figdata.json
"""
import argparse
import json
import os

import h5py
import numpy as np
import pandas as pd
import scipy.io as sio

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
SIGNALS = ["HR", "Resp", "Chin", "Leg", "Eye", "delta", "theta", "alpha", "sigma", "beta"]
STAGES = ["LS", "DS", "REM", "awake"]          # R2 order
MAT_STAGES = ["DS", "LS", "REM", "awake"]      # reproduction (GetLinksStrength) column-block order
VARIANTS = ["published", "tol_only", "win_only", "v2_calibrated", "v2_cal_lag3", "v2_cap120", "v2_cap600"]


def r(x, sig=4):
    """Round floats (recursively) to `sig` significant digits; NaN -> None."""
    if isinstance(x, dict):
        return {k: r(v, sig) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [r(v, sig) for v in x]
    if isinstance(x, np.ndarray):
        return r(x.tolist(), sig)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (float, np.floating)):
        x = float(x)
        if not np.isfinite(x):
            return None
        return float(f"{x:.{sig}g}")
    return x


def replication(msc):
    m = sio.loadmat(f"{msc}/01 repreduce NCOM/links_strength2000Sur.mat")
    real, sur_st = m["real_links"], m["sur_links"]
    sur_port = np.load(os.path.join(RES, "r1_sur_links_port.npy"))
    thr = json.load(open(os.path.join(RES, "r1_threshold.json")))
    pv = json.load(open(os.path.join(RES, "r1_pytds_vs_port.json")))
    rank = {}
    for k, st in enumerate(MAT_STAGES):
        sl = slice(35 * k, 35 * (k + 1))
        mr = real[:, sl].mean(0)
        order = np.argsort(-mr)
        rank[st] = {"real": mr[order], "stored_sur": np.sort(sur_st[:, sl].mean(0))[::-1],
                    "port_sur": np.sort(sur_port[:, sl].mean(0))[::-1],
                    "stored_sur_p95": np.percentile(sur_st[:, sl], 95, axis=0)[order],
                    "port_sur_p95": np.percentile(sur_port[:, sl], 95, axis=0)[order],
                    "real_subject_sd": real[:, sl].std(0, ddof=1)[order]}
    grid = np.round(1 + 0.1 * np.arange(191), 1)
    return {"stage_order": MAT_STAGES, "rank_curves": rank,
            "threshold_curves": {"grid": grid, "stored": m["signi_over"].ravel(),
                                 "port_surrogates_all_data_clean": thr["port_surrogates"]["signi_over"]},
            "thresholds": {"stored": float(m["thres"].ravel()[0]),
                           "stage_cells_with_port_surrogates": thr["thesis_stage_cells_mean_real"]["thres_with_sig_from_port_surrogates"],
                           "all_data_clean_pooling": thr["port_surrogates"]["thres"], "thesis_text": 9},
            "real_links_max_abs_diff": thr["real_links_max_abs_diff_vs_stored"],
            "surrogate_compare": thr["surrogate_compare"],
            "pytds_vs_port": pv}


def scores():
    return pd.read_csv(os.path.join(RES, "r2_scores.csv"))


def stage_matrices(sc):
    out = {}
    s = sc[sc.rule == "strict"]
    for v in ("published", "v2_calibrated"):
        out[v] = {}
        for st in STAGES:
            x = s[(s.variant == v) & (s.stage == st)].groupby("pair")[["real", "excess"]].mean()
            M = {m: np.full((10, 10), np.nan) for m in ("real", "excess")}
            for pair, row in x.iterrows():
                i, j = map(int, pair.split("-"))
                for m in M:
                    M[m][i, j] = M[m][j, i] = row[m]
            n = int(s[(s.variant == v) & (s.stage == st)].dropna(subset=["real"]).subject.nunique())
            out[v][st] = {"real": M["real"], "excess": M["excess"], "n_subjects": n}
    return out


def connectivity(sc):
    out = {}
    t = sc[sc.thesis_link == 1]
    for v in VARIANTS:
        out[v] = {}
        for rule in ("strict", "p80"):
            c = t[(t.variant == v) & (t.rule == rule)].groupby(["subject", "stage"]).excess.mean().unstack("stage")
            c = c.reindex(range(35))
            out[v][rule] = {st: c[st].tolist() for st in STAGES}
    return out


def effects():
    e = pd.read_csv(os.path.join(RES, "r2_effects.csv"))
    g = pd.read_csv(os.path.join(RES, "r2_gain.csv"))
    return {"effects": e[e.measure == "excess"].to_dict("records"),
            "gain": g[g.measure == "excess"].to_dict("records"),
            "summary": json.load(open(os.path.join(RES, "r2_summary.json")))["H5_decision_primary(excess,strict,v2_calibrated)"]}


def network():
    n = pd.read_csv(os.path.join(RES, "r2_network.csv"))
    agg = n.groupby(["variant", "stage"]).agg(mean_auc=("auc", "mean"), n_sig=("significant", "sum")).reset_index()
    cells = {v: {st: n[(n.variant == v) & (n.stage == st)].set_index("link").auc.to_dict() for st in STAGES}
             for v in n.variant.unique()}
    return {"by_stage": agg.to_dict("records"), "cells": cells,
            "H6": json.load(open(os.path.join(RES, "r2_summary.json")))["H6"]}


def mechanism():
    p = pd.read_csv(os.path.join(RES, "r2_params.csv"))
    summ = json.load(open(os.path.join(RES, "r2_summary.json")))
    raw = p[p.variant == "published"][[f"raw_bartlett_{s}" for s in SIGNALS]]
    med = raw.median()
    indep = {L: {s: float(L / (med[f"raw_bartlett_{s}"] / 30.0)) for s in SIGNALS} for L in (60, 120, 300, 600)}
    params = {v: {k: p[p.variant == v][k].describe()[["min", "50%", "max"]].tolist()
                  for k in ("window", "tolerance", "max_lag")} for v in VARIANTS}
    return {"median_raw_bartlett_window": {s: float(med[f"raw_bartlett_{s}"]) for s in SIGNALS},
            "raw_bartlett_window_quartiles": {s: raw[f"raw_bartlett_{s}"].quantile([.25, .5, .75]).tolist() for s in SIGNALS},
            "independent_samples_per_window": indep, "params_min_median_max": params,
            "n_defined": summ["n_defined"]}


def apnea(msc):
    r3 = json.load(open(os.path.join(RES, "r3_apnea_stored.json")))
    with h5py.File(f"{msc}/02 sleep apnea analysis/TDS imp/ranks_sick&healthy.mat") as f:
        sick, healthy = 100 * f["sick_ranks"][()].T, 100 * f["healthy_ranks"][()].T
    st = ["REM", "DS", "LS", "awake"]
    return {"stage_order": st, "threshold": 7.8,
            "ranks": {s: {"sick": sick[:, k], "healthy": healthy[:, k]} for k, s in enumerate(st)},
            "tests": r3["sick_vs_healthy"], "between_stages": r3["between_stages"],
            "full_night": r3["full_night"]}


def garmin():
    out = {}
    for part in ("main", "heldout"):
        s = pd.read_csv(os.path.join(RES, "garmin", f"summary_{part}.csv"))
        s["n_detect"] = (s.speed_hr_detect * s.n_runs).round().astype(int)
        s["n_xrun_fp"] = (s.xrun_fpr * s.n_xrun).round().astype(int)
        out[part] = s.to_dict("records")
    out["expectations"] = json.load(open(os.path.join(RES, "garmin", "expectations_check.json")))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--msc", required=True)
    ap.add_argument("--agg-out", required=True)
    a = ap.parse_args()
    sc = scores()
    out = {"signals": SIGNALS, "stages": STAGES, "replication": replication(a.msc),
           "stage_matrices": stage_matrices(sc), "connectivity": connectivity(sc), "effects": effects(),
           "network": network(), "mechanism": mechanism(), "apnea": apnea(a.msc), "garmin": garmin()}
    os.makedirs(os.path.dirname(os.path.abspath(a.agg_out)), exist_ok=True)
    json.dump(r(out), open(a.agg_out, "w"), separators=(",", ":"))
    print(a.agg_out, os.path.getsize(a.agg_out) // 1024, "KB")


if __name__ == "__main__":
    main()
