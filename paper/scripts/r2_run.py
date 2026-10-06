"""R2: thesis-TDS vs aTDS on D1 (35 subjects). See paper/ANALYSIS_PLAN.md, section R2.

Per subject x variant x pair: stage-wise TDS score (strict and p80 rules), circular-shift chance
(5 shifts) and, for the 35 thesis links with the strict rule and the 5 main variants, 50 null draws.
Usage: python paper/scripts/r2_run.py [--subjects 0 1] [--workers 3] [--outdir DIR]
Inputs are read from the npz in the scratchpad; no signal data is written to the repo.
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from itertools import combinations
from multiprocessing import Pool

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "experiments", "adaptive_groundtruth"))
sys.path.insert(0, ROOT)

from adaptive import adaptive_params, bartlett_factor  # noqa: E402
from common import (calibrate_tolerance, published_params, stable_seed,  # noqa: E402,F401
                    variant_params, window_labels)

DATA = "/tmp/claude-0/-home-user/b218c801-5d78-507e-9550-a8f2a0e7c418/scratchpad/ncom35.npz"
STAGES = ("LS", "DS", "REM", "awake")
MAIN = ("published", "win_only", "tol_only", "v2_calibrated", "v2_cal_lag3")
VARIANTS = MAIN + ("v2_cap120", "v2_cap600")
SIGNALS = ("HR", "Resp", "Chin", "Leg", "Eye", "delta", "theta", "alpha", "sigma", "beta")
PAIRS = list(combinations(range(10), 2))
THESIS = [(i, j) for (i, j) in PAIRS if not (i >= 5 and j >= 5)]
assert len(PAIRS) == 45 and len(THESIS) == 35
N_CHANCE, N_NULL = 5, 50


def window_stage(centres, st, W):
    """Stage index (0..3 per STAGES) of each window if >= 80 % of its samples carry it, else -1.
    Same logic as sleep_edf.py:window_stage, for per-sample integer stages (-1 = unlabelled)."""
    out = np.full(len(centres), -1)
    for i, c in enumerate(centres):
        s = int(c) - W // 2
        seg = st[max(0, s):s + W]
        if len(seg) < W:
            continue
        for k in range(4):
            if np.mean(seg == k) >= 0.8:
                out[i] = k
                break
    return out


def stage_masks(ws):
    """(strict, p80) boolean masks [4, n_windows]."""
    n = len(ws)
    p80 = np.stack([ws == k for k in range(4)])
    strict = np.zeros_like(p80)
    for k in range(4):
        m = p80[k]
        if n >= 5:
            strict[k, 2:n - 2] = m[:-4] & m[1:-3] & m[2:-2] & m[3:-1] & m[4:]
    return strict, p80


def score(lbl, mask):
    """[4] stage scores (NaN if < 5 counted windows)."""
    return np.array([100.0 * lbl[m].mean() if m.sum() >= 5 else np.nan for m in mask])


def run_subject(args):
    idx, = args
    d = np.load(DATA, allow_pickle=True)
    name = str(d["names"][idx])
    X = np.asarray(d["X_" + name], float)
    S = np.asarray(d["S_" + name]).astype(str)
    T = min(len(X), len(S))
    X, S = X[:T], S[:T]
    code = {s: k for k, s in enumerate(STAGES)}
    st = np.array([code.get(s, -1) for s in S])
    sig = [X[:, k] for k in range(10)]
    seed = stable_seed("ncom35", idx)

    vp = variant_params(sig, fs=1.0, max_window_s=300, seed=seed)
    pars = {v: vp[v] for v in MAIN}
    for cap in (120, 600):
        p, info = adaptive_params(sig, window_rule="bartlett", tolerance="calibrated",
                                  base=published_params(1.0), max_window=cap, n_shifts=10,
                                  max_pairs=10, seed=seed)
        pars[f"v2_cap{cap}"] = (p, {"window_clipped": info["window_clipped"],
                                    "calibration_failed": info["calibration_failed"]})
    raw_bart = [round(30 * bartlett_factor(x)) for x in sig]

    prow = []
    for v in VARIANTS:
        p, info = pars[v]
        r = {"subject": idx, "variant": v, "window": p.window, "step": p.overlap,
             "max_lag": p.max_lag, "tolerance": p.tolerance,
             "calibration_failed": info.get("calibration_failed", np.nan),
             "window_clipped": info.get("window_clipped", np.nan)}
        for k, s in enumerate(SIGNALS):
            r["raw_bartlett_" + s] = raw_bart[k]
        prow.append(r)

    rows = []
    null = np.full((len(MAIN), len(THESIS), 4, N_NULL), np.nan)
    for vi, v in enumerate(VARIANTS):
        p = pars[v][0]
        W = p.window
        lo = max(2 * W, 5 * p.max_lag)
        if 2 * lo >= T:
            lo = T // 3
        for (i, j) in PAIRS:
            pair = f"{i}-{j}"
            is_thesis = (i, j) in THESIS
            c, lbl, _ = window_labels(sig[i], sig[j], p)
            lbl = np.asarray(lbl, float)
            ws = window_stage(c, st, W)
            strict, p80 = stage_masks(ws)
            rng = np.random.default_rng(stable_seed("ncom35", idx, pair, v, "chance"))
            ch = []
            for _ in range(N_CHANCE):
                sh = int(rng.integers(lo, T - lo + 1))
                _, l2, _ = window_labels(sig[i], np.roll(sig[j], sh), p)
                ch.append((np.asarray(l2, float)))
            for rule, mask in (("strict", strict), ("p80", p80)):
                real = score(lbl, mask)
                chs = np.array([score(l2, mask) for l2 in ch])
                with np.errstate(all="ignore"):
                    import warnings
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        chance = np.nanmean(chs, axis=0)
                for k in range(4):
                    rows.append((idx, v, pair, int(is_thesis), STAGES[k], rule, int(mask[k].sum()),
                                 real[k], chance[k], real[k] - chance[k]))
            if is_thesis and v in MAIN:
                rng = np.random.default_rng(stable_seed("ncom35", idx, pair, v, "null"))
                li = THESIS.index((i, j))
                for n in range(N_NULL):
                    sh = int(rng.integers(lo, T - lo + 1))
                    _, l2, _ = window_labels(sig[i], np.roll(sig[j], sh), p)
                    null[MAIN.index(v), li, :, n] = score(np.asarray(l2, float), strict)
    return idx, rows, prow, null


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subjects", type=int, nargs="*", default=None)
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--outdir", default=os.path.join(ROOT, "paper", "results"))
    a = ap.parse_args()
    subs = a.subjects if a.subjects else list(range(35))
    os.makedirs(a.outdir, exist_ok=True)
    t0 = time.time()
    rows, prow, nulls = [], [], {}
    with Pool(a.workers) as pool:
        for idx, r, pr, nl in pool.imap_unordered(run_subject, [(s,) for s in subs]):
            rows += r
            prow += pr
            nulls[idx] = nl
            print(f"subject {idx} done, {time.time() - t0:.0f}s elapsed", flush=True)
    cols = ["subject", "variant", "pair", "thesis_link", "stage", "rule", "n_windows", "real", "chance", "excess"]
    pd.DataFrame(rows, columns=cols).sort_values(["subject", "variant", "pair", "rule", "stage"]) \
        .to_csv(os.path.join(a.outdir, "r2_scores.csv"), index=False)
    pd.DataFrame(prow).sort_values(["subject", "variant"]).to_csv(os.path.join(a.outdir, "r2_params.csv"), index=False)
    order = sorted(nulls)
    arrs = {v: np.stack([nulls[s][vi] for s in order]) for vi, v in enumerate(MAIN)}
    np.savez(os.path.join(a.outdir, "r2_null.npz"), subjects=np.array(order),
             link_order=np.array([f"{i}-{j}" for i, j in THESIS]), stage_order=np.array(STAGES),
             signals=np.array(SIGNALS),
             layout=np.array("per variant: [subject, link, stage, draw]; strict rule, 50 circular-shift null draws"),
             **arrs)
    print(f"total {time.time() - t0:.0f}s for {len(subs)} subjects", flush=True)


if __name__ == "__main__":
    main()
