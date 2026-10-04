"""Shared helpers for the adaptive-TDS synthetic experiments."""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.dirname(__file__))

from adaptive import adaptive_params  # noqa: E402
from pyTDS.params import TDSParams  # noqa: E402


def auc(pos, neg) -> float:
    """Mann-Whitney AUC, P(pos > neg) with ties counted 0.5."""
    pos = np.asarray(pos, dtype=float)
    neg = np.asarray(neg, dtype=float)
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    d = pos[:, None] - neg[None, :]
    return float((np.sum(d > 0) + 0.5 * np.sum(d == 0)) / d.size)


def variant_params(name: str, system):
    """Return (params, info) for a variant, given the pooled system signals."""
    if name == "fixed_default":
        p = TDSParams()
        return p, {"window": p.window, "step": p.overlap, "max_lag": p.max_lag,
                   "tolerance": p.tolerance}
    if name == "adaptive_fixedtol":
        return adaptive_params(system, tolerance="fixed")
    if name == "adaptive_calibrated":
        return adaptive_params(system, tolerance="calibrated", n_shifts=5, max_pairs=20)
    raise ValueError(name)


def oracle_grid() -> list[TDSParams]:
    grid = []
    for L in (30, 60, 120, 240, 480, 960, 1920):
        ml = (L - 1) // 2
        for tol in (0, 1, 2, 4, 8, 16):
            if tol <= ml:
                grid.append(TDSParams(window=L, overlap=L // 2, max_lag=ml, tolerance=tol))
    return grid


def delay_metrics(results, true_tau):
    """
    results: per coupled trial, either a tds() dict or a precomputed median
    (float) of the stable taus. Returns (hit_rate, median_rel_error).
    """
    meds = []
    for r in results:
        if isinstance(r, dict):
            st = np.asarray(r["stable_taus"])
            meds.append(float(np.median(st)) if len(st) else float("nan"))
        else:
            meds.append(float(r))
    meds = np.array(meds)
    tol = max(1.0, 0.1 * abs(true_tau))
    ok = np.abs(meds - true_tau) <= tol  # NaN -> False
    hit = float(np.mean(ok)) if len(meds) else float("nan")
    valid = ~np.isnan(meds)
    rel = (float(np.median(np.abs(meds[valid] - true_tau) / abs(true_tau)))
           if valid.any() else float("nan"))
    return hit, rel
