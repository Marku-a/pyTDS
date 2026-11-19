"""aTDS rule v2.1 (Algorithm 1) shared by the real-data experiments. Uses only this repo's pyTDS."""
from __future__ import annotations

import os
import sys
from collections import Counter
from dataclasses import replace
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "atds_brief", "experiments", "atds_code"))
sys.path.insert(0, ROOT)

from adaptive import (  # noqa: E402,F401
    acf, adaptive_params, bartlett_factor, dominant_period)
from pyTDS.core import stable_label, tds_score, time_delay_interaction  # noqa: E402
from pyTDS.params import TDSParams  # noqa: E402

PUBLISHED = TDSParams(window=60, overlap=30, max_lag=30, tolerance=1)
WORKERS = 4


def window_params(signals, d_max=None, cap=None) -> TDSParams:
    """Bartlett window rule (30 independent samples), min window 61, lag keeps the rhythm cap; tolerance fixed."""
    min_w = max(61, 2 * d_max + 1 if d_max else 61)
    p, _ = adaptive_params(signals, window_rule="bartlett", tolerance="fixed", min_window=min_w, max_window=cap)
    return p


def calib(p, pairs, n_shift=20, z=2.0, seed=0):
    """Largest tol with mean + z*SE of the chance score <= 5 %, rotations pooled over all training pairs.

    Copied from atds_brief/experiments/exp4.py; at most 50 pairs are used (random, seeded).
    """
    pairs = list(pairs)
    if len(pairs) > 50:
        idx = np.random.default_rng(seed).choice(len(pairs), 50, replace=False)
        pairs = [pairs[i] for i in sorted(idx)]
    rng = np.random.default_rng(seed); taus = []
    for x, y in pairs:
        n = len(x); lo = max(p.window, n // 10)
        for _ in range(n_shift): taus.append(time_delay_interaction(x, np.roll(y, int(rng.integers(lo, n - lo))), p)[0])
    tol, fail = 0, False
    for t in range(0, max(1, p.max_lag // 4) + 1):
        s = np.array([tds_score(stable_label(tt, replace(p, tolerance=t))) for tt in taus])
        ub = s.mean() + z * s.std(ddof=1) / np.sqrt(len(s))
        if ub <= 5.0: tol = t
        else:
            if t == 0: fail = True
            break
    return replace(p, tolerance=tol), fail


def delay_of(tau, labels):
    """Most frequent tau0 among stable windows (None if fewer than 3 stable windows)."""
    tau = np.asarray(tau); st = tau[np.asarray(labels).astype(bool)]
    if len(st) < 3:
        return None
    return float(Counter(st.tolist()).most_common(1)[0][0])


# ---- rotation null for a single pair (shared by Part B and Part F)
_CTX: dict = {}


def _null_worker(shift):
    x, y, p = _CTX["x"], _CTX["y"], _CTX["p"]
    tau = time_delay_interaction(x, np.roll(y, int(shift)), p)[0]
    return tds_score(stable_label(tau, p))


def rotation_null(x, y, p, n_null=200, seed=0, procs=WORKERS):
    """TDS scores of n_null circular rotations of y (shift in [T//10, T - T//10]); default anchor."""
    T = len(x)
    shifts = np.random.default_rng(seed).integers(T // 10, T - T // 10 + 1, n_null)
    _CTX.update(x=x, y=y, p=p)
    with Pool(procs) as pool:
        return np.array(pool.map(_null_worker, shifts, chunksize=5))


def null_stats(real, null):
    return {"null_mean": float(null.mean()), "null95": float(np.percentile(null, 95)),
            "p": float((1 + np.sum(null >= real)) / (1 + len(null)))}


def interp_nan(x):
    """Linear interpolation of NaNs (edges held constant)."""
    x = np.asarray(x, float).copy()
    bad = ~np.isfinite(x)
    if bad.any() and not bad.all():
        idx = np.arange(len(x))
        x[bad] = np.interp(idx[bad], idx[~bad], x[~bad])
    return x
