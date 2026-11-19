"""Shared helpers for the aTDS run-3 experiments (pyTDS only for TDS; baselines in numpy).

Needs a checkout of branch claude/practical-mccarthy-n2mec6 (see atds_brief/README.md); its path is read from the
environment variable ATDS_CODE. exp4.py / exp3.py / tests.py are imported from atds_brief/experiments; they read their
checkout path from sys.argv, so argv is set temporarily during the import.
"""
from __future__ import annotations

import csv
import json
import os
import sys
import time
from dataclasses import replace

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
BD = os.path.join(ROOT, "atds_brief", "experiments")
CODE = os.environ.get("ATDS_CODE")
if not CODE:
    raise SystemExit("set ATDS_CODE to a checkout of origin/claude/practical-mccarthy-n2mec6")
WORKERS = 4
RESULTS = os.path.join(HERE, "results")

sys.path.insert(0, os.path.join(ROOT, "experiments", "atds_real"))
import atds21 as A  # noqa: E402,F401  (puts pyTDS and atds_code on sys.path)

_argv = sys.argv
sys.argv = [_argv[0], CODE, BD]
sys.path.insert(0, BD)
import exp4 as E4  # noqa: E402
sys.argv = _argv
import adaptive as AD  # noqa: E402
from exp4 import (FIXED, MINW, N, ar1, ar2, boot_auc, boot_diff90, coupled_pair, found, osc_pair,  # noqa: E402,F401
                  pair, run_methods, score, summarize, wilson, WGRID, TGRID, auc, mode_tau)
from pyTDS.core import stable_label, tds_score, time_delay_interaction  # noqa: E402,F401
from pyTDS.params import TDSParams  # noqa: E402,F401


# ------------------------------------------------------------------ parametrised aTDS (all constants exposed)
def calib_var(p, train, n_shift=20, z=2.0, budget=5.0, seed=0):
    """exp4.calib with the chance-score budget exposed (budget=5, z=2, n_shift=20 reproduces exp4.calib)."""
    rng = np.random.default_rng(seed); taus = []
    for x, y in train:
        n = len(x); lo = max(p.window, n // 10)
        for _ in range(n_shift):
            taus.append(time_delay_interaction(x, np.roll(y, int(rng.integers(lo, n - lo))), p)[0])
    tol, fail = 0, False
    for t in range(0, max(1, p.max_lag // 4) + 1):
        s = np.array([tds_score(stable_label(tt, replace(p, tolerance=t))) for tt in taus])
        ub = s.mean() + z * s.std(ddof=1) / np.sqrt(len(s))
        if ub <= budget:
            tol = t
        else:
            if t == 0:
                fail = True
            break
    return replace(p, tolerance=tol), fail


def atds_var(train, n_eff=30.0, budget=5.0, z=2.0, floor=61, n_shift=20, cap=None, min_peak=0.3):
    """aTDS v2.1 with every design constant exposed (defaults = exp4.atds2)."""
    orig = AD.dominant_period
    if min_peak != 0.3:
        AD.dominant_period = lambda x, max_lag=None, min_peak_=min_peak: orig(x, max_lag, min_peak_)
    try:
        x, y = train[0]
        p0, _ = AD.adaptive_params([x, y], window_rule="bartlett", n_eff=n_eff, tolerance="fixed",
                                   min_window=floor, max_window=cap)
    finally:
        AD.dominant_period = orig
    return calib_var(p0, train, n_shift=n_shift, z=z, budget=budget)


def pdict(p):
    return dict(window=p.window, step=p.overlap, max_lag=p.max_lag, tol=p.tolerance)


# ------------------------------------------------------------------ non-TDS baselines (numpy)
LAGMAX = 1000


def _zs(x):
    x = np.asarray(x, float)
    return (x - x.mean()) / x.std()


def _lag_slice(n, m=LAGMAX):
    return np.arange(-m, m + 1)


def xcorr_full(x, y, m=LAGMAX):
    """Normalised linear cross-correlation of z-scored x, y, lags -m..m (positive lag: y lags x).

    Returns (score = max |C|, delay in the TDS convention = -argmax lag, C)."""
    x, y = _zs(x), _zs(y); n = len(x); nf = 1 << int(np.ceil(np.log2(2 * n)))
    c = np.fft.irfft(np.conj(np.fft.rfft(x, nf)) * np.fft.rfft(y, nf), nf)
    lags = _lag_slice(n, m); C = c[lags % nf] / (n - np.abs(lags))
    k = int(np.argmax(np.abs(C)))
    return float(abs(C[k])), -int(lags[k]), C


def gccphat_full(x, y, m=LAGMAX, eps=1e-12):
    """GCC-PHAT of z-scored x, y, lags -m..m; score = max |R|, delay = -argmax lag (TDS convention)."""
    x, y = _zs(x), _zs(y); n = len(x); nf = 1 << int(np.ceil(np.log2(2 * n)))
    G = np.conj(np.fft.rfft(x, nf)) * np.fft.rfft(y, nf)
    r = np.fft.irfft(G / (np.abs(G) + eps), nf)
    lags = _lag_slice(n, m); R = r[lags % nf]
    k = int(np.argmax(np.abs(R)))
    return float(abs(R[k])), -int(lags[k]), R


def hit_global(delay, d):
    """Delay found for a global method (same tolerance as the TDS hit rule, tau = -d)."""
    return bool(abs(delay + d) <= max(1, 0.1 * abs(d)))


def first_below_e(x):
    """First lag at which the autocorrelation falls below 1/e."""
    r = AD.acf(x, len(x) // 4); k = np.nonzero(r < np.exp(-1.0))[0]
    return float(k[0]) if len(k) else float(len(r) - 1)


# ------------------------------------------------------------------ stats / io
def auc_row(pos, neg, B=1000, seed=0):
    return dict(auc=round(auc(pos, neg), 3), auc_ci=boot_auc(pos, neg, B, seed))


def hit_row(hits):
    k, n = int(sum(hits)), len(hits)
    return dict(hit=round(k / n, 3), hit_ci=wilson(k, n), n=n)


def write_json(name, obj):
    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, name), "w") as f:
        json.dump(obj, f, indent=1, default=str)


def write_csv(name, rows):
    os.makedirs(RESULTS, exist_ok=True)
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(os.path.join(RESULTS, name), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)


def flat(d, prefix=""):
    """Flatten nested dict; lists stay as strings."""
    out = {}
    for k, v in d.items():
        if isinstance(v, dict):
            out.update(flat(v, f"{prefix}{k}_"))
        else:
            out[f"{prefix}{k}"] = json.dumps(v) if isinstance(v, (list, tuple)) else v
    return out


def smoke_n(default):
    """--smoke on the command line shrinks the sample (output goes to results_smoke)."""
    return (4, True) if "--smoke" in sys.argv else (default, False)


def set_results_for_smoke(smoke):
    global RESULTS
    if smoke:
        RESULTS = os.path.join(HERE, "results_smoke")
