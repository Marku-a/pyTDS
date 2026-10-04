"""
Shared core for the ground-truth study of adaptive TDS.

Every dataset script builds, per subject, a list of equal-length 1-D signals
sampled at `fs` Hz, then calls `variant_params` once per subject and
`pair_test` (or `window_labels`) per pair. Keeping the variants and the null
here guarantees every dataset is judged the same way.

Variants (an ablation of what v2-calibrated changes):

  default_samples  pyTDS defaults in raw samples (60 / 30 / 30 / +-1)
  published        the published TDS setting in SECONDS (Bashan 2012:
                   60 s window, 30 s step, +-30 s lag, +-1 s tolerance)
  win_only         v2 window rule; step and max lag follow the window;
                   tolerance stays the published +-1 s
  tol_only         published window / step / lag; tolerance calibrated on a
                   circular-shift null (5 % chance budget)
  v2_calibrated    the full adaptive rule (window + lag + calibrated tolerance)
  v2_cal_lag3      v2_calibrated with max lag capped at window / 3, tolerance
                   re-calibrated: is the gain from the wide lag search?
"""

from __future__ import annotations

import os
import sys
import zlib
from dataclasses import replace
from itertools import combinations

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.dirname(__file__))

from adaptive import adaptive_params  # noqa: E402
from pyTDS.core import stable_label, tds, tds_score, time_delay_interaction  # noqa: E402
from pyTDS.params import TDSParams  # noqa: E402

VARIANTS = ("default_samples", "published", "win_only", "tol_only", "v2_calibrated", "v2_cal_lag3")
ALPHA = 0.05


def stable_seed(*parts) -> int:
    """Seed that is identical across processes (Python's hash() is salted)."""
    return zlib.crc32("|".join(map(str, parts)).encode())


def published_params(fs: float) -> TDSParams:
    """Bashan et al. 2012 setting, in seconds, converted to samples at fs."""
    return TDSParams(window=max(8, round(60 * fs)), overlap=max(1, round(30 * fs)),
                     max_lag=max(1, round(30 * fs)), tolerance=max(1, round(1 * fs)))


def calibrate_tolerance(signals, params: TDSParams, alpha_stab: float = ALPHA, n_shifts: int = 10,
                        max_pairs: int = 10, seed: int = 0) -> tuple[TDSParams, dict]:
    """Largest tolerance whose chance TDS score stays <= alpha_stab, with params' window and lag.

    Same null and rule as adaptive.adaptive_params(tolerance="calibrated"), but the
    window and max lag are left as given (for the tol_only and lag3 ablations).
    """
    sigs = [np.asarray(s, float) for s in signals]
    T = min(len(s) for s in sigs)
    sigs = [s[:T] for s in sigs]
    L = params.window
    rng = np.random.default_rng(seed)
    pairs = list(combinations(range(len(sigs)), 2))
    if len(pairs) > max_pairs:
        pairs = [pairs[i] for i in rng.choice(len(pairs), max_pairs, replace=False)]
    lo, hi = max(L, T // 10), T - max(L, T // 10)
    null_taus = []
    for i, j in pairs:
        for _ in range(n_shifts):
            sh = int(rng.integers(lo, max(lo + 1, hi)))
            tau, _, _ = time_delay_interaction(sigs[i], np.roll(sigs[j], sh), params)
            null_taus.append(tau)
    eff_lag = min(params.max_lag, (L - 1) // 2)
    chosen, by_tol = 0, {}
    for tol in range(0, max(1, eff_lag // 4) + 1):
        p = replace(params, tolerance=tol)
        s = float(np.mean([tds_score(stable_label(t, p)) for t in null_taus])) if null_taus else 0.0
        by_tol[tol] = s
        if s <= 100.0 * alpha_stab:
            chosen = tol
        else:
            break  # chance score is monotone in tolerance
    return replace(params, tolerance=chosen), {"calibration_failed": by_tol.get(0, 0.0) > 100.0 * alpha_stab}


def variant_params(signals, fs: float, max_window_s: float, seed: int = 0) -> dict[str, tuple[TDSParams, dict]]:
    """variant -> (params, info) for one subject's system of signals.

    `max_window_s` is the user's time-resolution cap (the longest window accepted),
    fixed per dataset BEFORE looking at results.
    """
    pub = published_params(fs)
    cap = max(20, round(max_window_s * fs))
    out = {"default_samples": (TDSParams(), {}), "published": (pub, {})}

    p, info = adaptive_params(signals, window_rule="bartlett", tolerance="fixed", base=pub, max_window=cap)
    out["win_only"] = (p, {"window_clipped": info["window_clipped"]})

    p, info = calibrate_tolerance(signals, pub, seed=seed)
    out["tol_only"] = (p, info)

    p, info = adaptive_params(signals, window_rule="bartlett", tolerance="calibrated", base=pub,
                              max_window=cap, n_shifts=10, max_pairs=10, seed=seed)
    out["v2_calibrated"] = (p, {"window_clipped": info["window_clipped"],
                                "calibration_failed": info["calibration_failed"]})

    lag3 = replace(p, max_lag=max(1, p.window // 3))
    p3, info3 = calibrate_tolerance(signals, lag3, seed=seed)
    out["v2_cal_lag3"] = (p3, info3)
    return out


def params_row(params: TDSParams, fs: float) -> dict:
    """Parameters in samples and seconds, for reporting."""
    return {"window": params.window, "step": params.overlap, "max_lag": params.max_lag,
            "tolerance": params.tolerance, "window_s": params.window / fs,
            "max_lag_s": params.max_lag / fs, "tolerance_s": params.tolerance / fs}


def _summary(tau: np.ndarray, params: TDSParams) -> tuple[float, float]:
    """(score %, median stable tau in samples; NaN if none). pyTDS sign: tau < 0 means x leads y."""
    lbl = stable_label(tau, params)
    st = tau[lbl == 1]
    return tds_score(lbl), (float(np.median(st)) if len(st) else float("nan"))


def pair_test(x, y, params: TDSParams, n_null: int = 100, seed: int = 0) -> dict:
    """TDS score of (x, y) against a circular-shift null of y.

    Returns score (%), chance95 (95th pct of null scores), p (one-sided,
    (1 + #null >= score) / (1 + n_null)), lag_samples (median stable tau,
    pyTDS sign: negative = x leads y), n_windows.
    """
    x, y = np.asarray(x, float), np.asarray(y, float)
    n = len(x)
    tau, _, _ = time_delay_interaction(x, y, params)
    score, lag = _summary(tau, params)
    lo = max(2 * params.window, 5 * params.max_lag)
    if 2 * lo >= n:
        lo = n // 3
    rng = np.random.default_rng(seed)
    null = np.empty(n_null)
    for k in range(n_null):
        t0, _, _ = time_delay_interaction(x, np.roll(y, int(rng.integers(lo, n - lo + 1))), params)
        null[k] = _summary(t0, params)[0]
    return {"score": score, "chance95": float(np.percentile(null, 95)) if n_null else float("nan"),
            "p": float((1 + np.sum(null >= score)) / (1 + n_null)) if n_null else float("nan"),
            "lag_samples": lag, "n_windows": int(len(tau))}


def window_labels(x, y, params: TDSParams) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per-window (centre sample index, stable 0/1, tau samples) for state-wise scoring."""
    r = tds(np.asarray(x, float), np.asarray(y, float), params)
    return np.asarray(r["t_vec"]), np.asarray(r["stbl_lbl"]), np.asarray(r["tau"])


# ---------------------------------------------------------------- statistics

def auc(pos, neg) -> float:
    """Mann-Whitney AUC, P(pos > neg), ties 0.5."""
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    pos, neg = pos[np.isfinite(pos)], neg[np.isfinite(neg)]
    if not len(pos) or not len(neg):
        return float("nan")
    d = pos[:, None] - neg[None, :]
    return float((np.sum(d > 0) + 0.5 * np.sum(d == 0)) / d.size)


def paired_effect(a, b) -> dict:
    """Paired comparison a vs b (e.g. light vs deep sleep per subject).

    Returns n, mean difference, Cohen's dz, share of subjects with a > b, and an
    exact two-sided sign-test p (ties dropped).
    """
    from math import comb

    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    d = a[ok] - b[ok]
    n = len(d)
    if n == 0:
        return {"n": 0, "mean_diff": float("nan"), "dz": float("nan"), "frac_pos": float("nan"), "sign_p": float("nan")}
    pos, neg = int(np.sum(d > 0)), int(np.sum(d < 0))
    m = pos + neg
    k = min(pos, neg)
    sign_p = min(1.0, 2 * sum(comb(m, i) for i in range(k + 1)) / 2 ** m) if m else 1.0
    sd = float(np.std(d, ddof=1)) if n > 1 else float("nan")
    return {"n": n, "mean_diff": float(np.mean(d)), "dz": float(np.mean(d) / sd) if sd and sd > 0 else float("nan"),
            "frac_pos": pos / n, "sign_p": sign_p}


def group_effect(a, b) -> dict:
    """Unpaired comparison a vs b (e.g. young vs elderly): n's, means, Cohen's d, AUC P(a > b)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return {"n_a": len(a), "n_b": len(b), "d": float("nan"), "auc": float("nan")}
    sp = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    return {"n_a": len(a), "n_b": len(b), "mean_a": float(a.mean()), "mean_b": float(b.mean()),
            "d": float((a.mean() - b.mean()) / sp) if sp > 0 else float("nan"), "auc": auc(a, b)}
