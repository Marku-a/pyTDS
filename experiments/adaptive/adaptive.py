"""
System-adaptive TDS parameters (EXPERIMENTAL — not part of the pyTDS API).

Idea: instead of fixed sample counts (window=60, step=30, max_lag=30,
tolerance=1), derive ONE parameter set for a whole system of signals from
the signals' own time scales:

1. Each signal's decorrelation time tau_c (first lag where its
   autocorrelation drops below 1/e).
2. System time scale = a quantile of the tau_c values (default: the
   slowest signal, quantile=1.0), so every pair has enough independent
   samples per window.
3. window = cycles * system time scale; step = window // 2;
   max_lag = (window - 1) // 2, optionally capped at half the shortest
   detected oscillation period (beyond that, xcorr peaks alias).
4. tolerance: either kept fixed, or calibrated on a coupling-free null.
   The null circularly shifts one signal of each pair, which destroys the
   coupling but keeps each signal's own dynamics. The rule then picks the
   largest tolerance whose chance TDS score stays <= alpha_stab.

All choices use only per-signal properties or a coupling-free null, so
they cannot be tuned towards the coupling being measured.
"""

from __future__ import annotations

import os
import sys
from dataclasses import replace
from itertools import combinations

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from pyTDS.core import stable_label, tds, tds_score, time_delay_interaction  # noqa: E402
from pyTDS.params import TDSParams  # noqa: E402


def _as_list(signals) -> list[np.ndarray]:
    if isinstance(signals, np.ndarray) and signals.ndim == 2:
        return [np.asarray(signals[:, i], dtype=float) for i in range(signals.shape[1])]
    return [np.asarray(s, dtype=float) for s in signals]


def acf(x: np.ndarray, max_lag: int | None = None) -> np.ndarray:
    """Normalised autocorrelation r[0..max_lag] via FFT (r[0] = 1)."""
    x = np.asarray(x, dtype=float)
    n = len(x)
    if max_lag is None:
        max_lag = n // 4
    x = x - x.mean()
    if not np.any(x):
        return np.zeros(max_lag + 1)
    nfft = 1 << int(np.ceil(np.log2(2 * n)))
    X = np.fft.rfft(x, nfft)
    r = np.fft.irfft(X * np.conj(X), nfft)[: max_lag + 1]
    return r / r[0]


def decorrelation_time(x: np.ndarray, max_lag: int | None = None) -> float:
    """First lag (linearly interpolated) where the ACF drops below 1/e."""
    r = acf(x, max_lag)
    thr = np.exp(-1.0)
    below = np.nonzero(r < thr)[0]
    if len(below) == 0:
        return float(len(r) - 1)
    k = int(below[0])
    if k == 0:
        return 0.0
    # interpolate between k-1 (>= thr) and k (< thr)
    return float(k - 1 + (r[k - 1] - thr) / (r[k - 1] - r[k]))


def dominant_period(x: np.ndarray, max_lag: int | None = None, min_peak: float = 0.3) -> int | None:
    """Lag of the highest ACF peak after its first zero crossing, if > min_peak."""
    r = acf(x, max_lag)
    neg = np.nonzero(r <= 0)[0]
    if len(neg) == 0:
        return None
    z = int(neg[0])
    tail = r[z:]
    if len(tail) < 3:
        return None
    k = int(np.argmax(tail))
    if tail[k] < min_peak or k in (0, len(tail) - 1):
        return None
    return z + k


def bartlett_factor(x: np.ndarray, max_lag: int | None = None) -> float:
    """
    Variance inflation of a correlation estimate due to autocorrelation:
    B = 1 + 2 * sum_k r(k)^2, summed until the ACF first crosses zero.
    A window of L samples holds about L / B independent samples
    (Bartlett 1946). Unlike the 1/e decorrelation time, the sum sees the
    whole ACF, so slow structure (e.g. bursts riding on spikes) counts.
    """
    r = acf(x, max_lag)
    neg = np.nonzero(r <= 0)[0]
    K = int(neg[0]) if len(neg) else len(r)
    return float(1.0 + 2.0 * np.sum(r[1:K] ** 2))


def adaptive_params(
    signals,
    cycles: float = 10.0,
    window_rule: str = "decorr",
    n_eff: float = 30.0,
    quantile: float = 1.0,
    min_window: int = 20,
    max_window_frac: float = 0.1,
    cap_lag_by_period: bool = True,
    tolerance: str = "fixed",
    alpha_stab: float = 0.05,
    n_shifts: int = 20,
    max_pairs: int = 50,
    seed: int = 0,
    base: TDSParams | None = None,
) -> tuple[TDSParams, dict]:
    """
    Derive one TDSParams for a whole system of signals.

    Parameters
    ----------
    signals : T x N array or list of N equal-length 1-D arrays.
    cycles : window length in units of the system decorrelation time
        (window_rule="decorr", the v1 rule).
    window_rule : "decorr" (v1: window = cycles * decorrelation time) or
        "bartlett" (v2: window = n_eff * Bartlett factor, i.e. the window
        holds about n_eff independent samples for every signal).
    n_eff : target independent samples per window (window_rule="bartlett").
    quantile : which per-signal decorrelation time defines the system
        scale (1.0 = slowest signal, 0.5 = median signal).
    min_window, max_window_frac : window is clipped to
        [min_window, max(min_window, max_window_frac * T)].
    cap_lag_by_period : cap max_lag at half the shortest detected period.
    tolerance : "fixed" keeps base.tolerance; "calibrated" chooses it on a
        circular-shift null (see module docstring).
    alpha_stab : target chance TDS score (fraction) for calibration.
    n_shifts, max_pairs, seed : null-calibration effort / randomness.
    base : params to start from (stability rule, anchors). Default TDSParams().

    Returns
    -------
    params : TDSParams
    info : dict with the intermediate quantities (for reporting).
    """
    base = base or TDSParams()
    sigs = _as_list(signals)
    T = min(len(s) for s in sigs)
    sigs = [s[:T] for s in sigs]

    tau_c = np.array([decorrelation_time(s) for s in sigs])
    periods = [dominant_period(s) for s in sigs]
    tau_sys = float(np.quantile(tau_c, quantile))

    hi = max(min_window, int(max_window_frac * T))
    bart = np.array([bartlett_factor(s) for s in sigs])
    if window_rule == "decorr":
        L_raw = round(cycles * max(tau_sys, 1.0))
    elif window_rule == "bartlett":
        L_raw = round(n_eff * float(np.quantile(bart, quantile)))
    else:
        raise ValueError("window_rule must be 'decorr' or 'bartlett'")
    L = int(np.clip(L_raw, min_window, hi))
    step = max(1, L // 2)
    max_lag = (L - 1) // 2
    detected = [p for p in periods if p is not None]
    if cap_lag_by_period and detected:
        max_lag = max(1, min(max_lag, min(detected) // 2))

    params = replace(base, window=L, overlap=step, max_lag=max_lag)
    info = {
        "T": T,
        "tau_c": tau_c.tolist(),
        "periods": periods,
        "tau_sys": tau_sys,
        "bartlett": bart.tolist(),
        "window_rule": window_rule,
        "window": L,
        "step": step,
        "max_lag": max_lag,
        "window_clipped": L != L_raw,
    }

    if tolerance == "fixed":
        info["tolerance"] = params.tolerance
        return params, info
    if tolerance != "calibrated":
        raise ValueError("tolerance must be 'fixed' or 'calibrated'")

    rng = np.random.default_rng(seed)
    pairs = list(combinations(range(len(sigs)), 2))
    if len(pairs) > max_pairs:
        idx = rng.choice(len(pairs), max_pairs, replace=False)
        pairs = [pairs[i] for i in idx]
    lo_shift, hi_shift = max(L, T // 10), T - max(L, T // 10)
    null_taus = []
    for i, j in pairs:
        for _ in range(n_shifts):
            sh = int(rng.integers(lo_shift, max(lo_shift + 1, hi_shift)))
            tau, _, _ = time_delay_interaction(sigs[i], np.roll(sigs[j], sh), params)
            null_taus.append(tau)

    candidates = list(range(0, max(1, max_lag // 4) + 1))
    null_by_tol = {}
    chosen = 0
    for tol in candidates:
        p_tol = replace(params, tolerance=tol)
        s = float(np.mean([tds_score(stable_label(t, p_tol)) for t in null_taus]))
        null_by_tol[tol] = s
        if s <= 100.0 * alpha_stab:
            chosen = tol
        else:
            break  # null score is monotone in tolerance
    params = replace(params, tolerance=chosen)
    info["tolerance"] = chosen
    info["null_score_by_tolerance"] = null_by_tol
    return params, info


def tds_adaptive(s1, s2, system=None, **kwargs) -> tuple[dict, TDSParams, dict]:
    """
    Run TDS on (s1, s2) with system-adaptive params.

    If `system` (all signals of the system) is given, params come from it;
    otherwise from the pair alone (per-pair adaptation).
    """
    params, info = adaptive_params(system if system is not None else [s1, s2], **kwargs)
    return tds(s1, s2, params), params, info
