"""
Core TDS algorithm — 3 steps:
  1. time_delay_interaction() : sliding-window PBC cross-correlation → τ₀(t)
  2. stable_label()           : label each τ₀ as stable or unstable
  3. tds_score()              : fraction of stable points × 100
  4. tds()                    : full pipeline combining all three
"""

import numpy as np
from .params import TDSParams, VALID_STABILITY_ANCHORS
from .utils import zscore, pbc_xcorr


def time_delay_interaction(
    s1: np.ndarray,
    s2: np.ndarray,
    params: TDSParams = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Step 1: compute the time delay τ₀ between two signals over time using
    a sliding window with periodic-boundary-condition cross-correlation.

    Parameters
    ----------
    s1, s2 : np.ndarray
        Input signals of equal length (1-D).
    params : TDSParams, optional
        Algorithm parameters. Uses defaults if not provided.

    Returns
    -------
    tau : np.ndarray
        Time delay at each window (samples). Negative τ means ``s1`` leads
        ``s2``; positive τ means ``s2`` leads ``s1`` (see ``pbc_xcorr``).
    t_vec : np.ndarray
        Centre-time of each window (sample index).
    cmax : np.ndarray
        Peak cross-correlation value at each window (signed).
    """
    if params is None:
        params = TDSParams()

    s1 = np.asarray(s1, dtype=float)
    s2 = np.asarray(s2, dtype=float)

    N = len(s1)
    L = params.window
    step = params.overlap
    max_lag = params.max_lag

    # All valid window start positions (works for any L and step)
    starts = list(range(0, N - L + 1, step))
    n_windows = len(starts)

    tau = np.full(n_windows, np.nan)
    t_vec = np.full(n_windows, np.nan)
    cmax = np.full(n_windows, np.nan)

    for idx, start in enumerate(starts):
        end = start + L

        seg1 = zscore(s1[start:end])
        seg2 = zscore(s2[start:end])

        if params.window_anchor == "end":
            t_vec[idx] = end
        else:
            t_vec[idx] = start + L // 2

        # Flat segment → skip (leave as NaN → treated as unstable)
        if np.all(seg1 == 0) or np.all(seg2 == 0):
            tau[idx] = max_lag + 1    # sentinel: outside valid range
            cmax[idx] = 0.0
            continue

        lags, C = pbc_xcorr(seg1, seg2, max_lag)

        abs_C = np.abs(C)
        best = np.argmax(abs_C)
        tau[idx] = lags[best]
        cmax[idx] = C[best]        # signed: preserves direction of coupling

    # Trim to actual computed windows
    valid = ~np.isnan(t_vec)
    return tau[valid], t_vec[valid], cmax[valid]


def stable_label(
    tau: np.ndarray,
    params: TDSParams = None,
) -> np.ndarray:
    """
    Step 2: label each τ₀ value as stable (1) or unstable (0).

    For each sliding window of length stability_window, every unique τ₀ is
    tried as a candidate delay. If at least stability_min points fall within
    ±tolerance of it, those points are marked stable. Output length = len(tau).

    Which points within a winning window get labeled is controlled by
    ``params.stability_anchor`` (see ``TDSParams.stability_anchor``):

    - ``"any"`` (default): non-causal / look-ahead. A label at index ``j``
      may be set by any window start in ``[j - stability_window + 1, j]``,
      so it depends on future tau values up to ``j + stability_window - 1``.
      Matches the published offline definition (Bashan et al. 2012). Not
      safe for walk-forward backtests or real-time use.
    - ``"end"``: causal. Only index ``i + stability_window - 1`` of each
      window is writable, and only if that point is within ``tolerance`` of
      the winning candidate. Guarantee: ``stable_label(tau, p)[:k] ==
      stable_label(tau[:k], p)`` for every ``k``. The first
      ``stability_window - 1`` labels are always 0. Labels are always a
      subset of the ``"any"`` labels.

    Parameters
    ----------
    tau : np.ndarray
        Time delay series from time_delay_interaction().
    params : TDSParams, optional

    Returns
    -------
    stbl_lbl : np.ndarray
        Binary array, same length as tau. 1 = stable, 0 = unstable.
    """
    if params is None:
        params = TDSParams()

    anchor = params.stability_anchor
    if anchor not in VALID_STABILITY_ANCHORS:
        raise ValueError(
            f"stability_anchor must be one of {VALID_STABILITY_ANCHORS}, "
            f"got {anchor!r}"
        )

    N = len(tau)
    win = params.stability_window      # default 5
    min_stable = params.stability_min  # default 4
    tol = params.tolerance             # default 1
    max_lag = params.max_lag           # sentinel threshold

    stbl_lbl = np.zeros(N, dtype=int)

    for i in range(N - win + 1):
        seg = tau[i: i + win]

        # Skip windows containing NaN or sentinel values
        if np.any(np.isnan(seg)) or np.any(np.abs(seg) > max_lag):
            continue

        for d in np.unique(seg):
            in_range = np.abs(seg - d) <= tol
            if np.sum(in_range) >= min_stable:
                if anchor == "end":
                    if in_range[-1]:
                        stbl_lbl[i + win - 1] = 1
                else:
                    indices = np.where(in_range)[0] + i
                    stbl_lbl[indices] = 1
                break

    return stbl_lbl


def tds_score(stbl_lbl: np.ndarray) -> float:
    """
    Step 3: compute the TDS score as the percentage of stable points.

    Parameters
    ----------
    stbl_lbl : np.ndarray
        Binary stability label array from stable_label().

    Returns
    -------
    float
        TDS score in [0.0, 100.0].
    """
    if len(stbl_lbl) == 0:
        return 0.0
    return float(np.mean(stbl_lbl) * 100)


def tds(
    s1: np.ndarray,
    s2: np.ndarray,
    params: TDSParams = None,
) -> dict:
    """
    Full TDS pipeline: runs all three steps and returns everything.

    Parameters
    ----------
    s1, s2 : np.ndarray
        Input signals of equal length.
    params : TDSParams, optional

    Returns
    -------
    dict with keys:
        'score'     : float   — TDS score in [0, 100]
        'tau'       : ndarray — time delay series (τ < 0 → s1 leads s2,
                                τ > 0 → s2 leads s1)
        't_vec'     : ndarray — window centre times
        'cmax'      : ndarray — peak cross-correlation per window
        'stbl_lbl'  : ndarray — binary stable/unstable labels
        'stable_taus': ndarray — τ₀ values at stable points only

    'stbl_lbl' and 'score' honour params.stability_anchor (see
    TDSParams.stability_anchor and stable_label()).
    """
    if params is None:
        params = TDSParams()

    tau, t_vec, cmax = time_delay_interaction(s1, s2, params)
    stbl_lbl = stable_label(tau, params)
    score = tds_score(stbl_lbl)
    stable_taus = tau[stbl_lbl == 1]

    return {
        "score": score,
        "tau": tau,
        "t_vec": t_vec,
        "cmax": cmax,
        "stbl_lbl": stbl_lbl,
        "stable_taus": stable_taus,
    }
