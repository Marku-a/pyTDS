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

    Sign convention
    ----------------
    τ < 0 means s1 leads s2 (s1's pattern appears first, s2 follows).
    τ > 0 means s2 leads s1. Swapping s1 and s2 flips the sign of τ.
    See ``pbc_xcorr`` for the underlying convention.
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

    Order-independent rule: for each sliding window of length
    ``stability_window``, every candidate delay ``d`` present in the window
    (not just the first qualifying one) is evaluated. A candidate qualifies
    if at least ``stability_min`` of the points in the window fall within
    ``±tolerance`` of ``d``. Point ``j`` is labeled stable iff *some* window
    containing ``j`` has a qualifying candidate whose in-range set includes
    ``j``. Output length = len(tau).

    Which points within a winning window get labeled is controlled by
    ``params.stability_anchor`` (see ``TDSParams.stability_anchor``):

    - ``"any"`` (default): non-causal / look-ahead. A label at index ``j``
      may be set by any window start in ``[j - stability_window + 1, j]``,
      so it depends on future tau values up to ``j + stability_window - 1``.
      Union over all qualifying candidates' in-range points in that window.
      Matches the published offline definition (Bashan et al. 2012). Not
      safe for walk-forward backtests or real-time use.
    - ``"end"``: causal. Only index ``i + stability_window - 1`` of each
      window is writable, and only if ANY qualifying candidate in that
      window has that last point within tolerance. Guarantee:
      ``stable_label(tau, p)[:k] == stable_label(tau[:k], p)`` for every
      ``k``. The first ``stability_window - 1`` labels are always 0. Labels
      are always a subset of the ``"any"`` labels.

    Deviation note: versions <=0.1.0 stopped at the first qualifying
    candidate in ascending order (``np.unique`` order) and broke out of the
    loop, so results depended on candidate order and were biased toward
    negative lags (e.g. tau=[0,0,0,-1,1] and tau=[0,0,0,1,-1] produced
    different last labels). The ``break`` was present since the initial
    port; no MATLAB reference source is available in this repo to confirm
    whether it originated there. The published definition (Bashan et al.
    2012) does not specify a candidate evaluation order, so this
    implementation now evaluates every candidate and is order-independent.
    New labels (>=0.2.0) are always a superset of the old (<=0.1.0) labels.

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

        # Vectorized: pairwise distance between every point in the window
        # and every other point, used as candidate delays.
        dist = np.abs(seg[:, None] - seg[None, :])
        in_range_matrix = dist <= tol            # rows = candidate d = seg[c]
        qualifies = in_range_matrix.sum(axis=1) >= min_stable

        if not np.any(qualifies):
            continue

        if anchor == "end":
            if np.any(in_range_matrix[qualifies, -1]):
                stbl_lbl[i + win - 1] = 1
        else:
            union = np.any(in_range_matrix[qualifies], axis=0)
            indices = np.where(union)[0] + i
            stbl_lbl[indices] = 1

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

    Sign convention: 'tau' < 0 means s1 leads s2; 'tau' > 0 means s2 leads
    s1. Swapping s1 and s2 flips the sign of 'tau' (see
    time_delay_interaction() / pbc_xcorr()).
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
