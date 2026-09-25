"""
Low-level signal processing utilities for TDS.
"""

import numpy as np


def zscore(x: np.ndarray) -> np.ndarray:
    """
    Normalize a signal segment to zero mean and unit standard deviation.

    If the signal is flat (std == 0), returns an array of zeros to avoid
    division by zero — a flat segment produces NaN cross-correlations otherwise.

    Parameters
    ----------
    x : np.ndarray
        1-D signal segment.

    Returns
    -------
    np.ndarray
        Normalized segment.
    """
    std = np.std(x)
    if std == 0:
        return np.zeros_like(x, dtype=float)
    return (x - np.mean(x)) / std


def pbc_xcorr(seg1: np.ndarray, seg2: np.ndarray, max_lag: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Circular (periodic boundary condition) cross-correlation via FFT.

    This is the method used in the original NCOM paper (Bashan et al. 2012)
    and Ronny Bartsch's PLoS ONE 2015 implementation. Every lag uses all N
    data points, unlike MATLAB's xcorr which loses points at large lags.

    Parameters
    ----------
    seg1 : np.ndarray
        First signal segment (already z-scored).
    seg2 : np.ndarray
        Second signal segment (already z-scored), same length as seg1.
    max_lag : int
        Maximum lag to return (returns lags from -max_lag to +max_lag),
        subject to the aliasing limit described below.

    Returns
    -------
    lags : np.ndarray
        Array of lag values: [-m, ..., 0, ..., +m] where
        ``m = min(max_lag, (N - 1) // 2)``.
    C : np.ndarray
        Normalized cross-correlation values at each lag.

    Sign convention
    ----------------
    C[k] = sum_i seg1[(i + k) % N] * seg2[i]. The peak occurs at k < 0 when
    seg1 LEADS seg2 (seg2[i] ~= seg1[i - D] gives a peak at k = -D); the
    peak occurs at k > 0 when seg2 leads seg1. Swapping the two arguments
    flips the sign of the detected lag.

    Aliasing
    --------
    Under periodic boundary conditions a circular shift by k is
    indistinguishable from a shift by k - N, so lags with |k| > (N - 1) // 2
    alias onto the opposite sign (e.g. for even N, +N/2 and -N/2 are the
    exact same circular shift). To keep the returned lags unambiguous,
    ``pbc_xcorr`` only returns |k| <= (N - 1) // 2, even if ``max_lag`` is
    larger. Callers that rely on a fixed ``2 * max_lag + 1``-length output
    should account for this.
    """
    N = len(seg1)

    # Effective max lag: cap at (N - 1) // 2 to avoid returning the same
    # circular shift twice under two different signs (see "Aliasing" above).
    m = min(max_lag, (N - 1) // 2)

    # FFT-based circular cross-correlation
    # C[k] = sum_i seg1[(i + k) % N] * seg2[i]  (periodic shift)
    # ifft(F1 * conj(F2)) shifts seg1, not seg2 → negative k means seg1 leads.
    xc_full = np.fft.ifft(np.fft.fft(seg1) * np.conj(np.fft.fft(seg2))).real

    # Normalize to correlation coefficient in [-1, 1]
    # (equivalent to dividing by N * std1 * std2, but both are already z-scored so std≈1)
    xc_full /= N

    # Rearrange: FFT output has positive lags first, then negative lags
    # Positive lags: xc_full[0 .. m]
    # Negative lags: xc_full[N-m .. N-1]  (wrap-around)
    pos = xc_full[:m + 1]           # lags 0 .. +m
    neg = xc_full[N - m: N]         # lags -m .. -1

    lags = np.concatenate([np.arange(-m, 0), np.arange(0, m + 1)])
    C = np.concatenate([neg, pos])

    return lags, C
