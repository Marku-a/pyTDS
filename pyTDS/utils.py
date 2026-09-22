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
        Maximum lag to return (returns lags from -max_lag to +max_lag).

    Returns
    -------
    lags : np.ndarray
        Array of lag values: [-max_lag, ..., 0, ..., +max_lag].
    C : np.ndarray
        Normalized cross-correlation values at each lag.
    """
    N = len(seg1)

    # FFT-based circular cross-correlation
    # C[k] = sum_i seg1[i] * seg2[(i + k) % N]  (periodic shift)
    xc_full = np.fft.ifft(np.fft.fft(seg1) * np.conj(np.fft.fft(seg2))).real

    # Normalize to correlation coefficient in [-1, 1]
    # (equivalent to dividing by N * std1 * std2, but both are already z-scored so std≈1)
    xc_full /= N

    # Rearrange: FFT output has positive lags first, then negative lags
    # Positive lags: xc_full[0 .. max_lag]
    # Negative lags: xc_full[N-max_lag .. N-1]  (wrap-around)
    pos = xc_full[:max_lag + 1]           # lags 0 .. +max_lag
    neg = xc_full[N - max_lag: N]         # lags -max_lag .. -1

    lags = np.concatenate([np.arange(-max_lag, 0), np.arange(0, max_lag + 1)])
    C = np.concatenate([neg, pos])

    return lags, C
