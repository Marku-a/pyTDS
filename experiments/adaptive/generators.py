"""Synthetic generators with known coupling delay (experimental)."""

from __future__ import annotations

import numpy as np
from scipy.signal import lfilter


def ar1(n: int, tc: float, rng: np.random.Generator) -> np.ndarray:
    """z-scored AR(1) with phi = exp(-1/tc); burn-in int(10*tc)+100 discarded."""
    phi = np.exp(-1.0 / tc)
    burn = int(10 * tc) + 100
    e = rng.standard_normal(n + burn)
    x = lfilter([1.0], [1.0, -phi], e)[burn:]
    return (x - x.mean()) / x.std()


def coupled_pair(n, tc_x, tc_y, delay, c, rng, block_len=None):
    """
    y = c*m*roll(x, delay) + sqrt(1-c^2)*ar1(n, tc_y); m = coupling mask.
    m is all ones, or (block_len) alternating blocks OFF/ON/OFF/... .
    Returns x, y, m. TDS convention: stable tau = -delay.
    """
    x = ar1(n, tc_x, rng)
    if block_len:
        m = ((np.arange(n) // block_len) % 2 == 1).astype(float)
    else:
        m = np.ones(n)
    y = c * m * np.roll(x, delay) + np.sqrt(1 - c ** 2) * ar1(n, tc_y, rng)
    return x, y, m
