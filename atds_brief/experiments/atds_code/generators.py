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


def _rulkov_f(x, y, alpha):
    if x <= 0:
        return alpha / (1.0 - x) + y
    elif x < alpha + y:
        return alpha + y
    return -1.0


def rulkov_pair(n, delay, g, rng, mu=0.001, alpha=(4.5, 4.1), sigma=(0.01, -0.01),
                beta=1.0, sigma_e=1.0, extra_transient=1000):
    """
    Copy of examples/rulkov_validation.coupled_rulkov_with_delay with randomised
    initial conditions: x0 ~ U(-1.5,-0.5), y0 ~ U(-4.5,-3.5) per oscillator.
    The x0 draw is reused in the second pass (original resets x0=-1); y0 is reset to
    the first-pass median as in the original. Discards 50 + extra_transient steps.
    """
    disc = 50 + extra_transient
    n_total = n + disc
    x = [np.zeros(n_total), np.zeros(n_total)]
    y = [np.zeros(n_total), np.zeros(n_total)]
    x0 = rng.uniform(-1.5, -0.5, 2)
    y0 = rng.uniform(-4.5, -3.5, 2)
    for k in range(2):
        x[k][0], y[k][0] = x0[k], y0[k]
    for _ in range(2):
        for i in range(1, n_total):
            for k in range(2):
                o = 1 - k
                if delay >= i:
                    b = s = 0.0
                else:
                    b = g * beta * (x[o][i - 1 - delay] - x[k][i - 1])
                    s = g * sigma_e * (x[o][i - 1 - delay] - x[k][i - 1])
                x[k][i] = _rulkov_f(x[k][i - 1], y[k][i - 1] + b, alpha[k])
                y[k][i] = y[k][i - 1] - mu * (x[k][i - 1] + 1) + mu * sigma[k] + mu * s
        for k in range(2):
            y[k][0] = np.median(y[k])
            x[k][0] = x0[k]
    return x[0][disc:], x[1][disc:]
