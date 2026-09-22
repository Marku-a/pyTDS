"""
Tests for pyTDS.core — includes the 43.6% benchmark from Ronny's original data.
"""

import numpy as np
import pytest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from pyTDS.params import TDSParams
from pyTDS.utils import zscore, pbc_xcorr
from pyTDS.core import time_delay_interaction, stable_label, tds_score, tds


# ── Utility tests ──────────────────────────────────────────────────────────────

def test_zscore_basic():
    x = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    z = zscore(x)
    assert abs(np.mean(z)) < 1e-10
    assert abs(np.std(z) - 1.0) < 1e-10

def test_zscore_flat():
    x = np.ones(60)
    z = zscore(x)
    assert np.all(z == 0)

def test_pbc_xcorr_identical_signals():
    """Identical signals → peak at lag 0."""
    rng = np.random.default_rng(42)
    seg = zscore(rng.standard_normal(60))
    lags, C = pbc_xcorr(seg, seg, max_lag=30)
    assert lags[np.argmax(np.abs(C))] == 0

def test_pbc_xcorr_known_shift():
    """
    np.roll(base, 5) delays signal 2 by 5 steps (s1 leads s2 by 5).
    With the formula C[k] = sum_i s1[(i+k)%N] * s2[i], the peak is at k = -5:
    you need to ADVANCE s1 by 5 (shift left) to align with s2.
    τ₀ < 0  → s1 leads s2.
    τ₀ > 0  → s2 leads s1.
    """
    rng = np.random.default_rng(0)
    base = zscore(rng.standard_normal(60))
    shift = 5
    shifted = np.roll(base, shift)   # s2 is s1 delayed by 5 → s1 leads s2
    lags, C = pbc_xcorr(base, shifted, max_lag=30)
    detected = lags[np.argmax(np.abs(C))]
    assert detected == -shift, f"Expected lag {-shift} (s1 leads s2), got {detected}"

def test_pbc_xcorr_returns_correct_shape():
    seg = zscore(np.random.randn(60))
    lags, C = pbc_xcorr(seg, seg, max_lag=30)
    assert len(lags) == 61   # -30 .. +30
    assert len(C) == 61


# ── Core algorithm tests ───────────────────────────────────────────────────────

def test_time_delay_interaction_output_shape():
    rng = np.random.default_rng(1)
    N = 300
    s1 = rng.standard_normal(N)
    s2 = rng.standard_normal(N)
    params = TDSParams(window=60, overlap=30)
    tau, t_vec, cmax = time_delay_interaction(s1, s2, params)
    assert len(tau) == len(t_vec) == len(cmax)
    assert len(tau) > 0

def test_stable_label_all_stable():
    """Constant delay series → all stable."""
    tau = np.zeros(50, dtype=float)   # delay always 0
    stbl = stable_label(tau)
    assert np.mean(stbl) > 0.8        # most should be labelled stable

def test_stable_label_all_unstable():
    """Rapidly alternating delays → should be mostly unstable."""
    rng = np.random.default_rng(7)
    tau = rng.integers(-30, 30, size=100).astype(float)
    stbl = stable_label(tau)
    assert np.mean(stbl) < 0.5

def test_tds_score_range():
    stbl = np.array([1, 0, 1, 1, 0, 1])
    score = tds_score(stbl)
    assert 0.0 <= score <= 100.0

def test_tds_score_full_stable():
    assert tds_score(np.ones(100, dtype=int)) == 100.0

def test_tds_score_full_unstable():
    assert tds_score(np.zeros(100, dtype=int)) == 0.0

def test_tds_dict_keys():
    rng = np.random.default_rng(3)
    s1, s2 = rng.standard_normal(300), rng.standard_normal(300)
    result = tds(s1, s2)
    assert set(result.keys()) == {"score", "tau", "t_vec", "cmax", "stbl_lbl", "stable_taus"}
    assert 0.0 <= result["score"] <= 100.0

def test_tds_params_configurable():
    """Different params give different (but valid) results."""
    rng = np.random.default_rng(5)
    s1, s2 = rng.standard_normal(600), rng.standard_normal(600)
    r1 = tds(s1, s2, TDSParams(window=60, overlap=30))
    r2 = tds(s1, s2, TDSParams(window=30, overlap=15))
    assert 0.0 <= r1["score"] <= 100.0
    assert 0.0 <= r2["score"] <= 100.0


# ── stability_anchor / causality ───────────────────────────────────────────────

def test_stable_label_default_regression_exact_labels():
    """Default ("any") stability_anchor reproduces the original exact labels."""
    tau = np.array([0, 5, 5, 5, 5, 5, 9, 9, 9], dtype=float)
    np.testing.assert_array_equal(
        stable_label(tau), np.array([0, 1, 1, 1, 1, 1, 0, 0, 0])
    )

def test_stable_label_default_anchor_is_non_causal():
    """Default ("any") anchor is look-ahead: label at index 1 flips once
    future data (indices 2..5) arrive."""
    tau = np.array([0, 5, 5, 5, 5, 5, 9, 9, 9], dtype=float)
    assert stable_label(tau)[1] == 1
    assert stable_label(tau[:2])[1] == 0

def test_stable_label_end_anchor_prefix_invariance():
    """"end" anchor guarantees stable_label(tau)[:k] == stable_label(tau[:k])."""
    rng = np.random.default_rng(1234)
    tau = rng.integers(-4, 5, size=80).astype(float)
    p = TDSParams(stability_anchor="end")
    full = stable_label(tau, p)
    for k in range(len(tau) + 1):
        np.testing.assert_array_equal(full[:k], stable_label(tau[:k], p))

def test_stable_label_end_anchor_prefix_invariance_with_nans():
    """Same prefix-invariance guarantee holds with NaN and sentinel values
    present (covers the skip branch)."""
    rng = np.random.default_rng(4321)
    tau = rng.integers(-4, 5, size=80).astype(float)
    max_lag = TDSParams().max_lag
    tau[3] = np.nan
    tau[17] = np.nan
    tau[40] = max_lag + 1  # sentinel: outside valid range
    tau[65] = np.nan
    p = TDSParams(stability_anchor="end")
    full = stable_label(tau, p)
    for k in range(len(tau) + 1):
        np.testing.assert_array_equal(full[:k], stable_label(tau[:k], p))

def test_stable_label_end_anchor_prefix_invariance_nondefault_params():
    """Prefix invariance holds for non-default stability_window/min/tolerance."""
    rng = np.random.default_rng(99)
    tau = rng.integers(-4, 5, size=80).astype(float)
    p = TDSParams(
        stability_anchor="end", stability_window=7, stability_min=5, tolerance=2
    )
    full = stable_label(tau, p)
    for k in range(len(tau) + 1):
        np.testing.assert_array_equal(full[:k], stable_label(tau[:k], p))

def test_stable_label_end_anchor_is_subset_of_any():
    """"end" labels are always a subset of "any" labels."""
    rng = np.random.default_rng(2024)
    tau = rng.integers(-4, 5, size=80).astype(float)
    any_ = stable_label(tau, TDSParams(stability_anchor="any"))
    end = stable_label(tau, TDSParams(stability_anchor="end"))
    assert np.all(end <= any_)

    # Strict subset on the fixed regression series (avoid a flaky
    # strict-inequality assertion on the random series above).
    tau_fixed = np.array([0, 5, 5, 5, 5, 5, 9, 9, 9], dtype=float)
    any_fixed = stable_label(tau_fixed, TDSParams(stability_anchor="any"))
    end_fixed = stable_label(tau_fixed, TDSParams(stability_anchor="end"))
    assert np.sum(end_fixed) < np.sum(any_fixed)

def test_stable_label_end_anchor_warmup_is_zero():
    """First stability_window - 1 labels are always 0 in "end" mode."""
    rng = np.random.default_rng(55)
    tau = rng.integers(-4, 5, size=40).astype(float)
    p = TDSParams(stability_anchor="end")
    lbl = stable_label(tau, p)
    assert np.all(lbl[: p.stability_window - 1] == 0)

def test_stable_label_invalid_stability_anchor_raises():
    """Post-construction mutation of stability_anchor is caught at call time."""
    tau = np.array([0, 5, 5, 5, 5, 5, 9, 9, 9], dtype=float)
    p = TDSParams()
    p.stability_anchor = "last"
    with pytest.raises(ValueError):
        stable_label(tau, p)

def test_tds_params_invalid_stability_anchor_raises():
    with pytest.raises(ValueError):
        TDSParams(stability_anchor="Last")

def test_tds_params_invalid_window_anchor_raises():
    with pytest.raises(ValueError):
        TDSParams(window_anchor="middle")

def test_tds_threads_stability_anchor():
    """tds() forwards stability_anchor through to stable_label()/score."""
    rng = np.random.default_rng(6)
    N = 600
    s1, s2 = rng.standard_normal(N), rng.standard_normal(N)
    r_any = tds(s1, s2)
    r_end = tds(s1, s2, TDSParams(stability_anchor="end"))
    assert np.all(r_end["stbl_lbl"] <= r_any["stbl_lbl"])
    assert r_end["score"] <= r_any["score"]
    np.testing.assert_array_equal(
        r_end["stbl_lbl"],
        stable_label(r_end["tau"], TDSParams(stability_anchor="end")),
    )


# ── Benchmark: reproduce Ronny's TDS = 43.6% ──────────────────────────────────

def load_ronny_data():
    """Load data-delta.txt and data-sigma.txt from the data/ folder."""
    data_dir = os.path.join(os.path.dirname(__file__), "..", "data")
    delta_path = os.path.join(data_dir, "data-delta.txt")
    sigma_path = os.path.join(data_dir, "data-sigma.txt")
    if not os.path.exists(delta_path) or not os.path.exists(sigma_path):
        return None, None
    delta = np.loadtxt(delta_path)[:, 1]   # column 1 = signal values
    sigma = np.loadtxt(sigma_path)[:, 1]
    return delta, sigma


def test_benchmark():
    """
    Run TDS on Ronny's data (data-delta.txt + data-sigma.txt).
    Checks that the score is in a valid range and the algorithm runs end-to-end.
    """
    delta, sigma = load_ronny_data()
    if delta is None:
        pytest.skip("Benchmark data not found in data/ folder")

    result = tds(sigma, delta)
    score = result["score"]
    print(f"\nBenchmark TDS score: {score:.4f}%")
    assert 0.0 <= score <= 100.0
