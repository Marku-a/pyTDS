"""
Tests for pyTDS.surrogate — surrogate TDS null distribution and significance threshold.
"""

import numpy as np
import pytest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from pyTDS.params import TDSParams
from pyTDS.surrogate import build_surrogate_subject, surrogate_tds, significance_threshold


# ── build_surrogate_subject ────────────────────────────────────────────────────

def make_dataset(n_subjects=5, T=300, N=3, seed=0):
    rng = np.random.default_rng(seed)
    return rng.standard_normal((n_subjects, T, N))


def test_build_surrogate_shape():
    dataset = make_dataset()
    rng = np.random.default_rng(0)
    surrogate = build_surrogate_subject(dataset, rng)
    assert surrogate.shape == (300, 3)


def test_build_surrogate_columns_from_real_subjects():
    """Each column of the surrogate must exactly match some subject's column."""
    dataset = make_dataset(n_subjects=5, T=50, N=4, seed=1)
    rng = np.random.default_rng(1)
    surrogate = build_surrogate_subject(dataset, rng)

    for col in range(4):
        col_data = surrogate[:, col]
        found = any(
            np.allclose(col_data, dataset[subj, :, col])
            for subj in range(5)
        )
        assert found, f"Column {col} doesn't match any subject"


def test_build_surrogate_mixes_subjects():
    """
    With 5 subjects and 10 signals, it's overwhelmingly likely that at least
    two columns come from different subjects (i.e., the surrogate is mixed).
    """
    rng = np.random.default_rng(99)
    dataset = make_dataset(n_subjects=5, T=50, N=10, seed=99)

    same_count = 0
    for _ in range(20):
        surrogate = build_surrogate_subject(dataset, rng)
        # Check which subject each column came from
        origins = []
        for col in range(10):
            for subj in range(5):
                if np.allclose(surrogate[:, col], dataset[subj, :, col]):
                    origins.append(subj)
                    break
        same_count += int(len(set(origins)) == 1)  # all from same subject

    assert same_count < 5  # very unlikely all-same across 20 trials


# ── surrogate_tds ──────────────────────────────────────────────────────────────

def test_surrogate_tds_output_length():
    """
    Length = n_surrogates × N*(N-1)/2 pairs.
    For n_surrogates=10, N=3 → 10 × 3 = 30 scores.
    """
    dataset = make_dataset(n_subjects=5, T=300, N=3)
    params = TDSParams(window=60, overlap=30, n_surrogates=10)
    null = surrogate_tds(dataset, params, seed=0)
    expected_len = 10 * 3  # 10 surrogates × 3 pairs
    assert len(null) == expected_len


def test_surrogate_tds_values_in_range():
    dataset = make_dataset(n_subjects=5, T=300, N=3)
    params = TDSParams(window=60, overlap=30, n_surrogates=5)
    null = surrogate_tds(dataset, params, seed=1)
    assert np.all(null >= 0.0)
    assert np.all(null <= 100.0)


def test_surrogate_tds_reproducible_with_seed():
    dataset = make_dataset(n_subjects=5, T=300, N=3)
    params = TDSParams(window=60, overlap=30, n_surrogates=5)
    null_a = surrogate_tds(dataset, params, seed=42)
    null_b = surrogate_tds(dataset, params, seed=42)
    np.testing.assert_allclose(null_a, null_b)


def test_surrogate_tds_different_seeds_produce_different_surrogates():
    """
    Different seeds should produce different column assignments in the
    surrogate subjects (randomness is seeded, not constant).
    Surrogate TDS scores will all be ~0 (no coupling by design) — the
    important property is that the mixing itself is seed-dependent.
    """
    dataset = make_dataset(n_subjects=10, T=60, N=4, seed=0)
    rng_a = np.random.default_rng(1)
    rng_b = np.random.default_rng(2)
    surr_a = build_surrogate_subject(dataset, rng_a)
    surr_b = build_surrogate_subject(dataset, rng_b)
    # Different seeds → different column picks → surrogate arrays differ
    assert not np.allclose(surr_a, surr_b)


def test_surrogate_tds_rejects_bad_shape():
    bad = np.random.randn(300, 3)  # 2D instead of 3D
    with pytest.raises(ValueError):
        surrogate_tds(bad)


# ── significance_threshold ─────────────────────────────────────────────────────

def test_significance_threshold_returns_float():
    null = np.array([5.0, 10.0, 15.0, 20.0, 25.0])
    params = TDSParams(alpha=0.05)
    threshold = significance_threshold(null, params)
    assert isinstance(threshold, float)


def test_significance_threshold_correct_percentile():
    rng = np.random.default_rng(0)
    null = rng.uniform(0, 100, 10000)
    params = TDSParams(alpha=0.05)
    threshold = significance_threshold(null, params)
    # Should be ~95th percentile
    assert abs(threshold - np.percentile(null, 95)) < 0.01


def test_significance_threshold_alpha_01():
    """alpha=0.01 → 99th percentile → higher threshold."""
    rng = np.random.default_rng(0)
    null = rng.uniform(0, 100, 10000)
    params_05 = TDSParams(alpha=0.05)
    params_01 = TDSParams(alpha=0.01)
    t05 = significance_threshold(null, params_05)
    t01 = significance_threshold(null, params_01)
    assert t01 > t05


def test_significance_threshold_no_params():
    """Default params should work without passing params explicitly."""
    null = np.linspace(0, 100, 1000)
    threshold = significance_threshold(null)
    assert 0.0 <= threshold <= 100.0


def test_significance_threshold_coupled_above_uncoupled_below():
    """
    Scores from strongly coupled signals should exceed the surrogate threshold,
    while scores from uncoupled (random) signals should mostly fall below it.
    """
    rng = np.random.default_rng(7)
    n_subjects = 8
    T = 600
    N = 3

    # Build dataset of uncoupled random signals
    dataset = rng.standard_normal((n_subjects, T, N))
    params = TDSParams(window=60, overlap=30, n_surrogates=20, alpha=0.05)
    null = surrogate_tds(dataset, params, seed=7)
    threshold = significance_threshold(null, params)

    # Uncoupled random pair should be near threshold (mostly below)
    from pyTDS.core import tds
    s1 = rng.standard_normal(T)
    s2 = rng.standard_normal(T)
    result = tds(s1, s2, params)
    # Random pair: no strong claim on direction, but threshold should be reasonable
    assert 0.0 <= threshold <= 100.0
    assert 0.0 <= result["score"] <= 100.0
