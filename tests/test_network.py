"""
Tests for pyTDS.network — TDS matrix, thresholding, and networkx conversion.
"""

import numpy as np
import pytest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from pyTDS.params import TDSParams
from pyTDS.network import tds_matrix, fix_symmetry, apply_threshold, to_networkx


# ── Helpers ────────────────────────────────────────────────────────────────────

def make_signals(n_signals=4, T=300, seed=0):
    rng = np.random.default_rng(seed)
    return rng.standard_normal((T, n_signals))


# ── fix_symmetry ───────────────────────────────────────────────────────────────

def test_fix_symmetry_zeros_diagonal():
    mat = np.ones((4, 4)) * 5.0
    out = fix_symmetry(mat)
    assert np.all(np.diag(out) == 0.0)


def test_fix_symmetry_symmetric():
    rng = np.random.default_rng(1)
    mat = rng.uniform(0, 100, (5, 5))
    out = fix_symmetry(mat)
    np.testing.assert_allclose(out, out.T)


def test_fix_symmetry_preserves_symmetric_values():
    mat = np.array([[0, 30, 50],
                    [30, 0, 70],
                    [50, 70, 0]], dtype=float)
    out = fix_symmetry(mat)
    np.testing.assert_allclose(out, mat)


def test_fix_symmetry_averages_asymmetric():
    mat = np.array([[0, 10, 0],
                    [30, 0, 0],
                    [0, 0, 0]], dtype=float)
    out = fix_symmetry(mat)
    assert out[0, 1] == pytest.approx(20.0)
    assert out[1, 0] == pytest.approx(20.0)


# ── apply_threshold ────────────────────────────────────────────────────────────

def test_apply_threshold_basic():
    mat = np.array([[0, 80, 5],
                    [80, 0, 60],
                    [5, 60, 0]], dtype=float)
    adj = apply_threshold(mat, threshold=20.0)
    assert adj[0, 1] == 1
    assert adj[1, 2] == 1
    assert adj[0, 2] == 0
    assert np.all(np.diag(adj) == 0)


def test_apply_threshold_zeros_diagonal():
    mat = np.full((3, 3), 100.0)
    adj = apply_threshold(mat, threshold=50.0)
    assert np.all(np.diag(adj) == 0)


def test_apply_threshold_all_zeros():
    mat = np.zeros((4, 4))
    adj = apply_threshold(mat, threshold=5.0)
    assert np.all(adj == 0)


def test_apply_threshold_all_ones():
    mat = np.full((4, 4), 100.0)
    np.fill_diagonal(mat, 0)
    adj = apply_threshold(mat, threshold=5.0)
    expected_off_diag = np.ones((4, 4)) - np.eye(4)
    np.testing.assert_array_equal(adj, expected_off_diag.astype(int))


# ── tds_matrix ─────────────────────────────────────────────────────────────────

def test_tds_matrix_shape():
    signals = make_signals(n_signals=4, T=300)
    params = TDSParams(window=60, overlap=30)
    mat = tds_matrix(signals, params)
    assert mat.shape == (4, 4)


def test_tds_matrix_diagonal_zero():
    signals = make_signals(n_signals=3, T=300)
    params = TDSParams(window=60, overlap=30)
    mat = tds_matrix(signals, params)
    assert np.all(np.diag(mat) == 0.0)


def test_tds_matrix_symmetric():
    signals = make_signals(n_signals=3, T=300)
    params = TDSParams(window=60, overlap=30)
    mat = tds_matrix(signals, params)
    np.testing.assert_allclose(mat, mat.T)


def test_tds_matrix_values_in_range():
    signals = make_signals(n_signals=3, T=300)
    params = TDSParams(window=60, overlap=30)
    mat = tds_matrix(signals, params)
    assert np.all(mat >= 0.0)
    assert np.all(mat <= 100.0)


def test_tds_matrix_rejects_1d():
    s = np.random.randn(300)
    with pytest.raises(ValueError):
        tds_matrix(s)


def test_tds_matrix_coupled_pair_higher():
    """
    A known-coupled pair should have a higher TDS score than random noise.
    """
    rng = np.random.default_rng(42)
    T = 600
    base = np.sin(np.linspace(0, 20 * np.pi, T))
    shift = 5
    coupled = np.roll(base, shift) + 0.05 * rng.standard_normal(T)
    noise = rng.standard_normal(T)

    # Signals: [coupled_s1, coupled_s2, noise]
    signals = np.column_stack([base, coupled, noise])
    params = TDSParams(window=60, overlap=30)
    mat = tds_matrix(signals, params)

    # (0,1) should be higher than (0,2) and (1,2)
    assert mat[0, 1] > mat[0, 2]
    assert mat[0, 1] > mat[1, 2]


# ── to_networkx ────────────────────────────────────────────────────────────────

def test_to_networkx_node_count():
    pytest.importorskip("networkx")
    adj = np.array([[0, 1, 0],
                    [1, 0, 1],
                    [0, 1, 0]], dtype=int)
    G = to_networkx(adj)
    assert len(G.nodes) == 3


def test_to_networkx_edge_count():
    pytest.importorskip("networkx")
    adj = np.array([[0, 1, 0],
                    [1, 0, 1],
                    [0, 1, 0]], dtype=int)
    G = to_networkx(adj)
    assert len(G.edges) == 2


def test_to_networkx_custom_labels():
    pytest.importorskip("networkx")
    adj = np.array([[0, 1], [1, 0]], dtype=int)
    G = to_networkx(adj, labels=["HR", "SpO2"])
    assert "HR" in G.nodes
    assert "SpO2" in G.nodes
    assert G.has_edge("HR", "SpO2")


def test_to_networkx_edge_weights():
    pytest.importorskip("networkx")
    adj = np.array([[0, 1], [1, 0]], dtype=int)
    tds = np.array([[0.0, 75.0], [75.0, 0.0]])
    G = to_networkx(adj, tds_mat=tds)
    edge_data = G.get_edge_data(*list(G.edges)[0])
    assert edge_data["weight"] == pytest.approx(75.0)


def test_to_networkx_no_edges():
    pytest.importorskip("networkx")
    adj = np.zeros((3, 3), dtype=int)
    G = to_networkx(adj)
    assert len(G.edges) == 0
