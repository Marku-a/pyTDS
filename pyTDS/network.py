"""
Network module — build physiological interaction networks from TDS matrices.

Workflow:
  signals (T x N)
    -> tds_matrix()        : N x N TDS scores for every signal pair
    -> fix_symmetry()      : symmetrize, zero diagonal
    -> apply_threshold()   : binary adjacency matrix
    -> to_networkx()       : networkx Graph for analysis / plotting
"""

from __future__ import annotations
import numpy as np
from itertools import combinations
from .params import TDSParams
from .core import tds


def tds_matrix(
    signals: np.ndarray,
    params: TDSParams = None,
    verbose: bool = False,
) -> np.ndarray:
    """
    Compute the full N×N TDS matrix for a dataset of N signals.

    Parameters
    ----------
    signals : np.ndarray
        Shape (T, N) — T time points, N signals (columns).
    params : TDSParams, optional
    verbose : bool
        Print progress for large N.

    Returns
    -------
    np.ndarray
        Shape (N, N). Element [i, j] = TDS score between signal i and j.
        Diagonal is 0. Matrix is symmetrized (average of [i,j] and [j,i]).
    """
    if params is None:
        params = TDSParams()

    signals = np.asarray(signals, dtype=float)
    if signals.ndim == 1:
        raise ValueError("signals must be 2-D (T x N). Got 1-D array.")

    T, N = signals.shape
    mat = np.zeros((N, N), dtype=float)

    pairs = list(combinations(range(N), 2))
    total = len(pairs)

    for k, (i, j) in enumerate(pairs):
        if verbose:
            print(f"  TDS pair {k+1}/{total}  ({i},{j})", end="\r")
        result = tds(signals[:, i], signals[:, j], params)
        mat[i, j] = result["score"]
        mat[j, i] = result["score"]

    if verbose:
        print()

    return mat


def fix_symmetry(mat: np.ndarray) -> np.ndarray:
    """
    Symmetrize a TDS matrix and zero the diagonal.

    For a matrix computed with tds_matrix() the values are already symmetric,
    but if built from asymmetric estimates (e.g., directional measures) this
    averages the two directions.

    Parameters
    ----------
    mat : np.ndarray
        Shape (N, N).

    Returns
    -------
    np.ndarray
        Symmetrized matrix with zeroed diagonal.
    """
    mat = np.asarray(mat, dtype=float).copy()
    mat = (mat + mat.T) / 2.0
    np.fill_diagonal(mat, 0.0)
    return mat


def apply_threshold(
    tds_mat: np.ndarray,
    threshold: float,
) -> np.ndarray:
    """
    Apply a significance threshold to a TDS matrix → binary adjacency matrix.

    Parameters
    ----------
    tds_mat : np.ndarray
        Shape (N, N) — TDS score matrix.
    threshold : float
        Links with TDS > threshold are kept. Typically ~7% for physiology
        (Bashan et al. 2012, Nature Communications).

    Returns
    -------
    np.ndarray
        Shape (N, N) — binary adjacency matrix (0 or 1). Symmetric.
    """
    adj = (tds_mat > threshold).astype(int)
    np.fill_diagonal(adj, 0)
    return adj


def to_networkx(
    adj_mat: np.ndarray,
    tds_mat: np.ndarray = None,
    labels: list[str] = None,
):
    """
    Convert an adjacency matrix to a networkx Graph.

    Parameters
    ----------
    adj_mat : np.ndarray
        Shape (N, N) — binary adjacency matrix from apply_threshold().
    tds_mat : np.ndarray, optional
        Shape (N, N) — TDS scores used as edge weights.
    labels : list[str], optional
        Node labels. Defaults to '0', '1', ..., 'N-1'.

    Returns
    -------
    networkx.Graph
        Nodes are labelled, edges carry 'weight' = TDS score (if provided).
    """
    try:
        import networkx as nx
    except ImportError:
        raise ImportError("networkx is required: pip install networkx")

    N = adj_mat.shape[0]
    if labels is None:
        labels = [str(i) for i in range(N)]

    G = nx.Graph()
    G.add_nodes_from(labels)

    for i, j in combinations(range(N), 2):
        if adj_mat[i, j]:
            weight = float(tds_mat[i, j]) if tds_mat is not None else 1.0
            G.add_edge(labels[i], labels[j], weight=weight)

    return G
