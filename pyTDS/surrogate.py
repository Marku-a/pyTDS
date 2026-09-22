"""
Surrogate testing module — dataset-level null distribution for TDS significance.

The surrogate approach (from Bashan et al. 2012, Nature Communications):
  - Build "fake subjects" by taking each signal from a different real subject
  - This preserves each signal's real statistical properties (spectrum,
    autocorrelation, amplitude) while destroying inter-signal coupling
  - Run TDS on all signal pairs of each fake subject
  - Pool all TDS scores → null distribution
  - Threshold = (1 - alpha) percentile of null distribution

Why this is better than simple time-series shuffling:
  - Shuffling destroys the autocorrelation structure of signals
  - Cross-subject mixing keeps real physiological signal statistics intact
  - The null represents "what TDS would you get from real but unrelated signals?"
"""

from __future__ import annotations
import numpy as np
from itertools import combinations
from .params import TDSParams
from .core import tds


def build_surrogate_subject(
    dataset: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    """
    Build one surrogate subject by drawing each signal from a different
    randomly chosen real subject.

    Parameters
    ----------
    dataset : np.ndarray
        Shape (n_subjects, T, N) — subjects x time x signals.
    rng : np.random.Generator
        Random number generator for reproducibility.

    Returns
    -------
    np.ndarray
        Shape (T, N) — surrogate subject with mixed signals.
        Signal column j comes from a randomly chosen subject.
    """
    n_subjects, T, N = dataset.shape
    surrogate = np.empty((T, N), dtype=float)

    for j in range(N):
        subj_idx = rng.integers(0, n_subjects)
        surrogate[:, j] = dataset[subj_idx, :, j]

    return surrogate


def surrogate_tds(
    dataset: np.ndarray,
    params: TDSParams = None,
    seed: int = None,
    verbose: bool = False,
) -> np.ndarray:
    """
    Build a null TDS distribution using cross-subject surrogate mixing.

    Parameters
    ----------
    dataset : np.ndarray
        Shape (n_subjects, T, N) — subjects x time x signals.
        All subjects must have the same T and N.
    params : TDSParams, optional
        Uses params.n_surrogates for the number of surrogate subjects.
    seed : int, optional
        Random seed for reproducibility.
    verbose : bool
        Print progress.

    Returns
    -------
    np.ndarray
        1-D array of pooled null TDS scores from all surrogate pairs
        and all surrogate subjects. Length = n_surrogates * N*(N-1)/2.
    """
    if params is None:
        params = TDSParams()

    dataset = np.asarray(dataset, dtype=float)
    if dataset.ndim != 3:
        raise ValueError("dataset must be 3-D: (n_subjects, T, N)")

    n_subjects, T, N = dataset.shape
    n_surrogates = params.n_surrogates
    pairs = list(combinations(range(N), 2))
    rng = np.random.default_rng(seed)

    null_scores = []

    for s in range(n_surrogates):
        if verbose:
            pct = (s + 1) / n_surrogates * 100
            bar_len = 30
            filled  = int(bar_len * (s + 1) / n_surrogates)
            bar     = "#" * filled + "-" * (bar_len - filled)
            print(f"\r  [{bar}] {pct:5.1f}%  ({s+1}/{n_surrogates})", end="", flush=True)

        surrogate = build_surrogate_subject(dataset, rng)

        for i, j in pairs:
            result = tds(surrogate[:, i], surrogate[:, j], params)
            null_scores.append(result["score"])

    if verbose:
        print()

    return np.array(null_scores)


def significance_threshold(
    null_scores: np.ndarray,
    params: TDSParams = None,
) -> float:
    """
    Compute the significance threshold from a null TDS distribution.

    A link is considered significant if its TDS score exceeds this threshold.
    In Bashan et al. 2012, the threshold was ~7% TDS (all links above it
    were statistically significant at p < 0.001).

    Parameters
    ----------
    null_scores : np.ndarray
        Pooled null TDS scores from surrogate_tds().
    params : TDSParams, optional
        Uses params.alpha for the significance level (default 0.05).

    Returns
    -------
    float
        The (1 - alpha) percentile of null_scores.
    """
    if params is None:
        params = TDSParams()

    return float(np.percentile(null_scores, (1.0 - params.alpha) * 100.0))
