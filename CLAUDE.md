# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Project Does

**pyTDS** is a Python implementation of the Time Delay Stability (TDS) algorithm — a method for detecting whether two time-series signals maintain a stable fixed time delay, indicating stable coupling. Developed for physiological network analysis (EEG, heart rate, SpO2) based on Bashan et al. (Nature Communications 2012) and Bartsch et al. (PLoS ONE 2015).

## Commands

```bash
# Install for development
pip install -e ".[dev]"

# Run all tests (46 total)
python -m pytest tests/

# Run a single test file
python -m pytest tests/test_core.py

# Run a single test
python -m pytest tests/test_core.py::test_tds_score_basic

# Run tests with coverage
python -m pytest tests/ --cov=pyTDS

# Generate Excel test dashboard report (timestamped .xlsx in "tests reports/")
python run_dashboard.py
```

No linting or formatting tools are configured.

## Architecture

The algorithm is split into three stages and five modules:

**Stage 1 — Signal processing** (`utils.py`, `core.py`):
- `pbc_xcorr()` computes periodic-boundary-condition (circular) cross-correlation via FFT. This preserves all N time points (unlike linear cross-correlation) and is the core numerical primitive.
- `time_delay_interaction()` slides a window over two signals, calling `pbc_xcorr` in each window to extract the dominant time delay τ₀(t).
- `stable_label()` classifies each τ₀ as stable if ≥4 of 5 consecutive points agree within ±`tolerance` samples.
- `tds_score()` returns the percentage of stable-labeled windows.
- `tds()` orchestrates all three steps and returns a dict: `{score, tau, stable_taus, unstable_taus}`.

**Stage 2 — Network analysis** (`network.py`):
- `tds_matrix(signals)` runs pairwise `tds()` across all N signals (T×N input → N×N score matrix).
- `fix_symmetry()` symmetrizes the matrix and zeroes the diagonal.
- `apply_threshold()` binarizes scores into an adjacency matrix.
- `to_networkx()` converts to a `networkx.Graph` for graph analysis.

**Stage 3 — Significance testing** (`surrogate.py`):
- `build_surrogate_subject()` mixes signals across subjects to destroy inter-signal coupling while preserving each signal's temporal structure (autocorrelation/spectrum).
- `surrogate_tds()` pools null TDS scores from many surrogate subjects.
- `significance_threshold()` derives a percentile-based cutoff for `apply_threshold()`.

**Configuration** (`params.py`):
- `TDSParams` dataclass is the single source of truth for all algorithm parameters: `window`, `overlap`, `max_lag`, `tolerance`, `n_surrogates`, `alpha`. It is passed through all functions.

**Public API** (`__init__.py`):
- Exports `tds`, `TDSParams`, and the network/surrogate/viz submodules.

### Data Flow

```
Two signals (1D)  →  tds(s1, s2, params)  →  {score, tau, ...}

Signals (T×N)  →  tds_matrix()  →  N×N matrix
                                        ↓
                               fix_symmetry()
                                        ↓
                               apply_threshold(threshold)  ← significance_threshold() from surrogates
                                        ↓
                               to_networkx()  →  Graph
```

## Key Conventions

- **Python ≥ 3.10** — `|` union syntax is used in type hints throughout.
- All core functions accept `params: TDSParams | None = None`; a default `TDSParams()` is constructed internally when `None` is passed.
- The `data/` directory contains benchmark EEG signals (`data-delta.txt`, `data-sigma.txt`) used in `test_core.py` to validate against the known 43.6% TDS score from the original MATLAB implementation.
- `examples/rulkov_validation.py` is the algorithmic ground-truth benchmark using Rulkov chaotic oscillators with known coupling delays.
