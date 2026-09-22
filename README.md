# pyTDS — Time Delay Stability for Python

A clean, general-purpose Python implementation of the **Time Delay Stability (TDS)** algorithm.

TDS measures whether two signals are **stably coupled** over time — i.e., do they consistently maintain a fixed time delay between them? It was originally developed by Ronny Bartsch (PLoS ONE 2015) and applied to physiological network analysis in Bashan et al. (Nature Communications 2012).

---

## What is TDS?

Given two time-series signals, TDS answers:
> *"Are these signals consistently coupled with a fixed time delay, or is their relationship random and unstable?"*

A **high TDS score** means the signals maintain a stable lag — they are tightly coupled.  
A **low TDS score** means the relationship is erratic or non-existent.

### The Algorithm (3 steps)

**Step 1 — Sliding-window cross-correlation → τ₀(t)**

For each window of length `L`, compute the circular (PBC) cross-correlation between the two segments. The lag at the peak correlation is τ₀ — the dominant time delay at that moment.

**Step 2 — Stable labeling**

For each group of 5 consecutive τ₀ values, check if at least 4 out of 5 agree on roughly the same lag (within ±1 sample). If yes → those points are **stable** (the coupling is consistent). If no → **unstable**.

**Step 3 — TDS score**

```
TDS = (number of stable points / total points) × 100   →   [0%, 100%]
```

---

## Installation

```bash
git clone https://github.com/Marku-a/pyTDS.git
cd pyTDS
pip install -e .
```

**Dependencies:** `numpy`, `scipy`, `matplotlib`, `pandas`, `networkx`

---

## Quick Start

```python
import pyTDS
from pyTDS.params import TDSParams

# Two signals of equal length
result = pyTDS.tds(s1, s2)

print(result["score"])        # TDS score in [0, 100]
print(result["tau"])          # time delay series τ₀(t)
print(result["stable_taus"])  # τ₀ values at stable points only
```

### With custom parameters

```python
params = TDSParams(
    window=60,            # sliding window size (samples)
    overlap=30,           # step between windows (samples)
    max_lag=30,           # maximum lag to search ±30
    stability_window=5,   # points assessed together
    stability_min=4,      # min points that must agree
    tolerance=1,          # allowed deviation ±1 sample
    stability_anchor="any", # "any" (default, look-ahead) or "end" (causal)
)

result = pyTDS.tds(s1, s2, params)
```

---

## Full API Reference

### `TDSParams` — Configuration

All algorithm parameters in one place. Pass a single instance through all functions.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `window` | `60` | Sliding window size (samples) |
| `overlap` | `30` | Step between windows — `window/2` = 50% overlap |
| `max_lag` | `30` | Maximum cross-correlation lag searched (±samples) |
| `stability_window` | `5` | Number of consecutive τ₀ points assessed for stability |
| `stability_min` | `4` | Minimum points within tolerance to call stable (out of `stability_window`) |
| `tolerance` | `1` | Allowed deviation in τ₀ to still count as "same lag" (±samples) |
| `window_anchor` | `"center"` | Where to stamp each `t_vec` timestamp: `"center"` = window midpoint (`start + L//2`); `"end"` = last sample of window — the earliest moment the result could be known in real time |
| `stability_anchor` | `"any"` | Which points inside a stable window get labeled: `"any"` (default) = original Bashan et al. definition, **uses look-ahead**; `"end"` = causal, only the window's last point — safe for real-time / walk-forward use |
| `n_surrogates` | `1000` | Number of surrogate subjects for null distribution |
| `alpha` | `0.05` | Significance level — threshold = (1−α) percentile of null scores |

---

### `pyTDS.core` — Core Algorithm

#### `tds(s1, s2, params=None) → dict`

Full TDS pipeline. Returns a dictionary with:

| Key | Type | Description |
|-----|------|-------------|
| `score` | `float` | TDS score ∈ [0, 100] |
| `tau` | `ndarray` | Time delay series τ₀(t) |
| `t_vec` | `ndarray` | Window centre timestamps |
| `cmax` | `ndarray` | Peak cross-correlation per window |
| `stbl_lbl` | `ndarray` | Binary stable (1) / unstable (0) labels |
| `stable_taus` | `ndarray` | τ₀ values at stable points only |

#### `time_delay_interaction(s1, s2, params=None) → (tau, t_vec, cmax)`

Step 1 only. Returns the raw time delay series without stability labeling.

#### `stable_label(tau, params=None) → ndarray`

Step 2 only. Takes a τ₀ series and returns binary stability labels.

**Warning — look-ahead by default.** With the default `stability_anchor="any"`,
labeling is non-causal: the label at index `j` depends on τ₀ values up to
`j + stability_window - 1` (4 future windows at the default
`stability_window=5`) — only the final element of the series is causal. For
real-time or walk-forward/backtest use, pass `TDSParams(stability_anchor="end")`,
which guarantees `stable_label(tau)[:k] == stable_label(tau[:k])`, zeroes the
first `stability_window - 1` labels (warm-up), and yields scores no higher
than `"any"` (its labels are always a subset). Note that `window_anchor="end"` alone
does **not** make labels causal — both flags are needed together.

#### `tds_score(stbl_lbl) → float`

Step 3 only. Computes the TDS score from a binary label array.

---

### `pyTDS.network` — Physiological Networks

Build interaction networks from a set of N signals.

#### `tds_matrix(signals, params=None, verbose=False) → ndarray`

Compute the full N×N TDS matrix for a `(T, N)` array of signals.  
Element `[i, j]` is the TDS score between signal i and signal j.

```python
import numpy as np
from pyTDS.network import tds_matrix, fix_symmetry, apply_threshold, to_networkx

signals = np.column_stack([hr, spo2, thorax, abdomen])  # shape (T, 4)
mat = tds_matrix(signals, params)
```

#### `fix_symmetry(mat) → ndarray`

Symmetrize a TDS matrix: averages `(mat + mat.T) / 2` and zeros the diagonal.

#### `apply_threshold(tds_mat, threshold) → ndarray`

Apply a significance threshold → binary adjacency matrix.  
Links with `TDS > threshold` are kept (1), others removed (0).

```python
adj = apply_threshold(mat, threshold=7.0)   # ~7% is typical for physiological networks
```

#### `to_networkx(adj_mat, tds_mat=None, labels=None) → Graph`

Convert an adjacency matrix to a `networkx.Graph`.  
Edge weights = TDS scores (if `tds_mat` is provided).

```python
labels = ["HR", "SpO2", "Thorax", "Abdomen"]
G = to_networkx(adj, tds_mat=mat, labels=labels)
```

---

### `pyTDS.surrogate` — Significance Testing

Build a null TDS distribution using cross-subject surrogate mixing.

#### How it works

Instead of shuffling individual signals (which destroys autocorrelation structure), the surrogate method builds "fake subjects" where **each signal comes from a different real subject**. This preserves the real statistical properties of every signal while destroying inter-signal coupling. Running TDS on all signal pairs of many such fake subjects gives the null distribution.

A real link is significant only if its TDS score exceeds the (1−α) percentile of this null.

#### `surrogate_tds(dataset, params=None, seed=None, verbose=False) → ndarray`

```python
# dataset shape: (n_subjects, T, N)
# Returns 1-D array of pooled null TDS scores
null_scores = surrogate_tds(dataset, params, seed=42)
```

#### `significance_threshold(null_scores, params=None) → float`

```python
threshold = significance_threshold(null_scores, params)
# e.g. ~7% for physiological networks (Bashan et al. 2012)
```

---

### `pyTDS.viz` — Visualization

#### `plot_tau_series(tau, t_vec, stbl_lbl=None, title="", ax=None)`

Plot τ₀(t) with stable points in red and unstable in blue.  
Reproduces Figure 1b style from Bashan et al. 2012.

#### `plot_tds_matrix(mat, labels=None, threshold=None, title="", ax=None, cmap="hot")`

Plot an N×N TDS matrix as a heatmap.  
Reproduces Figure 2 style from Bashan et al. 2012.

#### `plot_network(graph, tds_mat=None, labels=None, threshold=None, title="", ax=None)`

Plot a circular physiological network diagram.  
Edge thickness and colour represent TDS link strength.

---

## Complete Example — Physiological Network

```python
import numpy as np
import pyTDS
from pyTDS.params import TDSParams
from pyTDS.network import tds_matrix, fix_symmetry, apply_threshold, to_networkx
from pyTDS.surrogate import surrogate_tds, significance_threshold
from pyTDS.viz import plot_tds_matrix, plot_network

# signals shape: (T, N) — one column per physiological signal
signals = np.load("patient_signals.npy")
labels  = ["HR", "SpO2", "Thorax", "Abdomen", "EEG-delta"]

params = TDSParams(window=60, overlap=30, max_lag=30)

# Step 1: compute TDS matrix
mat = tds_matrix(signals, params, verbose=True)
mat = fix_symmetry(mat)

# Step 2: surrogate significance threshold
# dataset shape: (n_subjects, T, N)
dataset = np.load("all_subjects.npy")
null    = surrogate_tds(dataset, params, seed=0)
thresh  = significance_threshold(null, params)
print(f"Significance threshold: {thresh:.1f}%")

# Step 3: build network
adj = apply_threshold(mat, threshold=thresh)
G   = to_networkx(adj, tds_mat=mat, labels=labels)

# Step 4: visualize
plot_tds_matrix(mat, labels=labels, threshold=thresh)
plot_network(G, tds_mat=mat, labels=labels, threshold=thresh)
```

---

## Validation — Rulkov Oscillators

The `examples/rulkov_validation.py` script validates TDS on **coupled Rulkov map oscillators** — a synthetic chaotic neuron model where the ground-truth time delay is known exactly.

```bash
python examples/rulkov_validation.py
```

Results:

| Test | TDS Score | Result |
|------|-----------|--------|
| Coupled Rulkov (g=0.12, delay=16) | ~18% | Stable coupling detected at correct lag |
| Uncoupled Rulkov (g=0) | ~2.5% | Below threshold — no significant coupling |
| White noise | 0.0% | Exactly zero — true null |

The detected lag matches the true delay of 16 samples.

---

## Running Tests

```bash
# Terminal output only
python -m pytest tests/

# Full test run + Excel dashboard saved to "tests reports/"
python run_dashboard.py
```

The dashboard saves a timestamped `.xlsx` file showing every test — green for passed, red for failed with the failure reason.

---

## Project Structure

```
pyTDS/
├── pyTDS/
│   ├── __init__.py       # public API exports
│   ├── params.py         # TDSParams dataclass
│   ├── core.py           # core algorithm (steps 1–3)
│   ├── utils.py          # pbc_xcorr(), zscore()
│   ├── network.py        # tds_matrix(), threshold, networkx
│   ├── surrogate.py      # null distribution via cross-subject mixing
│   └── viz.py            # plotting helpers
├── tests/
│   ├── test_core.py      # 14 tests — utils + core algorithm
│   ├── test_network.py   # 19 tests — network module
│   └── test_surrogate.py # 13 tests — surrogate + significance
├── examples/
│   ├── rulkov_validation.py   # synthetic ground-truth validation
│   └── two_signals.py         # minimal two-signal demo
├── data/
│   ├── data-delta.txt    # benchmark EEG signal (delta band)
│   └── data-sigma.txt    # benchmark EEG signal (sigma band)
├── run_dashboard.py      # run all tests → Excel report
└── pyproject.toml
```

---

## References

- Bashan, A. et al. (2012). *Network physiology reveals relations between network topology and physiological function.* Nature Communications 3, 702.
- Bartsch, R.P. et al. (2015). *Phase transitions in physiologic coupling.* PNAS 109, 10181.
- Rulkov, N.F. (2001). *Regularization of synchronized chaotic bursts.* Physical Review Letters 86, 183.
- Markuza, A. (2024). *MSc Thesis — TDS applied to polysomnography data for sleep apnea classification.*
