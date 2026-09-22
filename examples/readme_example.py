"""
README example — physiological network from synthetic polysomnography signals.

Generates a realistic 5-signal dataset (1 Hz, 1 hour) and runs the full
TDS pipeline: matrix → surrogate threshold → network → plots.

Signals modelled after Bashan et al. (2012, Nature Communications):

  Thorax   — frequency-modulated breathing oscillator (0.2–0.3 Hz), variable
              amplitude; the primary coupling driver for all other signals.
  Abdomen  — same breathing rhythm, slightly different amplitude and phase offset
              (~2–4 s lead).  Thorax–Abdomen coupling is the strongest link.
  HR       — respiratory sinus arrhythmia (RSA, ~5 s lag to Thorax) +
              Mayer-wave sympathetic oscillation (~0.1 Hz) + slow drift.
              Realistic range ~50–90 bpm.
  SpO2     — baseline ~97 %, small amplitude (~0.5 %), lags Thorax by ~15 s
              (lung ventilation → blood oxygenation transport delay).
  EEG-delta— delta-band power envelope (slow, <0.1 Hz), weakly coupled to
              breathing rhythm + coloured-noise baseline (AR(1)).

Usage
-----
    python examples/readme_example.py
"""

import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pyTDS
from pyTDS.params import TDSParams
from pyTDS.network import tds_matrix, fix_symmetry, apply_threshold, to_networkx
from pyTDS.surrogate import surrogate_tds, significance_threshold
from pyTDS.viz import plot_tds_matrix, plot_network
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))

# ── Progress helpers ───────────────────────────────────────────────────────────

def step(label, total=5):
    step.n = getattr(step, "n", 0) + 1
    step.t = time.perf_counter()
    print(f"\n[{step.n}/{total}] {label} ...", flush=True)

def done():
    print(f"      done  ({time.perf_counter() - step.t:.1f}s)", flush=True)

def progress_bar(current, total, label=""):
    bar_len = 30
    filled  = int(bar_len * current / total)
    bar     = "#" * filled + "-" * (bar_len - filled)
    pct     = current / total * 100
    print(f"\r  [{bar}] {pct:5.1f}%  ({current}/{total}{(' ' + label) if label else ''})",
          end="", flush=True)

# ── Synthetic data generation ──────────────────────────────────────────────────

T = 3600   # 1 hour at 1 Hz
t = np.arange(T)


def make_patient(seed):
    """Generate one synthetic polysomnography patient (T=3600, N=5)."""
    rng = np.random.default_rng(seed)

    def zscore(x):
        return (x - x.mean()) / (x.std() + 1e-9)

    # Thorax: frequency-modulated breathing oscillator
    base_freq  = rng.uniform(0.20, 0.28)           # 12–17 breaths/min
    freq_drift = np.cumsum(rng.standard_normal(T) * 0.00008)
    freq       = np.clip(base_freq + freq_drift, 0.15, 0.35)
    phase      = np.cumsum(2 * np.pi * freq)
    amp_env    = 1.0 + 0.25 * np.sin(2 * np.pi * 0.004 * t + rng.uniform(0, 2*np.pi))
    amp_env   += 0.10 * rng.standard_normal(T)
    amp_env    = np.clip(amp_env, 0.4, 1.8)
    thorax     = amp_env * np.sin(phase) + 0.08 * rng.standard_normal(T)

    # Abdomen: same rhythm, per-patient amplitude ratio and phase offset
    abd_phase_offset = rng.uniform(0.15, 0.45)     # ~2–4 s lead vs Thorax
    abd_amp_ratio    = rng.uniform(0.70, 1.00)
    abdomen = abd_amp_ratio * amp_env * np.sin(phase + abd_phase_offset)
    abdomen += 0.10 * rng.standard_normal(T)

    # HR: RSA + Mayer waves + slow circadian drift
    rsa_delay  = rng.integers(4, 8)                # 4–7 s
    rsa_amp    = rng.uniform(2.5, 4.5)             # ±3–5 bpm
    mayer_freq = rng.uniform(0.08, 0.12)           # ~0.1 Hz
    mayer_amp  = rng.uniform(1.5, 3.0)
    slow_drift = 4.0 * np.sin(2 * np.pi * 0.002 * t + rng.uniform(0, 2*np.pi))
    hr_mean    = rng.uniform(58, 72)
    hr = (hr_mean
          + rsa_amp   * np.roll(thorax, rsa_delay)
          + mayer_amp * np.sin(2 * np.pi * mayer_freq * t + rng.uniform(0, 2*np.pi))
          + slow_drift
          + 0.4 * rng.standard_normal(T))
    hr = np.clip(hr, 40, 105)

    # SpO2: stable baseline, small breathing-coupled dips
    spo2_delay = rng.integers(12, 18)              # ~15 s lung-to-pulse lag
    spo2_mean  = rng.uniform(96.0, 98.5)
    spo2_amp   = rng.uniform(0.20, 0.50)
    spo2_slow  = 0.4 * np.sin(2 * np.pi * 0.003 * t + rng.uniform(0, 2*np.pi))
    spo2 = (spo2_mean
            - spo2_amp * np.roll(thorax, spo2_delay)
            + spo2_slow
            + 0.08 * rng.standard_normal(T))
    spo2 = np.clip(spo2, 88, 100)

    # EEG-delta: AR(1) coloured noise + weak breathing coupling
    eeg_delay = rng.integers(6, 12)
    ar_noise  = np.zeros(T)
    ar_coef   = rng.uniform(0.92, 0.97)
    for i in range(1, T):
        ar_noise[i] = ar_coef * ar_noise[i - 1] + rng.standard_normal()
    eeg_slow  = 0.6 * np.sin(2 * np.pi * 0.018 * t + rng.uniform(0, 2*np.pi))
    eeg_delta = (0.35 * np.roll(thorax, eeg_delay)
                 + eeg_slow
                 + 0.3 * ar_noise
                 + 0.15 * rng.standard_normal(T))

    return np.column_stack([
        zscore(hr),
        zscore(spo2),
        zscore(thorax),
        zscore(abdomen),
        zscore(eeg_delta),
    ])


# ── Generate data ──────────────────────────────────────────────────────────────

step("Generating synthetic data")

patient_signals = make_patient(seed=42)
np.save(os.path.join(HERE, "patient_signals.npy"), patient_signals)

n_subjects = 30
subjects   = []
for s in range(n_subjects):
    progress_bar(s + 1, n_subjects, "subjects")
    subjects.append(make_patient(seed=s))
print()

all_subjects = np.stack(subjects)
np.save(os.path.join(HERE, "all_subjects.npy"), all_subjects)
done()

# ── README example ─────────────────────────────────────────────────────────────

signals = np.load(os.path.join(HERE, "patient_signals.npy"))
labels  = ["HR", "SpO2", "Thorax", "Abdomen", "EEG-delta"]
params  = TDSParams(window=60, overlap=30, max_lag=30)

# Step 1: compute TDS matrix
step("Computing TDS matrix")
mat = tds_matrix(signals, params, verbose=True)
mat = fix_symmetry(mat)
done()

# Step 2: surrogate significance threshold
step("Surrogate significance threshold")
dataset = np.load(os.path.join(HERE, "all_subjects.npy"))
null    = surrogate_tds(dataset, params, seed=0, verbose=True)
thresh  = significance_threshold(null, params)
done()
print(f"      threshold = {thresh:.1f}%")

# Step 3: build network
step("Building network")
adj = apply_threshold(mat, threshold=thresh)
G   = to_networkx(adj, tds_mat=mat, labels=labels)
done()

# Step 4: visualize
step("Saving plots")
fig_mat, _ = plot_tds_matrix(mat, labels=labels, threshold=thresh)
fig_net, _ = plot_network(G, tds_mat=mat, labels=labels, threshold=thresh)
fig_mat.savefig(os.path.join(HERE, "readme_tds_matrix.png"), dpi=150, bbox_inches="tight")
fig_net.savefig(os.path.join(HERE, "readme_tds_network.png"), dpi=150, bbox_inches="tight")
done()
print("      readme_tds_matrix.png, readme_tds_network.png")

plt.show()
