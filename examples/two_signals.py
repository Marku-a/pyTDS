"""
Minimal two-signal TDS demo
============================
Shows the simplest possible usage of pyTDS on two synthetic signals
with a known time delay.

Scenario
--------
  s1 = noisy sine wave
  s2 = s1 delayed by 8 samples + extra noise

Expected result:
  - TDS score clearly above 0% (signals are coupled)
  - Dominant detected lag near -8 (s1 leads s2)
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import matplotlib.pyplot as plt
import pyTDS
from pyTDS.params import TDSParams
from pyTDS.viz import plot_tau_series

# ── Generate synthetic signals ─────────────────────────────────────────────────
rng = np.random.default_rng(0)

T = 2000
t = np.arange(T)

TRUE_DELAY = 8   # s1 leads s2 by 8 samples → expected detected lag = -8

s1 = np.sin(2 * np.pi * t / 60) + 0.3 * rng.standard_normal(T)
s2 = np.roll(s1, TRUE_DELAY) + 0.3 * rng.standard_normal(T)  # delay s2 by 8

# ── Run TDS ────────────────────────────────────────────────────────────────────
params = TDSParams(window=60, overlap=30, max_lag=30)
result = pyTDS.tds(s1, s2, params)

score        = result["score"]
tau          = result["tau"]
t_vec        = result["t_vec"]
stbl_lbl     = result["stbl_lbl"]
stable_taus  = result["stable_taus"]

detected_lag = int(np.round(np.median(stable_taus))) if len(stable_taus) > 0 else None

print("=" * 50)
print("  pyTDS — Two-Signal Minimal Demo")
print("=" * 50)
print(f"  True delay      : {TRUE_DELAY}  (s1 leads s2)")
print(f"  Expected lag    : -{TRUE_DELAY}  (negative = s1 leads)")
print(f"  TDS score       : {score:.1f}%")
print(f"  Detected lag    : {detected_lag}  (median of stable taus)")
print("=" * 50)

# ── Plot ───────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 4))
fig.suptitle("pyTDS — Two-Signal Demo", fontsize=13, fontweight='bold')

# Raw signals
axes[0].plot(t[:200], s1[:200], color='steelblue', lw=1.0, label='s1')
axes[0].plot(t[:200], s2[:200], color='crimson',   lw=1.0, alpha=0.8, label=f's2  (delayed by {TRUE_DELAY})')
axes[0].set_title("Signals (first 200 samples)")
axes[0].set_xlabel("Time (samples)")
axes[0].legend(fontsize=9)
axes[0].grid(alpha=0.3)

# tau(t) series
plot_tau_series(
    tau, t_vec, stbl_lbl,
    title=f"Time Delay tau(t) — TDS = {score:.1f}%",
    ax=axes[1]
)
axes[1].axhline(-TRUE_DELAY, color='green', linestyle='--', lw=1.5,
                label=f'True lag = -{TRUE_DELAY}')
axes[1].legend(fontsize=9)

plt.tight_layout()
out_path = os.path.join(os.path.dirname(__file__), "two_signals.png")
plt.savefig(out_path, dpi=150, bbox_inches='tight')
print(f"\nPlot saved to: {out_path}")
plt.show()
