"""
Rulkov Oscillator Validation
============================
Validates the pyTDS core algorithm on synthetic coupled chaotic oscillators
where the ground-truth time delay is known exactly.

The Rulkov map (Rulkov 2001) is a 2D discrete-time model of a spiking/bursting
neuron. Two Rulkov oscillators are coupled with a known delay `tau_true`.
TDS should:
  1. Detect stable coupling at approximately tau_true
  2. Give a high TDS score (strongly coupled)
  3. Give a near-zero TDS score when the oscillators are uncoupled (g=0)

Parameters match the thesis:
  n=6000, mu=0.001, alpha=[4.5, 4.1], g=0.12, tau_true=16

Reference:
  Asaf Markuza MSc thesis, section "TDS on Rulkov Oscillator"
  Rulkov, N.F. (2001). Regularization of synchronized chaotic bursts.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import matplotlib.pyplot as plt
import pyTDS
from pyTDS.params import TDSParams


# ── Rulkov map ─────────────────────────────────────────────────────────────────

def rulkov_f(x, y, alpha):
    """Rulkov map fast variable update."""
    if x <= 0:
        return alpha / (1.0 - x) + y
    elif x < alpha + y:
        return alpha + y
    else:
        return -1.0


def coupled_rulkov_with_delay(n, mu, alpha, sigma, beta, sigma_e, g, delay):
    """
    Generate two coupled Rulkov map oscillators with a known time delay.

    Parameters
    ----------
    n : int         Number of time steps (extra 50 discarded as transient)
    mu : float      Slow variable time scale
    alpha : list    [alpha1, alpha2] — nonlinearity parameters
    sigma : list    [sigma1, sigma2] — DC input
    beta : float    Coupling coefficient
    sigma_e : float Noise coupling coefficient
    g : float       Coupling strength
    delay : int     Coupling delay in time steps (ground truth tau)

    Returns
    -------
    x1, x2 : np.ndarray  Fast variable time series (the "spikes")
    """
    n_total = n + 50   # add transient
    x1 = np.zeros(n_total)
    y1 = np.zeros(n_total)
    x2 = np.zeros(n_total)
    y2 = np.zeros(n_total)

    x1[0] = -1.0
    x2[0] = -1.0
    y1[0] = -4.0
    y2[0] = -4.0

    # Two-pass: first pass to find stable initial conditions
    for _ in range(2):
        for i in range(1, n_total):
            # Coupling from osc2 → osc1 (with delay)
            if delay >= i:
                b1 = 0.0
                s1 = 0.0
            else:
                b1 = g * beta * (x2[i - 1 - delay] - x1[i - 1])
                s1 = g * sigma_e * (x2[i - 1 - delay] - x1[i - 1])

            x1[i] = rulkov_f(x1[i - 1], y1[i - 1] + b1, alpha[0])
            y1[i] = y1[i - 1] - mu * (x1[i - 1] + 1) + mu * sigma[0] + mu * s1

            # Coupling from osc1 → osc2 (with delay)
            if delay >= i:
                b2 = 0.0
                s2 = 0.0
            else:
                b2 = g * beta * (x1[i - 1 - delay] - x2[i - 1])
                s2 = g * sigma_e * (x1[i - 1 - delay] - x2[i - 1])

            x2[i] = rulkov_f(x2[i - 1], y2[i - 1] + b2, alpha[1])
            y2[i] = y2[i - 1] - mu * (x2[i - 1] + 1) + mu * sigma[1] + mu * s2

        # Second pass: reset with median initial conditions
        y1[0] = np.median(y1)
        y2[0] = np.median(y2)
        x1[0] = -1.0
        x2[0] = -1.0

    # Discard transient
    return x1[50:], x2[50:]


# ── Main validation ─────────────────────────────────────────────────────────────

def run_validation():
    print("=" * 55)
    print("  pyTDS — Rulkov Oscillator Validation")
    print("=" * 55)

    # Parameters from Asaf's thesis
    n = 6000
    mu = 0.001
    alpha = [4.5, 4.1]
    sigma = [0.01, -0.01]
    beta = 1.0
    sigma_e = 1.0
    g_coupled = 0.12
    tau_true = 16

    params = TDSParams(window=60, overlap=30, max_lag=30)

    # ── Test 1: Coupled oscillators ──────────────────────────────────────────
    print(f"\n[1] Coupled  (g={g_coupled}, delay={tau_true})")
    x1, x2 = coupled_rulkov_with_delay(
        n, mu, alpha, sigma, beta, sigma_e, g_coupled, tau_true
    )
    result_coupled = pyTDS.tds(x1, x2, params)
    score_c = result_coupled["score"]
    tau_stable = result_coupled["stable_taus"]
    detected_delay = int(np.round(np.median(tau_stable))) if len(tau_stable) > 0 else None

    print(f"   TDS score    : {score_c:.1f}%")
    print(f"   True delay   : {tau_true}")
    print(f"   Detected lag : {detected_delay}  (median of stable taus)")

    # ── Test 2: Uncoupled oscillators (g=0) ──────────────────────────────────
    # Same Rulkov map, no coupling. Both oscillators are structurally similar
    # (same alpha range, same burst rhythm) so they can accidentally share a
    # lag for a few windows — gives a small non-zero score (~2-3%).
    # This is NOT a bug; it shows TDS is not fooled much by similar-looking
    # signals, but it also shows the score won't be exactly 0 in this case.
    print(f"\n[2] Uncoupled Rulkov (g=0)")
    x1r, x2r = coupled_rulkov_with_delay(
        n, mu, alpha, sigma, beta, sigma_e, g=0.0, delay=tau_true
    )
    result_rulkov_unc = pyTDS.tds(x1r, x2r, params)
    score_r = result_rulkov_unc["score"]
    print(f"   TDS score    : {score_r:.1f}%  (small but non-zero — same oscillator type)")

    # ── Test 3: White noise (true null) ──────────────────────────────────────
    # Pure white noise has zero structure. The lag jumps to a completely
    # random value each window → stability criterion never met → TDS = 0.0%.
    print(f"\n[3] White noise (true null)")
    rng = np.random.default_rng(0)
    x1n = rng.standard_normal(n)
    x2n = rng.standard_normal(n)
    result_noise = pyTDS.tds(x1n, x2n, params)
    score_n = result_noise["score"]
    print(f"   TDS score    : {score_n:.1f}%  (expected exactly 0.0%)")

    # ── Summary ──────────────────────────────────────────────────────────────
    print("\n" + "-" * 55)
    print("  VALIDATION SUMMARY")
    print("-" * 55)
    # Bidirectional coupling: either oscillator can lead → check |lag| ~ tau_true
    coupled_pass  = score_c > 10.0
    rulkov_unc_pass = score_r < 10.0   # below surrogate threshold (~7%)
    noise_pass    = score_n == 0.0
    delay_pass    = (
        detected_delay is not None
        and abs(abs(detected_delay) - tau_true) <= 3
    )

    print(f"  Coupled TDS > 10%            : {'PASS' if coupled_pass else 'FAIL'}  ({score_c:.1f}%)")
    print(f"  Uncoupled Rulkov TDS < 10%   : {'PASS' if rulkov_unc_pass else 'FAIL'}  ({score_r:.1f}%)")
    print(f"  White noise TDS = 0%         : {'PASS' if noise_pass else 'FAIL'}  ({score_n:.1f}%)")
    print(f"  |Detected lag| near {tau_true}   : {'PASS' if delay_pass else 'FAIL'}  ({detected_delay})")
    print("-" * 55)
    all_pass = coupled_pass and rulkov_unc_pass and noise_pass and delay_pass
    print(f"  Overall: {'ALL TESTS PASSED' if all_pass else 'SOME TESTS FAILED'}")
    print("=" * 55)

    # ── Plot ─────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(3, 2, figsize=(14, 12))
    fig.suptitle("pyTDS — Rulkov Oscillator Validation", fontsize=13, fontweight='bold')

    from pyTDS.viz import plot_tau_series
    t = np.arange(n)

    # Row 1: Coupled
    axes[0, 0].plot(t[:300], x1[:300], color='steelblue', lw=0.8, label='Osc 1')
    axes[0, 0].plot(t[:300], x2[:300], color='crimson',   lw=0.8, alpha=0.7, label='Osc 2')
    axes[0, 0].set_title(f"[1] Coupled signals (g={g_coupled}, delay={tau_true})")
    axes[0, 0].set_xlabel("Time step")
    axes[0, 0].legend(fontsize=8)
    axes[0, 0].grid(alpha=0.3)

    plot_tau_series(
        result_coupled["tau"], result_coupled["t_vec"],
        result_coupled["stbl_lbl"],
        title=f"[1] Coupled: tau(t) — TDS={score_c:.1f}%",
        ax=axes[0, 1]
    )

    # Row 2: Uncoupled Rulkov
    axes[1, 0].plot(t[:300], x1r[:300], color='steelblue', lw=0.8, label='Osc 1')
    axes[1, 0].plot(t[:300], x2r[:300], color='crimson',   lw=0.8, alpha=0.7, label='Osc 2')
    axes[1, 0].set_title("[2] Uncoupled Rulkov (g=0)")
    axes[1, 0].set_xlabel("Time step")
    axes[1, 0].legend(fontsize=8)
    axes[1, 0].grid(alpha=0.3)

    plot_tau_series(
        result_rulkov_unc["tau"], result_rulkov_unc["t_vec"],
        result_rulkov_unc["stbl_lbl"],
        title=f"[2] Uncoupled Rulkov: tau(t) — TDS={score_r:.1f}%",
        ax=axes[1, 1]
    )

    # Row 3: White noise
    axes[2, 0].plot(t[:300], x1n[:300], color='steelblue', lw=0.8, label='Noise 1')
    axes[2, 0].plot(t[:300], x2n[:300], color='crimson',   lw=0.8, alpha=0.7, label='Noise 2')
    axes[2, 0].set_title("[3] White noise (true null)")
    axes[2, 0].set_xlabel("Time step")
    axes[2, 0].legend(fontsize=8)
    axes[2, 0].grid(alpha=0.3)

    plot_tau_series(
        result_noise["tau"], result_noise["t_vec"],
        result_noise["stbl_lbl"],
        title=f"[3] White noise: tau(t) — TDS={score_n:.1f}%",
        ax=axes[2, 1]
    )

    plt.tight_layout()
    out_path = os.path.join(os.path.dirname(__file__), "rulkov_validation.png")
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f"\nPlot saved to: {out_path}")
    plt.show()

    return all_pass


if __name__ == "__main__":
    success = run_validation()
    sys.exit(0 if success else 1)
