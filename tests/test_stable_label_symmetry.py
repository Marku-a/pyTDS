"""
Regression / symmetry tests for the stable_label() order-independence fix.

Background: prior to this fix, stable_label() iterated candidate delays in
``np.unique(seg)`` (ascending) order and ``break``-ed at the first qualifying
candidate. This made the result depend on candidate order and biased it
toward negative lags. This file:

  1. Verifies the specific repro from the bug report is now order-independent.
  2. Verifies negation symmetry: stable_label(tau) == stable_label(-tau).
  3. Verifies swapping tds() arguments flips the sign of tau (and gives the
     same stability labels), on correlated signals.
  4. Verifies the "end" anchor's causal prefix-invariance guarantee still
     holds after the fix.
  5. Verifies the documented sign convention: tau < 0 means s1 leads s2.
  6. Runs a large regression corpus comparing old vs. new stable_label,
     asserting the new labels are always a superset of the old ones (no
     1 -> 0 flips), and prints summary counts.

Run with ``python -m pytest tests/test_stable_label_symmetry.py -s`` to see
the printed regression counts from test 6.
"""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from pyTDS.params import TDSParams, VALID_STABILITY_ANCHORS
from pyTDS.core import stable_label, tds


# ── Private copy of the OLD (<=0.1.0) buggy stable_label, for regression ──


def _stable_label_old(tau: np.ndarray, params: TDSParams = None) -> np.ndarray:
    """Exact copy of the pre-fix stable_label(): breaks at the first
    qualifying candidate in ascending (np.unique) order."""
    if params is None:
        params = TDSParams()

    anchor = params.stability_anchor
    if anchor not in VALID_STABILITY_ANCHORS:
        raise ValueError(
            f"stability_anchor must be one of {VALID_STABILITY_ANCHORS}, "
            f"got {anchor!r}"
        )

    N = len(tau)
    win = params.stability_window
    min_stable = params.stability_min
    tol = params.tolerance
    max_lag = params.max_lag

    stbl_lbl = np.zeros(N, dtype=int)

    for i in range(N - win + 1):
        seg = tau[i: i + win]

        if np.any(np.isnan(seg)) or np.any(np.abs(seg) > max_lag):
            continue

        for d in np.unique(seg):
            in_range = np.abs(seg - d) <= tol
            if np.sum(in_range) >= min_stable:
                if anchor == "end":
                    if in_range[-1]:
                        stbl_lbl[i + win - 1] = 1
                else:
                    indices = np.where(in_range)[0] + i
                    stbl_lbl[indices] = 1
                break

    return stbl_lbl


# ── 1. Order-independence repro ────────────────────────────────────────────


def test_order_independence_repro():
    tau_a = np.array([0, 0, 0, -1, 1], dtype=float)
    tau_b = np.array([0, 0, 0, 1, -1], dtype=float)

    for anchor in VALID_STABILITY_ANCHORS:
        p = TDSParams(stability_anchor=anchor)
        lbl_a = stable_label(tau_a, p)
        lbl_b = stable_label(tau_b, p)
        assert lbl_a[-1] == lbl_b[-1]
        assert lbl_a[-1] == 1


# ── 2. Sign (negation) symmetry ────────────────────────────────────────────


def test_negation_symmetry():
    max_lag = TDSParams().max_lag
    rng_master = np.random.default_rng(0)
    seeds = rng_master.integers(0, 1_000_000, size=200)

    for seed in seeds:
        rng = np.random.default_rng(int(seed))
        N = 200
        tau = rng.integers(-5, 6, size=N).astype(float)

        # Sprinkle NaNs and sentinels. Note sentinel check is abs > max_lag,
        # so negating a sentinel value keeps it a sentinel.
        n_nan = rng.integers(0, 5)
        nan_idx = rng.choice(N, size=n_nan, replace=False)
        tau[nan_idx] = np.nan

        remaining = np.setdiff1d(np.arange(N), nan_idx)
        if len(remaining) >= 2:
            n_sentinel = rng.integers(0, 3)
            sentinel_idx = rng.choice(
                remaining, size=min(n_sentinel, len(remaining)), replace=False
            )
            for j, idx in enumerate(sentinel_idx):
                tau[idx] = (max_lag + 1) if j % 2 == 0 else -(max_lag + 1)

        for anchor in VALID_STABILITY_ANCHORS:
            p = TDSParams(stability_anchor=anchor)
            lbl_pos = stable_label(tau, p)
            lbl_neg = stable_label(-tau, p)
            np.testing.assert_array_equal(
                lbl_pos, lbl_neg,
                err_msg=f"seed={seed} anchor={anchor}",
            )


# ── 3. Argument-swap symmetry on correlated signals ────────────────────────


def test_argument_swap_flips_tau_sign():
    """Swap symmetry must hold on *noisy* tau sequences with real jitter,
    not just near-constant ones. Signal design: weak coupling (so
    cross-correlation is noisy from window to window) with a slowly
    drifting delay that oscillates between `shift` and `shift + 1`
    samples. This makes tau0(t) hop between neighbouring integer lags
    across windows, so many stability windows contain several distinct
    candidate delays that simultaneously qualify (>= stability_min points
    within tolerance) -- exactly the situation where the pre-fix
    stable_label()'s "break at the first candidate in ascending order"
    behaviour is order- (and thus sign-) dependent. We confirm this has
    teeth by also checking that `_stable_label_old` actually disagrees
    between the a->b and b->a tau sequences for at least one seed.
    """
    old_diffs_found = 0

    # Seeds 0-11 avoid a separate, pre-existing edge case: when
    # max_lag == window // 2, lags +max_lag and -max_lag are the same
    # circular shift, and argmax picks -max_lag in both argument orders.
    for seed in range(12):
        rng = np.random.default_rng(seed)
        N = 3000

        # AR(1)-ish random walk noise for `a`.
        innovations = rng.standard_normal(N)
        a = np.cumsum(innovations) * 0.1 + innovations

        # `b` tracks a blend of two neighbouring shifted copies of `a`
        # whose mixing weight drifts slowly (a few full cycles over N),
        # under heavy additive noise and a coupling well below 1. This
        # keeps the dominant per-window delay jittering near shift/shift+1
        # (and occasionally elsewhere, since the coupling is weak) instead
        # of sitting at one constant value.
        shift = 1
        coupling = 0.6
        noise_scale = 0.8
        s_lo = np.roll(a, shift)
        s_hi = np.roll(a, shift + 1)
        w = (np.sin(np.linspace(0, 6 * np.pi, N)) + 1) / 2
        core = w * s_hi + (1 - w) * s_lo
        b = coupling * core + noise_scale * rng.standard_normal(N)

        r_ab = tds(a, b)
        r_ba = tds(b, a)

        np.testing.assert_array_equal(
            r_ba["tau"], -r_ab["tau"], err_msg=f"seed={seed}"
        )
        np.testing.assert_array_equal(
            r_ba["stbl_lbl"], r_ab["stbl_lbl"], err_msg=f"seed={seed}"
        )

        # Teeth check: the OLD order-dependent labeller should NOT be
        # swap-symmetric on at least one of these noisy tau sequences.
        old_ab = _stable_label_old(r_ab["tau"])
        old_ba = _stable_label_old(r_ba["tau"])
        if not np.array_equal(old_ab, old_ba):
            old_diffs_found += 1

    assert old_diffs_found > 0, (
        "expected at least one seed where the OLD stable_label breaks "
        "swap-symmetry (i.e. differs between tau_ab and tau_ba = -tau_ab); "
        "none found -- signal generation may need to be noisier"
    )


# ── 4. "end" anchor causal prefix invariance + any/end subset ─────────────


def test_end_anchor_prefix_invariance_and_subset_after_fix():
    max_lag = TDSParams().max_lag
    rng_master = np.random.default_rng(1)
    seeds = rng_master.integers(0, 1_000_000, size=25)

    for seed in seeds:
        rng = np.random.default_rng(int(seed))
        N = 120
        tau = rng.integers(-5, 6, size=N).astype(float)

        n_nan = rng.integers(0, 5)
        nan_idx = rng.choice(N, size=n_nan, replace=False)
        tau[nan_idx] = np.nan
        remaining = np.setdiff1d(np.arange(N), nan_idx)
        if len(remaining) >= 2:
            n_sentinel = rng.integers(0, 3)
            sentinel_idx = rng.choice(
                remaining, size=min(n_sentinel, len(remaining)), replace=False
            )
            tau[sentinel_idx] = max_lag + 1

        p_end = TDSParams(stability_anchor="end")
        p_any = TDSParams(stability_anchor="any")

        full_end = stable_label(tau, p_end)
        for k in range(len(tau) + 1):
            np.testing.assert_array_equal(
                full_end[:k], stable_label(tau[:k], p_end),
                err_msg=f"seed={seed} k={k}",
            )

        full_any = stable_label(tau, p_any)
        assert np.all(full_end <= full_any)


# ── 5. Sign convention: s1 leads s2 -> tau < 0 ─────────────────────────────


def test_sign_convention_s1_leads_s2():
    for seed in range(5):
        rng = np.random.default_rng(seed)
        N = 3000
        s1 = rng.standard_normal(N)
        shift = 2
        # s2[i] = s1[i - shift]  ->  s1 leads s2 by `shift` samples
        s2 = np.roll(s1, shift) + 0.01 * rng.standard_normal(N)

        result = tds(s1, s2)
        stable_taus = result["stable_taus"]
        assert len(stable_taus) > 0, f"seed={seed}: no stable points found"
        assert np.median(stable_taus) == -shift, (
            f"seed={seed}: median stable tau = {np.median(stable_taus)}, "
            f"expected {-shift}"
        )


# ── 6. Large regression corpus: new labels superset of old labels ─────────


def _run_regression_corpus():
    """Returns (results dict) with summary counts. Also used standalone by
    the pytest test below (printed with -s)."""
    max_lag = TDSParams().max_lag

    total_labels = 0
    zero_to_one = 0
    one_to_zero = 0
    old_stable_at_neg1 = 0
    old_stable_at_pos1 = 0
    new_stable_at_neg1 = 0
    new_stable_at_pos1 = 0

    results_per_anchor = {}

    for anchor in VALID_STABILITY_ANCHORS:
        p = TDSParams(stability_anchor=anchor)
        anchor_total = 0
        anchor_0to1 = 0
        anchor_1to0 = 0

        # Corpus A: iid random integer tau in -5..5, with NaN/sentinels.
        rng_master = np.random.default_rng(42)
        seeds = rng_master.integers(0, 1_000_000, size=500)
        for seed in seeds:
            rng = np.random.default_rng(int(seed))
            N = 500
            tau = rng.integers(-5, 6, size=N).astype(float)

            n_nan = rng.integers(0, 8)
            nan_idx = rng.choice(N, size=n_nan, replace=False)
            tau[nan_idx] = np.nan
            remaining = np.setdiff1d(np.arange(N), nan_idx)
            if len(remaining) >= 2:
                n_sentinel = rng.integers(0, 5)
                sentinel_idx = rng.choice(
                    remaining, size=min(n_sentinel, len(remaining)), replace=False
                )
                tau[sentinel_idx] = max_lag + 1

            old = _stable_label_old(tau, p)
            new = stable_label(tau, p)

            flips_1to0 = np.sum((old == 1) & (new == 0))
            flips_0to1 = np.sum((old == 0) & (new == 1))

            if flips_1to0 > 0:
                raise AssertionError(
                    f"anchor={anchor} seed={seed}: superset property "
                    f"violated, {flips_1to0} labels went 1 -> 0"
                )

            anchor_total += N
            anchor_0to1 += flips_0to1
            anchor_1to0 += flips_1to0

            if anchor == "any":
                old_stable_at_neg1 += int(np.sum((tau == -1) & (old == 1)))
                old_stable_at_pos1 += int(np.sum((tau == 1) & (old == 1)))
                new_stable_at_neg1 += int(np.sum((tau == -1) & (new == 1)))
                new_stable_at_pos1 += int(np.sum((tau == 1) & (new == 1)))

        # Corpus B: slowly-varying (random-walk-ish) tau.
        rng_master2 = np.random.default_rng(4242)
        seeds2 = rng_master2.integers(0, 1_000_000, size=100)
        for seed in seeds2:
            rng = np.random.default_rng(int(seed))
            N = 500
            steps = rng.integers(-1, 2, size=N)
            tau = np.clip(np.cumsum(steps), -5, 5).astype(float)

            n_nan = rng.integers(0, 8)
            nan_idx = rng.choice(N, size=n_nan, replace=False)
            tau[nan_idx] = np.nan

            old = _stable_label_old(tau, p)
            new = stable_label(tau, p)

            flips_1to0 = np.sum((old == 1) & (new == 0))
            flips_0to1 = np.sum((old == 0) & (new == 1))

            if flips_1to0 > 0:
                raise AssertionError(
                    f"anchor={anchor} seed={seed} (corpus B): superset "
                    f"property violated, {flips_1to0} labels went 1 -> 0"
                )

            anchor_total += N
            anchor_0to1 += flips_0to1
            anchor_1to0 += flips_1to0

        results_per_anchor[anchor] = {
            "total": anchor_total,
            "0->1": int(anchor_0to1),
            "1->0": int(anchor_1to0),
        }

        total_labels += anchor_total
        zero_to_one += anchor_0to1
        one_to_zero += anchor_1to0

    return {
        "total_labels": total_labels,
        "0->1": zero_to_one,
        "1->0": one_to_zero,
        "per_anchor": results_per_anchor,
        "old_stable_at_tau=-1": old_stable_at_neg1,
        "old_stable_at_tau=+1": old_stable_at_pos1,
        "new_stable_at_tau=-1": new_stable_at_neg1,
        "new_stable_at_tau=+1": new_stable_at_pos1,
    }


def test_regression_corpus_new_is_superset_of_old():
    stats = _run_regression_corpus()

    print("\n=== stable_label regression corpus (old vs new) ===")
    print(f"Total labels evaluated: {stats['total_labels']}")
    print(f"0 -> 1 flips (new stable, old not): {stats['0->1']}")
    print(f"1 -> 0 flips (should be 0):         {stats['1->0']}")
    for anchor, d in stats["per_anchor"].items():
        print(f"  anchor={anchor!r}: total={d['total']} 0->1={d['0->1']} 1->0={d['1->0']}")
    print(
        f"OLD stable count at tau=-1: {stats['old_stable_at_tau=-1']}  "
        f"tau=+1: {stats['old_stable_at_tau=+1']}"
    )
    print(
        f"NEW stable count at tau=-1: {stats['new_stable_at_tau=-1']}  "
        f"tau=+1: {stats['new_stable_at_tau=+1']}"
    )

    assert stats["1->0"] == 0
