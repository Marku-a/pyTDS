# Changelog

All notable changes to this project are documented in this file.

## [0.2.1] - 2026-09-25

### Fixed

- **`pbc_xcorr()` aliased lag at `max_lag >= window / 2` (the defaults).**
  Under periodic boundary conditions a circular shift by `k` is
  indistinguishable from a shift by `k - N`, so `pbc_xcorr()` was
  returning the same aliased `±N/2` circular shift twice (once as `+max_lag`
  and once as `-max_lag`) whenever `max_lag >= N/2`, which is true for the
  default parameters (`window=60`, `max_lag=30`). `time_delay_interaction()`
  then picked the peak with `np.argmax(np.abs(C))`, which ties on these
  duplicate values and always resolves to the first (most negative) one
  regardless of argument order, breaking the documented antisymmetry
  `tds(b, a)['tau'] == -tds(a, b)['tau']`. Lags returned by `pbc_xcorr()`
  are now limited to `|k| <= (window - 1) // 2`, so with the default
  parameters `tau = ±30` can no longer be reported (effective max lag is
  now 29). Results may change in windows whose cross-correlation peak was
  at `±30`.

## [0.2.0] - 2026-09-25

### Fixed

- **`stable_label()` order-dependence.** The stability check iterated
  candidate delays in ascending (`np.unique`) order and stopped at the
  first candidate that qualified (>= `stability_min` points within
  `tolerance`), instead of evaluating every candidate in the window. This
  made the result depend on the arbitrary order of `np.unique`, and
  because that order is ascending, it biased labeling toward
  **negative-lag** candidates — e.g. `tau=[0,0,0,-1,1]` and
  `tau=[0,0,0,1,-1]` produced different labels for the same last point.
  `stable_label()` now evaluates every candidate delay in each window and
  is order-independent: a point is stable iff *some* window containing it
  has a qualifying candidate covering it. This bug (the unconditional
  `break`) was present since the initial port from the original
  implementation; no MATLAB reference source is available in this
  repository to confirm whether it originated there, and the published
  definition (Bashan et al., Nature Communications 2012) does not specify
  a candidate evaluation order.

### Documented

- **Sign convention of `τ` (no computation change).** An inline comment in
  `pbc_xcorr()` had the cross-correlation formula backwards (it read
  `seg1[i] * seg2[(i + k) % N]`, the code actually computes
  `seg1[(i + k) % N] * seg2[i]`). This was a documentation-only bug — the
  underlying FFT computation was never changed. `pbc_xcorr()`,
  `time_delay_interaction()`, and `tds()` now explicitly document: `τ < 0`
  means `s1` leads `s2`; `τ > 0` means `s2` leads `s1`; swapping the two
  input signals flips the sign of `τ`.

### Upgrade notes

- **TDS scores, stable labels, and network/surrogate significance
  thresholds computed with pyTDS <= 0.1.0 may change** after upgrading.
  Because the new `stable_label()` rule is a superset of the old
  one (see the regression tests in
  `tests/test_stable_label_symmetry.py`), **scores can only increase, never
  decrease**, relative to the old (buggy) behaviour.
- The `τ` sign convention is now correctly documented. Previously an inline comment in
  `pbc_xcorr()` stated the opposite convention (while the actual
  numerical behaviour was unchanged). Any downstream interpretation of
  lead/lag relationships based on the old (incorrect) docstring wording
  should be re-checked against the actual sign convention now documented
  in `pbc_xcorr()`, `time_delay_interaction()`, and `tds()`.

## [0.1.0]

Initial versioned release.
