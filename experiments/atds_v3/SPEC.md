# aTDS run-3 experiments (answer the independent reviewers)

All TDS computations use **this repository's pyTDS only** (`pyTDS.core.time_delay_interaction`, `stable_label`, `tds_score`, `pyTDS.params.TDSParams`). Never TDSpy, TDSTOCK or any other TDS code; do not modify `pyTDS/`. Non-TDS baselines (cross-correlation, GCC-PHAT) are implemented here in numpy. Fixed seeds, 4 worker processes. Private sleep data are read from a path given on the command line; only aggregates are committed.

Reuse:
- aTDS rule: `experiments/atds_real/atds21.py` (`window_params`, `calib`, `PUBLISHED`, `delay_of`).
- Slow-down benchmark: `atds_brief/experiments/exp4.py` `U1` and helpers (`tests.pair`, `run_methods`, `summarize`, `boot_diff90`, `WGRID`, `TGRID`). `exp4.py` expects a checkout of branch `claude/practical-mccarthy-n2mec6` (for `experiments/adaptive/` and `tests.py`/`exp3.py` imports): create one with `git worktree add <scratch>/atds-code origin/claude/practical-mccarthy-n2mec6` outside the repo, as in `atds_brief/README.md`. Benchmark: `tc = 2r`, `d = 5r`, coupling `c = 0.3`, records of `N = 20000`, aTDS set on 5 training pairs, 100 coupled + 100 uncoupled scored pairs per r, seeds as in U1.
- Cohort code: `experiments/atds_real/ncom_cohort.py` (loading, stage assignment, cross-subject null, effects).

Outputs go to `experiments/atds_v3/results/` (JSON + CSV summary per experiment) with a `README.md` describing each file.

## E1 — Baselines beyond TDS settings (slow-down benchmark, r in {1, 2, 4, 8, 16})

Score the same 100 + 100 test pairs per r (same seeds as U1) with:

1. `xcorr_full`: full-record normalised cross-correlation of the z-scored pair over lags `|tau| <= 1000`; score = max |C|; delay = argmax.
2. `gccphat_full`: full-record GCC-PHAT (cross-spectrum divided by its magnitude, small epsilon), same lag range; score = peak value; delay = argmax.
3. `rescaled_est`: the published setting rescaled by an **estimated** factor `r_hat = max(1, round(tc_hat / 2))`, where `tc_hat` is the median over the 5 training pairs (both signals) of the first lag at which the autocorrelation falls below 1/e. Params: window 60 r_hat, step 30 r_hat, lag 30 r_hat, tolerance r_hat. Scored with pyTDS.
4. For reference, recompute `published`, `atds` and `rescaled` (true r) exactly as U1 (or load from `atds_brief/results/exp4.json` U1 and check the recomputation matches for one r).

Metrics per method and r: AUC with bootstrap 95 % CI, delay found (|delay - true| <= max(1, 0.1 d); for TDS methods use `delay_of` / the U1 hit rule unchanged; for global methods the argmax), Wilson CI; paired AUC difference aTDS minus each baseline with 95 % CI (bootstrap as `boot_diff90`).

## E2 — Sensitivity of the design constants (slow-down benchmark, r in {1, 4, 8}, 100 + 100 pairs)

Vary one constant at a time around the default (n_eff 30, budget 5 %, z 2, floor 61, 20 rotations): n_eff in {15, 30, 45, 60}; budget in {1, 5, 10} %; z in {0, 2, 3}; floor in {31, 61, 121}; rotations in {10, 20, 40}. Rhythm threshold in {0.2, 0.3, 0.4} on the rhythm case of `exp4.U6` (AR(2), period about 40, delays 8 and 30; 30 + 30 pairs). Report window, tolerance, AUC (CI), delay found, mean chance score on uncoupled pairs.

## E3 — Harder synthetic benchmarks (40 + 40 pairs each unless stated; 20000 samples)

1. **Regime switching** (the sleep-stage case): blocks of 2000 samples alternating regime A (x memory tc 8; coupled y = 0.6 x(t-20) + noise) and regime B (tc 32; uncoupled, c = 0). Controls: same regimes with c = 0 throughout. Settings: published; aTDS (whole record); aTDS capped at 1000 (half a block); aTDS per regime (window and tolerance set separately from the regime-A and regime-B segments of the training pairs, then each regime scored with its own parameters). Score per regime with the 80 % purity rule of the cohort code. Metrics: AUC regime A coupled vs control (regime A), false-stability rate in regime B (mean score of regime-B windows in coupled pairs), delay found in regime A.
2. **Nonlinear coupling**, tc 8, d 20: (a) monotone y = 0.6 tanh(2 x(t-d)) / std + noise; (b) non-monotone y = 0.6 (x(t-d)^2 - 1)/sqrt(2) + noise, for which linear correlation is zero by construction. 100 + 100 pairs; published vs aTDS vs xcorr_full; AUC and delay found.
3. **Drifting delay**: tc 8, c 0.6, delay drifting linearly from 10 to 30 samples over the record. 40 coupled pairs: published vs aTDS vs aTDS capped at 1000; report score, and the median absolute error between tau0 of stable windows and the true delay at the window centre.

## E4 — Stage-aware aTDS on the 35-subject sleep cohort

Using `ncom_cohort.py` machinery (same data path, stage assignment, 100-draw cross-subject stage-matched null, effects with 5 000 bootstrap resamples):

1. Report the distribution (median, IQR) of stage-episode lengths per stage (consecutive runs of the same label).
2. Setting `atds_stage`: window cap `L = clip(round(0.5 * median episode length over all W/LS/DS/REM episodes), 61, 300)`; lag `(L-1)//2`; **tolerance per stage**: for stage s, the largest t in `0..max(1, lag//4)` such that, over 50 random cross-subject stage-matched surrogate pairs (pooled over links), mean + 2 SE of the stage-s TDS is <= 5 %.
3. Score `published`, `atds_cap300` (reuse) and `atds_stage`; report per-stage means, the effects table (mean link TDS and number of significant links; LS-DS, W-DS, REM-DS, LS-REM) and paired gains vs published, plus the number of subjects per contrast.

## E5 — More real records in one field (hydrology)

Using `experiments/atds_real/fields.py` loaders/manifest: add airGR catchments `L0123001` (daily) and `L0123002` (daily; from the same airGR package data, `https://raw.githubusercontent.com/cran/airGR/master/data/L0123002.rda`), next to the hourly `L0123003` already run. Expected delay: report the domain expectation from the manifest where given (L0123001: 0-3 days) and, separately labelled, the full-record cross-correlation peak. Settings published and aTDS; 200 rotations; same outputs as Part F.

## Acceptance

- `python -m pytest tests/ -q` passes; `git diff --stat pyTDS/` empty; no raw private data committed.
- Each experiment writes JSON + CSV; the README lists wall time and the headline numbers.
- Report back the headline tables for E1–E5 and any deviation from this spec with the reason.
