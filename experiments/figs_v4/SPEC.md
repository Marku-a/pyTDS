# SPEC — extra figures for the simple aTDS paper (Part 2 simulations)

Goal: four new PNG figures that show, as graphs, results that the paper `paper/simple/atds_simple.tex` currently reports only in
tables/text. Compare **published (60/30/±30/±1)** vs **aTDS** only (do NOT plot `rescaled`, `rescaled_est`, oracle, GCC, BH or `atds_bh`).

Rules: TDS only via this repo's pyTDS (directly or through `experiments/atds_real/atds21.py`, `experiments/atds_v3/common.py`,
`experiments/atds_v3/e3_synthetic.py`, `atds_brief/experiments/exp4.py`). Do not modify `pyTDS/` or any existing experiment script or
result file. New code only in `experiments/figs_v4/`. Max 2 worker processes; keep memory < 6 GB (the container has 16 GB shared).
Style: copy rcParams, colours (published `#2a78d6`, aTDS `#eb6834`, ink `#262624`, muted `#8a8a85`) and the serif look from
`experiments/atds_real/plot_real_summary.py`; width 6.8 in; dpi 200; no titles overlapping (check by viewing the PNG).
Output PNGs to `paper/simple/figures/`.

## 1. `plot_classes.py` → `fig_signal_classes.png` (from result files only)
Source: `atds_brief/results/exp4.json` (U5 fGn H 0.6/0.75/0.9; U6 rhythm d 8/30; U7 or U8 fast signal d=50 with/without d_max;
U10 changing dynamics; U9 mixed / fast-slow) and `atds_brief/results/exp5.json` V1 (Rulkov d 16/32). Find the right keys by reading
`exp4.py`/`exp5.py`. Rows = the rows of table `tab:signals` in the paper (same order). Panel (a): AUC, horizontal dot plot, published
and aTDS per row, with 95 % CI whiskers where the JSON has them; vertical line at 0.5. Panel (b): "delay found" (%) per row, same layout.
Where a row has two aTDS variants (fast d=50: default and d_max=60; U10: one window), plot the default aTDS as the main dot and the
variant as a hollow marker, labelled in the legend. Every plotted number must equal the number in `tab:signals`; print a check table
(row, published AUC, aTDS AUC, found %) and list any mismatch.

## 2. `plot_networks.py` → `fig_networks.png` (from result files only)
(a) Per-link detection power and false-positive rate vs slow-down r ∈ {2,4,8} for published and aTDS with Wilson 95 % CIs, from
`experiments/atds_v3/results/x6_item5_power.json` (keys `r2|published`, `r2|atds`, ...). Power solid lines, FPR dashed, a dotted line at 5 %.
(b) Network test from `atds_brief/results/exp6.json` W1: per design (disjoint, common driver): true links found out of 30 and false links
in coupled and uncoupled networks, for aTDS (and published if W1 has it). Bars with counts written on them.
Print the numbers plotted.

## 3. `plot_hard.py` → `fig_hard_cases.png` (from result files only)
Three panels. (a) Regime switching, from `experiments/atds_v3/results/e3_synthetic.json` and `x6_item2_margin.json`: for published, aTDS,
aTDS with boundary margin: % stable windows in coupled regime A (coupled pairs), false stability in uncoupled regime B (coupled pairs) and
in controls; horizontal dotted line at 5 % budget; values written on bars. (b) Drifting delay (e3 `drifting_delay`): % stable windows and
median absolute delay error (samples) for published vs aTDS (two small grouped bars or bar + annotation). (c) Nonlinear coupling
(e3 `nonlinear_monotone`, `nonlinear_non_monotone`): AUC published vs aTDS for monotone and non-monotone. Print the numbers plotted; they
must match the paper text in `\subsection{Coupling that switches, drifts or is nonlinear}` (published 5.6 / aTDS 98.9 % regime A;
false stability 0.9 / 8.1 %; margin 3.8 %, controls 5.5 %; drifting 8.5 / 86.5 %, errors 0.30 / 0.25; AUC 0.997 / 1.0 and 0.803 / 0.747).

## 4. `plot_sim_traces.py` → `fig_sim_traces.png` (new computation, one seed each)
Same look as `paper/simple/figures/fig_real_traces.png` (script `experiments/atds_real/plot_real_traces.py`): rows = cases, left column
published, right column aTDS; each panel scatters τ0 per window (grey unstable, coloured stable), title "published: TDS xx %" /
"aTDS: TDS xx %", dashed line = true delay in pyTDS sign convention (τ0 = −d when x leads), shaded band = lag range, time on x-axis.
Rows:
  (i) slow-down benchmark r = 4 (U1 design: tc = 2r, d = 5r, c = 0.3, N = 20000), first coupled test pair; aTDS parameters exactly as U1
      (set on its 5 training pairs) — reuse exp4/E1 code.
  (ii) regime switching (E3 design), first coupled pair; shade the coupled regime-A blocks lightly; true delay line only over A blocks;
      aTDS parameters as in E3 (whole-record aTDS).
  (iii) drifting delay (E3 design), first pair; true delay drawn as the drifting line.
Use the same generators/seeds as the original experiments (import them; do not re-implement). Print each panel's score and the params.
Sanity: the printed aTDS/published params must equal those reported in the result JSONs for the same experiment (e.g. regime-switching
aTDS window 516, tolerance 29).

## Acceptance
- `python -m pytest tests/ -q` passes; `git diff --stat pyTDS/` empty; no existing file modified except new PNGs in `paper/simple/figures/`.
- Each script runs with `python experiments/figs_v4/<script>.py` from the repo root.
- Report: the check tables printed by each script, wall times, any mismatch with the paper numbers, and view each PNG yourself
  (Read the image) to confirm no overlapping text.
Out of scope: editing the paper .tex, committing.
