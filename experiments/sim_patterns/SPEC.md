# SPEC — simulated delay patterns: does aTDS recover the true delay?

Rules: TDS only via this repo's pyTDS (`pyTDS.core.time_delay_interaction`, `stable_label`, `tds_score`) and the aTDS rule in
`experiments/atds_real/atds21.py` (`window_params`, `calib`, `PUBLISHED`, `delay_of`). Do not modify `pyTDS/` or any existing file.
New files only in `experiments/sim_patterns/` and the PNG in `paper/simple/figures/`. Fixed seeds. At most 3 worker processes; memory < 6 GB.
Naming: the setting `PUBLISHED` (60/30/±30/±1) is called **"fixed TDS"** in all labels, legends, JSON keys (`fixed`) and printouts.
Style: copy rcParams/colours from `experiments/figs_v4/_style.py` (fixed TDS `#2a78d6`, aTDS `#eb6834`).

## Generator (`patterns.py`)
x: AR(1), correlation time tc = 8 (phi = exp(-1/tc)), N = 20000, z-scored. y(t) = c * x(t - d(t)) + sqrt(1-c^2) * eta(t), eta an independent
AR(1) with the same tc, z-scored; c = 0.6; d(t) rounded to integer samples, x indexed with clipping at 0 (first max(d) samples of y are pure eta).
Uncoupled control: same but c = 0. Patterns d(t), t = 0..N-1:
1. `fixed`: d = 20.
2. `step`: d = 10 for t < N/2, 30 after.
3. `drift`: linear from 10 to 30.
4. `harmonic`: d = 20 + 10 sin(2 pi t / 10000)  (two full periods).
5. `onoff`: d = 20 in blocks of 2000 samples alternating coupled (c = 0.6, first block coupled) / uncoupled (c = 0).
Seeds: pattern k (1..5), pair j: training pairs j = 0..4 use seed 1000*k + j; scored coupled pairs j = 0..39 seed 1000*k + 100 + j;
controls seed 1000*k + 500 + j.

## Settings
- fixed TDS = `atds21.PUBLISHED`.
- aTDS: `window_params` on the 5 training pairs' signals (all 10 signals), then `calib(p, training_pairs)`; one parameter set per pattern.

## Metrics per pattern and setting (40 coupled + 40 controls)
For each window w: centre t_w = w*step + L/2 (check pyTDS's window anchor and use the same convention; `window_anchor='center'`);
true tau_w = -d(t_w) (pyTDS sign: negative when x leads); for `onoff`, true tau only for windows lying >= 80 % inside a coupled block.
- `stable_pct`: mean % stable windows (coupled pairs; for onoff: inside coupled blocks).
- `track_pct`: % of stable windows with |tau0 - true tau| <= 3 samples (coupled pairs, pooled).
- `mae`: median |tau0 - true tau| over stable windows.
- `false_stable_pct`: for controls, mean % stable windows; for onoff also mean % stable windows in uncoupled blocks of coupled pairs.
- `auc`: coupled vs controls on the whole-record TDS score.
Save `results/patterns.json` and `results/patterns.csv` (+ the aTDS params per pattern).

## Figure (`plot_patterns.py` → `paper/simple/figures/fig_delay_patterns.png`)
Rows = the 5 patterns (first coupled scored pair, seed 1000*k+100), columns fixed TDS | aTDS, same look as
`paper/simple/figures/fig_real_traces.png` (script `experiments/atds_real/plot_real_traces.py`): tau0 per window, grey unstable, coloured
stable, shaded lag range, the true delay drawn as a black dashed line/curve (-d(t)); onoff: shade coupled blocks light grey and draw the true
line only there. Row titles: "Fixed delay (20)", "Delay that jumps (10 → 30)", "Drifting delay (10 → 30)", "Harmonic delay (20 ± 10)",
"Coupling on/off (delay 20)". Panel titles "fixed TDS: TDS xx %" / "aTDS: TDS xx %". y-limits about [-60, 15] (adjust so lines are
visible; the aTDS lag band will be clipped — fine). x axis "time (10^3 samples)". Width 6.8 in, dpi 200; check no overlapping text by viewing it.

## Acceptance
- `python experiments/sim_patterns/patterns.py` then `python experiments/sim_patterns/plot_patterns.py` from repo root.
- `python -m pytest tests/ -q` passes; `git diff --stat pyTDS/` empty.
- Report the metrics table (pattern × setting: stable_pct, track_pct, mae, false_stable_pct, auc, aTDS L/lag/tol), wall time, and view the PNG.
Out of scope: paper .tex, committing.
