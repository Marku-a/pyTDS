# aTDS run-3 experiments

Spec: `SPEC.md`. All TDS computations use this repository's pyTDS; baselines (cross-correlation, GCC-PHAT) are numpy. Fixed seeds, 4 workers. Only aggregates are stored (the private sleep cohort is read from a path on the command line).

## Setup

```
git worktree add <scratch>/atds-code origin/claude/practical-mccarthy-n2mec6    # outside the repo
export ATDS_CODE=<scratch>/atds-code
python experiments/atds_v3/e1_baselines.py          # also e2_sensitivity.py, e3_synthetic.py (add --smoke for a tiny run)
python experiments/atds_v3/e4_cohort_stage.py <All_data_clean.mat> experiments/atds_v3/results
python experiments/atds_v3/e5_hydrology.py <fields_data_dir> experiments/atds_v3/results   # needs raw/L0123001-3.rda
```

`common.py` imports `atds_brief/experiments/exp4.py` (U1 generators, `atds2`, `boot_diff90`, ...), adds the parametrised aTDS (`atds_var`) and the numpy baselines (`xcorr_full`, `gccphat_full`; delay reported in the TDS convention, tau = -d).

## Results (`results/`)

| Experiment | Files | Wall time |
|---|---|---|
| E1 baselines | `e1_baselines.json/.csv` (AUC + CI, delay found + Wilson CI, aTDS-minus-baseline AUC CIs, parameters, check against stored exp4 U1: identical) | 55 s |
| E2 sensitivity | `e2_sensitivity.json/.csv` (one row per constant value x case) | 146 s |
| E3 synthetic | `e3_synthetic.json/.csv` (regime switching, nonlinear monotone / non-monotone, drifting delay) | 51 s |
| E4 cohort stage | `e4_cohort_stage.json` (all), `e4_episodes.csv`, `e4_params.csv`, `e4_stage_means.csv`, `e4_effects.csv`, `e4_subject_stage.csv` | 427 s |
| E5 hydrology | `e5_hydrology.json/.csv` | 14 s |

`results/logs/` has the console logs and `times.txt` (seconds).

## Headline numbers

E1 (AUC; delay found in brackets), published / rescaled(true r) / aTDS / rescaled_est / xcorr_full / GCC-PHAT:
r=1 .997 / .997 / .995 / .997 / 1.0 / 1.0; r=2 .888 / .993 / .979 / .993 / 1.0 / 1.0; r=4 .590 / .955 / .923 / .955 / 1.0 / 1.0;
r=8 .541 / .806 / .877 / .806 / 1.0 / 1.0; r=16 .509 / .688 / .645 / .691 / 1.0 / 1.0.
The full-record global methods are perfect on this stationary benchmark, so aTDS does not beat them (aTDS-minus-xcorr AUC at r=16: -0.355, CI -0.43 to -0.29). aTDS beats published at every r >= 2 and beats true-r rescaling at r=8 only.

E4: episode lengths (s, median [IQR]): W 90 [60-210], LS 420 [120-1118], DS 480 [150-1050], REM 480 [180-900]; all 300. `atds_stage`: L=150, lag 74, tolerance W 6, LS 6, DS 6, REM 4.
Effects tables: `e4_effects.csv`.

E5 and the remaining tables: see the CSV files.

## Headline tables (all experiments)

**E2 sensitivity** (AUC at r = 1 / 4 / 8; default n_eff 30, budget 5 %, z 2, floor 61, 20 rotations: .995 / .923 / .877):
n_eff 15 .997/.717/.684, 45 1.0/.946/.899, 60 1.0/.991/.954 (chance 5.0-5.4 % at r 4, 8); budget 1 % 1.0/.920/.884, 10 % .981/.874/.851; z 0 .995/.911/.861, z 3 .995/.923/.881; floor 31 = 61, floor 121 1.0/.923/.877; rotations 10 .995/.921/.888, 40 .995/.923/.862. Rhythm threshold 0.2 and 0.3: lag capped at +-20, delay 30 not found (hit 0); 0.4: no cap, lag +-158, tol 16, both delays found.

**E3 harder benchmarks**: regime switching (regime-A AUC; regime-B false stability coupled / control): published .944, 0.89 / 0.81 %; aTDS 1.0, 8.12 / 5.78 %; aTDS per regime 1.0, 9.00 / 4.50 %. Nonlinear monotone AUC: published .997, aTDS 1.0, xcorr 1.0; non-monotone: published .803, aTDS .747, xcorr .341. Drifting delay: score published 8.5 %, aTDS 86.5 %; median |tau0 - true| 0.30 vs 0.25 samples.

**E4 stage-aware cohort** (window 150 s, tolerance W/LS/DS 6 s, REM 4 s): per-stage mean link TDS W/LS/DS/REM published 16.9/19.8/13.7/11.4, cap300 40.2/46.9/35.8/36.6, stage 31.8/36.4/27.0/22.6. Gains atds_stage minus published, mean link TDS: LS-DS +0.28 [-0.01, 0.57], W-DS -0.05, REM-DS -0.13, LS-REM +0.26 (all CIs include 0); significant links: W-DS -1.02 [-1.48, -0.68], others include 0. n subjects: 34 (LS-DS, REM-DS), 33 (W-DS), 35 (LS-REM).

**E5 catchments** (score, p, delay; xcorr peak): L0123001 daily published 59.1 %, .005, 1 d; aTDS 97.6 %, .005, 2 d; xcorr +2 d (r .301). L0123002 daily published 32.4 %, .005, 1 d; aTDS 21.1 %, .104, 162 d; largest |C| is a negative edge artefact (-30 d, r -.069), positive peak +1 d (r .054). L0123003 hourly published 1.2 %, .02, 13 h; aTDS 79.5 %, .005, 8 h; xcorr +8 h (r .499).

**Deviations from the spec**: global methods use the unbiased normalisation 1/(N-|k|) and report score = max |C| (sign-agnostic); E3 noise weighting sqrt(1-c^2) as in coupled_pair; per-regime calibration on concatenated regime blocks; E4 surrogate seeds differ from the 100-draw null seeds; E5 reports both the |C| peak (with edge flag) and the positive peak; the rhythm threshold is varied by patching adaptive.dominant_period(min_peak) temporarily.
