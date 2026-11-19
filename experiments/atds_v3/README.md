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
