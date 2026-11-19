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

## X6 follow-ups (`SPEC_X6.md`, `x6_followups.py`)

```
python experiments/atds_v3/x6_followups.py 1 <fields_data_dir>   # items 1 and 2 (and the real records); also: 3 <fields_data_dir>, 4, 5, all  (--smoke for a tiny run)
```
`rescaled_est` = published setting x r_hat, r_hat = max(1, round(tc_hat / 2)), tc_hat = first lag where the ACF of the slower signal falls below 1/e (median over the 5 training pairs; the record itself for real data). No window cap is applied to it. Wall times (s, 4 workers): item 1 + 2 + real records 133, item 3 8, item 4 41, item 5 1122 (`results/logs/times.txt`).

| Item | Files |
|---|---|
| 1 synthetic (E3 regime, E3 drift, exp4 U9, exp4 U6) | `x6_item1_synthetic.json/.csv` (includes a check of the re-implemented U9 / U6 against the stored exp4 results: identical AUC) |
| 1 real records (200 rotations, published / aTDS / rescaled_est) | `x6_item1_real.json/.csv` |
| 2 boundary margin | `x6_item2_margin.json/.csv` |
| 3 L0123002 diagnosis, capped reruns, stable-window counts | `x6_item3_L0123002.json`, `x6_item3_stable_windows.csv` |
| 4 n_eff estimators | `x6_item4_neff_estimators.json/.csv` |
| 5 per-link power | `x6_item5_power.json/.csv` |

Headline numbers.
- **Item 1**. Regime switching (A AUC; B false stability coupled/control %): published .944, 0.89/0.81; aTDS 1.0, 8.12/5.78; aTDS per regime 1.0, 9.00/4.50; rescaled_est (r_hat 8, window 480, tol 8) 1.0, 1.96/0.27. Drift (score %; median |err|): published 8.5, 0.30; aTDS 86.5, 0.25; rescaled_est (r_hat 4) 84.6, 0.24. U9 slow pair AUC/hit: published .438/0, per-pair aTDS .983/.97, rescaled_est (r_hat 18) .988/.97; fast pair 1.0/1.0 for per-pair aTDS and rescaled_est (r_hat 1); fast-slow pair: published .828/.57, aTDS 1.0/1.0, rescaled_est (r_hat 16) 1.0/1.0. U6 d=8 / d=30 (AUC, hit): published .974/.43 and .462/0; aTDS 1.0/1.0 and .995/0 (reported 11); rescaled_est (r_hat 4) 1.0/1.0 and 1.0/1.0 (reported 30). Real records (score %, p, delay samples): see `x6_item1_real.csv`; rescaled_est: L0123003 65.2 / .005 / 14 h; L0123001 93.0 / .005 / 2 d; L0123002 0 / 1.0 / none (0 stable windows, r_hat 19, window 1140); energy 100 / .005 / 2; gasoline = published (r_hat 1).
- **Item 2**: with the one-window margin the aTDS windows (516) leave 8 / 10 windows per pair in A / B; aTDS per regime has no regime-B window (its B window is 1000 > block / 3), so B false stability is not reportable there; regime-A AUC stays 1.0 for all aTDS-type settings and .942 for published. B false stability with margin (coupled/control %): published 0.89/0.77, aTDS 3.75/5.50, rescaled_est 0.38/0.19.
- **Item 3**: L0123002 flow has Bartlett factor 38.8 (precipitation 1.8), window before the cap 1165, capped by the 10 % rule at 1059, 19 windows, 4 stable, tolerance 9, rotation null mean 2.26 %, 95th percentile 21.05 % (= the real score, p = .104). The flow ACF is strongly seasonal (0.62 at lag 365): the first-zero-crossing Bartlett sum ignores it, windows end up about 3 years long. Capped at 365 d: 0 stable windows (p 1.0); at 120 d: score 10.3 %, null95 7.4 %, p .015, delay 1 d (the domain-plausible one).
- **Item 4** (AUC; delay found): r=1 .995/1.0/1.0 (a/b/c), r=4 .923/.938/.991, r=8 .877/.867/.961 (windows 563/736/1187 at r=8). Chance 3.3-4.9 % everywhere.
- **Item 5** (power / FPR, 50 + 50 pairs, Wilson 95 % CIs in the CSV): r=2 published .52/.02, rescaled .98/.02, rescaled_est .96/.00, aTDS .94/.00; r=4 .08/.04, .82/.00, .82/.00, .66/.06; r=8 .02/.04, .48/.08, .48/.08, .48/.04.

Deviations: item 4 reads "r = 1, 4, 8, 100 + 100 pairs" as r in {1, 4, 8} with 100 coupled + 100 uncoupled pairs; the margin variant requires >=1 window per regime (the 80 % purity rule >=3); rescaled_est in U9 is per pair (plus a system-wide variant `rescaled_est_system`).
