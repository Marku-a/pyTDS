# aTDS paper — analysis plan (written before the main runs)

Branch `claude/new-session-1at9bf` (from `claude/pytds-adaptive-parameters-760t3h`, dc2efdb).
Hypotheses H*/E* below are fixed before running; results are reported whatever they show.

## Data actually available (checked 2026-10-05)

| id | what | source | status |
|---|---|---|---|
| D1 | 35 subjects x 10 signals (HR, Resp, Chin, Leg, Eye, EEG delta/theta/alpha/sigma/beta), 1 Hz, per-second sleep stage (LS/REM/DS/awake/NA) | thesis repo `01 repreduce NCOM/surrogate test/All_data_clean.mat` (MATLAB table objects; decoded by `scripts/extract_ncom.py`) — the data of the thesis's reproduction of Bashan et al. 2012 (folder name `NCOMM-2012-data` in `makeDataForReconFig2.m`) | **inputs available** → full rerun with TDS and aTDS |
| D1-ref | thesis MATLAB outputs on D1: `real_links` (35x140), `sur_links` (2000x140), `p`, `thres`, `signi_over` | `01 repreduce NCOM/links_strength2000Sur.mat` | replication target |
| D2 | sleep-apnea cohort (Charité 2016-19): 28 patients, full-night 15x15 TDS matrices with AHI in file name; group/stage rank arrays; thresholds | `02 sleep apnea analysis/TDS_mats.mat`, `ranks_sick&healthy.mat`, `Thres.mat`, `TDS imp/stages_cells_sick&healthy.mat`, `TDS imp/patients_properties.xlsx` | **outputs only** (the 1 Hz "Embla Tables" are not available) → figures re-drawn from stored outputs; **aTDS cannot be run**, said so in the paper |
| D3 | 41 Garmin runs, 5 channels at 1 Hz | GarminTDS `run-data:share/garmin_runs.npz` | inputs available |

No raw or per-subject data are committed (pyTDS is public). Extracted arrays live outside git.
Committed results are aggregates (per-link/stage means, effect sizes), plus per-subject summary
numbers under anonymous indices only.

## Pilot already done
`scripts/thesis_port.py` (faithful port of the MATLAB TDS: linear `xcorr(..,30,'coeff')`, max |xc|,
60/30 windows, `NL=floor(2N/L-1)`, ±1 / 4-of-5, window time = sample (k-1)*30+30) reproduces
`real_links` of subjects 1-3 **exactly** (max |diff| = 0.0 over 140 values each).

## R1 Thesis replication on D1 (thesis setting)
1. `real_links` for all 35 subjects with the port; report max |diff| vs D1-ref.
2. 2000 surrogate subjects with the thesis method (`GetSurPatient.m`: 10 random distinct subjects,
   signal k from subject k, truncated to common length, then per stage the samples of each signal in
   that stage are concatenated and truncated to the shortest; stages LS, DS, REM, awake in that order,
   and the procedure is done twice ("add more Data") and concatenated — coder must port the full 96
   lines exactly). Seeds differ from MATLAB, so compare **distributions**, not values: per link x stage
   mean / 95th pct of ours vs `sur_links`; KS distance.
3. Thesis threshold: Welch? → MATLAB `ttest2` default = equal-variance two-sample t; p<1e-3 per link x
   stage; threshold = first t in 1:0.1:20 with 100 % of links with mean real > t significant. Report ours
   vs stored `thres`.
4. pyTDS (circular `pbc_xcorr`) at the thesis setting (`TDSParams(window=60, overlap=30, max_lag=30,
   tolerance=1)`, both `stability_anchor` values) vs the port: per-link-stage mean diff and correlation —
   explains pyTDS vs MATLAB differences (also the 43.6 % benchmark question).
5. Figures re-drawn: thesis Fig 11 (rank plot real vs 2000 surrogates per stage), Fig 12 (threshold
   curve), Bashan-Fig-2-style 10x10 stage matrices (`results/mats.jpg` in thesis repo).

## R2 thesis-TDS vs aTDS on D1
Variants (reuse `experiments/adaptive_groundtruth/common.py:variant_params`, fs=1, cap 300 s as in
the Sleep-EDF study): `published` (= thesis setting, 60/30/±30/±1, pyTDS core), `win_only`,
`tol_only`, `v2_calibrated` (aTDS), `v2_cal_lag3`. Plus `thesis_port` for R1 continuity.
Params are derived per subject from that subject's 10 signals.

Stage labelling for windows of any length: the Sleep-EDF rule (`sleep_edf.py:window_stage`, ≥80 % of
the window's samples in one stage, else unlabelled; stage score NaN if <5 labelled windows). For the
`thesis_port` row, also the thesis rule (centre sample) for continuity.

Outcome per subject x stage: **connectivity = mean over the 35 thesis links of (real − chance)**, where
chance = same pair with one signal circularly shifted (5 shifts, as in `sleep_edf.py`). Raw TDS % is NOT
compared across methods (different windows/tolerances → different chance levels). Secondary: all 45
pairs; raw scores.

Pre-registered hypotheses (direction from Bashan et al. 2012, Nat Commun, PMC3518900: connectivity
low in DS and REM, high in LS and wake):
- H1 LS > DS, H2 awake > DS, H3 LS > REM, H4 awake > REM (paired over subjects).
- H5 (aTDS claim) the paired effect size dz for H1–H4 is larger with `v2_calibrated` than `published`.
Statistics: paired dz, 10 000 subject-bootstrap 95 % CI (reuse `common.paired_effect`), Wilcoxon
signed-rank p; gain = dz(aTDS) − dz(published) with paired-bootstrap CI. Report n subjects with a
defined score per stage per variant (long windows can leave short stages empty — expected cost).

Network-level (thesis style, each method calibrated on its own null): 500 thesis-method surrogate
subjects per variant (params re-derived on each surrogate system, as aTDS only uses per-signal and
coupling-free information); per link x stage AUC real vs surrogate; number of links significant
(t-test p<1e-3, as thesis) and the thesis threshold per method. Expectation H6: aTDS separates real from
surrogate at least as well (mean AUC ≥ published).

Mechanism: chosen window/tolerance/max-lag per subject, `calibration_failed` count, Bartlett
independent samples per 60 s window per signal (reuse `independent_samples.py` logic).

## R3 sleep-apnea cohort (D2, stored outputs only)
Re-draw from stored outputs: Fig 20 (stage ranks healthy vs apnea, KS + Wilcoxon p), Fig 24 (link
appearance %), Fig 26 (degree PDFs), Fig 27 (sum of degrees vs AHI bins) where the stored data allow;
any figure that needs inputs not stored is listed as not reproducible. Group split AHI ≤ 10 / > 10
from file names (thesis p.13). No aTDS here.

## R4 Garmin (D3) — expectations fixed now, from the previous session's README
Rerun `GarminTDS/experiments/adaptive_compare/compare.py` at a pinned commit, unchanged.
Note: this is the same data the previous session used (export up to 2026-10-03), so R4 is a
**reproduction**, not a new test; any new run since then would be the only out-of-sample evidence.
- E1 first 30 usable runs, speed→HR detections: aTDS ≥ frozen (expected ≈ 8 vs 4); McNemar n.s.
- E2 cross-run false positives ≤ 2/30 for both.
- E3 aTDS window = 240 s cap on every run; tolerance ≈ ±14 s; lag IQR across runs wider than frozen.
- E4 held-out 4 runs: aTDS ≥ frozen.
- E5 removing the wide lag (lag = window/3) loses most of the aTDS gain.

## Out of scope / labelled untested
Public apnea cohort (e.g. PhysioNet `ucddb`, reachable on the AWS mirror) — not run; future work.
