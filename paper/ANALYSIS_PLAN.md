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

## Facts from the approved thesis PDF (`WirtingThesis/Thesis-Asaf Markuza_מאושר להדפסה.pdf`; English text identical to `Thesis-Asaf Markuza.pdf`)
- D1 "was taken from the same patient[s] used in [7]" = Bashan et al. 2012 Nat Commun (thesis p.22, Fig 11 caption).
- Thesis Fig 12 text: threshold 9 % TDS on D1 (2000 surrogates); stored `thres` = 8.6 (105/140 links significant).
- Apnea cohort: 7.8 % threshold (sleep-stage surrogates, Fig 17), 14 + 14 subjects (p.13).

## Pilot already done
`scripts/thesis_port.py` (faithful port of the MATLAB TDS: linear `xcorr(..,30,'coeff')`, max |xc|,
60/30 windows, `NL=floor(2N/L-1)`, ±1 / 4-of-5, window time = sample (k-1)*30+30) reproduces
`real_links` of **all 35 subjects exactly** (4900 values, max |diff| = 0.0; `scripts/r1_replicate_real.py`, `results/r1_real_links_check.json`).

## R1 Thesis replication on D1 (thesis setting)
1. `real_links` for all 35 subjects with the port; report max |diff| vs D1-ref.
2. 2000 surrogate subjects with the full thesis method — three parts, all to be ported:
   (a) `GetSurPatient.m` pass 1: 10 distinct random subjects (`randperm(35,10)`), signal k from subject k,
   all truncated to the common length; per stage in order LS, DS, REM, awake, each signal's in-stage samples
   are concatenated and truncated to the shortest; stage labels are taken from the 10th signal's subject;
   (b) pass 2 ("add more Data"): same again with a new draw, of which only the **first half**
   (`arr_data(1:end*0.5,:)`) is appended; (c) `randomizeSurData.m`: 5 rounds, each cutting random
   blocks (`randi` limits) out and re-concatenating them (labels move with the data) until < 5 % remains,
   which is appended. Seeds differ from MATLAB, so compare **distributions**, not values: per link x stage
   mean and 95th pct of ours vs stored `sur_links`; KS distance.
3. Thesis threshold (`looking_for_threshold.m`): MATLAB `ttest2` default = **equal-variance** two-sample t,
   real_links (35) vs sur_links per link x stage, significant if p < 1e-3; "mean real" is the **pooled**
   mean over windows of all subjects (`Stages_cell.mat` pooling, as `makeDataForReconFig2.m` builds it),
   NOT `real_links.mean(0)` (they differ by up to 4.18 points); threshold = first t in 1:0.1:20 where
   100 % of links with pooled mean > t are significant. Report ours vs stored `thres` (8.6) and the
   thesis text (9 %).
4. pyTDS (circular `pbc_xcorr`) at the thesis setting (`TDSParams(window=60, overlap=30, max_lag=30,
   tolerance=1)`, both `stability_anchor` values) vs the port: per-link-stage mean diff and correlation —
   explains pyTDS vs MATLAB differences (also the 43.6 % benchmark question).
5. Figures re-drawn: thesis Fig 11 (rank plot real vs 2000 surrogates per stage), Fig 12 (threshold
   curve), Bashan-Fig-2-style 10x10 stage matrices (`results/mats.jpg` in thesis repo).

## R2 thesis-TDS vs aTDS on D1 (revised after review, before any R2 run)
Variants (reuse `experiments/adaptive_groundtruth/common.py:variant_params`, fs=1, cap 300 s as in
the Sleep-EDF study): `published` (= thesis setting, 60/30/±30/±1, pyTDS core), `win_only`,
`tol_only`, `v2_calibrated` (aTDS), `v2_cal_lag3`. Params are derived per subject from that subject's
10 signals. Cap sensitivity (subject-level only): `v2_calibrated` with caps 120 s and 600 s.

Known before running (review on 3 subjects): the Bartlett rule asks for 5 500–14 000 s windows (EEG
delta), so the aTDS window is **set by the 300-s cap**, tolerance ≈ ±14–15, max lag ±149;
calibration did not fail. The paper must say that here "aTDS" = 300-s window + calibrated tolerance.

Stage labelling (same for every variant): a window is "pure k" if ≥ 80 % of its samples are stage k
(`sleep_edf.py:window_stage`). **Primary rule (strict):** window i counts for stage k only if windows
i−2..i+2 are all pure k, so the whole centred stability group lies in stage k (a stability label is
computed from 5 consecutive windows; at 300/150 s these span 900 s and otherwise mix stages). Sensitivity:
80 %-only rule. Stage score NaN if < 5 counted windows.

Outcome per subject x stage: **connectivity = mean over the 35 thesis links of (real − chance)**, chance =
same pair with one signal circularly shifted (5 shifts, shift in [max(2W, 5·lag), T−that]); reviewer
measured the subject-level chance SD at 0.1–0.6 points vs effects of 5–45. Raw TDS % is NOT compared
across methods. Secondary: all 45 pairs; raw scores.

Hypotheses (direction from Bashan et al. 2012, PMC3518900: connectivity low in DS and REM, high in LS and
wake). For `published` they are a **replication**, not blind: the stored thesis outputs already show
them (dz LS>DS 1.85, awake>DS 1.11, LS>REM 1.30, awake>REM 0.86, computed by the reviewer from real_links).
- H1 LS > DS, H2 awake > DS, H3 LS > REM, H4 awake > REM (paired). Wilcoxon signed-rank, Holm over the
  4 within each variant.
- H5 (aTDS claim): dz(v2_calibrated) > dz(published). dz of both and the gain are computed on the
  **same subjects** (defined for both variants in that contrast); gain CI = 10 000 paired subject-bootstrap.
  **Decision rule:** H5 supported if the gain CI excludes 0 (positive) in ≥ 3 of 4 contrasts; refuted if
  any gain CI is entirely < 0 in ≥ 2 contrasts; otherwise inconclusive. H2 and H4 are expected to be
  **underpowered for aTDS** (awake undefined in ~14/35 subjects at 300 s).
- Statistics code: new (bootstrap + Wilcoxon), modelled on `bootstrap_effects.py`.

Network-level (each method against its own null; replaces the thesis surrogate for R2 because the
thesis surrogates cut stages into ~354-s pieces, leaving 300-s windows almost no pure awake/DS windows):
per subject and link, **50 circular-shift null draws** scored with the same stage windows; per link x stage
AUC = P(real subject score > null score) over 35 real vs 35x50 null; link significant if one-sided
Mann-Whitney p < 1e-3 (thesis alpha). H6: mean AUC over links x stages, aTDS ≥ published; also count
of significant links per stage per variant.

Mechanism: chosen window/tolerance/max-lag per subject, `calibration_failed` count, raw (uncapped)
Bartlett window, independent samples per 60-s window per signal (reuse `independent_samples.py` logic),
n subjects with defined score per stage per variant.

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
Public apnea cohort (e.g. PhysioNet `ucddb`, reachable on the AWS mirror) — not run (author's decision 2026-10-05: NCOM data only); future work.
aTDS on the Charité apnea cohort — impossible here (inputs unavailable).
