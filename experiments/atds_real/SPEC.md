# Real-data tests of aTDS v2.1 (pyTDS only)

Goal: add real-world evidence to the aTDS paper. Two parts:

- **Part A — NCOM cohort** (the user's MSc data: 35 overnight recordings, 10 channels at 1 Hz, per-second sleep stages) and **Part B — benchmark recording M001602** (same 10 channels; source of the 43.6 % TDS benchmark).
- **Part F — other fields**: the 8 public datasets with known delays from branch `claude/practical-mccarthy-n2mec6` (`experiments/adaptive/data_fields/MANIFEST.json`).

**Hard rules**

1. Use only this repository's `pyTDS` package (`pyTDS.core.time_delay_interaction`, `stable_label`, `tds_score`, `pyTDS.params.TDSParams`). Never TDSpy, MATLAB ports, TDSTOCK code or any other TDS implementation. Do not modify `pyTDS/`.
2. The NCOM/M001602 data are private research recordings. Read them from the path given on the command line. Never copy raw signals into the repository; commit only derived numbers (scores, effect sizes, parameters).
3. All random seeds fixed. Use `multiprocessing` with 4 workers.

## Shared module `experiments/atds_real/atds21.py`

aTDS rule v2.1 exactly as in the paper (Algorithm 1), built on `atds_brief/experiments/atds_code/adaptive.py` (import it by adding that folder to `sys.path`):

- `PUBLISHED = TDSParams(window=60, overlap=30, max_lag=30, tolerance=1)`.
- `window_params(signals, d_max=None, cap=None) -> TDSParams`: `adaptive_params(signals, window_rule="bartlett", tolerance="fixed", min_window=max(61, 2*d_max+1 if d_max else 61), max_window=cap)`; keep its rhythm cap on the lag.
- `calib(p, pairs, n_shift=20, z=2.0, seed=0) -> (TDSParams, calibration_failed)`: copy verbatim from `atds_brief/experiments/exp4.py` `calib` (rotation shift in `[max(L, n//10), n - that)`, tolerance scan `0..max(1, max_lag//4)`, keep largest t with mean + 2 SE <= 5 %). Use at most 50 pairs (random, seeded) when more are given.
- `delay_of(tau, labels)`: the most frequent tau0 among stable windows (None if fewer than 3 stable windows).

## Part A — `experiments/atds_real/ncom_cohort.py <path to All_data_clean.mat> <out_dir>`

Data: load with `matio.load_from_mat` (package `mat-io`). Struct `p_data_stucrt` with `p_name` (35 ids), `data` (35 pandas DataFrames, columns `HR, Resp, Chin, Leg, Eye, delta, theta, alpha, sigma, beta`, 1 Hz), `sleep_stages` (per-second labels `awake`, `LS`, `DS`, `REM`; same length as data). Map `awake -> W`; drop seconds with other labels from stage assignment. Interpolate NaNs linearly; report NaN counts per subject/channel.

Settings (one parameter set per setting for the whole cohort):

1. `published` — `PUBLISHED`.
2. `atds` — window: per subject compute `max_i B_i` over the 10 channels, take the **median over subjects**, `L = clip(round(30*that), 61, 0.1*median length)`, lag `(L-1)//2`, step `L//2`. Tolerance: `calib` on rotations pooled over 50 random (subject, channel pair) combinations.
3. `atds_cap300` — as `atds` but window capped at 300 s before calibration.

For every subject, every one of the 45 channel pairs, every setting: `tau = time_delay_interaction(x, y, p)[0]`, `lab = stable_label(tau, p)` (default anchor). Window k covers samples `[k*step, k*step+L)`. **Assign a window to a stage if >= 80 % of its seconds carry that stage**; otherwise exclude it. Per-stage TDS = 100 * stable / assigned windows; NaN if fewer than 3 assigned windows.

Null (cross-subject, stage-matched, as in Bashan et al. 2012): for each pair (i, j), draw 100 subject pairs (a, b), a != b (seeded); take channel i of a and channel j of b, truncate both to the shorter length; compute per-stage TDS with **stage labels of subject a**. Per-link, per-stage, per-setting threshold = 95th percentile of the null. A link is significant in a subject/stage if its TDS > threshold.

Outputs (JSON + CSV in `<out_dir>`, `experiments/atds_real/results/` by default):

- params per setting (window, step, max_lag, tolerance, calibration_failed, B summary).
- per subject x stage x setting: mean TDS over 45 links, number of significant links, number of assigned windows.
- per stage x setting: 10x10 mean TDS matrix and fraction of subjects in which each link is significant.
- effects: for each setting and each contrast `LS-DS, W-DS, REM-DS, LS-REM`, on (i) mean link TDS and (ii) number of significant links: paired Cohen dz, 95 % bootstrap CI (5 000 subject resamples, seed 0), fraction of subjects with positive difference; and the paired gain `dz(atds) - dz(published)` and `dz(atds_cap300) - dz(published)` with bootstrap CI.

## Part B — `experiments/atds_real/benchmark_m001602.py <dir with M001602.dat and M001602_sleepstages_conv.dat> <out_dir>`

`M001602.dat`: whitespace table, 10 columns in the order HR, Resp, Chin, Leg, Eye, delta, theta, alpha, sigma, beta, 1 Hz, NaN allowed (interpolate). Stages file: `time label_code label` per 30 s.

1. Published setting with `stability_anchor="end"` on delta vs sigma: report TDS (MATLAB reference 43.6 %); also with the default anchor.
2. Check the M001602 table inside the cohort file equals this file (Pearson r per channel over the overlapping part).
3. aTDS (window from the 10 channels of this subject, calibration on this subject's 45 pairs) on delta vs sigma: TDS, rotation null (200 rotations, p, null95), delay.

## Part F — `experiments/atds_real/fields.py <data_dir> <out_dir>`

1. Copy `MANIFEST.json` from `origin/claude/practical-mccarthy-n2mec6:experiments/adaptive/data_fields/MANIFEST.json` into `experiments/atds_real/data_fields/` (commit the manifest only; `.gitignore` the data).
2. Download each source URL and rebuild the 8 CSVs following each entry's `preprocessing_notes` (use `pyreadr` for `.rda`). Reuse the loaders `d_*` and helpers `match`, `phys` from `origin/claude/practical-mccarthy-n2mec6:experiments/adaptive/exp_fields.py` (copy them; keep the same detrending/differencing).
3. For each dataset and settings `published` and `atds` (window from both signals; calibration on the record itself, 20 rotations): TDS score, rotation null (200 rotations: mean, 95th pct, p = (1+k)/(1+N)), number of windows, delay = `delay_of` converted to physical units, whether it matches the expected delay (`match`), parameters in samples and physical units, calibration_failed.
4. Output `fields_v21.json` and `fields_v21_summary.csv`.

## Acceptance

- `python -m pytest tests/` still passes (46 tests); `pyTDS/` unchanged (`git diff --stat pyTDS/` empty).
- All three scripts run end to end; report wall time and the key numbers (Part A effects table, Part B 43.x % reproduction, Part F summary table).
- No raw signals committed (`git status` shows only scripts, manifest, results, README).

## Out of scope

Editing the paper; choosing which field datasets appear in the paper; any pyTDS change; the sleep-apnea cohort.
