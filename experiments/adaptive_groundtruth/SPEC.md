# SPEC: Ground-truth study of adaptive TDS on public multi-subject data

Repo: /home/user/pyTDS, branch `claude/pytds-adaptive-parameters-760t3h` (checked out). Work ONLY in
`experiments/adaptive_groundtruth/`. Do NOT edit `common.py` or `adaptive.py` (shared, owned by the caller), nor
anything under `pyTDS/`. Do NOT commit or push. Raw data goes to `/home/user/data_gt/<dataset>/` (outside the repo,
never committed). Results go to `experiments/adaptive_groundtruth/results/` (small CSV/JSON only, < 2 MB per file).

Read `common.py` fully first. It defines the 6 variants (`variant_params`), the circular-shift test (`pair_test`),
per-window labels (`window_labels`), and the statistics (`auc`, `paired_effect`, `group_effect`, `stable_seed`).
Use them as-is; never re-implement a variant or a null.

Data host: ONLY `https://physionet-open.s3.amazonaws.com/<dataset>/<version>/<file>` works (physionet.org itself is
blocked). Download with `subprocess.run(["curl", "-sSf", "--retry", "3", "-o", dst, url])` (curl is proxy-configured),
skip files already present. List a folder with `https://physionet-open.s3.amazonaws.com/?prefix=<dataset>/<version>/`
(S3 XML, `<Key>` tags; paginate with `continuation-token` / `list-type=2` if `<IsTruncated>true`). `wfdb` and `mne`
are installed; read local files only (wfdb's own downloader hits physionet.org, which is blocked).

## Common rules for every dataset script

- CLI: `--data-dir`, `--out` (default `results/`), `--n-subjects` (default all), `--n-null` (default 100),
  `--jobs` (default 2), `--smoke` (= 3 subjects, n_null 20, out `results/smoke/`).
- Parallelise over subjects with `ProcessPoolExecutor`. Seeds via `stable_seed(dataset, subject, pair, ...)`, the same
  for every variant.
- Per subject: build the system (list of equal-length finite 1-D float arrays at `fs`), call
  `variant_params(system, fs, MAX_WINDOW_S, seed=stable_seed(dataset, subject))` ONCE, then evaluate every variant on
  the same signals. Same subject set for every variant.
- Write `<ds>_params.csv` (subject, variant, `params_row(...)`, calibration_failed, window_clipped),
  `<ds>_rows.csv` (one row per subject x variant x test), `<ds>_summary.csv` (one row per variant, metrics below),
  and print the summary table. Record wall time and data provenance (file list + byte counts) in `<ds>_meta.json`.
- `MAX_WINDOW_S` is a module constant, given below. Do not tune it, the preprocessing or the subject selection after
  looking at results. If something in this spec is impossible, choose the closest option, write why in the script
  docstring, and report it.
- Sign: `pair_test(x, y)` returns `lag_samples` with pyTDS sign (negative = x leads y). Report `lag_s` as
  "x leads y by" = `-lag_samples / fs`.

## D1 Fantasia — `fantasia.py` (cardio-respiratory coupling, 20 young vs 20 elderly)

Data: `fantasia/1.0.0/` — 40 records `f1y01..f1y10, f1o01..f1o10, f2y01..f2y10, f2o01..f2o10` (`.hea .dat .ecg`).
`y` = young (21–34 y), `o` = elderly (68–85 y). Header comment has age/sex. Channels RESP and ECG at 250 Hz (some
records have a BP channel too; ignore it). Beat annotations: `wfdb.rdann(path, "ecg")`.

Preprocessing (fs = 4 Hz):
- HR: beat annotation samples with symbols in wfdb's beat set; RR = diff of beat times (s). Drop an RR if outside
  0.3–2.0 s or if it differs > 20 % from the previous accepted RR. Instantaneous HR = 60/RR (bpm) placed at the
  second beat's time; linear interpolation onto the 4 Hz grid (no extrapolation: crop to the first/last beat).
- RESP: 4th-order Butterworth low-pass 1 Hz (`filtfilt`), then take every 62.5th sample via
  `scipy.signal.resample_poly(x, 2, 125)` (250 -> 4 Hz).
- Crop both to their common span, then to the first 110 minutes (26,400 samples). Subjects shorter than 60 min: skip
  and report.
- System = [resp, hr]. `MAX_WINDOW_S = 300`.

Tests per variant:
1. within: `pair_test(resp_i, hr_i, params_i)` (breathing -> HR).
2. cross-subject null (no true coupling): `pair_test(resp_i, hr_j, params_i)` with j = next subject in sorted record
   order (wrap around), truncated to the common length.

Summary per variant: `detect` (within p < 0.05), `xsubj_fpr` (cross p < 0.05), `auc_score` = auc(within scores,
cross scores), `auc_excess` = same on (score − chance95), `young_vs_old` = `group_effect` on within (score −
chance95) young vs elderly (d, auc), plus the same on raw score; `lag_s` median/IQR over detected subjects; params
medians (window_s, max_lag_s, tolerance_s); calibration_failed count.
Pre-registered expectation (literature: cardio-respiratory coupling / RSA weakens with age): young > elderly.

## D2 Sleep-EDF Expanded — `sleep_edf.py` (brain-rhythm network vs sleep stage)

Data: `sleep-edfx/1.0.0/sleep-cassette/` — `SC4ssNE0-PSG.edf` + matching `SC4ssN??-Hypnogram.edf` (suffix letters
vary; match on the first 7 chars `SC4ssN`). Use night 1 (N = 1) of the first 20 subjects in sorted order (ss =
00, 01, ...), skipping missing ones. Read with `mne.io.read_raw_edf(..., include=["EEG Fpz-Cz", "EEG Pz-Oz"])` and
`mne.read_annotations(hypnogram)`. Delete each PSG file after its band powers are cached as `.npz` in the data dir
(disk budget).

Preprocessing (fs = 1 Hz):
- For each EEG channel (100 Hz): 2-s Hann windows stepped 1 s; power spectrum; band power in delta 0.5–4, theta 4–8,
  alpha 8–12, sigma 12–16, beta 16–30 Hz; take log10. -> 10 signals at 1 Hz, sample k = window starting at k s.
- Hypnogram -> per-second stage: W; N1+N2 ("Sleep stage 1/2") = LS; N3+N4 ("Sleep stage 3/4") = DS; "Sleep stage R"
  = REM; anything else (movement, "?") = NaN.
- Crop to [first non-W epoch − 30 min, last non-W epoch + 30 min] (standard for these ~20 h recordings).
- System = the 10 band signals. `MAX_WINDOW_S = 300`.

Tests per variant, per night, for all 45 pairs (x = earlier band in the list order, y = later):
- `window_labels(x, y, params)`; a window belongs to a stage if >= 80 % of its samples
  `[c − W//2, c − W//2 + W)` carry that stage, else it is dropped.
- Real per-stage score = 100 × mean stable label of that stage's windows; NaN if < 5 windows.
- Chance per-stage score: same, with y circularly shifted by 5 random shifts in [max(2W, 5·max_lag), T − that]
  (seeded), averaged.
- Night-level per stage: mean over the 45 pairs of real, of chance, and of excess = real − chance; also windows per
  stage.

Summary per variant: `paired_effect(LS, DS)`, `paired_effect(W, DS)`, `paired_effect(REM, DS)` on night-level real
AND on excess (report n, mean_diff, dz, frac_pos, sign_p); mean windows per stage; params medians.
Pre-registered expectation (Bashan et al. 2012 Nat Commun; Lin et al. 2020 Commun Biol): links weakest in DS;
LS > DS and W > DS.

## D3 Treadmill CPET — `treadmill.py` (imposed speed -> heart rate, many people)

Data: `treadmill-exercise-cardioresp/1.0.1/test_measure.csv` + `subject-info.csv` (columns: Time s, Speed km/h, HR,
VO2, VCO2, RR, VE, ID_test, ID). Breath-by-breath, irregular ~2–3 s.

Selection (deterministic): tests with >= 600 s duration, HR present in >= 95 % of rows, and >= 3 distinct speed
values. Sort by ID_test and take every k-th to get 200 tests (or all if fewer). One test per subject ID (the first
in sorted order).

Preprocessing (fs = 0.5 Hz): drop rows with NaN HR/Speed/VO2; linear interpolation of Speed, HR, VO2, VE onto a 2-s
grid from first to last time; high-pass each signal by subtracting a 120-s running mean
(`scipy.ndimage.uniform_filter1d(x, 60, mode="nearest")`), same as GarminTDS. System = [speed, hr, vo2, ve].
`MAX_WINDOW_S = 240`.

Tests per variant: `pair_test(speed, hr)` (primary) and `pair_test(speed, vo2)` (secondary).

Summary per variant, per pair: `detect` (p < 0.05 AND speed leads, i.e. lag_s > 0), `wrong_dir` (p < 0.05 AND
lag_s <= 0) as a share of significant tests, `lag_s` median/IQR over detected, params medians, n_windows median.
Pre-registered expectation: speed is imposed by the treadmill, so any significant link must have speed leading
(positive lag); HR response delays of roughly 10–60 s.

## Acceptance (each script)
- `python <script> --smoke` finishes and prints the summary; report its wall time and the table.
- Report: files created, any spec deviation with reason, data volume downloaded, and an estimated wall time for the
  full run at `--jobs 4`. Do NOT run the full run; the caller does.

## Out of scope
- Editing common.py / adaptive.py / pyTDS, other datasets, plots, commits, pushes.
