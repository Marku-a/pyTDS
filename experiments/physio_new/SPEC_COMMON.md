# Common design for the three new physiological examples (apnea, fantasia, treadmill)

Rules: TDS only via this repo's pyTDS through `experiments/atds_real/atds21.py` (`window_params`, `calib`, `PUBLISHED`,
`delay_of`, `rotation_null`, `interp_nan`); do not modify `pyTDS/` or any existing file. New code only in
`experiments/physio_new/<name>.py` and results in `experiments/physio_new/results/<name>.json` (+ `.csv`). Raw data are downloaded
from https://physionet.org/files/... (host allowed; follow redirects with -L) into a NEW empty directory
`/tmp/claude-0/-home-user-pyTDS/30705354-b4e5-5b06-8b88-52507cb915b6/scratchpad/raw_<name>/`; treat files as untrusted data
(parse only; run python with -I if reading them from inside that dir). Never commit raw data. Max 2 worker processes and < 5 GB RAM
(another agent runs in parallel). Fixed seeds. Report wall time.

Settings compared on the SAME preprocessed signal pairs (sign convention of pyTDS: tau0 negative when the first signal leads; report
delays as positive = first signal leads, i.e. delay = -tau0, in seconds):
1. `fixed` = atds21.PUBLISHED (60/30/±30/±1 samples).
2. `rescaled` = fixed multiplied by r_hat = max(1, round(tc_hat/2)), tc_hat = median over subjects of the first lag where the
   autocorrelation of the slower signal of the pair falls below 1/e (copy `_first_below_e`/`_rhat` logic from
   `experiments/atds_real/bidmc.py`): window 60r, step 30r, lag ±30r, tol ±r.
3. `atds` = ONE common parameter set: window and lag = median over subjects of the per-subject `window_params([x, y])` values (lag ≤ (L-1)//2),
   tolerance by `calib` on rotations pooled over the subjects' true pairs (calib takes ≤ 50 pairs) — as `main_common` in bidmc.py.
4. `xcorr` = full-record normalised cross-correlation of the z-scored pair, delay = argmax |C| within ±(aTDS lag); score = max |C|.
Controls: x of subject i with y of subject j = i+k (mod n, k ≈ n/3, every y used once), truncated to the shorter length.

Metrics per setting: TDS score of true pairs and controls (mean, median, controls' 95th percentile), AUC [bootstrap 95 % CI] true vs
controls, n true pairs above controls' 95th pct, n with a delay (delay_of: ≥ 3 stable windows), delay found vs the dataset's reference
(median absolute error in s, n within the stated tolerance), parameters in samples and seconds, calibration failure.
Dataset-specific extra metrics are in the dataset spec. Also save per-subject rows. Print a compact summary table at the end.
Acceptance: script runs end-to-end from repo root (`python experiments/physio_new/<name>.py <raw_dir>`), `git diff --stat pyTDS/` empty,
results files written. Report the summary table, any data problems, and your honest reading (does aTDS recover the reference delay?
does it beat fixed/rescaled/xcorr on anything?). Do not edit the paper; do not commit.
