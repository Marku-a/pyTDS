# SPEC — rename "published" to "fixed TDS" in all paper figures; new real-data figures

Rules: TDS only via this repo's pyTDS / `experiments/atds_real/atds21.py`. Do not modify `pyTDS/`. You MAY edit the plotting scripts
listed below (label strings and the new figures only — no change to any computation or result file). Max 2 processes, memory < 6 GB
(a previous run hit 14 GB by loading `data/data-*.txt` without `[:, 1]`). Output PNGs to `paper/simple/figures/` (overwrite).
Style: as `experiments/figs_v4/_style.py`; fixed TDS `#2a78d6`, aTDS `#eb6834`.

## A. Relabel
In every figure used by `paper/simple/atds_simple.tex` (`grep includegraphics`), every visible occurrence of "published" (legend,
title, axis, annotation; any case) becomes "fixed TDS" (e.g. "published 60/±30/±1" → "fixed TDS 60/±30/±1", "published: TDS 2 %" →
"fixed TDS: TDS 2 %"). Find each figure's generating script (search the repo: `paper/make_figures.py`, `experiments/atds_v3/plot_v3.py`,
`experiments/atds_real/plot_real.py`, `plot_real_traces.py`, `plot_real_summary.py`, `experiments/figs_v4/*.py`, `atds_brief/...`;
copies were made from `paper/paperorchestra/final_run3/figures/`). Regenerate each and confirm by viewing that nothing else changed.
Skip `fig_networks.png`, `fig_hard_cases.png`, `fig_sim_traces.png`, `fig_fields_significance.png` (being removed from the paper).
If a figure's script needs data not available (e.g. private cohort raw data), check whether it only needs committed result files; if it
truly cannot be regenerated, report it — do not edit PNGs by hand.
Some scripts need `ATDS_CODE=/tmp/claude-0/-home-user-pyTDS/30705354-b4e5-5b06-8b88-52507cb915b6/scratchpad/atds-code`;
field data dir: `/tmp/claude-0/-home-user-pyTDS/30705354-b4e5-5b06-8b88-52507cb915b6/scratchpad/fields_data`;
BIDMC CSVs: `/tmp/claude-0/-home-user-pyTDS/30705354-b4e5-5b06-8b88-52507cb915b6/scratchpad/bidmc_raw`.

## B. Split `fig_real_traces.png` (edit `experiments/atds_real/plot_real_traces.py`)
1. `fig_physio_traces.png`: rows EEG benchmark (as now) and ECG→finger pulse patient 6 (as now), columns fixed TDS | aTDS.
2. `fig_fields_traces.png`: **aTDS only**, one panel per record, 5 rows × 1 column (or 3×2 grid with the 6th cell used for (C)), records:
   rainfall→river hourly (L0123003; expected 8 h), rainfall→river daily L0123001 (`hydrology_airgr_L0123001_daily.csv`; expected 0–3 days,
   loader as in `experiments/atds_v3/e5_hydrology.py`), daily L0123002 (expected 0–3 days), temperature→electricity (expected 1 h),
   crude oil→gasoline weekly changes (expected 1–2 weeks). Use the aTDS parameters exactly as in the result files
   (`experiments/atds_real/results/fields_v21.json`, `experiments/atds_v3/results/e5_hydrology.json`; print and assert they match).
   Each panel: tau0 per window in physical units (grey unstable, orange stable), the expected delay as a shaded band (range) or dashed line
   (point), y zoomed around the expected delay (show the whole aTDS lag range only if needed to show where the estimates are), title
   "<record>: aTDS delay <found> (expected <x>), TDS <score> %, p = <p>" using the result-file values.
## C. `fig_fields_summary.png`: estimated vs expected delay (aTDS only)
One panel: for each of the 5 records, the expected delay (range as a vertical bar, in hours, log y axis or per-record normalised — choose
what is readable) and the aTDS delay found as an orange dot; mark non-significant results (p > 0.05) hollow. Values from the result files.
## D. `fig_bidmc_delay.png` (from `experiments/atds_real/results/bidmc_ecg_ppg.json` only)
Scatter per patient: x = reference delay (`xcorr_delay`/125 s, the full-record alignment within half a beat), y = delay found
(fixed TDS blue, aTDS orange; only patients with a delay), both wrapped to [-rr/2, rr/2) relative to the reference as in `bidmc.py`'s
`wrap` (i.e. plot reference + wrap(found − reference, rr)); identity line; text with counts within 0.04 s (fixed TDS 22/41, aTDS 40/53).
Second panel: histogram of (found − reference) in seconds for both settings.

## Acceptance
- Every regenerated PNG viewed; no "published" left in any figure the paper uses: report the list figure → script → status.
- `python -m pytest tests/ -q` passes; `git diff --stat pyTDS/` empty; no result JSON/CSV changed (`git status`).
- Report printed values for B2, C, D. Out of scope: paper .tex, committing.
