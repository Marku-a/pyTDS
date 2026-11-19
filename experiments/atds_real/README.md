# Real-data tests of aTDS v2.1 (pyTDS only)

Spec: `SPEC.md`. Uses only this repository's `pyTDS` (not modified). Shared rule: `atds21.py`.

- `ncom_cohort.py <All_data_clean.mat> <out_dir>` (Part A): 35-subject NCOM cohort, settings `published`, `atds`, `atds_cap300`,
  stage-matched cross-subject null (100 draws/pair). Optional smoke-test flags: `--n-subjects N --n-null N`.
- `benchmark_m001602.py <data dir> <out_dir>` (Part B): M001602 checks (the 43.6 % MATLAB benchmark itself is reproduced on `data/data-*.txt`, which are a different series from the M001602.dat columns); the cohort `.mat` in the same dir is used for the consistency check.
- `fields.py <data_dir> <out_dir>` (Part F): downloads sources and rebuilds CSVs in `<data_dir>` (needs `pyreadr`).
  **Scope reduced by user decision: only 3 of the 8 manifest datasets are run** (hydrology hourly L0123003, energy vic_elec summer,
  economics WTI -> retail gasoline weekly). Gas furnace, climate, daily hydrology, Beijing PM2.5 and Jena are not downloaded or run.
  `data_fields/MANIFEST.json` is committed; the data are git-ignored.

Seeds are fixed, 4 worker processes. The NCOM/M001602 recordings are private: pass their path on the command line; only derived
numbers are written to `results/`.
