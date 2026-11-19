"""Part B: benchmark recording M001602 (source of the 43.6 % TDS benchmark).

Usage: python benchmark_m001602.py <dir with *M001602.dat and *M001602_sleepstages_conv.dat> <out_dir>
(the cohort *All_data_clean.mat is auto-detected in the same dir for the consistency check)
Private data are read from the given path only; only derived numbers are written.
"""
from __future__ import annotations

import glob
import json
import os
import sys
import time
from dataclasses import replace
from itertools import combinations

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import atds21 as A  # noqa: E402
from pyTDS.core import stable_label, time_delay_interaction  # noqa: E402
from pyTDS.core import tds as pytds  # noqa: E402

CH = ["HR", "Resp", "Chin", "Leg", "Eye", "delta", "theta", "alpha", "sigma", "beta"]
HERE = os.path.dirname(os.path.abspath(__file__))


def load_cohort_subject(mat_path, pid="M001602"):
    import matio
    d = matio.load_from_mat(mat_path)
    r = d["p_data_stucrt"][0, 0]
    names = [str(np.asarray(n).ravel()[0]) for n in r[0].ravel()]
    return r[1].ravel()[names.index(pid)]


def main(data_dir, out_dir, mat_path=None):
    t0 = time.time()
    os.makedirs(out_dir, exist_ok=True)
    f_dat = [f for f in glob.glob(os.path.join(data_dir, "*M001602.dat"))]
    f_st = glob.glob(os.path.join(data_dir, "*M001602_sleepstages_conv.dat"))
    assert len(f_dat) == 1 and len(f_st) == 1, (f_dat, f_st)
    # the file holds the 10-column table followed by unrelated trailing text/numbers: keep the leading 10-field rows
    rows = []
    with open(f_dat[0]) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) != 10:
                break
            rows.append([float(v) for v in parts])
    raw = np.array(rows)
    assert raw.shape[1] == 10, raw.shape
    nan_counts = {c: int(np.isnan(raw[:, i]).sum()) for i, c in enumerate(CH)}
    sig = np.column_stack([A.interp_nan(raw[:, i]) for i in range(10)])
    st = pd.read_csv(f_st[0], sep=r"\s+", header=None, names=["t", "code", "label"])
    out = {"n_samples": int(len(sig)), "nan_counts": nan_counts,
           "stage_epochs_30s": st["label"].value_counts().to_dict()}
    d, s = sig[:, CH.index("delta")], sig[:, CH.index("sigma")]

    # 1. published, end anchor and default anchor
    pe = replace(A.PUBLISHED, stability_anchor="end")
    out["published_end_anchor"] = {"tds": float(pytds(d, s, pe)["score"])}
    out["published_default_anchor"] = {"tds": float(pytds(d, s, A.PUBLISHED)["score"])}

    # The 43.6 % MATLAB benchmark belongs to data/data-delta.txt and data/data-sigma.txt in this
    # repository, which are NOT the delta/sigma columns of M001602.dat (different length and scale,
    # r ~ 0 at any shift). Report the benchmark on those files separately.
    root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    td = np.loadtxt(os.path.join(root, "data", "data-delta.txt"))
    ts = np.loadtxt(os.path.join(root, "data", "data-sigma.txt"))
    td, ts = td[:, 1], ts[:, 1]  # column 0 = time, column 1 = signal
    out["repo_txt_benchmark"] = {
        "tds_end_anchor": float(pytds(td, ts, pe)["score"]), "matlab_reference": 43.6052,
        "n_samples": int(len(td)),
        "note": "data/data-*.txt (the 43.6 % benchmark) differ from the M001602.dat delta/sigma "
                "columns: lengths 29941 vs 29286, values ~1e7 vs ~1e4, correlation ~0 at shifts up to 800 s."}

    # 2. consistency with the cohort-file copy
    mats = glob.glob(os.path.join(data_dir, "*All_data_clean.mat"))
    mat_path = mat_path or (mats[0] if mats else None)
    if mat_path:
        c = load_cohort_subject(mat_path)
        cm = c.values.astype(float)
        cmi = [A.interp_nan(cm[:, i]) for i in range(10)]

        def r_at(i, lag):  # cohort sample k == file sample k + lag
            n = min(len(cm), len(raw) - lag)
            return float(np.corrcoef(cmi[i][:n], sig[lag:lag + n, i])[0, 1]), int(n)

        best_lag = max(range(0, 200), key=lambda g: r_at(CH.index("delta"), g)[0])
        out["cohort_vs_file"] = {
            "len_cohort": int(len(cm)), "len_file": int(len(raw)), "best_lag_file_minus_cohort": int(best_lag),
            "n_overlap_at_best_lag": r_at(0, best_lag)[1],
            "pearson_r_at_best_lag": {ch: r_at(i, best_lag)[0] for i, ch in enumerate(CH)},
            "pearson_r_at_lag0": {ch: r_at(i, 0)[0] for i, ch in enumerate(CH)}}

    # 3. aTDS on delta vs sigma
    p0 = A.window_params([sig[:, i] for i in range(10)])
    pairs = [(sig[:, i], sig[:, j]) for i, j in combinations(range(10), 2)]
    p, fail = A.calib(p0, pairs)
    tau, _, _ = time_delay_interaction(d, s, p)
    lab = stable_label(tau, p)
    score = float(100.0 * lab.mean())
    null = A.rotation_null(d, s, p, 200, seed=0)
    out["atds"] = {"params": {"window": p.window, "step": p.overlap, "max_lag": p.max_lag, "tolerance": p.tolerance,
                              "calibration_failed": bool(fail)},
                   "tds": score, "n_windows": int(len(tau)), **A.null_stats(score, null),
                   "delay_samples": A.delay_of(tau, lab)}
    out["wall_seconds"] = time.time() - t0
    with open(os.path.join(out_dir, "benchmark_m001602.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    # optional 3rd arg: path to the cohort .mat (for the consistency check)
    main(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else None)
