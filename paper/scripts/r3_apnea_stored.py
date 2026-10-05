"""R3: thesis sleep-apnea results re-computed from STORED thesis outputs only (inputs unavailable).

Sources (thesis repo, 02 sleep apnea analysis):
- TDS imp/ranks_sick&healthy.mat : sick_ranks / healthy_ranks, MATLAB 105 x 4 (stages REM, DS, LS, awake), each column the
  105 upper-triangle link strengths (fraction) of the group-pooled stage map, SORTED descending
  (construct_sick_healthy_nets.m:134-186).
- TDS imp/stages_cells_sick&healthy.mat : per group, per stage, 15 x 15 cells of pooled stable labels.
- TDS_mats.mat : 28 per-patient full-night 15 x 15 TDS matrices; AHI in the file name (SL032_AHI_2.9.mat).
Reproduces results_tests.m (Fig 20-22: signrank + kstest2 between groups and between stages) and checks the
ranks against the pooled stage cells. Adds, as a labelled supplement, a link-matched test (same link in both
groups) because the thesis pairs links by rank. Writes ../results/r3_apnea_stored.json.

usage: python r3_apnea_stored.py <msc-project dir>
"""
import json
import os
import re
import sys

import h5py
import numpy as np
from scipy.stats import ks_2samp, wilcoxon, spearmanr

STAGES = ["REM", "DS", "LS", "awake"]
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "r3_apnea_stored.json")


def cell_means(f, name):
    """name -> 4 x 15 x 15 mean of pooled stable labels (stage order as stored)."""
    top = f[name][()]
    out = np.full((4, 15, 15), np.nan)
    for s in range(4):
        c = f[top.flat[s]][()]
        for i in range(15):
            for j in range(15):
                v = f[c[i, j]][()].ravel()  # HDF5 stores MATLAB arrays transposed; c[i,j] is MATLAB (j,i)
                out[s, j, i] = v.mean() if v.size and v.dtype.kind == "f" else np.nan
    return out


def main(repo):
    d = f"{repo}/02 sleep apnea analysis"
    with h5py.File(f"{d}/TDS imp/ranks_sick&healthy.mat") as f:  # the copy in TDS imp/ (where results_tests.m runs); the one in the parent dir is an older version
        sick, healthy = f["sick_ranks"][()].T, f["healthy_ranks"][()].T  # 105 x 4
    with h5py.File(f"{d}/TDS imp/stages_cells_sick&healthy.mat") as f:
        maps = {"sick": cell_means(f, "stages_cell_sick"), "healthy": cell_means(f, "stages_cell_healthy")}
    iu = np.triu_indices(15, 1)
    res = {"n_links": int(sick.shape[0])}

    chk = {}
    for g, ranks in (("sick", sick), ("healthy", healthy)):
        mine = np.stack([np.sort(maps[g][s][iu])[::-1] for s in range(4)], 1)
        chk[g] = float(np.nanmax(np.abs(mine - ranks)))
    res["ranks_vs_stage_cells_max_abs_diff"] = chk

    between = {}
    for s, st in enumerate(STAGES):
        w = wilcoxon(sick[:, s], healthy[:, s])  # thesis: signrank on rank-sorted columns
        k = ks_2samp(sick[:, s], healthy[:, s])
        a, b = maps["sick"][s][iu], maps["healthy"][s][iu]
        ok = np.isfinite(a) & np.isfinite(b)
        wl = wilcoxon(a[ok], b[ok])
        between[st] = {"thesis_signrank_p": float(w.pvalue), "thesis_ks_p": float(k.pvalue),
                       "median_sick_pct": float(100 * np.median(sick[:, s])),
                       "median_healthy_pct": float(100 * np.median(healthy[:, s])),
                       "supp_link_matched_signrank_p": float(wl.pvalue),
                       "supp_link_matched_frac_sick_gt_healthy": float(np.mean(a[ok] > b[ok]))}
    res["sick_vs_healthy"] = between

    within = {}
    for g, ranks in (("sick", sick), ("healthy", healthy)):
        within[g] = {f"{STAGES[i]}|{STAGES[j]}": {"signrank_p": float(wilcoxon(ranks[:, i], ranks[:, j]).pvalue),
                                                   "ks_p": float(ks_2samp(ranks[:, i], ranks[:, j]).pvalue)}
                     for i in range(4) for j in range(i + 1, 4)}
    res["between_stages"] = within

    # full-night per-patient networks (stored TDS_mats) vs AHI — NOT the thesis's per-stage Fig 27
    with h5py.File(f"{d}/TDS_mats.mat") as f:
        R = f["#refs#"]
        mats = names = None
        for k in R:
            o = R[k]
            if isinstance(o, h5py.Dataset) and o.attrs.get("MATLAB_class") == b"cell" and o.shape == (1, 28):
                els = [f[r] for r in o[0]]
                if els[0].attrs.get("MATLAB_class") == b"double" and els[0].shape == (15, 15):
                    mats = np.stack([e[()].T for e in els])
                elif els[0].attrs.get("MATLAB_class") == b"char":
                    s = ["".join(map(chr, e[()].ravel())) for e in els]
                    if "AHI" in s[0]:
                        names = s
    ahi = np.array([float(re.search(r"AHI_([\d.]+)\.mat", n).group(1)) for n in names])
    thr = 7.8
    deg_sum = np.array([int(((m[iu] > thr)).sum()) for m in mats])  # number of links = sum of degrees / 2
    rho = spearmanr(ahi, deg_sum)
    res["full_night"] = {"n_patients": int(len(ahi)), "n_ahi_gt_10": int((ahi > 10).sum()),
                         "n_ahi_le_10": int((ahi <= 10).sum()), "threshold_pct": thr,
                         "spearman_links_vs_ahi": float(rho.statistic), "spearman_p": float(rho.pvalue),
                         "median_links_sick": float(np.median(deg_sum[ahi > 10])),
                         "median_links_healthy": float(np.median(deg_sum[ahi <= 10])),
                         "per_patient_ahi_links": [[float(a), int(n)] for a, n in sorted(zip(ahi, deg_sum))]}
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(res, open(OUT, "w"), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k != "between_stages"}, indent=1)[:3000])


if __name__ == "__main__":
    main(sys.argv[1])
