"""R1.3: thesis significance threshold, ported from msc-project '01 repreduce NCOM/surrogate test/looking_for_threshold.m'
(+ pooled stage cells of '01 repreduce NCOM/makeDataForReconFig2.m', section 'Divide to sleep stages').

 (a) pooled mean_real (140, order DS,LS,REM,awake x 35 links): stable labels of all 35 subjects pooled per link x stage
     (thesis_port.links_counts, centre-sample stage rule), x100; compared to the stored mean_real.
 (b) MATLAB ttest2(real_links, sur_links) = equal-variance two-sample t per column (NaNs omitted), significant if p < 1e-3.
 (c) signi_over(thres) for thres = 1:0.1:20 and thres = first thres with signi_over == 100.
Run with the stored sur_links and with our ported surrogates (r1_sur_links_port.npy); compares the two surrogate
distributions per column (mean, 95th pct, two-sample KS statistic).

Usage: python r1_threshold.py DATA.npz REPO [SUR_PORT.npy]
Writes ../results/r1_threshold.json and ../results/r1_surrogate_compare.csv
"""
import csv, json, os, sys
from concurrent.futures import ProcessPoolExecutor
import numpy as np, scipy.io as sio
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__))
from thesis_port import links_counts, links_strength, STAGES

DATA, REPO = sys.argv[1], sys.argv[2]
SUR = sys.argv[3] if len(sys.argv) > 3 else os.path.join(os.path.dirname(__file__), "..", "results", "r1_sur_links_port.npy")
OUT = os.path.join(os.path.dirname(__file__), "..", "results")
THRES = 1 + 0.1 * np.arange(191)                      # 1:0.1:20


def one(i):
    z = np.load(DATA); n = z["names"][i]
    X, S = z["X_" + n], z["S_" + n]
    return links_counts(X, S), links_strength(X, S).flatten(order="F")


def ttest(real, sur):
    return stats.ttest_ind(real, sur, equal_var=True, nan_policy="omit")[1]


def threshold(mean_real, sig):
    over = np.array([100 * np.sum(mean_real > t) / 140 for t in THRES])
    with np.errstate(invalid="ignore", divide="ignore"):
        so = np.array([100 * np.sum((mean_real > t) & sig) / np.sum(mean_real > t) for t in THRES])
    k = np.flatnonzero(so == 100)
    return so, (float(THRES[k[0]]) if len(k) else None)


def analyse(name, real, sur, mean_real, ref):
    p = ttest(real, sur)
    sig = p < 1e-3
    so, thr = threshold(mean_real, sig)
    r = {"surrogates": name, "n_sur": int(sur.shape[0]), "n_significant": int(sig.sum()), "thres": thr,
         "n_links_over_thres": None if thr is None else int((mean_real > thr).sum())}
    if ref is not None:
        r["max_abs_diff_p"] = float(np.max(np.abs(p - ref["p"])))
        r["max_rel_diff_p"] = float(np.max(np.abs(p - ref["p"])[ref["p"] > 0] / ref["p"][ref["p"] > 0]))
        r["significant_links_mismatch"] = int((sig != ref["sig"]).sum())
        r["max_abs_diff_signi_over"] = float(np.nanmax(np.abs(so - ref["so"])))
        r["signi_over_nan_equal"] = bool(np.array_equal(np.isnan(so), np.isnan(ref["so"])))
        r["thres_stored"] = ref["thres"]
        r["thres_diff"] = None if thr is None else abs(thr - ref["thres"])
    return r, p, sig, so


def stored_stage_cells_mean_real():
    """mean_real as looking_for_threshold.m builds it: from the thesis's own 01 repreduce NCOM/Stages_cell.mat
    (pooled stable labels written by makeDataForReconFig2.m from the raw .dat files and their own stage files),
    stage order DS, LS, REM, awake, links in maps(good_links') order."""
    sc = sio.loadmat(f"{REPO}/01 repreduce NCOM/Stages_cell.mat")["stages_cell"]
    g = np.triu(np.ones((10, 10)), 1).astype(bool)
    g[5:, 5:] = False
    order = [(r, c) for c in range(10) for r in range(10) if g.T[r, c]]
    out = []
    for k in range(4):
        cell = sc[k, 0]
        out += [100 * np.asarray(cell[r, c], float).ravel().mean() for r, c in order]
    return np.array(out)


if __name__ == "__main__":
    m = sio.loadmat(f"{REPO}/01 repreduce NCOM/links_strength2000Sur.mat")
    real_st, sur_st, mean_st = m["real_links"], m["sur_links"], m["mean_real"].ravel()
    ref = {"p": m["p"].ravel(), "sig": m["significant_links"].ravel().astype(bool), "so": m["signi_over"].ravel(),
           "thres": float(m["thres"].ravel()[0])}
    with ProcessPoolExecutor(4) as ex:
        res = list(ex.map(one, range(35)))
    cnt = np.sum([r[0] for r in res], axis=0)                     # 35 x 4 x 2 pooled
    real = np.array([r[1] for r in res])
    pooled = 100 * cnt[:, :, 0] / cnt[:, :, 1]                    # 35 links x 4 stages (STAGES order = DS,LS,REM,awake)
    mean_real = pooled.flatten(order="F")
    out = {"pooled_mean_real": {"max_abs_diff_vs_stored": float(np.max(np.abs(mean_real - mean_st))),
                                 "stages_order": list(STAGES), "mean_real": mean_real.tolist()},
           "real_links_max_abs_diff_vs_stored": float(np.nanmax(np.abs(real - real_st))),
           "mean_of_real_links_vs_pooled_max_abs_diff": float(np.max(np.abs(real.mean(0) - mean_real)))}
    r, p_st, sig_st, so_st = analyse("stored sur_links", real, sur_st, mean_real, ref)
    out["stored_surrogates"] = r
    out["stored_surrogates"]["p"] = p_st.tolist()
    out["stored_surrogates"]["signi_over"] = so_st.tolist()
    # diagnostics: which input explains the difference to the stored thres? (stored mean_real and/or stored p)
    so_ss, thr_ss = threshold(mean_st, ref["sig"])
    out["diagnostic_stored_mean_real_stored_p"] = {"thres": thr_ss, "max_abs_diff_signi_over": float(np.nanmax(np.abs(so_ss - ref["so"])))}
    out["diagnostic_stored_p_recomputed_from_stored_arrays"] = {
        "max_abs_diff_p": float(np.max(np.abs(p_st - ref["p"]))), "n_significant_stored_p": int(ref["sig"].sum()),
        "n_significant_recomputed": int(sig_st.sum()), "significant_mismatch": int((sig_st != ref["sig"]).sum())}
    out["diagnostic_stored_mean_real_recomputed_sig"] = {"thres": threshold(mean_st, sig_st)[1]}
    out["thres_grid"] = THRES.tolist()
    out["thres_thesis_text"] = 9
    # The thesis's mean_real came from its own Stages_cell.mat (raw-.dat pipeline), not from All_data_clean.
    mr_cells = stored_stage_cells_mean_real()
    out["thesis_stage_cells_mean_real"] = {
        "max_abs_diff_vs_stored_mean_real": float(np.nanmax(np.abs(mr_cells - mean_st))),
        "thres_with_stored_sig": threshold(mr_cells, ref["sig"])[1],
        "thres_with_sig_recomputed_from_stored_sur": threshold(mr_cells, sig_st)[1]}
    if os.path.exists(SUR):
        sur = np.load(SUR)
        r, p_p, sig_p, so_p = analyse("ported surrogates", real, sur, mean_real, ref)
        r["p"] = p_p.tolist(); r["signi_over"] = so_p.tolist()
        r["significant_links_mismatch_vs_stored"] = int((sig_p != sig_st).sum())
        out["port_surrogates"] = r
        out["diagnostic_stored_mean_real_port_sig"] = {"thres": threshold(mean_st, sig_p)[1]}
        out["thesis_stage_cells_mean_real"]["thres_with_sig_from_port_surrogates"] = threshold(mr_cells, sig_p)[1]
        rows, S = [], ["stage", "link"]
        for c in range(140):
            a, b = sur[:, c], sur_st[:, c]
            ks = stats.ks_2samp(a, b).statistic
            rows.append([STAGES[c // 35], c % 35, c, a.mean(), b.mean(), np.percentile(a, 95), np.percentile(b, 95), ks])
        arr = np.array([[x[3], x[4], x[5], x[6], x[7]] for x in rows])
        d = lambda i, j: np.abs(arr[:, i] - arr[:, j])
        out["surrogate_compare"] = {
            "mean_abs_diff_of_means": {"median": float(np.median(d(0, 1))), "max": float(d(0, 1).max())},
            "abs_diff_p95": {"median": float(np.median(d(2, 3))), "max": float(d(2, 3).max())},
            "ks_statistic": {"median": float(np.median(arr[:, 4])), "max": float(arr[:, 4].max())},
            "mean_of_means_ours": float(arr[:, 0].mean()), "mean_of_means_stored": float(arr[:, 1].mean())}
        with open(os.path.join(OUT, "r1_surrogate_compare.csv"), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["stage", "link_idx", "col", "mean_ours", "mean_stored", "p95_ours", "p95_stored", "ks"])
            w.writerows(rows)
    json.dump(out, open(os.path.join(OUT, "r1_threshold.json"), "w"), indent=1)
    for k in ("pooled_mean_real", "real_links_max_abs_diff_vs_stored", "mean_of_real_links_vs_pooled_max_abs_diff"):
        print(k, out[k] if k != "pooled_mean_real" else out[k]["max_abs_diff_vs_stored"])
    for k in ("stored_surrogates", "port_surrogates"):
        if k in out:
            print(k, {a: b for a, b in out[k].items() if a not in ("p", "signi_over")})
    print(out.get("surrogate_compare"))
