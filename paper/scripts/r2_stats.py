"""R2 statistics (reads only paper/results/r2_scores.csv, r2_params.csv, r2_null.npz). See ANALYSIS_PLAN.md."""
from __future__ import annotations

import argparse
import json
import os
import warnings

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
STAGES = ("LS", "DS", "REM", "awake")
MAIN = ("published", "win_only", "tol_only", "v2_calibrated", "v2_cal_lag3")
VARIANTS = MAIN + ("v2_cap120", "v2_cap600")
RULES = ("strict", "p80")
MEASURES = ("excess", "excess_all45", "real", "real_all45")
HYP = {"H1": ("LS", "DS"), "H2": ("awake", "DS"), "H3": ("LS", "REM"), "H4": ("awake", "REM")}
NBOOT = 10000


def holm(p):
    p = np.asarray(p, float)
    out = np.full(len(p), np.nan)
    ok = np.flatnonzero(np.isfinite(p))
    order = ok[np.argsort(p[ok])]
    run, m = 0.0, len(order)
    for r, i in enumerate(order):
        run = max(run, (m - r) * p[i])
        out[i] = min(1.0, run)
    return out


def dz_of(d):
    return float(np.mean(d) / np.std(d, ddof=1)) if len(d) > 2 and np.std(d, ddof=1) > 0 else np.nan


def boot_ci(d, rng):
    """95 % CI of mean and dz, subject bootstrap."""
    n = len(d)
    ix = rng.integers(0, n, (NBOOT, n))
    b = d[ix]
    m = b.mean(1)
    sd = b.std(1, ddof=1)
    with np.errstate(all="ignore"):
        z = np.where(sd > 0, m / sd, np.nan)
    return np.nanpercentile(m, [2.5, 97.5]), np.nanpercentile(z, [2.5, 97.5])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--resdir", default=os.path.join(HERE, "..", "results"))
    a = ap.parse_args()
    R = a.resdir
    sc = pd.read_csv(os.path.join(R, "r2_scores.csv"))
    pr = pd.read_csv(os.path.join(R, "r2_params.csv"))
    nz = np.load(os.path.join(R, "r2_null.npz"))
    subs = sorted(sc.subject.unique())
    ns = len(subs)

    # connectivity[measure][(variant, rule)] -> DataFrame subject x stage
    def agg(df, col):
        g = df.groupby(["subject", "variant", "rule", "stage"])[col].apply(lambda s: s.mean() if s.notna().any() else np.nan)
        return g
    thesis = sc[sc.thesis_link == 1]
    conn = {"excess": agg(thesis, "excess"), "real": agg(thesis, "real"),
            "excess_all45": agg(sc, "excess"), "real_all45": agg(sc, "real")}

    def get(meas, v, rule):
        g = conn[meas]
        return g.xs((v, rule), level=("variant", "rule")).unstack("stage").reindex(index=subs, columns=list(STAGES))

    # ---- n defined per stage x variant x rule
    ndef = []
    for v in VARIANTS:
        for rule in RULES:
            c = get("excess", v, rule)
            for s in STAGES:
                ndef.append({"variant": v, "rule": rule, "stage": s, "n_defined": int(c[s].notna().sum())})
    ndef = pd.DataFrame(ndef)

    # ---- H1-H4
    rng = np.random.default_rng(0)
    eff = []
    diffs = {}
    for meas in MEASURES:
        for v in VARIANTS:
            for rule in RULES:
                c = get(meas, v, rule)
                rows = []
                for h, (s1, s2) in HYP.items():
                    d = (c[s1] - c[s2]).dropna()
                    diffs[(meas, v, rule, h)] = d
                    n = len(d)
                    r = {"measure": meas, "variant": v, "rule": rule, "hyp": h, "contrast": f"{s1}>{s2}", "n": n,
                         "mean_diff": np.nan, "sd": np.nan, "dz": np.nan, "mean_ci_lo": np.nan, "mean_ci_hi": np.nan,
                         "dz_ci_lo": np.nan, "dz_ci_hi": np.nan, "wilcoxon_p": np.nan}
                    if n >= 3:
                        dv = d.values
                        r["mean_diff"], r["sd"], r["dz"] = float(dv.mean()), float(dv.std(ddof=1)), dz_of(dv)
                        (r["mean_ci_lo"], r["mean_ci_hi"]), (r["dz_ci_lo"], r["dz_ci_hi"]) = boot_ci(dv, rng)
                        if np.any(dv != 0):
                            r["wilcoxon_p"] = float(stats.wilcoxon(dv, alternative="two-sided").pvalue)
                    rows.append(r)
                hp = holm([r["wilcoxon_p"] for r in rows])
                for r, p in zip(rows, hp):
                    r["p_holm"] = p
                eff += rows
    eff = pd.DataFrame(eff)
    eff.to_csv(os.path.join(R, "r2_effects.csv"), index=False)

    # ---- H5 gain
    gain = []
    rng = np.random.default_rng(1)
    for meas in MEASURES:
        for rule in RULES:
            for v in VARIANTS:
                if v == "published":
                    continue
                for h in HYP:
                    dp, dv = diffs[(meas, "published", rule, h)], diffs[(meas, v, rule, h)]
                    both = dp.index.intersection(dv.index)
                    r = {"measure": meas, "rule": rule, "variant": v, "hyp": h, "n_both": len(both),
                         "dz_pub": np.nan, "dz_var": np.nan, "gain": np.nan, "gain_ci_lo": np.nan, "gain_ci_hi": np.nan}
                    if len(both) >= 3:
                        x, y = dp[both].values, dv[both].values
                        r["dz_pub"], r["dz_var"] = dz_of(x), dz_of(y)
                        r["gain"] = r["dz_var"] - r["dz_pub"]
                        ix = rng.integers(0, len(both), (NBOOT, len(both)))
                        bx, by = x[ix], y[ix]
                        with np.errstate(all="ignore"):
                            g = by.mean(1) / by.std(1, ddof=1) - bx.mean(1) / bx.std(1, ddof=1)
                        r["gain_ci_lo"], r["gain_ci_hi"] = np.nanpercentile(g, [2.5, 97.5])
                    gain.append(r)
    gain = pd.DataFrame(gain)
    gain.to_csv(os.path.join(R, "r2_gain.csv"), index=False)

    def decide(g):
        pos = int((g.gain_ci_lo > 0).sum())
        neg = int((g.gain_ci_hi < 0).sum())
        dec = "supported" if pos >= 3 else ("refuted" if neg >= 2 else "inconclusive")
        return {"decision": dec, "n_contrasts_ci_above_0": pos, "n_contrasts_ci_below_0": neg}
    h5 = {}
    for meas in MEASURES:
        for rule in RULES:
            for v in VARIANTS:
                if v == "published":
                    continue
                g = gain[(gain.measure == meas) & (gain.rule == rule) & (gain.variant == v)]
                h5[f"{meas}|{rule}|{v}"] = decide(g)

    # ---- network
    net = []
    for v in MAIN:
        A = nz[v]  # [subj, link, stage, draw]
        subj_all = list(nz["subjects"])
        links = list(nz["link_order"])
        s_v = sc[(sc.variant == v) & (sc.rule == "strict") & (sc.thesis_link == 1)]
        for li, lk in enumerate(links):
            for si, s in enumerate(STAGES):
                rr = s_v[(s_v.pair == lk) & (s_v.stage == s)].set_index("subject").real.reindex(subj_all).values
                real = rr[np.isfinite(rr)]
                nul = A[:, li, si, :].ravel()
                nul = nul[np.isfinite(nul)]
                auc = p = np.nan
                if len(real) >= 3 and len(nul) >= 10:
                    u = stats.mannwhitneyu(real, nul, alternative="greater", method="asymptotic")
                    auc, p = float(u.statistic / (len(real) * len(nul))), float(u.pvalue)
                net.append({"variant": v, "link": lk, "stage": s, "n_real": len(real), "n_null": len(nul),
                            "auc": auc, "p": p, "significant": bool(np.isfinite(p) and p < 1e-3)})
    net = pd.DataFrame(net)
    net.to_csv(os.path.join(R, "r2_network.csv"), index=False)
    netsum = {}
    for v in MAIN:
        n = net[net.variant == v]
        netsum[v] = {"mean_auc_all": float(n.auc.mean()), "n_cells_defined": int(n.auc.notna().sum()),
                     "n_significant_total": int(n.significant.sum()),
                     "per_stage": {s: {"n_significant": int(n[n.stage == s].significant.sum()),
                                       "mean_auc": float(n[n.stage == s].auc.mean()),
                                       "n_cells_defined": int(n[n.stage == s].auc.notna().sum())} for s in STAGES}}
    a_pub = net[net.variant == "published"].set_index(["link", "stage"]).auc
    a_v2 = net[net.variant == "v2_calibrated"].set_index(["link", "stage"]).auc
    both = a_pub.notna() & a_v2.notna()
    h6 = {"mean_auc_published_all": float(a_pub.mean()), "mean_auc_v2_calibrated_all": float(a_v2.mean()),
          "mean_auc_published_common_cells": float(a_pub[both].mean()),
          "mean_auc_v2_calibrated_common_cells": float(a_v2[both].mean()),
          "n_common_cells": int(both.sum()),
          "v2_ge_published": bool(a_v2[both].mean() >= a_pub[both].mean()) if both.any() else None}

    # ---- params
    psum = {}
    for v in VARIANTS:
        p = pr[pr.variant == v]
        psum[v] = {k: {"median": float(p[k].median()), "min": float(p[k].min()), "max": float(p[k].max())}
                   for k in ("window", "tolerance", "max_lag")}
        psum[v]["calibration_failed_count"] = int(p.calibration_failed.fillna(0).astype(bool).sum())
        psum[v]["window_clipped_count"] = int(p.window_clipped.fillna(0).astype(bool).sum())
    bart = {c[len("raw_bartlett_"):]: float(pr[pr.variant == "published"][c].median())
            for c in pr.columns if c.startswith("raw_bartlett_")}

    # ---- NaN accounting
    nan_ct = {}
    for meas in ("excess", "real"):
        for v in VARIANTS:
            for rule in RULES:
                c = get(meas, v, rule)
                nan_ct[f"{meas}|{v}|{rule}"] = {s: int(c[s].isna().sum()) for s in STAGES}

    out = {"n_subjects": ns, "H5_decision_primary(excess,strict,v2_calibrated)": h5["excess|strict|v2_calibrated"],
           "H5_all": h5, "H6": h6, "network_summary": netsum,
           "n_defined": {f"{r.variant}|{r.rule}|{r.stage}": r.n_defined for r in ndef.itertuples()},
           "params_summary": psum, "median_raw_bartlett_window_by_signal": bart,
           "n_subjects_undefined_connectivity": nan_ct}
    with open(os.path.join(R, "r2_summary.json"), "w") as f:
        json.dump(out, f, indent=1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for v in ("published", "v2_calibrated"):
            print(eff[(eff.variant == v) & (eff.rule == "strict") & (eff.measure == "excess")].round(4).to_string(index=False))
        print(gain[(gain.variant == "v2_calibrated") & (gain.rule == "strict") & (gain.measure == "excess")].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
