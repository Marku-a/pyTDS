"""Part A: aTDS v2.1 vs published TDS parameters on the NCOM overnight cohort (35 subjects, 10 channels, 1 Hz).

Usage: python ncom_cohort.py <All_data_clean.mat> [<out_dir>] [--n-subjects N] [--n-null N]
Private data are read from the given path only; only derived numbers are written.
"""
from __future__ import annotations

import csv
import json
import os
import sys
import time
from dataclasses import replace
from itertools import combinations
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import atds21 as A  # noqa: E402
from pyTDS.core import stable_label, time_delay_interaction  # noqa: E402
from pyTDS.params import TDSParams  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CH = ["HR", "Resp", "Chin", "Leg", "Eye", "delta", "theta", "alpha", "sigma", "beta"]
STAGES = ["W", "LS", "DS", "REM"]
STAGE_MAP = {"awake": 0, "LS": 1, "DS": 2, "REM": 3}
SETTINGS = ["published", "atds", "atds_cap300"]
PAIRS = list(combinations(range(10), 2))
CONTRASTS = [("LS", "DS"), ("W", "DS"), ("REM", "DS"), ("LS", "REM")]
FRAC = 0.8
N_NULL = 100
N_BOOT = 5000
SEED = 0
_G: dict = {}


# ------------------------------------------------------------------ data
def _label(e):
    e = np.asarray(e).ravel()
    return str(e[0]) if len(e) else ""


def load_cohort(path, n_subjects=None):
    import matio
    r = matio.load_from_mat(path)["p_data_stucrt"][0, 0]
    ids = [str(np.asarray(n).ravel()[0]) for n in r[0].ravel()]
    sigs, stages, nan_counts = [], [], {}
    for pid, df, st in zip(ids, r[1].ravel(), r[3].ravel()):
        x = df[CH].values.astype(float)
        nan_counts[pid] = {c: int(np.isnan(x[:, i]).sum()) for i, c in enumerate(CH)}
        sigs.append(np.column_stack([A.interp_nan(x[:, i]) for i in range(10)]))
        lab = [_label(e) for e in st.ravel()]
        assert len(lab) == len(x), (pid, len(lab), len(x))
        stages.append(np.array([STAGE_MAP.get(s, -1) for s in lab], dtype=int))
    if n_subjects:
        ids, sigs, stages = ids[:n_subjects], sigs[:n_subjects], stages[:n_subjects]
    return ids, sigs, stages, nan_counts


# ------------------------------------------------------------------ parameters
def cohort_params(sigs, calib_seed=SEED):
    """One parameter set per setting for the whole cohort."""
    B = np.array([[A.bartlett_factor(x[:, i]) for i in range(10)] for x in sigs])
    bmax = B.max(axis=1)
    med_len = float(np.median([len(x) for x in sigs]))
    L0 = int(np.clip(round(30 * float(np.median(bmax))), 61, 0.1 * med_len))
    rng = np.random.default_rng(calib_seed)
    combos = [(s, k) for s in range(len(sigs)) for k in range(len(PAIRS))]
    sel = rng.choice(len(combos), min(50, len(combos)), replace=False)
    cal_pairs = [(sigs[combos[c][0]][:, PAIRS[combos[c][1]][0]], sigs[combos[c][0]][:, PAIRS[combos[c][1]][1]])
                 for c in sel]
    params, info = {"published": A.PUBLISHED}, {"published": {"calibration_failed": False}}
    for name, L in (("atds", L0), ("atds_cap300", min(L0, 300))):
        p0 = TDSParams(window=L, overlap=L // 2, max_lag=(L - 1) // 2, tolerance=1)
        p, fail = A.calib(p0, cal_pairs)
        params[name] = p
        info[name] = {"calibration_failed": bool(fail)}
    summ = {"B_max_per_subject_median": float(np.median(bmax)), "B_max_min": float(bmax.min()),
            "B_max_max": float(bmax.max()), "B_per_channel_median": dict(zip(CH, np.median(B, axis=0).tolist())),
            "median_length": med_len, "uncapped_window": L0}
    for n, p in params.items():
        info[n].update(window=p.window, step=p.overlap, max_lag=p.max_lag, tolerance=p.tolerance, B_summary=summ)
    return params, info


# ------------------------------------------------------------------ stage assignment / per-stage TDS
def assign_windows(stage, p, nwin):
    """Stage code of each window (>= 80 % of its seconds in one stage) else -1."""
    out = np.full(nwin, -1)
    starts = np.arange(nwin) * p.overlap
    for s in range(4):
        cs = np.concatenate([[0], np.cumsum(stage == s)])
        ok = (cs[starts + p.window] - cs[starts]) >= FRAC * p.window
        out[ok] = s
    return out


def stage_tds(lab, asg):
    r = np.full(4, np.nan)
    for s in range(4):
        m = asg == s
        if m.sum() >= 3:
            r[s] = 100.0 * lab[m].sum() / m.sum()
    return r


def pair_stage_tds(x, y, p, stage):
    tau = time_delay_interaction(x, y, p)[0]
    lab = stable_label(tau, p)
    return stage_tds(lab, assign_windows(stage, p, len(tau))), len(tau)


# ------------------------------------------------------------------ workers
def _subject_worker(s):
    sig, stage = _G["sigs"][s], _G["stages"][s]
    res = np.full((3, 45, 4), np.nan)
    nass = np.zeros((3, 4), dtype=int)
    for si, name in enumerate(SETTINGS):
        p = _G["params"][name]
        for k, (i, j) in enumerate(PAIRS):
            res[si, k], nw = pair_stage_tds(sig[:, i], sig[:, j], p, stage)
        asg = assign_windows(stage, p, nw)
        nass[si] = [(asg == q).sum() for q in range(4)]
    return res, nass


def _null_worker(k):
    sigs, stages, params = _G["sigs"], _G["stages"], _G["params"]
    i, j = PAIRS[k]
    rng = np.random.default_rng([SEED, k])
    S = len(sigs)
    out = np.full((_G["n_null"], 3, 4), np.nan)
    for d in range(_G["n_null"]):
        a = int(rng.integers(S)); b = int(rng.integers(S - 1)); b += b >= a
        n = min(len(sigs[a]), len(sigs[b]))
        x, y, st = sigs[a][:n, i], sigs[b][:n, j], stages[a][:n]
        for si, name in enumerate(SETTINGS):
            out[d, si] = pair_stage_tds(x, y, params[name], st)[0]
    with np.errstate(all="ignore"):
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return np.nanpercentile(out, 95, axis=0), np.nanmean(out, axis=0)


# ------------------------------------------------------------------ statistics
def dz(d):
    return float(np.mean(d) / np.std(d, ddof=1))


def boot_effect(d, rng_seed=SEED):
    d = np.asarray(d, float)
    n = len(d)
    idx = np.random.default_rng(rng_seed).integers(0, n, (N_BOOT, n))
    bd = d[idx]
    b = bd.mean(1) / bd.std(1, ddof=1)
    return {"n": int(n), "dz": dz(d), "ci95": [float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))],
            "frac_positive": float(np.mean(d > 0))}


def boot_gain(da, dp, rng_seed=SEED):
    da, dp = np.asarray(da, float), np.asarray(dp, float)
    n = len(da)
    idx = np.random.default_rng(rng_seed).integers(0, n, (N_BOOT, n))
    ba, bp = da[idx], dp[idx]
    b = ba.mean(1) / ba.std(1, ddof=1) - bp.mean(1) / bp.std(1, ddof=1)
    return {"n": int(n), "gain": dz(da) - dz(dp), "ci95": [float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))]}


def effects(mean_tds, nsig):
    """mean_tds, nsig: arrays (settings, subjects, stages)."""
    out = {}
    for metric, arr in (("mean_link_tds", mean_tds), ("n_significant_links", nsig)):
        out[metric] = {}
        for si, name in enumerate(SETTINGS):
            out[metric][name] = {}
            for a, b in CONTRASTS:
                d = arr[si][:, STAGES.index(a)] - arr[si][:, STAGES.index(b)]
                d = d[np.isfinite(d)]
                out[metric][name][f"{a}-{b}"] = boot_effect(d)
        out[metric]["gain"] = {}
        for name in ("atds", "atds_cap300"):
            out[metric]["gain"][name] = {}
            for a, b in CONTRASTS:
                ia, ib = STAGES.index(a), STAGES.index(b)
                d1 = arr[SETTINGS.index(name)][:, ia] - arr[SETTINGS.index(name)][:, ib]
                d0 = arr[0][:, ia] - arr[0][:, ib]
                ok = np.isfinite(d1) & np.isfinite(d0)
                out[metric]["gain"][name][f"{a}-{b}"] = boot_gain(d1[ok], d0[ok])
    return out


# ------------------------------------------------------------------ main
def main(path, out_dir, n_subjects=None, n_null=N_NULL):
    t0 = time.time()
    os.makedirs(out_dir, exist_ok=True)
    ids, sigs, stages, nan_counts = load_cohort(path, n_subjects)
    S = len(ids)
    params, pinfo = cohort_params(sigs)
    print("params", {k: (p.window, p.overlap, p.max_lag, p.tolerance) for k, p in params.items()}, flush=True)
    _G.update(sigs=sigs, stages=stages, params=params, n_null=n_null)
    with Pool(A.WORKERS) as pool:
        real = pool.map(_subject_worker, range(S))
        print(f"real done {time.time() - t0:.0f}s", flush=True)
        nulls = pool.map(_null_worker, range(len(PAIRS)))
    print(f"null done {time.time() - t0:.0f}s", flush=True)
    R = np.stack([r[0] for r in real])             # (S, 3, 45, 4)
    NA = np.stack([r[1] for r in real])            # (S, 3, 4)
    TH = np.stack([n[0] for n in nulls], axis=2)   # (3, 4, 45)
    NM = np.stack([n[1] for n in nulls], axis=2)
    sig = np.zeros_like(R, dtype=bool)
    with np.errstate(invalid="ignore"):
        sig = R > TH.transpose(0, 2, 1)[None]      # (S, 3, 45, 4)
    sig &= np.isfinite(R)
    mean_tds = np.full((3, S, 4), np.nan)
    nsig = np.zeros((3, S, 4))
    for si in range(3):
        for q in range(4):
            for s in range(S):
                v = R[s, si, :, q]
                if np.isfinite(v).any():
                    mean_tds[si, s, q] = np.nanmean(v)
                    nsig[si, s, q] = sig[s, si, :, q].sum()
                else:
                    nsig[si, s, q] = np.nan  # stage not scorable (< 3 assigned windows)

    rows = [dict(subject=ids[s], stage=STAGES[q], setting=SETTINGS[si], mean_link_tds=mean_tds[si, s, q],
                 n_significant_links=nsig[si, s, q], n_assigned_windows=int(NA[s, si, q]))
            for s in range(S) for si in range(3) for q in range(4)]
    matrices = {}
    mrows = []
    for si, name in enumerate(SETTINGS):
        matrices[name] = {}
        for q, st in enumerate(STAGES):
            M = np.full((10, 10), np.nan); F = np.full((10, 10), np.nan)
            for k, (i, j) in enumerate(PAIRS):
                v = R[:, si, k, q]; ok = np.isfinite(v)
                if ok.any():
                    M[i, j] = M[j, i] = v[ok].mean()
                    F[i, j] = F[j, i] = sig[:, si, k, q][ok].mean()
                    mrows.append(dict(setting=name, stage=st, ch_i=CH[i], ch_j=CH[j], mean_tds=M[i, j],
                                      frac_subjects_significant=F[i, j], n_subjects=int(ok.sum()),
                                      null95=TH[si, q, k], null_mean=NM[si, q, k]))
            matrices[name][st] = {"mean_tds": M.tolist(), "frac_significant": F.tolist()}
    eff = effects(mean_tds, nsig)
    out = {"channels": CH, "subjects": ids, "n_subjects": S, "n_null_draws": n_null,
           "params": pinfo, "nan_counts": nan_counts, "matrices": matrices, "effects": eff,
           "subject_stage": rows, "wall_seconds": time.time() - t0}
    with open(os.path.join(out_dir, "ncom_results.json"), "w") as f:
        json.dump(out, f, indent=1)

    def wcsv(name, rs):
        with open(os.path.join(out_dir, name), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rs[0])); w.writeheader(); w.writerows(rs)
    wcsv("ncom_subject_stage.csv", rows)
    wcsv("ncom_link_matrices.csv", mrows)
    erows = []
    for metric, d in eff.items():
        for name, cd in d.items():
            if name == "gain":
                for n2, c2 in cd.items():
                    for con, e in c2.items():
                        erows.append(dict(metric=metric, setting=f"gain_{n2}_minus_published", contrast=con, n=e["n"],
                                          value=e["gain"], ci_lo=e["ci95"][0], ci_hi=e["ci95"][1], frac_positive=""))
            else:
                for con, e in cd.items():
                    erows.append(dict(metric=metric, setting=name, contrast=con, n=e["n"], value=e["dz"],
                                      ci_lo=e["ci95"][0], ci_hi=e["ci95"][1], frac_positive=e["frac_positive"]))
    wcsv("ncom_effects.csv", erows)
    wcsv("ncom_params.csv", [dict(setting=k, **{a: b for a, b in v.items() if a != "B_summary"}) for k, v in pinfo.items()])
    for r in erows:
        print(r["metric"], r["setting"], r["contrast"], r["n"], round(r["value"], 3), round(r["ci_lo"], 3), round(r["ci_hi"], 3), r["frac_positive"])
    print(f"total {time.time() - t0:.0f}s")


if __name__ == "__main__":
    args = sys.argv[1:]
    opt = {}
    for flag in ("--n-subjects", "--n-null"):
        if flag in args:
            k = args.index(flag); opt[flag] = int(args[k + 1]); del args[k:k + 2]
    main(args[0], args[1] if len(args) > 1 else os.path.join(HERE, "results"),
         opt.get("--n-subjects"), opt.get("--n-null", N_NULL))
