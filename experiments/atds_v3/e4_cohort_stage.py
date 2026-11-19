"""E4: stage-aware aTDS on the 35-subject sleep cohort (reuses experiments/atds_real/ncom_cohort.py machinery).

Usage: ATDS_CODE=<checkout> python e4_cohort_stage.py <All_data_clean.mat> [<out_dir>] [--n-subjects N] [--n-null N]
Private data are read from the given path only; only aggregates are written.

Settings: published, atds_cap300 (as in ncom_cohort.py) and atds_stage: window cap L = clip(round(0.5 * median stage-episode
length), 61, 300), lag (L-1)//2, step L//2, and one tolerance per stage chosen on 50 cross-subject stage-matched surrogate pairs.
"""
import json
import os
import sys
import time
import warnings
from dataclasses import replace
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402
import ncom_cohort as NC  # noqa: E402

SETTINGS = ["published", "atds_cap300", "atds_stage"]
N_SURR = 50
_orig_pair_stage_tds = NC.pair_stage_tds


class StageParams:
    """TDS parameters with one tolerance per stage (index = stage code 0..3 = W, LS, DS, REM)."""

    def __init__(self, base, tols):
        self.base, self.tols = base, list(tols)
        self.window, self.overlap, self.max_lag, self.tolerance = base.window, base.overlap, base.max_lag, base.tolerance


def pair_stage_tds(x, y, p, stage):
    if not isinstance(p, StageParams):
        return _orig_pair_stage_tds(x, y, p, stage)
    tau = C.time_delay_interaction(x, y, p.base)[0]
    asg = NC.assign_windows(stage, p.base, len(tau))
    r = np.full(4, np.nan)
    for s in range(4):
        m = asg == s
        if m.sum() >= 3:
            r[s] = 100.0 * C.stable_label(tau, replace(p.base, tolerance=p.tols[s]))[m].mean()
    return r, len(tau)


# ------------------------------------------------------------------ episode lengths
def episodes(stages):
    """Lengths (samples) of runs of the same label, per stage code 0..3 (label -1 is dropped)."""
    out = {s: [] for s in range(4)}
    for st in stages:
        b = np.nonzero(np.diff(st) != 0)[0] + 1
        starts = np.concatenate([[0], b]); ends = np.concatenate([b, [len(st)]])
        for a, e in zip(starts, ends):
            if st[a] >= 0:
                out[int(st[a])].append(int(e - a))
    return out


def ep_stats(v):
    v = np.asarray(v)
    q = np.percentile(v, [25, 75])
    return dict(n_episodes=int(len(v)), median=float(np.median(v)), q25=float(q[0]), q75=float(q[1]), mean=float(v.mean()),
                max=int(v.max()))


# ------------------------------------------------------------------ per-stage tolerance
def _surr(k):
    sigs, stages, p = NC._G["sigs"], NC._G["stages"], NC._G["base"]
    rng = np.random.default_rng([777, k])
    S = len(sigs)
    i, j = NC.PAIRS[int(rng.integers(len(NC.PAIRS)))]
    a = int(rng.integers(S)); b = int(rng.integers(S - 1)); b += b >= a
    n = min(len(sigs[a]), len(sigs[b]))
    tau = C.time_delay_interaction(sigs[a][:n, i], sigs[b][:n, j], p)[0]
    asg = NC.assign_windows(stages[a][:n], p, len(tau))
    res = np.full((max(1, p.max_lag // 4) + 1, 4), np.nan)
    for t in range(res.shape[0]):
        lab = C.stable_label(tau, replace(p, tolerance=t))
        res[t] = NC.stage_tds(lab, asg)
    return res


def stage_tolerances(sigs, stages, base, procs=C.WORKERS):
    NC._G.update(sigs=sigs, stages=stages, base=base)
    with Pool(procs) as pool:
        R = np.stack(pool.map(_surr, range(N_SURR)))             # (draws, tol, stage)
    tols, info = [], []
    for s in range(4):
        tol, fail, rows = 0, False, []
        for t in range(R.shape[1]):
            v = R[:, t, s]; v = v[np.isfinite(v)]
            ub = v.mean() + 2 * v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else np.nan
            rows.append(dict(tol=t, n=int(len(v)), mean=float(v.mean()) if len(v) else None, ub=float(ub)))
        for r_ in rows:
            if r_["ub"] <= 5.0:
                tol = r_["tol"]
            else:
                fail = r_["tol"] == 0
                break
        tols.append(tol); info.append(dict(stage=NC.STAGES[s], tolerance=tol, calibration_failed=bool(fail), scan=rows))
    return tols, info


# ------------------------------------------------------------------ effects (as ncom_cohort.effects, settings list differs)
def effects(mean_tds, nsig):
    out = {}
    for metric, arr in (("mean_link_tds", mean_tds), ("n_significant_links", nsig)):
        out[metric] = {}
        for si, name in enumerate(SETTINGS):
            out[metric][name] = {}
            for a, b in NC.CONTRASTS:
                d = arr[si][:, NC.STAGES.index(a)] - arr[si][:, NC.STAGES.index(b)]
                out[metric][name][f"{a}-{b}"] = NC.boot_effect(d[np.isfinite(d)])
        out[metric]["gain"] = {}
        for name in SETTINGS[1:]:
            out[metric]["gain"][name] = {}
            for a, b in NC.CONTRASTS:
                ia, ib = NC.STAGES.index(a), NC.STAGES.index(b)
                d1 = arr[SETTINGS.index(name)][:, ia] - arr[SETTINGS.index(name)][:, ib]
                d0 = arr[0][:, ia] - arr[0][:, ib]
                ok = np.isfinite(d1) & np.isfinite(d0)
                out[metric]["gain"][name][f"{a}-{b}"] = NC.boot_gain(d1[ok], d0[ok])
    return out


def main(path, out_dir, n_subjects=None, n_null=NC.N_NULL):
    t0 = time.time()
    C.RESULTS = out_dir
    ids, sigs, stages, nan_counts = NC.load_cohort(path, n_subjects)
    S = len(ids)
    # 1. episode lengths
    ep = [episodes([st[:]]) for st in stages]
    eps = {s: [l for e in ep for l in e[s]] for s in range(4)}
    allep = [l for s in range(4) for l in eps[s]]
    ep_summary = {NC.STAGES[s]: ep_stats(eps[s]) for s in range(4)}
    ep_summary["all"] = ep_stats(allep)
    print("episodes", {k: (v["n_episodes"], v["median"], v["q25"], v["q75"]) for k, v in ep_summary.items()}, flush=True)
    # 2. parameters
    base_params, pinfo = NC.cohort_params(sigs)
    L = int(np.clip(round(0.5 * float(np.median(allep))), 61, 300))
    base = C.TDSParams(window=L, overlap=L // 2, max_lag=(L - 1) // 2, tolerance=1)
    tols, tinfo = stage_tolerances(sigs, stages, base)
    params = {"published": base_params["published"], "atds_cap300": base_params["atds_cap300"],
              "atds_stage": StageParams(base, tols)}
    pinfo = {k: pinfo[k] for k in ("published", "atds_cap300")}
    for k in pinfo:
        pinfo[k].pop("B_summary", None)
    pinfo["atds_stage"] = dict(window=L, step=L // 2, max_lag=base.max_lag, tolerance_by_stage=dict(zip(NC.STAGES, tols)),
                               calibration_failed_by_stage={t["stage"]: t["calibration_failed"] for t in tinfo},
                               median_episode_length=float(np.median(allep)))
    print("params", {k: (p.window, p.overlap, p.max_lag, getattr(p, "tols", p.tolerance)) for k, p in params.items()}, flush=True)
    # 3. scoring (ncom_cohort workers; settings and per-stage scoring patched)
    NC.SETTINGS = SETTINGS
    NC.pair_stage_tds = pair_stage_tds
    NC._G.update(sigs=sigs, stages=stages, params=params, n_null=n_null)
    with Pool(C.WORKERS) as pool:
        real = pool.map(NC._subject_worker, range(S))
        print(f"real done {time.time() - t0:.0f}s", flush=True)
        nulls = pool.map(NC._null_worker, range(len(NC.PAIRS)))
    print(f"null done {time.time() - t0:.0f}s", flush=True)
    R = np.stack([r[0] for r in real]); NA = np.stack([r[1] for r in real])
    TH = np.stack([n[0] for n in nulls], axis=2)
    with np.errstate(invalid="ignore"):
        sig = R > TH.transpose(0, 2, 1)[None]
    sig &= np.isfinite(R)
    mean_tds = np.full((3, S, 4), np.nan); nsig = np.zeros((3, S, 4))
    for si in range(3):
        for q in range(4):
            for s in range(S):
                v = R[s, si, :, q]
                if np.isfinite(v).any():
                    mean_tds[si, s, q] = np.nanmean(v); nsig[si, s, q] = sig[s, si, :, q].sum()
                else:
                    nsig[si, s, q] = np.nan
    eff = effects(mean_tds, nsig)
    stage_means = []
    for si, nm in enumerate(SETTINGS):
        for q, st in enumerate(NC.STAGES):
            ok = np.isfinite(mean_tds[si, :, q])
            stage_means.append(dict(setting=nm, stage=st, n_subjects=int(ok.sum()), mean_link_tds=float(mean_tds[si, ok, q].mean()),
                                    sd_across_subjects=float(mean_tds[si, ok, q].std(ddof=1)),
                                    mean_n_significant_links=float(np.nanmean(nsig[si, :, q])),
                                    mean_assigned_windows=float(NA[:, si, q].mean())))
    rows = [dict(subject=ids[s], stage=NC.STAGES[q], setting=SETTINGS[si], mean_link_tds=mean_tds[si, s, q],
                 n_significant_links=nsig[si, s, q], n_assigned_windows=int(NA[s, si, q]))
            for s in range(S) for si in range(3) for q in range(4)]
    out = dict(n_subjects=S, n_null_draws=n_null, n_boot=NC.N_BOOT, episodes=ep_summary, params=pinfo, tolerance_scan=tinfo,
               stage_means=stage_means, effects=eff, subject_stage=rows, wall_seconds=round(time.time() - t0))
    C.write_json("e4_cohort_stage.json", out)
    C.write_csv("e4_episodes.csv", [dict(stage=k, **v) for k, v in ep_summary.items()])
    C.write_csv("e4_params.csv", [dict(setting=k, **C.flat(v)) for k, v in pinfo.items()])
    C.write_csv("e4_stage_means.csv", stage_means)
    C.write_csv("e4_subject_stage.csv", rows)
    erows = []
    for metric, d in eff.items():
        for name, cd in d.items():
            if name == "gain":
                for n2, c2 in cd.items():
                    for con, e in c2.items():
                        erows.append(dict(metric=metric, setting=f"gain_{n2}_minus_published", contrast=con, n=e["n"], value=e["gain"],
                                          ci_lo=e["ci95"][0], ci_hi=e["ci95"][1], frac_positive=""))
            else:
                for con, e in cd.items():
                    erows.append(dict(metric=metric, setting=name, contrast=con, n=e["n"], value=e["dz"], ci_lo=e["ci95"][0],
                                      ci_hi=e["ci95"][1], frac_positive=e["frac_positive"]))
    C.write_csv("e4_effects.csv", erows)
    for r in stage_means:
        print({k: (round(v, 2) if isinstance(v, float) else v) for k, v in r.items()})
    for r in erows:
        print(r["metric"], r["setting"], r["contrast"], r["n"], round(r["value"], 3), round(r["ci_lo"], 3), round(r["ci_hi"], 3), r["frac_positive"])
    print(f"total {time.time() - t0:.0f}s")


if __name__ == "__main__":
    args = sys.argv[1:]
    opt = {}
    for flag in ("--n-subjects", "--n-null"):
        if flag in args:
            k = args.index(flag); opt[flag] = int(args[k + 1]); del args[k:k + 2]
    warnings.simplefilter("ignore")
    main(args[0], args[1] if len(args) > 1 else os.path.join(C.HERE, "results"), opt.get("--n-subjects"), opt.get("--n-null", NC.N_NULL))
