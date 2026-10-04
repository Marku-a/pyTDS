"""
Synthetic ground-truth test: do system-adaptive TDS params beat fixed defaults?

Usage: python experiments/adaptive/exp_synthetic.py --exp E1|E1c|E1d|E1b|E2|all [--quick]

System-level rule: for each condition the params are computed ONCE from the
pooled signals (x and y) of the first 3 coupled and first 3 uncoupled trials
and reused for all trials. The "oracle" picks, per condition, the grid param
with best AUC -- it peeks at the labels, so it is an OPTIMISTIC UPPER BOUND.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from itertools import combinations
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "examples"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from adaptive import adaptive_params, decorrelation_time  # noqa: E402
from common import auc, delay_metrics, oracle_grid, variant_params  # noqa: E402
from generators import ar1, coupled_pair  # noqa: E402
from pyTDS.core import stable_label, tds_score, time_delay_interaction  # noqa: E402
from pyTDS.params import TDSParams  # noqa: E402

RES_DIR = os.path.join(HERE, "results")
FIG_DIR = os.path.join(ROOT, "adaptive-TDS reports", "figures")
VARIANTS = ["fixed_default", "adaptive_fixedtol", "adaptive_calibrated"]
EXP_ID = {"E1": 1, "E1c": 2, "E1d": 3, "E1b": 4, "E2": 5}
GRID = oracle_grid()
COLORS = {"fixed_default": "tab:gray", "adaptive_fixedtol": "tab:blue",
          "adaptive_calibrated": "tab:green", "oracle": "tab:red"}


# ---------------------------------------------------------------- generation
def gen_trial(spec: dict, seed_key: list):
    """Return (x, y, mask) for a trial spec; deterministic given seed_key."""
    rng = np.random.default_rng(seed_key)
    if spec["kind"] == "ar":
        return coupled_pair(spec["n"], spec["tc_x"], spec["tc_y"], spec["delay"],
                            spec["c"], rng, spec.get("block_len"))
    if spec["kind"] == "rulkov":
        from rulkov_validation import coupled_rulkov_with_delay
        x1, x2 = coupled_rulkov_with_delay(
            spec["n"], 0.001, [4.5, 4.1], [0.01, -0.01], 1.0, 1.0, spec["g"], spec["delay"])
        # Generator is deterministic -> trial variety via per-trial observation noise
        x1 = x1 + rng.normal(0, 0.01 * x1.std(), len(x1))
        x2 = x2 + rng.normal(0, 0.01 * x2.std(), len(x2))
        return x1, x2, np.ones(len(x1))
    raise ValueError(spec["kind"])


def seed_key(eid, ci, coupled, t):
    return [EXP_ID[eid], ci, int(coupled), t]


def sys_signals(spec_c, spec_u, eid, ci):
    sig = []
    for flag, spec in ((1, spec_c), (0, spec_u)):
        for t in range(3):
            x, y, _ = gen_trial(spec, seed_key(eid, ci, flag, t))
            sig += [x, y]
    return sig


# ------------------------------------------------------------------ workers
def _taus(x, y, plist, cache):
    out = []
    for p in plist:
        k = (p.window, p.overlap, p.max_lag, p.window_anchor)
        if k not in cache:
            cache[k] = time_delay_interaction(x, y, p)
        out.append(cache[k])
    return out


def worker_score(args):
    """Score + median stable tau for each (key, params); tau cached per (L, step, max_lag)."""
    spec, sk, named = args
    x, y, _ = gen_trial(spec, sk)
    plist = [p for _, p in named] + GRID
    cache = {}
    res = {}
    for i, ((tau, _, _), p) in enumerate(zip(_taus(x, y, plist, cache), plist)):
        lbl = stable_label(tau, p)
        st = tau[lbl == 1]
        key = named[i][0] if i < len(named) else f"grid{i - len(named)}"
        res[key] = (tds_score(lbl), float(np.median(st)) if len(st) else float("nan"))
    return res


def worker_episode(args):
    spec, sk, named = args
    x, y, m = gen_trial(spec, sk)
    out = {}
    for name, p in named:
        tau, t_vec, _ = time_delay_interaction(x, y, p)
        lbl = stable_label(tau, p)
        truth = m[t_vec.astype(int)] > 0
        tpr = float(np.mean(lbl[truth] == 1)) if truth.any() else float("nan")
        tnr = float(np.mean(lbl[~truth] == 0)) if (~truth).any() else float("nan")
        out[name] = (tds_score(lbl), 0.5 * (tpr + tnr))
    return out


def gen_network(seed):
    rng = np.random.default_rng([EXP_ID["E1b"], seed])
    z = lambda a: (a - a.mean()) / a.std()
    n = 30000
    tcs = [2] * 4 + [40] * 4
    own = [ar1(n, tc, rng) for tc in tcs]
    s = [None] * 8
    for i in (0, 2, 4, 6):
        s[i] = own[i]
    s[1] = z(0.6 * np.roll(s[0], 5) + 0.8 * own[1])
    s[3] = z(0.6 * np.roll(s[2], 3) + 0.8 * own[3])
    s[7] = z(0.6 * np.roll(s[6], 60) + 0.8 * own[7])
    s[5] = z(0.5 * np.roll(s[4], 30) + 0.5 * np.roll(s[1], 10) + 0.7 * own[5])
    return [z(a) for a in s]


# (i, j) pair with i<j -> true tau (s1=si, s2=sj; sj = roll(si, d) -> tau=-d)
NET_EDGES = {(0, 1): -5, (2, 3): -3, (4, 5): -30, (6, 7): -60, (1, 5): -10}


def worker_network(seed):
    sig = gen_network(seed)
    pairs = list(combinations(range(8), 2))
    kw = dict(tolerance="calibrated", n_shifts=5, max_pairs=20)
    p_q1, info_q1 = adaptive_params(sig, quantile=1.0, **kw)
    p_q05, _ = adaptive_params(sig, quantile=0.5, **kw)
    tc = np.array(info_q1["tau_c"])
    gm = float(np.sqrt(tc.min() * tc.max()))
    groups = [[i for i in range(8) if tc[i] <= gm], [i for i in range(8) if tc[i] > gm]]
    gparams = {}
    for g, idx in enumerate(groups):
        gparams[g] = adaptive_params([sig[i] for i in idx], **kw)[0] if len(idx) >= 2 else p_q1
    grp_of = {i: g for g, idx in enumerate(groups) for i in idx}
    from pyTDS.core import tds
    variants = ["fixed_default", "system_q1", "system_q05", "perpair", "grouped"]
    res = {v: {} for v in variants}
    used = {"fixed_default": TDSParams(), "system_q1": p_q1, "system_q05": p_q05,
            "grouped": {f"g{g}": gparams[g] for g in gparams}}
    for (i, j) in pairs:
        pp = {"fixed_default": TDSParams(), "system_q1": p_q1, "system_q05": p_q05,
              "perpair": adaptive_params([sig[i], sig[j]], tolerance="calibrated", n_shifts=5)[0],
              "grouped": gparams[grp_of[i]] if grp_of[i] == grp_of[j] else p_q1}
        for v in variants:
            r = tds(sig[i], sig[j], pp[v])
            st = r["stable_taus"]
            res[v][f"{i}-{j}"] = (r["score"], float(np.median(st)) if len(st) else float("nan"))
            if v == "perpair":
                res.setdefault("_perpair_params", {})[f"{i}-{j}"] = [pp[v].window, pp[v].tolerance]
    pdesc = {}
    for v, p in used.items():
        if isinstance(p, dict):
            pdesc[v] = {k: (q.window, q.overlap, q.max_lag, q.tolerance) for k, q in p.items()}
        else:
            pdesc[v] = (p.window, p.overlap, p.max_lag, p.tolerance)
    return {"seed": seed, "res": res, "params": pdesc, "tau_c": tc.tolist(), "groups": groups}


# ---------------------------------------------------------------- summaries
def pdict(p):
    return {"window": p.window, "step": p.overlap, "max_lag": p.max_lag, "tolerance": p.tolerance}


def nan_ms(a):
    a = np.asarray(a, dtype=float)
    return float(np.nanmean(a)), float(np.nanstd(a))


def run_detection(pool, eid, ci, cname, spec_c, spec_u, ntr, true_tau):
    """Generic coupled-vs-uncoupled condition with 3 variants + oracle."""
    system = sys_signals(spec_c, spec_u, eid, ci)
    named, infos = [], {}
    for v in VARIANTS:
        p, info = variant_params(v, system)
        named.append((v, p))
        infos[v] = {k: info.get(k) for k in ("window", "step", "max_lag", "tolerance", "tau_sys")}
    jobs = [(spec_c, seed_key(eid, ci, 1, t), named) for t in range(ntr)] + \
           [(spec_u, seed_key(eid, ci, 0, t), named) for t in range(ntr)]
    out = pool.map(worker_score, jobs, chunksize=1)
    rc, ru = out[:ntr], out[ntr:]

    def stats(key):
        sc = [r[key][0] for r in rc]
        su = [r[key][0] for r in ru]
        hit, rel = delay_metrics([r[key][1] for r in rc], true_tau)
        return {"auc": auc(sc, su), "coupled_mean": float(np.mean(sc)), "coupled_sd": float(np.std(sc)),
                "uncoupled_mean": float(np.mean(su)), "uncoupled_sd": float(np.std(su)),
                "hit_rate": hit, "median_rel_err": rel,
                "coupled_scores": sc, "uncoupled_scores": su,
                "coupled_median_tau": [r[key][1] for r in rc]}

    variants = {}
    for v in VARIANTS:
        variants[v] = stats(v)
        variants[v]["params"] = infos[v]
    best, bkey = None, None
    for gi, gp in enumerate(GRID):
        s = stats(f"grid{gi}")
        k = (s["auc"], s["coupled_mean"])
        if best is None or k > best:
            best, bkey, bs = k, gi, s
    bs["params"] = pdict(GRID[bkey])
    bs["note"] = "OPTIMISTIC UPPER BOUND: grid param chosen using the labels"
    variants["oracle"] = bs
    return {"condition": cname, "true_tau": true_tau, "variants": variants}


def rows_from(eid, conds):
    rows = []
    for c in conds:
        for v, s in c["variants"].items():
            p = s["params"]
            rows.append({"exp": eid, "condition": c["condition"], "variant": v,
                         "auc": s["auc"], "coupled_mean": s["coupled_mean"],
                         "coupled_sd": s["coupled_sd"], "uncoupled_mean": s["uncoupled_mean"],
                         "uncoupled_sd": s["uncoupled_sd"], "hit_rate": s["hit_rate"],
                         "median_rel_err": s["median_rel_err"], "window": p["window"],
                         "step": p["step"], "max_lag": p["max_lag"], "tolerance": p["tolerance"]})
    return rows


def save(eid, payload, rows, quick):
    os.makedirs(RES_DIR, exist_ok=True)
    suf = "_quick" if quick else ""
    with open(os.path.join(RES_DIR, f"synthetic_{eid}{suf}.json"), "w") as f:
        json.dump(payload, f, indent=1, default=float)
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(os.path.join(RES_DIR, f"synthetic_{eid}{suf}_summary.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def savefig(fig, name, quick):
    os.makedirs(FIG_DIR, exist_ok=True)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, name + ("_quick" if quick else "") + ".png"), dpi=130)
    plt.close(fig)


def get(conds, name, v, key):
    for c in conds:
        if c["condition"] == name:
            return c["variants"][v][key]


def bar_groups(ax, labels, series, ylabel, title, errs=None):
    w = 0.8 / len(series)
    x = np.arange(len(labels))
    for k, (name, vals) in enumerate(series.items()):
        ax.bar(x + k * w - 0.4 + w / 2, vals, w, label=name, color=COLORS.get(name),
               yerr=None if errs is None else errs[name], capsize=2)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(fontsize=8)


ALLV = VARIANTS + ["oracle"]


# -------------------------------------------------------------- experiments
def run_E1(pool, quick):
    ntr = 2 if quick else 20
    rs, cs = [1, 2, 4, 8, 16], [0.3, 0.6]
    conds, ci = [], 0
    for c in cs:
        for r in rs:
            sc = dict(kind="ar", n=30000, tc_x=2 * r, tc_y=2 * r, delay=5 * r, c=c)
            su = dict(sc, c=0.0)
            conds.append(run_detection(pool, "E1", ci, f"r={r},c={c}", sc, su, ntr, -5 * r))
            conds[-1].update(r=r, c=c)
            ci += 1
            print("  E1", conds[-1]["condition"], flush=True)
    save("E1", {"conditions": conds, "note": "oracle = optimistic upper bound"}, rows_from("E1", conds), quick)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    for ax, c in zip(axes, cs):
        for v in ALLV:
            ax.plot(rs, [get(conds, f"r={r},c={c}", v, "auc") for r in rs], "o-", label=v, color=COLORS[v])
        ax.set_xscale("log", base=2)
        ax.set_xlabel("time-scale factor r (tc=2r, delay=5r)")
        ax.set_ylabel("AUC (coupled vs uncoupled)")
        ax.set_title(f"E1: AUC vs scale, c={c}")
        ax.legend(fontsize=8)
    savefig(fig, "E1_auc_vs_scale", quick)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    for ax, c in zip(axes, cs):
        for v in ALLV:
            ax.plot(rs, [get(conds, f"r={r},c={c}", v, "hit_rate") for r in rs], "o-", label=v, color=COLORS[v])
        ax.set_xscale("log", base=2)
        ax.set_xlabel("time-scale factor r (tc=2r, delay=5r)")
        ax.set_ylabel("delay hit rate")
        ax.set_title(f"E1: delay hit rate vs scale, c={c}")
        ax.legend(fontsize=8)
    savefig(fig, "E1_hit_vs_scale", quick)
    return conds


def run_E1c(pool, quick):
    ntr = 2 if quick else 20
    conds = []
    for ci, (tx, ty) in enumerate([(2, 2), (2, 30), (30, 2), (30, 30)]):
        sc = dict(kind="ar", n=30000, tc_x=tx, tc_y=ty, delay=10, c=0.5)
        conds.append(run_detection(pool, "E1c", ci, f"tcx={tx},tcy={ty}", sc, dict(sc, c=0.0), ntr, -10))
        print("  E1c", conds[-1]["condition"], flush=True)
    save("E1c", {"conditions": conds, "note": "oracle = optimistic upper bound"}, rows_from("E1c", conds), quick)
    labels = [c["condition"] for c in conds]
    fig, ax = plt.subplots(figsize=(9, 4))
    bar_groups(ax, labels, {v: [c["variants"][v]["auc"] for c in conds] for v in ALLV},
               "AUC", "E1c: mixed time scales, AUC (oracle = optimistic upper bound)")
    ax.set_xlabel("(tc_x, tc_y)")
    ax.set_ylim(0, 1.05)
    savefig(fig, "E1c_auc", quick)
    return conds


def run_E1d(pool, quick):
    ntr = 2 if quick else 10
    conds = []
    for ci, r in enumerate([1, 4]):
        sc = dict(kind="ar", n=30000, tc_x=2 * r, tc_y=2 * r, delay=5 * r, c=0.6, block_len=2000)
        su = dict(sc, c=0.0)
        system = sys_signals(sc, su, "E1d", ci)
        named, vs = [], {}
        for v in VARIANTS:
            p, _ = variant_params(v, system)
            named.append((v, p))
        out = pool.map(worker_episode, [(sc, seed_key("E1d", ci, 1, t), named) for t in range(ntr)], chunksize=1)
        for v, p in named:
            bacc = [o[v][1] for o in out]
            sco = [o[v][0] for o in out]
            vs[v] = {"balanced_acc_mean": float(np.mean(bacc)), "balanced_acc_sd": float(np.std(bacc)),
                     "score_mean": float(np.mean(sco)), "score_sd": float(np.std(sco)),
                     "balanced_acc": bacc, "scores": sco, "params": pdict(p)}
        conds.append({"condition": f"r={r}", "variants": vs})
        print("  E1d", f"r={r}", flush=True)
    rows = []
    for c in conds:
        for v, s in c["variants"].items():
            rows.append({"exp": "E1d", "condition": c["condition"], "variant": v,
                         "balanced_acc_mean": s["balanced_acc_mean"], "balanced_acc_sd": s["balanced_acc_sd"],
                         "score_mean": s["score_mean"], "score_sd": s["score_sd"], **s["params"]})
    save("E1d", {"conditions": conds}, rows, quick)
    fig, ax = plt.subplots(figsize=(7, 4))
    bar_groups(ax, [c["condition"] for c in conds],
               {v: [c["variants"][v]["balanced_acc_mean"] for c in conds] for v in VARIANTS},
               "balanced accuracy (window label vs coupling mask)", "E1d: episodic coupling, c=0.6",
               errs={v: [c["variants"][v]["balanced_acc_sd"] for c in conds] for v in VARIANTS})
    ax.set_xlabel("time-scale factor r")
    ax.axhline(0.5, color="k", ls=":", lw=0.8)
    savefig(fig, "E1d_balanced_accuracy", quick)
    return conds


def run_E1b(pool, quick):
    nseed = 2 if quick else 10
    outs = pool.map(worker_network, range(nseed), chunksize=1)
    variants = ["fixed_default", "system_q1", "system_q05", "perpair", "grouped"]
    summ = {}
    for v in variants:
        aucs, hits = [], []
        for o in outs:
            sc = {k: val[0] for k, val in o["res"][v].items()}
            pos = [sc[f"{i}-{j}"] for (i, j) in NET_EDGES]
            neg = [s for k, s in sc.items() if tuple(map(int, k.split("-"))) not in NET_EDGES]
            aucs.append(auc(pos, neg))
            meds = [o["res"][v][f"{i}-{j}"][1] for (i, j) in NET_EDGES]
            hits += [float(abs(m - t) <= max(1.0, 0.1 * abs(t))) if not np.isnan(m) else 0.0
                     for m, t in zip(meds, NET_EDGES.values())]
        summ[v] = {"edge_auc_mean": float(np.mean(aucs)), "edge_auc_sd": float(np.std(aucs)),
                   "edge_auc": aucs, "edge_hit_rate": float(np.mean(hits)),
                   "params_seed0": outs[0]["params"].get(v, "per-pair (see _perpair_params)")}
    pp = {k: v for k, v in outs[0]["res"].get("_perpair_params", {}).items()}
    for o in outs:
        o["res"].pop("_perpair_params", None)
    save("E1b", {"summary": summ, "per_seed": outs, "true_edges_tau": {f"{i}-{j}": t for (i, j), t in NET_EDGES.items()},
                 "perpair_params_seed0 [window,tol]": pp},
         [{"exp": "E1b", "condition": "network", "variant": v, "edge_auc_mean": s["edge_auc_mean"],
           "edge_auc_sd": s["edge_auc_sd"], "edge_hit_rate": s["edge_hit_rate"],
           "params_seed0": json.dumps(s["params_seed0"])} for v, s in summ.items()], quick)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(variants, [summ[v]["edge_auc_mean"] for v in variants],
           yerr=[summ[v]["edge_auc_sd"] for v in variants], capsize=3,
           color=["tab:gray", "tab:blue", "tab:cyan", "tab:orange", "tab:green"], label="mean +- sd over seeds")
    ax.set_ylabel("edge AUC (5 edges vs 23 non-edges)")
    ax.set_xlabel("variant")
    ax.set_title("E1b: mixed-scale network, edge detection")
    ax.legend(fontsize=8)
    savefig(fig, "E1b_edge_auc", quick)
    return summ


def run_E2(pool, quick):
    ntr = 2 if quick else 10
    conds = []
    for ci, d in enumerate([8, 16, 32]):
        sc = dict(kind="rulkov", n=6000, g=0.12, delay=d)
        # empirically (default tds, see report) stable tau = +delay for x1, x2 as returned
        conds.append(run_detection(pool, "E2", ci, f"delay={d}", sc, dict(sc, g=0.0), ntr, +d))
        conds[-1]["delay"] = d
        print("  E2", conds[-1]["condition"], flush=True)
    save("E2", {"conditions": conds, "variation_method": "generator is deterministic; per-trial "
                "observation noise N(0, 0.01*std) added to x1,x2 with per-trial rng seed; true tau=+delay "
                "(verified empirically with default tds: delay 16 -> median stable tau +17)",
                "note": "oracle = optimistic upper bound"}, rows_from("E2", conds), quick)
    ds = [8, 16, 32]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for ax, key, yl in zip(axes, ["auc", "hit_rate"], ["AUC (g=0.12 vs g=0)", "delay hit rate"]):
        for v in ALLV:
            ax.plot(ds, [get(conds, f"delay={d}", v, key) for d in ds], "o-", label=v, color=COLORS[v])
        ax.set_xlabel("true delay (samples)")
        ax.set_ylabel(yl)
        ax.set_title(f"E2 Rulkov: {yl} vs delay")
        ax.legend(fontsize=8)
    savefig(fig, "E2_auc_hit_vs_delay", quick)
    return conds


RUNNERS = {"E1": run_E1, "E1c": run_E1c, "E1d": run_E1d, "E1b": run_E1b, "E2": run_E2}


def print_summary(eid, res):
    print(f"\n=== {eid} ===")
    if eid == "E1b":
        for v, s in res.items():
            print(f"{v:14s} AUC {s['edge_auc_mean']:.3f}+-{s['edge_auc_sd']:.3f} hit {s['edge_hit_rate']:.2f} params {s['params_seed0']}")
        return
    for c in res:
        for v, s in c["variants"].items():
            p = s["params"]
            ptxt = f"L={p['window']} step={p['step']} lag={p['max_lag']} tol={p['tolerance']}"
            if eid == "E1d":
                print(f"{c['condition']:16s} {v:20s} bacc {s['balanced_acc_mean']:.3f} score {s['score_mean']:.1f} {ptxt}")
            else:
                print(f"{c['condition']:16s} {v:20s} AUC {s['auc']:.3f} c {s['coupled_mean']:.1f} u {s['uncoupled_mean']:.1f} "
                      f"hit {s['hit_rate']:.2f} {ptxt}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exp", default="all", choices=list(RUNNERS) + ["all"])
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    names = list(RUNNERS) if a.exp == "all" else [a.exp]
    t0 = time.time()
    with Pool(4) as pool:
        for e in names:
            t1 = time.time()
            res = RUNNERS[e](pool, a.quick)
            print_summary(e, res)
            print(f"[{e}] {time.time() - t1:.0f}s", flush=True)
    print(f"total {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
