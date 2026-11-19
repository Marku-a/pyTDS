"""Simulated delay patterns: does aTDS recover the true delay?  See SPEC.md.

    python experiments/sim_patterns/patterns.py

Uses this repo's pyTDS (through atds21). "fixed TDS" = atds21.PUBLISHED (60/30/+-30/+-1).
"""
import csv
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "atds_real"))
import atds21 as A  # noqa: E402

N, TC, C = 20000, 8.0, 0.6
PATTERNS = ["fixed", "step", "drift", "harmonic", "onoff"]
BLOCK = 2000
OUT = os.path.join(HERE, "results")


def delay(kind, n=N):
    t = np.arange(n)
    if kind in ("fixed", "onoff"):
        d = np.full(n, 20.0)
    elif kind == "step":
        d = np.where(t < n / 2, 10.0, 30.0)
    elif kind == "drift":
        d = 10 + 20 * t / (n - 1)
    elif kind == "harmonic":
        d = 20 + 10 * np.sin(2 * np.pi * t / 10000)
    else:
        raise ValueError(kind)
    return np.rint(d).astype(int)


def coupling(kind, c, n=N):
    """Per-sample coupling c(t)."""
    if kind == "onoff":
        return np.where((np.arange(n) // BLOCK) % 2 == 0, c, 0.0)
    return np.full(n, float(c))


def ar1(rng, n=N, tc=TC):
    phi = np.exp(-1.0 / tc)
    e = rng.standard_normal(n + 500)
    x = np.zeros(n + 500)
    for i in range(1, n + 500):
        x[i] = phi * x[i - 1] + e[i]
    x = x[500:]
    return (x - x.mean()) / x.std()


def make_pair(kind, c, seed, n=N):
    rng = np.random.default_rng(seed)
    x = ar1(rng, n); eta = ar1(rng, n)
    d = delay(kind, n); cc = coupling(kind, c, n)
    idx = np.clip(np.arange(n) - d, 0, None)
    y = cc * x[idx] + np.sqrt(1 - cc ** 2) * eta
    return x, y


def seeds(k, role, j):
    return 1000 * k + {"train": 0, "coupled": 100, "control": 500}[role] + j


def window_truth(kind, p, n_win):
    """Centre-time, true tau, and fractional coupled share of each window (pyTDS: centre = start + L//2)."""
    L, step = p.window, p.overlap
    starts = np.arange(n_win) * step
    t = starts + L // 2
    d = delay(kind)
    tau_true = -d[np.minimum(t, N - 1)].astype(float)
    if kind == "onoff":
        cs = np.concatenate([[0], np.cumsum(coupling(kind, 1.0) > 0)])
        frac = (cs[starts + L] - cs[starts]) / L
    else:
        frac = np.ones(n_win)
    return t, tau_true, frac


def analyse(args):
    kind, role, j, p = args
    k = PATTERNS.index(kind) + 1
    x, y = make_pair(kind, C if role == "coupled" else 0.0, seeds(k, role, j))
    tau = A.time_delay_interaction(x, y, p)[0]
    lab = np.asarray(A.stable_label(tau, p)).astype(bool)
    return tau, lab, float(A.tds_score(lab))


def auc(pos, neg):
    pos, neg = np.asarray(pos), np.asarray(neg)
    return float(((pos[:, None] > neg[None, :]).sum() + 0.5 * (pos[:, None] == neg[None, :]).sum()) / (len(pos) * len(neg)))


def evaluate(kind, p, pool):
    jobs = [(kind, "coupled", j, p) for j in range(40)] + [(kind, "control", j, p) for j in range(40)]
    res = pool.map(analyse, jobs, chunksize=4)
    cpl, ctl = res[:40], res[40:]
    n_win = len(cpl[0][0])
    _, tau_true, frac = window_truth(kind, p, n_win)
    inb = frac >= 0.8            # windows inside coupled blocks (all windows unless onoff)
    outb = frac <= 0.2           # windows inside uncoupled blocks (onoff only)
    stable_pct, err_all, false_blocks = [], [], []
    n_st = 0
    for tau, lab, _ in cpl:
        stable_pct.append(100.0 * lab[inb].mean())
        m = lab & inb
        n_st += int(m.sum())
        err_all.append(np.abs(tau[m] - tau_true[m]))
        if kind == "onoff":
            false_blocks.append(100.0 * lab[outb].mean())
    err = np.concatenate(err_all)
    out = dict(
        stable_pct=float(np.mean(stable_pct)),
        track_pct=float(100.0 * np.mean(err <= 3)) if len(err) else float("nan"),
        mae=float(np.median(err)) if len(err) else float("nan"),
        false_stable_pct=float(np.mean([100.0 * r[1].mean() for r in ctl])),
        auc=auc([r[2] for r in cpl], [r[2] for r in ctl]),
        n_stable_windows=n_st,
    )
    if kind == "onoff":
        out["false_stable_uncoupled_blocks_pct"] = float(np.mean(false_blocks))
    return out


def main():
    t0 = time.time()
    results, rows = {}, []
    with Pool(3) as pool:
        for k, kind in enumerate(PATTERNS, 1):
            train = [make_pair(kind, C, seeds(k, "train", j)) for j in range(5)]
            sigs = [s for pr in train for s in pr]
            pa = A.window_params(sigs)
            pa, fail = A.calib(pa, train)
            setting = {"fixed": A.PUBLISHED, "atds": pa}
            results[kind] = {"atds_params": dict(window=pa.window, overlap=pa.overlap, max_lag=pa.max_lag,
                                                 tolerance=pa.tolerance, calib_fail=bool(fail))}
            for name, p in setting.items():
                m = evaluate(kind, p, pool)
                results[kind][name] = m
                rows.append(dict(pattern=kind, setting="fixed TDS" if name == "fixed" else "aTDS",
                                 window=p.window, max_lag=p.max_lag, tol=p.tolerance, **m))
                print(kind, name, f"L={p.window} lag={p.max_lag} tol={p.tolerance}",
                      {a: round(b, 2) for a, b in m.items()}, flush=True)
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "patterns.json"), "w") as f:
        json.dump(results, f, indent=2)
    keys = ["pattern", "setting", "window", "max_lag", "tol", "stable_pct", "track_pct", "mae", "false_stable_pct",
            "false_stable_uncoupled_blocks_pct", "auc", "n_stable_windows"]
    with open(os.path.join(OUT, "patterns.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, keys, restval=""); w.writeheader(); w.writerows(rows)
    print(f"wall time {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
