"""E3: harder synthetic benchmarks (regime switching, nonlinear coupling, drifting delay), 20000 samples.

Usage: ATDS_CODE=<checkout> python e3_synthetic.py [--smoke]
"""
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

NS = 20000
_G = {}


# ------------------------------------------------------------------ generators
def regime_pair(seed, c, n=NS, blk=2000, d=20):
    """Blocks alternate regime A (tc 8, coupled with strength c, delay d) and B (tc 32, uncoupled)."""
    rng = np.random.default_rng(seed)
    A = ((np.arange(n) // blk) % 2 == 0)
    x = np.where(A, C.ar1(n, 8, rng), C.ar1(n, 32, rng))
    noise = np.where(A, C.ar1(n, 8, rng), C.ar1(n, 32, rng))
    m = c * A
    y = m * np.roll(x, d) + np.sqrt(1 - m ** 2) * noise
    return x, y, A


def nl_pair(seed, kind, c, n=NS, d=20):
    rng = np.random.default_rng(seed)
    x = C.ar1(n, 8, rng); xd = np.roll(x, d)
    f = np.tanh(2 * xd) if kind == "monotone" else (xd ** 2 - 1) / np.sqrt(2)
    if kind == "monotone":
        f = f / f.std()
    return x, c * f + np.sqrt(1 - c * c) * C.ar1(n, 8, rng)


def drift_pair(seed, c=0.6, n=NS):
    rng = np.random.default_rng(seed)
    x = C.ar1(n, 8, rng)
    dly = 10 + 20 * np.arange(n) / n
    y = c * x[(np.arange(n) - np.round(dly).astype(int)) % n] + np.sqrt(1 - c * c) * C.ar1(n, 8, rng)
    return x, y


# ------------------------------------------------------------------ regime scoring (80 % purity rule)
def regime_score(x, y, A, p, d=20):
    tau = C.time_delay_interaction(x, y, p)[0]
    lab = C.stable_label(tau, p).astype(bool)
    starts = np.arange(len(tau)) * p.overlap
    cs = np.concatenate([[0], np.cumsum(A)])
    fa = (cs[starts + p.window] - cs[starts]) / p.window
    inA, inB = fa >= 0.8, fa <= 0.2
    out = {}
    for nm, m in (("A", inA), ("B", inB)):
        out[nm] = 100.0 * lab[m].mean() if m.sum() >= 3 else np.nan
    sa = lab & inA
    out["found"] = int(sa.sum() >= 3 and abs(np.median(tau[sa]) + d) <= max(1, 0.1 * d))
    return out


def _setup_regime(_):
    train = [regime_pair(1000000 + j, 0.6) for j in range(5)]
    segA = [(x[A], y[A]) for x, y, A in train]; segB = [(x[~A], y[~A]) for x, y, A in train]
    tr = [(x, y) for x, y, _ in train]
    return {"published": C.FIXED, "atds": C.E4.atds2(tr)[0], "atds_cap1000": C.E4.atds2(tr, cap=1000)[0],
            "regime_A": C.E4.atds2(segA)[0], "regime_B": C.E4.atds2(segB)[0]}


def _regime_task(k):
    P = _G["regime"]
    xc, yc, A = regime_pair(2000000 + k, 0.6)
    xn, yn, _ = regime_pair(3000000 + k, 0.0)
    out = {}
    for nm in ("published", "atds", "atds_cap1000"):
        out[nm] = (regime_score(xc, yc, A, P[nm]), regime_score(xn, yn, A, P[nm]))
    ca, na = regime_score(xc, yc, A, P["regime_A"]), regime_score(xn, yn, A, P["regime_A"])
    cb, nb = regime_score(xc, yc, A, P["regime_B"]), regime_score(xn, yn, A, P["regime_B"])
    out["atds_per_regime"] = (dict(A=ca["A"], B=cb["B"], found=ca["found"]), dict(A=na["A"], B=nb["B"], found=na["found"]))
    return out


# ------------------------------------------------------------------ nonlinear
def _setup_nl(kind):
    train = [nl_pair(4000000 + 100 * (kind == "non_monotone") + j, kind, 0.6) for j in range(5)]
    return {"published": C.FIXED, "atds": C.E4.atds2(train)[0]}


def _nl_task(args):
    kind, k = args
    o = 100000 * (kind == "non_monotone")
    xp, yp = nl_pair(5000000 + o + k, kind, 0.6); xn, yn = nl_pair(6000000 + o + k, kind, 0.0)
    out = {}
    for nm, p in _G["nl"][kind].items():
        s, tau, lab = C.score(xp, yp, p)
        out[nm] = (s, C.score(xn, yn, p)[0], int(C.found(tau, lab, 20)))
    sp, dl, _ = C.xcorr_full(xp, yp)
    out["xcorr_full"] = (sp, C.xcorr_full(xn, yn)[0], int(C.hit_global(dl, 20)))
    return kind, out


# ------------------------------------------------------------------ drift
def _setup_drift(_):
    tr = [drift_pair(7000000 + j) for j in range(5)]
    return {"published": C.FIXED, "atds": C.E4.atds2(tr)[0], "atds_cap1000": C.E4.atds2(tr, cap=1000)[0]}


def _drift_task(k):
    x, y = drift_pair(8000000 + k)
    out = {}
    for nm, p in _G["drift"].items():
        tau, t, _ = C.time_delay_interaction(x, y, p); lab = C.stable_label(tau, p).astype(bool)
        true = 10 + 20 * t / NS
        out[nm] = (float(C.tds_score(lab)), np.abs(tau[lab] + true[lab]))
    return out


def main():
    n, smoke = C.smoke_n(40)
    C.set_results_for_smoke(smoke)
    nn = 4 if smoke else 100
    t0 = time.time(); res = {}; rows = []
    _G["regime"] = _setup_regime(0)
    _G["nl"] = {k: _setup_nl(k) for k in ("monotone", "non_monotone")}
    _G["drift"] = _setup_drift(0)
    print(f"setup {time.time() - t0:.0f}s", flush=True)
    with Pool(C.WORKERS) as pool:
        rg = pool.map(_regime_task, range(n))
        print(f"regime {time.time() - t0:.0f}s", flush=True)
        nl = pool.map(_nl_task, [(kd, k) for kd in ("monotone", "non_monotone") for k in range(nn)], chunksize=4)
        print(f"nonlinear {time.time() - t0:.0f}s", flush=True)
        dr = pool.map(_drift_task, range(n))

    # 1. regime switching
    E = {}
    for nm in ("published", "atds", "atds_cap1000", "atds_per_regime"):
        cA = [r[nm][0]["A"] for r in rg]; nA = [r[nm][1]["A"] for r in rg]
        cB = [r[nm][0]["B"] for r in rg]; nB = [r[nm][1]["B"] for r in rg]
        ok = [i for i in range(n) if np.isfinite(cA[i]) and np.isfinite(nA[i])]
        pp = _G["regime"].get(nm)
        pr = ({"A": C.pdict(_G["regime"]["regime_A"]), "B": C.pdict(_G["regime"]["regime_B"])} if nm == "atds_per_regime"
              else C.pdict(pp))
        E[nm] = dict(params=pr, **C.auc_row([cA[i] for i in ok], [nA[i] for i in ok]),
                     regimeA_coupled_mean=round(float(np.nanmean(cA)), 2), regimeA_control_mean=round(float(np.nanmean(nA)), 2),
                     false_stability_B_coupled_mean=round(float(np.nanmean(cB)), 2),
                     false_stability_B_control_mean=round(float(np.nanmean(nB)), 2),
                     **C.hit_row([r[nm][0]["found"] for r in rg]))
        rows.append(dict(benchmark="regime_switching", method=nm, **C.flat({k: v for k, v in E[nm].items()})))
    res["regime_switching"] = E

    # 2. nonlinear
    for kind in ("monotone", "non_monotone"):
        E = {}
        sel = [o for kd, o in nl if kd == kind]
        for nm in sel[0]:
            pos, neg, hit = zip(*[o[nm] for o in sel])
            E[nm] = dict(params=C.pdict(_G["nl"][kind][nm]) if nm in _G["nl"][kind] else None,
                         **C.auc_row(pos, neg), **C.hit_row(hit), chance=round(float(np.mean(neg)), 3))
            rows.append(dict(benchmark=f"nonlinear_{kind}", method=nm, **C.flat(E[nm])))
        res[f"nonlinear_{kind}"] = E

    # 3. drifting delay
    E = {}
    for nm in _G["drift"]:
        sc = [o[nm][0] for o in dr]; err = [o[nm][1] for o in dr]; allerr = np.concatenate(err)
        E[nm] = dict(params=C.pdict(_G["drift"][nm]), score_mean=round(float(np.mean(sc)), 2), score_sd=round(float(np.std(sc, ddof=1)), 2),
                     n_stable_windows=int(len(allerr)),
                     median_abs_err_pooled=round(float(np.median(allerr)), 2) if len(allerr) else None,
                     iqr_abs_err=[round(float(v), 2) for v in np.percentile(allerr, [25, 75])] if len(allerr) else None,
                     median_of_pair_medians=round(float(np.median([np.median(e) for e in err if len(e)])), 2))
        rows.append(dict(benchmark="drifting_delay", method=nm, **C.flat(E[nm])))
    res["drifting_delay"] = E
    res["wall_seconds"] = round(time.time() - t0); res["n_pairs"] = dict(regime=n, drift=n, nonlinear=nn)
    C.write_json("e3_synthetic.json", res); C.write_csv("e3_synthetic.csv", rows)
    for b, v in res.items():
        if isinstance(v, dict) and b != "n_pairs":
            print(b)
            for m, e in v.items():
                print("  ", m, {k: x for k, x in e.items()})
    print(f"total {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
