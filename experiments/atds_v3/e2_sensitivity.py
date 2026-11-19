"""E2: sensitivity of the aTDS design constants, one at a time around the default.

Slow-down benchmark (tc = 2r, d = 5r, c = 0.3, N = 20000), r in {1, 4, 8}, 100 + 100 pairs, seeds as exp4.U1;
rhythm threshold on the exp4.U6 rhythm case (AR(2), period about 40, delays 8 and 30, 30 + 30 pairs).
Usage: ATDS_CODE=<checkout> python e2_sensitivity.py [--smoke]
"""
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

DEFAULT = dict(n_eff=30.0, budget=5.0, z=2.0, floor=61, n_shift=20)
GRID = [("n_eff", v) for v in (15, 30, 45, 60)] + [("budget", v) for v in (1.0, 5.0, 10.0)] + \
       [("z", v) for v in (0.0, 2.0, 3.0)] + [("floor", v) for v in (31, 61, 121)] + [("n_shift", v) for v in (10, 20, 40)]
RS = [1, 4, 8]
_G = {}


def _key(kw):
    return tuple(sorted(kw.items()))


def _cfg(name, v):
    kw = dict(DEFAULT); kw[name] = v
    return kw


def _setup(args):
    r, kw = args
    tc, d, c = 2 * r, 5 * r, 0.3
    train = [C.pair(800000 + r * 10 + j, tc, d, c) for j in range(5)]
    p, fail = C.atds_var(train, **dict(kw))
    return (r, _key(kw)), (p, fail)


def _score(args):
    r, k = args
    tc, d, c = 2 * r, 5 * r, 0.3
    xp, yp = C.pair(60000 + r * 1000 + k, tc, d, c)
    xn, yn = C.pair(70000 + r * 1000 + k, tc, d, 0.0)
    out = {}
    for key, (p, _) in _G["P"].items():
        if key[0] == r:
            s, tau, lab = C.score(xp, yp, p)
            out[key] = (s, C.score(xn, yn, p)[0], int(C.found(tau, lab, d)))
    return r, out


def _setup_rhythm(args):
    d, thr = args
    train = [C.osc_pair(700000 + d + j, d, 0.6) for j in range(5)]
    p, fail = C.atds_var(train, min_peak=thr)
    return (d, thr), (p, fail)


def _score_rhythm(args):
    d, k = args
    xp, yp = C.osc_pair(710000 + d * 100 + k, d, 0.6)
    xn, yn = C.osc_pair(720000 + d * 100 + k, d, 0.0)
    out = {}
    for key, (p, _) in _G["R"].items():
        if key[0] == d:
            s, tau, lab = C.score(xp, yp, p)
            out[key] = (s, C.score(xn, yn, p)[0], int(C.found(tau, lab, d)))
    return d, out


def _rows(label, value, case, p, fail, pos, neg, hit):
    return dict(constant=label, value=value, case=case, window=p.window, step=p.overlap, max_lag=p.max_lag,
                tolerance=p.tolerance, calib_failed=fail, **C.auc_row(pos, neg), **C.hit_row(hit),
                chance_mean=round(float(np.mean(neg)), 3), coupled_mean=round(float(np.mean(pos)), 2))


def main():
    n, smoke = C.smoke_n(100)
    C.set_results_for_smoke(smoke)
    rs = [1, 4] if smoke else RS
    t0 = time.time()
    uniq = {}
    for name, v in GRID:
        uniq[_key(_cfg(name, v))] = _cfg(name, v)
    with Pool(C.WORKERS) as pool:
        _G["P"] = dict(pool.map(_setup, [(r, kw) for r in rs for kw in uniq.values()]))
        print(f"setup {time.time() - t0:.0f}s ({len(_G['P'])} parameter sets)", flush=True)
    with Pool(C.WORKERS) as pool:
        res = pool.map(_score, [(r, k) for r in rs for k in range(n)], chunksize=2)
    rows = []
    for name, v in GRID:
        key = _key(_cfg(name, v))
        for r in rs:
            pos, neg, hit = zip(*[o[(r, key)] for rr, o in res if rr == r])
            p, fail = _G["P"][(r, key)]
            rows.append(_rows(name, v, f"slowdown_r{r}", p, fail, pos, neg, hit))
    # rhythm threshold
    nr = 4 if smoke else 30
    thrs = (0.2, 0.3, 0.4)
    with Pool(C.WORKERS) as pool:
        _G["R"] = dict(pool.map(_setup_rhythm, [(d, t) for d in (8, 30) for t in thrs]))
    with Pool(C.WORKERS) as pool:
        resr = pool.map(_score_rhythm, [(d, k) for d in (8, 30) for k in range(nr)])
    for t in thrs:
        for d in (8, 30):
            pos, neg, hit = zip(*[o[(d, t)] for dd, o in resr if dd == d])
            p, fail = _G["R"][(d, t)]
            rows.append(_rows("rhythm_threshold", t, f"rhythm_d{d}", p, fail, pos, neg, hit))
    C.write_json("e2_sensitivity.json", dict(default=DEFAULT, n_pairs=n, rows=rows, wall_seconds=round(time.time() - t0)))
    C.write_csv("e2_sensitivity.csv", rows)
    for rw in rows:
        print(f"{rw['constant']:17s} {rw['value']:>5} {rw['case']:13s} L={rw['window']:4d} tol={rw['tolerance']:3d} auc {rw['auc']} {rw['auc_ci']} "
              f"hit {rw['hit']} chance {rw['chance_mean']}")
    print(f"total {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
