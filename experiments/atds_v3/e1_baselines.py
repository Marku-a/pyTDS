"""E1: baselines beyond TDS settings on the slow-down benchmark (tc = 2r, d = 5r, c = 0.3, N = 20000).

Usage: ATDS_CODE=<checkout of claude/practical-mccarthy-n2mec6> python e1_baselines.py [--smoke]
Same seeds as exp4.U1; adds xcorr_full, gccphat_full and rescaled_est (r_hat from the training pairs).
"""
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

RS = [1, 2, 4, 8, 16]
_G = {}


def _setup(r):
    tc, d, c = 2 * r, 5 * r, 0.3
    train = [C.pair(800000 + r * 10 + j, tc, d, c) for j in range(5)]
    ad, fail = C.E4.atds2(train)
    tc_hat = float(np.median([C.first_below_e(s) for xy in train for s in xy]))
    r_hat = max(1, int(round(tc_hat / 2)))
    M = {"published": C.FIXED,
         "rescaled": C.TDSParams(window=60 * r, overlap=30 * r, max_lag=30 * r, tolerance=r),
         "atds": ad,
         "rescaled_est": C.TDSParams(window=60 * r_hat, overlap=30 * r_hat, max_lag=30 * r_hat, tolerance=r_hat)}
    return r, dict(M=M, fail=fail, tc_hat=tc_hat, r_hat=r_hat)


def _score_pair(args):
    r, k = args
    tc, d, c = 2 * r, 5 * r, 0.3
    xp, yp = C.pair(60000 + r * 1000 + k, tc, d, c)
    xn, yn = C.pair(70000 + r * 1000 + k, tc, d, 0.0)
    out = {}
    for name, p in _G[r]["M"].items():
        s, tau, lab = C.score(xp, yp, p)
        out[name] = (s, C.score(xn, yn, p)[0], int(C.found(tau, lab, d)), int((lab == 1).sum() >= 3))
    for name, f in (("xcorr_full", C.xcorr_full), ("gccphat_full", C.gccphat_full)):
        sp, dl, _ = f(xp, yp)
        out[name] = (sp, f(xn, yn)[0], int(C.hit_global(dl, d)), 1)
    return r, k, out


def main():
    n, smoke = C.smoke_n(100)
    C.set_results_for_smoke(smoke)
    rs = [1, 4] if smoke else RS
    t0 = time.time()
    with Pool(C.WORKERS) as pool:
        for r, g in pool.map(_setup, rs):
            _G[r] = g
        print(f"setup {time.time() - t0:.0f}s", {r: (g["r_hat"], {k: C.pdict(p) for k, p in g["M"].items()}) for r, g in _G.items()}, flush=True)
    with Pool(C.WORKERS) as pool:   # fork after _G is filled
        res = pool.map(_score_pair, [(r, k) for r in rs for k in range(n)], chunksize=2)
    out, rows = {}, []
    for r in rs:
        d = 5 * r
        per = {m: dict(pos=[], neg=[], hit=[], cond=[]) for m in res[0][2]}
        for rr, k, o in res:
            if rr == r:
                for m, (sp, sn, h, cd) in o.items():
                    per[m]["pos"].append(sp); per[m]["neg"].append(sn); per[m]["hit"].append(h); per[m]["cond"].append(cd)
        g = _G[r]
        entry = dict(r=r, n=n, d=d, calib_failed=g["fail"], tc_hat=g["tc_hat"], r_hat=g["r_hat"],
                     params={m: C.pdict(p) for m, p in g["M"].items()}, methods={}, diff_atds_minus={})
        for m, v in per.items():
            entry["methods"][m] = dict(**C.auc_row(v["pos"], v["neg"]), **C.hit_row(v["hit"]),
                                       coupled=round(float(np.mean(v["pos"])), 3), chance=round(float(np.mean(v["neg"])), 3))
            rows.append(dict(r=r, method=m, **C.flat(entry["methods"][m]), **{f"p_{k}": v2 for k, v2 in
                                                                                (C.pdict(g["M"][m]) if m in g["M"] else {}).items()}))
        for m in per:
            if m != "atds":
                c90, c95 = C.boot_diff90(per["atds"]["pos"], per["atds"]["neg"], per[m]["pos"], per[m]["neg"])
                dd = round(entry["methods"]["atds"]["auc"] - entry["methods"][m]["auc"], 3)
                entry["diff_atds_minus"][m] = dict(diff=dd, ci95=c95)
                rows.append(dict(r=r, method=f"diff_atds_minus_{m}", auc=dd, auc_ci=str(c95)))
        out[str(r)] = entry
    # consistency with the stored exp4 U1 (published = "fixed")
    ref = {x["r"]: x for x in __import__("json").load(open(os.path.join(C.ROOT, "atds_brief", "results", "exp4.json")))["U1"]}
    out["check_vs_exp4_U1"] = {str(r): {m: dict(now=out[str(r)]["methods"][m]["auc"], stored=ref[r][s]["auc"])
                                        for m, s in (("published", "fixed"), ("rescaled", "rescaled"), ("atds", "atds"))}
                               for r in rs if r in ref}
    out["wall_seconds"] = round(time.time() - t0)
    C.write_json("e1_baselines.json", out)
    C.write_csv("e1_baselines.csv", rows)
    for r in rs:
        e = out[str(r)]
        print(f"r={r} r_hat={e['r_hat']}")
        for m, v in e["methods"].items():
            print(f"  {m:14s} auc {v['auc']} {v['auc_ci']} hit {v['hit']} {v['hit_ci']} chance {v['chance']}")
        print("  diffs", e["diff_atds_minus"])
    print("check", out["check_vs_exp4_U1"]); print(f"total {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
