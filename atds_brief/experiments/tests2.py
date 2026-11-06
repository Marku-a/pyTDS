import sys, json
from dataclasses import replace
from multiprocessing import Pool
import numpy as np
sys.argv_ = sys.argv
sys.path.insert(0, sys.argv[1] + "/experiments/adaptive"); sys.path.insert(0, sys.argv[1])
sys.path.insert(0, sys.argv[4])
from tests import pair, score, auc, mode_tau, N
from adaptive import adaptive_params

def lagB(d):
    out = []
    for k in range(5):
        x, y = pair(20000 + d * 10 + k, 16, d, 0.6)
        ad, _ = adaptive_params([x, y], window_rule="bartlett", tolerance="calibrated")
        row = dict(d=d, k=k, window=ad.window, tol=ad.tolerance)
        for name, p in [("lag30", replace(ad, max_lag=30)), ("atds", ad)]:
            s, tau, lab = score(x, y, p); m, _ = mode_tau(tau, lab)
            row[name] = dict(found=-m, score=round(s, 1), max_lag=p.max_lag)
        out.append(row)
    return out

def tolC(r):
    tc, d, c = 2 * r, 5 * r, 0.3
    ad, _ = adaptive_params(list(pair(800000 + r, tc, d, c)), window_rule="bartlett", tolerance="calibrated")
    res = dict(r=r, window=ad.window, max_lag=ad.max_lag, tol=ad.tolerance)
    for name, p in [("tol1", replace(ad, tolerance=1)), ("atds", ad)]:
        pos, neg, hit = [], [], 0
        for k in range(20):
            s, tau, lab = score(*pair(60000 + r * 100 + k, tc, d, c), p); pos.append(s)
            if (lab == 1).sum() >= 3 and abs(np.median(tau[lab == 1]) + d) <= max(1, 0.1 * d): hit += 1
            neg.append(score(*pair(70000 + r * 100 + k, tc, d, 0.0), p)[0])
        res[name] = dict(auc=round(auc(pos, neg), 3), hit=round(hit / 20, 2), chance=round(float(np.mean(neg)), 2),
                         coupled=round(float(np.mean(pos)), 1))
    return res

if __name__ == "__main__":
    with Pool(4) as P:
        B = sum(P.map(lagB, [5, 10, 20, 30, 40, 60, 80, 120, 160]), [])
        C = P.map(tolC, [1, 2, 4, 8, 16])
    json.dump(dict(B2=B, C2=C), open(sys.argv[2], "w"))
    for d in sorted(set(b["d"] for b in B)):
        rr = [b for b in B if b["d"] == d]
        print(d, [(b["lag30"]["found"], b["lag30"]["score"], b["atds"]["found"], b["atds"]["score"]) for b in rr])
    [print(c) for c in C]
