"""Rule v2.1 + optional longest plausible delay d_max: window >= max(61, 2*d_max + 1)."""
import sys, json, time
from multiprocessing import Pool
import numpy as np
W, BD = sys.argv[1], sys.argv[2]
sys.path.insert(0, W + "/experiments/adaptive"); sys.path.insert(0, W); sys.path.insert(0, BD)
from tests import pair, score, auc
from exp4 import calib, summarize, run_methods, FIXED, MINW
from adaptive import adaptive_params

def atds_opt(train, d_max=None):
    x, y = train[0]; mw = MINW if d_max is None else max(MINW, 2 * d_max + 1)
    p0, _ = adaptive_params([x, y], window_rule="bartlett", tolerance="fixed", min_window=mw)
    return calib(p0, train)

def W1(d):
    tc, c, n = 1, 0.5, 100
    train = [pair(1600000 + d * 10 + j, tc, d, c) for j in range(5)]
    M = {"fixed": FIXED, "atds_default": atds_opt(train)[0], "atds_dmax60": atds_opt(train, d_max=60)[0]}
    res = run_methods(M, lambda k: pair(1610000 + d * 1000 + k, tc, d, c), lambda k: pair(1620000 + d * 1000 + k, tc, d, 0.0), n, d)
    return dict(d=d, tc=tc, c=c, params={k: dict(window=p.window, max_lag=p.max_lag, tol=p.tolerance) for k, p in M.items()},
                **{k: summarize(v["pos"], v["neg"], v["hit"], v["cond"], n) for k, v in res.items()})

if __name__ == "__main__":
    with Pool(3) as P:
        t0 = time.time(); out = {"W1": P.map(W1, [20, 50])}
    json.dump(out, open(sys.argv[3], "w")); print("done", round(time.time() - t0), "s")
