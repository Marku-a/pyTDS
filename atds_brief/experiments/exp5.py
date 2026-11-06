"""Round-3 experiments: Rulkov with the Methods hit rule + offset; networks with 999 rotations, incl. a common driver."""
import sys, json, time
from multiprocessing import Pool
import numpy as np
W, BD = sys.argv[1], sys.argv[2]
sys.path.insert(0, W + "/experiments/adaptive"); sys.path.insert(0, W); sys.path.insert(0, BD)
from tests import score, auc
from exp3 import found
from exp4 import atds2, calib, rulkov_oneway, MINW, FIXED
from adaptive import adaptive_params
from generators import ar1

def V1(d):
    n, nn = 6000, 20
    train = [rulkov_oneway(n, d, 0.12, np.random.default_rng(800000 + d * 10 + j)) for j in range(3)]; ad, fail = atds2(train)
    out = dict(d=d, calib_failed=fail)
    for name, p in [("fixed", FIXED), ("atds", ad)]:
        pos, neg, hit, med = [], [], 0, []
        for k in range(nn):
            xp, yp = rulkov_oneway(n, d, 0.12, np.random.default_rng(810000 + d * 100 + k)); xn, yn = rulkov_oneway(n, d, 0.0, np.random.default_rng(820000 + d * 100 + k))
            s, tau, lab = score(xp, yp, p); pos.append(s); hit += int(found(tau, lab, d))
            if (lab == 1).sum() >= 3: med.append(-float(np.median(tau[lab == 1])))
            neg.append(score(xn, yn, p)[0])
        out[name] = dict(auc=round(auc(pos, neg), 3), hit=round(hit / nn, 2), n_with_stable=len(med), median_reported=round(float(np.median(med)), 1) if med else None,
                         median_offset=round(float(np.median(med)) - d, 1) if med else None, window=p.window, max_lag=p.max_lag, tol=p.tolerance, chance=round(float(np.mean(neg)), 2))
    return out

LINKS = {"disjoint": [(0, 1, 10), (2, 3, 20), (4, 5, 15)], "driver": [(0, 1, 10), (0, 2, 20), (3, 4, 15)]}
def network(seed, kind, coupled, n=10000):
    rng = np.random.default_rng(seed); s = [ar1(n, 8, rng) for _ in range(6)]; c = 0.4 if coupled else 0.0
    for a, b, d in LINKS[kind]: s[b] = c * np.roll(s[a], d) + np.sqrt(1 - c * c) * s[b]
    return s
def V2(args):
    seed, kind, coupled = args; s = network(seed, kind, coupled); n = len(s[0])
    p0, _ = adaptive_params(s, window_rule="bartlett", tolerance="fixed", min_window=MINW)
    pairs = [(i, j) for i in range(6) for j in range(i + 1, 6)]
    p, _ = calib(p0, [(s[i], s[j]) for i, j in pairs], n_shift=4)
    rng = np.random.default_rng(seed + 7); pv = []
    for i, j in pairs:
        real = score(s[i], s[j], p)[0]; lo = max(p.window, n // 10)
        null = [score(s[i], np.roll(s[j], int(rng.integers(lo, n - lo))), p)[0] for _ in range(999)]
        pv.append((1 + sum(v >= real for v in null)) / 1000)
    pv = np.array(pv); m = len(pv); order = np.argsort(pv); thr = 0.05 * np.arange(1, m + 1) / m
    ok = pv[order] <= thr; kmax = np.nonzero(ok)[0].max() + 1 if ok.any() else 0; sig = np.zeros(m, bool); sig[order[:kmax]] = True
    true = np.array([coupled and any((i, j) == (a, b) for a, b, _ in LINKS[kind]) for i, j in pairs])
    indirect = np.array([coupled and kind == "driver" and (i, j) == (1, 2) for i, j in pairs])
    other = ~true & ~indirect
    return dict(seed=seed, kind=kind, coupled=coupled, window=p.window, tol=p.tolerance, true_found=int((sig & true).sum()), n_true=int(true.sum()),
                indirect_flagged=int((sig & indirect).sum()), indirect_p=float(pv[indirect][0]) if indirect.any() else None,
                false_other=int((sig & other).sum()), n_other=int(other.sum()), raw_fp_other=int(((pv <= 0.05) & other).sum()), n_sig=int(sig.sum()))

TASKS = {"V1": (V1, [16, 32]), "V2": (V2, [(s, k, c) for k in ("disjoint", "driver") for c in (True, False) for s in range(10)])}
if __name__ == "__main__":
    out = {}
    with Pool(4) as P:
        for t in sys.argv[4].split(","):
            t0 = time.time(); f, a = TASKS[t]; out[t] = P.map(f, a); json.dump(out, open(sys.argv[3], "w")); print(t, round(time.time() - t0), "s", flush=True)
