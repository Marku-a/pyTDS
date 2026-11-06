import sys, json, time
from dataclasses import replace
from multiprocessing import Pool
import numpy as np
W = sys.argv[1]
sys.path.insert(0, W + "/experiments/adaptive"); sys.path.insert(0, W)
from adaptive import adaptive_params
from generators import coupled_pair
from pyTDS.params import TDSParams
from pyTDS.core import time_delay_interaction, stable_label, tds_score
N = 20000
FIXED = TDSParams()

def score(x, y, p):
    tau, _, _ = time_delay_interaction(x, y, p)
    lab = stable_label(tau, p)
    return float(tds_score(lab)), tau, lab

def auc(pos, neg):
    pos, neg = np.asarray(pos), np.asarray(neg)
    return float(((pos[:, None] > neg[None, :]).sum() + 0.5 * (pos[:, None] == neg[None, :]).sum()) / (len(pos) * len(neg)))

def pair(seed, tc, delay, c):
    return coupled_pair(N, tc, tc, delay, c, np.random.default_rng(seed))[:2]

def calibrate(p, train, seed, n_shift=20):
    rng = np.random.default_rng(seed); taus = []
    for k in range(n_shift):
        x, y = train[k % len(train)]
        lo = max(p.window, N // 10)
        taus.append(time_delay_interaction(x, np.roll(y, int(rng.integers(lo, N - lo))), p)[0])
    chosen = 0
    for tol in range(0, max(1, p.max_lag // 4) + 1):
        s = np.mean([tds_score(stable_label(t, replace(p, tolerance=tol))) for t in taus])
        if s <= 5.0: chosen = tol
        else: break
    return replace(p, tolerance=chosen)

# ---- Test A: window sweep (tolerance calibrated per window), aTDS choice marked
def testA(args):
    r, L = args
    tc, d, c = 2 * r, 5 * r, 0.3
    train = [pair(900000 + r * 10 + k, tc, d, c) for k in range(2)]
    p = TDSParams(window=L, overlap=L // 2, max_lag=(L - 1) // 2)
    p = calibrate(p, train, 7)
    pos = [score(*pair(1000 * r + k, tc, d, c), p)[0] for k in range(20)]
    neg = [score(*pair(5000 * r + k, tc, d, 0.0), p)[0] for k in range(20)]
    return dict(r=r, L=L, tol=p.tolerance, auc=round(auc(pos, neg), 3))

# ---- Test B: delay recovery, slow signals tc=16, delays beyond the fixed +-30
def mode_tau(tau, lab):
    t = tau[lab == 1] if (lab == 1).sum() >= 3 else tau
    v, cnt = np.unique(t, return_counts=True); return int(v[np.argmax(cnt)]), bool((lab == 1).sum() >= 3)
def testB(d):
    out = []
    for k in range(5):
        x, y = pair(20000 + d * 10 + k, 16, d, 0.6)
        ad, _ = adaptive_params([x, y], window_rule="bartlett", tolerance="calibrated")
        row = dict(d=d, k=k)
        for name, p in [("fixed", FIXED), ("atds", ad)]:
            s, tau, lab = score(x, y, p)
            m, st = mode_tau(tau, lab)
            row[name] = dict(found=-m, score=round(s, 1), from_stable=st, max_lag=p.max_lag)
        out.append(row)
    return out

# ---- Test C: chance level on uncoupled pairs across time scales
def testC(tc):
    train = [pair(30000 + tc * 10 + k, tc, 0, 0.0) for k in range(2)]
    ad, info = adaptive_params(list(train[0]), window_rule="bartlett", tolerance="calibrated")
    res = {}
    for name, p in [("fixed", FIXED), ("atds", ad)]:
        sc = [score(*pair(40000 + tc * 100 + k, tc, 0, 0.0), p)[0] for k in range(30)]
        res[name] = dict(mean=round(float(np.mean(sc)), 2), p95=round(float(np.percentile(sc, 95)), 2))
    return dict(tc=tc, window=ad.window, tol=ad.tolerance, **res)

# ---- Test D: overall detection AUC vs slowness, all three rules on
def testD(r):
    tc, d, c = 2 * r, 5 * r, 0.3
    train = pair(800000 + r, tc, d, c)
    ad, _ = adaptive_params(list(train), window_rule="bartlett", tolerance="calibrated")
    res = dict(r=r, window=ad.window, max_lag=ad.max_lag, tol=ad.tolerance)
    for name, p in [("fixed", FIXED), ("atds", ad)]:
        pos, neg, hit = [], [], 0
        for k in range(20):
            s, tau, lab = score(*pair(60000 + r * 100 + k, tc, d, c), p); pos.append(s)
            if (lab == 1).sum() >= 3 and abs(np.median(tau[lab == 1]) + d) <= max(1, 0.1 * d): hit += 1
            neg.append(score(*pair(70000 + r * 100 + k, tc, d, 0.0), p)[0])
        res[name] = dict(auc=round(auc(pos, neg), 3), hit=round(hit / 20, 2))
    return res

if __name__ == "__main__":
    t0 = time.time()
    with Pool(int(sys.argv[3])) as P:
        # aTDS windows for test A markers
        marks = {}
        for r in (2, 8):
            ad, _ = adaptive_params(list(pair(900000 + r * 10, 2 * r, 5 * r, 0.3)), window_rule="bartlett", tolerance="calibrated")
            marks[r] = ad.window
        A = P.map(testA, [(r, L) for r in (2, 8) for L in (30, 60, 120, 240, 480, 960, 1920)])
        print("A", time.time() - t0, flush=True)
        B = sum(P.map(testB, [5, 10, 20, 30, 40, 60, 80, 120, 160]), [])
        print("B", time.time() - t0, flush=True)
        C = P.map(testC, [1, 2, 4, 8, 16, 32, 64])
        print("C", time.time() - t0, flush=True)
        D = P.map(testD, [1, 2, 4, 8, 16])
        print("D", time.time() - t0, flush=True)
    json.dump(dict(A=A, A_marks=marks, B=B, C=C, D=D), open(sys.argv[2], "w"))
    print(marks); [print(a) for a in A]
    for d in sorted(set(b["d"] for b in B)):
        rows = [b for b in B if b["d"] == d]
        print(d, [(b["fixed"]["found"], b["atds"]["found"], b["fixed"]["score"], b["atds"]["score"]) for b in rows])
    [print(c) for c in C]; [print(x) for x in D]
