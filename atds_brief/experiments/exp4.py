"""Round-2 experiments (rule v2.1: min window 61; tolerance from pooled rotations, mean + 2 SE <= 5 %)."""
import sys, json, time
from dataclasses import replace
from multiprocessing import Pool
import numpy as np
W, BD = sys.argv[1], sys.argv[2]
sys.path.insert(0, W + "/experiments/adaptive"); sys.path.insert(0, W); sys.path.insert(0, BD)
from tests import pair, score, auc, mode_tau, N
from exp3 import boot_auc, wilson, found, fgn_pair, osc_pair, pair_bart, WGRID, TGRID
from exp3 import fgn, ar2
from adaptive import adaptive_params, acf, bartlett_factor
from generators import ar1, coupled_pair, _rulkov_f
from pyTDS.params import TDSParams
from pyTDS.core import time_delay_interaction, stable_label, tds_score
FIXED = TDSParams(); MINW = 61

def calib(p, train, n_shift=20, z=2.0, seed=0):
    """Largest tol with mean + z*SE of the chance score <= 5 %, rotations pooled over all training pairs."""
    rng = np.random.default_rng(seed); taus = []
    for x, y in train:
        n = len(x); lo = max(p.window, n // 10)
        for _ in range(n_shift): taus.append(time_delay_interaction(x, np.roll(y, int(rng.integers(lo, n - lo))), p)[0])
    tol, fail = 0, False
    for t in range(0, max(1, p.max_lag // 4) + 1):
        s = np.array([tds_score(stable_label(tt, replace(p, tolerance=t))) for tt in taus])
        ub = s.mean() + z * s.std(ddof=1) / np.sqrt(len(s))
        if ub <= 5.0: tol = t
        else:
            if t == 0: fail = True
            break
    return replace(p, tolerance=tol), fail

def atds2(train, cap=None, pairwise=False):
    """aTDS v2.1: window from the first training pair (max B, or pairwise B_xy), min 61; tolerance pooled over all training pairs."""
    x, y = train[0]
    if pairwise:
        Bxy = pair_bart(x, y); n = len(x); L = int(np.clip(round(30 * Bxy), MINW, max(MINW, 0.1 * n)))
        p0 = TDSParams(window=L, overlap=L // 2, max_lag=(L - 1) // 2)
    else:
        p0, _ = adaptive_params([x, y], window_rule="bartlett", tolerance="fixed", min_window=MINW, max_window=cap)
    return calib(p0, train)

def boot_diff90(pa, na, pb, nb, B=1000, seed=1):
    rng = np.random.default_rng(seed); v = []; pa, na, pb, nb = map(np.asarray, (pa, na, pb, nb))
    for _ in range(B):
        i, j = rng.integers(0, len(pa), len(pa)), rng.integers(0, len(na), len(na)); v.append(auc(pa[i], na[j]) - auc(pb[i], nb[j]))
    return [round(float(np.percentile(v, 5)), 3), round(float(np.percentile(v, 95)), 3)], [round(float(np.percentile(v, 2.5)), 3), round(float(np.percentile(v, 97.5)), 3)]

def summarize(pos, neg, hits, conds, n):
    k = sum(conds); hc = sum(h for h, c in zip(hits, conds) if c)
    return dict(auc=round(auc(pos, neg), 3), auc_ci=boot_auc(pos, neg), hit=round(sum(hits) / n, 2), hit_ci=wilson(sum(hits), n),
                hit_cond=round(hc / k, 2) if k else None, n_cond=k, coupled=round(float(np.mean(pos)), 2), chance=round(float(np.mean(neg)), 2),
                chance_p95=round(float(np.percentile(neg, 95)), 2))

def run_methods(M, gen_pos, gen_neg, n, d):
    res = {k: dict(pos=[], neg=[], hit=[], cond=[]) for k in M}
    for k in range(n):
        xp, yp = gen_pos(k); xn, yn = gen_neg(k)
        for name, p in M.items():
            s, tau, lab = score(xp, yp, p); res[name]["pos"].append(s); res[name]["hit"].append(int(found(tau, lab, d))); res[name]["cond"].append(int((lab == 1).sum() >= 3))
            res[name]["neg"].append(score(xn, yn, p)[0])
    return res

# U1: overall test (Fig 6, 7)
def U1(r, n=100):
    tc, d, c = 2 * r, 5 * r, 0.3
    train = [pair(800000 + r * 10 + j, tc, d, c) for j in range(5)]
    ad, fail = atds2(train)
    resc = TDSParams(window=60 * r, overlap=30 * r, max_lag=30 * r, tolerance=r)
    trP = [pair(910000 + r * 100 + k, tc, d, c) for k in range(20)]; trN = [pair(920000 + r * 100 + k, tc, d, 0.0) for k in range(20)]
    best = (-1, None)
    for L in WGRID:
        p0 = TDSParams(window=L, overlap=L // 2, max_lag=(L - 1) // 2)
        tp = [time_delay_interaction(x, y, p0)[0] for x, y in trP]; tn = [time_delay_interaction(x, y, p0)[0] for x, y in trN]
        for tol in TGRID:
            if tol > p0.max_lag: continue
            pt = replace(p0, tolerance=tol); a = auc([tds_score(stable_label(t, pt)) for t in tp], [tds_score(stable_label(t, pt)) for t in tn])
            if a > best[0]: best = (a, pt)
    M = {"fixed": FIXED, "rescaled": resc, "atds": ad, "atds_tol1": replace(ad, tolerance=1), "oracle": best[1]}
    res = run_methods(M, lambda k: pair(60000 + r * 1000 + k, tc, d, c), lambda k: pair(70000 + r * 1000 + k, tc, d, 0.0), n, d)
    out = dict(r=r, n=n, calib_failed=fail, params={k: dict(window=p.window, max_lag=p.max_lag, tol=p.tolerance) for k, p in M.items()})
    for name, v in res.items(): out[name] = summarize(v["pos"], v["neg"], v["hit"], v["cond"], n)
    for other in ("fixed", "rescaled", "oracle", "atds_tol1"):
        c90, c95 = boot_diff90(res["atds"]["pos"], res["atds"]["neg"], res[other]["pos"], res[other]["neg"])
        out["diff_" + other] = dict(ci90=c90, ci95=c95)
    return out

# U2: window sweep (Fig 3), 100 per class
def U2(args):
    r, L = args; tc, d, c = 2 * r, 5 * r, 0.3
    train = [pair(900000 + r * 10 + k, tc, d, c) for k in range(5)]
    p, _ = calib(TDSParams(window=L, overlap=L // 2, max_lag=(L - 1) // 2), train)
    pos = [score(*pair(1000 * r + k, tc, d, c), p)[0] for k in range(100)]; neg = [score(*pair(5000 * r + k, tc, d, 0.0), p)[0] for k in range(100)]
    return dict(r=r, L=L, tol=p.tolerance, auc=round(auc(pos, neg), 3), auc_ci=boot_auc(pos, neg))
def U2marks(r):
    tc, d, c = 2 * r, 5 * r, 0.3; train = [pair(900000 + r * 10 + k, tc, d, c) for k in range(5)]
    p, _ = atds2(train)
    pos = [score(*pair(1000 * r + k, tc, d, c), p)[0] for k in range(100)]; neg = [score(*pair(5000 * r + k, tc, d, 0.0), p)[0] for k in range(100)]
    return dict(r=r, L=p.window, tol=p.tolerance, auc=round(auc(pos, neg), 3), auc_ci=boot_auc(pos, neg))

# U3: block coupling, fixed +-2 criterion, three block lengths
def U3(args):
    neff, bl = args; tc, d, c, n = 8, 20, 0.6, 40000
    x0, y0, _ = coupled_pair(n, tc, tc, d, c, np.random.default_rng(5), block_len=bl)
    B = bartlett_factor(x0); L = max(20, int(round(neff * B)))
    p, _ = calib(TDSParams(window=L, overlap=L // 2, max_lag=(L - 1) // 2), [(x0, y0)])
    ba = []
    for k in range(20):
        x, y, m = coupled_pair(n, tc, tc, d, c, np.random.default_rng(100 + k), block_len=bl)
        tau, _, _ = time_delay_interaction(x, y, p); lab = stable_label(tau, p).astype(bool)
        starts = np.arange(0, n - L + 1, p.overlap)[:len(tau)]; truth = np.array([m[s:s + L].mean() >= 0.5 for s in starts])
        corr = lab & (np.abs(tau + d) <= 2)
        ba.append((corr[truth].mean() + (1 - lab[~truth].mean())) / 2)
    return dict(neff=neff, block=bl, B=round(B, 2), window=L, tol=p.tolerance, bal_acc=round(float(np.mean(ba)), 3), sd=round(float(np.std(ba)), 3))

# U4: lag test (Fig 4) with v2.1
def U4(d):
    train = [pair(30000 + d * 10 + j, 16, d, 0.6) for j in range(5)]; ad, _ = atds2(train); rows = []
    for k in range(10):
        x, y = pair(20000 + d * 10 + k, 16, d, 0.6); row = dict(d=d, k=k, window=ad.window, tol=ad.tolerance)
        for name, p in [("lag30", replace(ad, max_lag=30)), ("atds", ad)]:
            s, tau, lab = score(x, y, p); m, _ = mode_tau(tau, lab); row[name] = dict(found=-m, score=round(s, 1), max_lag=p.max_lag)
        rows.append(row)
    return rows

# U5: long memory, weaker coupling
def U5(H):
    d, c, n = 20, 0.2, 40
    train = [fgn_pair(500000 + int(H * 100) + j, H, d, c) for j in range(5)]; ad, fail = atds2(train)
    old, _ = adaptive_params(list(train[0]), window_rule="bartlett", tolerance="calibrated")
    M = {"fixed": FIXED, "atds": ad, "atds_v2_min20": old}
    res = run_methods(M, lambda k: fgn_pair(510000 + int(H * 1000) + k, H, d, c), lambda k: fgn_pair(520000 + int(H * 1000) + k, H, d, 0.0), n, d)
    return dict(H=H, calib_failed=fail, params={k: dict(window=p.window, max_lag=p.max_lag, tol=p.tolerance) for k, p in M.items()},
                **{k: summarize(v["pos"], v["neg"], v["hit"], v["cond"], n) for k, v in res.items()})

# U6: rhythm, with sign of C at tau0
def U6(d):
    c, n = 0.6, 30; train = [osc_pair(700000 + d + j, d, c) for j in range(5)]; ad, fail = atds2(train)
    info = adaptive_params(list(train[0]), window_rule="bartlett", tolerance="fixed")[1]
    out = dict(d=d, period=info["periods"], calib_failed=fail)
    for name, p in [("fixed", FIXED), ("atds", ad)]:
        pos, neg, hit, rep, sgn = [], [], 0, [], []
        for k in range(n):
            x, y = osc_pair(710000 + d * 100 + k, d, c); s, tau, lab = score(x, y, p); pos.append(s); hit += int(found(tau, lab, d)); rep.append(-mode_tau(tau, lab)[0])
            # sign of the correlation at the reported lag, in the first stable window
            if (lab == 1).any():
                i = int(np.nonzero(lab == 1)[0][0]); st = i * p.overlap; a, b = x[st:st + p.window], y[st:st + p.window]
                a = (a - a.mean()) / a.std(); b = (b - b.mean()) / b.std(); kk = int(tau[i]); sgn.append(float(np.mean(np.roll(a, -kk) * b)))
            neg.append(score(*osc_pair(720000 + d * 100 + k, d, 0.0), p)[0])
        v, cn = np.unique(rep, return_counts=True)
        out[name] = dict(auc=round(auc(pos, neg), 3), hit=round(hit / n, 2), chance=round(float(np.mean(neg)), 2), reported_mode=int(v[np.argmax(cn)]),
                         share_negative_C=round(float(np.mean(np.array(sgn) < 0)), 2) if sgn else None, window=p.window, max_lag=p.max_lag, tol=p.tolerance)
    return out

# U7: Rulkov, one-way coupling (neuron 0 drives neuron 1)
def rulkov_oneway(n, delay, g, rng, mu=0.001, alpha=(4.5, 4.1), sigma=(0.01, -0.01), beta=1.0, sigma_e=1.0, extra=1000):
    disc = 50 + extra; T = n + disc; x = [np.zeros(T), np.zeros(T)]; y = [np.zeros(T), np.zeros(T)]
    x0 = rng.uniform(-1.5, -0.5, 2); y0 = rng.uniform(-4.5, -3.5, 2)
    for k in range(2): x[k][0], y[k][0] = x0[k], y0[k]
    for _ in range(2):
        for i in range(1, T):
            for k in range(2):
                if k == 1 and delay < i: b = g * beta * (x[0][i - 1 - delay] - x[1][i - 1]); s = g * sigma_e * (x[0][i - 1 - delay] - x[1][i - 1])
                else: b = s = 0.0
                x[k][i] = _rulkov_f(x[k][i - 1], y[k][i - 1] + b, alpha[k]); y[k][i] = y[k][i - 1] - mu * (x[k][i - 1] + 1) + mu * sigma[k] + mu * s
        for k in range(2): y[k][0] = np.median(y[k]); x[k][0] = x0[k]
    return x[0][disc:], x[1][disc:]
def U7(d):
    n, nn = 6000, 20
    train = [rulkov_oneway(n, d, 0.12, np.random.default_rng(800000 + d * 10 + j)) for j in range(3)]; ad, fail = atds2(train)
    out = dict(d=d, calib_failed=fail, periods=adaptive_params(list(train[0]), window_rule="bartlett", tolerance="fixed")[1]["periods"])
    for name, p in [("fixed", FIXED), ("atds", ad)]:
        pos, neg, rep, hit = [], [], [], 0
        for k in range(nn):
            xp, yp = rulkov_oneway(n, d, 0.12, np.random.default_rng(810000 + d * 100 + k)); xn, yn = rulkov_oneway(n, d, 0.0, np.random.default_rng(820000 + d * 100 + k))
            s, tau, lab = score(xp, yp, p); pos.append(s); m, _ = mode_tau(tau, lab); rep.append(-m); hit += int(-m == d or abs(-m - d) <= max(1, 0.1 * d))
            neg.append(score(xn, yn, p)[0])
        v, cn = np.unique(rep, return_counts=True)
        out[name] = dict(auc=round(auc(pos, neg), 3), hit=round(hit / nn, 2), reported_mode=int(v[np.argmax(cn)]), window=p.window, max_lag=p.max_lag, tol=p.tolerance, chance=round(float(np.mean(neg)), 2))
    return out

# U8: EEG, parameters from first half, scored on second half
def U8(_):
    a = np.loadtxt("/home/user/pyTDS/data/data-delta.txt")[:, 1]; b = np.loadtxt("/home/user/pyTDS/data/data-sigma.txt")[:, 1]
    n = min(len(a), len(b)); h = n // 2; A1, B1, A2, B2 = a[:h], b[:h], a[h:n], b[h:n]
    out = dict(n_half=h)
    variants = [("fixed", FIXED), ("fixed_end", replace(FIXED, stability_anchor="end")), ("atds", atds2([(A1, B1)])[0]), ("atds_cap300", atds2([(A1, B1)], cap=300)[0])]
    for name, p in variants:
        s, tau, lab = score(A2, B2, p); rng = np.random.default_rng(9); null = []
        for _ in range(100):
            lo = max(p.window, h // 10); null.append(score(A2, np.roll(B2, int(rng.integers(lo, (n - h) - lo))), p)[0])
        st = tau[lab == 1]
        out[name] = dict(window=p.window, max_lag=p.max_lag, tol=p.tolerance, score=round(s, 1), chance=round(float(np.mean(null)), 1), chance_p95=round(float(np.percentile(null, 95)), 1),
                         p_emp=round((1 + sum(v >= s for v in null)) / 101, 3), n_windows=len(tau), median_stable_tau=float(np.median(st)) if len(st) else None)
    return out

# U9: mixed speeds: (a) system with separate slow and fast pairs, (b) a coupled fast-slow pair
def U9(_):
    n, nt = N, 40
    def system(seed, coupled):
        rng = np.random.default_rng(seed); cc = 0.4 if coupled else 0.0
        xs, ys, _ = coupled_pair(n, 32, 32, 40, cc, rng); xf, yf, _ = coupled_pair(n, 2, 2, 5, cc, rng); return xs, ys, xf, yf
    def fs_pair(seed, c):
        rng = np.random.default_rng(seed); x = ar1(n, 32, rng); e = ar1(n, 2, rng); return x, c * np.roll(x, 10) + np.sqrt(1 - c * c) * e
    out = {}
    trS = [system(990000 + j, True) for j in range(5)]
    sysw, _ = adaptive_params(list(trS[0]), window_rule="bartlett", tolerance="fixed", min_window=MINW); sysw, _ = calib(sysw, [(s[0], s[1]) for s in trS] + [(s[2], s[3]) for s in trS])
    rules = {"system_wide": {"slow": sysw, "fast": sysw},
             "per_pair_maxB": {"slow": atds2([(s[0], s[1]) for s in trS])[0], "fast": atds2([(s[2], s[3]) for s in trS])[0]},
             "pairwise_Bxy": {"slow": atds2([(s[0], s[1]) for s in trS], pairwise=True)[0], "fast": atds2([(s[2], s[3]) for s in trS], pairwise=True)[0]},
             "fixed": {"slow": FIXED, "fast": FIXED}}
    for rn, P in rules.items():
        o = {}
        for which, (i, j, d) in {"slow": (0, 1, 40), "fast": (2, 3, 5)}.items():
            p = P[which]; res = run_methods({"m": p}, lambda k: system(1000000 + k, True)[i:j + 1:j - i] if False else (system(1000000 + k, True)[i], system(1000000 + k, True)[j]),
                                            lambda k: (system(2000000 + k, False)[i], system(2000000 + k, False)[j]), nt, d)["m"]
            o[which] = dict(window=p.window, max_lag=p.max_lag, tol=p.tolerance, **summarize(res["pos"], res["neg"], res["hit"], res["cond"], nt))
        out[rn] = o
    trF = [fs_pair(995000 + j, 0.5) for j in range(5)]
    M = {"per_pair_maxB": atds2(trF)[0], "pairwise_Bxy": atds2(trF, pairwise=True)[0], "fixed": FIXED}
    res = run_methods(M, lambda k: fs_pair(1100000 + k, 0.5), lambda k: fs_pair(1200000 + k, 0.0), nt, 10)
    out["fast_slow_pair"] = {k: dict(window=M[k].window, max_lag=M[k].max_lag, tol=M[k].tolerance, **summarize(v["pos"], v["neg"], v["hit"], v["cond"], nt)) for k, v in res.items()}
    out["fast_slow_B"] = dict(B_x=round(bartlett_factor(trF[0][0]), 2), B_y=round(bartlett_factor(trF[0][1]), 2), B_xy=round(pair_bart(*trF[0]), 2))
    return out

# U10: nonstationary: memory changes halfway (tc 4 -> 16), delay 10, coupling 0.4
def ns_pair(seed, c, n=N):
    rng = np.random.default_rng(seed); h = n // 2
    x = np.concatenate([ar1(h, 4, rng), ar1(n - h, 16, rng)]); e = np.concatenate([ar1(h, 4, rng), ar1(n - h, 16, rng)])
    return x, c * np.roll(x, 10) + np.sqrt(1 - c * c) * e
def U10(_):
    n, nt, h = N, 40, N // 2; train = [ns_pair(1300000 + j, 0.4) for j in range(5)]
    whole, _ = atds2(train); first, _ = atds2([(x[:h], y[:h]) for x, y in train]); second, _ = atds2([(x[h:], y[h:]) for x, y in train])
    out = dict(params=dict(whole=whole.window, whole_tol=whole.tolerance, first=first.window, first_tol=first.tolerance, second=second.window, second_tol=second.tolerance))
    for half, sl in (("first", slice(0, h)), ("second", slice(h, n))):
        M = {"fixed": FIXED, "atds_whole": whole, "atds_per_segment": first if half == "first" else second}
        res = run_methods(M, lambda k: tuple(v[sl] for v in ns_pair(1400000 + k, 0.4)), lambda k: tuple(v[sl] for v in ns_pair(1500000 + k, 0.0)), nt, 10)
        out[half] = {k: summarize(v["pos"], v["neg"], v["hit"], v["cond"], nt) for k, v in res.items()}
    return out

# U11: network with per-link rotation test + Benjamini-Hochberg FDR
def network(seed, coupled, n=10000):
    rng = np.random.default_rng(seed); s = [ar1(n, 8, rng) for _ in range(6)]; c = 0.4 if coupled else 0.0
    for a, b, d in ((0, 1, 10), (2, 3, 20), (4, 5, 15)): s[b] = c * np.roll(s[a], d) + np.sqrt(1 - c * c) * s[b]
    return s
def U11(args):
    seed, coupled = args; s = network(seed, coupled); n = len(s[0])
    p0, _ = adaptive_params(s, window_rule="bartlett", tolerance="fixed", min_window=MINW)
    pairs = [(i, j) for i in range(6) for j in range(i + 1, 6)]
    p, _ = calib(p0, [(s[i], s[j]) for i, j in pairs], n_shift=4)
    rng = np.random.default_rng(seed + 7); pv = []
    for i, j in pairs:
        real = score(s[i], s[j], p)[0]; lo = max(p.window, n // 10)
        null = [score(s[i], np.roll(s[j], int(rng.integers(lo, n - lo))), p)[0] for _ in range(99)]
        pv.append((1 + sum(v >= real for v in null)) / 100)
    pv = np.array(pv); order = np.argsort(pv); m = len(pv); thr = 0.05 * (np.arange(1, m + 1)) / m
    passed = pv[order] <= thr; kmax = np.nonzero(passed)[0].max() + 1 if passed.any() else 0
    sig = np.zeros(m, bool); sig[order[:kmax]] = True
    true = np.array([(i, j) in ((0, 1), (2, 3), (4, 5)) and coupled for i, j in pairs])
    return dict(seed=seed, coupled=coupled, window=p.window, tol=p.tolerance, true_found=int((sig & true).sum()), false_links=int((sig & ~true).sum()),
                n_sig=int(sig.sum()), raw_p05_false=int(((pv <= 0.05) & ~true).sum()))

# U12: Figure 1/5 example with v2.1 (training pairs separate)
def U12(_):
    x, y, _ = coupled_pair(30000, 8, 8, 20, 0.6, np.random.default_rng(7))
    ad, _ = atds2([tuple(coupled_pair(30000, 8, 8, 20, 0.6, np.random.default_rng(8 + 100 * j))[:2]) for j in range(5)])
    ex = {}
    for name, p in [("fixed", FIXED), ("atds", ad)]:
        tau, t, _ = time_delay_interaction(x, y, p); lab = stable_label(tau, p); m = t <= 12000
        ex[name] = dict(params=dict(window=p.window, step=p.overlap, max_lag=p.max_lag, tol=p.tolerance), score=round(float(tds_score(lab)), 1),
                        near_true=round(float(np.mean(np.abs(tau + 20) <= 2)) * 100, 1), t=t[m].astype(int).tolist(), tau=tau[m].astype(int).tolist(), stable=lab[m].astype(int).tolist())
    r2 = np.random.default_rng(3); T_ = len(x); lo, hi = max(ad.window, T_ // 10), T_ - max(ad.window, T_ // 10)
    nt = [time_delay_interaction(x, np.roll(y, int(r2.integers(lo, hi))), ad)[0] for _ in range(20)]; rt = time_delay_interaction(x, y, ad)[0]
    curve = []
    for tol in range(31):
        sc = np.array([tds_score(stable_label(t, replace(ad, tolerance=tol))) for t in nt])
        curve.append(dict(tol=tol, null=round(float(sc.mean()), 2), ub=round(float(sc.mean() + 2 * sc.std(ddof=1) / np.sqrt(len(sc))), 2), real=round(float(tds_score(stable_label(rt, replace(ad, tolerance=tol)))), 2)))
    return dict(example=ex, example_info=dict(window=ad.window, tol=ad.tolerance), tol_curve=curve)

TASKS = {"U12": (U12, [0]), "U8": (U8, [0]), "U6": (U6, [8, 30]), "U7": (U7, [16, 32]), "U5": (U5, [0.6, 0.75, 0.9]), "U4": (U4, [5, 10, 20, 30, 40, 60, 80, 120, 160]),
         "U10": (U10, [0]), "U9": (U9, [0]), "U3": (U3, [(ne, bl) for bl in (1000, 2000, 4000) for ne in (5, 10, 15, 20, 30, 45, 60, 90, 120)]),
         "U11": (U11, [(s, c) for c in (True, False) for s in range(10)]), "U2m": (U2marks, [2, 8]), "U2": (U2, [(r, L) for r in (2, 8) for L in WGRID]), "U1": (U1, [1, 2, 4, 8, 16])}
if __name__ == "__main__":
    which = list(TASKS) if sys.argv[4] == "all" else sys.argv[4].split(",")
    try: out = json.load(open(sys.argv[3]))
    except Exception: out = {}
    with Pool(4) as P:
        for t in which:
            t0 = time.time(); f, args = TASKS[t]; r = P.map(f, args); out[t] = sum(r, []) if t == "U4" else r
            json.dump(out, open(sys.argv[3], "w")); print(t, "done", round(time.time() - t0), "s", flush=True)
