"""Experiments answering the referee report. Usage: python -I exp3.py <worktree> <builddir> <out.json> <task|all>"""
import sys, json, time
from dataclasses import replace
from multiprocessing import Pool
import numpy as np
from scipy.signal import lfilter
W, BD = sys.argv[1], sys.argv[2]
sys.path.insert(0, W + "/experiments/adaptive"); sys.path.insert(0, W); sys.path.insert(0, BD)
from tests import pair, score, auc, calibrate, mode_tau, N
from adaptive import adaptive_params, acf, bartlett_factor
from generators import ar1, coupled_pair, rulkov_pair
from pyTDS.params import TDSParams
from pyTDS.core import time_delay_interaction, stable_label, tds_score
FIXED = TDSParams()
ATDS = dict(window_rule="bartlett", tolerance="calibrated")

def boot_auc(pos, neg, B=1000, seed=0):
    rng = np.random.default_rng(seed); pos, neg = np.asarray(pos), np.asarray(neg); v = []
    for _ in range(B): v.append(auc(pos[rng.integers(0, len(pos), len(pos))], neg[rng.integers(0, len(neg), len(neg))]))
    return [round(float(np.percentile(v, 2.5)), 3), round(float(np.percentile(v, 97.5)), 3)]
def boot_diff(posA, negA, posB, negB, B=1000, seed=1):
    rng = np.random.default_rng(seed); v = []; n1, n0 = len(posA), len(negA)
    for _ in range(B):
        i, j = rng.integers(0, n1, n1), rng.integers(0, n0, n0)
        v.append(auc(np.asarray(posA)[i], np.asarray(negA)[j]) - auc(np.asarray(posB)[i], np.asarray(negB)[j]))
    return [round(float(np.percentile(v, 2.5)), 3), round(float(np.percentile(v, 97.5)), 3)]
def wilson(k, n):
    if n == 0: return [0, 0]
    p, z = k / n, 1.96; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d; h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [round(c - h, 3), round(c + h, 3)]
def found(tau, lab, d):
    return (lab == 1).sum() >= 3 and abs(np.median(tau[lab == 1]) + d) <= max(1, 0.1 * abs(d))

# ---------- T1: overall test with fair baselines, n = 100 per class, CIs ----------
WGRID = [30, 60, 120, 240, 480, 960, 1920]
TGRID = [0, 1, 2, 3, 5, 8, 12, 18, 27, 40, 60, 90, 135]
def T1(r, n=100):
    tc, d, c = 2 * r, 5 * r, 0.3
    ad, _ = adaptive_params(list(pair(800000 + r, tc, d, c)), **ATDS)
    resc = TDSParams(window=60 * r, overlap=30 * r, max_lag=30 * r, tolerance=r)
    # oracle: grid chosen by AUC on separate training pairs
    trP = [pair(910000 + r * 100 + k, tc, d, c) for k in range(20)]; trN = [pair(920000 + r * 100 + k, tc, d, 0.0) for k in range(20)]
    best = (-1, None)
    for L in WGRID:
        p0 = TDSParams(window=L, overlap=L // 2, max_lag=(L - 1) // 2)
        tp = [time_delay_interaction(x, y, p0)[0] for x, y in trP]; tn = [time_delay_interaction(x, y, p0)[0] for x, y in trN]
        for tol in TGRID:
            if tol > p0.max_lag: continue
            pt = replace(p0, tolerance=tol)
            a = auc([tds_score(stable_label(t, pt)) for t in tp], [tds_score(stable_label(t, pt)) for t in tn])
            if a > best[0]: best = (a, pt)
    orc = best[1]
    M = {"fixed": FIXED, "rescaled": resc, "atds": ad, "atds_tol1": replace(ad, tolerance=1), "oracle": orc}
    res = {k: dict(pos=[], neg=[], hit=0) for k in M}
    for k in range(n):
        xp, yp = pair(60000 + r * 1000 + k, tc, d, c); xn, yn = pair(70000 + r * 1000 + k, tc, d, 0.0)
        for name, p in M.items():
            s, tau, lab = score(xp, yp, p); res[name]["pos"].append(s); res[name]["hit"] += int(found(tau, lab, d))
            res[name]["neg"].append(score(xn, yn, p)[0])
    out = dict(r=r, n=n, params={k: dict(window=p.window, max_lag=p.max_lag, tol=p.tolerance) for k, p in M.items()})
    for name, v in res.items():
        out[name] = dict(auc=round(auc(v["pos"], v["neg"]), 3), auc_ci=boot_auc(v["pos"], v["neg"]), hit=round(v["hit"] / n, 2),
                         hit_ci=wilson(v["hit"], n), coupled=round(float(np.mean(v["pos"])), 2), chance=round(float(np.mean(v["neg"])), 2),
                         chance_p95=round(float(np.percentile(v["neg"], 95)), 2))
    for other in ("fixed", "rescaled", "oracle", "atds_tol1"):
        out["diff_atds_minus_" + other] = boot_diff(res["atds"]["pos"], res["atds"]["neg"], res[other]["pos"], res[other]["neg"])
    return out

# ---------- T2: lag isolated, parameters from a separate training pair ----------
def T2(d):
    ad, _ = adaptive_params(list(pair(30000 + d * 10, 16, d, 0.6)), **ATDS)
    rows = []
    for k in range(10):
        x, y = pair(20000 + d * 10 + k, 16, d, 0.6)
        row = dict(d=d, k=k, window=ad.window, tol=ad.tolerance)
        for name, p in [("lag30", replace(ad, max_lag=30)), ("atds", ad)]:
            s, tau, lab = score(x, y, p); m, _ = mode_tau(tau, lab); row[name] = dict(found=-m, score=round(s, 1), max_lag=p.max_lag)
        rows.append(row)
    return rows

# ---------- T3: block coupling, time-localised detection vs independent samples per window ----------
def T3(args):
    neff, c = args
    tc, d, n, bl = 8, 20, 40000, 2000
    rng0 = np.random.default_rng(5)
    x0, y0, _ = coupled_pair(n, tc, tc, d, c, rng0, block_len=bl)
    B = bartlett_factor(x0); L = max(20, int(round(neff * B)))
    p = TDSParams(window=L, overlap=L // 2, max_lag=(L - 1) // 2)
    rng = np.random.default_rng(6); taus = []
    for _ in range(20):
        lo = max(L, n // 10); taus.append(time_delay_interaction(x0, np.roll(y0, int(rng.integers(lo, n - lo))), p)[0])
    tol = 0
    for t_ in range(0, max(1, p.max_lag // 4) + 1):
        if np.mean([tds_score(stable_label(tt, replace(p, tolerance=t_))) for tt in taus]) <= 5: tol = t_
        else: break
    p = replace(p, tolerance=tol); ba = []
    for k in range(10):
        x, y, m = coupled_pair(n, tc, tc, d, c, np.random.default_rng(100 + k), block_len=bl)
        tau, tv, _ = time_delay_interaction(x, y, p); lab = stable_label(tau, p)
        starts = np.arange(0, n - L + 1, p.overlap)[:len(tau)]
        truth = np.array([m[s:s + L].mean() >= 0.5 for s in starts])
        corr = lab.astype(bool) & (np.abs(tau + d) <= max(tol, 1))  # stable AND at the true delay
        tpr = corr[truth].mean() if truth.any() else 0; tnr = 1 - lab.astype(bool)[~truth].mean() if (~truth).any() else 0
        ba.append((tpr + tnr) / 2)
    return dict(neff=neff, c=c, B=round(B, 2), window=L, tol=tol, bal_acc=round(float(np.mean(ba)), 3), sd=round(float(np.std(ba)), 3))

# ---------- T4: long-range-correlated signals (fGn) ----------
def fgn(n, H, rng):
    k = np.arange(n + 1); g = 0.5 * ((k + 1) ** (2 * H) - 2 * k ** (2 * H) + np.abs(k - 1) ** (2 * H))
    row = np.concatenate([g, g[-2:0:-1]]); lam = np.fft.fft(row).real; lam[lam < 0] = 0
    m = len(row); w = rng.standard_normal(m) + 1j * rng.standard_normal(m)
    z = np.fft.fft(np.sqrt(lam / m) * w)[:n].real; return (z - z.mean()) / z.std()
def fgn_pair(seed, H, d, c, n=N):
    rng = np.random.default_rng(seed); x = fgn(n, H, rng); e = fgn(n, H, rng)
    return x, c * np.roll(x, d) + np.sqrt(1 - c * c) * e
def T4(H):
    d, c, n = 20, 0.5, 40
    ad, info = adaptive_params(list(fgn_pair(500000 + int(H * 100), H, d, c)), **ATDS)
    res = {}
    for name, p in [("fixed", FIXED), ("atds", ad)]:
        pos, neg, hit = [], [], 0
        for k in range(n):
            s, tau, lab = score(*fgn_pair(510000 + int(H * 1000) + k, H, d, c), p); pos.append(s); hit += int(found(tau, lab, d))
            neg.append(score(*fgn_pair(520000 + int(H * 1000) + k, H, d, 0.0), p)[0])
        res[name] = dict(auc=round(auc(pos, neg), 3), auc_ci=boot_auc(pos, neg), hit=round(hit / n, 2), chance=round(float(np.mean(neg)), 2),
                         chance_p95=round(float(np.percentile(neg, 95)), 2), coupled=round(float(np.mean(pos)), 2))
    Bn = {}
    for nn in (5000, 20000, 80000):
        Bn[nn] = round(float(np.median([bartlett_factor(fgn(nn, H, np.random.default_rng(600000 + j))) for j in range(5)])), 2)
    return dict(H=H, window=ad.window, max_lag=ad.max_lag, tol=ad.tolerance, calib_failed=bool(info.get("calibration_failed")), B_vs_length=Bn, **res)

# ---------- T5: rhythmic signals (AR(2), period 40) ----------
def ar2(n, P, rho, rng):
    a1, a2 = 2 * rho * np.cos(2 * np.pi / P), -rho * rho
    x = lfilter([1.0], [1.0, -a1, -a2], rng.standard_normal(n + 500))[500:]; return (x - x.mean()) / x.std()
def osc_pair(seed, d, c, n=N, P=40):
    rng = np.random.default_rng(seed); x = ar2(n, P, 0.97, rng); e = ar2(n, P, 0.97, rng)
    return x, c * np.roll(x, d) + np.sqrt(1 - c * c) * e
def T5(d):
    c, n = 0.6, 30
    ad, info = adaptive_params(list(osc_pair(700000 + d, d, c)), **ATDS)
    res = {}
    for name, p in [("fixed", FIXED), ("atds", ad)]:
        pos, neg, hit, rep = [], [], 0, []
        for k in range(n):
            s, tau, lab = score(*osc_pair(710000 + d * 100 + k, d, c), p); pos.append(s); hit += int(found(tau, lab, d)); rep.append(-mode_tau(tau, lab)[0])
            neg.append(score(*osc_pair(720000 + d * 100 + k, d, 0.0), p)[0])
        vals, cnt = np.unique(rep, return_counts=True)
        res[name] = dict(auc=round(auc(pos, neg), 3), hit=round(hit / n, 2), chance=round(float(np.mean(neg)), 2), coupled=round(float(np.mean(pos)), 2),
                         most_reported=int(vals[np.argmax(cnt)]), max_lag=p.max_lag, tol=p.tolerance, window=p.window)
    return dict(d=d, period_detected=[p for p in info["periods"]], calib_failed=bool(info.get("calibration_failed")), **res)

# ---------- T6: Rulkov neurons ----------
def T6(d):
    n, nn = 6000, 20
    rng = np.random.default_rng(800000 + d); ad, info = adaptive_params(list(rulkov_pair(n, d, 0.12, rng)), **ATDS)
    res = {}
    pos_s = {k: [] for k in ("fixed", "atds")}; neg_s = {k: [] for k in ("fixed", "atds")}; hit = {k: 0 for k in ("fixed", "atds")}; rep = {k: [] for k in ("fixed", "atds")}
    for k in range(nn):
        xp, yp = rulkov_pair(n, d, 0.12, np.random.default_rng(810000 + d * 100 + k)); xn, yn = rulkov_pair(n, d, 0.0, np.random.default_rng(820000 + d * 100 + k))
        for name, p in [("fixed", FIXED), ("atds", ad)]:
            s, tau, lab = score(xp, yp, p); pos_s[name].append(s); hit[name] += int(found(tau, lab, -d)); rep[name].append(int(mode_tau(tau, lab)[0]))
            neg_s[name].append(score(xn, yn, p)[0])
    for name, p in [("fixed", FIXED), ("atds", ad)]:
        vals, cnt = np.unique(rep[name], return_counts=True)
        res[name] = dict(auc=round(auc(pos_s[name], neg_s[name]), 3), hit=round(hit[name] / nn, 2), most_reported=int(vals[np.argmax(cnt)]),
                         window=p.window, max_lag=p.max_lag, tol=p.tolerance, chance=round(float(np.mean(neg_s[name])), 2))
    return dict(d=d, periods=info["periods"], calib_failed=bool(info.get("calibration_failed")), **res)

# ---------- T7: EEG benchmark pair (delta, sigma) ----------
def T7(_):
    a = np.loadtxt("/home/user/pyTDS/data/data-delta.txt")[:, 1]; b = np.loadtxt("/home/user/pyTDS/data/data-sigma.txt")[:, 1]
    n = min(len(a), len(b)); a, b = a[:n], b[:n]
    ad, info = adaptive_params([a, b], **ATDS); out = dict(n=n, B=[round(v, 1) for v in info["bartlett"]])
    for name, p in [("fixed", FIXED), ("fixed_end", replace(FIXED, stability_anchor="end")), ("atds", ad), ("atds_cap300", adaptive_params([a, b], max_window=300, **ATDS)[0])]:
        s, tau, lab = score(a, b, p); rng = np.random.default_rng(9); null = []
        for _ in range(100):
            lo = max(p.window, n // 10); null.append(score(a, np.roll(b, int(rng.integers(lo, n - lo))), p)[0])
        st = tau[lab == 1]
        out[name] = dict(window=p.window, max_lag=p.max_lag, tol=p.tolerance, score=round(s, 1), chance=round(float(np.mean(null)), 1),
                         chance_p95=round(float(np.percentile(null, 95)), 1), excess=round(s - float(np.mean(null)), 1),
                         p_emp=round((1 + sum(v >= s for v in null)) / 101, 3), n_windows=len(tau), median_stable_tau=float(np.median(st)) if len(st) else None)
    return out

# ---------- T8: mixed-speed system ----------
def pair_bart(x, y, maxk=None):
    rx, ry = acf(x, maxk), acf(y, maxk); pr = rx * ry; neg = np.nonzero(pr <= 0)[0]; K = int(neg[0]) if len(neg) else len(pr)
    return float(1 + 2 * np.sum(pr[1:K]))
def T8(_):
    n, ntr = N, 40
    def system(seed, coupled):
        rng = np.random.default_rng(seed); cc = 0.4 if coupled else 0.0
        xs, ys, _ = coupled_pair(n, 32, 32, 40, cc, rng); xf, yf, _ = coupled_pair(n, 2, 2, 5, cc, rng); return xs, ys, xf, yf
    tr = system(990000, True)
    sysP, _ = adaptive_params(list(tr), **ATDS)
    def pair_rule(x, y, pairwise):
        if pairwise:
            Bxy = pair_bart(x, y); L = int(np.clip(round(30 * Bxy), 20, 0.1 * n))
            p = TDSParams(window=L, overlap=L // 2, max_lag=(L - 1) // 2)
            return calibrate(p, [(x, y)], 3)
        return adaptive_params([x, y], **ATDS)[0]
    rules = {"system_wide": {"slow": sysP, "fast": sysP}, "per_pair_maxB": {"slow": pair_rule(tr[0], tr[1], False), "fast": pair_rule(tr[2], tr[3], False)},
             "pairwise_Bxy": {"slow": pair_rule(tr[0], tr[1], True), "fast": pair_rule(tr[2], tr[3], True)}, "fixed": {"slow": FIXED, "fast": FIXED}}
    # cross pair (slow x with fast x): never coupled -> check false stability under pairwise rule
    cross = pair_rule(tr[0], tr[2], True)
    out = {"cross_pair_pairwise": dict(window=cross.window, tol=cross.tolerance)}
    for rn, P in rules.items():
        o = {}
        for which, (i, j, d) in {"slow": (0, 1, 40), "fast": (2, 3, 5)}.items():
            p = P[which]; pos, neg, hit = [], [], 0
            for k in range(ntr):
                sp = system(1000000 + k, True); sn = system(2000000 + k, False)
                s, tau, lab = score(sp[i], sp[j], p); pos.append(s); hit += int(found(tau, lab, d)); neg.append(score(sn[i], sn[j], p)[0])
            o[which] = dict(window=p.window, max_lag=p.max_lag, tol=p.tolerance, auc=round(auc(pos, neg), 3), auc_ci=boot_auc(pos, neg), hit=round(hit / ntr, 2), chance=round(float(np.mean(neg)), 2))
        out[rn] = o
    xs = [system(3000000 + k, False) for k in range(20)]
    out["cross_pair_chance_pairwise"] = round(float(np.mean([score(s[0], s[2], cross)[0] for s in xs])), 2)
    return out

# ---------- T9: Figure 1/5 example with parameters from a separate training pair; Fig 3 AUC at the aTDS windows ----------
def T9(_):
    rng = np.random.default_rng(7); x, y, _ = coupled_pair(30000, 8, 8, 20, 0.6, rng)
    ad, info = adaptive_params(list(coupled_pair(30000, 8, 8, 20, 0.6, np.random.default_rng(8))[:2]), **ATDS)
    ex = {}
    for name, p in [("fixed", FIXED), ("atds", ad)]:
        tau, t, _ = time_delay_interaction(x, y, p); lab = stable_label(tau, p); m = t <= 12000
        ex[name] = dict(params=dict(window=p.window, step=p.overlap, max_lag=p.max_lag, tol=p.tolerance), score=round(float(tds_score(lab)), 1),
                        near_true=round(float(np.mean(np.abs(tau + 20) <= 2)) * 100, 1), t=t[m].astype(int).tolist(), tau=tau[m].astype(int).tolist(), stable=lab[m].astype(int).tolist())
    T_ = len(x); r2 = np.random.default_rng(3); lo, hi = max(ad.window, T_ // 10), T_ - max(ad.window, T_ // 10)
    nt = [time_delay_interaction(x, np.roll(y, int(r2.integers(lo, hi))), ad)[0] for _ in range(20)]; rt = time_delay_interaction(x, y, ad)[0]
    curve = [dict(tol=tol, null=round(float(np.mean([tds_score(stable_label(t, replace(ad, tolerance=tol))) for t in nt])), 2),
                  real=round(float(tds_score(stable_label(rt, replace(ad, tolerance=tol)))), 2)) for tol in range(31)]
    # Fig 3: AUC exactly at the aTDS windows
    A3 = {}
    for r, L in ((2, 121), (8, 515)):
        tc, d, c = 2 * r, 5 * r, 0.3; train = [pair(900000 + r * 10 + k, tc, d, c) for k in range(2)]
        p = calibrate(TDSParams(window=L, overlap=L // 2, max_lag=(L - 1) // 2), train, 7)
        pos = [score(*pair(1000 * r + k, tc, d, c), p)[0] for k in range(20)]; neg = [score(*pair(5000 * r + k, tc, d, 0.0), p)[0] for k in range(20)]
        A3[r] = dict(L=L, tol=p.tolerance, auc=round(auc(pos, neg), 3))
    return dict(example=ex, example_info=dict(window=info["window"], tol=info["tolerance"]), tol_curve=curve, fig3_at_atds=A3)

TASKS = {"T1": (T1, [1, 2, 4, 8, 16]), "T2": (T2, [5, 10, 20, 30, 40, 60, 80, 120, 160]),
         "T3": (T3, [(ne, c) for c in (0.3, 0.6) for ne in (5, 10, 15, 20, 30, 45, 60, 90, 120)]),
         "T4": (T4, [0.6, 0.75, 0.9]), "T5": (T5, [8, 30]), "T6": (T6, [16, 32]), "T7": (T7, [0]), "T8": (T8, [0]), "T9": (T9, [0])}
if __name__ == "__main__":
    which = list(TASKS) if sys.argv[4] == "all" else sys.argv[4].split(",")
    try: out = json.load(open(sys.argv[3]))
    except Exception: out = {}
    with Pool(4) as P:
        for t in which:
            t0 = time.time(); f, args = TASKS[t]; r = P.map(f, args)
            out[t] = sum(r, []) if t == "T2" else r
            json.dump(out, open(sys.argv[3], "w")); print(t, "done", round(time.time() - t0), "s", flush=True)
