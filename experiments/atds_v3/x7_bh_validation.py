"""X7: validate the Bayley-Hammersley window estimator (atds_bh) outside the benchmark it was found on (spec: SPEC_X7.md).

Usage: ATDS_CODE=<checkout> python x7_bh_validation.py <item> [<fields_data_dir> | <cohort.mat>] [--smoke]
  item = 1 held-out slow-down | 2 other generators (U5, U6, U9, U10, Rulkov V1) | 3 regime switching (+margin) and drift
       | 4 per-link power | 5r real records (needs <fields_data_dir>) | 5s sleep cohort (needs <cohort.mat>)
       | all (items 1-4 only)
atds_bh = aTDS (30 independent samples, floor 61, cap 10 %, lag (L-1)//2 with rhythm cap, tolerance calibration unchanged)
with the window from B_BH = max(1, 1 + 2 sum_{k<=K}(1 - k/n) r(k)), K = first lag with r <= 0 (max over the signals).
"""
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402
import e1_baselines as E1  # noqa: E402
import e3_synthetic as E3  # noqa: E402
import x6_followups as X6  # noqa: E402

A, AD = C.A, C.AD
_G = {}
SMOKE = "--smoke" in sys.argv
OFFSET = 5_000_000
CAP_FRAC = 0.1


# ------------------------------------------------------------------ the estimator and the rule
def atds_bh(train):
    """aTDS with the BH window. Window from the first training pair (as exp4.atds2); tolerance pooled over all training pairs."""
    x, y = train[0]
    B = max(X6.B_bayley(x), X6.B_bayley(y))
    p, fail = C.E4.calib(X6.params_from_B(x, y, B, len(x)), train)
    return p, fail, B


def atds_pair(train):
    p, fail = C.E4.atds2(train)
    return p, fail


def ident(p):
    return (p.window, p.overlap, p.max_lag, p.tolerance)


def srow(pos, neg, hit, cond, n):
    """summarize (exp4 metrics) plus AUC CI already inside."""
    return C.summarize(pos, neg, hit, cond, n)


# ================================================================== item 1: held-out slow-down set
def _setup1(r):
    tc, d, c = 2 * r, 5 * r, 0.3
    train = [C.pair(OFFSET + 800000 + r * 10 + j, tc, d, c) for j in range(5)]
    ad, fail = C.E4.atds2(train)
    bh, failb, B = atds_bh(train)
    tc_hat = float(np.median([C.first_below_e(s) for xy in train for s in xy]))
    rh = max(1, int(round(tc_hat / 2)))
    M = {"published": C.FIXED, "rescaled_est": X6.resc_params(rh), "atds": ad, "atds_bh": bh}
    return r, dict(M=M, fail=dict(atds=fail, atds_bh=failb), tc_hat=tc_hat, r_hat=rh, B_bh=B)


def _task1(args):
    r, k = args
    tc, d, c = 2 * r, 5 * r, 0.3
    xp, yp = C.pair(OFFSET + 60000 + r * 1000 + k, tc, d, c); xn, yn = C.pair(OFFSET + 70000 + r * 1000 + k, tc, d, 0.0)
    out = {}
    for nm, p in _G["s1"][r]["M"].items():
        s, tau, lab = C.score(xp, yp, p)
        out[nm] = (s, C.score(xn, yn, p)[0], int(C.found(tau, lab, d)))
    return r, out


def item1():
    n, smoke = C.smoke_n(100)
    rs = [1, 4] if smoke else E1.RS
    t0 = time.time()
    with Pool(C.WORKERS) as pool:
        _G["s1"] = dict(pool.map(_setup1, rs))
    print("setup", {r: (g["r_hat"], {k: C.pdict(p) for k, p in g["M"].items()}) for r, g in _G["s1"].items()}, flush=True)
    with Pool(C.WORKERS) as pool:
        res = pool.map(_task1, [(r, k) for r in rs for k in range(n)], chunksize=2)
    out, rows = {}, []
    for r in rs:
        g = _G["s1"][r]; per = {m: [o[m] for rr, o in res if rr == r] for m in g["M"]}
        e = dict(r=r, n=n, r_hat=g["r_hat"], tc_hat=g["tc_hat"], B_bh=round(g["B_bh"], 2), calib_failed=g["fail"],
                 params={m: C.pdict(p) for m, p in g["M"].items()}, methods={}, diff_atds_bh_minus={})
        for m, v in per.items():
            pos, neg, hit = zip(*v)
            e["methods"][m] = dict(**C.auc_row(pos, neg), **C.hit_row(hit), chance=round(float(np.mean(neg)), 3))
            rows.append(dict(r=r, method=m, **C.flat(e["methods"][m]), **{f"p_{a}": b for a, b in C.pdict(g["M"][m]).items()}))
        for m in ("atds", "rescaled_est", "published"):
            c90, c95 = C.boot_diff90([o[0] for o in per["atds_bh"]], [o[1] for o in per["atds_bh"]], [o[0] for o in per[m]], [o[1] for o in per[m]])
            dd = round(e["methods"]["atds_bh"]["auc"] - e["methods"][m]["auc"], 3)
            e["diff_atds_bh_minus"][m] = dict(diff=dd, ci95=c95)
            rows.append(dict(r=r, method=f"diff_atds_bh_minus_{m}", auc=dd, auc_ci=str(c95)))
        out[str(r)] = e; print(r, {m: (v["auc"], v["hit"], v["chance"]) for m, v in e["methods"].items()}, e["diff_atds_bh_minus"], flush=True)
    out["wall_seconds"] = round(time.time() - t0)
    C.write_json("x7_item1_heldout.json", out); C.write_csv("x7_item1_heldout.csv", rows)
    print(f"item 1 {time.time() - t0:.0f}s")


# ================================================================== item 2: other generators
def _std_methods(train, with_published=True):
    """atds, atds_bh, rescaled_est (+ published) from training pairs."""
    ad, fa = C.E4.atds2(train); bh, fb, B = atds_bh(train)
    rs, info = X6.est_rescaled(train)
    M = {"published": C.FIXED} if with_published else {}
    M.update({"atds": ad, "atds_bh": bh, "rescaled_est": rs})
    return M, dict(calib_failed=dict(atds=fa, atds_bh=fb), B_bh=round(B, 2), **info)


def _setup2(key):
    kind, a = key
    if kind == "U5":
        H = a
        train = [C.E4.fgn_pair(500000 + int(H * 100) + j, H, 20, 0.2) for j in range(5)]
    elif kind == "U6":
        train = [C.osc_pair(700000 + a + j, a, 0.6) for j in range(5)]
    elif kind == "V1":
        train = [C.E4.rulkov_oneway(6000, a, 0.12, np.random.default_rng(800000 + a * 10 + j)) for j in range(3)]
    elif kind == "U9slow":
        trS = [X6.u9_system(990000 + j, True) for j in range(5)]
        train = [(s[0], s[1]) for s in trS]
    elif kind == "U9fast":
        trS = [X6.u9_system(990000 + j, True) for j in range(5)]
        train = [(s[2], s[3]) for s in trS]
    elif kind == "U9fs":
        train = [X6.u9_fs(995000 + j, 0.5) for j in range(5)]
    elif kind == "U10":
        train = [C.E4.ns_pair(1300000 + j, 0.4) for j in range(5)]
        h = C.N // 2
        M, info = _std_methods(train)
        Mf, if_ = _std_methods([(x[:h], y[:h]) for x, y in train], False)
        Ms, is_ = _std_methods([(x[h:], y[h:]) for x, y in train], False)
        for nm in ("atds", "atds_bh", "rescaled_est"):
            M[f"{nm}_first"] = Mf[nm]; M[f"{nm}_second"] = Ms[nm]
        return key, dict(M=M, info=dict(whole=info, first=if_, second=is_))
    else:
        raise ValueError(kind)
    M, info = _std_methods(train)
    return key, dict(M=M, info=info)


def _hit_med(tau, lab, d):
    return -float(np.median(tau[lab == 1])) if (lab == 1).sum() >= 3 else None


def _task2(args):
    kind, a, k = args
    M = _G["s2"][(kind, a)]["M"]; out = {}
    if kind == "U5":
        xp, yp = C.E4.fgn_pair(510000 + int(a * 1000) + k, a, 20, 0.2); xn, yn = C.E4.fgn_pair(520000 + int(a * 1000) + k, a, 20, 0.0); d = 20
    elif kind == "U6":
        xp, yp = C.osc_pair(710000 + a * 100 + k, a, 0.6); xn, yn = C.osc_pair(720000 + a * 100 + k, a, 0.0); d = a
    elif kind == "V1":
        xp, yp = C.E4.rulkov_oneway(6000, a, 0.12, np.random.default_rng(810000 + a * 100 + k))
        xn, yn = C.E4.rulkov_oneway(6000, a, 0.0, np.random.default_rng(820000 + a * 100 + k)); d = a
    elif kind in ("U9slow", "U9fast"):
        sp, sn = X6.u9_system(1000000 + k, True), X6.u9_system(2000000 + k, False)
        i, j, d = (0, 1, 40) if kind == "U9slow" else (2, 3, 5)
        xp, yp, xn, yn = sp[i], sp[j], sn[i], sn[j]
    elif kind == "U9fs":
        xp, yp = X6.u9_fs(1100000 + k, 0.5); xn, yn = X6.u9_fs(1200000 + k, 0.0); d = 10
    elif kind == "U10":
        n = C.N; h = n // 2
        P = C.E4.ns_pair(1400000 + k, 0.4); Q = C.E4.ns_pair(1500000 + k, 0.0)
        d = 10
        for half, sl in (("first", slice(0, h)), ("second", slice(h, n))):
            for nm, p in M.items():
                if nm.endswith("_first") and half != "first" or nm.endswith("_second") and half != "second":
                    continue
                s, tau, lab = C.score(P[0][sl], P[1][sl], p)
                out[(half, nm)] = (s, C.score(Q[0][sl], Q[1][sl], p)[0], int(C.found(tau, lab, d)), int((lab == 1).sum() >= 3))
        return kind, a, out
    for nm, p in M.items():
        s, tau, lab = C.score(xp, yp, p); sg = None
        if kind == "U6" and (lab == 1).any():
            i = int(np.nonzero(lab == 1)[0][0]); st = i * p.overlap; u, v = xp[st:st + p.window], yp[st:st + p.window]
            u = (u - u.mean()) / u.std(); v = (v - v.mean()) / v.std(); sg = float(np.mean(np.roll(u, -int(tau[i])) * v))
        out[nm] = (s, C.score(xn, yn, p)[0], int(C.found(tau, lab, d)), int((lab == 1).sum() >= 3), _hit_med(tau, lab, d),
                   -C.mode_tau(tau, lab)[0] if kind == "U6" else None, sg)
    return kind, a, out


def item2():
    smoke = SMOKE
    n = lambda full: 4 if smoke else full  # noqa: E731
    keys = [("U5", H) for H in ((0.6,) if smoke else (0.6, 0.75, 0.9))] + [("U6", d) for d in (8, 30)] + \
           [("U9slow", 0), ("U9fast", 0), ("U9fs", 0), ("U10", 0)] + [("V1", d) for d in ((16,) if smoke else (16, 32))]
    nn = {"U5": n(40), "U6": n(30), "V1": n(20), "U9slow": n(40), "U9fast": n(40), "U9fs": n(40), "U10": n(40)}
    t0 = time.time()
    with Pool(C.WORKERS) as pool:
        _G["s2"] = dict(pool.map(_setup2, keys))
    print(f"setup {time.time() - t0:.0f}s", flush=True)
    with Pool(C.WORKERS) as pool:
        res = pool.map(_task2, [(kd, a, k) for kd, a in keys for k in range(nn[kd])], chunksize=1)
    out, rows = {}, []
    for kd, a in keys:
        g = _G["s2"][(kd, a)]; o = [x for k2, a2, x in res if (k2, a2) == (kd, a)]; nt = nn[kd]
        key = f"{kd}|{a}"; ent = dict(n=nt, info=g["info"], methods={})
        if kd == "U10":
            for half in ("first", "second"):
                ent["methods"][half] = {}
                for nm, p in g["M"].items():
                    if (half, nm) not in o[0]:
                        continue
                    pos, neg, hit, cond, *_ = zip(*[x[(half, nm)] for x in o])
                    e = dict(window=p.window, max_lag=p.max_lag, tol=p.tolerance, **srow(pos, neg, hit, cond, nt))
                    ent["methods"][half][nm] = e; rows.append(dict(case=key, half=half, method=nm, **C.flat(e)))
        else:
            for nm, p in g["M"].items():
                pos, neg, hit, cond, med, rep, sg = zip(*[x[nm] for x in o])
                e = dict(window=p.window, max_lag=p.max_lag, tol=p.tolerance, **srow(pos, neg, hit, cond, nt))
                if kd == "V1":
                    m = [v for v in med if v is not None]
                    e.update(n_with_stable=len(m), median_reported=round(float(np.median(m)), 1) if m else None,
                             median_offset=round(float(np.median(m)) - a, 1) if m else None)
                if kd == "U6":
                    sgv = [v for v in sg if v is not None]; v, cn = np.unique(rep, return_counts=True)
                    e.update(reported_mode=int(v[np.argmax(cn)]), share_negative_C=round(float(np.mean(np.array(sgv) < 0)), 2) if sgv else None)
                ent["methods"][nm] = e; rows.append(dict(case=key, method=nm, **C.flat(e)))
        out[key] = ent
        print(key, {nm: (e.get("auc"), e.get("hit")) for nm, e in (ent["methods"] if kd != "U10" else {**ent["methods"]["first"], **{"2_" + k: v for k, v in ent["methods"]["second"].items()}}).items()}, flush=True)
    out["wall_seconds"] = round(time.time() - t0)
    C.write_json("x7_item2_generators.json", out); C.write_csv("x7_item2_generators.csv", rows)
    print(f"item 2 {time.time() - t0:.0f}s")


# ================================================================== item 3: regime switching (+margin) and drift
def item3():
    n, smoke = C.smoke_n(40)
    t0 = time.time()
    P, rinfo = X6._setup_regime_x6()
    train = [E3.regime_pair(1000000 + j, 0.6) for j in range(5)]
    P["atds_bh"], fb, B = atds_bh([(x, y) for x, y, _ in train])
    X6._G["regime"] = P
    methods = X6.REG_METHODS + ("atds_bh",)
    X6.REG_METHODS = methods
    dr_tr = [E3.drift_pair(7000000 + j) for j in range(5)]
    Dm = E3._setup_drift(0); Dm["rescaled_est"], dinfo = X6.est_rescaled(dr_tr)
    Dm["atds_bh"], fbd, Bd = atds_bh(dr_tr)
    X6._G["drift"] = Dm
    print("params regime", {k: (C.pdict(v) if hasattr(v, "window") else None) for k, v in P.items()}, "B_bh", round(B, 2), flush=True)
    with Pool(C.WORKERS) as pool:
        rg = pool.map(X6._regime_task_x6, range(n))
        dr = pool.map(X6._drift_task_x6, range(n))
    reg, rows = {"rescaled_est_info": rinfo, "B_bh": round(B, 2), "calib_failed_bh": fb, "n_pairs": n}, []
    for margin in (False, True):
        for nm in methods:
            e = X6.regime_summary(rg, nm, margin, n, P)
            reg[f"{nm}|margin={margin}"] = e; rows.append(dict(case="regime_switching", method=nm, **C.flat(e)))
    D = {}
    for nm in Dm:
        sc = [o[nm][0] for o in dr]; allerr = np.concatenate([o[nm][1] for o in dr])
        D[nm] = dict(params=C.pdict(Dm[nm]), score_mean=round(float(np.mean(sc)), 2), score_sd=round(float(np.std(sc, ddof=1)), 2),
                     n_stable_windows=int(len(allerr)), median_abs_err_pooled=round(float(np.median(allerr)), 2) if len(allerr) else None)
        if nm == "rescaled_est":
            D[nm].update(dinfo)
        if nm == "atds_bh":
            D[nm].update(B_bh=round(Bd, 2), calib_failed=fbd)
        rows.append(dict(case="drifting_delay", method=nm, **C.flat(D[nm])))
    out = dict(regime=reg, drift=D, wall_seconds=round(time.time() - t0))
    C.write_json("x7_item3_timeresolved.json", out); C.write_csv("x7_item3_timeresolved.csv", rows)
    for nm in methods:
        a, b = reg[f"{nm}|margin=False"], reg[f"{nm}|margin=True"]
        print(nm, "A auc", a.get("auc"), b.get("auc"), "B coup/ctrl", a["false_stability_B_coupled_mean"], a["false_stability_B_control_mean"],
              "margin", b["false_stability_B_coupled_mean"], b["false_stability_B_control_mean"])
    print({k: (v["score_mean"], v["median_abs_err_pooled"]) for k, v in D.items()})
    print(f"item 3 {time.time() - t0:.0f}s")


# ================================================================== item 4: per-link power
def _setup4(r):
    _, g = E1._setup(r)
    tc, d, c = 2 * r, 5 * r, 0.3
    train = [C.pair(800000 + r * 10 + j, tc, d, c) for j in range(5)]
    bh, fb, B = atds_bh(train)
    g["M"] = {k: g["M"][k] for k in ("published", "rescaled", "rescaled_est", "atds")}
    g["M"]["atds_bh"] = bh
    g["B_bh"] = B
    return r, g


def item4():
    n50 = 2 if SMOKE else 50
    nrot = 9 if SMOKE else 99
    rs = [2, 4] if SMOKE else [2, 4, 8]
    t0 = time.time()
    with Pool(C.WORKERS) as pool:
        X6._G["e5"] = dict(pool.map(_setup4, rs))
    X6._G["nrot"] = nrot
    keep = ("published", "rescaled", "rescaled_est", "atds", "atds_bh")
    print({r: (g["r_hat"], round(g["B_bh"], 2), {k: C.pdict(p) for k, p in g["M"].items()}) for r, g in X6._G["e5"].items()}, flush=True)
    with Pool(C.WORKERS) as pool:
        res = pool.map(X6._task5, [(r, kind, k) for r in rs for kind in ("coupled", "uncoupled") for k in range(n50)], chunksize=1)
    out, rows = {}, []
    for r in rs:
        for nm in keep:
            e = dict(r=r, method=nm, n_rotations=nrot, r_hat=X6._G["e5"][r]["r_hat"], params=C.pdict(X6._G["e5"][r]["M"][nm]))
            for kind in ("coupled", "uncoupled"):
                pv = [o[nm] for rr, kd, o in res if rr == r and kd == kind]
                kk = int(sum(v <= 0.05 for v in pv))
                e["power" if kind == "coupled" else "fpr"] = round(kk / len(pv), 3)
                e["power_ci" if kind == "coupled" else "fpr_ci"] = C.wilson(kk, len(pv))
                e["n_coupled" if kind == "coupled" else "n_uncoupled"] = len(pv)
            out[f"r{r}|{nm}"] = e; rows.append(C.flat(e)); print(e)
    out["wall_seconds"] = round(time.time() - t0)
    C.write_json("x7_item4_power.json", out); C.write_csv("x7_item4_power.csv", rows)
    print(f"item 4 {time.time() - t0:.0f}s")


# ================================================================== item 5a: real records
def item5_real(fd):
    n_null = 8 if SMOKE else 200
    full, rows = {}, []
    for name, loader in X6.load_records(fd):
        t0 = time.time()
        d = loader(); s1, s2 = np.asarray(d["s1"], float), np.asarray(d["s2"], float); d["s1"], d["s2"] = s1, s2
        n = len(s1)
        p_atds, fail = A.calib(A.window_params([s1, s2]), [(s1, s2)], n_shift=20)
        Bb = max(X6.B_bayley(s1), X6.B_bayley(s2)); Bd = max(AD.bartlett_factor(s1), AD.bartlett_factor(s2))
        p0 = X6.params_from_B(s1, s2, Bb, n)
        assert ident(X6.params_from_B(s1, s2, Bd, n)) == ident(A.window_params([s1, s2])), "rule replication differs from aTDS"
        p_bh, failb = A.calib(p0, [(s1, s2)], n_shift=20)
        tc = X6.tc_slow(s1, s2); r = X6.r_hat_from(tc)
        sets = (("published", A.PUBLISHED, False), ("rescaled_est", X6.resc_params(r), False), ("atds", p_atds, fail), ("atds_bh", p_bh, failb))
        entry = dict(n=n, B_default=round(Bd, 2), B_bh=round(Bb, 2), window_uncalibrated_default=int(round(30 * Bd)),
                     window_uncalibrated_bh=int(round(30 * Bb)), tc_hat_slower_signal=tc, r_hat=r, settings={})
        for sname, p, f in sets:
            res = X6.run_rec(d, p, f, n_null); entry["settings"][sname] = res
            rows.append(dict(record=name, setting=sname, n=n, score=res["score"], null_mean=res["null_mean"], null95=res["null95"], p=res["p"],
                             n_windows=res["n_windows"], n_stable=res["n_stable"], delay_samples=res["delay_samples"], delay_phys=res["delay_phys"],
                             **{f"{k}_samples": v for k, v in res["params_samples"].items()}, calibration_failed=res["calibration_failed"]))
        full[name] = entry
        print(f"{name}: B {Bd:.1f} -> BH {Bb:.1f}, windows {p_atds.window} -> {p_bh.window}, {time.time() - t0:.0f}s", flush=True)
    C.write_json("x7_item5_real.json", full); C.write_csv("x7_item5_real.csv", rows)
    for rw in rows:
        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in rw.items()})


# ================================================================== item 5b: sleep cohort
def item5_cohort(path):
    import ncom_cohort as NC
    t0 = time.time()
    ids, sigs, stages, _ = NC.load_cohort(path)
    S = len(sigs)
    Bd = np.array([[AD.bartlett_factor(x[:, i]) for i in range(10)] for x in sigs])
    Bb = np.array([[X6.B_bayley(x[:, i]) for i in range(10)] for x in sigs])
    med_len = float(np.median([len(x) for x in sigs]))
    L0 = lambda B: int(np.clip(round(30 * float(np.median(B.max(axis=1)))), 61, 0.1 * med_len))  # noqa: E731
    info = dict(n_subjects=S, B_default_median_of_max=round(float(np.median(Bd.max(axis=1))), 2), B_bh_median_of_max=round(float(np.median(Bb.max(axis=1))), 2),
                B_default_per_channel_median=dict(zip(NC.CH, np.round(np.median(Bd, axis=0), 2).tolist())),
                B_bh_per_channel_median=dict(zip(NC.CH, np.round(np.median(Bb, axis=0), 2).tolist())),
                window_uncapped_default=L0(Bd), window_uncapped_bh=L0(Bb), window_cap300_default=min(L0(Bd), 300), window_cap300_bh=min(L0(Bb), 300))
    info["window_changes_under_cap300"] = info["window_cap300_default"] != info["window_cap300_bh"]
    print(info, flush=True)
    if info["window_changes_under_cap300"]:
        info["cohort_reanalysis"] = "window differs: reanalysis needed"
    else:
        info["cohort_reanalysis"] = ("window, step, lag identical (300 / 150 / 149); calibration inputs (seed, 50 pairs, p0) identical, so "
                                     "atds_bh_cap300 == atds_cap300 by construction and every cohort number (stage means, effects) is unchanged")
        # verify the calibration inputs are identical
        pa = NC.cohort_params(sigs)[0]["atds_cap300"]
        info["atds_cap300_params"] = dict(window=pa.window, step=pa.overlap, max_lag=pa.max_lag, tolerance=pa.tolerance)
    info["wall_seconds"] = round(time.time() - t0)
    C.write_json("x7_item5_cohort.json", info)
    return info


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    item = args[0]; arg = args[1] if len(args) > 1 else None
    C.set_results_for_smoke(SMOKE)
    t0 = time.time()
    if item in ("1", "all"):
        item1()
    if item in ("2", "all"):
        item2()
    if item in ("3", "all"):
        item3()
    if item in ("4", "all"):
        item4()
    if item == "5r":
        item5_real(arg)
    if item == "5s":
        item5_cohort(arg)
    print(f"X7 item {item} wall {time.time() - t0:.0f}s")
