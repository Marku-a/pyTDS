"""X6: follow-up experiments requested by the run-3 reviewers (spec: SPEC_X6.md).

Usage: ATDS_CODE=<checkout> python x6_followups.py <item> [<fields_data_dir>] [--smoke]
  item = 1 | 2 | 3 | 4 | 5 | all      (items 1 and 3 need <fields_data_dir> with the raw/*.rda and CSVs of E5 / Part F)
Item 1 parts: regime switching + drifting delay (E3), U9 and U6 (exp4), real records (200 rotations).
Item 2 recomputes the regime-resolved scores with and without a one-window boundary margin.
Item 3 diagnoses L0123002 and tabulates stable-window counts (reads the item-1 real-record results).
Items 4 and 5: slow-down benchmark (window rule variants; per-link rotation tests).
"""
import json
import os
import sys
import time
from dataclasses import replace
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402
import e1_baselines as E1  # noqa: E402
import e3_synthetic as E3  # noqa: E402
import e5_hydrology as E5  # noqa: E402
import fields as F  # noqa: E402

A = C.A
AD = C.AD
_G = {}
SMOKE = "--smoke" in sys.argv
BLK = 2000


# ------------------------------------------------------------------ shared helpers
def tc_slow(x, y):
    """Decorrelation time (first lag below 1/e) of the slower of two signals."""
    return max(C.first_below_e(x), C.first_below_e(y))


def r_hat_from(tc_hat):
    return max(1, int(round(tc_hat / 2)))


def resc_params(r):
    """Published setting (60 / 30 / 30 / 1) multiplied by r."""
    return C.TDSParams(window=60 * r, overlap=30 * r, max_lag=30 * r, tolerance=r)


def est_rescaled(train_pairs):
    """rescaled_est from training pairs: tc_hat = median over pairs of the slower signal's 1/e lag."""
    tc = float(np.median([tc_slow(x, y) for x, y in train_pairs]))
    r = r_hat_from(tc)
    return resc_params(r), dict(tc_hat=tc, r_hat=r)


def fnum(v, nd=3):
    return None if v is None or (isinstance(v, float) and not np.isfinite(v)) else round(float(v), nd)


# ================================================================== item 1a / 2: regime switching
def regime_eval(x, y, Amask, p, margin, d=20):
    tau = C.time_delay_interaction(x, y, p)[0]
    lab = C.stable_label(tau, p).astype(bool)
    W = p.window
    st = np.arange(len(tau)) * p.overlap
    if margin:   # entirely inside one block and >= one window length from both block edges
        b = st // BLK
        ok = ((st + W - 1) // BLK == b) & (st - b * BLK >= W) & (BLK * (b + 1) - (st + W) >= W)
        inA, inB = ok & (b % 2 == 0), ok & (b % 2 == 1)
        need = 1
    else:        # 80 % purity rule of E3
        cs = np.concatenate([[0], np.cumsum(Amask)])
        fa = (cs[st + W] - cs[st]) / W
        inA, inB = fa >= 0.8, fa <= 0.2
        need = 3
    out = {}
    for nm, m in (("A", inA), ("B", inB)):
        out[nm] = 100.0 * lab[m].mean() if m.sum() >= need else np.nan
        out["n" + nm] = int(m.sum())
    sa = lab & inA
    out["found"] = int(sa.sum() >= 3 and abs(np.median(tau[sa]) + d) <= max(1, 0.1 * d))
    return out


def _setup_regime_x6():
    P = E3._setup_regime(0)
    train = [E3.regime_pair(1000000 + j, 0.6) for j in range(5)]
    P["rescaled_est"], info = est_rescaled([(x, y) for x, y, _ in train])
    return P, info


REG_METHODS = ("published", "atds", "atds_cap1000", "atds_per_regime", "rescaled_est")


def _regime_task_x6(k):
    P = _G["regime"]
    xc, yc, Am = E3.regime_pair(2000000 + k, 0.6)
    xn, yn, _ = E3.regime_pair(3000000 + k, 0.0)
    out = {}
    for margin in (False, True):
        for nm in REG_METHODS:
            if nm == "atds_per_regime":
                res = []
                for x, y in ((xc, yc), (xn, yn)):
                    a = regime_eval(x, y, Am, P["regime_A"], margin); b = regime_eval(x, y, Am, P["regime_B"], margin)
                    res.append(dict(A=a["A"], B=b["B"], nA=a["nA"], nB=b["nB"], found=a["found"]))
                out[(nm, margin)] = tuple(res)
            else:
                out[(nm, margin)] = (regime_eval(xc, yc, Am, P[nm], margin), regime_eval(xn, yn, Am, P[nm], margin))
    return out


def _nm(v):
    v = np.asarray(v, float)
    return float(np.nanmean(v)) if np.isfinite(v).any() else None


def regime_summary(res, nm, margin, n, P):
    cA = [r[(nm, margin)][0]["A"] for r in res]; nA = [r[(nm, margin)][1]["A"] for r in res]
    cB = [r[(nm, margin)][0]["B"] for r in res]; nB = [r[(nm, margin)][1]["B"] for r in res]
    ok = [i for i in range(n) if np.isfinite(cA[i]) and np.isfinite(nA[i])]
    pr = ({"A": C.pdict(P["regime_A"]), "B": C.pdict(P["regime_B"])} if nm == "atds_per_regime" else C.pdict(P[nm]))
    e = dict(params=pr, margin=margin, n_pairs_with_A_windows=len(ok),
             mean_windows_A=fnum(np.mean([r[(nm, margin)][0]["nA"] for r in res]), 1),
             mean_windows_B=fnum(np.mean([r[(nm, margin)][0]["nB"] for r in res]), 1))
    if ok:
        e.update(C.auc_row([cA[i] for i in ok], [nA[i] for i in ok]))
    else:
        e.update(auc=None, auc_ci=None)
    e.update(regimeA_coupled_mean=fnum(_nm(cA), 2), regimeA_control_mean=fnum(_nm(nA), 2),
             false_stability_B_coupled_mean=fnum(_nm(cB), 2), false_stability_B_control_mean=fnum(_nm(nB), 2),
             n_pairs_with_B_windows=int(np.isfinite(cB).sum()),
             **C.hit_row([r[(nm, margin)][0]["found"] for r in res]))
    if not ok:
        e["note"] = "no window inside a regime-A block after the margin: not reportable"
    return e


# ================================================================== item 1b: drifting delay
def _drift_task_x6(k):
    x, y = E3.drift_pair(8000000 + k)
    out = {}
    for nm, p in _G["drift"].items():
        tau, t, _ = C.time_delay_interaction(x, y, p); lab = C.stable_label(tau, p).astype(bool)
        true = 10 + 20 * t / E3.NS
        out[nm] = (float(C.tds_score(lab)), np.abs(tau[lab] + true[lab]))
    return out


# ================================================================== item 1c: exp4 U9
def u9_system(seed, coupled, n=C.N):
    rng = np.random.default_rng(seed); cc = 0.4 if coupled else 0.0
    xs, ys, _ = C.coupled_pair(n, 32, 32, 40, cc, rng); xf, yf, _ = C.coupled_pair(n, 2, 2, 5, cc, rng)
    return xs, ys, xf, yf


def u9_fs(seed, c, n=C.N):
    rng = np.random.default_rng(seed); x = C.ar1(n, 32, rng); e = C.ar1(n, 2, rng)
    return x, c * np.roll(x, 10) + np.sqrt(1 - c * c) * e


def _setup_u9():
    E = C.E4
    trS = [u9_system(990000 + j, True) for j in range(5)]
    sysw, _ = AD.adaptive_params(list(trS[0]), window_rule="bartlett", tolerance="fixed", min_window=C.MINW)
    sysw, _ = E.calib(sysw, [(s[0], s[1]) for s in trS] + [(s[2], s[3]) for s in trS])
    slow_tr = [(s[0], s[1]) for s in trS]; fast_tr = [(s[2], s[3]) for s in trS]
    rs, ris = est_rescaled(slow_tr); rf, rif = est_rescaled(fast_tr)
    tc_sys = float(np.median([max(tc_slow(s[0], s[1]), tc_slow(s[2], s[3])) for s in trS]))
    rsys = resc_params(r_hat_from(tc_sys))
    rules = {"system_wide": {"slow": sysw, "fast": sysw},
             "per_pair_maxB": {"slow": E.atds2(slow_tr)[0], "fast": E.atds2(fast_tr)[0]},
             "pairwise_Bxy": {"slow": E.atds2(slow_tr, pairwise=True)[0], "fast": E.atds2(fast_tr, pairwise=True)[0]},
             "fixed": {"slow": C.FIXED, "fast": C.FIXED},
             "rescaled_est": {"slow": rs, "fast": rf},
             "rescaled_est_system": {"slow": rsys, "fast": rsys}}
    info = dict(rescaled_est_slow=ris, rescaled_est_fast=rif, rescaled_est_system=dict(tc_hat=tc_sys, r_hat=r_hat_from(tc_sys)))
    trF = [u9_fs(995000 + j, 0.5) for j in range(5)]
    rF, riF = est_rescaled(trF)
    M = {"per_pair_maxB": E.atds2(trF)[0], "pairwise_Bxy": E.atds2(trF, pairwise=True)[0], "fixed": C.FIXED, "rescaled_est": rF}
    info["rescaled_est_fast_slow"] = riF
    return rules, M, info


def _u9_task(k):
    rules, M = _G["u9"]
    sp, sn = u9_system(1000000 + k, True), u9_system(2000000 + k, False)
    out = {}
    for rn, P in rules.items():
        for which, (i, j, d) in {"slow": (0, 1, 40), "fast": (2, 3, 5)}.items():
            p = P[which]
            s, tau, lab = C.score(sp[i], sp[j], p)
            out[(rn, which)] = (s, C.score(sn[i], sn[j], p)[0], int(C.found(tau, lab, d)), int((lab == 1).sum() >= 3))
    xp, yp = u9_fs(1100000 + k, 0.5); xn, yn = u9_fs(1200000 + k, 0.0)
    for nm, p in M.items():
        s, tau, lab = C.score(xp, yp, p)
        out[("fast_slow_pair", nm)] = (s, C.score(xn, yn, p)[0], int(C.found(tau, lab, 10)), int((lab == 1).sum() >= 3))
    return out


# ================================================================== item 1d: exp4 U6
def _setup_u6(d):
    train = [C.osc_pair(700000 + d + j, d, 0.6) for j in range(5)]
    ad, fail = C.E4.atds2(train)
    resc, info = est_rescaled(train)
    return d, dict(M={"fixed": C.FIXED, "atds": ad, "rescaled_est": resc}, fail=fail, info=info)


def _u6_task(args):
    d, k = args
    x, y = C.osc_pair(710000 + d * 100 + k, d, 0.6); xn, yn = C.osc_pair(720000 + d * 100 + k, d, 0.0)
    out = {}
    for nm, p in _G["u6"][d]["M"].items():
        s, tau, lab = C.score(x, y, p); sg = None
        if (lab == 1).any():
            i = int(np.nonzero(lab == 1)[0][0]); st = i * p.overlap; a, b = x[st:st + p.window], y[st:st + p.window]
            a = (a - a.mean()) / a.std(); b = (b - b.mean()) / b.std(); sg = float(np.mean(np.roll(a, -int(tau[i])) * b))
        out[nm] = (s, C.score(xn, yn, p)[0], int(C.found(tau, lab, d)), -C.mode_tau(tau, lab)[0], sg)
    return d, out


# ================================================================== item 1e: real records
def run_rec(d, p, fail, n_null):
    """Same output as fields.run_setting but with a configurable number of rotations."""
    s1, s2, dt = d["s1"], d["s2"], d["dt"]
    tau = C.time_delay_interaction(s1, s2, p)[0]
    lab = C.stable_label(tau, p)
    score = float(C.tds_score(lab))
    null = A.rotation_null(s1, s2, p, n_null, seed=0)
    dl = A.delay_of(tau, lab)
    return {"score": score, "n_windows": int(len(tau)), "n_stable": int(lab.sum()), **A.null_stats(score, null),
            "n_null": n_null, "delay_samples": dl, "delay_phys": None if dl is None else F.phys(dl, dt),
            "delay_after_driver_samples": None if dl is None else -dl,
            "params_samples": {"window": p.window, "step": p.overlap, "max_lag": p.max_lag, "tolerance": p.tolerance},
            "params_phys": {"window": F.phys(p.window, dt), "step": F.phys(p.overlap, dt), "max_lag": F.phys(p.max_lag, dt),
                            "tolerance": F.phys(p.tolerance, dt) if p.tolerance else "0"},
            "calibration_failed": bool(fail)}


def load_records(fd):
    F.FD = fd
    for code in ("L0123001", "L0123002"):
        E5.build_daily(fd, code)
    return [("L0123003_hourly", F.d_hydro_hourly), ("L0123001_daily", E5.load_daily("L0123001")),
            ("L0123002_daily", E5.load_daily("L0123002")), ("energy_temp_electricity", F.d_energy),
            ("gasoline_weekly_diff", F.d_gasoline)]


def part_real(fd):
    n_null = 8 if SMOKE else 200
    full, rows = {}, []
    for name, loader in load_records(fd):
        t0 = time.time()
        d = loader(); s1, s2 = np.asarray(d["s1"], float), np.asarray(d["s2"], float); d["s1"], d["s2"] = s1, s2
        p_atds, fail = A.calib(A.window_params([s1, s2]), [(s1, s2)], n_shift=20)
        tc = tc_slow(s1, s2); r = r_hat_from(tc)
        sets = (("published", A.PUBLISHED, False), ("atds", p_atds, fail), ("rescaled_est", resc_params(r), False))
        entry = dict(n=len(s1), sampling_interval_s=d["dt"], tc_hat_slower_signal=tc, r_hat=r, settings={})
        for sname, p, f in sets:
            res = run_rec(d, p, f, n_null)
            entry["settings"][sname] = res
            rows.append(dict(record=name, setting=sname, n=len(s1), score=res["score"], null_mean=res["null_mean"], null95=res["null95"],
                             p=res["p"], n_windows=res["n_windows"], n_stable=res["n_stable"], delay_samples=res["delay_samples"],
                             delay_phys=res["delay_phys"], r_hat=r if sname == "rescaled_est" else None,
                             **{f"{k}_samples": v for k, v in res["params_samples"].items()}, calibration_failed=res["calibration_failed"]))
        full[name] = entry
        print(f"{name}: n={len(s1)} r_hat={r} done in {time.time() - t0:.0f}s", flush=True)
    C.write_json("x6_item1_real.json", full); C.write_csv("x6_item1_real.csv", rows)
    for rw in rows:
        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in rw.items()})
    return full


# ================================================================== item 1 / 2 driver
def item1_item2(fd):
    n, smoke = C.smoke_n(40)
    t0 = time.time(); res1 = {}; rows1 = []
    # ---- regime switching (items 1 and 2)
    _G["regime"], rinfo = _setup_regime_x6()
    _G["drift"] = E3._setup_drift(0)
    train = [E3.drift_pair(7000000 + j) for j in range(5)]
    _G["drift"]["rescaled_est"], dinfo = est_rescaled(train)
    with Pool(C.WORKERS) as pool:
        rg = pool.map(_regime_task_x6, range(n))
        dr = pool.map(_drift_task_x6, range(n))
    P = _G["regime"]
    reg = {"rescaled_est_info": rinfo, "n_pairs": n}; rows2 = []; item1_reg = {}
    for margin in (False, True):
        for nm in REG_METHODS:
            e = regime_summary(rg, nm, margin, n, P)
            reg[f"{nm}|margin={margin}"] = e
            rows2.append(dict(method=nm, **C.flat(e)))
            if not margin:
                item1_reg[nm] = e
    C.write_json("x6_item2_margin.json", dict(blocks=BLK, rule="windows entirely inside a block and >= one window length from both block edges "
                                                          "(>=1 window per regime required); without margin: 80 % purity rule (>=3 windows)",
                                              results=reg, wall_seconds=round(time.time() - t0)))
    C.write_csv("x6_item2_margin.csv", rows2)
    res1["regime_switching"] = item1_reg
    rows1 += [dict(case="regime_switching", method=nm, **C.flat(e)) for nm, e in item1_reg.items()]
    print(f"regime {time.time() - t0:.0f}s")
    # ---- drifting delay
    E = {}
    for nm in _G["drift"]:
        sc = [o[nm][0] for o in dr]; err = [o[nm][1] for o in dr]; allerr = np.concatenate(err)
        E[nm] = dict(params=C.pdict(_G["drift"][nm]), score_mean=round(float(np.mean(sc)), 2), score_sd=round(float(np.std(sc, ddof=1)), 2),
                     n_stable_windows=int(len(allerr)),
                     median_abs_err_pooled=round(float(np.median(allerr)), 2) if len(allerr) else None,
                     median_of_pair_medians=round(float(np.median([np.median(e) for e in err if len(e)])), 2) if any(len(e) for e in err) else None)
        if nm == "rescaled_est":
            E[nm].update(dinfo)
        rows1.append(dict(case="drifting_delay", method=nm, **C.flat(E[nm])))
    res1["drifting_delay"] = E
    print(f"drift {time.time() - t0:.0f}s")
    # ---- U9
    rules, M, uinfo = _setup_u9(); _G["u9"] = (rules, M)
    nt = 4 if smoke else 40
    with Pool(C.WORKERS) as pool:
        u9r = pool.map(_u9_task, range(nt))
    U9 = dict(info=uinfo)
    for key in u9r[0]:
        rn, which = key
        pos, neg, hit, cond = zip(*[o[key] for o in u9r])
        p = (rules[rn][which] if which in ("slow", "fast") else M[which])
        e = dict(window=p.window, max_lag=p.max_lag, tol=p.tolerance, **C.summarize(pos, neg, hit, cond, nt))
        if which in ("slow", "fast"):
            U9.setdefault(rn, {})[which] = e
        else:
            U9.setdefault("fast_slow_pair", {})[which] = e
        rows1.append(dict(case=f"U9_{rn}" if which in ("slow", "fast") else "U9_fast_slow_pair", method=which, **C.flat(e)))
    res1["U9"] = U9
    print(f"U9 {time.time() - t0:.0f}s")
    # ---- U6
    nr = 4 if smoke else 30
    with Pool(C.WORKERS) as pool:
        _G["u6"] = dict(pool.map(_setup_u6, [8, 30]))
    with Pool(C.WORKERS) as pool:
        u6r = pool.map(_u6_task, [(d, k) for d in (8, 30) for k in range(nr)])
    U6 = {}
    for d in (8, 30):
        o = [x for dd, x in u6r if dd == d]; g = _G["u6"][d]
        U6[f"d{d}"] = dict(rescaled_est_info=g["info"], calib_failed=g["fail"])
        for nm, p in g["M"].items():
            pos, neg, hit, rep, sg = zip(*[x[nm] for x in o]); sg = [v for v in sg if v is not None]
            v, cn = np.unique(rep, return_counts=True)
            e = dict(auc=round(C.auc(pos, neg), 3), auc_ci=C.boot_auc(pos, neg), hit=round(sum(hit) / nr, 2), chance=round(float(np.mean(neg)), 2),
                     reported_mode=int(v[np.argmax(cn)]), share_negative_C=round(float(np.mean(np.array(sg) < 0)), 2) if sg else None,
                     window=p.window, max_lag=p.max_lag, tol=p.tolerance)
            U6[f"d{d}"][nm] = e
            rows1.append(dict(case=f"U6_d{d}", method=nm, **C.flat(e)))
    res1["U6"] = U6
    print(f"U6 {time.time() - t0:.0f}s")
    res1["wall_seconds_synthetic"] = round(time.time() - t0)
    C.write_json("x6_item1_synthetic.json", res1); C.write_csv("x6_item1_synthetic.csv", rows1)
    # consistency check of the replicated U9 / U6 against the stored exp4 results
    try:
        ref = json.load(open(os.path.join(C.ROOT, "atds_brief", "results", "exp4.json")))
        chk = {}
        if not smoke:
            r9 = ref["U9"][0]
            for rn in ("system_wide", "per_pair_maxB", "pairwise_Bxy", "fixed"):
                chk[f"U9_{rn}"] = {w: dict(now=U9[rn][w]["auc"], stored=r9[rn][w]["auc"]) for w in ("slow", "fast")}
            for d in (8, 30):
                rr = [x for x in ref["U6"] if x["d"] == d][0]
                chk[f"U6_d{d}"] = {nm: dict(now=U6[f"d{d}"][nm]["auc"], stored=rr[nm2]["auc"]) for nm, nm2 in (("fixed", "fixed"), ("atds", "atds"))}
        res1["check_vs_exp4"] = chk
        C.write_json("x6_item1_synthetic.json", res1)
        print("check vs stored exp4:", chk)
    except Exception as ex:  # noqa: BLE001
        print("check failed:", ex)
    # ---- real records
    if fd:
        part_real(fd)
    return res1


# ================================================================== item 3
def item3(fd):
    n_null = 8 if SMOKE else 200
    t0 = time.time()
    loaders = dict(load_records(fd))
    d = loaders["L0123002_daily"](); s1, s2 = np.asarray(d["s1"], float), np.asarray(d["s2"], float); d["s1"], d["s2"] = s1, s2
    n = len(s1)
    Bx, By = AD.bartlett_factor(s1), AD.bartlett_factor(s2)
    out = dict(n=n, bartlett_precip=round(Bx, 2), bartlett_flow=round(By, 2), tc_1e_precip=C.first_below_e(s1), tc_1e_flow=C.first_below_e(s2),
               window_before_cap=int(round(30 * max(Bx, By))), cap_10pct=max(61, int(0.1 * n)), runs={})
    out["acf_flow"] = {str(k): round(float(v), 3) for k, v in zip((1, 10, 30, 100, 365), AD.acf(s2, 400)[[1, 10, 30, 100, 365]])}
    for label, cap in (("atds_default", None), ("atds_cap365", 365), ("atds_cap120", 120)):
        p, fail = A.calib(A.window_params([s1, s2], cap=cap), [(s1, s2)], n_shift=20)
        res = run_rec(d, p, fail, n_null)
        out["runs"][label] = dict(window=p.window, step=p.overlap, max_lag=p.max_lag, tolerance=p.tolerance, n_windows=res["n_windows"],
                                  n_stable=res["n_stable"], score=res["score"], null_mean=res["null_mean"], null95=res["null95"], p=res["p"],
                                  delay_samples=res["delay_samples"], calibration_failed=fail)
        print(label, out["runs"][label], flush=True)
    res = run_rec(d, A.PUBLISHED, False, n_null)
    out["runs"]["published"] = dict(window=60, step=30, max_lag=30, tolerance=1, n_windows=res["n_windows"], n_stable=res["n_stable"], score=res["score"],
                                    null_mean=res["null_mean"], null95=res["null95"], p=res["p"], delay_samples=res["delay_samples"])
    # stable windows behind the delay estimate, every real record / setting (from the item-1 real run)
    rp = os.path.join(C.RESULTS, "x6_item1_real.json")
    counts = []
    if os.path.exists(rp):
        for rec, e in json.load(open(rp)).items():
            for sname, r in e["settings"].items():
                counts.append(dict(record=rec, setting=sname, n_windows=r["n_windows"], n_stable=r["n_stable"], delay_samples=r["delay_samples"],
                                   delay_estimated=r["delay_samples"] is not None, score=r["score"], p=r["p"]))
    out["stable_window_counts"] = counts
    out["wall_seconds"] = round(time.time() - t0)
    C.write_json("x6_item3_L0123002.json", out)
    C.write_csv("x6_item3_stable_windows.csv", counts or [dict(note="run item 1 first")])
    print(json.dumps(out, indent=1, default=str)); print(f"total {time.time() - t0:.0f}s")


# ================================================================== item 4: n_eff estimators
def B_default(x):
    return AD.bartlett_factor(x)


def B_all_lags(x):
    r = AD.acf(x); return float(1.0 + 2.0 * np.sum(r[1:] ** 2))      # lags up to T/4, negative lobes add


def B_bayley(x):
    r = AD.acf(x); n = len(x); z = np.nonzero(r <= 0)[0]; K = int(z[0]) if len(z) else len(r) - 1
    k = np.arange(1, K + 1)
    return float(max(1.0, 1.0 + 2.0 * np.sum((1 - k / n) * r[1:K + 1])))


ESTS = {"a_bartlett_first_zero": B_default, "b_bartlett_all_lags": B_all_lags, "c_bayley_hammersley": B_bayley}


def params_from_B(x, y, B, n):
    L = int(np.clip(round(30 * B), 61, max(61, int(0.1 * n))))
    ml = (L - 1) // 2
    det = [q for q in (AD.dominant_period(s) for s in (x, y)) if q is not None]
    if det:
        ml = max(1, min(ml, min(det) // 2))
    return C.TDSParams(window=L, overlap=max(1, L // 2), max_lag=ml)


def _setup4(args):
    r, en = args
    tc, d, c = 2 * r, 5 * r, 0.3
    train = [C.pair(800000 + r * 10 + j, tc, d, c) for j in range(5)]
    x, y = train[0]
    B = max(ESTS[en](x), ESTS[en](y))
    p0 = params_from_B(x, y, B, len(x))
    p, fail = C.E4.calib(p0, train)
    if en.startswith("a_"):   # must equal the stored aTDS rule
        pa, _ = C.E4.atds2(train)
        assert (pa.window, pa.overlap, pa.max_lag, pa.tolerance) == (p.window, p.overlap, p.max_lag, p.tolerance), (pa, p)
    return (r, en), dict(p=p, fail=fail, B=B)


def _task4(args):
    r, k = args
    tc, d, c = 2 * r, 5 * r, 0.3
    xp, yp = C.pair(60000 + r * 1000 + k, tc, d, c); xn, yn = C.pair(70000 + r * 1000 + k, tc, d, 0.0)
    out = {}
    for en in ESTS:
        p = _G["e4"][(r, en)]["p"]
        s, tau, lab = C.score(xp, yp, p)
        out[en] = (s, C.score(xn, yn, p)[0], int(C.found(tau, lab, d)))
    return r, out


def item4():
    n, smoke = C.smoke_n(100)
    rs = [1, 4] if smoke else [1, 4, 8]
    t0 = time.time()
    with Pool(C.WORKERS) as pool:
        _G["e4"] = dict(pool.map(_setup4, [(r, en) for r in rs for en in ESTS]))
    print(f"setup {time.time() - t0:.0f}s", {k: (v["B"], C.pdict(v["p"])) for k, v in _G["e4"].items()}, flush=True)
    with Pool(C.WORKERS) as pool:
        res = pool.map(_task4, [(r, k) for r in rs for k in range(n)], chunksize=2)
    rows, out = [], {}
    for r in rs:
        for en in ESTS:
            pos, neg, hit = zip(*[o[en] for rr, o in res if rr == r])
            g = _G["e4"][(r, en)]
            e = dict(r=r, estimator=en, B=round(g["B"], 2), window=g["p"].window, max_lag=g["p"].max_lag, tolerance=g["p"].tolerance,
                     calib_failed=g["fail"], **C.auc_row(pos, neg), **C.hit_row(hit), chance=round(float(np.mean(neg)), 3))
            out[f"r{r}|{en}"] = e; rows.append(C.flat(e))
            print(e)
    out["wall_seconds"] = round(time.time() - t0)
    C.write_json("x6_item4_neff_estimators.json", out); C.write_csv("x6_item4_neff_estimators.csv", rows)
    print(f"total {time.time() - t0:.0f}s")


# ================================================================== item 5: per-link power
def _task5(args):
    r, kind, k = args
    tc, d = 2 * r, 5 * r
    cpl = kind == "coupled"
    x, y = C.pair((60000 if cpl else 70000) + r * 1000 + k, tc, d, 0.3 if cpl else 0.0)
    nrot = _G["nrot"]; n = len(x); out = {}
    for mi, (nm, p) in enumerate(_G["e5"][r]["M"].items()):
        rng = np.random.default_rng(9000000 + r * 100000 + (0 if cpl else 50000) + k * 10 + mi)
        real = C.score(x, y, p)[0]; lo = max(p.window, n // 10)
        null = [C.score(x, np.roll(y, int(rng.integers(lo, n - lo))), p)[0] for _ in range(nrot)]
        out[nm] = (1 + sum(v >= real for v in null)) / (nrot + 1)
    return r, kind, out


def item5():
    n50 = 2 if SMOKE else 50
    nrot = 9 if SMOKE else 99
    rs = [2, 4] if SMOKE else [2, 4, 8]
    t0 = time.time()
    with Pool(C.WORKERS) as pool:
        _G["e5"] = dict(pool.map(E1._setup, rs))
    _G["nrot"] = nrot
    keep = ("published", "rescaled", "rescaled_est", "atds")
    for r in rs:
        _G["e5"][r]["M"] = {k: _G["e5"][r]["M"][k] for k in keep}
    print(f"setup {time.time() - t0:.0f}s", {r: (g["r_hat"], {k: C.pdict(p) for k, p in g["M"].items()}) for r, g in _G["e5"].items()}, flush=True)
    with Pool(C.WORKERS) as pool:
        res = pool.map(_task5, [(r, kind, k) for r in rs for kind in ("coupled", "uncoupled") for k in range(n50)], chunksize=1)
    out, rows = {}, []
    for r in rs:
        for nm in keep:
            e = dict(r=r, method=nm, n_rotations=nrot, r_hat=_G["e5"][r]["r_hat"], params=C.pdict(_G["e5"][r]["M"][nm]))
            for kind in ("coupled", "uncoupled"):
                pv = [o[nm] for rr, kd, o in res if rr == r and kd == kind]
                kk = int(sum(v <= 0.05 for v in pv))
                e["power" if kind == "coupled" else "fpr"] = round(kk / len(pv), 3)
                e["power_ci" if kind == "coupled" else "fpr_ci"] = C.wilson(kk, len(pv))
                e["n_coupled" if kind == "coupled" else "n_uncoupled"] = len(pv)
            out[f"r{r}|{nm}"] = e; rows.append(C.flat(e)); print(e)
    out["wall_seconds"] = round(time.time() - t0)
    C.write_json("x6_item5_power.json", out); C.write_csv("x6_item5_power.csv", rows)
    print(f"total {time.time() - t0:.0f}s")


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    item = args[0]; fd = args[1] if len(args) > 1 else None
    C.set_results_for_smoke(SMOKE)
    t0 = time.time()
    if item in ("1", "all"):
        item1_item2(fd)
    if item == "2":
        item1_item2(None)    # regime part only matters; also refreshes the item-1 synthetic outputs
    if item in ("3", "all"):
        item3(fd)
    if item in ("4", "all"):
        item4()
    if item in ("5", "all"):
        item5()
    print(f"X6 item {item} wall {time.time() - t0:.0f}s")
