"""
Real-data experiments R1-R4: system-adaptive TDS params vs fixed defaults.

Usage: python experiments/adaptive/exp_real.py --exp R1|R2|R3|R4|all
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from dataclasses import asdict
from multiprocessing import Pool

import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import butter, filtfilt, find_peaks, resample_poly

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

from adaptive import adaptive_params  # noqa: E402
from pyTDS.core import tds  # noqa: E402
from pyTDS.params import TDSParams  # noqa: E402

DATA = os.path.join(ROOT, "data")
EXT = os.path.join(HERE, "data_external")
RES = os.path.join(HERE, "results")
FIG = os.path.join(ROOT, "adaptive-TDS reports", "figures")
N_NULL = 200
VARIANTS = ["fixed_default", "adaptive_fixedtol", "adaptive_calibrated"]
COLORS = {"fixed_default": "tab:gray", "adaptive_fixedtol": "tab:blue", "adaptive_calibrated": "tab:red"}
RECS = {"resting5": "bio_resting_5min_100hz.csv", "resting8": "bio_resting_8min_100hz.csv",
        "event": "bio_eventrelated_100hz.csv"}

_CTX: dict = {}


# ---------------------------------------------------------------- helpers
def z(x):
    x = np.asarray(x, float)
    return (x - x.mean()) / x.std()


def make_params(system, variant, **kw):
    """Return (params, info) for a variant; system = list of equal-length signals."""
    if variant == "fixed_default":
        p = TDSParams()
        return p, {"window": p.window, "step": p.overlap, "max_lag": p.max_lag, "tolerance": p.tolerance}
    tol = "fixed" if variant == "adaptive_fixedtol" else "calibrated"
    return adaptive_params(system, tolerance=tol, **kw)


def _null_worker(sh):
    r = tds(_CTX["s1"], np.roll(_CTX["s2"], int(sh)), _CTX["p"])
    return r["score"]


def null_scores(s1, s2, params, procs=4):
    T = len(s1)
    rng = np.random.default_rng(0)
    shifts = rng.integers(T // 10, T - T // 10 + 1, N_NULL)
    _CTX.update(s1=s1, s2=s2, p=params)
    with Pool(procs) as pool:
        return np.array(pool.map(_null_worker, shifts, chunksize=5))


def null_stats(real, null):
    sd = null.std(ddof=1)
    return {"null_mean": float(null.mean()), "null_sd": float(sd),
            "null_p95": float(np.percentile(null, 95)),
            "p_emp": float((1 + np.sum(null >= real)) / (1 + len(null))),
            "z": float((real - null.mean()) / sd) if sd > 0 else None}


def run_pair(s1, s2, params, null=True):
    r = tds(s1, s2, params)
    st = np.asarray(r["stable_taus"])
    out = {"score": float(r["score"]),
           "median_stable_tau_samples": float(np.median(st)) if len(st) else None,
           "n_windows": int(len(r["tau"])), "n_stable": int(len(st))}
    if null:
        out.update(null_stats(out["score"], null_scores(s1, s2, params)))
    return out, r


def pinfo(params, info):
    d = {"window": params.window, "step": params.overlap, "max_lag": params.max_lag,
         "tolerance": params.tolerance}
    d["info"] = info
    return d


def bandpass(x, fs, lo=0.5, hi=20.0):
    b, a = butter(4, [lo / (fs / 2), hi / (fs / 2)], btype="band")
    return filtfilt(b, a, x)


def jdump(obj, name):
    os.makedirs(RES, exist_ok=True)

    def h(o):
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        return str(o)

    with open(os.path.join(RES, f"real_{name}.json"), "w") as f:
        json.dump(obj, f, indent=1, default=h)


def write_csv(rows, name):
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(os.path.join(RES, f"real_{name}_summary.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def table(rows, cols):
    print("  " + " | ".join(f"{c:>14}" for c in cols))
    for r in rows:
        print("  " + " | ".join(f"{_fmt(r.get(c)):>14}" for c in cols))


def _fmt(v):
    if isinstance(v, float):
        return f"{v:.3g}"
    return str(v)


def plt_():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(FIG, exist_ok=True)
    return plt


def flat(base, res, params, info):
    row = dict(base)
    row.update({k: v for k, v in res.items()})
    row.update(window=params.window, step=params.overlap, max_lag=params.max_lag, tolerance=params.tolerance)
    return row


# ---------------------------------------------------------------- R1
def load_eeg():
    d = np.loadtxt(os.path.join(DATA, "data-delta.txt"))[:, 1]
    s = np.loadtxt(os.path.join(DATA, "data-sigma.txt"))[:, 1]
    return d, s


def exp_R1():
    d, s = load_eeg()
    out, rows = {}, []
    configs = [(v, v, {}) for v in VARIANTS]
    configs += [(f"adaptive_calibrated_cycles{c}", "adaptive_calibrated", {"cycles": c}) for c in (5, 20)]
    configs.append(("fixed_default_anchor_end", "fixed_default", {"stability_anchor": "end"}))
    for name, var, kw in configs:
        if name == "fixed_default_anchor_end":  # MATLAB-style rule: reproduces ~43.6%
            p, info = TDSParams(stability_anchor="end"), {"note": "stability_anchor=end"}
        else:
            p, info = make_params([d, s], var, **kw)
        res, _ = run_pair(d, s, p)
        out[name] = {"result": res, "params": pinfo(p, info)}
        rows.append(flat({"variant": name}, res, p, info))
        print(f"  R1 {name}: score={res['score']:.2f} null95={res['null_p95']:.2f} p={res['p_emp']:.4f}")
    jdump(out, "R1")
    write_csv(rows, "R1")
    table(rows, ["variant", "score", "null_mean", "null_p95", "p_emp", "z", "median_stable_tau_samples",
                 "window", "step", "max_lag", "tolerance"])
    plt = plt_()
    fig, ax = plt.subplots(figsize=(8, 4.5))
    x = np.arange(len(rows))
    ax.bar(x - 0.2, [r["score"] for r in rows], 0.4, label="real TDS score")
    ax.bar(x + 0.2, [r["null_p95"] for r in rows], 0.4, label="null 95th pct", color="tab:orange")
    ax.axhline(43.6, ls="--", c="k", lw=1, label="MATLAB benchmark 43.6%")
    ax.set_xticks(x)
    ax.set_xticklabels([r["variant"].replace("adaptive_", "adapt_\n") for r in rows], fontsize=8)
    ax.set_ylabel("TDS score (%)")
    ax.set_title("R1: EEG delta vs sigma, real vs circular-shift null (200 shifts)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "R1_eeg_score_vs_null.png"), dpi=130)
    plt.close(fig)
    return out


# ---------------------------------------------------------------- R2
def detect_r(ecg, fs):
    pk, _ = find_peaks(ecg, distance=int(0.4 * fs), height=np.percentile(ecg, 98) * 0.5)
    return pk


def exp_R2():
    df = np.genfromtxt(os.path.join(EXT, RECS["resting5"]), delimiter=",", names=True)
    fs0 = 100
    ecg_f = bandpass(df["ECG"], fs0)
    ppg_f = bandpass(df["PPG"], fs0)
    ecg, ppg = z(ecg_f), z(ppg_f)
    rp = detect_r(ecg, fs0)
    pp, _ = find_peaks(ppg, distance=40, prominence=0.5 * np.std(ppg))
    hr = 60 * fs0 / np.diff(rp)
    delays = []
    for r in rp:
        c = pp[(pp >= r + 0.1 * fs0) & (pp <= r + 0.6 * fs0)]
        if len(c):
            delays.append((c[0] - r) * 1000 / fs0)
    delays = np.array(delays)
    pat = float(np.median(delays))
    q1, q3 = np.percentile(delays, [25, 75])
    sanity = {"n_rpeaks": int(len(rp)), "n_ppg_peaks": int(len(pp)), "hr_median_bpm": float(np.median(hr)),
              "hr_p5_bpm": float(np.percentile(hr, 5)), "hr_p95_bpm": float(np.percentile(hr, 95)),
              "n_pat_pairs": int(len(delays)), "PAT_ref_ms": pat, "PAT_q1_ms": float(q1), "PAT_q3_ms": float(q3)}
    print("  R2 sanity:", sanity)
    out = {"sanity": sanity, "rates": {}}
    rows, perwin = [], {}
    for fs, q in ((100, 1), (50, 2), (25, 4)):
        e = ecg if q == 1 else z(resample_poly(ecg_f, 1, q))
        p_ = ppg if q == 1 else z(resample_poly(ppg_f, 1, q))
        out["rates"][fs] = {}
        for v in VARIANTS:
            par, info = make_params([e, p_], v)
            res, r = run_pair(e, p_, par)
            st = np.asarray(r["stable_taus"])
            ms = st * 1000.0 / fs
            res["median_stable_tau_ms"] = float(np.median(ms)) if len(ms) else None
            res["frac_within_50ms_of_-PAT"] = float(np.mean(np.abs(ms + pat) <= 50)) if len(ms) else None
            out["rates"][fs][v] = {"result": res, "params": pinfo(par, info)}
            rows.append(flat({"fs": fs, "variant": v}, res, par, info))
            if fs == 100:
                tau = np.asarray(r["tau"]) * 1000.0 / fs
                t = (np.arange(len(tau)) * par.overlap + par.window / 2) / fs
                perwin[v] = (t, tau, np.asarray(r["tau"]) )
            print(f"  R2 fs={fs} {v}: score={res['score']:.1f} p={res['p_emp']:.3f} med_tau_ms={res['median_stable_tau_ms']}")
    jdump(out, "R2")
    write_csv(rows, "R2")
    table(rows, ["fs", "variant", "score", "null_p95", "p_emp", "median_stable_tau_ms",
                 "frac_within_50ms_of_-PAT", "window", "step", "max_lag", "tolerance"])
    plt = plt_()
    fig, axs = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
    for ax, v in zip(axs, VARIANTS):
        t, tau, _ = perwin[v]
        ax.plot(t, tau, ".", ms=3, color=COLORS[v], label="per-window tau")
        ax.axhline(-pat, c="k", ls="--", label=f"-PAT_ref = {-pat:.0f} ms")
        ax.set_ylabel("tau (ms)")
        pr = out["rates"][100][v]
        ax.set_title(f"{v} (score {pr['result']['score']:.1f}%, window {pr['params']['window']}, "
                     f"max_lag {pr['params']['max_lag']}, tol {pr['params']['tolerance']})", fontsize=9)
        ax.legend(loc="upper right", fontsize=8)
    axs[-1].set_xlabel("time (s)")
    fig.suptitle("R2: ECG->PPG per-window delay at 100 Hz (s1=ECG, s2=PPG)")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "R2_pat_tau_over_time.png"), dpi=130)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for v in VARIANTS:
        ys = [out["rates"][fs][v]["result"]["median_stable_tau_ms"] for fs in (100, 50, 25)]
        ax.plot([100, 50, 25], [np.nan if y is None else y for y in ys], "o-", color=COLORS[v], label=v)
    ax.axhline(-pat, c="k", ls="--", label="-PAT_ref")
    ax.set_xlabel("sampling rate (Hz)")
    ax.set_ylabel("median stable tau (ms)")
    ax.set_title("R2: sampling-rate invariance of the estimated delay")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "R2_rate_invariance.png"), dpi=130)
    plt.close(fig)
    return out


# ---------------------------------------------------------------- R3
def hr_rsp(name):
    df = np.genfromtxt(os.path.join(EXT, RECS[name]), delimiter=",", names=True)
    fs = 100
    ecg = z(bandpass(df["ECG"], fs))
    rp = detect_r(ecg, fs)
    tr = rp / fs
    hr = 60.0 / np.diff(tr)
    f = interp1d(tr[1:], hr, kind="cubic")
    k = np.arange(int(np.ceil(tr[1] * 4)), int(np.floor(tr[-1] * 4)) + 1)
    hr4 = f(k / 4.0)
    b, a = butter(4, 1.0 / (fs / 2))
    rsp4 = resample_poly(filtfilt(b, a, df["RSP"]), 1, 25)  # sample j at j/4 s
    k = k[k < len(rsp4)]
    hr4, rsp4 = hr4[: len(k)], rsp4[k]
    return z(hr4), z(rsp4), {"n_rpeaks": int(len(rp)), "hr_median_bpm": float(np.median(hr)),
                             "hr_min": float(hr.min()), "hr_max": float(hr.max()), "n4hz": int(len(k))}


def exp_R3():
    sig, san = {}, {}
    for n in RECS:
        h, r, s = hr_rsp(n)
        sig[n] = (h, r)
        san[n] = s
        print("  R3 sanity", n, s)
    T = min(len(x) for p in sig.values() for x in p)
    system = [x[:T] for n in RECS for x in sig[n]]
    names = list(RECS)
    out = {"sanity": san, "variants": {}}
    rows = []
    for v in VARIANTS:
        par, info = make_params(system, v)
        out["variants"][v] = {"params": pinfo(par, info), "true": {}, "cross": {}}
        for n in names:
            h, r = sig[n]
            res, _ = run_pair(h, r, par)
            mt = res["median_stable_tau_samples"]
            res["median_stable_tau_s"] = None if mt is None else mt / 4.0
            out["variants"][v]["true"][n] = res
            rows.append(flat({"variant": v, "pair": f"true:{n}"}, res, par, info))
        for a in names:
            for b in names:
                if a == b:
                    continue
                m = min(len(sig[a][0]), len(sig[b][1]))
                res, _ = run_pair(sig[a][0][:m], sig[b][1][:m], par, null=False)
                out["variants"][v]["cross"][f"HR_{a}-RSP_{b}"] = res
                rows.append(flat({"variant": v, "pair": f"cross:HR_{a}-RSP_{b}"}, res, par, info))
        tr = [out["variants"][v]["true"][n]["score"] for n in names]
        cr = [c["score"] for c in out["variants"][v]["cross"].values()]
        out["variants"][v]["summary"] = {"true_mean": float(np.mean(tr)), "cross_mean": float(np.mean(cr)),
                                         "cross_max": float(np.max(cr))}
        print(f"  R3 {v}: true={np.round(tr, 1)} cross mean={np.mean(cr):.1f} max={np.max(cr):.1f}")
    jdump(out, "R3")
    write_csv(rows, "R3")
    table([r for r in rows if r["pair"].startswith("true")],
          ["variant", "pair", "score", "null_p95", "p_emp", "median_stable_tau_s", "window", "step", "max_lag", "tolerance"])
    plt = plt_()
    fig, ax = plt.subplots(figsize=(9, 4.5))
    x = np.arange(len(VARIANTS))
    for i, n in enumerate(names):
        ax.bar(x + (i - 1) * 0.12 - 0.2, [out["variants"][v]["true"][n]["score"] for v in VARIANTS], 0.12,
               label=f"true {n}")
    ax.bar(x + 0.3, [out["variants"][v]["summary"]["cross_mean"] for v in VARIANTS], 0.2, color="k",
           alpha=0.6, label="cross-recording mean (6)")
    ax.scatter(np.repeat(x + 0.3, 6), [c["score"] for v in VARIANTS for c in out["variants"][v]["cross"].values()],
               c="w", edgecolors="k", s=14, zorder=3, label="cross-recording pairs")
    ax.set_xticks(x)
    ax.set_xticklabels(VARIANTS, fontsize=8)
    ax.set_ylabel("TDS score (%)")
    ax.set_title("R3: HR vs respiration, true pairs vs cross-recording null")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "R3_hr_rsp_vs_cross.png"), dpi=130)
    plt.close(fig)
    return out


# ---------------------------------------------------------------- R4
def exp_R4():
    d, s = load_eeg()
    d, s = z(d), z(s)
    factors = [("down2", 1, 2), ("orig", 1, 1), ("up2", 2, 1), ("up4", 4, 1)]
    out, rows = {}, []
    for lab, up, dn in factors:
        a = resample_poly(d, up, dn) if (up, dn) != (1, 1) else d
        b = resample_poly(s, up, dn) if (up, dn) != (1, 1) else s
        a, b = z(a), z(b)
        fs = up / dn
        out[lab] = {}
        for v in VARIANTS:
            par, info = make_params([a, b], v)
            res, _ = run_pair(a, b, par, null=False)
            mt = res["median_stable_tau_samples"]
            res["median_stable_tau_s"] = None if mt is None else mt / fs
            out[lab][v] = {"result": res, "params": pinfo(par, info)}
            rows.append(flat({"factor": fs, "variant": v}, res, par, info))
            print(f"  R4 {lab} {v}: score={res['score']:.1f} tau_s={res['median_stable_tau_s']} win={par.window}")
    jdump(out, "R4")
    write_csv(rows, "R4")
    table(rows, ["factor", "variant", "score", "median_stable_tau_s", "window", "step", "max_lag", "tolerance"])
    plt = plt_()
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for v in VARIANTS:
        pts = sorted((r["factor"], r["score"]) for r in rows if r["variant"] == v)
        ax.plot([p[0] for p in pts], [p[1] for p in pts], "o-", color=COLORS[v], label=v)
    ax.set_xscale("log", base=2)
    ax.set_xlabel("sampling factor (x original 1 Hz)")
    ax.set_ylabel("TDS score (%)")
    ax.set_title("R4: EEG delta vs sigma, score vs resampling factor")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "R4_rate_invariance.png"), dpi=130)
    plt.close(fig)
    return out


EXPS = {"R1": exp_R1, "R2": exp_R2, "R3": exp_R3, "R4": exp_R4}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exp", default="all", choices=["R1", "R2", "R3", "R4", "all"])
    a = ap.parse_args()
    for k in (EXPS if a.exp == "all" else [a.exp]):
        t0 = time.time()
        print(f"== {k} ==")
        EXPS[k]()
        print(f"== {k} done in {time.time() - t0:.0f}s ==")


if __name__ == "__main__":
    main()
