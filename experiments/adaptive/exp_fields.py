"""
Field datasets experiment: fixed default TDS parameters vs system-adaptive parameters.

Run: python experiments/adaptive/exp_fields.py
Outputs: results/fields.json, results/fields_summary.csv, results/fields_excerpts.json
"""
from __future__ import annotations

import csv
import json
import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import exp_real as ER  # noqa: E402
from adaptive import adaptive_params  # noqa: E402
from pyTDS.core import tds, time_delay_interaction  # noqa: E402
from pyTDS.params import TDSParams  # noqa: E402

FD = os.path.join(HERE, "data_fields")
RES = ER.RES
DAY = 86400.0
MONTH = 30.4375 * DAY


def _csv(name, **kw):
    return pd.read_csv(os.path.join(FD, name), **kw)


# ------------------------------------------------------------ dataset loaders
# Each returns dict(name, s1, s2, dt (seconds), names (leader, follower), expected, prep, note)
def d_furnace():
    df = _csv("engineering_gas_furnace.csv")
    return dict(s1=df.gas_rate.values, s2=df.co2_pct.values, dt=9.0, names=("gas feed rate", "CO2"),
                expected=("delay", 3), prep="As is (n=296, below the 500-sample minimum, so small and noisy).")


def d_climate():
    df = _csv("climate_nino34_gistemp.csv", parse_dates=["date"])[["NINO34_ANOM", "gistemp_loti"]].dropna()
    ma = df.rolling(121, center=True).mean()
    d = (df - ma).dropna()
    return dict(s1=d.NINO34_ANOM.values, s2=d.gistemp_loti.values, dt=MONTH, names=("Nino3.4", "global temperature"),
                expected=("delay", 3),
                prep="Dropped the 1 NaN month; subtracted a centred 121-month moving average from each series "
                     "(removes trend) and dropped the 120 edge months.")


def d_hydro_hourly():
    df = _csv("hydrology_airgr_L0123003_hourly.csv")[["precip_mm", "q_mm"]].interpolate(limit_direction="both")
    return dict(s1=df.precip_mm.values, s2=df.q_mm.values, dt=3600.0, names=("precipitation", "river flow"),
                expected=("delay", 8), prep="Precip and flow as is; NaNs filled linearly (none present).")


def d_hydro_daily():
    df = _csv("hydrology_airgr_L0123001_daily.csv")[["precip_mm", "q_mm"]]
    nn = int(df.q_mm.isna().sum())
    df = df.interpolate(limit_direction="both")
    return dict(s1=df.precip_mm.values, s2=df.q_mm.values, dt=DAY, names=("precipitation", "river flow"),
                expected=("delay", 2), prep=f"Precip and flow as is; {nn} NaNs in flow filled linearly.")


def d_energy():
    df = _csv("energy_vic_elec.csv")
    t = pd.to_datetime(df.Time, utc=True).dt.tz_convert("Australia/Melbourne")
    df = pd.DataFrame({"T": df.Temperature.values, "D": df.Demand.values,
                       "month": t.dt.month.values, "wk": (t.dt.weekday.values >= 5).astype(int),
                       "hh": (t.dt.hour.values * 2 + t.dt.minute.values // 30)})
    df = df[df.month.isin([12, 1, 2])].reset_index(drop=True)
    for c in ("T", "D"):
        df[c] = df[c] - df.groupby(["wk", "hh"])[c].transform("mean")
    return dict(s1=df["T"].values, s2=df["D"].values, dt=1800.0, names=("temperature", "electricity demand"),
                expected=("delay", 2),
                prep="Kept Dec/Jan/Feb only (local Melbourne time); subtracted the mean profile by (weekday/weekend, "
                     "half-hour of day) from both. The summers are concatenated into one series (TDS needs a "
                     "contiguous series), so there are 2 artificial junctions (Feb->Dec) and 1 jump at the record start.")


def d_beijing():
    df = _csv("atmosphere_beijing_pm25.csv", na_values=["NA"])
    nn = int(df["pm2.5"].isna().sum())
    pm = df["pm2.5"].interpolate(limit_direction="both")
    return dict(s1=df.Iws.values, s2=np.log1p(pm.values), dt=3600.0, names=("cumulated wind speed", "log PM2.5"),
                expected=("sign", None),
                prep=f"{nn} PM2.5 NaNs interpolated linearly (leading ones back-filled); follower = log1p(PM2.5), "
                     "leader = Iws (cumulated wind speed).")


def d_jena():
    df = _csv("weather_jena_climate_hourly.csv")
    return dict(s1=df["T (degC)"].values, s2=df["rh (%)"].values, dt=3600.0, names=("temperature", "relative humidity"),
                expected=("delay", 0), prep="T vs rh as is (zero-lag control, strong daily rhythm).", n_null=100)


def d_gasoline():
    df = _csv("economics_wti_vs_us_retail_gasoline_weekly.csv")[["wti_usd_bbl", "retail_gasoline_usd_gal"]]
    df = df.interpolate(limit_direction="both").diff().dropna()
    return dict(s1=df.wti_usd_bbl.values, s2=df.retail_gasoline_usd_gal.values, dt=7 * DAY,
                names=("crude oil price", "retail gasoline price"), expected=("delay", 1.5),
                prep="Weekly first differences of both prices (6 NaN retail weeks interpolated linearly before differencing).")


def d_eeg(_=None):
    d, s = ER.load_eeg()
    return dict(s1=d, s2=s, dt=1.0, names=("delta power", "sigma power"), expected=("delay", 0),
                prep="Sleep EEG delta and sigma band series as is (column 1 of data-*.txt, 1 s); leader unknown, "
                     "expected tau 0 (sign not assessed).", eeg=True)


def d_hrrsp():
    h, r, s = ER.hr_rsp("resting8")
    return dict(s1=h, s2=r, dt=0.25, names=("heart rate", "respiration"), expected=("unknown", None),
                prep="NeuroKit resting 8 min recording: R3 preprocessing from exp_real (ECG band-pass, R-peak HR "
                     "cubic-interpolated to 4 Hz, RSP low-passed and decimated to 4 Hz), z-scored; s1=HR, s2=RSP.")


DATASETS = [
    ("gas_furnace", d_furnace), ("climate_nino34_gistemp", d_climate), ("hydrology_hourly", d_hydro_hourly),
    ("hydrology_daily", d_hydro_daily), ("energy_vic_summer", d_energy), ("beijing_pm25", d_beijing),
    ("jena_T_rh", d_jena), ("gasoline_weekly_diff", d_gasoline), ("sleep_eeg_delta_sigma", d_eeg),
    ("neurokit_hr_rsp", d_hrrsp),
]


# ------------------------------------------------------------ helpers
def phys(x, dt):
    """Format samples*dt (seconds) in a readable unit."""
    s = abs(x) * dt
    sign = "-" if x < 0 else ""
    for unit, u in (("years", 365.25 * DAY), ("months", MONTH), ("weeks", 7 * DAY), ("days", DAY),
                    ("h", 3600.0), ("min", 60.0), ("s", 1.0)):
        if s >= u * (0.95 if unit != "s" else 0):
            # months only when native step is a month; weeks only when native step is a week
            if unit == "months" and dt < 0.9 * MONTH:
                continue
            if unit == "weeks" and dt < 0.9 * 7 * DAY:
                continue
            if unit == "years" and dt < 0.9 * MONTH:
                continue
            v = s / u
            return f"{sign}{v:.3g} {unit}"
    return f"{sign}{s:.3g} s"


def why_sentence(info, names, dt, params):
    B = np.array(info["bartlett"])
    k = int(np.argmax(B))
    L = info["window"]
    raw = 30 * B[k]
    s = (f"Slowest signal: {names[k]}, 1 independent sample per {phys(B[k], dt)} (B={B[k]:.3g}) "
         f"-> window = 30x{phys(B[k], dt)} = {phys(raw, dt)}")
    if info["window_clipped"]:
        s += f", clipped to {phys(L, dt)} ({L} samples)"
    pers = [p for p in info["periods"] if p is not None]
    base_cap = (L - 1) // 2
    if pers:
        pm = min(pers)
        s += f"; rhythm detected ({phys(pm, dt)}) -> max lag " + (
            f"capped at {phys(info['max_lag'], dt)}" if info["max_lag"] < base_cap else f"{phys(info['max_lag'], dt)} (cap not binding)")
    else:
        s += f"; no rhythm detected -> max lag = half window = {phys(info['max_lag'], dt)}"
    tol = info["tolerance"]
    if info.get("calibration_failed"):
        s += (f"; tolerance +-{phys(tol, dt) if tol else '0'}: even tolerance 0 gives fake stability above 5% "
              f"(calibration failed, null {info['null_score_by_tolerance'][0]:.1f}%).")
    else:
        s += f"; tolerance +-{phys(tol, dt) if tol else '0'} keeps fake stability <= 5%."
    return s


def match(exp, tau, dt):
    kind, val = exp
    if tau is None:
        return {"match": None, "sign_ok": None}
    if kind == "delay":
        tol = max(1.0, 0.5 * val)
        ok = abs(tau - (-val)) <= tol
        sign_ok = None if val == 0 else bool(tau < 0)
        return {"match": bool(ok), "sign_ok": sign_ok}
    if kind == "sign":
        return {"match": None, "sign_ok": bool(tau < 0)}
    return {"match": None, "sign_ok": None}


def zs(x):
    x = np.asarray(x, float)
    return (x - x.mean()) / x.std()


def run_variant(d, variant, params, info):
    s1, s2 = d["s1"], d["s2"]
    dt = d["dt"]
    ER.N_NULL = d.get("n_null", 200)
    r = tds(s1, s2, params)
    st = np.asarray(r["stable_taus"])
    med = float(np.median(st)) if len(st) else None
    nul = ER.null_stats(float(r["score"]), ER.null_scores(s1, s2, params))
    out = {"score": float(r["score"]), "n_windows": int(len(r["tau"])), "n_stable": int(len(st)),
           "null_mean": nul["null_mean"], "null95": nul["null_p95"], "p": nul["p_emp"], "z": nul["z"],
           "n_null": ER.N_NULL,
           "median_stable_tau_samples": med,
           "median_stable_tau_phys": None if med is None else phys(med, dt),
           "median_stable_tau_seconds": None if med is None else med * dt}
    out.update(match(d["expected"], med, dt))
    out["params_samples"] = {"window": params.window, "step": params.overlap, "max_lag": params.max_lag,
                             "tolerance": params.tolerance}
    out["params_phys"] = {"window": phys(params.window, dt), "step": phys(params.overlap, dt),
                          "max_lag": phys(params.max_lag, dt), "tolerance": phys(params.tolerance, dt) if params.tolerance else "0"}
    if variant != "fixed_default":
        out["info"] = {k: info.get(k) for k in ("bartlett", "tau_c", "periods", "window_clipped", "calibration_failed",
                                                "null_score_by_tolerance", "T")}
        out["why"] = why_sentence(info, d["names"], dt, params)
    else:
        out["why"] = "Fixed default TDSParams (window=60, step=30, max_lag=30, tolerance=1 sample), independent of the signals."
    return out, r


def excerpt(d, results):
    """Excerpt of ~300-600 samples where adaptive windows are stable; tau series of both variants."""
    T = len(d["s1"])
    n = min(600, T)
    ad = results["adaptive"]
    tau_a, t_a, _ = time_delay_interaction(d["s1"], d["s2"], ad["_params"])
    from pyTDS.core import stable_label
    lab = np.asarray(stable_label(tau_a, ad["_params"])).astype(bool)
    tc = np.asarray(t_a)[lab]
    best, start = -1, 0
    for s0 in range(0, T - n + 1, max(1, n // 10)):
        c = int(np.sum((tc >= s0) & (tc < s0 + n)))
        if c > best:
            best, start = c, s0
    seg = slice(start, start + n)
    out = {"dt_seconds": d["dt"], "names": list(d["names"]), "segment_start": int(start), "segment_n": int(n),
           "s1_z": zs(d["s1"][seg]).round(4).tolist(), "s2_z": zs(d["s2"][seg]).round(4).tolist(),
           "segment_stable_adaptive_windows": best, "variants": {}}
    for v, r in results.items():
        tau, tv, _ = time_delay_interaction(d["s1"], d["s2"], r["_params"])
        tau, tv = np.asarray(tau), np.asarray(tv)
        stab = np.asarray(stable_label(tau, r["_params"])).astype(bool)
        half = r["_params"].window / 2
        m = (tv + half > start) & (tv - half < start + n)  # windows overlapping the segment
        out["variants"][v] = {"window": r["_params"].window,
                              "centres": (tv[m] - start).tolist(), "tau": tau[m].tolist(),
                              "stable": stab[m].tolist()}
    return out


def main():
    os.makedirs(RES, exist_ok=True)
    full, rows, exc = {}, [], {}
    t00 = time.time()
    for name, loader in DATASETS:
        t0 = time.time()
        d = loader()
        s1, s2 = np.asarray(d["s1"], float), np.asarray(d["s2"], float)
        d["s1"], d["s2"] = s1, s2
        n = len(s1)
        variants = {"fixed_default": (TDSParams(), {})}
        variants["adaptive"] = adaptive_params([s1, s2], window_rule="bartlett", tolerance="calibrated")
        if d.get("eeg"):
            variants["adaptive_cap5min"] = adaptive_params([s1, s2], window_rule="bartlett", tolerance="calibrated",
                                                           max_window=300)
        entry = {"n": n, "sampling_interval_s": d["dt"], "leader": d["names"][0], "follower": d["names"][1],
                 "expected": {"kind": d["expected"][0], "value_samples": d["expected"][1]},
                 "preprocessing": d["prep"], "variants": {}}
        keep = {}
        for v, (p, info) in variants.items():
            res, _ = run_variant(d, v, p, info)
            entry["variants"][v] = res
            keep[v] = {"_params": p}
            rows.append({"dataset": name, "variant": v, "n": n, "score": res["score"], "null_mean": res["null_mean"],
                         "null95": res["null95"], "p": res["p"], "n_null": res["n_null"], "n_windows": res["n_windows"],
                         "tau_samples": res["median_stable_tau_samples"], "tau_phys": res["median_stable_tau_phys"],
                         "match": res["match"], "sign_ok": res["sign_ok"],
                         **{f"{k}_samples": x for k, x in res["params_samples"].items()},
                         **{f"{k}_phys": x for k, x in res["params_phys"].items()},
                         "calibration_failed": res.get("info", {}).get("calibration_failed"),
                         "window_clipped": res.get("info", {}).get("window_clipped"),
                         "why": res["why"]})
        full[name] = entry
        exc[name] = excerpt(d, keep)
        print(f"{name}: n={n} done in {time.time() - t0:.0f}s", flush=True)

    def h(o):
        if isinstance(o, np.integer):
            return int(o)
        if isinstance(o, np.floating):
            return float(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, (np.bool_,)):
            return bool(o)
        return str(o)

    json.dump(full, open(os.path.join(RES, "fields.json"), "w"), indent=1, default=h)
    json.dump(exc, open(os.path.join(RES, "fields_excerpts.json"), "w"), default=h)
    keys = list(rows[0])
    with open(os.path.join(RES, "fields_summary.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    cols = ["dataset", "variant", "n", "score", "null95", "p", "tau_phys", "match", "sign_ok", "window_phys",
            "tolerance_phys", "calibration_failed"]
    ER.table(rows, cols)
    print(f"total {time.time() - t00:.0f}s")


if __name__ == "__main__":
    main()
