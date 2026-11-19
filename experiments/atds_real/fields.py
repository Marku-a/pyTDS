"""Part F: aTDS v2.1 vs published parameters on public datasets with known delays.

Scope (user decision): only three of the eight manifest datasets are run:
hydrology hourly (airGR L0123003), energy (vic_elec, summer/diurnal-profile removed), economics (WTI -> retail gasoline).

Usage: python fields.py <data_dir> <out_dir>
Downloads the sources listed in data_fields/MANIFEST.json into <data_dir>/raw and rebuilds the CSVs in <data_dir>.
Loaders d_*, helpers match/phys are copied from exp_fields.py on branch claude/practical-mccarthy-n2mec6.
"""
from __future__ import annotations

import csv
import json
import os
import sys
import time
import urllib.request

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import atds21 as A  # noqa: E402
from pyTDS.core import stable_label, tds_score, time_delay_interaction  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DAY = 86400.0
MONTH = 30.4375 * DAY
R = "https://raw.githubusercontent.com"
SOURCES = {
    "L0123003.rda": f"{R}/cran/airGR/master/data/L0123003.rda",
    "vic_elec.csv": f"{R}/vincentarelbundock/Rdatasets/master/csv/tsibbledata/vic_elec.csv",
    "wti-weekly.csv": f"{R}/datasets/oil-prices/main/data/wti-weekly.csv",
    "Retail.csv": f"{R}/mihirpat98/Gas-Price-Prediction-Time-Series/13bdd4671572326737b484b1570bec37b9921877/data_raw/Retail.csv",
}


def fetch(raw_dir):
    os.makedirs(raw_dir, exist_ok=True)
    for name, url in SOURCES.items():
        dst = os.path.join(raw_dir, name)
        if not os.path.exists(dst):
            urllib.request.urlretrieve(url, dst)


def build_csvs(fd):
    raw = os.path.join(fd, "raw")
    fetch(raw)
    import pyreadr
    b = pyreadr.read_r(os.path.join(raw, "L0123003.rda"))["BasinObs"]
    pd.DataFrame({"date": pd.to_datetime(b.DatesR), "precip_mm": b.P, "temp_c": b["T"], "pet_mm": b.E,
                  "q_ls": b.Qls, "q_mm": b.Qmm}).to_csv(os.path.join(fd, "hydrology_airgr_L0123003_hourly.csv"), index=False)
    v = pd.read_csv(os.path.join(raw, "vic_elec.csv"))
    v[["Time", "Demand", "Temperature"]].to_csv(os.path.join(fd, "energy_vic_elec.csv"), index=False)
    w = pd.read_csv(os.path.join(raw, "wti-weekly.csv"), parse_dates=["Date"]).rename(columns={"Date": "wti_date", "Price": "wti_usd_bbl"})
    rt = pd.read_csv(os.path.join(raw, "Retail.csv"), skiprows=3, header=None, names=["date", "retail_gasoline_usd_gal"],
                     usecols=[0, 1])
    rt["date"] = pd.to_datetime(rt["date"], errors="coerce")
    rt = rt.dropna(subset=["date"]).sort_values("date")
    m = pd.merge_asof(rt, w.sort_values("wti_date"), left_on="date", right_on="wti_date", direction="backward",
                      tolerance=pd.Timedelta(days=10))
    m.to_csv(os.path.join(fd, "economics_wti_vs_us_retail_gasoline_weekly.csv"), index=False)


FD = None


def _csv(name, **kw):
    return pd.read_csv(os.path.join(FD, name), **kw)


# ---- loaders copied from exp_fields.py
def d_hydro_hourly():
    df = _csv("hydrology_airgr_L0123003_hourly.csv")[["precip_mm", "q_mm"]].interpolate(limit_direction="both")
    return dict(s1=df.precip_mm.values, s2=df.q_mm.values, dt=3600.0, names=("precipitation", "river flow"),
                expected=("delay", 8), prep="Precip and flow as is; NaNs filled linearly (none present).")


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
                     "half-hour of day) from both. The summers are concatenated into one series.")


def d_gasoline():
    df = _csv("economics_wti_vs_us_retail_gasoline_weekly.csv")[["wti_usd_bbl", "retail_gasoline_usd_gal"]]
    df = df.interpolate(limit_direction="both").diff().dropna()
    return dict(s1=df.wti_usd_bbl.values, s2=df.retail_gasoline_usd_gal.values, dt=7 * DAY,
                names=("crude oil price", "retail gasoline price"), expected=("delay", 1.5),
                prep="Weekly first differences of both prices (NaN retail weeks interpolated linearly before differencing).")


DATASETS = [("hydrology_hourly", d_hydro_hourly), ("energy_vic_summer", d_energy), ("gasoline_weekly_diff", d_gasoline)]


def phys(x, dt):
    """Format samples*dt (seconds) in a readable unit."""
    s = abs(x) * dt
    sign = "-" if x < 0 else ""
    for unit, u in (("years", 365.25 * DAY), ("months", MONTH), ("weeks", 7 * DAY), ("days", DAY),
                    ("h", 3600.0), ("min", 60.0), ("s", 1.0)):
        if s >= u * (0.95 if unit != "s" else 0):
            if unit == "months" and dt < 0.9 * MONTH:
                continue
            if unit == "weeks" and dt < 0.9 * 7 * DAY:
                continue
            if unit == "years" and dt < 0.9 * MONTH:
                continue
            v = s / u
            return f"{sign}{v:.3g} {unit}"
    return f"{sign}{s:.3g} s"


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


def run_setting(d, name, p, fail):
    s1, s2, dt = d["s1"], d["s2"], d["dt"]
    tau = time_delay_interaction(s1, s2, p)[0]
    lab = stable_label(tau, p)
    score = float(tds_score(lab))
    null = A.rotation_null(s1, s2, p, 200, seed=0)
    dl = A.delay_of(tau, lab)
    out = {"score": score, "n_windows": int(len(tau)), "n_stable": int(lab.sum()), **A.null_stats(score, null),
           "n_null": 200, "delay_samples": dl, "delay_phys": None if dl is None else phys(dl, dt),
           "expected_delay_samples": d["expected"][1], **match(d["expected"], dl, dt),
           "params_samples": {"window": p.window, "step": p.overlap, "max_lag": p.max_lag, "tolerance": p.tolerance},
           "params_phys": {"window": phys(p.window, dt), "step": phys(p.overlap, dt), "max_lag": phys(p.max_lag, dt),
                           "tolerance": phys(p.tolerance, dt) if p.tolerance else "0"},
           "calibration_failed": bool(fail)}
    return out


def main(data_dir, out_dir):
    global FD
    FD = data_dir
    os.makedirs(out_dir, exist_ok=True)
    build_csvs(data_dir)
    manifest = {m["file"]: m for m in json.load(open(os.path.join(HERE, "data_fields", "MANIFEST.json")))}
    full, rows = {}, []
    t00 = time.time()
    for name, loader in DATASETS:
        t0 = time.time()
        d = loader()
        s1, s2 = np.asarray(d["s1"], float), np.asarray(d["s2"], float)
        d["s1"], d["s2"] = s1, s2
        p_atds = A.window_params([s1, s2])
        p_atds, fail = A.calib(p_atds, [(s1, s2)], n_shift=20)
        entry = {"n": len(s1), "sampling_interval_s": d["dt"], "leader": d["names"][0], "follower": d["names"][1],
                 "expected": {"kind": d["expected"][0], "value_samples": d["expected"][1]},
                 "preprocessing": d["prep"], "settings": {}}
        for sname, p, f in (("published", A.PUBLISHED, False), ("atds", p_atds, fail)):
            res = run_setting(d, sname, p, f)
            entry["settings"][sname] = res
            rows.append({"dataset": name, "setting": sname, "n": len(s1), "score": res["score"], "null_mean": res["null_mean"],
                         "null95": res["null95"], "p": res["p"], "n_windows": res["n_windows"],
                         "delay_samples": res["delay_samples"], "delay_phys": res["delay_phys"],
                         "expected_delay_samples": res["expected_delay_samples"], "match": res["match"],
                         **{f"{k}_samples": v for k, v in res["params_samples"].items()},
                         **{f"{k}_phys": v for k, v in res["params_phys"].items()},
                         "calibration_failed": res["calibration_failed"]})
        full[name] = entry
        print(f"{name}: n={len(s1)} done in {time.time() - t0:.0f}s", flush=True)
    full["_wall_seconds"] = time.time() - t00
    json.dump(full, open(os.path.join(out_dir, "fields_v21.json"), "w"), indent=1, default=str)
    with open(os.path.join(out_dir, "fields_v21_summary.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    for r in rows:
        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()})
    print(f"total {time.time() - t00:.0f}s")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "results"))
