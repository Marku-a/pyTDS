"""E5: more real records in one field. airGR catchments L0123001 and L0123002 (daily) next to L0123003 (hourly).

Usage: ATDS_CODE=<checkout> python e5_hydrology.py <fields_data_dir> [<out_dir>]
<fields_data_dir>/raw must contain L0123001.rda, L0123002.rda, L0123003.rda (from
https://raw.githubusercontent.com/cran/airGR/master/data/<name>.rda) and the hourly CSV built by experiments/atds_real/fields.py.
Settings published and aTDS, 200 rotations (as Part F); the domain expectation (manifest) and the full-record
cross-correlation peak are reported separately.
"""
import os
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402
import fields as F  # noqa: E402

A = C.A
XLAG = {"daily": 30, "hourly": 72}   # cross-correlation search range in samples (30 days / 3 days)


def build_daily(fd, code):
    import pyreadr
    b = pyreadr.read_r(os.path.join(fd, "raw", f"{code}.rda"))["BasinObs"]
    out = os.path.join(fd, f"hydrology_airgr_{code}_daily.csv")
    pd.DataFrame({"date": pd.to_datetime(b.DatesR), "precip_mm": b.P, "temp_c": b["T"], "pet_mm": b.E,
                  "q_ls": b.Qls, "q_mm": b.Qmm}).to_csv(out, index=False)
    return out


def load_daily(code):
    def f():
        df = F._csv(f"hydrology_airgr_{code}_daily.csv")[["precip_mm", "q_mm"]]
        nn = int(df.q_mm.isna().sum())
        df = df.interpolate(limit_direction="both")
        return dict(s1=df.precip_mm.values, s2=df.q_mm.values, dt=F.DAY, names=("precipitation", "river flow"),
                    expected=("none", 0), prep=f"Precip and flow as is; {nn} NaNs in flow filled linearly.")
    return f


def main(fd, out_dir):
    F.FD = fd
    os.makedirs(out_dir, exist_ok=True)
    C.RESULTS = out_dir
    for code in ("L0123001", "L0123002"):
        build_daily(fd, code)
    # domain expectation from the manifest where given (delay of the flow after the rain, in samples); None = not given
    DATA = [("L0123001_daily", load_daily("L0123001"), "daily", "manifest: 0-3 days", (0.0, 3.0)),
            ("L0123002_daily", load_daily("L0123002"), "daily", "not given in the manifest", None),
            ("L0123003_hourly", F.d_hydro_hourly, "hourly", "manifest: positive, several hours (Part F used 8 h)", None)]
    t00 = time.time(); full, rows = {}, []
    for name, loader, kind, exp_text, exp_range in DATA:
        t0 = time.time()
        d = loader(); s1, s2 = np.asarray(d["s1"], float), np.asarray(d["s2"], float); d["s1"], d["s2"] = s1, s2
        sc, dl_x, _ = C.xcorr_full(s1, s2, XLAG[kind]); lag_x = -dl_x          # positive = flow lags rain
        entry = dict(n=len(s1), sampling_interval_s=d["dt"], preprocessing=d["prep"], domain_expectation=exp_text,
                     xcorr_full=dict(peak_lag_samples=lag_x, peak_lag_phys=F.phys(lag_x, d["dt"]), peak_abs_r=round(sc, 3),
                                     search_range_samples=XLAG[kind], peak_at_search_edge=bool(abs(lag_x) == XLAG[kind])), settings={})
        p_atds, fail = A.calib(A.window_params([s1, s2]), [(s1, s2)], n_shift=20)
        for sname, p, f in (("published", A.PUBLISHED, False), ("atds", p_atds, fail)):
            res = F.run_setting(d, sname, p, f)
            for k in ("expected_delay_samples", "match", "sign_ok"):
                res.pop(k, None)
            dl = res["delay_samples"]                                              # TDS convention: tau = -delay
            res["delay_after_rain_samples"] = None if dl is None else -dl
            res["in_domain_range"] = (None if exp_range is None or dl is None else bool(exp_range[0] <= -dl <= exp_range[1]))
            res["matches_xcorr_peak_within_1_or_10pct"] = (None if dl is None else bool(abs(-dl - lag_x) <= max(1, 0.1 * abs(lag_x))))
            entry["settings"][sname] = res
            rows.append(dict(catchment=name, setting=sname, n=len(s1), score=res["score"], null_mean=res["null_mean"], null95=res["null95"],
                             p=res["p"], n_windows=res["n_windows"], delay_samples=res["delay_samples"], delay_phys=res["delay_phys"],
                             domain_expectation=exp_text, in_domain_range=res["in_domain_range"],
                             xcorr_peak_lag_samples=lag_x, xcorr_peak_lag_phys=entry["xcorr_full"]["peak_lag_phys"],
                             matches_xcorr_peak=res["matches_xcorr_peak_within_1_or_10pct"],
                             **{f"{k}_samples": v for k, v in res["params_samples"].items()},
                             **{f"{k}_phys": v for k, v in res["params_phys"].items()}, calibration_failed=res["calibration_failed"]))
        full[name] = entry
        print(f"{name}: n={len(s1)} done in {time.time() - t0:.0f}s", flush=True)
    full["_wall_seconds"] = round(time.time() - t00)
    C.write_json("e5_hydrology.json", full); C.write_csv("e5_hydrology.csv", rows)
    for r in rows:
        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()})
    print(f"total {time.time() - t00:.0f}s")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else C.RESULTS)
