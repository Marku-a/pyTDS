"""Data for the aTDS evidence reader (all reports in one page): aggregates only.

Sources (all committed):
- pyTDS branch claude/practical-mccarthy-n2mec6, experiments/adaptive/results/*_summary.csv (Reports 02-07)
- this branch: experiments/adaptive_groundtruth/results/ (Sleep-EDF, Fantasia, treadmill; ground-truth study)
- this branch: paper/figures/figdata.json (Bashan-2012 cohort, apnea, Garmin)

usage: python reader_data.py --out ../figures/reader.json
"""
import argparse
import io
import json
import os
import subprocess

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PRIOR = "origin/claude/practical-mccarthy-n2mec6"
GT = os.path.join(ROOT, "experiments", "adaptive_groundtruth", "results")
from figdata import r  # noqa: E402


def prior_csv(name):
    txt = subprocess.check_output(["git", "-C", ROOT, "show", f"{PRIOR}:experiments/adaptive/results/{name}"], text=True)
    return pd.read_csv(io.StringIO(txt))


def synthetic():
    out = {}
    e1 = prior_csv("synthetic_E1_summary.csv")
    e1 = e1[e1.condition.str.endswith("c=0.3")].copy()
    e1["r"] = e1.condition.str.extract(r"r=(\d+)").astype(int)
    out["E1_weak_coupling_vs_slowness"] = {v: e1[e1.variant == v].sort_values("r")[["r", "auc", "hit_rate", "window", "tolerance"]].to_dict("list")
                                           for v in ("fixed_default", "v2_calibrated", "oracle")}
    e2 = prior_csv("synthetic_E2_summary.csv")
    e2["delay"] = e2.condition.str.extract(r"delay=(\d+)").astype(int)
    out["E2_rulkov_vs_delay"] = {v: e2[e2.variant == v].sort_values("delay")[["delay", "auc", "hit_rate", "window", "tolerance"]].to_dict("list")
                                 for v in ("fixed_default", "v2_calibrated", "oracle")}
    e1c = prior_csv("synthetic_E1c_summary.csv")
    out["E1c_mixed_speeds"] = e1c[e1c.variant.isin(["fixed_default", "v2_calibrated", "oracle"])][["condition", "variant", "auc", "hit_rate"]].to_dict("records")
    e1d = prior_csv("synthetic_E1d_summary.csv")
    out["E1d_on_off"] = e1d[e1d.variant.isin(["fixed_default", "v2_calibrated"])][["condition", "variant", "balanced_acc_mean"]].to_dict("records")
    r4 = prior_csv("real_R4_summary.csv")
    out["R4_sampling_rate"] = {v: r4[r4.variant == v].sort_values("factor")[["factor", "score", "window", "tolerance"]].to_dict("list")
                               for v in ("fixed_default", "v2_calibrated")}
    r2 = prior_csv("real_R2_summary.csv")
    out["R2_ecg_ppg"] = r2[r2.variant.isin(["fixed_default", "v2_calibrated"])][["fs", "variant", "score", "null_p95", "median_stable_tau_ms"]].to_dict("records")
    r3 = prior_csv("real_R3_summary.csv")
    out["R3_hr_resp"] = r3[r3.variant.isin(["fixed_default", "v2_calibrated"]) & r3.pair.str.startswith("true")][["variant", "pair", "score", "null_p95", "window", "tolerance"]].to_dict("records")
    return out


def fields():
    f = prior_csv("fields_summary.csv")
    # Report 07 ran sleep EEG with a 5-min cap; use that row as "adaptive" there.
    f = f[~((f.dataset == "sleep_eeg_delta_sigma") & (f.variant == "adaptive"))].copy()
    f.loc[f.variant == "adaptive_cap5min", "variant"] = "adaptive"
    f = f[f.variant.isin(["fixed_default", "adaptive"])]
    return f[["dataset", "variant", "n", "score", "null95", "p", "tau_phys", "match", "window_phys", "tolerance_phys", "calibration_failed"]].to_dict("records")


def groundtruth():
    e = pd.read_csv(os.path.join(GT, "effects_bootstrap.csv"))
    e = e[e.variant.isin(["published", "tol_only", "win_only", "v2_calibrated"])]
    ind = pd.read_csv(os.path.join(GT, "independent_samples.csv"))
    fs = pd.read_csv(os.path.join(GT, "fantasia_summary.csv"))
    return {"effects": e.to_dict("records"), "independent_samples": ind.to_dict("records"),
            "fantasia": fs[fs.variant.isin(["published", "v2_calibrated"])][["variant", "n_subjects", "detect", "xsubj_fpr", "lag_s_median"]].to_dict("records")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = {"synthetic": synthetic(), "fields": fields(), "groundtruth": groundtruth()}
    json.dump(r(out), open(a.out, "w"), separators=(",", ":"))
    print(a.out, os.path.getsize(a.out) // 1024, "KB")


if __name__ == "__main__":
    main()
