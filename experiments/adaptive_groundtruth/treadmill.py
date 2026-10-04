"""
D3 Treadmill CPET (imposed speed -> heart rate / VO2), see SPEC.md.

Data: physionet treadmill-exercise-cardioresp 1.0.1 (test_measure.csv + subject-info.csv).
Raw data goes to --data-dir (default /home/user/data_gt/treadmill), results to --out.

Notes / closest-option choices:
- The CSV header is lower-case ("time"), the spec says "Time s"; columns are matched by the real header.
- ID_test is sorted as a string ("sorted order"). One test per subject ID (first in sorted order) is applied
  BEFORE the every-k-th thinning, so the thinned list really has 200 distinct subjects.
- VE is interpolated from its own finite rows (the spec only drops NaN HR/Speed/VO2); a test with no finite
  VE is skipped and reported.
- Duplicate time stamps are dropped (first kept) so the linear interpolation is well defined.
- With --n-subjects N the same every-k-th rule is used to pick N tests from the qualifying list.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from scipy.ndimage import uniform_filter1d

from common import VARIANTS, pair_test, params_row, stable_seed, variant_params

DATASET = "treadmill"
FS = 0.5
MAX_WINDOW_S = 240
BASE = "https://physionet-open.s3.amazonaws.com/treadmill-exercise-cardioresp/1.0.1/"
FILES = ("test_measure.csv", "subject-info.csv")
PAIRS = {"speed_hr": 1, "speed_vo2": 2}  # name -> index of y in system [speed, hr, vo2, ve]
HERE = os.path.dirname(os.path.abspath(__file__))


def download(data_dir: str) -> dict:
    os.makedirs(data_dir, exist_ok=True)
    prov = {}
    for f in FILES:
        dst = os.path.join(data_dir, f)
        if not os.path.exists(dst):
            subprocess.run(["curl", "-sSf", "--retry", "3", "-o", dst, BASE + f], check=True)
        prov[f] = os.path.getsize(dst)
    return prov


def select_tests(df: pd.DataFrame, n: int | None) -> tuple[list[str], dict]:
    tcol = next(c for c in df.columns if c.lower().startswith("time"))
    scol = next(c for c in df.columns if c.lower().startswith("speed"))
    keep = []
    for tid, g in df.groupby("ID_test", sort=True):
        t = g[tcol].to_numpy(float)
        if len(g) < 2 or t.max() - t.min() < 600:
            continue
        if g["HR"].notna().mean() < 0.95:
            continue
        if g[scol].dropna().nunique() < 3:
            continue
        keep.append((str(tid), g["ID"].iloc[0]))
    n_qual = len(keep)
    seen, uniq = set(), []
    for tid, sid in sorted(keep):  # string sort of ID_test
        if sid not in seen:
            seen.add(sid)
            uniq.append(tid)
    target = 200 if n is None else n
    if len(uniq) > target:
        idx = np.unique(np.floor(np.arange(target) * len(uniq) / target).astype(int))
        uniq = [uniq[i] for i in idx]
    return uniq, {"n_tests": int(df["ID_test"].nunique()), "n_qualifying": n_qual, "n_selected": len(uniq)}


def preprocess(g: pd.DataFrame) -> list[np.ndarray] | None:
    tcol = next(c for c in g.columns if c.lower().startswith("time"))
    scol = next(c for c in g.columns if c.lower().startswith("speed"))
    g = g.sort_values(tcol).drop_duplicates(tcol)
    g = g[g["HR"].notna() & g[scol].notna() & g["VO2"].notna()]
    if len(g) < 10:
        return None
    t = g[tcol].to_numpy(float)
    grid = np.arange(t[0], t[-1] + 1e-9, 1 / FS)
    sig = []
    for col in (scol, "HR", "VO2", "VE"):
        ok = g[col].notna().to_numpy()
        if ok.sum() < 10:
            return None
        x = np.interp(grid, t[ok], g[col].to_numpy(float)[ok])
        sig.append(x - uniform_filter1d(x, 60, mode="nearest"))
    if len(grid) < 300 or not all(np.all(np.isfinite(s)) for s in sig):
        return None
    return sig


def run_subject(args):
    tid, sig, n_null = args
    pv = variant_params(sig, FS, MAX_WINDOW_S, seed=stable_seed(DATASET, tid))
    prow, rows = [], []
    for v in VARIANTS:
        p, info = pv[v]
        prow.append({"subject": tid, "variant": v, **params_row(p, FS),
                     "calibration_failed": bool(info.get("calibration_failed", False)),
                     "window_clipped": bool(info.get("window_clipped", False))})
        for name, j in PAIRS.items():
            r = pair_test(sig[0], sig[j], p, n_null=n_null, seed=stable_seed(DATASET, tid, name))
            rows.append({"subject": tid, "variant": v, "test": name, **r, "lag_s": -r["lag_samples"] / FS})
    return tid, prow, rows


def summarise(rows: pd.DataFrame, params: pd.DataFrame) -> pd.DataFrame:
    out = []
    for v in VARIANTS:
        pm = params[params.variant == v]
        for name in PAIRS:
            r = rows[(rows.variant == v) & (rows.test == name)]
            sig = r[r.p < 0.05]
            det = sig[sig.lag_s > 0]
            wrong = sig[sig.lag_s <= 0]  # lag NaN (no stable window) cannot be p<0.05 in practice; counts as <=0 only if finite
            wrong = wrong[np.isfinite(wrong.lag_s)]
            out.append({
                "variant": v, "test": name, "n": len(r),
                "detect": len(det) / len(r) if len(r) else np.nan,
                "n_sig": len(sig), "wrong_dir": len(wrong) / len(sig) if len(sig) else np.nan,
                "lag_s_median": det.lag_s.median() if len(det) else np.nan,
                "lag_s_q25": det.lag_s.quantile(.25) if len(det) else np.nan,
                "lag_s_q75": det.lag_s.quantile(.75) if len(det) else np.nan,
                "window_s_med": pm.window_s.median(), "max_lag_s_med": pm.max_lag_s.median(),
                "tolerance_s_med": pm.tolerance_s.median(), "n_windows_med": r.n_windows.median(),
                "calibration_failed": int(pm.calibration_failed.sum()),
            })
    return pd.DataFrame(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/home/user/data_gt/treadmill")
    ap.add_argument("--out", default=os.path.join(HERE, "results"))
    ap.add_argument("--n-subjects", type=int, default=None)
    ap.add_argument("--n-null", type=int, default=100)
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--smoke", action="store_true")
    a = ap.parse_args()
    if a.smoke:
        a.n_subjects, a.n_null, a.out = 3, 20, os.path.join(HERE, "results", "smoke")
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    prov = download(a.data_dir)
    df = pd.read_csv(os.path.join(a.data_dir, "test_measure.csv"))
    ids, sel = select_tests(df, a.n_subjects)
    print(f"tests {sel['n_tests']}, qualifying {sel['n_qualifying']}, selected {sel['n_selected']}", flush=True)
    groups = {tid: g for tid, g in df[df.ID_test.astype(str).isin(ids)].groupby(df.ID_test.astype(str))}
    jobs, skipped = [], []
    for tid in ids:
        s = preprocess(groups[tid])
        if s is None:
            skipped.append(tid)
        else:
            jobs.append((tid, s, a.n_null))
    if skipped:
        print("skipped (preprocessing failed):", skipped)
    prows, rows = [], []
    with ProcessPoolExecutor(a.jobs) as ex:
        for tid, pr, r in ex.map(run_subject, jobs):
            prows += pr
            rows += r
            print("done", tid, f"{time.time() - t0:.0f}s", flush=True)
    params, rows = pd.DataFrame(prows), pd.DataFrame(rows)
    summ = summarise(rows, params)
    params.to_csv(os.path.join(a.out, f"{DATASET}_params.csv"), index=False)
    rows.to_csv(os.path.join(a.out, f"{DATASET}_rows.csv"), index=False)
    summ.to_csv(os.path.join(a.out, f"{DATASET}_summary.csv"), index=False)
    wall = time.time() - t0
    with open(os.path.join(a.out, f"{DATASET}_meta.json"), "w") as f:
        json.dump({"dataset": DATASET, "files_bytes": prov, "selection": sel, "subjects": [j[0] for j in jobs],
                   "skipped": skipped, "n_null": a.n_null, "fs": FS, "max_window_s": MAX_WINDOW_S,
                   "jobs": a.jobs, "wall_s": wall}, f, indent=1)
    with pd.option_context("display.width", 250, "display.max_columns", 30, "display.float_format", "{:.3f}".format):
        print(summ.to_string(index=False))
    print(f"wall {wall:.1f}s")


if __name__ == "__main__":
    main()
