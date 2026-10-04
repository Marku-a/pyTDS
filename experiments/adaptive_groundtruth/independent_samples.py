"""Mechanism check: how many independent samples does one TDS window hold?

n_eff = window / B, with B the Bartlett variance-inflation factor of the slowest
signal in the system (adaptive.bartlett_factor; the same quantity the v2 rule sizes
its window with). Reported for the published 60-s window and for the v2 window, per
dataset. Writes results/independent_samples.csv.
"""

import glob
import os

import numpy as np
import pandas as pd

import fantasia
import sleep_edf
import treadmill
from adaptive import bartlett_factor
from common import published_params

rows = []
for p in sorted(glob.glob("/home/user/data_gt/sleep-edfx/SC4*.npz")):
    z = np.load(p)
    sigs, _, _ = sleep_edf.build_system({"bp": z["bp"], "stage": z["stage"]})
    rows.append(("sleep_edf", os.path.basename(p)[:6], 1.0, max(bartlett_factor(s) for s in sigs)))

beats = fantasia.beat_symbols()
for rec in sorted({os.path.basename(p)[:5] for p in glob.glob("/home/user/data_gt/fantasia/f*.hea")}):
    data, _, _ = fantasia.preprocess("/home/user/data_gt/fantasia", rec, beats)
    if data is not None:
        rows.append(("fantasia", rec, 4.0, max(bartlett_factor(s) for s in data)))

df = pd.read_csv("/home/user/data_gt/treadmill/test_measure.csv")
ids, _ = treadmill.select_tests(df, None)
for tid in ids:
    sig = treadmill.preprocess(df[df["ID_test"] == tid])
    if sig is not None:
        rows.append(("treadmill", tid, 0.5, max(bartlett_factor(s) for s in sig)))

out = pd.DataFrame(rows, columns=["dataset", "subject", "fs", "bartlett_slowest"])
prm = pd.concat([pd.read_csv(f"results/{d}_params.csv").assign(dataset=d) for d in ("sleep_edf", "fantasia", "treadmill")])
v2 = prm[prm.variant == "v2_calibrated"].groupby("dataset").window_s.median()
summ = []
for d, g in out.groupby("dataset"):
    fs = g.fs.iloc[0]
    B = g.bartlett_slowest.median()
    pub = published_params(fs).window
    summ.append({"dataset": d, "n_subjects": len(g), "bartlett_median": B,
                 "published_window_s": pub / fs, "n_eff_published": pub / B,
                 "v2_window_s": v2[d], "n_eff_v2": v2[d] * fs / B})
summ = pd.DataFrame(summ)
summ.to_csv("results/independent_samples.csv", index=False)
print(summ.round(1).to_string(index=False))
