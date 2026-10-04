"""Control: is the bigger sleep-stage contrast of 300-s windows just window SELECTION?

A 300-s window counts for a stage only if >= 80 % of it is that stage, so long-window
variants only see consolidated stage episodes. Here the published 60-s setting is scored
twice per night: (a) the normal 80 % rule on its own 60 s, (b) only windows whose
surrounding 300 s also pass the 80 % rule. If (b) closes the gap to v2_calibrated, the
gain is selection; if not, it is the method. Writes results/sleep_selection_control.csv.
"""

import glob
import os
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations

import numpy as np
import pandas as pd

from common import published_params, window_labels
from sleep_edf import STAGES, build_system, stage_scores, window_stage

DATA = "/home/user/data_gt/sleep-edfx"


def run(path):
    z = np.load(path)
    sigs, st, _ = build_system({"bp": z["bp"], "stage": z["stage"]})
    p = published_params(1.0)
    out = []
    for i, j in combinations(range(len(sigs)), 2):
        c, lbl, _ = window_labels(sigs[i], sigs[j], p)
        for rule, W in (("own_60s", p.window), ("within_300s", 300)):
            out.append([os.path.basename(path)[:6], rule, *stage_scores(lbl, window_stage(c, st, W))])
    return out


if __name__ == "__main__":
    paths = sorted(glob.glob(f"{DATA}/SC4*.npz"))
    with ProcessPoolExecutor(4) as ex:
        rows = [r for rr in ex.map(run, paths) for r in rr]
    df = pd.DataFrame(rows, columns=["subject", "rule", *STAGES])
    night = df.groupby(["subject", "rule"])[list(STAGES)].mean().reset_index()
    night.to_csv("results/sleep_selection_control.csv", index=False)
    for rule, g in night.groupby("rule"):
        for s in ("LS", "W", "REM"):
            d = (g[s] - g["DS"]).dropna()
            print(f"{rule:12s} {s}>DS  n={len(d)}  mean={d.mean():6.2f}  dz={d.mean() / d.std(ddof=1):.2f}")
