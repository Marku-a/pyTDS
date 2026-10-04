"""Bootstrap 95 % CIs for each variant's ground-truth effect size, and for its
difference from the published fixed setting (resampling subjects, paired).

Fantasia: Cohen's d, young vs elderly, on within-subject excess (score - chance95).
Sleep-EDF: Cohen's dz of night-level excess, LS vs DS, W vs DS, REM vs DS.
Writes results/effects_bootstrap.csv.
"""

import numpy as np
import pandas as pd

B, REF = 5000, "published"
rng = np.random.default_rng(0)


def d_unpaired(a, b):
    sp = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    return (a.mean() - b.mean()) / sp


def dz(x):
    return x.mean() / x.std(ddof=1)


rows = []
f = pd.read_csv("results/fantasia_rows.csv")
f = f[f.test == "within"].pivot(index="subject", columns="variant", values="excess")
young = np.array(["y" in s for s in f.index])
fy, fo = f[young], f[~young]
for v in f.columns:
    est = d_unpaired(fy[v].values, fo[v].values)
    bs, bd = [], []
    for _ in range(B):
        a = fy.sample(len(fy), replace=True, random_state=rng.integers(1 << 31))
        b = fo.sample(len(fo), replace=True, random_state=rng.integers(1 << 31))
        x = d_unpaired(a[v].values, b[v].values)
        bs.append(x)
        bd.append(x - d_unpaired(a[REF].values, b[REF].values))
    rows.append({"dataset": "fantasia", "effect": "young>old d", "variant": v, "estimate": est,
                 "ci_lo": np.percentile(bs, 2.5), "ci_hi": np.percentile(bs, 97.5),
                 "diff_vs_published": est - d_unpaired(fy[REF].values, fo[REF].values),
                 "diff_ci_lo": np.percentile(bd, 2.5), "diff_ci_hi": np.percentile(bd, 97.5)})

s = pd.read_csv("results/sleep_edf_rows.csv")
night = s.groupby(["subject", "variant", "stage"]).excess.mean().unstack("stage")
for st in ("LS", "W", "REM"):
    diff = (night[st] - night["DS"]).unstack("variant").dropna()  # nights usable in every variant
    for v in diff.columns:
        x = diff[v].values
        bs, bd = [], []
        for _ in range(B):
            i = rng.integers(0, len(x), len(x))
            bs.append(dz(x[i]))
            bd.append(dz(x[i]) - dz(diff[REF].values[i]))
        rows.append({"dataset": "sleep_edf", "effect": f"{st}>DS dz", "variant": v, "n": len(x),
                     "estimate": dz(x), "ci_lo": np.percentile(bs, 2.5), "ci_hi": np.percentile(bs, 97.5),
                     "diff_vs_published": dz(x) - dz(diff[REF].values),
                     "diff_ci_lo": np.percentile(bd, 2.5), "diff_ci_hi": np.percentile(bd, 97.5)})

out = pd.DataFrame(rows)
out.to_csv("results/effects_bootstrap.csv", index=False)
pd.set_option("display.width", 200)
print(out.round(2).to_string(index=False))
