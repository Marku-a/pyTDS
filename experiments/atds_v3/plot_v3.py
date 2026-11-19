"""Run-3 figure: (a) slow-down AUC for all methods incl. full-record baselines; (b) regime switching.

    python experiments/atds_v3/plot_v3.py <out_png>
"""
import ast
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

RES = Path(__file__).resolve().parent / "results"
INK, MUTED = "#262624", "#8a8a85"
plt.rcParams.update({
    "font.family": "serif", "font.size": 8.5, "axes.titlesize": 9, "legend.fontsize": 7.2,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
    "grid.color": "#e6e5e0", "grid.linewidth": 0.6, "axes.edgecolor": MUTED,
    "xtick.color": MUTED, "ytick.color": MUTED, "legend.frameon": False,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})
SER = [  # method, label, colour, marker, linestyle
    ("published", "published 60/±30/±1", "#2a78d6", "o", "-"),
    ("rescaled", "rescaled by true r", "#1baf7a", "^", "-"),
    ("rescaled_est", "rescaled by estimated r", "#1baf7a", "v", ":"),
    ("atds", "aTDS", "#eb6834", "s", "-"),
    ("xcorr_full", "full-record xcorr / GCC-PHAT", "#6b6b66", "D", "--"),
]


def ci(v):
    return ast.literal_eval(v) if isinstance(v, str) else v


e1 = pd.read_csv(RES / "e1_baselines.csv")
e3 = pd.read_csv(RES / "e3_synthetic.csv")
fig, (a, b) = plt.subplots(1, 2, figsize=(6.6, 2.7), gridspec_kw={"width_ratios": [1.25, 1]})
for m, lab, c, mk, ls in SER:
    d = e1[e1.method == m].sort_values("r")
    lo = np.array([ci(x)[0] for x in d.auc_ci]); hi = np.array([ci(x)[1] for x in d.auc_ci])
    a.errorbar(d.r, d.auc, yerr=[d.auc - lo, hi - d.auc], fmt=mk, ls=ls, color=c, capsize=2, lw=1.3, ms=4, label=lab)
a.set_xscale("log", base=2)
a.set_xticks([1, 2, 4, 8, 16], ["1", "2", "4", "8", "16"])
a.set_xlabel("slow-down factor r")
a.set_ylabel("AUC")
a.set_ylim(0.4, 1.04)
a.set_title("(a) stationary benchmark: all methods", loc="left")
a.legend(loc="upper center", bbox_to_anchor=(0.5, -0.3), ncol=2)

rs = e3[e3.benchmark == "regime_switching"].set_index("method")
names = [("published", "published"), ("atds", "aTDS"), ("atds_per_regime", "aTDS,\nper regime")]
x = np.arange(len(names)); w = 0.36
cA = [rs.loc[k, "regimeA_coupled_mean"] for k, _ in names]
cB = [rs.loc[k, "false_stability_B_coupled_mean"] for k, _ in names]
b.bar(x - w / 2, cA, w, color="#eb6834", label="coupled regime: stable windows")
b.bar(x + w / 2, cB, w, color="#2a78d6", label="uncoupled regime: false stability")
for xi, v in zip(x - w / 2, cA):
    b.text(xi, v + 1.5, f"{v:.0f}", ha="center", fontsize=7, color=INK)
for xi, v in zip(x + w / 2, cB):
    b.text(xi, v + 1.5, f"{v:.1f}", ha="center", fontsize=7, color=INK)
b.axhline(5, color=INK, ls=":", lw=0.8, label="5 % false-stability budget")
b.set_xticks(x, [n for _, n in names])
b.set_ylabel("score (%)")
b.set_ylim(0, 112)
b.set_title("(b) coupling that switches with the regime", loc="left")
b.legend(loc="upper center", bbox_to_anchor=(0.5, -0.3), ncol=1)
fig.tight_layout()
out = Path(sys.argv[1]); out.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(out, dpi=200)
print("wrote", out)
