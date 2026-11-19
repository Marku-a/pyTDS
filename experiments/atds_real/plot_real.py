"""Figure of the NCOM cohort results (per-stage network strength), from results/ only.

    python experiments/atds_real/plot_real.py <out_png>
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

RES = Path(__file__).resolve().parent / "results"
C_PUB, C_ATDS, INK, MUTED = "#2a78d6", "#eb6834", "#262624", "#8a8a85"
plt.rcParams.update({
    "font.family": "serif", "font.size": 8.5, "axes.titlesize": 9, "legend.fontsize": 7.5,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
    "grid.color": "#e6e5e0", "grid.linewidth": 0.6, "axes.edgecolor": MUTED,
    "xtick.color": MUTED, "ytick.color": MUTED, "legend.frameon": False,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})

d = pd.read_csv(RES / "ncom_subject_stage.csv")
stages = ["W", "LS", "REM", "DS"]
labels = {"W": "wake", "LS": "light", "REM": "REM", "DS": "deep"}
settings = [("published", C_PUB, "o", "fixed TDS 60/±30/±1"),
            ("atds_cap300", C_ATDS, "s", "aTDS, window capped at 300 s")]
fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.6))
for ax, col, title in [(axes[0], "mean_link_tds", "(a) mean TDS over 45 links (%)"),
                       (axes[1], "n_significant_links", "(b) significant links (of 45)")]:
    for k, (s, c, m, lab) in enumerate(settings):
        off = -0.15 if k == 0 else 0.15
        for i, st in enumerate(stages):
            v = d[(d.setting == s) & (d.stage == st)][col].dropna().to_numpy()
            x = i + off + np.random.default_rng(i + 10 * k).uniform(-0.06, 0.06, len(v))
            ax.scatter(x, v, s=6, color=c, alpha=0.35, linewidths=0)
            se = v.std(ddof=1) / np.sqrt(len(v))
            ax.errorbar(i + off, v.mean(), yerr=1.96 * se, fmt=m, color=c, ms=5, capsize=2.5, lw=1.4,
                        label=lab if i == 0 else None)
    ax.set_xticks(range(4), [labels[s] for s in stages])
    ax.set_title(title, loc="left")
axes[0].set_ylim(0, None)
h, lab = axes[0].get_legend_handles_labels()
fig.legend(h, lab, loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.07))
fig.tight_layout()
out = Path(sys.argv[1])
out.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(out, dpi=200)
print("wrote", out)
