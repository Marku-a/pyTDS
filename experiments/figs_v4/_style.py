"""Shared style for the figs_v4 figures (copied from experiments/atds_real/plot_real_summary.py / plot_real_traces.py)."""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
OUTDIR = os.path.join(ROOT, "paper", "simple", "figures")
INK, MUTED, GREY = "#262624", "#8a8a85", "#c9c8c2"
COL = {"published": "#2a78d6", "atds": "#eb6834"}
WIDTH, DPI = 6.8, 200
plt.rcParams.update({
    "font.family": "serif", "font.size": 8, "axes.titlesize": 8.5,
    "legend.fontsize": 7.2, "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": MUTED,
    "axes.grid": True, "grid.color": "#e6e5e0", "grid.linewidth": 0.6, "legend.frameon": False,
    "xtick.color": MUTED, "ytick.color": MUTED, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})


def save(fig, name):
    os.makedirs(OUTDIR, exist_ok=True)
    out = os.path.join(OUTDIR, name)
    fig.savefig(out, dpi=DPI)
    print("wrote", out)
