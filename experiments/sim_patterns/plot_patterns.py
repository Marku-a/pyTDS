"""Figure: tau0 per window on simulated delay patterns, fixed TDS vs aTDS (first coupled pair, seed 1000*k+100).

    python experiments/sim_patterns/plot_patterns.py   (run patterns.py first)
"""
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
import patterns as P  # noqa: E402

A = P.A
from dataclasses import replace  # noqa: E402

INK, MUTED, GREY = "#262624", "#8a8a85", "#c9c8c2"
COL = {"fixed": "#2a78d6", "atds": "#eb6834"}
plt.rcParams.update({
    "font.family": "serif", "font.size": 8, "axes.titlesize": 8.5,
    "legend.fontsize": 7.2, "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": MUTED,
    "axes.grid": True, "grid.color": "#e6e5e0", "grid.linewidth": 0.6, "legend.frameon": False,
    "xtick.color": MUTED, "ytick.color": MUTED, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})
TITLES = {"fixed": "Fixed delay (20)", "step": "Delay that jumps (10 → 30)",
          "drift": "Drifting delay (10 → 30)", "harmonic": "Harmonic delay (20 ± 10)",
          "onoff": "Coupling on/off (delay 20)"}
OUT = os.path.join(ROOT, "paper", "simple", "figures", "fig_delay_patterns.png")


def main():
    res = json.load(open(os.path.join(HERE, "results", "patterns.json")))
    fig, axes = plt.subplots(5, 2, figsize=(6.8, 9.6), sharey="row", gridspec_kw=dict(hspace=0.95, wspace=0.08))
    for k, kind in enumerate(P.PATTERNS, 1):
        x, y = P.make_pair(kind, P.C, P.seeds(k, "coupled", 0))
        a = res[kind]["atds_params"]
        pa = replace(A.PUBLISHED, window=a["window"], overlap=a["overlap"], max_lag=a["max_lag"], tolerance=a["tolerance"])
        d = P.delay(kind)
        for c, (name, p) in enumerate((("fixed", A.PUBLISHED), ("atds", pa))):
            ax = axes[k - 1, c]
            tau, tv, _ = A.time_delay_interaction(x, y, p)
            lab = np.asarray(A.stable_label(tau, p)).astype(bool)
            score = float(A.tds_score(lab))
            t = tv / 1000.0
            if kind == "onoff":
                for b in range(0, P.N // P.BLOCK, 2):
                    ax.axvspan(b * P.BLOCK / 1000, (b + 1) * P.BLOCK / 1000, color="#ececea", lw=0, zorder=0)
            ax.scatter(t[~lab], tau[~lab], s=2.5, color=GREY, lw=0, label="unstable window")
            ax.scatter(t[lab], tau[lab], s=3.5, color=COL[name], lw=0, label="stable window")
            ax.axhspan(-p.max_lag, p.max_lag, color=COL[name], alpha=0.06, lw=0)
            tt = np.arange(P.N)
            tr = -d.astype(float)
            if kind == "onoff":
                cpl = (tt // P.BLOCK) % 2 == 0
                for b in range(0, P.N // P.BLOCK, 2):
                    sl = slice(b * P.BLOCK, (b + 1) * P.BLOCK)
                    ax.plot(tt[sl] / 1000, tr[sl], color=INK, ls="--", lw=0.9, label="true delay" if b == 0 else None)
            else:
                ax.plot(tt / 1000, tr, color=INK, ls="--", lw=0.9, label="true delay")
            ax.set_ylim(-60, 15)
            ax.set_xlim(0, P.N / 1000)
            ax.set_xlabel("time (10$^3$ samples)")
            ax.set_ylabel("delay (samples)" if c == 0 else "")
            lbl = "fixed TDS" if name == "fixed" else "aTDS"
            ax.set_title(f"{lbl}: TDS {score:.0f} %", loc="left", color=COL[name])
        axes[k - 1, 0].annotate(TITLES[kind], xy=(0, 1.3), xycoords="axes fraction", fontsize=8.5, color=INK, fontweight="bold")
    h, l = axes[1, 1].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.02), markerscale=2.5)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, dpi=200)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
