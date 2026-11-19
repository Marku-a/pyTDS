"""Real-data figure: delay found in each window (tau0) under the published setting and aTDS.

    python experiments/atds_real/plot_real_traces.py <fields_data_dir> <bidmc_csv_dir> <bidmc_record> <out_png>

Uses this repo's pyTDS only (through atds21). Stable windows are coloured, unstable ones grey;
the dashed line is the expected delay (tau0 = -delay when the first signal leads).
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import atds21 as A  # noqa: E402
import bidmc  # noqa: E402
import fields as F  # noqa: E402

INK, MUTED, GREY = "#262624", "#8a8a85", "#c9c8c2"
COL = {"published": "#2a78d6", "atds": "#eb6834"}
plt.rcParams.update({
    "font.family": "serif", "font.size": 8, "axes.titlesize": 8.5,
    "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": MUTED,
    "xtick.color": MUTED, "ytick.color": MUTED, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})


def eeg_benchmark():
    root = os.path.abspath(os.path.join(HERE, "..", ".."))
    x = np.loadtxt(os.path.join(root, "data", "data-delta.txt"))[:, 1]
    y = np.loadtxt(os.path.join(root, "data", "data-sigma.txt"))[:, 1]
    h = len(x) // 2  # parameters set on the first half of the night, scored on the second (as in the paper)
    return dict(s1=x[h:], s2=y[h:], train=(x[:h], y[:h]), dt=1.0, names=("EEG delta", "EEG sigma"), expected=None,
                ref=0.0, unit=("min", 60.0), tunit=("min", 60.0), ylim=8,
                title="Sleep EEG power, delta vs sigma (1 Hz; TDS benchmark, second half of the night)")


def rows(fd, bdir, rec):
    F.FD = fd
    h = F.d_hydro_hourly(); h.update(unit=("h", 3600.0), tunit=("days", 86400.0), ylim=40,
                                     title="Rainfall to river flow (hourly; expected delay 8 h)")
    e = F.d_energy(); e.update(unit=("h", 3600.0), tunit=("days", 86400.0), ylim=12,
                               title="Temperature to electricity demand (30 min, three summers joined; expected delay 1 h)")
    b = bidmc.pair(bdir, rec)
    b.update(expected=None, ref=bidmc.xcorr_delay(b["s1"], b["s2"], bidmc.reference(*bidmc.load(bdir, rec))[1]),
             unit=("s", 1.0), tunit=("s", 1.0), ylim=0.45,
             title=f"ECG to finger pulse, ICU patient {rec} of 53 (125 Hz; median aTDS score)")
    return [eeg_benchmark(), b, h, e]


def panel(ax, d, name, p):
    x, y = np.asarray(d["s1"], float), np.asarray(d["s2"], float)
    tau = A.time_delay_interaction(x, y, p)[0]
    lab = np.asarray(A.stable_label(tau, p)).astype(bool)
    score = float(A.tds_score(lab))
    u, us = d["unit"]
    tu, tus = d["tunit"]
    t = (np.arange(len(tau)) * p.overlap + p.window / 2) * d["dt"] / tus
    yd = np.asarray(tau) * d["dt"] / us
    ax.scatter(t[~lab], yd[~lab], s=2.5, color=GREY, lw=0, label="unstable window")
    ax.scatter(t[lab], yd[lab], s=3.5, color=COL[name], lw=0, label="stable window")
    lag = p.max_lag * d["dt"] / us
    ax.axhspan(-lag, lag, color=COL[name], alpha=0.06, lw=0)
    if d.get("expected") is not None:
        ax.axhline(-d["expected"][1] * d["dt"] / us, color=INK, ls="--", lw=0.8)
    elif d.get("ref") is not None:
        ax.axhline(d["ref"] * d["dt"] / us, color=INK, ls=":", lw=0.8)
    ax.set_ylim(-d["ylim"], d["ylim"])
    ax.set_xlabel(f"time ({tu})")
    lbl = "published" if name == "published" else "aTDS"
    ax.set_title(f"{lbl}: TDS {score:.0f} %", loc="left", color=COL[name])
    ax.set_ylabel(f"delay ({u})")
    print("  ", name, f"score {score:.1f}", flush=True)
    return score


def main(fd, bdir, rec, out):
    R = rows(fd, bdir, rec)
    fig, axes = plt.subplots(len(R), 2, figsize=(6.8, 1.9 * len(R)), sharey="row",
                             gridspec_kw=dict(hspace=1.05, wspace=0.08))
    for k, d in enumerate(R):
        x, y = np.asarray(d["s1"], float), np.asarray(d["s2"], float)
        tx, ty = d.get("train", (x, y))
        pa = A.window_params([tx, ty]); pa, _ = A.calib(pa, [(tx, ty)])
        print(d["title"], "aTDS", pa, flush=True)
        for c, (name, p) in enumerate((("published", A.PUBLISHED), ("atds", pa))):
            panel(axes[k, c], d, name, p)
        axes[k, 0].annotate(d["title"], xy=(0, 1.32), xycoords="axes fraction", fontsize=8.5, color=INK,
                            fontweight="bold")
        axes[k, 1].set_ylabel("")
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    fig.savefig(out, dpi=200)
    print("wrote", out)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4])
