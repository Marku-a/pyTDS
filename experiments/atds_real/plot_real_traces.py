"""Real-data figures: delay found in each window (tau0) under fixed TDS and aTDS.

    python experiments/atds_real/plot_real_traces.py <fields_data_dir> <bidmc_csv_dir> <bidmc_record> <out_dir>

Writes fig_physio_traces.png (EEG benchmark, ECG -> finger pulse; fixed TDS | aTDS), fig_fields_traces.png (aTDS only,
five field records) and fig_fields_summary.png (aTDS delay vs expected delay). Field panels use the aTDS parameters stored in
the result files (asserted to reproduce their score and delay).

Uses this repo's pyTDS only (through atds21). Stable windows are coloured, unstable ones grey;
the dashed line is the expected delay (tau0 = -delay when the first signal leads).
"""
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import atds21 as A  # noqa: E402
import bidmc  # noqa: E402
import fields as F  # noqa: E402
from pyTDS.params import TDSParams  # noqa: E402

INK, MUTED, GREY = "#262624", "#8a8a85", "#c9c8c2"
COL = {"published": "#2a78d6", "atds": "#eb6834"}
plt.rcParams.update({
    "font.family": "serif", "font.size": 8, "axes.titlesize": 8.5,
    "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": MUTED,
    "legend.frameon": False,
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


def rows(bdir, rec):
    b = bidmc.pair(bdir, rec)
    b.update(expected=None, ref=bidmc.xcorr_delay(b["s1"], b["s2"], bidmc.reference(*bidmc.load(bdir, rec))[1]),
             unit=("s", 1.0), tunit=("s", 1.0), ylim=0.45,
             title=f"ECG to finger pulse, ICU patient {rec} of 53 (125 Hz; median aTDS score)")
    return [eeg_benchmark(), b]


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
    lbl = "fixed TDS" if name == "published" else "aTDS"
    ax.set_title(f"{lbl}: TDS {score:.0f} %", loc="left", color=COL[name])
    ax.set_ylabel(f"delay ({u})")
    print("  ", name, f"score {score:.1f}", flush=True)
    return score


def physio(bdir, rec, out):
    R = rows(bdir, rec)
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
    fig.savefig(out, dpi=200)
    print("wrote", out)


# ---- field records, aTDS only; parameters and values from the result files
RES = os.path.join(HERE, "results")
V3RES = os.path.join(HERE, "..", "atds_v3", "results")
H, D, W = 3600.0, 86400.0, 7 * 86400.0


def _daily(fd, code):
    df = pd.read_csv(os.path.join(fd, f"hydrology_airgr_{code}_daily.csv"))[["precip_mm", "q_mm"]].interpolate(limit_direction="both")
    return dict(s1=df.precip_mm.values, s2=df.q_mm.values, dt=D)   # as e5_hydrology.load_daily


def field_records(fd):
    F.FD = fd
    f21 = json.load(open(os.path.join(RES, "fields_v21.json")))
    e5 = json.load(open(os.path.join(V3RES, "e5_hydrology.json")))
    # name, loader, result entry, unit (label, seconds), expected range in unit (lo, hi) or point, expected text, ylim (in unit)
    return [
        dict(name="Rainfall to river, hourly (L0123003)", d=F.d_hydro_hourly(), r=f21["hydrology_hourly"]["settings"]["atds"],
             unit=("h", H), exp=(8, 8), exp_txt="8 h", ylim=(-10, 50)),
        dict(name="Rainfall to river, daily (L0123001)", d=_daily(fd, "L0123001"), r=e5["L0123001_daily"]["settings"]["atds"],
             unit=("days", D), exp=(0, 3), exp_txt="0-3 days", ylim=(-12, 40)),
        dict(name="Rainfall to river, daily (L0123002)", d=_daily(fd, "L0123002"), r=e5["L0123002_daily"]["settings"]["atds"],
             unit=("days", D), exp=(1, 1), exp_txt="1 day, xcorr peak", ylim=(-250, 250)),
        dict(name="Temperature to electricity demand", d=F.d_energy(), r=f21["energy_vic_summer"]["settings"]["atds"],
             unit=("h", H), exp=(1, 1), exp_txt="1 h", ylim=(-10, 15)),
        dict(name="Crude oil to gasoline, weekly changes", d=F.d_gasoline(), r=f21["gasoline_weekly_diff"]["settings"]["atds"],
             unit=("weeks", W), exp=(1, 2), exp_txt="1-2 weeks", ylim=(-10, 12)),
    ]


def fields_traces(fd, out):
    recs = field_records(fd)
    fig, axes = plt.subplots(len(recs), 1, figsize=(6.8, 1.75 * len(recs)), gridspec_kw=dict(hspace=0.85))
    for ax, R in zip(axes, recs):
        r, d = R["r"], R["d"]
        ps = r["params_samples"]
        p = TDSParams(window=ps["window"], overlap=ps["step"], max_lag=ps["max_lag"], tolerance=ps["tolerance"])
        print(R["name"], "aTDS params", ps, flush=True)
        x, y = np.asarray(d["s1"], float), np.asarray(d["s2"], float)
        tau = np.asarray(A.time_delay_interaction(x, y, p)[0])
        lab = np.asarray(A.stable_label(tau, p)).astype(bool)
        score = float(A.tds_score(lab)); dl = A.delay_of(tau, lab)
        assert abs(score - r["score"]) < 1e-9 and dl == r["delay_samples"] and len(tau) == r["n_windows"], (R["name"], score, dl)
        u, us = R["unit"]
        t = (np.arange(len(tau)) * p.overlap + p.window / 2) * d["dt"] / (365.25 * 86400)
        yd = -tau * d["dt"] / us                      # delay of the second signal after the first (tau = -delay)
        lo, hi = R["exp"]
        if lo == hi:
            ax.axhline(lo, color=INK, ls="--", lw=0.8)
        else:
            ax.axhspan(lo, hi, color=INK, alpha=0.15, lw=0)
        ax.scatter(t[~lab], yd[~lab], s=6, color=GREY, lw=0, label="unstable window")
        ax.scatter(t[lab], yd[lab], s=8, color=COL["atds"], lw=0, label="stable window")
        ax.set_ylim(*R["ylim"])
        found = abs(dl) * d["dt"] / us
        ax.set_title(f"{R['name']}: aTDS delay {found:.3g} {u} (expected {R['exp_txt']}), TDS {r['score']:.0f} %, "
                     f"p = {r['p']:.3f}", loc="left", fontsize=7.6)
        ax.set_ylabel(f"delay ({u})")
        ax.set_xlabel("time in the record (years)")
        nout = int(np.sum((yd < R["ylim"][0]) | (yd > R["ylim"][1])))
        print(f"  score {score:.2f} delay_samples {dl} found {found:.3g} {u} p {r['p']:.4f} windows {len(tau)} outside view {nout}", flush=True)
    axes[0].legend(loc="upper right", ncol=2, fontsize=6.8)
    fig.savefig(out, dpi=200)
    print("wrote", out)
    return recs


def fields_summary(recs, out):
    short = ["hourly\nL0123003", "daily\nL0123001", "daily\nL0123002", "temperature\nto electricity", "oil to\ngasoline"]
    fig, ax = plt.subplots(figsize=(5.2, 3.0))
    for i, R in enumerate(recs):
        d, r = R["d"], R["r"]
        u, us = R["unit"]
        lo, hi = (v * us / H for v in R["exp"])
        found = abs(r["delay_samples"]) * d["dt"] / H
        sig = r["p"] <= 0.05
        print(f"C {R['name']}: expected {R['exp_txt']} = [{lo:g}, {hi:g}] h; found {found:g} h; p {r['p']:.4f}; significant {sig}")
        if hi > lo:
            ax.plot([i, i], [lo, hi], color=MUTED, lw=6, alpha=0.5, solid_capstyle="butt", label="expected delay" if i == 0 else None)
        else:
            ax.plot([i - 0.22, i + 0.22], [lo, lo], color=MUTED, lw=3, alpha=0.8, solid_capstyle="butt",
                    label="expected delay" if i == 0 else None)
        ax.scatter([i], [found], s=36, marker="o", facecolor=COL["atds"] if sig else "white", edgecolor=COL["atds"], lw=1.4, zorder=3,
                   label=("aTDS delay found" if i == 0 else None) if sig else ("aTDS, not significant (p > 0.05)" if i == 2 else None))
    ax.set_yscale("symlog", linthresh=1)
    ax.set_ylim(-0.3, 6000)
    ax.set_xticks(range(len(recs)), short, fontsize=7)
    ax.set_ylabel("delay (hours, symlog)")
    ax.set_xlim(-0.6, len(recs) - 0.4)
    ax.legend(loc="upper left", fontsize=6.8)
    ax.set_title("Delay found by aTDS against the expected delay", loc="left")
    fig.savefig(out, dpi=200)
    print("wrote", out)


def main(fd, bdir, rec, outdir):
    os.makedirs(outdir, exist_ok=True)
    physio(bdir, rec, os.path.join(outdir, "fig_physio_traces.png"))
    recs = fields_traces(fd, os.path.join(outdir, "fig_fields_traces.png"))
    fields_summary(recs, os.path.join(outdir, "fig_fields_summary.png"))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4])
