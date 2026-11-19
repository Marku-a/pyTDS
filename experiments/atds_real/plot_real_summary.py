"""Real-data summary figures (from committed result files only; no raw data needed).

    python experiments/atds_real/plot_real_summary.py <out_dir>

fig_bidmc.png:   ICU ECG -> finger pulse, 53 patients: (a) scores of true pairs vs cross-patient controls,
                 (b) published vs aTDS score per patient.
fig_fields_significance.png: score vs the 95th percentile of 200 circular rotations, per record and setting,
                 with the delay found.
"""
import json
import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
INK, MUTED = "#262624", "#8a8a85"
COL = {"published": "#2a78d6", "atds": "#eb6834"}
LAB = {"published": "fixed TDS", "atds": "aTDS"}
plt.rcParams.update({
    "font.family": "serif", "font.size": 8, "axes.titlesize": 8.5, "legend.fontsize": 7.2,
    "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": MUTED, "axes.grid": True,
    "grid.color": "#e6e5e0", "grid.linewidth": 0.6, "xtick.color": MUTED, "ytick.color": MUTED,
    "legend.frameon": False, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})


def fig_bidmc(out):
    B = json.load(open(HERE / "results" / "bidmc_ecg_ppg.json"))
    R, S = B["records"], B["summary"]
    c = sorted([r for r in R if r["coupled"]], key=lambda r: r["ecg"])
    u = [r for r in R if not r["coupled"]]
    fig, (a, b) = plt.subplots(1, 2, figsize=(6.8, 2.7), gridspec_kw={"width_ratios": [1.25, 1], "wspace": 0.38})
    rng = np.random.default_rng(0)
    for k, m in enumerate(("published", "atds")):
        for j, (rows, kind) in enumerate(((c, "same patient"), (u, "different patients"))):
            x0 = k * 2.4 + j
            v = np.array([r[f"{m}_score"] for r in rows])
            a.scatter(x0 + rng.uniform(-0.22, 0.22, v.size), v, s=7, lw=0,
                      color=COL[m] if j == 0 else MUTED, alpha=0.85)
            a.plot([x0 - 0.3, x0 + 0.3], [np.median(v)] * 2, color=INK, lw=1.2)
        a.plot([k * 2.4 - 0.45, k * 2.4 + 1.45], [S[m]["control_95pct"]] * 2, color=INK, ls=":", lw=0.8)
        a.text(k * 2.4 + 0.5, 104, f"{LAB[m]}\nAUC {S[m]['auc']:.2f}  [{S[m]['auc_ci95'][0]:.2f}, {S[m]['auc_ci95'][1]:.2f}]\n"
               f"{S[m]['n_coupled_above_control95']}/53 above controls",
               ha="center", va="bottom", fontsize=7.2, color=COL[m])
    a.set_xticks([0, 1, 2.4, 3.4], ["same\npatient", "different\npatients"] * 2, fontsize=7)
    a.set_ylabel("TDS score (%)")
    a.set_ylim(-3, 135); a.set_yticks([0, 25, 50, 75, 100])
    a.set_title("(a) ECG vs finger pulse: true pairs and controls", loc="left")
    p = np.array([r["published_score"] for r in c]); q = np.array([r["atds_score"] for r in c])
    b.plot([0, 100], [0, 100], color=MUTED, lw=0.8, ls="--")
    b.scatter(p, q, s=10, color=COL["atds"], lw=0)
    b.set_xlabel("fixed TDS score (%)"); b.set_ylabel("aTDS score (%)")
    b.set_xlim(-3, 103); b.set_ylim(-3, 103)
    b.set_title(f"(b) per patient (aTDS higher in {int(np.sum(q > p))}/53)", loc="left")
    fig.savefig(out, dpi=200)
    print("wrote", out)


def fig_fields(out):
    F = json.load(open(HERE / "results" / "fields_v21.json"))
    E = json.load(open(ROOT / "atds_v3" / "results" / "e5_hydrology.json"))
    recs = [  # label, settings dict, sampling interval (h), expected
        ("River, hourly\n(L0123003)", F["hydrology_hourly"]["settings"], 1.0, "8 h"),
        ("River, daily\n(L0123001)", E["L0123001_daily"]["settings"], 24.0, "0–3 d"),
        ("River, daily\n(L0123002)", E["L0123002_daily"]["settings"], 24.0, "0–3 d"),
        ("Temperature →\nelectricity", F["energy_vic_summer"]["settings"], 0.5, "1 h"),
        ("Crude oil →\ngasoline", F["gasoline_weekly_diff"]["settings"], 168.0, "1–2 wk"),
    ]

    def fmt(dl, dt_h):
        if dl is None:
            return "–"
        h = abs(dl) * dt_h
        return f"{h:.0f} h" if h < 24 else (f"{h / 24:.0f} d" if h < 24 * 14 else f"{h / 168:.0f} wk")

    fig, ax = plt.subplots(figsize=(6.8, 2.8))
    w = 0.36
    for i, (lab, st, dt_h, exp) in enumerate(recs):
        for k, m in enumerate(("published", "atds")):
            r = st[m]; x = i + (k - 0.5) * w
            sig = r["p"] <= 0.05
            ax.bar(x, r["score"], w * 0.92, color=COL[m], alpha=1.0 if sig else 0.35,
                   label=LAB[m] if i == 0 else None)
            ax.plot([x - w / 2, x + w / 2], [r["null95"]] * 2, color=INK, lw=1.3,
                    label="95th pct of rotations (chance)" if i == 0 and k == 0 else None)
            ax.text(x, r["score"] + 2, fmt(r["delay_samples"], dt_h) + ("" if sig else "\nn.s."),
                    ha="center", va="bottom", fontsize=6.6, color=COL[m])
        ax.text(i, -24, f"expected {exp}", ha="center", fontsize=6.8, color=MUTED)
    ax.set_xticks(range(len(recs)), [r[0] for r in recs], fontsize=7.2)
    ax.set_ylabel("TDS score (%)"); ax.set_ylim(0, 118)
    ax.set_title("Score, chance level and delay found (pale bar: not significant, p > 0.05)", loc="left")
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, 1.0), ncol=3)
    fig.savefig(out, dpi=200)
    print("wrote", out)


if __name__ == "__main__":
    o = Path(sys.argv[1]); os.makedirs(o, exist_ok=True)
    fig_bidmc(o / "fig_bidmc.png")
    fig_fields(o / "fig_fields_significance.png")
