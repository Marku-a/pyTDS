"""fig_hard_cases.png: regime switching, drifting delay, nonlinear coupling (published vs aTDS). Result files only.

    python experiments/figs_v4/plot_hard.py
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _style import COL, INK, MUTED, ROOT, WIDTH, plt, save  # noqa: E402

RS = os.path.join(ROOT, "experiments", "atds_v3", "results")
E3 = json.load(open(os.path.join(RS, "e3_synthetic.json")))
MG = json.load(open(os.path.join(RS, "x6_item2_margin.json")))["results"]

# ---- (a) regime switching
reg = E3["regime_switching"]
groups = [("published", reg["published"], "published"), ("aTDS", reg["atds"], "atds"), ("aTDS + margin", MG["atds|margin=True"], "atds")]
bars = [("regime A, coupled pairs", "regimeA_coupled_mean"), ("regime B, coupled pairs", "false_stability_B_coupled_mean"),
        ("regime B, controls", "false_stability_B_control_mean")]
print("(a) regime switching: % stable windows")
print(f"{'':16s}" + "".join(f"{b[0]:>26s}" for b in bars))
for g, d, _ in groups:
    print(f"{g:16s}" + "".join(f"{d[k]:26.2f}" for _, k in bars))
paper = {"published": (5.6, 0.9, None), "aTDS": (98.9, 8.1, None), "aTDS + margin": (None, 3.8, 5.5)}
bad = []
for g, d, _ in groups:
    for (nm, k), pv in zip(bars, paper[g]):
        if pv is not None and abs(round(d[k], 1) - pv) > 0.051:
            bad.append((g, nm, round(d[k], 1), pv))

# ---- (b), (c)
dr = E3["drifting_delay"]
nm_, nl_ = E3["nonlinear_monotone"], E3["nonlinear_non_monotone"]
print("(b) drifting delay: score %, median |error| (samples):",
      {m: (dr[m]["score_mean"], dr[m]["median_abs_err_pooled"]) for m in ("published", "atds")})
print("(c) AUC monotone / non-monotone:", {m: (nm_[m]["auc"], nl_[m]["auc"]) for m in ("published", "atds")})
for lab, got, exp in (("drift score pub", dr["published"]["score_mean"], 8.5), ("drift score aTDS", dr["atds"]["score_mean"], 86.5),
                      ("drift err pub", dr["published"]["median_abs_err_pooled"], 0.30), ("drift err aTDS", dr["atds"]["median_abs_err_pooled"], 0.25),
                      ("AUC mono pub", nm_["published"]["auc"], 0.997), ("AUC mono aTDS", nm_["atds"]["auc"], 1.0),
                      ("AUC nonmono pub", nl_["published"]["auc"], 0.803), ("AUC nonmono aTDS", nl_["atds"]["auc"], 0.747)):
    tol = 0.051 if "score" in lab else 0.0051
    if abs(got - exp) > tol:
        bad.append((lab, got, exp))
print("mismatches vs paper text:", bad if bad else "none")

fig = plt.figure(figsize=(WIDTH, 3.1))
gs = fig.add_gridspec(1, 3, width_ratios=[2.2, 1.0, 1.0], wspace=0.5)
a = fig.add_subplot(gs[0]); b = fig.add_subplot(gs[1]); c = fig.add_subplot(gs[2])

w = 0.26
shade = [1.0, 0.62, 0.3]
for i, (g, d, m) in enumerate(groups):
    for j, (_, k) in enumerate(bars):
        x = i + (j - 1) * w
        v = d[k]
        a.bar(x, v, w * 0.92, color=COL[m], alpha=shade[j], edgecolor=COL[m], lw=0.6)
        a.text(x, v + 1.5, f"{v:.1f}", ha="center", va="bottom", fontsize=6.4, color=INK, rotation=90)
a.axhline(5, color=INK, ls=":", lw=0.8)
a.set_xticks(range(3), [g[0] for g in groups], fontsize=7.4)
a.set_ylim(0, 150); a.set_yticks([0, 25, 50, 75, 100]); a.set_xlim(-0.55, 2.6)
a.set_ylabel("stable windows (%)"); a.grid(axis="x", visible=False)
a.set_title("(a) regime switching", loc="left")
h = [plt.Rectangle((0, 0), 1, 1, color=INK, alpha=s, label=l) for s, l in zip(shade, ("regime A, coupled pairs", "regime B, coupled pairs", "regime B, controls"))]
h.append(plt.Line2D([], [], color=INK, ls=":", lw=0.8, label="5 % budget"))
a.legend(handles=h, loc="upper left", ncol=2, fontsize=6.4, bbox_to_anchor=(-0.02, 1.02), columnspacing=1.0, handlelength=1.2, handletextpad=0.4)

# (b) drifting delay: bars = % stable windows; annotation = median abs error
for i, m in enumerate(("published", "atds")):
    v = dr[m]["score_mean"]
    b.bar(i, v, 0.62, color=COL[m])
    b.text(i, v + 2, f"{v:.1f} %", ha="center", va="bottom", fontsize=7.2, color=INK)
b.set_xticks([0, 1], [f"{n}\nerr. {dr[m]['median_abs_err_pooled']:.2f}" for n, m in (("published", "published"), ("aTDS", "atds"))], fontsize=6.8); b.set_ylim(0, 105)
b.set_ylabel("stable windows (%)"); b.grid(axis="x", visible=False)
b.set_title("(b) drifting delay", loc="left")
b.set_xlabel("(median |delay error|, samples)", fontsize=6.8, color=MUTED)

# (c) AUC
for i, (kind, d) in enumerate((("monotone", nm_), ("non-monotone", nl_))):
    for j, m in enumerate(("published", "atds")):
        x = i + (j - 0.5) * 0.36
        c.bar(x, d[m]["auc"], 0.33, color=COL[m])
        c.text(x, d[m]["auc"] + 0.012, f"{d[m]['auc']:.3f}".rstrip("0").rstrip(".") if d[m]["auc"] < 1 else "1.0", ha="center", va="bottom", fontsize=6.8, color=INK, rotation=90)
c.axhline(0.5, color=INK, ls=":", lw=0.8)
c.set_xticks([0, 1], ["monotone", "non-\nmonotone"], fontsize=7.4); c.set_ylim(0.4, 1.12); c.set_xlim(-0.6, 1.6)
c.set_ylabel("AUC"); c.grid(axis="x", visible=False)
c.set_title("(c) nonlinear", loc="left")
h = [plt.Rectangle((0, 0), 1, 1, color=COL["published"], label="published"), plt.Rectangle((0, 0), 1, 1, color=COL["atds"], label="aTDS")]
fig.legend(handles=h, loc="lower center", ncol=2, bbox_to_anchor=(0.72, -0.12), fontsize=7)
save(fig, "fig_hard_cases.png")
