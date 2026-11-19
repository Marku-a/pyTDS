"""fig_signal_classes.png: AUC and 'delay found' per signal class of tab:signals (published vs aTDS), from result files only.

    python experiments/figs_v4/plot_classes.py
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _style import COL, INK, MUTED, ROOT, WIDTH, plt, save  # noqa: E402

R = os.path.join(ROOT, "atds_brief", "results")
E4 = json.load(open(os.path.join(R, "exp4.json")))
E5 = json.load(open(os.path.join(R, "exp5.json")))
E6 = json.load(open(os.path.join(R, "exp6.json")))


def mk(d):
    """(auc, auc_ci or None, found %, found_ci or None) from a summarize()-style dict."""
    ac, hc = d.get("auc_ci"), d.get("hit_ci")
    return dict(auc=d["auc"], auc_ci=ac, found=100 * d["hit"], found_ci=None if hc is None else [100 * max(0, hc[0]), 100 * hc[1]])


# each row: label, group, published, aTDS, variant (or None), variant label, expected table values
# expected = (pub AUC, atds AUC, pub found %, atds found %); None where the table does not state it
rows = []
for h, e in zip((0.6, 0.75, 0.9), E4["U5"]):
    rows.append(dict(label=f"fGn $H$ = {h}", pub=mk(e["fixed"]), atds=mk(e["atds"]), var=None))
exp_tab = {  # tab:signals numbers: (pub AUC, atds AUC, pub found, atds found)
    "fGn $H$ = 0.6": (0.60, 0.59, 28, 23), "fGn $H$ = 0.75": (0.70, 0.65, 33, 23), "fGn $H$ = 0.9": (0.46, 0.82, 15, 85),
    "Rhythm, $d$ = 8": (0.97, 1.00, 43, 100), "Rhythm, $d$ = 30": (0.46, 0.995, 0, 0),
    "Rulkov, $d$ = 16": (0.73, 1.00, None, None), "Rulkov, $d$ = 32": (0.79, 1.00, None, None),
    "Fast ($\\tau_c$=1), $d$ = 50": (0.56, 0.56, 0, 0),
    "Memory change, 1st half": (0.98, 1.00, None, None), "Memory change, 2nd half": (0.71, 0.96, 45, 90),
    "Mixed: slow pair": (0.44, 0.99, None, None), "Mixed: fast pair": (1.00, 1.00, None, None),
    "Slow $x$ drives fast $y$": (0.83, 1.00, 57, 100),
}
for d, e in zip((8, 30), E4["U6"]):
    rows.append(dict(label=f"Rhythm, $d$ = {d}", pub=mk(dict(auc=e["fixed"]["auc"], hit=e["fixed"]["hit"])),
                     atds=mk(dict(auc=e["atds"]["auc"], hit=e["atds"]["hit"])), var=None))
for e in E5["V1"]:
    rows.append(dict(label=f"Rulkov, $d$ = {e['d']}", pub=mk(e["fixed"]), atds=mk(e["atds"]), var=None,
                     note=f"median delay {e['fixed']['median_reported']:g} / {e['atds']['median_reported']:g}"))
w = [e for e in E6["W1"] if e["d"] == 50][0]
rows.append(dict(label="Fast ($\\tau_c$=1), $d$ = 50", pub=mk(w["fixed"]), atds=mk(w["atds_default"]),
                 var=mk(w["atds_dmax60"]), vlabel="aTDS variant ($d_{\\max}$ = 60; per half; joint $B_{xy}$)"))
u10 = E4["U10"][0]
for half, nm in (("first", "1st"), ("second", "2nd")):
    rows.append(dict(label=f"Memory change, {nm} half", pub=mk(u10[half]["fixed"]), atds=mk(u10[half]["atds_whole"]),
                     var=mk(u10[half]["atds_per_segment"]), vlabel="aTDS variant ($d_{\\max}$ = 60; per half; joint $B_{xy}$)"))
u9 = E4["U9"][0]
for w_, nm in (("slow", "slow pair"), ("fast", "fast pair")):
    rows.append(dict(label=f"Mixed: {nm}", pub=mk(u9["fixed"][w_]), atds=mk(u9["system_wide"][w_]), var=None))
fs = u9["fast_slow_pair"]
rows.append(dict(label="Slow $x$ drives fast $y$", pub=mk(fs["fixed"]), atds=mk(fs["per_pair_maxB"]), var=mk(fs["pairwise_Bxy"]),
                 vlabel="aTDS variant ($d_{\\max}$ = 60; per half; joint $B_{xy}$)"))
# order = table order
order = list(exp_tab)
rows = sorted(rows, key=lambda r: order.index(r["label"]))
assert [r["label"] for r in rows] == order, [r["label"] for r in rows]

# ---- check table
print(f"{'row':28s} {'pub AUC':>8s} {'aTDS AUC':>9s} {'pub found':>10s} {'aTDS found':>11s}  table  -> mismatch")
bad = []
for r in rows:
    pa, aa, pf, af = exp_tab[r["label"]]
    got = (r["pub"]["auc"], r["atds"]["auc"], r["pub"]["found"], r["atds"]["found"])
    flags = []
    for nm, g, e, tol in zip(("pub AUC", "aTDS AUC", "pub found", "aTDS found"), got, (pa, aa, pf, af), (0.0051, 0.0051, 0.51, 0.51)):
        if e is not None and abs(g - e) > tol:
            flags.append(f"{nm}: json {g} vs table {e}")
    bad += [(r["label"], f) for f in flags]
    extra = (" [" + r["note"] + "]") if r.get("note") else ""
    var = f" variant AUC {r['var']['auc']} found {r['var']['found']:.0f}" if r["var"] else ""
    print(f"{r['label'].replace('$', ''):28s} {got[0]:8.3f} {got[1]:9.3f} {got[2]:9.0f}% {got[3]:10.0f}%  {'ok' if not flags else 'MISMATCH'}{extra}{var}")
print("mismatches:", bad if bad else "none")

# ---- figure
n = len(rows)
ypos = np.arange(n)[::-1].astype(float)
fig, axes = plt.subplots(1, 2, figsize=(WIDTH, 0.33 * n + 1.3), sharey=True, gridspec_kw=dict(wspace=0.06))
for ax, key, ci, xl, vl in ((axes[0], "auc", "auc_ci", "AUC (coupled vs uncoupled)", 0.5), (axes[1], "found", "found_ci", "delay found (% of coupled pairs)", None)):
    for y, r in zip(ypos, rows):
        hv = r["var"] is not None
        for name, dy in (("pub", 0.24 if hv else 0.16), ("atds", 0.0 if hv else -0.16)):
            d = r[name]; c = COL["published" if name == "pub" else "atds"]
            ax.plot(d[key], y + dy, "o", color=c, ms=4.2, zorder=3)
            if d[ci] is not None:
                ax.plot(d[ci], [y + dy] * 2, color=c, lw=1.1, zorder=2)
        if r["var"] is not None:
            v = r["var"]
            ax.plot(v[key], y - 0.24, "o", mfc="white", mec=COL["atds"], mew=1.2, ms=5.2, zorder=4)
            if v[ci] is not None and v[ci][0] != v[ci][1]:
                ax.plot(v[ci], [y - 0.24] * 2, color=COL["atds"], lw=1.1, alpha=0.5, zorder=2)
    if vl is not None:
        ax.axvline(vl, color=INK, ls=":", lw=0.8)
    ax.set_xlabel(xl)
    ax.set_ylim(-0.7, n - 0.3)
    ax.grid(axis="y", visible=False)
axes[0].set_xlim(0.3, 1.04); axes[1].set_xlim(-6, 106)
axes[0].set_yticks(ypos, [r["label"] for r in rows], fontsize=7.4)
axes[0].tick_params(axis="y", length=0)
axes[1].tick_params(axis="y", length=0)
axes[0].set_title("(a) AUC", loc="left"); axes[1].set_title("(b) delay found", loc="left")
h = [plt.Line2D([], [], marker="o", color=COL["published"], ls="", ms=4.5, label="published (60 / $\\pm$30 / $\\pm$1)"),
     plt.Line2D([], [], marker="o", color=COL["atds"], ls="", ms=4.5, label="aTDS"),
     plt.Line2D([], [], marker="o", mfc="white", mec=COL["atds"], mew=1.2, color=COL["atds"], ls="", ms=5.2,
                label="aTDS variant: $d_{\\max}$ = 60 (fast); per-half window (memory); joint $B_{xy}$ (slow$\\to$fast)")]
fig.legend(handles=h, loc="lower center", ncol=1, bbox_to_anchor=(0.5, -0.1), fontsize=7)
save(fig, "fig_signal_classes.png")
