"""fig_networks.png: (a) per-link power / false-positive rate vs slow-down r; (b) network test with FDR control. Result files only.

    python experiments/figs_v4/plot_networks.py

Note: the spec names exp6.json W1 for (b), but exp6.json W1 is the d_max experiment; the network experiment is
atds_brief/results/exp5.json V2 (disjoint and common-driver designs, 999 rotations, BH at q = 0.05). V2 has aTDS only.
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _style import COL, INK, MUTED, ROOT, WIDTH, plt, save  # noqa: E402

P = json.load(open(os.path.join(ROOT, "experiments", "atds_v3", "results", "x6_item5_power.json")))
V2 = json.load(open(os.path.join(ROOT, "atds_brief", "results", "exp5.json")))["V2"]
rs = (2, 4, 8)

print("(a) per-link test, 50 coupled + 50 uncoupled pairs, 99 rotations")
print(f"{'r':>2} {'method':10s} {'power':>6s} {'95% CI':>14s} {'FPR':>6s} {'95% CI':>14s}")
for r in rs:
    for m in ("published", "atds"):
        e = P[f"r{r}|{m}"]
        print(f"{r:>2} {m:10s} {e['power']:6.2f} {str(e['power_ci']):>14s} {e['fpr']:6.2f} {str(e['fpr_ci']):>14s}")

fig, (a, b) = plt.subplots(1, 2, figsize=(WIDTH, 2.9), gridspec_kw=dict(width_ratios=[1, 1.1], wspace=0.3))
for m, dx in (("published", -0.04), ("atds", 0.04)):
    x = np.array(rs) * (2 ** dx)
    for key, ls, mk in (("power", "-", "o"), ("fpr", "--", "s")):
        v = np.array([P[f"r{r}|{m}"][key] for r in rs]) * 100
        ci = np.array([P[f"r{r}|{m}"][key + "_ci"] for r in rs]) * 100
        a.errorbar(x, v, yerr=[v - ci[:, 0], ci[:, 1] - v], color=COL[m], ls=ls, marker=mk, ms=4, lw=1.2, capsize=2, elinewidth=0.8)
a.axhline(5, color=INK, ls=":", lw=0.8)
a.text(5.7, 9, "5 %", ha="center", va="bottom", fontsize=7, color=INK)
a.set_xscale("log", base=2); a.set_xticks(rs, [str(r) for r in rs]); a.minorticks_off()
a.set_xlim(1.7, 9.5); a.set_ylim(-3, 103)
a.set_xlabel("slow-down $r$ (time-scale factor)"); a.set_ylabel("per-link rate (%)")
a.set_title("(a) per-link test", loc="left")
h = [plt.Line2D([], [], color=COL["published"], marker="o", ms=4, label="published"),
     plt.Line2D([], [], color=COL["atds"], marker="o", ms=4, label="aTDS"),
     plt.Line2D([], [], color=MUTED, ls="-", label="power (coupled)"),
     plt.Line2D([], [], color=MUTED, ls="--", marker="s", ms=4, label="false positive (uncoupled)")]
a.legend(handles=h, loc="upper center", ncol=2, fontsize=6.6, bbox_to_anchor=(0.5, -0.2), columnspacing=1.0, handlelength=1.8)

# (b) network test
kinds = (("disjoint", "disjoint links"), ("driver", "common driver"))
cats = (("true links found\n(of 30)", "true"), ("false links,\ncoupled nets", "fc"), ("false links,\nuncoupled nets", "fu"))
print("\n(b) network test (aTDS only; 10 coupled + 10 uncoupled networks per design; BH q = 0.05)")
cnt = {}
for k, _ in kinds:
    co = [e for e in V2 if e["kind"] == k and e["coupled"]]
    un = [e for e in V2 if e["kind"] == k and not e["coupled"]]
    assert len(co) == 10 and len(un) == 10
    cnt[k] = dict(true=sum(e["true_found"] for e in co), ntrue=sum(e["n_true"] for e in co),
                  fc=sum(e["false_other"] + e["indirect_flagged"] for e in co), fu=sum(e["false_other"] + e["indirect_flagged"] for e in un),
                  ind=sum(e["indirect_flagged"] for e in co))
    print(f"  {k:9s} true found {cnt[k]['true']}/{cnt[k]['ntrue']}, false links in coupled nets {cnt[k]['fc']} (indirect {cnt[k]['ind']}), in uncoupled nets {cnt[k]['fu']}")
w = 0.36
for j, (k, kl) in enumerate(kinds):
    c = cnt[k]
    for i, (key, col) in enumerate((("true", COL["atds"]), ("fc", MUTED), ("fu", "#c9c8c2"))):
        x = i + (j - 0.5) * w
        b.bar(x, c[key], w * 0.92, color=col, alpha=1.0 if j == 0 else 0.55, edgecolor=col if j else "none", lw=0.6)
        b.text(x, c[key] + 0.5, str(c[key]), ha="center", va="bottom", fontsize=7.5, color=INK)
b.set_xticks(range(3), [c[0] for c in cats], fontsize=6.8)
b.set_ylim(0, 33); b.set_ylabel("links (summed over 10 networks)")
b.grid(axis="x", visible=False)
b.set_title("(b) aTDS network test (BH, $q$ = 0.05)", loc="left")
h = [plt.Rectangle((0, 0), 1, 1, color=INK, alpha=0.85, label="disjoint links"),
     plt.Rectangle((0, 0), 1, 1, color=INK, alpha=0.35, label="common driver")]
b.legend(handles=h, loc="upper right", fontsize=6.8)
save(fig, "fig_networks.png")
