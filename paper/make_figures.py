"""Build every figure of the aTDS paper from the saved experiment results.

Reads only atds_brief/results/*.json (no experiment is re-run) and writes
PDF figures to paper/figures/.  Run from the repository root:

    python paper/make_figures.py
"""
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "atds_brief" / "results"
OUT = ROOT / "paper" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

# Validated categorical palette (fixed order, one hue per entity).
C_FIXED, C_ATDS, C_RESC, C_EXTRA = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
C_ORACLE, INK, MUTED = "#6b6b66", "#262624", "#8a8a85"

plt.rcParams.update({
    "font.family": "serif", "font.size": 8.5, "axes.titlesize": 9,
    "axes.labelsize": 8.5, "legend.fontsize": 7.5, "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": "#e6e5e0", "grid.linewidth": 0.6,
    "lines.linewidth": 1.6, "lines.markersize": 4.5, "legend.frameon": False,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})


def load(name):
    return json.loads((RES / f"{name}.json").read_text())


theory, exp4 = load("theory"), load("exp4")


def save(fig, name):
    fig.savefig(OUT / f"{name}.pdf")
    po = ROOT / "paper" / "paperorchestra" / "raw_materials" / "figures"
    po.mkdir(parents=True, exist_ok=True)
    fig.savefig(po / f"{name}.png", dpi=200)
    plt.close(fig)
    print("wrote", OUT / f"{name}.pdf")


# 1. Motivating example: same coupled pair, fixed vs aTDS parameters.
def fig_example():
    ex = exp4["U12"][0]["example"]
    fig, axes = plt.subplots(2, 1, figsize=(6.6, 3.4), sharex=True)
    for ax, key, col, title in [
        (axes[0], "fixed", C_FIXED, "(a) fixed TDS"),
        (axes[1], "atds", C_ATDS, "(b) aTDS"),
    ]:
        d = ex[key]
        t, tau, st = map(np.asarray, (d["t"], d["tau"], d["stable"]))
        m = t <= 12000
        t, tau, st = t[m], tau[m], st[m].astype(bool)
        ax.scatter(t[~st], tau[~st], s=10, facecolors="none", edgecolors=col,
                   linewidths=0.7, label="unstable window")
        ax.scatter(t[st], tau[st], s=10, color=col, label="stable window")
        ax.axhline(-20, color=INK, ls="--", lw=0.9)
        p = d["params"]
        ax.set_title(f"{title}: L = {p['window']}, $\\tau_{{\\max}}$ = ±{p['max_lag']}, "
                     f"tol = ±{p['tol']};  TDS = {d['score']:.1f} %", loc="left")
        ax.set_ylabel(r"$\tau_0$ (samples)")
        ax.set_ylim(-p["max_lag"] - 3, p["max_lag"] + 3)
    axes[0].legend(loc="lower right", bbox_to_anchor=(1.0, 1.12), ncol=2, borderaxespad=0)
    axes[1].text(12000, -45, "true delay (−20)", ha="right", va="top", fontsize=7.5)
    axes[1].set_xlabel("time (samples)")
    fig.tight_layout()
    save(fig, "fig_example_fixed_vs_atds")


# 2. Window rule: independent samples per window, and AUC vs window length.
def fig_window():
    fig, (a, b) = plt.subplots(1, 2, figsize=(6.6, 2.5))
    mem = theory["memory"]
    tc = [m["tc"] for m in mem]
    a.plot(tc, [m["indep_fixed"] for m in mem], "o-", color=C_FIXED, label="fixed L = 60")
    a.plot(tc, [m["indep_atds"] for m in mem], "s-", color=C_ATDS, label="aTDS L = 30·B")
    a.set_xscale("log", base=2)
    a.set_yscale("log")
    a.set_xlabel(r"correlation time $\tau_c$ (samples)")
    a.set_ylabel("independent samples per window, L/B")
    a.set_title("(a) information per window", loc="left")
    a.legend()

    marks = {2: ("faster: $\\tau_c$=4, d=10", "o", C_FIXED), 8: ("slower: $\\tau_c$=16, d=40", "s", C_RESC)}
    for r, (lab, mk, col) in marks.items():
        rows = [x for x in exp4["U2"] if x["r"] == r]
        L = np.array([x["L"] for x in rows])
        auc = np.array([x["auc"] for x in rows])
        ci = np.array([x["auc_ci"] for x in rows])
        b.errorbar(L, auc, yerr=[auc - ci[:, 0], ci[:, 1] - auc], fmt=mk + "-",
                   color=col, capsize=2, lw=1.4, label=lab)
    for x in exp4["U2m"]:
        b.axvline(x["L"], color=marks[x["r"]][2], ls=":", lw=1)
        b.plot(x["L"], x["auc"], "D", color=C_ATDS, ms=6, zorder=5,
               label="aTDS choice" if x["r"] == 2 else None)
    b.axvline(60, color=MUTED, ls="--", lw=1)
    b.text(58, 1.0, "L = 60", color=MUTED, fontsize=7, ha="right")
    b.set_xscale("log", base=2)
    b.set_xlabel("window L (samples)")
    b.set_ylabel("AUC (coupled vs uncoupled)")
    b.set_ylim(0.4, 1.03)
    b.set_title("(b) detection vs window", loc="left")
    b.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0)
    fig.tight_layout()
    save(fig, "fig_window_rule")


# 3. On/off coupling: balanced accuracy vs independent samples per window.
def fig_onoff():
    fig, ax = plt.subplots(figsize=(3.3, 2.5))
    for block, mk, col in [(1000, "o", C_FIXED), (2000, "s", C_RESC), (4000, "^", C_EXTRA)]:
        rows = [x for x in exp4["U3"] if x["block"] == block]
        ax.plot([x["neff"] for x in rows], [x["bal_acc"] for x in rows], mk + "-",
                color=col, label=f"blocks of {block}")
    ax.axvline(30, color=C_ATDS, ls="--", lw=1)
    ax.text(31, 0.5, "aTDS: 30", color=INK, fontsize=7)
    ax.set_xlabel("independent samples per window, L/B")
    ax.set_ylabel("balanced accuracy")
    ax.legend(loc="lower right")
    fig.tight_layout()
    save(fig, "fig_onoff_tradeoff")


# 4. Lag range: reported delay vs true delay, +-30 vs aTDS.
def fig_lag():
    rows = exp4["U4"]
    d = np.array([x["d"] for x in rows], float)
    fig, (a, b) = plt.subplots(1, 2, figsize=(6.6, 2.6))
    rng = np.random.default_rng(0)
    jit = rng.uniform(-2, 2, len(d))
    for key, mk, col, lab in [("lag30", "o", C_FIXED, r"lag $\pm 30$"),
                              ("atds", "s", C_ATDS, r"aTDS lag $(L-1)/2$")]:
        f = np.array([x[key]["found"] for x in rows], float)
        s = np.array([x[key]["score"] for x in rows], float)
        ok = np.abs(f - d) <= np.maximum(1, 0.1 * d)
        for axis, y in ((a, f), (b, s)):
            axis.scatter((d + jit)[ok], y[ok], marker=mk, s=14, color=col, label=lab)
            axis.scatter((d + jit)[~ok], y[~ok], marker=mk, s=14, facecolors="none",
                         edgecolors=col, linewidths=0.8)
    a.plot([0, 165], [0, 165], color=INK, ls="--", lw=0.8)
    a.set_xlabel("true delay d (samples)")
    a.set_ylabel("reported delay (samples)")
    a.set_title("(a) reported delay", loc="left")
    a.legend(loc="upper left")
    b.set_xlabel("true delay d (samples)")
    b.set_ylabel("TDS score (%)")
    b.set_title("(b) TDS score of the same runs", loc="left")
    fig.tight_layout()
    save(fig, "fig_lag_range")


# 5. Tolerance calibration on the example pair.
def fig_tol():
    tc = exp4["U12"][0]["tol_curve"]
    tol = [x["tol"] for x in tc]
    fig, (a, b) = plt.subplots(1, 2, figsize=(6.6, 2.4))
    a.plot(tol, [x["null"] for x in tc], "o-", color=C_FIXED, ms=3, label="chance score C(tol)")
    a.plot(tol, [x["ub"] for x in tc], "--", color=C_FIXED, label="upper bound C + 2 SE")
    a.axhline(5, color=INK, lw=0.8, ls=":")
    a.axvline(13, color=C_ATDS, lw=1)
    a.text(13.5, 1, "chosen ±13", fontsize=7)
    a.set_xlabel("tolerance (samples)")
    a.set_ylabel("score of rotated copies (%)")
    a.set_title("(a) coupling-free copies", loc="left")
    a.legend(loc="upper left")
    b.plot(tol, [x["real"] for x in tc], "s-", color=C_ATDS, ms=3)
    b.axvline(13, color=C_ATDS, lw=1)
    b.set_ylim(0, 100)
    b.set_xlabel("tolerance (samples)")
    b.set_ylabel("TDS score of real pair (%)")
    b.set_title("(b) the coupled pair", loc="left")
    fig.tight_layout()
    save(fig, "fig_tolerance_calibration")


# 6. Slow-down test: four settings.
def fig_slowdown():
    rows = exp4["U1"]
    r = np.array([x["r"] for x in rows])
    fig, (a, b) = plt.subplots(1, 2, figsize=(6.6, 2.6))
    for key, mk, col, lab, ls in [
        ("fixed", "o", C_FIXED, "fixed TDS 60/±30/±1", "-"),
        ("rescaled", "^", C_RESC, "rescaled by hand (knows r)", "-"),
        ("atds", "s", C_ATDS, "aTDS", "-"),
        ("oracle", "D", C_ORACLE, "oracle (labelled grid)", ":"),
    ]:
        auc = np.array([x[key]["auc"] for x in rows])
        ci = np.array([x[key]["auc_ci"] for x in rows])
        hit = np.array([x[key]["hit"] for x in rows])
        hci = np.array([x[key]["hit_ci"] for x in rows])
        a.errorbar(r, auc, yerr=[auc - ci[:, 0], ci[:, 1] - auc], fmt=mk, ls=ls,
                   color=col, capsize=2, lw=1.4, label=lab)
        b.errorbar(r, 100 * hit, yerr=[100 * (hit - hci[:, 0]), 100 * (hci[:, 1] - hit)],
                   fmt=mk, ls=ls, color=col, capsize=2, lw=1.4)
    for axis in (a, b):
        axis.set_xscale("log", base=2)
        axis.set_xticks(r, [str(v) for v in r])
        axis.set_xlabel("slow-down factor r")
    a.set_ylabel("AUC")
    a.set_ylim(0.4, 1.03)
    a.set_title("(a) detection", loc="left")
    h, lab = a.get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.06))
    b.set_ylabel("true delay found (%)")
    b.set_title("(b) delay recovery", loc="left")
    fig.tight_layout()
    save(fig, "fig_slowdown_baselines")


if __name__ == "__main__":
    fig_example()
    fig_window()
    fig_onoff()
    fig_lag()
    fig_tol()
    fig_slowdown()
