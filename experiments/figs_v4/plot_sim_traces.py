"""fig_sim_traces.png: delay found in each window (tau0) for three simulated cases, published (left) vs aTDS (right).

    ATDS_CODE=<checkout of claude/practical-mccarthy-n2mec6> python experiments/figs_v4/plot_sim_traces.py

Same generators, seeds and aTDS parameter-setting code as the original experiments (imported, not re-implemented):
  (i)   slow-down r = 4, atds_brief/experiments/exp4.py U1 (test pair seed 60000 + 4*1000 + 0, aTDS set on 5 training pairs)
  (ii)  regime switching, experiments/atds_v3/e3_synthetic.py (first coupled pair, seed 2000000; whole-record aTDS)
  (iii) drifting delay, e3_synthetic.py (first pair, seed 8000000)
Dashed line = true delay in pyTDS sign convention (tau0 = -d when x leads).
"""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "atds_v3"))
from _style import COL, GREY, INK, ROOT, WIDTH, plt, save  # noqa: E402

import common as C  # noqa: E402  (needs ATDS_CODE)
import e3_synthetic as E3  # noqa: E402

t0 = time.time()
NS = E3.NS


def trace(x, y, p):
    tau, t, _ = C.time_delay_interaction(x, y, p)
    lab = np.asarray(C.stable_label(tau, p)).astype(bool)
    return np.asarray(tau), np.asarray(t), lab, float(C.tds_score(lab))


def pd(p):
    return dict(window=p.window, step=p.overlap, max_lag=p.max_lag, tol=p.tolerance)


# ---- (i) slow-down r = 4 (exp4.U1 design)
r = 4; tc, d, c = 2 * r, 5 * r, 0.3
train = [C.pair(800000 + r * 10 + j, tc, d, c) for j in range(5)]
ad_i, fail = C.E4.atds2(train)
x_i, y_i = C.pair(60000 + r * 1000 + 0, tc, d, c)
assert len(x_i) == 20000 and x_i.ndim == 1 and y_i.ndim == 1
# ---- (ii) regime switching
P_reg = E3._setup_regime(0)
x_r, y_r, A_r = E3.regime_pair(2000000 + 0, 0.6)
# ---- (iii) drifting delay
P_dr = E3._setup_drift(0)
x_d, y_d = E3.drift_pair(8000000 + 0)
for v in (x_r, y_r, x_d, y_d):
    assert v.shape == (NS,)

cases = [
    dict(title="Slow-down $r$ = 4 ($\\tau_c$ = 8, delay 20, coupling 0.3)", x=x_i, y=y_i, P=C.FIXED, A=ad_i, ylim=140, kind="const", d=d),
    dict(title="Regime switching (A: coupled, delay 20; B: uncoupled; blocks of 2000)", x=x_r, y=y_r, P=C.FIXED, A=P_reg["atds"], ylim=290, kind="regime", blocks=A_r),
    dict(title="Drifting delay (10 to 30 samples, coupling 0.6)", x=x_d, y=y_d, P=C.FIXED, A=P_dr["atds"], ylim=140, kind="drift"),
]

# ---- sanity against result files
U1 = [e for e in json.load(open(os.path.join(ROOT, "atds_brief", "results", "exp4.json")))["U1"] if e["r"] == r][0]
e3 = json.load(open(os.path.join(ROOT, "experiments", "atds_v3", "results", "e3_synthetic.json")))
chk = [("(i) published", pd(C.FIXED), dict(window=60, step=30, max_lag=30, tol=1)),
       ("(i) aTDS", pd(ad_i), dict(window=U1["params"]["atds"]["window"], step=U1["params"]["atds"]["window"] // 2,
                                  max_lag=U1["params"]["atds"]["max_lag"], tol=U1["params"]["atds"]["tol"])),
       ("(ii) aTDS", pd(P_reg["atds"]), e3["regime_switching"]["atds"]["params"]),
       ("(iii) aTDS", pd(P_dr["atds"]), e3["drifting_delay"]["atds"]["params"])]
print("parameter sanity (computed vs result JSON):")
for nm, got, ref in chk:
    print(f"  {nm:14s} {got}  {'ok' if got == ref else 'MISMATCH vs ' + str(ref)}")

fig, axes = plt.subplots(3, 2, figsize=(WIDTH, 1.95 * 3), sharey="row", gridspec_kw=dict(hspace=1.05, wspace=0.08))
for k, cs in enumerate(cases):
    print(cs["title"].replace("$", ""))
    for j, (name, p) in enumerate((("published", cs["P"]), ("atds", cs["A"]))):
        ax = axes[k, j]
        tau, t, lab, sc = trace(cs["x"], cs["y"], p)
        tt = (np.arange(len(tau)) * p.overlap + p.window / 2) / 1000.0  # window centres
        if cs["kind"] == "regime":
            Ab = cs["blocks"]; edges = np.flatnonzero(np.diff(np.concatenate([[0], Ab.astype(int), [0]])))
            for s_, e_ in zip(edges[::2], edges[1::2]):
                ax.axvspan(s_ / 1000, e_ / 1000, color=INK, alpha=0.07, lw=0)
                ax.hlines(-20, s_ / 1000, e_ / 1000, color=INK, ls="--", lw=0.8)
        ax.scatter(tt[~lab], tau[~lab], s=2.5, color=GREY, lw=0, label="unstable window")
        ax.scatter(tt[lab], tau[lab], s=3.5, color=COL[name], lw=0, label="stable window")
        ax.axhspan(-p.max_lag, p.max_lag, color=COL[name], alpha=0.06, lw=0)
        if cs["kind"] == "const":
            ax.axhline(-cs["d"], color=INK, ls="--", lw=0.8)
        elif cs["kind"] == "drift":
            xs = np.linspace(0, NS, 200); ax.plot(xs / 1000, -(10 + 20 * xs / NS), color=INK, ls="--", lw=0.8)
        ax.set_ylim(-cs["ylim"], cs["ylim"]); ax.set_xlim(0, NS / 1000)
        ax.set_xlabel("time (10$^3$ samples)")
        ax.set_title(f"{'published' if name == 'published' else 'aTDS'}: TDS {sc:.0f} %", loc="left", color=COL[name])
        if j == 0:
            ax.set_ylabel("delay (samples)")
        print(f"   {name:9s} score {sc:.1f}  params {pd(p)}  stable windows {int(lab.sum())}/{len(lab)}")
    axes[k, 0].annotate(cs["title"], xy=(0, 1.32), xycoords="axes fraction", fontsize=8.5, color=INK, fontweight="bold")
save(fig, "fig_sim_traces.png")
print(f"wall {time.time() - t0:.0f}s")
