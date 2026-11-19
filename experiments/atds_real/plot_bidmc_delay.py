"""fig_bidmc_delay.png: delay found per ICU patient against the full-record reference delay (result file only).

    python experiments/atds_real/plot_bidmc_delay.py <out_dir>

Reads results/bidmc_ecg_ppg.json. Reference = xcorr_delay / 125 Hz (full-record alignment within half a beat); the delay found is
plotted as reference + wrap(found - reference, rr), i.e. wrapped to [-rr/2, rr/2) around the reference as in bidmc.py (stored as *_err_xcorr_s).
"""
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
FS = 125.0
INK, MUTED = "#262624", "#8a8a85"
COL = {"published": "#2a78d6", "atds": "#eb6834"}
LAB = {"published": "fixed TDS", "atds": "aTDS"}
plt.rcParams.update({
    "font.family": "serif", "font.size": 8, "axes.titlesize": 8.5, "legend.fontsize": 7.2,
    "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": MUTED, "axes.grid": True,
    "grid.color": "#e6e5e0", "grid.linewidth": 0.6, "legend.frameon": False,
    "xtick.color": MUTED, "ytick.color": MUTED, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})


def main(outdir):
    B = json.load(open(os.path.join(HERE, "results", "bidmc_ecg_ppg.json")))
    recs = [r for r in B["records"] if r["coupled"]]
    ref = np.array([r["xcorr_delay"] / FS for r in recs])
    fig, (a, b) = plt.subplots(1, 2, figsize=(6.8, 3.0), gridspec_kw={"width_ratios": [1.1, 1], "wspace": 0.3})
    lo, hi = float(ref.min()) - 0.05, float(ref.max()) + 0.05
    a.plot([lo, hi], [lo, hi], color=MUTED, lw=0.8, ls="--")
    errs, txt = {}, []
    for m in ("published", "atds"):
        e = np.array([np.nan if r[f"{m}_err_xcorr_s"] is None else r[f"{m}_err_xcorr_s"] for r in recs])
        ok = np.isfinite(e)
        n_in = int(np.sum(np.abs(e[ok]) <= 0.04))
        print(f"D {LAB[m]}: {ok.sum()} patients with a delay, {n_in} within 0.04 s of the reference; median err {np.median(e[ok]):.3f} s; "
              f"summary says {B['summary'][m]['n_within_0p04s_of_xcorr']}/{B['summary'][m]['n_with_delay']}")
        assert n_in == B["summary"][m]["n_within_0p04s_of_xcorr"] and int(ok.sum()) == B["summary"][m]["n_with_delay"]
        a.scatter(ref[ok], ref[ok] + e[ok], s=10, color=COL[m], lw=0, alpha=0.85, label=LAB[m])
        errs[m] = e[ok]
        txt.append(f"{LAB[m]}: {n_in}/{int(ok.sum())} within 0.04 s")
    a.set_xlabel("reference delay: full-record alignment (s)"); a.set_ylabel("delay found (s)")
    a.set_title("(a) delay found against reference", loc="left")
    a.text(0.03, 0.97, "\n".join(txt), transform=a.transAxes, va="top", fontsize=7.2, color=INK)
    a.legend(loc="lower right", markerscale=1.3)
    bins = np.arange(-0.34, 0.341, 0.02)
    for m in ("published", "atds"):
        b.hist(errs[m], bins=bins, color=COL[m], alpha=0.6, label=LAB[m])
    b.axvspan(-0.04, 0.04, color=MUTED, alpha=0.12, lw=0)
    b.set_xlabel("delay found minus reference (s)"); b.set_ylabel("patients")
    b.set_title("(b) difference to the reference", loc="left")
    b.legend(loc="upper right")
    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, "fig_bidmc_delay.png")
    fig.savefig(out, dpi=200)
    print("wrote", out)


if __name__ == "__main__":
    main(sys.argv[1])
