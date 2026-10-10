"""Narrated explainer video: how adaptive TDS (aTDS) works.

Voice: Piper TTS (pip `piper-tts`, voice en-us-lessac-medium from the rhasspy/piper v0.0.2 GitHub release).
Real numbers (scenes 4-7) come from the R5 worked example (examples.json, subject #29, Resp-theta);
scenes 2-3 use synthetic signals made here.

usage: python make_video.py --examples examples.json --voice en-us-lessac-medium.onnx --work DIR --out atds_explainer.mp4
"""
import argparse
import json
import os
import subprocess
import wave

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.animation import FFMpegWriter  # noqa: E402

W, H, FPS = 12.8, 7.2, 24
BG, INK, INK2, GRID = "#141517", "#eceae4", "#a9a8a0", "#2a2b2f"
TDS, ATDS, NULL, AQUA = "#3987e5", "#d95926", "#8a8c92", "#199e70"
plt.rcParams.update({"font.family": "DejaVu Sans", "text.color": INK, "axes.labelcolor": INK2, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.edgecolor": "#55565b", "figure.facecolor": BG, "axes.facecolor": BG,
                     "savefig.facecolor": BG, "font.size": 15})

SCENES = [
    ("title", "This is adaptive Time Delay Stability, or a T D S. In about two minutes: what it changes in classic TDS, and why."),
    ("tds", "Classic TDS slides a window along two signals. In each window it finds the delay at which the two signals line up best. "
            "If that delay stays nearly the same for four out of five windows in a row, the two signals are called coupled. "
            "The TDS score is the share of time they are coupled."),
    ("problem", "The settings are fixed: a sixty second window, and a delay that may wander by only one second. "
                "But slow signals remember their past for minutes. A sixty second window then holds less than one independent sample, "
                "so the delay it finds is mostly noise, and stability becomes a matter of luck."),
    ("memory", "So a T D S first measures each signal's memory. It takes the autocorrelation, and adds up its squares until the first zero crossing. "
               "This Bartlett factor says how many samples are worth one independent sample."),
    ("window", "Step two. The window is made long enough to hold about thirty independent samples of the slowest signal. "
               "Because that can be very long, you set a cap, which states the time resolution you need. Here, five minutes."),
    ("lag", "Step three. The window moves by half its length, and the delay search grows with it, up to half the window, "
            "but never beyond half of any rhythm in the data, where delays become ambiguous."),
    ("tolerance", "Step four, the tolerance. One signal is rotated in time by a large random amount. This keeps each signal's own shape but destroys any coupling. "
                  "a T D S then picks the widest tolerance at which these fake pairs still look coupled at most five percent of the time. "
                  "Here that is fifteen seconds."),
    ("summary", "The result is one set of parameters for the whole recording. It uses only each signal's own properties and a coupling free test, "
                "so it cannot be tuned toward the answer. Then TDS runs exactly as before."),
    ("when", "When does it help? When the fixed window is much shorter than the signals' memory, and the recording is long. "
             "What does it cost? Coarser timing, less precise delays, and short states or short recordings lose windows. "
             "That is adaptive TDS."),
]


def say(text, voice, path):
    subprocess.run(["python3", "-m", "piper", "-m", voice, "-f", path, "--sentence-silence", "0.25"], input=text.encode(), check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    with wave.open(path) as w:
        return w.getnframes() / w.getframerate()


def ar_pair(n, tau, delay, seed):
    rng = np.random.default_rng(seed)
    a = np.exp(-1 / tau)
    x = np.zeros(n)
    e = rng.standard_normal(n)
    for i in range(1, n):
        x[i] = a * x[i - 1] + e[i]
    y = np.roll(x, delay) + 0.6 * rng.standard_normal(n) * x.std()
    return (x - x.mean()) / x.std(), (y - y.mean()) / y.std()


def ease(t):
    t = np.clip(t, 0, 1)
    return t * t * (3 - 2 * t)


def heading(fig, text, sub=None):
    fig.text(0.06, 0.9, text, fontsize=30, fontweight="bold", color=INK)
    if sub:
        fig.text(0.06, 0.845, sub, fontsize=17, color=INK2)


class Scene:
    """Draws frame k of n into a fresh figure via draw(fig, p) with p in [0, 1]."""

    def __init__(self, draw):
        self.draw = draw


def s_title(fig, p):
    fig.text(0.5, 0.58, "Adaptive TDS", ha="center", fontsize=58, fontweight="bold", alpha=ease(p * 4))
    fig.text(0.5, 0.47, "how it chooses its own window, lag and tolerance", ha="center", fontsize=22, color=INK2, alpha=ease(p * 4 - 0.5))
    ax = fig.add_axes([0.2, 0.18, 0.6, 0.18])
    ax.axis("off")
    t = np.linspace(0, 6 * np.pi, 600)
    k = int(600 * ease(p * 1.5))
    ax.plot(t[:k], np.sin(t[:k]), color=TDS, lw=3)
    ax.plot(t[:k], np.sin(t[:k] - 0.9) - 2.6, color=ATDS, lw=3)
    ax.set_xlim(0, 6 * np.pi)
    ax.set_ylim(-4, 1.4)


X2, Y2 = ar_pair(900, 4, 6, 1)


def s_tds(fig, p):
    heading(fig, "Classic TDS", "slide a window · find the best delay τ₀ · stable if 4 of 5 agree")
    ax = fig.add_axes([0.07, 0.42, 0.88, 0.35])
    t = np.arange(900)
    ax.plot(t, X2 + 2.2, color=INK2, lw=1.2)
    ax.plot(t, Y2 - 2.2, color=AQUA, lw=1.2)
    ax.set_yticks([])
    ax.set_xlim(0, 900)
    ax.set_xlabel("time (s)")
    for s in ax.spines.values():
        s.set_visible(False)
    L, step = 60, 30
    nwin = int(1 + (900 - L) / step)
    k = int(nwin * ease(p * 1.15))
    ax.axvspan(k * step, k * step + L, color=TDS, alpha=0.25)
    ax.text(10, 4.1, "signal 1", color=INK2, fontsize=13)
    ax.text(10, -0.4, "signal 2", color=AQUA, fontsize=13)
    bx = fig.add_axes([0.07, 0.08, 0.88, 0.25])
    taus = []
    for i in range(nwin):
        a, b = X2[i * step:i * step + L], Y2[i * step:i * step + L]
        lags = np.arange(-30, 31)
        c = [np.corrcoef(a, np.roll(b, -l))[0, 1] for l in lags]
        taus.append(lags[int(np.argmax(np.abs(c)))])
    taus = np.array(taus)
    stable = np.zeros(nwin, bool)
    for i in range(nwin - 4):
        seg = taus[i:i + 5]
        for d in np.unique(seg):
            m = np.abs(seg - d) <= 1
            if m.sum() >= 4:
                stable[i:i + 5][m] = True
                break
    cx = np.arange(nwin) * step + L / 2
    bx.scatter(cx[:k + 1], taus[:k + 1], s=60, facecolors=np.where(stable[:k + 1], TDS, BG), edgecolors=TDS, linewidths=2, zorder=3)
    bx.set_xlim(0, 900)
    bx.set_ylim(-32, 32)
    bx.set_ylabel("τ₀ (s)")
    bx.axhline(6, color=TDS, lw=0.8, ls=":")
    bx.text(895, 10, "true delay", va="bottom", ha="right", fontsize=12, color=INK2)
    bx.grid(color=GRID)
    bx.text(0, 36, "filled = stable window", fontsize=13, color=INK2)


XF, _ = ar_pair(600, 1.5, 0, 3)
XS, _ = ar_pair(600, 40, 0, 4)


def s_problem(fig, p):
    heading(fig, "The problem with fixed settings", "60-s window, ±1-s tolerance, whatever the signal")
    for row, (x, name, indep, col) in enumerate([(XF, "fast signal (e.g. muscle activity)", "≈ 16 independent samples", AQUA),
                                                 (XS, "slow signal (e.g. EEG δ power, heart rate)", "< 1 independent sample", ATDS)]):
        ax = fig.add_axes([0.07, 0.47 - row * 0.33, 0.52, 0.24])
        ax.plot(np.arange(600), x, color=col, lw=1.4)
        ax.axvspan(270, 330, color=TDS, alpha=0.3)
        ax.set_yticks([])
        ax.set_xlim(0, 600)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_title(name, loc="left", fontsize=15, color=INK2)
        if p > 0.25 + row * 0.25:
            fig.text(0.64, 0.58 - row * 0.33, "in one 60-s window:", fontsize=15, color=INK2)
            fig.text(0.64, 0.53 - row * 0.33, indep, fontsize=21, fontweight="bold", color=col)
    if p > 0.75:
        fig.text(0.07, 0.06, "Real data, Bashan 2012 cohort: EEG δ 0.36, θ 0.82, heart rate 0.93 independent samples per minute",
                 fontsize=14, color=INK2)


def make_acf_scene(ex):
    d = ex["sleep"]["atds_derivation"]
    name, a = next(iter(d["acf"].items()))
    a = np.array(a)
    B = d["bartlett_factor"][name]
    zc = int(np.argmax(a <= 0)) if (a <= 0).any() else len(a)

    def s_memory(fig, p):
        heading(fig, "Step 1 · measure each signal's memory", f"autocorrelation of {name} (real night, subject #29)")
        ax = fig.add_axes([0.08, 0.17, 0.55, 0.58])
        k = max(2, int(len(a) * ease(p * 1.6)))
        lags = np.arange(len(a))
        ax.plot(lags[:k], a[:k], color=ATDS, lw=2.5)
        ax.axhline(0, color=INK2, lw=0.8)
        if p > 0.45:
            ax.fill_between(lags[:min(k, zc)], 0, a[:min(k, zc)], color=ATDS, alpha=0.25)
            if zc < len(a):
                ax.axvline(zc, color=INK2, ls=":", lw=1)
                ax.text(zc + 5, 0.8, "first zero\ncrossing", fontsize=13, color=INK2)
            else:
                ax.text(len(a) - 8, 0.85, "no zero crossing yet:\nthe sum continues\nbeyond 600 s", fontsize=13, color=INK2, ha="right")
        ax.set_xlim(0, len(a))
        ax.set_ylim(-0.4, 1.05)
        ax.set_xlabel("lag (s)")
        ax.set_ylabel("autocorrelation r(k)")
        ax.grid(color=GRID)
        if p > 0.6:
            fig.text(0.68, 0.62, "Bartlett factor", fontsize=18, color=INK2)
            fig.text(0.68, 0.54, "B = 1 + 2 Σ r(k)²", fontsize=26, fontweight="bold")
            fig.text(0.68, 0.44, f"{name}: B ≈ {B:.0f}", fontsize=21, color=ATDS)
            fig.text(0.68, 0.33, f"{B:.0f} samples ≈ one\nindependent sample", fontsize=16, color=INK2)
    return s_memory


def make_window_scene(ex):
    rw = ex["sleep"]["atds_derivation"]["raw_window_by_signal"]
    sig = sorted(rw, key=lambda s: -rw[s])

    def s_window(fig, p):
        heading(fig, "Step 2 · window = 30 independent samples", "of the slowest signal, clipped at the cap you choose")
        ax = fig.add_axes([0.08, 0.14, 0.86, 0.62])
        k = int(len(sig) * ease(p * 1.6))
        vals = [rw[s] for s in sig]
        ax.bar(range(k), vals[:k], color=[ATDS if i == 0 else NULL for i in range(k)], width=0.65)
        ax.set_yscale("log")
        ax.set_ylim(30, 10000)
        ax.set_xticks(range(len(sig)))
        ax.set_xticklabels(sig)
        ax.set_xlim(-0.6, len(sig) - 0.4)
        ax.set_ylabel("window for 30 independent samples (s)")
        ax.grid(axis="y", color=GRID)
        ax.axhline(60, color=TDS, lw=2, ls="--")
        fig.text(0.1, 0.79, "dashed: 60 s (classic TDS window)", color=TDS, fontsize=14)
        if p > 0.55:
            ax.axhline(300, color=ATDS, lw=2.5)
            ax.text(len(sig) - 0.5, 330, "cap: 300 s → aTDS window", ha="right", color=ATDS, fontsize=15, fontweight="bold")
    return s_window


def s_lag(fig, p):
    heading(fig, "Step 3 · step and lag range follow the window")
    ax = fig.add_axes([0.08, 0.25, 0.84, 0.45])
    ax.set_xlim(0, 900)
    ax.set_ylim(0, 3)
    ax.axis("off")
    g = ease(p * 1.4)
    for row, (L, col, lab) in enumerate([(60, TDS, "classic TDS: window 60 s · step 30 s · lag ±30 s"),
                                         (300, ATDS, "aTDS: window 300 s · step 150 s · lag ±149 s")]):
        y = 2.1 - row * 1.4
        for j in range(3):
            x0 = 120 + j * L / 2 * (1 + 2 * row * 0) + (j * L / 2 if row else 0) * 0
            x0 = 120 + j * L / 2
            ax.add_patch(plt.Rectangle((x0, y - 0.18 + 0.12 * j), L, 0.3, fill=False, ec=col, lw=2.5, alpha=1 if g > j / 3 else 0.15))
        if g > 0.6:
            c = 120 + L / 2
            ax.annotate("", xy=(c - L / 2, y + 0.55), xytext=(c + L / 2, y + 0.55), arrowprops=dict(arrowstyle="<->", color=col, lw=2))
        ax.text(120, y - 0.55, lab, fontsize=17, color=col)
    if p > 0.7:
        fig.text(0.08, 0.1, "Lag search: up to half the window, never past half of a detected rhythm.",
                 fontsize=16, color=INK2)


def make_tol_scene(ex):
    nb = ex["sleep"]["atds_derivation"]["null_score_by_tolerance"]
    tol = np.array(sorted(map(int, nb)))
    val = np.array([nb[str(t)] for t in tol])
    chosen = ex["sleep"]["atds_derivation"]["chosen_tolerance"]

    def s_tol(fig, p):
        heading(fig, "Step 4 · tolerance from a coupling-free test", "rotate one signal in time → coupling gone, own shape kept")
        ax = fig.add_axes([0.08, 0.14, 0.55, 0.6])
        k = max(2, int(len(tol) * ease(p * 1.5)))
        ax.plot(tol[:k], val[:k], color=ATDS, lw=2.5, marker="o")
        ax.axhline(5, color=INK2, ls="--")
        ax.text(0.3, 5.25, "5 % chance budget", color=INK2, fontsize=14)
        if k >= len(tol) - 1 and p > 0.65:
            i = list(tol).index(chosen)
            ax.scatter([chosen], [val[i]], s=260, color=ATDS, zorder=4, ec=INK, lw=2)
            ax.annotate(f"chosen: ±{chosen} s", (chosen, val[i]), xytext=(chosen - 7, val[i] + 1.3), color=ATDS, fontsize=17,
                        fontweight="bold", arrowprops=dict(arrowstyle="->", color=ATDS))
        ax.set_xlim(-0.5, tol.max() + 0.8)
        ax.set_ylim(0, 7)
        ax.set_xlabel("tolerance (± s)")
        ax.set_ylabel("fake-pair TDS score (%)")
        ax.grid(color=GRID)
        rx = fig.add_axes([0.68, 0.32, 0.27, 0.38])
        rx.axis("off")
        t = np.linspace(0, 1, 200)
        sh = int(200 * min(1, p * 2) * 0.45)
        rx.plot(t, np.sin(14 * t) * np.exp(-t) + 1.4, color=INK2, lw=2)
        rx.plot(t, np.roll(np.sin(14 * t - 1) * np.exp(-t), sh) - 0.8, color=AQUA, lw=2)
        rx.text(0, 2.6, "signal 1", color=INK2, fontsize=13)
        rx.text(0, 0.35, "signal 2, rotated", color=AQUA, fontsize=13)
        rx.set_ylim(-2.2, 3)
    return s_tol


def s_summary(fig, p):
    heading(fig, "What changes", "one parameter set per recording, from the signals themselves")
    rows = [("window", "60 s", "300 s"), ("step", "30 s", "150 s"), ("lag search", "±30 s", "±149 s"), ("tolerance", "±1 s", "±15 s"),
            ("rule", "4 of 5 windows", "4 of 5 windows")]
    fig.text(0.3, 0.72, "classic TDS", fontsize=20, color=TDS, fontweight="bold", ha="center")
    fig.text(0.62, 0.72, "aTDS (subject #29)", fontsize=20, color=ATDS, fontweight="bold", ha="center")
    for i, (a, b, c) in enumerate(rows):
        if p > i / 7:
            y = 0.63 - i * 0.075
            fig.text(0.08, y, a, fontsize=19, color=INK2)
            fig.text(0.3, y, b, fontsize=21, ha="center")
            fig.text(0.62, y, c, fontsize=21, ha="center", color=ATDS if b != c else INK)
    if p > 0.7:
        fig.text(0.08, 0.14, "Uses only each signal's own memory and a coupling-free test, so it cannot be tuned toward the answer.",
                 fontsize=16, color=INK2)


def s_when(fig, p):
    heading(fig, "When it helps · what it costs")
    good = ["fixed window ≪ the signals' memory", "slow physiology, rivers, climate, running", "long recordings",
            "different sampling rates give the same answer"]
    bad = ["coarser time resolution (set the cap)", "less precise delays (wider tolerance)", "short states / short records lose windows",
           "rhythmic raw waveforms: use beat-to-beat series"]
    fig.text(0.08, 0.72, "helps", fontsize=24, color=ATDS, fontweight="bold")
    fig.text(0.53, 0.72, "costs", fontsize=24, color=TDS, fontweight="bold")
    for i, (g, b) in enumerate(zip(good, bad)):
        if p > i / 6:
            fig.text(0.08, 0.63 - i * 0.09, "• " + g, fontsize=16)
        if p > (i + 1) / 6:
            fig.text(0.53, 0.63 - i * 0.09, "• " + b, fontsize=16)
    if p > 0.85:
        fig.text(0.08, 0.12, "Evidence: aTDS reports 01–07, ground-truth study, aTDS paper draft (pyTDS)", fontsize=13, color=INK2)


def render(draw, dur, path):
    n = max(2, int(round(dur * FPS)))
    fig = plt.figure(figsize=(W, H), dpi=100)
    writer = FFMpegWriter(fps=FPS, codec="libx264", extra_args=["-pix_fmt", "yuv420p", "-preset", "veryfast", "-crf", "20"])
    with writer.saving(fig, path, dpi=100):
        for k in range(n):
            fig.clf()
            draw(fig, k / (n - 1))
            writer.grab_frame(facecolor=BG)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    for k in ("--examples", "--voice", "--work", "--out"):
        ap.add_argument(k, required=True)
    a = ap.parse_args()
    ex = json.load(open(a.examples))
    draws = {"title": s_title, "tds": s_tds, "problem": s_problem, "memory": make_acf_scene(ex), "window": make_window_scene(ex),
             "lag": s_lag, "tolerance": make_tol_scene(ex), "summary": s_summary, "when": s_when}
    os.makedirs(a.work, exist_ok=True)
    parts = []
    for i, (key, text) in enumerate(SCENES):
        wav = os.path.join(a.work, f"{i:02d}_{key}.wav")
        dur = say(text, a.voice, wav) + 0.8
        vid = os.path.join(a.work, f"{i:02d}_{key}_v.mp4")
        render(draws[key], dur, vid)
        out = os.path.join(a.work, f"{i:02d}_{key}.mp4")
        subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-i", vid, "-i", wav, "-af", "adelay=300|300,apad", "-c:v", "copy",
                        "-c:a", "aac", "-b:a", "128k", "-ar", "44100", "-ac", "1", "-shortest", out], check=True)
        parts.append(out)
        print(f"{key}: {dur:.1f} s")
    lst = os.path.join(a.work, "list.txt")
    open(lst, "w").write("".join(f"file '{p}'\n" for p in parts))
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-f", "concat", "-safe", "0", "-i", lst, "-c", "copy", a.out], check=True)
    print(a.out)


if __name__ == "__main__":
    main()
