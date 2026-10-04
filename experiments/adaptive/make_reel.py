"""Render the adaptive-TDS Instagram reel (1080x1920, 30 fps, H.264, no audio).

Scene length = narration + 0.3 s lead-in + 0.9 s tail (min 4 s); visuals and subtitles are synchronised to the
sentence times in reel_audio/sentences.json (written by make_voice.py).

Run from anywhere:  python experiments/adaptive/make_reel.py
Output: "adaptive-TDS reports/adaptive_TDS_reel.mp4" plus one still PNG per scene
in "adaptive-TDS reports/figures/reel_frames/".
"""
import json
import subprocess
import sys
from functools import lru_cache
from multiprocessing import Pool
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Circle, FancyBboxPatch, Rectangle
from matplotlib.textpath import TextPath

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from pyTDS import TDSParams, tds  # noqa: E402

sys.path.insert(0, str(ROOT / "experiments/adaptive"))
import exp_fields as EF  # noqa: E402
from adaptive import acf  # noqa: E402

RES = ROOT / "experiments/adaptive/results"
OUT = ROOT / "adaptive-TDS reports/adaptive_TDS_reel.mp4"
FRAMES = ROOT / "adaptive-TDS reports/figures/reel_frames"
FIELDS = json.load(open(RES / "fields.json"))
EXC = json.load(open(RES / "fields_excerpts.json"))
SENT = json.load(open(ROOT / "experiments/adaptive/reel_audio/sentences.json"))
LEAD, TAIL = 0.3, 0.9  # narration starts 0.3 s into each scene; 0.9 s hold after it

W, H, FPS = 1080, 1920, 30
BG, TX, TX2, GRID = "#1a1a19", "#ffffff", "#c3c2b7", "#3a3a38"
BLUE, ORANGE, GRAY, AQUA = "#3987e5", "#d95926", "#8a8984", "#199e70"
STATUS = {"better": ("✓", AQUA), "mixed": ("≈", "#c98500"), "tie": ("≈", "#c98500"),
          "worse": ("✗", "#e66767")}
plt.rcParams["font.family"] = "DejaVu Sans"
PT = 0.72  # px -> pt at dpi 100
FIELD_DUR = 5.5

# ---------------------------------------------------------------- helpers
class Ctx:
    def __init__(self, fig, ax, sec, sid=None):
        self.fig, self.ax, self.t = fig, ax, sec  # t = seconds into scene
        sc = SENT.get(sid, {"sentences": []})["sentences"]
        self.S = [(s["start"] + LEAD, s["end"] + LEAD) for s in sc]  # sentence (start, end) in scene time
        self.chunks = [(k["start"] + LEAD, k["text"]) for s in sc for k in s["chunks"]]


def S(c, i):
    """Scene time at which sentence i starts being spoken."""
    return c.S[min(i, len(c.S) - 1)][0]


def E(c, i):
    """Scene time at which sentence i ends."""
    return c.S[min(i, len(c.S) - 1)][1]


def fade(c, t0, d=0.4):
    return float(np.clip((c.t - t0) / d, 0, 1))


def prog(c, t0, t1):
    return float(np.clip((c.t - t0) / (t1 - t0), 0, 1))


@lru_cache(maxsize=None)
def _w(s, px, bold):
    fp = FontProperties(family="DejaVu Sans", weight="bold" if bold else "normal")
    return TextPath((0, 0), s, size=px, prop=fp).get_extents().width


def wrap(s, px, bold, maxw):
    lines, cur = [], ""
    for word in s.split(" "):
        trial = (cur + " " + word).strip()
        if cur and _w(trial, px, bold) > maxw:
            lines.append(cur)
            cur = word
        else:
            cur = trial
    return lines + [cur]


def T(c, x, y, s, px, col=TX, bold=False, a=1.0, ha="left", maxw=None, lh=1.28):
    """Draw (wrapped) text, top-left anchored at y; return y below the block."""
    lines = wrap(s, px, bold, maxw) if maxw else [s]
    for i, ln in enumerate(lines):
        if a > 0:
            c.ax.text(x, y + i * px * lh, ln, fontsize=px * PT, color=col, alpha=a, ha=ha, va="top",
                      fontweight="bold" if bold else "normal")
    return y + len(lines) * px * lh


def box(c, x, y, w, h, fc="none", ec="none", a=1.0, r=20, lw=2):
    if a <= 0 or w <= 0 or h <= 0:
        return
    r = min(r, w / 2, h / 2)
    c.ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}", fc=fc,
                                  ec=ec, lw=lw * PT, alpha=a))


def sub_axes(c, x, y, w, h, a):
    ax = c.fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H], facecolor=BG)
    for sp in ax.spines.values():
        sp.set_color(GRID)
        sp.set_alpha(a)
    ax.set_xticks([])
    ax.set_yticks([])
    return ax


def lines_plot(ax, series, p, a, lw=2.5):
    """Draw each (x, y, colour) progressively left to right up to fraction p."""
    for x, y, col in series:
        k = max(2, int(len(x) * p))
        ax.plot(x[:k], y[:k], color=col, lw=lw * PT, alpha=a, solid_capstyle="round")


def status(s, col_alpha=1.0):
    return STATUS[s]


def readable(sec):
    """Seconds -> short human string."""
    sec = abs(sec)
    for lim, div, u in [(60, 1, "s"), (3600, 60, "min"), (86400, 3600, "h"), (7 * 86400, 86400, "days"),
                        (28 * 86400, 7 * 86400, "weeks"), (45 * 86400, 86400, "days"),
                        (730 * 86400, 30.44 * 86400, "months")]:
        if sec < lim:
            return f"{sec / div:.3g} {u}"
    return f"{sec / (365.25 * 86400):.3g} years"


def minus(s):
    return s.replace("-", "−")



def ease(p):
    return float(p * p * (3 - 2 * p))


DISPLAY = [("T D S", "TDS"), ("E E G", "EEG"), ("El Nino", "El Niño")]


def shown(s):
    for a, b in DISPLAY:
        s = s.replace(a, b)
    return s


def subtitle(c):
    """Current narration chunk as a 1-2 line caption at the bottom of the safe zone (y 1418-1540)."""
    cur = None
    for t0, txt in c.chunks:
        if c.t >= t0 - 0.05:
            cur, ct = shown(txt), t0
    if cur is None:
        return
    a = fade(c, ct - 0.05, 0.15)
    mw = 920
    lines = wrap(cur, 38, False, mw)
    if len(lines) > 2:
        print("WARNING: subtitle wraps to", len(lines), "lines:", cur, flush=True)
    elif len(lines) == 2 and len(lines[1].split()) == 1:  # avoid a one-word second line: balance the two lines
        while mw > 500 and len(wrap(cur, 38, False, mw - 20)) == 2:
            mw -= 20
    box(c, 40, 1418, 1000, 122, "#0d0d0c", a=0.92 * a, r=18)
    T(c, 540, 1479 - len(lines) * 38 * 1.25 / 2, cur, 38, TX, a=a, ha="center", maxw=mw, lh=1.25)


def plot_axes(c, x, y, w, h, a, xlim, ylim, xt=None, yt=None, xl="", yl="", fs=28):
    """Axes at pixel box (x, y, w, h) with real tick labels; returns None while invisible."""
    if a <= 0:
        return None
    ax = c.fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H], facecolor=BG)
    for sp in ax.spines.values():
        sp.set_color(GRID)
        sp.set_alpha(a)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xticks(xt if xt is not None else [])
    ax.set_yticks(yt if yt is not None else [])
    ax.tick_params(colors=TX2, labelsize=fs * PT, length=6, width=1.5 * PT)
    for lb in ax.get_xticklabels() + ax.get_yticklabels():
        lb.set_alpha(a)
    ax.set_xlabel(xl, color=TX2, fontsize=(fs + 4) * PT, labelpad=10, alpha=a)
    ax.set_ylabel(yl, color=TX2, fontsize=(fs + 4) * PT, labelpad=10, alpha=a)
    return ax


def at(ax, x, y, s, px, col=TX, a=1.0, ha="left", va="top", bold=False):
    """Text in axes-fraction coordinates."""
    if a > 0:
        ax.text(x, y, s, transform=ax.transAxes, fontsize=px * PT, color=col, alpha=a, ha=ha, va=va,
                fontweight="bold" if bold else "normal")


def hdr(c, text, t0=0.1, maxw=960):
    """Scene header; returns y below it."""
    y = 270
    for ln in text.split("\n"):
        y = T(c, 60, y, ln, 64, bold=True, a=fade(c, t0), maxw=maxw, lh=1.15)
    return y


# ---------------------------------------------------------------- hydrology data (real)
HYD = EF.d_hydro_hourly()
ACF_RAIN, ACF_RIVER = acf(HYD["s1"], 240), acf(HYD["s2"], 240)
HV = FIELDS["hydrology_hourly"]["variants"]["adaptive"]
HB = max(HV["info"]["bartlett"])  # 63.8 hours per independent sample
HPS, HPP = HV["params_samples"], HV["params_phys"]
NULL = {int(k): v for k, v in HV["info"]["null_score_by_tolerance"].items()}
TOL = HPS["tolerance"]
HEX = EXC["hydrology_hourly"]
RAIN_Z, RIVER_Z = np.array(HEX["s1_z"], float), np.array(HEX["s2_z"], float)


def band(x, top):
    """Normalise a series to [0,1] and place it in a band starting at `top`."""
    return top + (x - x.min()) / (x.max() - x.min())


# ---------------------------------------------------------------- scene 1
def scene_hook(c):
    T(c, W / 2, 560, "Can an algorithm tune itself?", 84, bold=True, a=fade(c, 0.3), ha="center", maxw=960, lh=1.2)
    T(c, W / 2, 800, "Adaptive TDS · tested on 9 fields of science", 40, TX2, a=fade(c, S(c, 1)), ha="center", maxw=960)
    ax = sub_axes(c, 60, 1100, 960, 260, 0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    x = np.linspace(0, 1, 400)
    f = lambda u: np.sin(2 * np.pi * 3 * u) + 0.5 * np.sin(2 * np.pi * 7.3 * u)
    p = prog(c, S(c, 1), S(c, 1) + 2.0)
    lines_plot(ax, [(x, f(x), BLUE), (x, f(x - 0.04), ORANGE)], p, fade(c, S(c, 1)))
    ax.set_ylim(-2, 2)


# ---------------------------------------------------------------- scene 2
def _demo():
    rng = np.random.default_rng(3)
    n = 400
    e = rng.standard_normal(n + 7)
    x = np.zeros(n + 7)
    for i in range(1, n + 7):
        x[i] = 0.9 * x[i - 1] + e[i]
    s1 = x[7:]
    s2 = x[:n] + 0.15 * x.std() * rng.standard_normal(n)  # s2[t] = s1[t-7] + noise
    return s1, s2, tds(s1, s2, TDSParams(window=60))


DEMO = _demo()


def scene_what(c):
    s1, s2, r = DEMO
    n = len(s1)
    T(c, 60, 270, "What TDS does", 64, bold=True, a=fade(c, 0.1))
    X0, PW, PY, PH = 60, 960, 380, 400
    ax = sub_axes(c, X0, PY, PW, PH, fade(c, S(c, 1) - 0.3, 0.3))
    idx = np.arange(n)
    # sentences: 0 title, 1 slides a window, 2 finds the delay, 3 counts stability, 4 the score
    lines_plot(ax, [(idx, s1 / s1.std() + 2.2, BLUE), (idx, s2 / s2.std() - 2.2, ORANGE)],
               prog(c, S(c, 1) - 0.3, S(c, 1) + 1.0), fade(c, S(c, 1) - 0.3, 0.3))
    ax.set_xlim(0, n)
    ax.set_ylim(-4.8, 4.8)
    q = prog(c, S(c, 1) + 0.8, S(c, 4) + 0.2)
    centre = 30 + q * (n - 60)
    on = c.t > S(c, 1) + 0.8
    if on:
        ax.add_patch(Rectangle((centre - 30, -4.7), 60, 9.4, fc="#ffffff", alpha=0.14, ec=TX, lw=1.5 * PT))
    k = int(np.argmin(abs(r["t_vec"] - centre)))
    if c.t > S(c, 2):
        T(c, 60, 840, f"delay in this window: {minus(str(int(r['tau'][k])))}", 44, bold=True, a=fade(c, S(c, 2), 0.3))
    nd = len(r["tau"])
    if c.t > S(c, 3):
        for i in range(nd):
            if r["t_vec"][i] <= centre + 1e-6:
                c.ax.add_patch(Circle((60 + 30 + i * 80, 990), 24, fc=AQUA if r["stbl_lbl"][i] else GRAY, ec="none"))
        T(c, 60, 1050, "aqua = stable    gray = unstable", 30, TX2, a=fade(c, S(c, 3), 0.3))
    if c.t > S(c, 4):
        sc = r["score"]
        a = fade(c, S(c, 4), 0.3)
        T(c, 60, 1120, f"{int(np.sum(r['stbl_lbl']))} of {nd} windows are stable", 40, TX2, a=a, maxw=960)
        T(c, 60, 1180, f"TDS score = {sc:.0f}%", 64, AQUA, True, a=a, maxw=960)


# ---------------------------------------------------------------- scene 3
def scene_problem(c):
    T(c, 60, 270, "Fixed settings break", 64, bold=True, a=fade(c, 0.1))
    T(c, 60, 370, "Same night of sleep EEG — only the sampling rate changed", 40, TX2, a=fade(c, S(c, 1)), maxw=960)
    vals, labs = [74.7, 58.1, 41.3, 26.8], ["0.5×", "1×", "2×", "4×"]
    # bars appear one by one while the narration runs (sentences: 0 fixed numbers, 1 one night, 2 same brain, 3 score drops)
    ts = [S(c, 1) + 0.6, S(c, 1) + 0.5 * (E(c, 1) - S(c, 1)) + 0.6, S(c, 2), S(c, 3)]
    base, top, bw = 1070, 560, 150
    for i, (v, lb) in enumerate(zip(vals, labs)):
        x = 90 + i * 235
        a = fade(c, ts[i])
        hh = (base - top) * v / 80 * prog(c, ts[i], ts[i] + 1.0)
        box(c, x, base - hh, bw, hh, GRAY, a=a, r=24)
        if hh > 5:
            T(c, x + bw / 2, base - hh - 56, f"{v}%", 44, bold=True, a=a, ha="center")
        T(c, x + bw / 2, base + 20, lb, 44, TX2, a=a, ha="center")
    c.ax.plot([60, 1020], [base, base], color=GRID, lw=2 * PT, alpha=fade(c, ts[0]))
    T(c, 60, 1160, "Original TDS uses fixed sample counts (window = 60 samples).", 40, a=fade(c, S(c, 0)), maxw=960)
    T(c, 60, 1285, "Change the recording device → different answer.", 40, TX, bold=True, a=fade(c, S(c, 3) + 1.2), maxw=960)


# ---------------------------------------------------------------- scene 4: how_intro
def scene_how_intro(c):
    y = hdr(c, "How it tunes itself —\n4 steps", S(c, 0))
    steps = ["1  Memory", "2  Window", "3  Max delay", "4  Tolerance"]
    for i, s in enumerate(steps):
        a = fade(c, S(c, 0) + 0.7 + 0.45 * i)
        x, yy = 60 + (i % 2) * 490, y + 30 + (i // 2) * 110
        box(c, x, yy, 470, 90, "#232322", AQUA, a=a, lw=3)
        T(c, x + 235, yy + 20, s, 44, TX, True, a=a, ha="center")
    py = y + 280
    a = fade(c, S(c, 1))
    ax = plot_axes(c, 60, py, 960, 1385 - py - 20, a, (0, 69), (-0.05, 2.55))
    if ax is not None:
        lo, hi = 432, 502
        rain, river = RAIN_Z[lo:hi], RIVER_Z[lo:hi]
        xs = np.arange(len(rain))
        pr = prog(c, S(c, 1) + 0.4, E(c, 2) - 0.3)
        ra, rv = band(rain, 1.45), band(river, 0.0)
        k = max(2, int(len(xs) * pr))
        ax.plot(xs[:k], ra[:k], color=BLUE, lw=4 * PT)
        ax.plot(xs[:k], rv[:k], color=ORANGE, lw=4 * PT)
        at(ax, 0.02, 0.97, "rain", 36, BLUE, bold=True, a=fade(c, S(c, 1) + 0.4))
        at(ax, 0.02, 0.43, "river flow", 36, ORANGE, bold=True, a=fade(c, S(c, 1) + 0.4))
        at(ax, 0.98, 0.02, "time (hours)  →", 30, TX2, ha="right", va="bottom", a=fade(c, S(c, 1) + 0.4))
        if c.t > S(c, 2) + 1.0:
            ap = fade(c, S(c, 2) + 1.0, 0.5)
            pr_i, pv_i = int(np.argmax(rain)), int(np.argmax(river))
            ax.annotate("", xy=(pv_i, rv[pv_i] + 0.04), xytext=(pr_i + 1, ra[pr_i] - 0.12),
                        arrowprops=dict(arrowstyle="-|>", color=AQUA, lw=4 * PT, alpha=ap, mutation_scale=30))
            ax.text(max(pr_i, pv_i) + 4, 1.2, "rain → river", color=AQUA, fontsize=40 * PT, fontweight="bold",
                    alpha=ap, ha="left", va="center")


# ---------------------------------------------------------------- scene 5: step 1 memory
def scene_step1(c):
    hdr(c, "STEP 1 · Memory", S(c, 0))
    # sentence 1: a signal compared with a shifted copy of itself
    a = fade(c, S(c, 1))
    ax = plot_axes(c, 60, 370, 960, 250, a, (0, 150), (-0.15, 2.6))
    if ax is not None:
        v = band(RIVER_Z[440:590], 0.0) * 0.9
        sh = 28
        p = prog(c, S(c, 1) + 0.3, S(c, 1) + 1.8)
        xs = np.arange(150)
        ax.plot(xs, v + 1.55, color=ORANGE, lw=4 * PT)
        ax.plot(xs + sh * ease(p), v + 0.0, color=GRAY, lw=4 * PT, ls=(0, (5, 3)))
        at(ax, 0.99, 0.97, "river flow", 32, ORANGE, bold=True, ha="right")
        at(ax, 0.99, 0.50, "copy of itself, shifted", 32, TX2, bold=True, ha="right")
    # ACF plot
    a = fade(c, S(c, 2) - 0.4)
    ax = plot_axes(c, 160, 700, 860, 560, a, (0, 240), (0, 1.02), xt=[0, 48, 96, 144, 192, 240], yt=[0, 0.5, 1],
                   xl="shift (hours)", yl="similarity to itself")
    if ax is not None:
        lags = np.arange(241)
        p1 = prog(c, S(c, 2), S(c, 2) + 2.0)
        p2 = prog(c, S(c, 4), S(c, 4) + 2.0)
        if p1 > 0:
            k = max(2, int(241 * p1))
            ax.plot(lags[:k], ACF_RAIN[:k], color=BLUE, lw=5 * PT)
        if p2 > 0:
            k = max(2, int(241 * p2))
            ax.plot(lags[:k], ACF_RIVER[:k], color=ORANGE, lw=5 * PT)
        at(ax, 0.30, 0.97, "rain forgets in hours", 38, BLUE, bold=True, a=fade(c, S(c, 2) + 0.6))
        at(ax, 0.30, 0.86, "river remembers for days", 38, ORANGE, bold=True, a=fade(c, S(c, 4) + 0.6))


# ---------------------------------------------------------------- scene 6: step 2 window
def scene_step2(c):
    hdr(c, "STEP 2 · Window", S(c, 0))
    n, pitch, sq = 64, 15, 12
    # sentence 1: 64 hourly points, neighbours are almost copies
    a1 = fade(c, S(c, 1))
    if a1 > 0:
        T(c, 60, 360, "slow signal: neighbouring points ≈ copies", 36, TX2, a=a1, maxw=960)
    merge = ease(prog(c, S(c, 2), S(c, 2) + 1.2))
    for i in range(n):
        ai = fade(c, S(c, 1) + (E(c, 1) - S(c, 1)) * 0.8 * i / n, 0.15)
        if ai <= 0 or merge >= 1:
            continue
        x0, y0 = 60 + i * pitch, 440
        x1, y1 = 465 + 150 * (i + 0.5) / n - sq / 2, 560 + 12 * ((i * 7) % 5)
        x, y = x0 + (x1 - x0) * merge, y0 + (y1 - y0) * merge
        c.ax.add_patch(Rectangle((x, y), sq, sq, fc=ORANGE, ec="none", alpha=ai))
    blk = fade(c, S(c, 2) + 1.1, 0.25)
    # slot geometry of the 30-block bar: 1 hour = 960/window px
    win = HPS["window"]
    sc = 960 / win
    slot = HB * sc
    by = 790
    t3 = S(c, 3)
    mv = ease(prog(c, t3, t3 + 0.8))
    if blk > 0 and c.t < t3 + 0.8:
        bx, bw_, bh = 465, 150, 62
        x = bx + (60 - bx) * mv
        y = 585 + (by - 585) * mv
        w_, h_ = bw_ + (slot - 2 - bw_) * mv, bh + (80 - bh) * mv
        box(c, x, y, w_, h_, ORANGE, a=blk, r=8)
        T(c, W / 2, 665, f"1 new point = {HB:.0f} h", 44, TX, True, a=blk * (1 - mv), ha="center")
    if c.t >= t3:
        n_show = int(np.clip(1 + (c.t - t3 - 0.8) / ((E(c, 3) - t3 - 0.8) / 29), 1, 30)) if c.t > t3 + 0.8 else 1
        for j in range(n_show):
            box(c, 60 + j * slot, by, slot - 2, 80, ORANGE, a=1.0 if j else fade(c, t3), r=4)
        T(c, W / 2, 890, f"30 new points = 30 × {HB:.0f} h", 44, TX, True, a=fade(c, t3 + 0.8), ha="center", maxw=960)
    t4 = S(c, 4)
    if c.t > t4:
        T(c, W / 2, 970, f"≈ {win / 24:.0f} days", 72, AQUA, True, a=fade(c, t4), ha="center")
        T(c, W / 2, 1065, f"window = {win} h = {HPP['window']}", 34, TX2, a=fade(c, t4 + 0.2), ha="center")
        a = fade(c, t4 + 0.9)
        fx = 60 + 60 * sc
        box(c, 60, 1170, max(fx - 60, 3), 50, GRAY, a=a, r=3)
        T(c, 60 + max(fx - 60, 3) + 20, 1175, "fixed: 60 h", 40, TX2, True, a=a)
        T(c, 60, 1235, "(same scale)", 30, TX2, a=a)


# ---------------------------------------------------------------- scene 7: step 3 delay
def pulse(u):
    """Heartbeat-like pulse train, period 1."""
    ph = (u % 1.0)
    return np.exp(-((ph - 0.25) / 0.035) ** 2) - 0.35 * np.exp(-((ph - 0.42) / 0.06) ** 2) + 0.45 * np.exp(-((ph - 0.7) / 0.07) ** 2)


def scene_step3(c):
    hdr(c, "STEP 3 · Max delay", S(c, 0))
    a = fade(c, S(c, 0) + 0.2)
    box(c, 60, 440, 960, 90, AQUA, a=a, r=14)
    T(c, W / 2, 462, f"window ≈ {HPS['window'] / 24:.0f} days", 44, "#ffffff", True, a=a, ha="center")
    a2 = fade(c, S(c, 0) + 1.2)
    if a2 > 0:
        box(c, 60, 440, 480, 90, "#ffffff", a=0.22 * a2, r=14)
        p = prog(c, S(c, 0) + 1.2, S(c, 0) + 2.0)
        ar = c.ax.annotate("", xy=(60 + 480 * p, 575), xytext=(60, 575),
                           arrowprops=dict(arrowstyle="<|-|>", color=TX, lw=4 * PT, mutation_scale=28))
        ar.set_alpha(a2)
        T(c, 60, 595, "max delay = ½ window", 40, TX, True, a=a2)
    T(c, 60, 655, f"≈ {HPS['max_lag'] / 24:.0f} days", 64, AQUA, True, a=fade(c, S(c, 1)))
    # sentence 2: rhythm
    a = fade(c, S(c, 2))
    ax = plot_axes(c, 60, 800, 960, 330, a, (0, 4.2), (-0.7, 1.35))
    if ax is not None:
        u = np.linspace(0, 4.2, 800)
        s = 0.5 + 0.5 * ease(prog(c, S(c, 2) + 0.8, E(c, 2) - 0.5))
        ax.plot(u, pulse(u - s), color=ORANGE, lw=9 * PT, alpha=0.6)
        ax.plot(u, pulse(u), color=BLUE, lw=3.5 * PT)
        at(ax, 0.01, 0.97, "rhythm (e.g. heartbeat)", 32, BLUE, bold=True)
        at(ax, 0.99, 0.97, f"copy shifted by {s:.1f} beat", 32, ORANGE, bold=True, ha="right")
    if c.t > E(c, 2) - 0.5:
        a = fade(c, E(c, 2) - 0.5)
        T(c, W / 2, 1160, "shift by 1 beat looks identical", 40, TX, True, a=a, ha="center", maxw=960)
        T(c, W / 2, 1220, "→ stop at ½ beat", 48, AQUA, True, a=a, ha="center")


# ---------------------------------------------------------------- scene 8: step 4 tolerance
def scene_step4(c):
    hdr(c, "STEP 4 · Tolerance", S(c, 0))
    a = fade(c, S(c, 1) - 0.3)
    ax = plot_axes(c, 60, 380, 960, 270, a, (0, 600), (-0.05, 2.55))
    if ax is not None:
        xs = np.arange(600)
        sl = 300 * ease(prog(c, S(c, 1) + 0.5, E(c, 1) - 0.3))
        ax.plot(xs, band(RAIN_Z, 1.45), color=BLUE, lw=3 * PT)
        ax.plot(xs, np.roll(band(RIVER_Z, 0.0), int(sl)), color=ORANGE, lw=3 * PT)
        at(ax, 0.01, 0.97, "rain", 32, BLUE, bold=True)
        at(ax, 0.01, 0.40, "river flow", 32, ORANGE, bold=True)
        if sl > 150:
            at(ax, 0.99, 0.40, "slid in time → no real link", 32, "#e66767", bold=True, ha="right", a=fade(c, S(c, 1) + 1.5))
    T(c, W / 2, 680, "any stability found now is fake", 44, "#e66767", True, a=fade(c, S(c, 2)), ha="center", maxw=960)
    # chart
    a = fade(c, S(c, 3) - 0.2)
    ax = plot_axes(c, 160, 790, 860, 480, a, (0, 72), (0, 8), xt=[0, 24, 48, 72], yt=[0, 2, 4, 6],
                   xl="tolerance (hours)", yl="fake stability (%)")
    if ax is not None:
        tol = np.arange(0, 73)
        val = np.array([NULL[int(t)] for t in tol])
        p = prog(c, S(c, 3) + 0.3, E(c, 3) - 0.2)
        k = max(2, int(len(tol) * p))
        ax.plot(tol[:k], val[:k], color=ORANGE, lw=5 * PT)
        if c.t > S(c, 3):
            al = fade(c, S(c, 3))
            ax.plot([0, 72], [5, 5], color=TX, lw=3 * PT, ls=(0, (6, 4)), alpha=al)
            ax.text(1, 5.25, "fake-stability limit 5%", color=TX, fontsize=34 * PT, alpha=al, va="bottom")
        if c.t > S(c, 4):
            al = fade(c, S(c, 4), 0.4)
            ax.plot([TOL], [NULL[TOL]], "o", color=AQUA, ms=26 * PT, alpha=al, zorder=5)
            ax.annotate(f"chosen: ±{round(TOL / 24)} days", xy=(TOL, NULL[TOL]), xytext=(50, 6.9),
                        color=AQUA, fontsize=44 * PT, fontweight="bold", alpha=al, ha="center", va="center",
                        arrowprops=dict(arrowstyle="-|>", color=AQUA, lw=4 * PT, alpha=al, mutation_scale=28, shrinkB=16))


# ---------------------------------------------------------------- scene 9: step 5 you
def scene_step5(c):
    y = hdr(c, "Your one choice:\ntime resolution", S(c, 0))
    cards = [("long window", "reliable, blurry timing", 560, S(c, 1)),
             ("short window", "sharp timing, noisier", 190, S(c, 1) + 0.5 * (E(c, 1) - S(c, 1)))]
    for i, (name, txt, wbar, t0) in enumerate(cards):
        a = fade(c, t0)
        yy = y + 40 + i * 400
        box(c, 60, yy, 960, 370, "#232322", GRID, a=a)
        T(c, 95, yy + 25, f"{name} →", 52, AQUA, True, a=a)
        T(c, 95, yy + 100, txt, 46, TX, a=a, maxw=890)
        ly = yy + 280
        if a > 0:
            c.ax.plot([95, 985], [ly, ly], color=GRID, lw=3 * PT, alpha=a)
            box(c, 540 - wbar / 2, ly - 28, wbar, 56, AQUA, a=0.55 * a, r=10)
            c.ax.add_patch(Circle((540, ly), 11, fc=ORANGE, ec="none", alpha=a))
        T(c, 95, ly + 38, "window", 28, TX2, a=a)
        T(c, 985, ly + 38, "when the link happens", 28, TX2, a=a, ha="right")


# ---------------------------------------------------------------- field scenes
ORDER = [
    ("hydrology_hourly", "Hydrology", "rain → river flow (hourly)", "better", "Fixed can't separate it from chance (1.2% vs 1.1%)"),
    ("energy_vic_summer", "Energy", "heat → electricity demand", "better", "Adaptive finds the 1-hour delay; fixed reports 0"),
    ("climate_nino34_gistemp", "Climate", "El Niño → global temperature", "mixed", "Stronger detection, but delay 1 mo vs expected ~3 mo"),
    ("beijing_pm25", "Air quality", "wind → PM2.5 pollution", "better", "Much stronger detection (57% vs chance 12%)"),
    ("gas_furnace", "Engineering", "gas feed → CO₂ (furnace)", "worse", "Only 296 points — too short for adaptive windows"),
    ("gasoline_weekly_diff", "Economics", "crude oil → petrol price", "tie", "Both detect it; neither resolves the 1–2 week lag"),
    ("jena_T_rh", "Weather", "temperature ↔ humidity (control)", "tie", "Both find 0 delay; adaptive flags its own calibration limit"),
    ("sleep_eeg_delta_sigma", "Neuroscience", "sleep EEG delta ↔ sigma", "better", "Same answer at any sampling rate (see scene 3)"),
    ("neurokit_hr_rsp", "Physiology", "heart rate ↔ breathing", "tie", "Adaptive picks near-original settings when they already fit"),
]
BOTTOM = 1400  # field-scene content stays above the subtitle box (y >= 1418)


def expected_text(d):
    e = d["expected"]
    if e["kind"] == "delay":
        if not e["value_samples"]:
            return "≈ 0"
        return "~" + readable(e["value_samples"] * d["sampling_interval_s"])
    return "unknown"


def scene_field(c, key):
    _, dom, rest, st, verdict = next(o for o in ORDER if o[0] == key)
    d, ex = FIELDS[key], EXC[key]
    cap = key == "sleep_eeg_delta_sigma"
    va = d["variants"]["adaptive_cap5min" if cap else "adaptive"]
    vf = d["variants"]["fixed_default"]
    info, pp, dt = va["info"], va["params_phys"], d["sampling_interval_s"]
    T(c, 60, 250, dom, 64, bold=True, a=fade(c, 0))
    T(c, 60, 322, rest, 40, TX, a=fade(c, 0), maxw=960)
    # legend
    al = fade(c, 0.1)
    lx = 60
    for nm, col, lab in [(d["leader"], BLUE, "leader"), (d["follower"], ORANGE, "follower")]:
        if al > 0:
            c.ax.plot([lx, lx + 44], [405, 405], color=col, lw=3 * PT, alpha=al)
        T(c, lx + 54, 389, nm, 30, TX2, a=al)
        lx += 54 + _w(nm, 30, False) + 40
    # chose panel: header + three parameter columns + explanation lines
    cols = [("Window", pp["window"]), ("Max delay", pp["max_lag"]), ("Tolerance", f"±{pp['tolerance']}")]
    bmax = max(info["bartlett"])
    win_txt = pp["window"]
    why = (f"Slowest signal ≈ 1 independent sample per {readable(bmax * dt)} → 30 of them = "
           f"{d['variants']['adaptive']['params_phys']['window'] if cap else win_txt}")
    if cap:
        why += f", capped to {win_txt}"
    elif info.get("window_clipped"):
        why += " (clipped)"
    extra = [(why, 28, TX2)]
    pers = [p for p in info["periods"] if p]
    if pers:
        extra.append((f"Rhythm detected: {readable(min(pers) * dt)} → delay capped", 28, TX2))
    if info.get("calibration_failed"):
        extra.append(("⚠ calibration warning: fake stability > 5% even at tolerance 0", 28, "#c98500"))
    ch_h = 16 + 48 + 84 + 8 + sum(len(wrap(s, px, False, 925)) * px * 1.25 for s, px, _ in extra) + 12
    # bottom-up layout
    vlines = wrap(verdict, 32, False, 910)
    bad_h = 16 + 50 + len(vlines) * 40 + 10
    bad_y = BOTTOM - bad_h
    del_y = bad_y - 10 - 84
    bar_y = del_y - 10 - 176
    ch_y = bar_y - 10 - ch_h
    plot_y, plot_h = 430, ch_y - 10 - 430
    # excerpt plot
    ax = sub_axes(c, 60, plot_y, 960, plot_h, fade(c, 0.1))
    s1, s2 = np.nan_to_num(np.array(ex["s1_z"], float)), np.nan_to_num(np.array(ex["s2_z"], float))
    idx = np.arange(len(s1))
    lines_plot(ax, [(idx, s1, BLUE), (idx, s2, ORANGE)], prog(c, 0.1, 1.6), fade(c, 0.1), lw=2)
    ax.set_xlim(0, len(s1))
    lo, hi = min(s1.min(), s2.min()), max(s1.max(), s2.max())
    ax.set_ylim(lo - 0.1 * (hi - lo), hi + 0.1 * (hi - lo))
    # chose panel
    a = fade(c, 0.6)
    box(c, 60, ch_y, 960, ch_h, "#232322", GRID, a=a)
    T(c, 90, ch_y + 16, "Adaptive chose:", 40, AQUA, True, a=a)
    if cap:
        T(c, 90 + _w("Adaptive chose:", 40, True) + 20, ch_y + 24, "with 5-min time resolution", 28, TX2, a=a)
    for i, (lab, val) in enumerate(cols):
        T(c, 90 + i * 300, ch_y + 16 + 48, lab, 26, TX2, a=a)
        T(c, 90 + i * 300, ch_y + 16 + 48 + 32, val, 36, TX, True, a=a)
    y = ch_y + 16 + 48 + 84 + 8
    for s, px, col in extra:
        y = T(c, 90, y, s, px, col, a=a, maxw=925, lh=1.25)
    # bars
    a = fade(c, 1.2)
    g = prog(c, 1.2, 2.2)
    X0, SC = 270, 5.4
    for i, (nm, v, col) in enumerate([("Fixed", vf, GRAY), ("Adaptive", va, AQUA)]):
        y = bar_y + 34 + i * 70
        T(c, 60, y + 6, nm, 34, TX2, a=a)
        box(c, X0, y, max(v["score"] * SC * g, 1), 54, col, a=a, r=14)
        if g > 0.99:
            T(c, X0 + v["score"] * SC + 14, y + 6, f"{v['score']:.1f}%", 36, TX, True, a=a)
        tx = X0 + v["null95"] * SC
        if a > 0:
            c.ax.plot([tx, tx], [y - 6, y + 60], color=TX, lw=1.5 * PT, alpha=a)
        if i == 0:
            T(c, max(tx, X0 + 45), bar_y + 4, "chance", 28, TX2, a=a, ha="center")
    # delay line
    a = fade(c, 1.8)
    T(c, 60, del_y, f"Delay found: fixed {minus(vf['median_stable_tau_phys'])} · adaptive {minus(va['median_stable_tau_phys'])}",
      36, a=a, maxw=960, lh=1.1)
    T(c, 60, del_y + 44, f"expected {expected_text(d)}", 36, TX2, a=a)
    # verdict badge
    a = fade(c, 2.3)
    icon, col = STATUS[st]
    box(c, 60, bad_y, 960, bad_h, "#232322", col, a=a, lw=3)
    T(c, 90, bad_y + 14, f"{icon} {st}", 44, col, True, a=a)
    T(c, 90, bad_y + 14 + 52, verdict, 32, TX, a=a, maxw=910, lh=1.25)


# ---------------------------------------------------------------- scoreboard / outro
def scene_board(c):
    T(c, 60, 260, "9 fields, honest scoreboard", 64, bold=True, a=fade(c, 0.1), maxw=960, lh=1.15)
    # sentences: 0 nine fields, 1 four better, 2 four ties or mixed, 3 one worse
    cur = None
    for i, grp in [(1, ("better",)), (2, ("tie", "mixed")), (3, ("worse",))]:
        if c.t >= S(c, i):
            cur = grp
    for i, (_, dom, _, st, _) in enumerate(ORDER):
        a = fade(c, S(c, 0) + 0.12 * i)
        if cur is not None and st not in cur:
            a *= 0.3
        y = 460 + i * 88
        icon, col = STATUS[st]
        T(c, 80, y, dom, 44, a=a)
        T(c, 1000, y, f"{icon} {st}", 44, col, True, a=a, ha="right")
        if a > 0:
            c.ax.plot([60, 1020], [y + 68, y + 68], color=GRID, lw=2 * PT, alpha=a)
    T(c, 60, 1290, "4 better · 4 tie/mixed · 1 worse", 48, AQUA, True, a=fade(c, S(c, 3) + 0.6), maxw=960)


def scene_outro(c):
    T(c, 60, 270, "Know the limits", 64, bold=True, a=fade(c, 0.1))
    pts = ["Very short data → keep fixed settings", "Long windows = coarser timing → set your time resolution",
           "Rhythmic raw waveforms → max delay = half a cycle"]
    y = 450
    for i, s in enumerate(pts):
        a = fade(c, S(c, i))
        c.ax.add_patch(Circle((78, y + 26), 11, fc=AQUA, ec="none", alpha=a))
        y = T(c, 110, y, s, 46, a=a, maxw=910) + 50
    T(c, 60, 1200, "pyTDS · adaptive branch (experimental)", 46, AQUA, True, a=fade(c, E(c, 2) - 0.8), maxw=960)


# ---------------------------------------------------------------- render
FUNCS = {"hook": scene_hook, "what": scene_what, "problem": scene_problem, "how_intro": scene_how_intro,
         "step1_memory": scene_step1, "step2_window": scene_step2, "step3_delay": scene_step3,
         "step4_tolerance": scene_step4, "step5_you": scene_step5, "board": scene_board, "outro": scene_outro}
IDS = [s["id"] for s in json.load(open(ROOT / "experiments/adaptive/reel_voiceover.json"))["scenes"]]
SCENES = [(f"{i + 1:02d}_{sid}", FUNCS.get(sid) or (lambda c, k=sid: scene_field(c, k)), sid) for i, sid in enumerate(IDS)]


def scene_durations():
    """Scene duration = narration + 0.3 s lead-in + 0.9 s tail, minimum 4 s."""
    nar = json.load(open(ROOT / "experiments/adaptive/reel_audio/durations.json"))
    return [max(4.0, nar[sid] + LEAD + TAIL) for _, _, sid in SCENES]


def settle_time(sid, dur):
    """After this many seconds nothing animates any more (last caption + longest animation)."""
    last = SENT[sid]["sentences"][-1]["chunks"][-1]["start"] + LEAD
    return min(dur, last + 3.2)


_FIG = None


def render(si, sec):
    global _FIG
    if _FIG is None:
        _FIG = plt.figure(figsize=(10.8, 19.2), dpi=100, facecolor=BG)
    fig = _FIG
    fig.clf()
    fig.patch.set_facecolor(BG)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(H, 0)
    ax.axis("off")
    c = Ctx(fig, ax, sec, SCENES[si][2])
    SCENES[si][1](c)
    subtitle(c)
    fig.canvas.draw()
    return np.ascontiguousarray(np.asarray(fig.canvas.buffer_rgba())[..., :3])


def _job(job):
    si, sec = job
    return render(si, sec).tobytes()


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--probe":  # --probe <scene index> <seconds> <out.png>
        plt.imsave(sys.argv[4], render(int(sys.argv[2]), float(sys.argv[3])))
        return
    FRAMES.mkdir(parents=True, exist_ok=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    for old in FRAMES.glob("*.png"):
        old.unlink()
    durs = scene_durations()
    json.dump(durs, open(ROOT / "experiments/adaptive/reel_audio/scene_durations.json", "w"))
    jobs, reps = [], []
    for si, ((_, _, sid), dur) in enumerate(zip(SCENES, durs)):
        n, st, frozen = int(round(dur * FPS)), settle_time(sid, dur), False
        for f in range(n):
            sec = f / FPS
            if sec >= st and frozen:
                reps[-1] += 1
            else:
                jobs.append((si, sec))
                reps.append(1)
                frozen = sec >= st
    cmd = ["ffmpeg", "-y", "-loglevel", "error", "-f", "rawvideo", "-pix_fmt", "rgb24", "-s", f"{W}x{H}", "-r", str(FPS),
           "-i", "-", "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "20", "-movflags", "+faststart", str(OUT)]
    pipe = subprocess.Popen(cmd, stdin=subprocess.PIPE)
    print(f"{len(jobs)} frames to render", flush=True)
    with Pool(4) as pool:
        for k, (buf, r) in enumerate(zip(pool.imap(_job, jobs, chunksize=4), reps)):
            for _ in range(r):
                pipe.stdin.write(buf)
            if k % 200 == 0:
                print(k, flush=True)
    pipe.stdin.close()
    pipe.wait()
    t0 = 0.0
    for si, ((name, _, sid), dur) in enumerate(zip(SCENES, durs)):
        plt.imsave(FRAMES / f"{name}.png", render(si, dur))
        print(f"{name}: {t0:.1f}s - {t0 + dur:.1f}s", flush=True)
        t0 += dur
    print("total", t0, "s ->", OUT)


if __name__ == "__main__":
    main()
