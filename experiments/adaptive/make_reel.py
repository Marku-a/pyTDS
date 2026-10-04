"""Render the adaptive-TDS Instagram reel (1080x1920, 30 fps, H.264, no audio).

Run from anywhere:  python experiments/adaptive/make_reel.py
Output: "adaptive-TDS reports/adaptive_TDS_reel.mp4" plus one still PNG per scene
in "adaptive-TDS reports/figures/reel_frames/".
"""
import json
import subprocess
import sys
from functools import lru_cache
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

RES = ROOT / "experiments/adaptive/results"
OUT = ROOT / "adaptive-TDS reports/adaptive_TDS_reel.mp4"
FRAMES = ROOT / "adaptive-TDS reports/figures/reel_frames"
FIELDS = json.load(open(RES / "fields.json"))
EXC = json.load(open(RES / "fields_excerpts.json"))

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
    def __init__(self, fig, ax, sec):
        self.fig, self.ax, self.t = fig, ax, sec  # t = seconds into scene


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


# ---------------------------------------------------------------- scene 1
def scene_hook(c):
    T(c, W / 2, 560, "Can an algorithm tune itself?", 84, bold=True, a=fade(c, 0.2), ha="center", maxw=960, lh=1.2)
    T(c, W / 2, 800, "Adaptive TDS · tested on 9 fields of science", 40, TX2, a=fade(c, 1.2), ha="center", maxw=960)
    ax = sub_axes(c, 60, 1100, 960, 260, 0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    x = np.linspace(0, 1, 400)
    f = lambda u: np.sin(2 * np.pi * 3 * u) + 0.5 * np.sin(2 * np.pi * 7.3 * u)
    p = prog(c, 1.6, 3.6)
    lines_plot(ax, [(x, f(x), BLUE), (x, f(x - 0.04), ORANGE)], p, fade(c, 1.6))
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
CAPS2 = ["TDS slides a window along two signals…", "…finds the delay where they match best…",
         "…and counts how often that delay stays stable = TDS score."]


def scene_what(c):
    s1, s2, r = DEMO
    n = len(s1)
    T(c, 60, 270, "What TDS does", 64, bold=True, a=fade(c, 0))
    X0, PW, PY, PH = 60, 960, 380, 400
    ax = sub_axes(c, X0, PY, PW, PH, fade(c, 0.2))
    idx = np.arange(n)
    lines_plot(ax, [(idx, s1 / s1.std() + 2.2, BLUE), (idx, s2 / s2.std() - 2.2, ORANGE)],
               prog(c, 0.3, 2.3), fade(c, 0.2))
    ax.set_xlim(0, n)
    ax.set_ylim(-4.8, 4.8)
    q = prog(c, 2.3, 8.7)
    centre = 30 + q * (n - 60)
    if c.t > 2.3:
        ax.add_patch(Rectangle((centre - 30, -4.7), 60, 9.4, fc="#ffffff", alpha=0.14, ec=TX, lw=1.5 * PT))
    k = int(np.argmin(abs(r["t_vec"] - centre)))
    if c.t > 2.3:
        T(c, 60, 840, f"delay in this window: {minus(str(int(r['tau'][k])))}", 44, bold=True)
    a = fade(c, 0.4)
    nd = len(r["tau"])
    for i in range(nd):
        if r["t_vec"][i] <= centre + 1e-6 and c.t > 2.3:
            c.ax.add_patch(Circle((60 + 30 + i * 80, 990), 24, fc=AQUA if r["stbl_lbl"][i] else GRAY, ec="none"))
    T(c, 60, 1050, "aqua = stable    gray = unstable", 30, TX2, a=a)
    for i, cap in enumerate(CAPS2):
        t0 = i * 3.0
        al = fade(c, t0) * (1 - fade(c, t0 + 3.0 - 0.3, 0.3) if i < 2 else 1)
        T(c, 60, 1160, cap, 52, a=al, maxw=960, bold=True)


# ---------------------------------------------------------------- scene 3
def scene_problem(c):
    T(c, 60, 270, "Fixed settings break", 64, bold=True, a=fade(c, 0))
    T(c, 60, 370, "Same night of sleep EEG — only the sampling rate changed", 40, TX2, a=fade(c, 0.3), maxw=960)
    vals, labs = [74.7, 58.1, 41.3, 26.8], ["0.5×", "1×", "2×", "4×"]
    base, top, bw = 1070, 560, 150
    for i, (v, lb) in enumerate(zip(vals, labs)):
        x = 90 + i * 235
        a = fade(c, 0.6 + 0.2 * i)
        hh = (base - top) * v / 80 * prog(c, 0.6 + 0.2 * i, 1.8 + 0.2 * i)
        box(c, x, base - hh, bw, hh, GRAY, a=a, r=24)
        if hh > 5:
            T(c, x + bw / 2, base - hh - 56, f"{v}%", 44, bold=True, a=a, ha="center")
        T(c, x + bw / 2, base + 20, lb, 44, TX2, a=a, ha="center")
    c.ax.plot([60, 1020], [base, base], color=GRID, lw=2 * PT, alpha=fade(c, 0.6))
    T(c, 60, 1190, "Original TDS uses fixed sample counts (window = 60 samples).", 40, a=fade(c, 2.6), maxw=960)
    T(c, 60, 1330, "Change the recording device → different answer.", 40, TX, bold=True, a=fade(c, 4.2), maxw=960)


# ---------------------------------------------------------------- scene 4
CARDS = [("WINDOW", "holds ~30 independent samples of the slowest signal (measured from how long each signal 'remembers' itself)"),
         ("MAX DELAY", "half the window, or half a rhythm cycle (longer delays are ambiguous)"),
         ("TOLERANCE", "widest setting that keeps fake stability ≤ 5% (tested on time-shifted copies with no coupling)")]


def scene_how(c):
    T(c, 60, 260, "How adaptive sets the parameters", 64, bold=True, a=fade(c, 0), maxw=960, lh=1.15)
    for i, (name, why) in enumerate(CARDS):
        a = fade(c, 0.4 + 3.0 * i)
        y = 440 + i * 285
        box(c, 60, y, 960, 270, "#232322", GRID, a=a)
        T(c, 95, y + 22, name, 52, AQUA, True, a=a)
        T(c, 95, y + 95, why, 40, TX, a=a, maxw=890, lh=1.22)
    T(c, 60, 1330, "You choose only one thing: your time resolution.", 52, AQUA, True, a=fade(c, 9.6), maxw=960)


# ---------------------------------------------------------------- scene 5
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
    T(c, 60, 372, "", 30)
    # legend
    al = fade(c, 0.1)
    for j, (nm, col) in enumerate([(d["leader"], BLUE), (d["follower"], ORANGE)]):
        pass
    lx = 60
    for nm, col, lab in [(d["leader"], BLUE, "leader"), (d["follower"], ORANGE, "follower")]:
        if al > 0:
            c.ax.plot([lx, lx + 44], [405, 405], color=col, lw=3 * PT, alpha=al)
        T(c, lx + 54, 389, nm, 30, TX2, a=al)
        lx += 54 + _w(nm, 30, False) + 40
    # chose panel text
    lines = [("Adaptive chose:", 40, TX, True), (f"Window {pp['window']}", 40, TX, False),
             (f"Max delay {pp['max_lag']}", 40, TX, False), (f"Tolerance ±{pp['tolerance']}", 40, TX, False)]
    bmax = max(info["bartlett"])
    win_txt = pp["window"]
    why = (f"Slowest signal ≈ 1 independent sample per {readable(bmax * dt)} → 30 of them = "
           f"{d['variants']['adaptive']['params_phys']['window'] if cap else win_txt}")
    if cap:
        why += f", capped to {win_txt}"
    elif info.get("window_clipped"):
        why += " (clipped)"
    extra = [(why, 30, TX2)]
    pers = [p for p in info["periods"] if p]
    if pers:
        extra.append((f"Rhythm detected: {readable(min(pers) * dt)} → delay capped", 30, TX2))
    if info.get("calibration_failed"):
        extra.append(("⚠ calibration warning: fake stability > 5% even at tolerance 0", 30, "#c98500"))
    ch_h = 20 + 50 + 3 * 52 + 8 + sum(len(wrap(s, px, False, 925)) * px * 1.28 for s, px, _ in extra) + 14
    # bottom-up layout
    vlines = wrap(verdict, 34, False, 880)
    bad_h = 24 + 56 + len(vlines) * 44 + 12
    bad_y = 1545 - bad_h
    del_y = bad_y - 12 - 84
    bar_y = del_y - 12 - 190
    ch_y = bar_y - 12 - ch_h
    plot_y, plot_h = 430, ch_y - 12 - 430
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
    y = ch_y + 18
    for i, (s, px, col, b) in enumerate(lines):
        T(c, 90, y, s, px, AQUA if i == 0 else col, b, a=a)
        if i == 0 and cap:
            T(c, 90 + _w(s, 40, True) + 20, y + 8, "with 5-min time resolution", 30, TX2, a=a)
        y += 50 if i == 0 else 52
    y += 8
    for s, px, col in extra:
        y = T(c, 90, y, s, px, col, a=a, maxw=925)
    # bars
    a = fade(c, 1.2)
    g = prog(c, 1.2, 2.2)
    X0, SC = 270, 5.4
    for i, (nm, v, col) in enumerate([("Fixed", vf, GRAY), ("Adaptive", va, AQUA)]):
        y = bar_y + 36 + i * 76
        T(c, 60, y + 6, nm, 34, TX2, a=a)
        box(c, X0, y, max(v["score"] * SC * g, 1), 56, col, a=a, r=14)
        if g > 0.99:
            T(c, X0 + v["score"] * SC + 14, y + 6, f"{v['score']:.1f}%", 36, TX, True, a=a)
        tx = X0 + v["null95"] * SC
        if a > 0:
            c.ax.plot([tx, tx], [y - 6, y + 62], color=TX, lw=1.5 * PT, alpha=a)
        if i == 0:
            T(c, max(tx, X0 + 45), bar_y - 2, "chance", 30, TX2, a=a, ha="center")
    # delay line
    a = fade(c, 1.8)
    T(c, 60, del_y, f"Delay found: fixed {minus(vf['median_stable_tau_phys'])} · adaptive {minus(va['median_stable_tau_phys'])}",
      36, a=a, maxw=960, lh=1.1)
    T(c, 60, del_y + 44, f"expected {expected_text(d)}", 36, TX2, a=a)
    # verdict badge
    a = fade(c, 2.3)
    icon, col = STATUS[st]
    box(c, 60, bad_y, 960, bad_h, "#232322", col, a=a, lw=3)
    T(c, 90, bad_y + 18, f"{icon} {st}", 48, col, True, a=a)
    T(c, 90, bad_y + 24 + 56, verdict, 34, TX, a=a, maxw=880, lh=1.29)


# ---------------------------------------------------------------- scene 6
def scene_board(c):
    T(c, 60, 260, "9 fields, honest scoreboard", 64, bold=True, a=fade(c, 0), maxw=960, lh=1.15)
    for i, (_, dom, _, st, _) in enumerate(ORDER):
        a = fade(c, 0.9 + 0.35 * i)
        y = 460 + i * 88
        icon, col = STATUS[st]
        T(c, 80, y, dom, 44, a=a)
        T(c, 1000, y, f"{icon} {st}", 44, col, True, a=a, ha="right")
        if a > 0:
            c.ax.plot([60, 1020], [y + 68, y + 68], color=GRID, lw=2 * PT, alpha=a)
    T(c, 60, 1290, "4 better · 4 tie/mixed · 1 worse", 48, AQUA, True, a=fade(c, 4.4), maxw=960)


# ---------------------------------------------------------------- scene 7
def scene_outro(c):
    T(c, 60, 270, "Know the limits", 64, bold=True, a=fade(c, 0))
    pts = ["Very short data → keep fixed settings", "Long windows = coarser timing → set your time resolution",
           "Rhythmic raw waveforms → max delay = half a cycle"]
    y = 450
    for i, s in enumerate(pts):
        a = fade(c, 0.4 + 1.0 * i)
        c.ax.add_patch(Circle((78, y + 26), 11, fc=AQUA, ec="none", alpha=a))
        y = T(c, 110, y, s, 46, a=a, maxw=910) + 50
    T(c, 60, 1200, "pyTDS · adaptive branch (experimental)", 46, AQUA, True, a=fade(c, 4.2), maxw=960)


# ---------------------------------------------------------------- render
SCENES = [("01_hook", scene_hook, 4, 3.6), ("02_what", scene_what, 9, 9), ("03_problem", scene_problem, 7, 4.7),
          ("04_how", scene_how, 12, 10.1)]
SCENES += [(f"05_{i + 1}_{o[0]}", (lambda c, k=o[0]: scene_field(c, k)), FIELD_DUR, 2.8) for i, o in enumerate(ORDER)]
SCENES += [("06_board", scene_board, 7, 4.9), ("07_outro", scene_outro, 7, 4.7)]


def main():
    FRAMES.mkdir(parents=True, exist_ok=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(10.8, 19.2), dpi=100, facecolor=BG)
    cmd = ["ffmpeg", "-y", "-loglevel", "error", "-f", "rawvideo", "-pix_fmt", "rgb24", "-s", f"{W}x{H}", "-r", str(FPS),
           "-i", "-", "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "20", "-movflags", "+faststart", str(OUT)]
    pipe = subprocess.Popen(cmd, stdin=subprocess.PIPE)

    def render(fn, sec, dur):
        fig.clf()
        fig.patch.set_facecolor(BG)
        ax = fig.add_axes([0, 0, 1, 1])
        ax.set_xlim(0, W)
        ax.set_ylim(H, 0)
        ax.axis("off")
        fn(Ctx(fig, ax, sec))
        fig.canvas.draw()
        return np.ascontiguousarray(np.asarray(fig.canvas.buffer_rgba())[..., :3])

    t0 = 0.0
    for name, fn, dur, settle in SCENES:
        n, last = int(round(dur * FPS)), None
        for f in range(n):
            sec = f / FPS
            if sec >= settle and last is not None and settle < dur:
                frame = last
            else:
                frame = render(fn, sec, dur)
                if sec >= settle:
                    last = frame
            pipe.stdin.write(frame.tobytes())
        final = render(fn, dur, dur)
        plt.imsave(FRAMES / f"{name}.png", final)
        print(f"{name}: {t0:.1f}s - {t0 + dur:.1f}s", flush=True)
        t0 += dur
    pipe.stdin.close()
    pipe.wait()
    print("total", t0, "s ->", OUT)


if __name__ == "__main__":
    main()
