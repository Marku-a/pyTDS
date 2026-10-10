"""3Blue1Brown-style explainer of adaptive TDS (Manim Community 0.22).

Env: VDATA=video_data.json (from manim_data.py), VOICE=<piper .onnx>, VCACHE=<dir for voice wavs>.
Render: manim -qm atds_3b1b.py S01Hook ... S11Recap ; then concatenate (see render.sh).
Narration voice: Piper TTS (en-us-lessac-medium). Every number shown comes from VDATA.
"""
import hashlib
import json
import os
import subprocess
import wave
from contextlib import contextmanager

import numpy as np
from manim import *  # noqa: F403

D = json.load(open(os.environ["VDATA"]))
VOICE = os.environ.get("VOICE", "")
CACHE = os.environ.get("VCACHE", "vo_cache")
os.makedirs(CACHE, exist_ok=True)

BG, INK, INK2, GRID = "#141517", "#eceae4", "#a9a8a0", "#3a3b40"
C_TDS, C_ATDS, C_AQUA, C_NULL, C_GOLD = "#3987e5", "#d95926", "#199e70", "#8a8c92", "#e3b341"
config.background_color = BG
Text.set_default(font="DejaVu Sans", color=INK)
MathTex.set_default(color=INK)


SENT_GAP, PARA_GAP, LENGTH_SCALE = 0.3, 1.0, float(os.environ.get("LENGTH_SCALE", "1.0"))


TTS = os.environ.get("TTS", "piper")          # "kokoro" (natural, default for the final video) or "piper"
KOKORO_DIR = os.environ.get("KOKORO_DIR", "")
KOKORO_VOICE = os.environ.get("KOKORO_VOICE", "af_heart")
_KOKORO = None


def _synth_sentence(text):
    tag = f"kokoro:{KOKORO_VOICE}" if TTS == "kokoro" else VOICE
    h = hashlib.md5((tag + str(LENGTH_SCALE) + text).encode()).hexdigest()[:14]
    p = os.path.join(CACHE, f"s_{h}.wav")
    if not os.path.exists(p):
        if TTS == "kokoro":
            global _KOKORO
            import soundfile as sf
            from kokoro_onnx import Kokoro
            if _KOKORO is None:
                _KOKORO = Kokoro(os.path.join(KOKORO_DIR, "kokoro-v1.0.onnx"), os.path.join(KOKORO_DIR, "voices-v1.0.bin"))
            audio, sr = _KOKORO.create(text, voice=KOKORO_VOICE, speed=1.0 / LENGTH_SCALE, lang="en-us")
            sf.write(p, audio, sr, subtype="PCM_16")
        else:
            subprocess.run(["/usr/bin/python3", "-m", "piper", "-m", VOICE, "-f", p, "--length-scale", str(LENGTH_SCALE)],
                           input=text.encode(), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    with wave.open(p) as w:
        return w.getparams(), w.readframes(w.getnframes())


def tts(text):
    """Synthesize sentence by sentence with SENT_GAP between them. Returns (wav path, duration, [(start, end, sentence)])."""
    import re
    sents = [x.strip() for x in re.split(r"(?<=[.!?])\s+", text) if x.strip()]
    h = hashlib.md5((TTS + KOKORO_VOICE + VOICE + str(LENGTH_SCALE) + str(SENT_GAP) + text).encode()).hexdigest()[:14]
    out = os.path.join(CACHE, f"b_{h}.wav")
    frames, cues, t, params = [], [], 0.0, None
    for i, snt in enumerate(sents):
        params, fr = _synth_sentence(snt)
        rate, width = params.framerate, params.sampwidth
        dur = len(fr) / (width * params.nchannels) / rate
        cues.append((t, t + dur, snt))
        frames.append(fr)
        t += dur
        if i < len(sents) - 1:
            gap = b"\x00" * int(SENT_GAP * rate) * width * params.nchannels
            frames.append(gap)
            t += SENT_GAP
    if not os.path.exists(out):
        with wave.open(out, "wb") as w:
            w.setnchannels(params.nchannels)
            w.setsampwidth(params.sampwidth)
            w.setframerate(params.framerate)
            w.writeframes(b"".join(frames))
    return out, t, cues


def axes(xr, yr, w, h, xlabel=None, ylabel=None, nums=True):
    ax = Axes(x_range=xr, y_range=yr, x_length=w, y_length=h, tips=False,
              axis_config={"color": GRID, "stroke_width": 2, "include_numbers": nums, "font_size": 20,
                           "decimal_number_config": {"color": INK2, "num_decimal_places": 0}})
    g = VGroup(ax)
    if xlabel:
        g.add(Text(xlabel, font_size=20, color=INK2).next_to(ax, DOWN, buff=0.3))
    if ylabel:
        g.add(Text(ylabel, font_size=20, color=INK2).rotate(PI / 2).next_to(ax, LEFT, buff=0.45))
    return ax, g


def curve(ax, ys, color, x0=0, width=2.5):
    xs = np.arange(len(ys)) + x0
    return ax.plot_line_graph(xs, np.asarray(ys), line_color=color, add_vertex_dots=False, stroke_width=width)


def chapter(n, title):
    t = Text(f"{n}  ·  {title}", font_size=24, color=INK2).to_corner(UL, buff=0.4)
    return t


class VO(Scene):
    def setup(self):
        self.cues = []

    @contextmanager
    def beat(self, text, pause=PARA_GAP):
        p, d, cues = tts(text)
        t0 = self.renderer.time
        self.add_sound(p)
        self.cues += [(t0 + a, t0 + b, c) for a, b, c in cues]
        yield d
        rest = t0 + d + pause - self.renderer.time
        if rest > 0.03:
            self.wait(rest)

    def intro(self, n, title):
        """Chapter card: big title, a silent pause, then it moves to the corner."""
        big = VGroup(Text(f"Chapter {n}", font_size=30, color=INK2), Text(title[0].upper() + title[1:], font_size=52)).arrange(DOWN, buff=0.3)
        self.play(FadeIn(big, shift=UP * 0.3), run_time=0.8)
        self.wait(1.2)
        small = chapter(n, title)
        self.play(ReplacementTransform(big, small), run_time=0.7)
        self.wait(0.3)

    def tear_down(self):
        os.makedirs("cues", exist_ok=True)
        json.dump({"duration": self.renderer.time, "cues": self.cues}, open(os.path.join("cues", type(self).__name__ + ".json"), "w"))
        super().tear_down()


# ------------------------------------------------------------------ 1 hook
class S01Hook(VO):
    def construct(self):
        t = np.linspace(0, 40, 400)
        breath = np.sin(2 * np.pi * t / 8)
        heart = 0.8 * np.sin(2 * np.pi * (t - 2.5) / 8) + 0.12 * np.sin(2 * np.pi * t / 0.9)
        ax, _ = axes([0, 40, 10], [-1.5, 1.5, 1], 11, 2.2, nums=False)
        ax2, _ = axes([0, 40, 10], [-1.5, 1.5, 1], 11, 2.2, nums=False)
        ax.shift(UP * 1.3)
        ax2.shift(DOWN * 1.5)
        cb = ax.plot_line_graph(t, breath, line_color=C_AQUA, add_vertex_dots=False, stroke_width=4)
        ch = ax2.plot_line_graph(t, heart, line_color=C_GOLD, add_vertex_dots=False, stroke_width=4)
        lb = Text("breathing", font_size=28, color=C_AQUA).next_to(ax, UP, buff=0.1).align_to(ax, LEFT)
        lh = Text("heart rate", font_size=28, color=C_GOLD).next_to(ax2, UP, buff=0.1).align_to(ax2, RIGHT)
        with self.beat("Here's something you've felt without noticing. Your heart and your breathing aren't independent. Breathe in, and your heart speeds up a little. But not instantly. It happens a moment later, with a delay."):
            self.play(Create(cb), FadeIn(lb), run_time=2.5)
            self.play(Create(ch), FadeIn(lh), run_time=2.5)
            p1 = ax.c2p(2, 1)
            p2 = ax2.c2p(4.5, 0.8)
            l1 = DashedLine(p1, [p1[0], ax2.c2p(0, -1.5)[1], 0], color=INK2)
            l2 = DashedLine([p2[0], ax.c2p(0, 1.5)[1], 0], p2, color=INK2)
            arr = DoubleArrow([p1[0], 0, 0], [p2[0], 0, 0], color=INK, buff=0, stroke_width=4)
            lab = MathTex(r"\text{delay}", font_size=36).next_to(arr, UP, buff=0.1)
            self.play(Create(l1), Create(l2), GrowFromCenter(arr), FadeIn(lab), run_time=1.5)
        self.play(FadeOut(VGroup(ax, ax2, cb, ch, lb, lh, l1, l2, arr, lab)), run_time=0.8)
        names = ["brain", "heart", "lungs", "eyes", "muscles"]
        cols = [C_TDS, C_GOLD, C_AQUA, INK2, C_ATDS]
        nodes = VGroup(*[VGroup(Circle(0.55, color=c, fill_opacity=0.15, stroke_width=4), Text(n, font_size=22))
                         for n, c in zip(names, cols)])
        for i, nd in enumerate(nodes):
            nd.move_to(2.1 * np.array([np.cos(PI / 2 + i * TAU / 5), np.sin(PI / 2 + i * TAU / 5), 0]) + DOWN * 0.9)
        links = VGroup(*[Line(nodes[i].get_center(), nodes[j].get_center(), color=INK2, stroke_width=3, buff=0.6)
                         for i, j in [(0, 1), (1, 2), (0, 3), (0, 4), (2, 4)]])
        title = Text("Time Delay Stability (TDS)", font_size=40).to_edge(UP)
        with self.beat("Physiologists want to know which systems in the body are linked like this, and when. One popular tool for that is called Time Delay Stability. T D S, for short."):
            self.play(LaggedStart(*[GrowFromCenter(n) for n in nodes], lag_ratio=0.15), run_time=2)
            self.play(LaggedStart(*[Create(l) for l in links], lag_ratio=0.2), run_time=2)
            self.play(Write(title), run_time=1.2)
        sub = Text("adaptive TDS", font_size=48, color=C_ATDS).next_to(title, DOWN, buff=0.25)
        with self.beat("In this video, we'll look at a version that picks its own settings. We'll call it adaptive T D S. But to get there, we first need to see how T D S itself thinks."):
            self.play(Write(sub), run_time=1.5)
            self.play(Indicate(sub, color=C_ATDS), run_time=1.2)
        self.play(FadeOut(Group(*self.mobjects)), run_time=0.8)


# ------------------------------------------------------------------ 2 cross-correlation
class S02XCorr(VO):
    def construct(self):
        F = D["fast"]
        s = F["win_start"]
        x = np.array(F["x"][s:s + 60])
        y = np.array(F["y"][s:s + 60])
        x = (x - x.mean()) / x.std()
        y = (y - y.mean()) / y.std()
        self.intro(1, "finding the delay")
        ax, g = axes([0, 60, 10], [-3, 3, 3], 10.5, 2.0, xlabel="time in the window (s)", nums=True)
        g.shift(UP * 1.6)
        cx = curve(ax, x, C_AQUA, width=3)
        lx = Text("signal 1", font_size=22, color=C_AQUA).next_to(ax, UL, buff=0.05).shift(DOWN * 0.3)
        k = ValueTracker(0)

        def shifted():
            sh = int(round(k.get_value()))
            ys = np.roll(y, -sh)
            return curve(ax, ys, C_GOLD, width=3)
        cy = always_redraw(shifted)
        ly = Text("signal 2, shifted", font_size=22, color=C_GOLD).next_to(lx, DOWN, buff=0.1).align_to(lx, LEFT)
        with self.beat("Okay. Take a one minute piece of two signals. Is the second one just a delayed copy of the first? Well, let's find out. We slide it in time."):
            self.play(FadeIn(g), Create(cx), FadeIn(lx), run_time=2)
            self.play(Create(cy), FadeIn(ly), run_time=1.5)
            self.play(k.animate.set_value(-12), run_time=1.5)
            self.play(k.animate.set_value(12), run_time=2)
            self.play(k.animate.set_value(-30), run_time=1.5)
        lags = np.array(F["xc_lag"])
        C = np.array(F["xc_C"])
        order = np.argsort(lags)
        lags, C = lags[order], C[order]
        bx, bg = axes([-30, 30, 10], [-1, 1, 0.5], 10.5, 2.3, xlabel="shift τ (s)", ylabel="C(τ)")
        bg.shift(DOWN * 1.9)

        def traced():
            v = k.get_value()
            m = lags <= v + 1e-9
            if m.sum() < 2:
                return VGroup()
            return bx.plot_line_graph(lags[m], C[m], line_color=C_TDS, add_vertex_dots=False, stroke_width=4)
        tr = always_redraw(traced)
        dot = always_redraw(lambda: Dot(bx.c2p(k.get_value(), np.interp(k.get_value(), lags, C)), color=C_TDS, radius=0.08))
        with self.beat("At every shift, we multiply the two curves point by point, and add it all up. When the shapes line up, most of those products are positive, so the sum gets big. This running score is called the cross correlation."):
            self.play(FadeIn(bg), run_time=1)
            self.add(tr, dot)
            self.play(k.animate.set_value(30), run_time=9, rate_func=linear)
            full = bx.plot_line_graph(lags, C, line_color=C_TDS, add_vertex_dots=False, stroke_width=4)
            self.remove(tr)
            self.add(full)
        peak = int(lags[np.argmax(np.abs(C))])
        with self.beat(f"And its peak tells us the delay. We call it tau zero. Here, the second signal follows the first by {peak} seconds."):
            self.play(k.animate.set_value(peak), run_time=2)
            pk = Dot(bx.c2p(peak, C[lags == peak][0]), color=C_ATDS, radius=0.13)
            lab = MathTex(rf"\tau_0 = {peak}\,\text{{s}}", font_size=44, color=C_ATDS).next_to(pk, UP, buff=0.2)
            self.play(GrowFromCenter(pk), Write(lab), run_time=1.2)
            self.play(Indicate(cy, color=C_ATDS), run_time=1.2)
        self.play(FadeOut(Group(*self.mobjects)), run_time=0.8)


# ------------------------------------------------------------------ 3 stability
class S03Stability(VO):
    def construct(self):
        F = D["fast"]
        self.intro(2, "stable delays")
        ax, g = axes([0, 900, 100], [-4, 4, 4], 11, 1.7, nums=False)
        g.shift(UP * 2.0)
        x, y = np.array(F["x"]), np.array(F["y"])
        c1 = curve(ax, x / 1.3 + 1.6, C_AQUA, width=1.6)
        c2 = curve(ax, y / 1.3 - 1.6, C_GOLD, width=1.6)
        T = F["tds"]
        cen, lag, st = np.array(T["centre"]), np.array(T["lag"]), np.array(T["stable"])
        bx, bg = axes([0, 900, 100], [-30, 30, 10], 11, 3.0, xlabel="time (s)", ylabel="τ₀ (s)")
        bg.shift(DOWN * 1.5)
        win = Rectangle(width=ax.c2p(60, 0)[0] - ax.c2p(0, 0)[0], height=1.9, color=C_TDS, fill_opacity=0.18, stroke_width=2)
        win.move_to(ax.c2p(30, 0))
        dots = [Dot(bx.c2p(c, l), radius=0.06, color=C_TDS) for c, l in zip(cen, lag)]
        with self.beat("Now, let's do that again. And again. Window after window, along the whole recording. Every window gives us one delay. One dot."):
            self.play(FadeIn(g), Create(c1), Create(c2), FadeIn(bg), run_time=2)
            self.add(win)
            for i in range(len(dots)):
                self.play(win.animate.move_to(ax.c2p(cen[i], 0)), FadeIn(dots[i], scale=2), run_time=0.17)
        true = DashedLine(bx.c2p(0, 6), bx.c2p(900, 6), color=INK2)
        tl = Text("same delay, window after window", font_size=22, color=INK2).next_to(bx.c2p(900, 6), UP, buff=0.12).shift(LEFT * 2.2)
        rng = np.random.default_rng(3)
        rnd = VGroup(*[Dot(bx.c2p(c, rng.integers(-30, 31)), radius=0.06, color=C_NULL) for c in cen])
        with self.beat("If two systems really are coupled, the dots line up at the same delay. If they're not, they just jump around."):
            self.play(Create(true), FadeIn(tl), run_time=1.5)
            self.wait(1)
            self.play(FadeIn(rnd), *[d.animate.set_opacity(0.15) for d in dots], run_time=1.5)
            self.wait(1.2)
            self.play(FadeOut(rnd), *[d.animate.set_opacity(1) for d in dots], FadeOut(tl), run_time=1)
        i0 = 5
        grp = VGroup(*dots[i0:i0 + 5])
        box = SurroundingRectangle(grp, color=C_ATDS, buff=0.15)
        band = Rectangle(width=bx.c2p(cen[i0 + 4], 0)[0] - bx.c2p(cen[i0], 0)[0] + 0.4,
                         height=bx.c2p(0, 7)[1] - bx.c2p(0, 5)[1], color=C_ATDS, fill_opacity=0.25, stroke_width=0)
        band.move_to(bx.c2p((cen[i0] + cen[i0 + 4]) / 2, 6))
        rule = Text("stable: ≥ 4 of 5 neighbouring delays within ± tolerance", font_size=24, color=C_ATDS).next_to(bg, DOWN, buff=0.1)
        with self.beat("T D S turns this into a simple rule. A window counts as stable when at least four out of five neighbouring delays agree, within a small tolerance. In the published method, that's plus or minus one second."):
            self.play(Create(box), run_time=1)
            self.play(FadeIn(band), Write(rule), run_time=2)
            self.wait(1)
        score = T["score"]
        ctr = ValueTracker(0)
        sc = always_redraw(lambda: Text(f"TDS score = {ctr.get_value():.0f} %", font_size=34, color=C_TDS).to_corner(UR, buff=0.4))
        with self.beat("Let's fill in the stable windows. The T D S score is just the share of windows that are stable."):
            self.add(sc)
            self.play(*[d.animate.set_fill(C_TDS, 1).scale(1.4) for d, s_ in zip(dots, st) if s_ == 1],
                      *[d.animate.set_fill(BG, 1) for d, s_ in zip(dots, st) if s_ == 0],
                      ctr.animate.set_value(score), FadeOut(box), FadeOut(band), run_time=2.5)
        self.play(FadeOut(Group(*self.mobjects)), run_time=0.8)


# ------------------------------------------------------------------ 4 hidden assumption
class S04Assumption(VO):
    def construct(self):
        self.intro(3, "the hidden assumption")
        card = VGroup(*[VGroup(Text(a, font_size=26, color=INK2), Text(b, font_size=40, color=C_TDS)).arrange(DOWN, buff=0.15)
                        for a, b in [("window", "60 s"), ("step", "30 s"), ("delay search", "±30 s"), ("tolerance", "±1 s")]]).arrange(RIGHT, buff=0.9)
        with self.beat("Now, those settings, a sixty second window, a thirty second step, plus or minus one second, were chosen for sleep recordings sampled once per second. And they quietly assume something. They assume that one minute holds plenty of independent information."):
            self.play(LaggedStart(*[FadeIn(c, shift=UP) for c in card], lag_ratio=0.25), run_time=2.5)
            q = Text("enough independent information in one minute?", font_size=30, color=C_GOLD).next_to(card, DOWN, buff=0.7)
            self.play(Write(q), run_time=2)
        self.play(FadeOut(card), FadeOut(q), run_time=0.6)
        M = D["memory"]
        rows = []
        for i, (key, col, name) in enumerate([("fast", C_AQUA, "fast signal"), ("slow", C_GOLD, "slow signal")]):
            ax, _ = axes([0, 60, 10], [-3, 3, 3], 8.5, 1.8, nums=False)
            ax.shift(UP * (1.3 - 2.9 * i) + LEFT * 1.5)
            xs = np.array(M[key]["x"][:60])
            c = curve(ax, xs, col, width=3)
            B = M[key]["B"]
            beads = VGroup(*[Dot(ax.c2p(t, xs[int(t)]), radius=0.1, color=INK) for t in np.arange(0, 60, B)])
            n = Text(f"≈ {60 / B:.0f} independent\nsamples / minute", font_size=24, color=col).next_to(ax, RIGHT, buff=0.4)
            lab = Text(name, font_size=24, color=col).next_to(ax, UP, buff=0.05).align_to(ax, LEFT)
            rows.append((ax, c, beads, n, lab))
        with self.beat("Look at a fast signal. Every few seconds, it's doing something new. Now look at a slow one. Neighbouring seconds have almost the same value. So, sixty points, but really only a couple of facts."):
            for ax, c, beads, n, lab in rows:
                self.play(Create(ax), Create(c), FadeIn(lab), run_time=1.4)
                self.play(LaggedStart(*[GrowFromCenter(b) for b in beads], lag_ratio=0.05), FadeIn(n), run_time=1.6)
        self.play(*[FadeOut(VGroup(*r)) for r in rows], run_time=0.6)
        S = D["slow"]
        T = S["tds"]
        cen, lag = np.array(T["centre"]), np.array(T["lag"])
        m = cen < 3600
        bx, bg = axes([0, 3600, 600], [-30, 30, 10], 11, 4.2, xlabel="time (s)", ylabel="τ₀ (s)")
        bg.shift(DOWN * 0.4)
        true = DashedLine(bx.c2p(0, 10), bx.c2p(3600, 10), color=C_GOLD, stroke_width=3)
        tl = Text("true delay: 10 s", font_size=22, color=C_GOLD).next_to(bx.c2p(3600, 10), UP, buff=0.1).shift(LEFT * 1.2)
        dots = VGroup(*[Dot(bx.c2p(c, l), radius=0.05, color=C_TDS) for c, l in zip(cen[m], lag[m])])
        sc = Text(f"TDS score {T['score']:.1f} %   ·   chance {S['null_mean']['tds']:.1f} %", font_size=28, color=C_TDS).to_corner(UR, buff=0.45)
        with self.beat("And with that little information, the delay in each window is close to random. Here are two slow signals that really are coupled, with a ten second delay. Watch what classic T D S does. The dots scatter. Its score is under three percent, barely above chance."):
            self.play(FadeIn(bg), Create(true), FadeIn(tl), run_time=1.5)
            self.play(LaggedStart(*[FadeIn(d, scale=2) for d in dots], lag_ratio=0.01), run_time=5)
            self.play(Write(sc), run_time=1.5)
        self.play(FadeOut(Group(*self.mobjects)), run_time=0.8)


# ------------------------------------------------------------------ 5 memory
class S05Memory(VO):
    def construct(self):
        self.intro(4, "step 1: measure memory")
        M = D["memory"]
        xs = np.array(M["slow"]["x"][:240])
        ax, _ = axes([0, 240, 60], [-3, 3, 3], 10.5, 1.6, nums=False)
        ax.shift(UP * 2.1)
        k = ValueTracker(0)
        c1 = curve(ax, xs, C_GOLD, width=2.5)
        c2 = always_redraw(lambda: curve(ax, xs[:240 - int(k.get_value())], C_AQUA, x0=int(k.get_value()), width=2.5).set_stroke(opacity=0.85))
        kl = always_redraw(lambda: MathTex(rf"k = {int(k.get_value())}\,\text{{s}}", font_size=34).next_to(ax, RIGHT, buff=0.2))
        with self.beat("So here's the idea behind adaptive T D S. Before anything else, measure memory. Put a signal next to a copy of itself, shifted by k seconds, and ask: how similar are they?"):
            self.play(Create(ax), Create(c1), run_time=1.5)
            self.add(c2, kl)
            self.play(k.animate.set_value(15), run_time=2)
            self.play(k.animate.set_value(0), run_time=1)
        a_s = np.array(M["slow"]["acf"])
        a_f = np.array(M["fast"]["acf"])
        bx, bg = axes([0, 120, 20], [-0.2, 1, 0.5], 7.5, 3.4, xlabel="shift k (s)", ylabel="similarity r(k)")
        bg.shift(DOWN * 1.2 + LEFT * 1.9)

        def tr():
            n = int(k.get_value()) + 1
            if n < 2:
                return VGroup()
            return bx.plot_line_graph(np.arange(n), a_s[:n], line_color=C_GOLD, add_vertex_dots=False, stroke_width=4)
        acf_s = always_redraw(tr)
        with self.beat("With no shift, they match perfectly. A correlation of one. As the shift grows, the match fades. That curve is the autocorrelation. And notice, a fast signal forgets within seconds."):
            self.play(FadeIn(bg), run_time=0.8)
            self.add(acf_s)
            self.play(k.animate.set_value(120), run_time=6, rate_func=linear)
            cf = bx.plot_line_graph(np.arange(121), a_f[:121], line_color=C_AQUA, add_vertex_dots=False, stroke_width=4)
            self.play(Create(cf), run_time=1.5)
            lf = Text("fast", font_size=24, color=C_AQUA).move_to(bx.c2p(14, 0.2))
            ls = Text("slow", font_size=24, color=C_GOLD).move_to(bx.c2p(60, 0.5))
            self.play(FadeIn(lf), FadeIn(ls), run_time=0.6)
        zc = int(np.argmax(a_s <= 0)) if (a_s <= 0).any() else 120
        bars = VGroup(*[Rectangle(width=bx.c2p(1, 0)[0] - bx.c2p(0, 0)[0], height=max(1e-3, bx.c2p(0, a_s[i] ** 2)[1] - bx.c2p(0, 0)[1]),
                                  color=C_GOLD, fill_opacity=0.45, stroke_width=0).move_to(bx.c2p(i + 0.5, a_s[i] ** 2 / 2))
                        for i in range(1, min(zc, 120))])
        f = MathTex(r"B = 1 + 2\sum_{k \ge 1} r(k)^2", font_size=46).to_edge(RIGHT, buff=0.4).shift(UP * 0.2)
        f2 = Text("summed until r first crosses zero", font_size=20, color=INK2).next_to(f, DOWN, buff=0.2)
        with self.beat("The Bartlett factor adds up the squared autocorrelation, until the curve first crosses zero. B equals one, plus twice the sum of r squared. Think of it as: how many samples in a row are worth one independent sample?"):
            self.play(Write(f), FadeIn(f2), run_time=2)
            self.play(LaggedStart(*[GrowFromEdge(b, DOWN) for b in bars], lag_ratio=0.01), run_time=3)
        bf = Text(f"fast:  B ≈ {M['fast']['B']:.0f}", font_size=30, color=C_AQUA).next_to(f2, DOWN, buff=0.5)
        bs = Text(f"slow:  B ≈ {M['slow']['B']:.0f}", font_size=30, color=C_GOLD).next_to(bf, DOWN, buff=0.2)
        with self.beat(f"For the fast signal, B is about {M['fast']['B']:.0f}. For the slow one, it's about {M['slow']['B']:.0f}. So {M['slow']['B']:.0f} seconds of that slow signal carry roughly one second's worth of new information."):
            self.play(FadeIn(bf, shift=UP), run_time=1)
            self.play(FadeIn(bs, shift=UP), run_time=1)
        self.play(FadeOut(Group(*self.mobjects)), run_time=0.8)


# ------------------------------------------------------------------ 6 window
class S06Window(VO):
    def construct(self):
        self.intro(5, "step 2: size the window")
        S = D["slow"]
        L = S["atds_params"]["window"]
        f = MathTex(r"L = 30 \times B", font_size=60).shift(UP * 2.3)
        f2 = MathTex(rf"= 30 \times {S['B']:.0f} \approx {L}\,\text{{s}}", font_size=48, color=C_ATDS).next_to(f, DOWN, buff=0.3)
        with self.beat("Step two. Make the window long enough to hold about thirty independent samples. So the window is thirty times B. And we take B from the slowest signal in the recording, so that every pair is measured with the same ruler."):
            self.play(Write(f), run_time=2)
            self.wait(2)
            self.play(Write(f2), run_time=2)
        ax = NumberLine(x_range=[0, 1500, 300], length=11, color=GRID, include_numbers=True, font_size=22,
                        decimal_number_config={"color": INK2, "num_decimal_places": 0}).shift(DOWN * 1.2)
        lbl = Text("seconds", font_size=20, color=INK2).next_to(ax, DOWN, buff=0.45)
        w = ValueTracker(60)
        bar = always_redraw(lambda: Rectangle(width=ax.n2p(w.get_value())[0] - ax.n2p(0)[0], height=0.6, color=C_ATDS, fill_opacity=0.5)
                            .move_to(ax.n2p(w.get_value() / 2) + UP * 0.6))
        b60 = Rectangle(width=ax.n2p(60)[0] - ax.n2p(0)[0], height=0.6, color=C_TDS, fill_opacity=0.7).move_to(ax.n2p(30) + UP * 1.4)
        t60 = Text("classic: 60 s", font_size=22, color=C_TDS).next_to(b60, UP, buff=0.1).align_to(b60, LEFT)
        tw = always_redraw(lambda: Text(f"aTDS: {w.get_value():.0f} s", font_size=22, color=C_ATDS).next_to(bar, RIGHT, buff=0.2))
        with self.beat(f"For our slow pair, that's a window of about {L / 60:.0f} minutes. Not one."):
            self.play(Create(ax), FadeIn(lbl), FadeIn(b60), FadeIn(t60), run_time=1)
            self.add(bar, tw)
            self.play(w.animate.set_value(L), run_time=3)
        cap = DashedLine(ax.n2p(300) + UP * 2.2, ax.n2p(300) + DOWN * 0.2, color=C_GOLD, stroke_width=4)
        ct = Text("cap you choose\n(5 min for sleep)", font_size=22, color=C_GOLD).next_to(cap, RIGHT, buff=0.15).align_to(cap, UP).shift(DOWN * 0.05)
        with self.beat("Of course, a longer window means coarser timing. So you set a cap: the longest window you're willing to accept. For the sleep data, we used five minutes."):
            self.play(Create(cap), FadeIn(ct), run_time=1.5)
            self.play(Indicate(ct, color=C_GOLD), run_time=1.2)
        self.play(FadeOut(Group(*self.mobjects)), run_time=0.8)


# ------------------------------------------------------------------ 7 lag + rhythm
class S07Lag(VO):
    def construct(self):
        self.intro(6, "step 3: the delay search")
        S = D["slow"]["atds_params"]
        nl = NumberLine(x_range=[-500, 500, 100], length=11, color=GRID, include_numbers=True, font_size=20,
                        decimal_number_config={"color": INK2, "num_decimal_places": 0}).shift(UP * 0.6)
        lab = Text("delays searched (s)", font_size=20, color=INK2).next_to(nl, DOWN, buff=0.45)
        r1 = Line(nl.n2p(-30), nl.n2p(30), color=C_TDS, stroke_width=14).shift(UP * 0.5)
        r2 = Line(nl.n2p(-S["max_lag"]), nl.n2p(S["max_lag"]), color=C_ATDS, stroke_width=14).shift(UP * 0.95)
        t1 = Text("classic: ±30 s", font_size=22, color=C_TDS).next_to(r1, RIGHT, buff=0.2)
        t2 = Text(f"aTDS: ±{S['max_lag']} s  (half the window)", font_size=22, color=C_ATDS).next_to(r2, UP, buff=0.12)
        with self.beat("Step three. The window moves forward by half its length, and the delay search grows with it, up to half the window. Slow systems can have long delays, and a fixed thirty second search would simply never find them."):
            self.play(Create(nl), FadeIn(lab), run_time=1)
            self.play(Create(r1), FadeIn(t1), run_time=1.2)
            self.play(Create(r2), FadeIn(t2), run_time=1.8)
        self.play(FadeOut(VGroup(nl, lab, r1, r2, t1, t2)), run_time=0.6)
        ax, _ = axes([0, 20, 5], [-1.5, 1.5, 1], 10, 1.8, nums=False)
        ax.shift(UP * 1.3)
        T = 4.0
        t = np.linspace(0, 20, 400)
        base = ax.plot_line_graph(t, np.sin(2 * PI * t / T), line_color=C_AQUA, add_vertex_dots=False, stroke_width=4)
        sh = ValueTracker(0)
        mov = always_redraw(lambda: ax.plot_line_graph(t, np.sin(2 * PI * (t - sh.get_value()) / T) * 0.95, line_color=C_GOLD,
                                                       add_vertex_dots=False, stroke_width=3))
        sl = always_redraw(lambda: MathTex(rf"\text{{shift}} = {sh.get_value():+.1f}\,\text{{s}}", font_size=36, color=C_GOLD)
                           .next_to(ax, DOWN, buff=0.3))
        per = Text("rhythm: repeats every T = 4 s", font_size=24, color=C_AQUA).next_to(ax, UP, buff=0.1)
        with self.beat("There's one exception, and that's rhythms. If a signal repeats every T seconds, then a delay of just over half a period looks exactly like a negative delay of just under half a period."):
            self.play(Create(base), FadeIn(per), run_time=1.2)
            self.add(mov, sl)
            self.play(sh.animate.set_value(2.6), run_time=2.5)
            self.wait(0.6)
            self.play(sh.animate.set_value(-1.4), run_time=0.01)
            self.wait(1.2)
            self.play(sh.animate.set_value(2.6), run_time=0.01)
            self.wait(1.2)
        f = MathTex(r"\text{search} \le \tfrac{T}{2}", font_size=54, color=C_ATDS).shift(DOWN * 2.2)
        with self.beat("See? The two overlays look the same. So the search stops at half of the shortest rhythm in the data."):
            self.play(Write(f), run_time=1.5)
        self.play(FadeOut(Group(*self.mobjects)), run_time=0.8)


# ------------------------------------------------------------------ 8 tolerance
class S08Tolerance(VO):
    def construct(self):
        self.intro(7, "step 4: the tolerance")
        rng = np.random.default_rng(7)
        bx, bg = axes([0, 30, 5], [-40, 40, 20], 9, 3.2, nums=False)
        bg.shift(UP * 1.0)
        lags = 10 + rng.normal(0, 4, 30)
        dots = VGroup(*[Dot(bx.c2p(i + 0.5, l), radius=0.07, color=C_ATDS) for i, l in enumerate(lags)])
        tol = ValueTracker(1)
        band = always_redraw(lambda: Rectangle(width=9, height=bx.c2p(0, 10 + tol.get_value())[1] - bx.c2p(0, 10 - tol.get_value())[1],
                                               color=C_ATDS, fill_opacity=0.2, stroke_width=0).move_to(bx.c2p(15, 10)))
        tl = always_redraw(lambda: MathTex(rf"\pm {tol.get_value():.0f}\,\text{{s}}", font_size=40, color=C_ATDS).next_to(bg, RIGHT, buff=0.3))
        with self.beat("Step four is the subtle one. The tolerance. With a long window, plus or minus one second is way too strict, because a real delay wanders a little. But if you make the tolerance too loose, chance alone starts to look like stability."):
            self.play(FadeIn(bg), LaggedStart(*[FadeIn(d) for d in dots], lag_ratio=0.03), run_time=2)
            self.add(band, tl)
            self.wait(1.5)
            self.play(tol.animate.set_value(8), run_time=2)
            self.play(tol.animate.set_value(35), run_time=2)
        self.play(FadeOut(VGroup(bg, dots, band, tl)), run_time=0.6)
        S = D["slow"]
        y = np.array(S["y"][:600])
        ax, _ = axes([0, 600, 100], [-3, 3, 3], 10.5, 1.6, nums=False)
        ax.shift(UP * 1.2)
        x = np.array(S["x"][:600])
        cx = curve(ax, x, C_AQUA, width=2)
        ax2, _ = axes([0, 600, 100], [-3, 3, 3], 10.5, 1.6, nums=False)
        ax2.shift(DOWN * 1.0)
        r = ValueTracker(0)
        cy = always_redraw(lambda: curve(ax2, np.roll(y, int(r.get_value())), C_GOLD, width=2))
        l1 = Text("signal 1", font_size=22, color=C_AQUA).next_to(ax, UP, buff=0.05).align_to(ax, LEFT)
        l2 = Text("signal 2, rotated", font_size=22, color=C_GOLD).next_to(ax2, UP, buff=0.05).align_to(ax2, LEFT)
        arrow = CurvedArrow(ax2.c2p(590, -2.6), ax2.c2p(10, -2.6), angle=-TAU / 6, color=INK2)
        with self.beat("So, adaptive T D S asks the data. It takes one signal and rotates it in time by a large random amount. Whatever falls off the end wraps around to the start. Each signal keeps its own shape and rhythm. But any real coupling between them is gone."):
            self.play(Create(cx), Create(ax), FadeIn(l1), run_time=1.2)
            self.add(cy)
            self.play(FadeIn(l2), Create(arrow), run_time=1)
            self.play(r.animate.set_value(260), run_time=4, rate_func=smooth)
        self.play(FadeOut(VGroup(ax, ax2, cx, cy, l1, l2, arrow)), run_time=0.6)
        nb = S["null_by_tol"]
        tt = np.array(sorted(map(int, nb)))
        vv = np.array([nb[str(i)] for i in tt])
        chosen = S["atds_params"]["tolerance"]
        cx_, cg = axes([0, int(tt.max()) + 5, 10], [0, max(8, float(vv.max()) + 1), 2], 9.5, 4.2, xlabel="tolerance (± s)",
                       ylabel="fake pairs looking coupled (%)")
        cg.shift(DOWN * 0.4)
        five = DashedLine(cx_.c2p(0, 5), cx_.c2p(tt.max() + 5, 5), color=INK2)
        ft = Text("5 % chance budget", font_size=22, color=INK2).next_to(cx_.c2p(2, 5), UP, buff=0.1).align_to(cx_.c2p(2, 5), LEFT)
        line = cx_.plot_line_graph(tt, vv, line_color=C_ATDS, add_vertex_dots=False, stroke_width=4)
        with self.beat("Then it runs T D S on these fake pairs, while slowly widening the tolerance. And it keeps the widest tolerance at which the fake pairs still look coupled no more than five percent of the time."):
            self.play(FadeIn(cg), Create(five), FadeIn(ft), run_time=1.5)
            self.play(Create(line), run_time=4, rate_func=linear)
        i = list(tt).index(chosen)
        pk = Dot(cx_.c2p(chosen, vv[i]), radius=0.14, color=C_ATDS)
        pl = Text(f"chosen: ±{chosen} s", font_size=30, color=C_ATDS).next_to(pk, LEFT, buff=0.4).shift(DOWN * 0.6)
        with self.beat(f"For our slow pair, that's plus or minus {chosen} seconds. That's wide. But it's earned. With it, chance stays at or below five percent."):
            self.play(GrowFromCenter(pk), Write(pl), run_time=1.5)
        self.play(FadeOut(Group(*self.mobjects)), run_time=0.8)


# ------------------------------------------------------------------ 9 example + ground truth
class S09Example(VO):
    def construct(self):
        self.intro(8, "putting it together")
        S = D["slow"]
        P = S["atds_params"]
        panels = []
        for i, (key, col, title, lim) in enumerate([("tds", C_TDS, "classic TDS  ·  60 s windows", 30),
                                                    ("atds", C_ATDS, f"aTDS  ·  {P['window']} s windows, ±{P['tolerance']} s", 60)]):
            T = S[key]
            bx, bg = axes([0, 14400, 3600], [-lim, lim, lim / 2], 5.6, 3.4, xlabel="time (h)", nums=False)
            bg.shift(LEFT * 3.2 + RIGHT * 6.4 * i + DOWN * 0.3)
            t = Text(title, font_size=22, color=col).next_to(bg, UP, buff=0.15)
            true = DashedLine(bx.c2p(0, 10), bx.c2p(14400, 10), color=C_GOLD)
            cen, lag, st = np.array(T["centre"]), np.array(T["lag"]), np.array(T["stable"])
            m = np.abs(lag) <= lim
            ds = VGroup(*[Dot(bx.c2p(c, l), radius=0.045 if key == "tds" else 0.09, color=col,
                              fill_opacity=1 if s_ else 0.25) for c, l, s_ in zip(cen[m], lag[m], st[m])])
            sc = Text(f"score {T['score']:.1f} %  ·  chance {S['null_mean'][key]:.1f} %", font_size=22, color=col).next_to(bg, DOWN, buff=0.2)
            yt = VGroup(*[Text(f"{v:+d} s" if v else "0", font_size=16, color=INK2).next_to(bx.c2p(0, v), LEFT, buff=0.1) for v in (-lim, 0, lim)])
            tl = Text("true delay 10 s", font_size=16, color=C_GOLD).next_to(bx.c2p(14400, 10), UP, buff=0.05).shift(LEFT * 0.8)
            bg.add(yt, tl)
            panels.append((bg, t, true, ds, sc))
        with self.beat(f"Now let's put it all together, on that same slow pair. Classic T D S, with one minute windows. The dots scatter. Adaptive T D S, with {P['window'] / 60:.0f} minute windows. And now, the dots lock onto the true ten second delay. Half of its windows are stable, against a chance level of zero."):
            for bg, t, true, ds, sc in panels:
                self.play(FadeIn(bg), FadeIn(t), Create(true), run_time=1)
                self.play(LaggedStart(*[FadeIn(d) for d in ds], lag_ratio=0.01), run_time=3)
                self.play(Write(sc), run_time=1)
        self.play(FadeOut(Group(*self.mobjects[1:])), run_time=0.6)
        E = D["e1"]
        r = np.array(E["fixed_default"]["r"])
        ax, g = axes([0, 4.5, 1], [0, 1, 0.25], 8.5, 4.4, ylabel="separation (AUC)", nums=False)
        g.add(Text("signals slower by ×", font_size=20, color=INK2).next_to(ax, DOWN, buff=0.6))
        g.shift(DOWN * 0.3 + LEFT * 1.3)
        xt = VGroup(*[Text(f"{v}×", font_size=20, color=INK2).next_to(ax.c2p(np.log2(v), 0), DOWN, buff=0.15) for v in r])
        yt = VGroup(*[Text(f"{v}", font_size=20, color=INK2).next_to(ax.c2p(0, v), LEFT, buff=0.15) for v in (0, 0.5, 1)])
        chance = DashedLine(ax.c2p(0, 0.5), ax.c2p(4.5, 0.5), color=GRID)
        lines = []
        for key, col, name, dash in [("oracle", C_NULL, "best possible", True), ("fixed_default", C_TDS, "classic TDS", False),
                                     ("v2_calibrated", C_ATDS, "adaptive TDS", False)]:
            ys = np.array(E[key]["auc"])
            ln = ax.plot_line_graph(np.log2(r), ys, line_color=col, add_vertex_dots=True, stroke_width=4,
                                    vertex_dot_radius=0.07)
            if dash:
                ln = DashedVMobject(ln["line_graph"], num_dashes=30)
            lab = Text(name, font_size=22, color=col).next_to(ax.c2p(np.log2(r[-1]), ys[-1]), RIGHT, buff=0.2)
            lines.append((ln, lab))
        with self.beat("That was just one example. But in earlier tests, with pairs whose coupling and delay were known, the pattern held. As the signals get slower, classic T D S falls to chance. Adaptive T D S stays close to the best possible setting. That's the one you'd pick if you already knew the answer."):
            self.play(FadeIn(g), FadeIn(xt), FadeIn(yt), Create(chance), run_time=1.2)
            for ln, lab in lines:
                self.play(Create(ln), FadeIn(lab), run_time=2.2)
        self.play(FadeOut(Group(*self.mobjects)), run_time=0.8)


# ------------------------------------------------------------------ 10 real data
class S10Real(VO):
    def construct(self):
        self.intro(9, "real data")
        C = D["sleep_stage_means"]
        stages = [("LS", "light"), ("awake", "wake"), ("REM", "REM"), ("DS", "deep")]
        grp = VGroup()
        for i, (v, col, name) in enumerate([("published", C_TDS, "classic TDS"), ("v2_calibrated", C_ATDS, "aTDS")]):
            mx = max(C[v].values())
            bars = VGroup()
            for j, (s, lab) in enumerate(stages):
                h = 3.0 * C[v][s] / mx
                b = Rectangle(width=0.9, height=h, color=col, fill_opacity=0.85, stroke_width=0)
                b.move_to(np.array([-5.4 + 6.4 * i + 1.25 * j, -2.2 + h / 2, 0]))
                t = Text(lab, font_size=20, color=INK2).next_to(b, DOWN, buff=0.12).set_y(-2.45)
                bars.add(VGroup(b, t))
            title = Text(name, font_size=26, color=col).move_to(np.array([-5.4 + 6.4 * i + 1.9, 1.45, 0]))
            grp.add(VGroup(bars, title))
        note = Text("Bashan 2012 cohort, 35 people · each method scaled to its own maximum", font_size=20, color=INK2).to_edge(DOWN, buff=0.25)
        with self.beat("What about real data? On sleep recordings from Bashan's study, thirty five people, both methods find the known result. The body's network is weakest in deep sleep. Here, adaptive T D S doesn't add much. And that makes sense, because the classic settings were designed for exactly this data."):
            for g in grp:
                self.play(FadeIn(g[1]), LaggedStart(*[GrowFromEdge(b[0], DOWN) for b in g[0]], lag_ratio=0.15),
                          FadeIn(VGroup(*[b[1] for b in g[0]])), run_time=2)
            self.play(FadeIn(note), run_time=0.8)
        self.play(FadeOut(Group(*self.mobjects[1:])), run_time=0.6)
        G = D["garmin"]
        keys = [k for k in G if k != "length_min"]
        bx, bg = axes([0, 1700, 300], [-120, 120, 60], 11, 4.2, xlabel="time in the run (s)", ylabel="speed leads HR by (s)")
        bg.shift(DOWN * 0.2)
        zero = DashedLine(bx.c2p(0, 0), bx.c2p(1700, 0), color=GRID)
        rows = []
        for k, col in zip(keys, (C_TDS, C_ATDS)):
            m = G[k]
            ds = VGroup(*[Dot(bx.c2p(c, l), radius=0.07 if col == C_TDS else 0.12, color=col, fill_opacity=1 if s_ else 0.2)
                          for c, l, s_ in zip(m["centre"], m["lag"], m["stable"]) if abs(l) <= 120])
            rows.append(ds)
        leg = VGroup(Text("app setting (120 s windows)", font_size=22, color=C_TDS),
                     Text("aTDS (240 s windows)", font_size=22, color=C_ATDS)).arrange(RIGHT, buff=0.8).to_edge(UP, buff=0.9)
        st = [l for l, s_ in zip(G[keys[1]]["lag"], G[keys[1]]["stable"]) if s_]
        with self.beat(f"And on one of your runs, the fixed app setting found no stable delay between speed and heart rate. Adaptive T D S found {len(st)} stable windows, with heart rate following speed by about {np.median(st):.0f} seconds."):
            self.play(FadeIn(bg), Create(zero), FadeIn(leg), run_time=1.2)
            self.play(LaggedStart(*[FadeIn(d) for d in rows[0]], lag_ratio=0.03), run_time=2)
            self.play(LaggedStart(*[FadeIn(d, scale=2) for d in rows[1]], lag_ratio=0.05), run_time=2)
            box = SurroundingRectangle(VGroup(*[d for d in rows[1] if d.get_fill_opacity() > 0.5]), color=C_ATDS, buff=0.15)
            self.play(Create(box), run_time=1)
        self.play(FadeOut(Group(*self.mobjects)), run_time=0.8)


# ------------------------------------------------------------------ 11 trade-offs + recap
class S11Recap(VO):
    def construct(self):
        self.intro(10, "the price, and a summary")
        nl = NumberLine(x_range=[0, 1800, 300], length=11, color=GRID, include_numbers=False).shift(UP * 0.3)
        segs = [(0, 700, "light", "#2b3542"), (700, 820, "wake", "#5a4a2a"), (820, 1800, "deep", "#34445a")]
        stripes = VGroup()
        for a, b, n, c in segs:
            r = Rectangle(width=nl.n2p(b)[0] - nl.n2p(a)[0], height=0.7, color=c, fill_opacity=1, stroke_width=0)
            r.move_to((nl.n2p(a) + nl.n2p(b)) / 2 + UP * 0.6)
            stripes.add(VGroup(r, Text(n, font_size=20).move_to(r)))
        w_short = VGroup(*[Rectangle(width=nl.n2p(60)[0] - nl.n2p(0)[0], height=0.3, color=C_TDS, stroke_width=2).move_to(nl.n2p(c) + DOWN * 0.3)
                           for c in range(30, 1800, 60)])
        w_long = VGroup(*[Rectangle(width=nl.n2p(300)[0] - nl.n2p(0)[0], height=0.3, color=C_ATDS, stroke_width=3).move_to(nl.n2p(c) + DOWN * 0.8)
                          for c in range(150, 1800, 300)])
        with self.beat("Of course, nothing is free. A long window blurs short events. A five minute window can't isolate a two minute stretch of being awake. And a wide tolerance blurs the exact delay."):
            self.play(FadeIn(stripes), run_time=1)
            self.play(Create(w_short), run_time=1.5)
            self.play(Create(w_long), run_time=1.5)
            self.play(Indicate(stripes[1], color=C_GOLD), run_time=1.5)
        self.play(FadeOut(VGroup(nl, stripes, w_short, w_long)), run_time=0.6)
        steps = VGroup(*[Text(s, font_size=30) for s in ["1.  Measure each signal's memory  (Bartlett factor B)",
                                                         "2.  Window = 30 × B of the slowest signal, up to your cap",
                                                         "3.  Step and delay search grow with the window",
                                                         "4.  Tolerance: widest with ≤ 5 % fake coupling"]]).arrange(DOWN, aligned_edge=LEFT, buff=0.45)
        steps.shift(DOWN * 0.2)
        for s, c in zip(steps, (C_GOLD, C_ATDS, C_ATDS, C_ATDS)):
            s[:2].set_color(c)
        with self.beat("So, let's recap. Adaptive T D S, in four steps. One: measure each signal's memory. Two: size the window to hold enough independent samples. Three: let the delay search grow with it. Four: calibrate the tolerance on pairs where the coupling has been destroyed. And then, choose the cap that matches the time resolution you need."):
            for s in steps:
                self.play(FadeIn(s, shift=RIGHT), run_time=1.2)
                self.wait(1.5)
        end = Text("adaptive Time Delay Stability", font_size=44, color=C_ATDS)
        with self.beat("And that's adaptive Time Delay Stability. Thanks for watching."):
            self.play(FadeOut(steps), run_time=0.6)
            self.play(Write(end), run_time=1.5)
        self.wait(1)
