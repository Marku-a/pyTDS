"""Synthesise an original calm ambient/lo-fi bed (numpy only) for the reel.

Output: reel_audio/music.wav (44.1 kHz stereo, length = video length, peak <= -1 dBFS).
Length is the sum of reel_audio/scene_durations.json (written by make_reel.py).
"""
import json
import wave
from pathlib import Path

import numpy as np
from scipy.signal import butter, sosfilt

SR, BPM = 44100, 90
A = Path(__file__).resolve().parent / "reel_audio"
rng = np.random.default_rng(7)
# Am - F - C - G ; pad voicings (Hz) and bass roots
CHORDS = [([220.00, 261.63, 329.63, 392.00], 55.00), ([174.61, 220.00, 261.63, 349.23], 43.65),
          ([196.00, 261.63, 329.63, 392.00], 65.41), ([196.00, 246.94, 293.66, 392.00], 49.00)]
BAR = 4 * 60 / BPM * 2  # 2 bars of 4/4 per chord = 5.33 s


def lp(x, fc, order=2):
    return sosfilt(butter(order, fc, "low", fs=SR, output="sos"), x, axis=-1)


def main():
    durs = json.load(open(A / "scene_durations.json"))
    total = round(sum(round(d * 30) / 30 for d in durs), 3)
    n = int(total * SR)
    t = np.arange(n) / SR
    pad = np.zeros((2, n))
    bass = np.zeros(n)
    nch = int(np.ceil(total / BAR)) + 1
    for k in range(nch):
        notes, root = CHORDS[k % 4]
        s0 = int(k * BAR * SR)
        L = min(int((BAR + 2.5) * SR), n - s0)  # overlap into next chord
        if L <= 0:
            break
        tt = np.arange(L) / SR
        env = np.minimum(1, tt / 1.6) * np.exp(-np.maximum(tt - BAR, 0) * 1.6)  # slow attack, release
        env *= 1 - 0.0 * tt
        for j, f in enumerate(notes):
            for d, side in [(-0.004, 0), (0.0, 0.5), (0.004, 1)]:  # detune, spread
                ff = f * (1 + d)
                ph = rng.uniform(0, 2 * np.pi)
                w = np.sin(2 * np.pi * ff * tt + ph) + 0.35 * (2 * ((ff * tt + ph / 6.283) % 1) - 1) * 0.5
                g = 0.05 * env * (1 + 0.15 * np.sin(2 * np.pi * 0.11 * tt + j))
                pad[0, s0:s0 + L] += w * g * (1 - side * 0.6)
                pad[1, s0:s0 + L] += w * g * (0.4 + side * 0.6)
        be = np.minimum(1, tt / 0.05) * np.exp(-np.maximum(tt - BAR, 0) * 4)
        bass[s0:s0 + L] += 0.30 * np.sin(2 * np.pi * root * tt) * be * (1 + 0.25 * np.sin(2 * np.pi * root * 2 * tt) * 0.3)
    pad = lp(pad, 1400)
    # muted hi-hat on 8ths + soft kick-like pulse on beats, after the hook
    beat = 60 / BPM
    hat = np.zeros(n)
    start = round(durs[0] * 30) / 30
    tk = start
    i = 0
    while tk < total - 1:
        k = int(tk * SR)
        m = min(int(0.09 * SR), n - k)
        nz = rng.standard_normal(m) * np.exp(-np.arange(m) / SR * 55)
        hat[k:k + m] += nz * (0.05 if i % 2 == 0 else 0.03)
        tk += beat / 2
        i += 1
    hat = lp(hat, 7000) - lp(hat, 2500)  # band-limited, muted
    pulse = np.zeros(n)
    tk = start
    while tk < total - 1:
        k = int(tk * SR)
        m = min(int(0.3 * SR), n - k)
        tt = np.arange(m) / SR
        pulse[k:k + m] += 0.18 * np.sin(2 * np.pi * (50 + 60 * np.exp(-tt * 30)) * tt) * np.exp(-tt * 12)
        tk += beat * 2
    ramp = np.clip((t - start) / 2.0, 0, 1)
    # swell at scene boundaries: filtered noise burst peaking at each boundary
    sw = np.zeros((2, n))
    tb = 0.0
    for d in durs[:-1]:
        tb += round(d * 30) / 30
        k0, k1 = int((tb - 1.2) * SR), int((tb + 0.4) * SR)
        env = np.concatenate([np.linspace(0, 1, int(1.2 * SR)) ** 2, np.linspace(1, 0, k1 - int(tb * SR))])
        sw[:, k0:k0 + len(env)] += (lp(rng.standard_normal((2, len(env))), 1800) * env) * 0.04
    mix = pad + bass * 0.9 + (hat * 0.7 + pulse * 0.6) * ramp + sw
    mix = mix.copy()
    mix[0] = mix[0]
    mix[1] = mix[1]
    mix *= np.minimum(1, t / 1.5) * np.minimum(1, (total - t) / 3.0).clip(0)  # 1.5 s in, 3 s out
    mix = np.tanh(mix * 1.2) / 1.2  # soft saturation, no hard clipping
    pk = np.abs(mix).max()
    mix *= 10 ** (-3 / 20) / pk  # peak -3 dBFS (well under -1)
    pcm = (mix.T * 32767).astype(np.int16)
    with wave.open(str(A / "music.wav"), "wb") as w:
        w.setnchannels(2)
        w.setsampwidth(2)
        w.setframerate(SR)
        w.writeframes(pcm.tobytes())
    print(f"music.wav {total:.2f}s peak {20 * np.log10(np.abs(mix).max()):.1f} dBFS rms {20 * np.log10(np.sqrt((mix ** 2).mean())):.1f} dBFS")


if __name__ == "__main__":
    main()
