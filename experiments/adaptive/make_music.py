"""Synthesise an original lo-fi / ambient-pop bed (numpy only) for the reel.

Built to be audible on phone speakers: nearly all energy sits above 250 Hz.
  * 86 BPM, chords Am-F-C-G (2 bars each)
  * electric-piano style FM arpeggio (chord tones 220-1000 Hz, slight stereo spread)
  * warm pad (content 300-1500 Hz)
  * soft kick (with a mid-range click so phones reproduce it) and closed hat from ~4 s
Output: reel_audio/music_part<N>.wav (44.1 kHz stereo 16-bit, length = that part's video length, peak -3 dBFS).
Length is the sum of reel_audio/scene_durations_part<N>.json (written by make_reel.py).
  python make_music.py --part 1|2
Levels in the final mix are set by mix_reel.py (LUFS targets), not here.
"""
import json
import sys
import wave
from pathlib import Path

import numpy as np
from scipy.signal import butter, sosfilt

SR, BPM = 44100, 86
A = Path(__file__).resolve().parent / "reel_audio"
PART = int(sys.argv[sys.argv.index("--part") + 1]) if "--part" in sys.argv else 1
rng = np.random.default_rng(7)
BEAT = 60 / BPM
BAR = 4 * BEAT
CHORD_LEN = 2 * BAR  # 5.58 s
# Am - F - C - G : arpeggio tones (Hz, 220-1000) and pad voicing (Hz)
ARP = [[220.00, 261.63, 329.63, 440.00, 523.25, 659.25],
       [220.00, 261.63, 349.23, 440.00, 523.25, 698.46],
       [261.63, 329.63, 392.00, 523.25, 659.25, 783.99],
       [246.94, 293.66, 392.00, 493.88, 587.33, 783.99]]
PAD = [[261.63, 329.63, 440.00, 659.25], [261.63, 349.23, 440.00, 523.25],
       [329.63, 392.00, 523.25, 659.25], [293.66, 392.00, 493.88, 587.33]]
PATTERN = [0, 2, 3, 4, 3, 2, 3, 5]  # 8th-note order over the 6 tones


def filt(x, kind, fc, order=2):
    return sosfilt(butter(order, fc, kind, fs=SR, output="sos"), x, axis=-1)


def epiano(f, dur=1.6):
    """FM electric piano note: 1:1 body with decaying index + short bell-like attack."""
    t = np.arange(int(dur * SR)) / SR
    idx = 1.6 * np.exp(-t * 4.5) + 0.3
    body = np.sin(2 * np.pi * f * t + idx * np.sin(2 * np.pi * f * t))
    bell = 0.35 * np.sin(2 * np.pi * f * 4 * t + 0.8 * np.sin(2 * np.pi * f * 7 * t)) * np.exp(-t * 14)
    env = np.minimum(1, t / 0.004) * np.exp(-t * 2.6)
    return (body + bell) * env


def add(buf, x, k, pan):
    n = min(len(x), buf.shape[1] - k)
    if n <= 0:
        return
    buf[0, k:k + n] += x[:n] * (1 - pan)
    buf[1, k:k + n] += x[:n] * pan


def main():
    durs = json.load(open(A / f"scene_durations_part{PART}.json"))
    total = round(sum(round(d * 30) / 30 for d in durs), 3)
    n = int(total * SR)
    t = np.arange(n) / SR
    arp = np.zeros((2, n))
    pad = np.zeros((2, n))
    drums = np.zeros((2, n))
    # --- arpeggio: 8th notes, slight pan alternation, velocity accents
    nch = int(np.ceil(total / CHORD_LEN)) + 1
    notes = {}
    for k in range(nch):
        for j in range(16):  # 16 eighths per chord
            tk = k * CHORD_LEN + j * BEAT / 2
            if tk >= total:
                break
            f = ARP[k % 4][PATTERN[j % 8]]
            if f not in notes:
                notes[f] = epiano(f)
            vel = (1.0 if j % 4 == 0 else 0.7 if j % 2 == 0 else 0.55) * rng.uniform(0.9, 1.0)
            pan = 0.5 + 0.22 * np.sin(j * 1.7 + k)  # slight stereo spread
            add(arp, notes[f] * vel, int(tk * SR), pan)
    # --- pad: detuned saws, slow attack, overlapping chords, band-limited 300-1500 Hz
    for k in range(nch):
        s0 = int(k * CHORD_LEN * SR)
        L = min(int((CHORD_LEN + 2.0) * SR), n - s0)
        if L <= 0:
            break
        tt = np.arange(L) / SR
        env = np.minimum(1, tt / 1.2) * np.exp(-np.maximum(tt - CHORD_LEN, 0) * 1.8)
        for j, f in enumerate(PAD[k % 4]):
            for d, pan in [(-0.003, 0.2), (0.0, 0.5), (0.003, 0.8)]:
                ff = f * (1 + d)
                ph = rng.uniform(0, 1)
                w = 2 * ((ff * tt + ph) % 1) - 1
                w = 0.8 * w + 0.5 * np.sin(2 * np.pi * ff * tt)
                g = 0.05 * env * (1 + 0.12 * np.sin(2 * np.pi * 0.13 * tt + j))
                pad[0, s0:s0 + L] += w * g * (1 - pan)
                pad[1, s0:s0 + L] += w * g * pan
    pad = filt(filt(pad, "low", 1500), "high", 300)
    # --- drums from ~4 s: soft kick on beats 1 and 3 (+ click), closed hat on off-8ths
    start = 4.0
    k_t = start
    bi = 0
    while k_t < total - 0.5:
        k = int(k_t * SR)
        m = min(int(0.35 * SR), n - k)
        tt = np.arange(m) / SR
        body = np.sin(2 * np.pi * np.cumsum(75 + 110 * np.exp(-tt * 30)) / SR) * np.exp(-tt * 11)
        body = np.tanh(2.5 * body)  # saturation adds harmonics phone speakers can reproduce
        click = filt(rng.standard_normal(m), "band", [900, 3500]) * np.exp(-tt * 90) * 0.6
        v = 1.0 if bi % 2 == 0 else 0.8
        drums[:, k:k + m] += ((body * 0.28 + click * 0.8) * v)[None, :]
        k_t += 2 * BEAT
        bi += 1
    h_t = start + BEAT / 2
    i = 0
    while h_t < total - 0.5:
        k = int(h_t * SR)
        m = min(int(0.07 * SR), n - k)
        nz = filt(rng.standard_normal(m), "high", 6000) * np.exp(-np.arange(m) / SR * 70)
        drums[0, k:k + m] += nz * 0.16 * (0.5 + 0.5 * (i % 2 == 0)) * 0.9
        drums[1, k:k + m] += nz * 0.16 * (0.5 + 0.5 * (i % 2 == 0)) * 1.1
        h_t += BEAT
        i += 1
    ramp = np.clip((t - start + 0.5) / 1.5, 0, 1)
    mix = arp * 0.62 + pad * 0.9 + drums * ramp
    mix *= np.minimum(1, t / 1.5) * np.clip((total - t) / 3.0, 0, 1)  # 1.5 s in, 3 s out
    mix = np.tanh(mix * 1.1) / 1.1
    mix *= 10 ** (-3 / 20) / np.abs(mix).max()
    pcm = (mix.T * 32767).astype(np.int16)
    with wave.open(str(A / f"music_part{PART}.wav"), "wb") as w:
        w.setnchannels(2)
        w.setsampwidth(2)
        w.setframerate(SR)
        w.writeframes(pcm.tobytes())
    mono = mix.mean(0)
    X = np.abs(np.fft.rfft(mono)) ** 2
    fr = np.fft.rfftfreq(n, 1 / SR)
    frac = X[fr > 250].sum() / X.sum()
    print(f"music_part{PART}.wav {total:.2f}s peak {20 * np.log10(np.abs(mix).max()):.1f} dBFS; energy above 250 Hz: {100 * frac:.1f}%")


if __name__ == "__main__":
    main()
