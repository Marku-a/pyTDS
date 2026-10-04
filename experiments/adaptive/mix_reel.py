"""Mix narration + music with numpy/ffmpeg and mux with the re-timed silent reel.

Output (python mix_reel.py --part 1|2): "adaptive-TDS reports/adaptive_TDS_reel_part1_how_it_works.mp4" / ..._part2_nine_fields.mp4 (AAC 192k, 44.1 kHz stereo, -14 LUFS, TP -1.5 dB).
Run after make_voice.py, make_reel.py, make_music.py.

Loudness targets (measured on the separate stems, before the final loudnorm):
  voice stem                  -16 LUFS
  music in gaps (no speech)   -23 LUFS
  music during speech         -27 LUFS (music ducked ~4 dB while the narrator speaks)
Final mix: ffmpeg loudnorm (two pass, linear) to -14 LUFS integrated, true peak -1.5 dB.
"""
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy.io import wavfile
from scipy.signal import resample_poly

from make_voice import lufs, read_wav

ROOT = Path(__file__).resolve().parents[2]
A = ROOT / "experiments/adaptive/reel_audio"
PART = int(sys.argv[sys.argv.index("--part") + 1]) if "--part" in sys.argv else 1
VIDEO = A / f"silent_part{PART}.mp4"
OUT = ROOT / ("adaptive-TDS reports/adaptive_TDS_reel_part1_how_it_works.mp4" if PART == 1
              else "adaptive-TDS reports/adaptive_TDS_reel_part2_nine_fields.mp4")
SR = 44100
LEAD = 0.3
VOICE_LUFS, GAP_LUFS, SPEECH_LUFS = -16.0, -23.0, -27.0
DOWN, UP = 0.25, 0.35  # duck ramp before speech starts / after it ends (s)


def measure(x, tmp):
    p = Path(tmp) / "m.wav"
    wavfile.write(p, SR, x.T.astype(np.float32))
    return lufs(p)


def ramp_env(n, intervals, depth_db):
    """Gain envelope: 1 outside speech, -depth_db inside, linear ramps (DOWN before, UP after)."""
    t = np.arange(n) / SR
    d = np.zeros(n)
    for a, b in intervals:
        d = np.maximum(d, np.clip((t - (a - DOWN)) / DOWN, 0, 1) * np.clip(((b + UP) - t) / UP, 0, 1))
    return 10 ** (-depth_db * d / 20)


def masks(n, intervals):
    """Boolean masks of the flat speech part and the flat gap part (ramps excluded, fade-in/out excluded)."""
    t = np.arange(n) / SR
    sp = np.zeros(n, bool)
    gap = np.ones(n, bool)
    for a, b in intervals:
        sp |= (t >= a + 0.05) & (t <= b - 0.05)
        gap &= ~((t >= a - DOWN - 0.05) & (t <= b + UP + 0.05))
    gap &= (t > 1.6) & (t < n / SR - 3.1)
    return sp, gap


def loudnorm(inp, out, tp=-1.5):
    base = f"loudnorm=I=-14:TP={tp}:LRA=11"
    r = subprocess.run(["ffmpeg", "-nostats", "-i", str(inp), "-af", base + ":print_format=json", "-f", "null", "-"],
                       capture_output=True, text=True, check=True)
    js = json.loads(re.search(r"\{[^{}]*\"input_i\"[^{}]*\}", r.stderr, re.S).group(0))
    f = (f"{base}:measured_I={js['input_i']}:measured_TP={js['input_tp']}:measured_LRA={js['input_lra']}"
         f":measured_thresh={js['input_thresh']}:offset={js['target_offset']}:linear=true,aresample={SR}")
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-i", str(inp), "-af", f, "-c:a", "pcm_s16le", str(out)], check=True)


def final_stats(path):
    r = subprocess.run(["ffmpeg", "-nostats", "-i", str(path), "-af", "ebur128=peak=true:framelog=quiet", "-f", "null", "-"],
                       capture_output=True, text=True)
    i = float(re.findall(r"I:\s+(-?[\d.]+) LUFS", r.stderr)[-1])
    tp = float(re.findall(r"Peak:\s+(-?[\d.]+) dBFS", r.stderr)[-1])
    return i, tp


def main():
    sc = json.load(open(A / "sentences.json"))
    allids = [s["id"] for s in json.load(open(ROOT / "experiments/adaptive/reel_voiceover.json"))["scenes"]]
    ids = json.load(open(ROOT / "experiments/adaptive/reel_voiceover.json"))["parts"][f"part{PART}"]
    durs = json.load(open(A / f"scene_durations_part{PART}.json"))
    starts, t0 = [], 0.0
    for d in durs:
        starts.append(t0)
        t0 += round(d * 30) / 30
    total = t0
    n = int(round(total * SR))
    voice = np.zeros(n)
    intervals = []
    for i, sid in enumerate(ids):
        a, sr = read_wav(A / f"voice_{allids.index(sid):02d}_{sid}.wav")
        a = resample_poly(a, SR, sr)
        k = int(round((starts[i] + LEAD) * SR))
        voice[k:k + len(a)] += a[: n - k]
        s = sc[sid]["sentences"]
        # speech span of this scene's narration (speech to speech incl. pauses between sentences)
        intervals.append((starts[i] + LEAD + s[0]["start"], starts[i] + LEAD + s[-1]["end"]))
    voice = np.stack([voice, voice])
    music = wavfile.read(A / f"music_part{PART}.wav")[1].T.astype(np.float64) / 32768
    music = music[:, :n] if music.shape[1] >= n else np.pad(music, ((0, 0), (0, n - music.shape[1])))
    with tempfile.TemporaryDirectory() as td:
        voice *= 10 ** ((VOICE_LUFS - measure(voice, td)) / 20)
        sp, gap = masks(n, intervals)
        # music level: gaps at GAP_LUFS (stem unducked there), then pick duck depth so speech lands on SPEECH_LUFS
        music *= 10 ** ((GAP_LUFS - measure(music[:, gap], td)) / 20)
        depth = GAP_LUFS - SPEECH_LUFS
        for _ in range(3):
            mus = music * ramp_env(n, intervals, depth)
            err = measure(mus[:, sp], td) - SPEECH_LUFS
            depth = float(np.clip(depth + err, 3.0, 5.5))
        mus = music * ramp_env(n, intervals, depth)
        m_whole, m_sp, m_gap = measure(mus, td), measure(mus[:, sp], td), measure(mus[:, gap], td)
        v_whole = measure(voice, td)
        mix = voice + mus
        pk = np.abs(mix).max()
        print(f"duck depth {depth:.1f} dB; stems peak: voice {20 * np.log10(np.abs(voice).max()):.1f} dBFS, mix {20 * np.log10(pk):.1f} dBFS")
        pre = Path(td) / "pre.wav"
        wavfile.write(pre, SR, mix.T.astype(np.float32))
        for tp in (-1.5, -2.0, -3.0):  # AAC encode can overshoot; tighten the limiter if so
            normed = Path(td) / "norm.wav"
            loudnorm(pre, normed, tp)
            subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-i", str(normed), "-i", str(VIDEO), "-map", "1:v", "-map", "0:a",
                            "-c:v", "copy", "-c:a", "aac", "-b:a", "192k", "-ar", str(SR), "-ac", "2", "-shortest",
                            "-movflags", "+faststart", str(OUT)], check=True)
            fi, ftp = final_stats(OUT)
            if ftp <= -1.0:
                break
    print(f"LUFS  music stem (whole) {m_whole:.1f} | music during speech {m_sp:.1f} | music in gaps {m_gap:.1f} | "
          f"voice stem {v_whole:.1f} | final mix {fi:.1f} (true peak {ftp:.1f} dBFS)")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
