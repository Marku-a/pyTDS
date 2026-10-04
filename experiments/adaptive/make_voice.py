"""Synthesise one narration WAV per reel scene with piper, then clean it with an ffmpeg chain.

Voice: en-us-libritts-high (CC BY 4.0), speaker/length_scale chosen by measurement (see CREDITS.md).
Voice files (reel_audio/en-us-libritts-high.onnx[.json]) come from
https://github.com/rhasspy/piper/releases/download/v0.0.2/voice-en-us-libritts-high.tar.gz
Each scene is split into sentences (0.45 s of silence inserted between them); sentences longer than
MAX_CHARS are further split at commas/colons (0.15 s) so on-screen subtitles fit two lines.
Writes reel_audio/voice_<NN>_<id>.wav (processed, ~-16 LUFS each), reel_audio/durations.json and
reel_audio/sentences.json ({id: {dur, sentences: [{text, start, end, chunks: [{text, start, end}]}]}}).
  python make_voice.py          synthesise + process
  python make_voice.py --eval   word error rate (pocketsphinx) of the processed WAVs vs the script
"""
import json
import re
import subprocess
import sys
import tempfile
import wave
from pathlib import Path

import numpy as np

D = Path(__file__).resolve().parent
A = D / "reel_audio"
SPEAKER, LENGTH_SCALE = "4535", 1.45
NOISE_SCALE, NOISE_W = 0.5, 0.6
GAP_SENT, GAP_CLAUSE, MAX_CHARS = 0.45, 0.15, 80
CHAIN = ("highpass=f=80,equalizer=f=3000:t=q:w=1:g=2.5,"
         "acompressor=threshold=-20dB:ratio=3:attack=10:release=120:makeup=2")
TARGET_LUFS = -16.0


def split_sentences(text):
    return [s for s in re.split(r"(?<=[.?!])\s+", text.strip()) if s]


def split_chunks(sent):
    """Split a long sentence at commas/colons into pieces of <= MAX_CHARS where possible."""
    if len(sent) <= MAX_CHARS:
        return [sent]
    parts = re.split(r"(?<=[,:])\s+", sent)
    out, cur = [], ""
    for p in parts:
        if cur and len(cur) + 1 + len(p) > MAX_CHARS:
            out.append(cur)
            cur = p
        else:
            cur = (cur + " " + p).strip()
    return out + [cur]


def lufs(path):
    """Integrated loudness (LUFS) of an audio file via ffmpeg ebur128."""
    r = subprocess.run(["ffmpeg", "-nostats", "-i", str(path), "-af", "ebur128=framelog=quiet", "-f", "null", "-"],
                       capture_output=True, text=True)
    m = re.findall(r"I:\s+(-?[\d.]+|-inf) LUFS", r.stderr)
    return float(m[-1])


def write_wav(path, a, sr):
    pcm = (np.clip(a, -1, 1) * 32767).astype(np.int16)
    with wave.open(str(path), "wb") as w:
        w.setnchannels(1)
        w.setsampwidth(2)
        w.setframerate(sr)
        w.writeframes(pcm.tobytes())


def read_wav(path):
    with wave.open(str(path), "rb") as w:
        sr = w.getframerate()
        return np.frombuffer(w.readframes(w.getnframes()), np.int16).astype(np.float64) / 32767, sr


def synth_scene(voice, cfg, say, sr):
    """Return (audio, sentences) with per-sentence and per-chunk times in seconds."""
    parts, sents, pos = [], [], 0
    ss = split_sentences(say)
    for si, s in enumerate(ss):
        chunks = split_chunks(s)
        sent = {"text": s, "start": pos / sr, "chunks": []}
        for ci, ch in enumerate(chunks):
            a = np.concatenate([c.audio_float_array for c in voice.synthesize(ch, cfg)])
            last_in_sent = ci == len(chunks) - 1
            gap = 0 if (last_in_sent and si == len(ss) - 1) else int((GAP_SENT if last_in_sent else GAP_CLAUSE) * sr)
            sent["chunks"].append({"text": ch, "start": pos / sr, "end": (pos + len(a)) / sr})
            parts += [a, np.zeros(gap)]
            pos += len(a) + gap
        sent["end"] = sent["chunks"][-1]["end"]
        sents.append(sent)
    return np.concatenate(parts), sents


def process(a, sr, tmp):
    """highpass + presence boost + light compression, then gain to TARGET_LUFS (peak-safe)."""
    raw, out = tmp / "raw.wav", tmp / "proc.wav"
    write_wav(raw, a * (0.89 / max(np.abs(a).max(), 1e-9)), sr)
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-i", str(raw), "-af", CHAIN, "-ar", str(sr), "-ac", "1",
                    "-c:a", "pcm_s16le", str(out)], check=True)
    p, _ = read_wav(out)
    g = 10 ** ((TARGET_LUFS - lufs(out)) / 20)
    p = p * g
    pk = np.abs(p).max()
    if pk > 0.97:  # gentle peak safety (compressed speech rarely needs it)
        p = p * (0.97 / pk)
    return p


def main():
    from piper import PiperVoice, SynthesisConfig
    voice = PiperVoice.load(str(A / "en-us-libritts-high.onnx"))
    sr = voice.config.sample_rate
    sid = json.load(open(A / "en-us-libritts-high.onnx.json"))["speaker_id_map"][SPEAKER]
    cfg = SynthesisConfig(speaker_id=sid, length_scale=LENGTH_SCALE, noise_scale=NOISE_SCALE, noise_w_scale=NOISE_W)
    durs, info = {}, {}
    with tempfile.TemporaryDirectory() as td:
        for i, s in enumerate(json.load(open(D / "reel_voiceover.json"))["scenes"]):
            a, sents = synth_scene(voice, cfg, s["say"], sr)
            p = process(a, sr, Path(td))
            write_wav(A / f"voice_{i:02d}_{s['id']}.wav", p, sr)
            durs[s["id"]] = len(p) / sr
            info[s["id"]] = {"dur": len(p) / sr, "sentences": sents}
            nw = len(s["say"].split())
            print(f"{s['id']}: {durs[s['id']]:.2f}s {nw / durs[s['id']]:.2f} words/s peak={np.abs(p).max():.3f}", flush=True)
    json.dump(durs, open(A / "durations.json", "w"), indent=1)
    json.dump(info, open(A / "sentences.json", "w"), indent=1)


def evaluate():
    from pocketsphinx import Decoder
    from scipy.signal import resample_poly
    dec = Decoder(samprate=16000)
    norm = lambda t: re.sub(r"[^a-z' ]", " ", t.lower()).split()
    rows = []
    for i, s in enumerate(json.load(open(D / "reel_voiceover.json"))["scenes"]):
        a, sr = read_wav(next(A.glob(f"voice_{i:02d}_*.wav")))
        a16 = resample_poly(a, 16000, sr)
        dec.start_utt()
        dec.process_raw((a16 * 32767).astype(np.int16).tobytes(), full_utt=True)
        dec.end_utt()
        r, h = norm(s["say"]), norm(dec.hyp().hypstr if dec.hyp() else "")
        d = np.arange(len(h) + 1)
        for x in range(1, len(r) + 1):
            p, d = d, np.zeros(len(h) + 1, int)
            d[0] = x
            for y in range(1, len(h) + 1):
                d[y] = min(p[y] + 1, d[y - 1] + 1, p[y - 1] + (r[x - 1] != h[y - 1]))
        rows.append(d[-1] / len(r))
        print(f"{s['id']}: WER {rows[-1]:.2f}")
    print(f"mean WER over {len(rows)} scenes: {np.mean(rows):.3f}")


if __name__ == "__main__":
    evaluate() if "--eval" in sys.argv else main()
