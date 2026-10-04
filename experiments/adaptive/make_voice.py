"""Synthesise one narration WAV per reel scene with piper (LibriTTS high, CC BY 4.0).

Voice files (reel_audio/en-us-libritts-high.onnx[.json]) come from
https://github.com/rhasspy/piper/releases/download/v0.0.2/voice-en-us-libritts-high.tar.gz
Writes reel_audio/voice_<NN>_<id>.wav and reel_audio/durations.json.
"""
import json
import wave
from pathlib import Path

import numpy as np
from piper import PiperVoice, SynthesisConfig

D = Path(__file__).resolve().parent
A = D / "reel_audio"
SPEAKER, LENGTH_SCALE = "4535", 1.05


def main():
    voice = PiperVoice.load(str(A / "en-us-libritts-high.onnx"))
    sid = json.load(open(A / "en-us-libritts-high.onnx.json"))["speaker_id_map"][SPEAKER]
    cfg = SynthesisConfig(speaker_id=sid, length_scale=LENGTH_SCALE)
    out = {}
    for i, s in enumerate(json.load(open(D / "reel_voiceover.json"))["scenes"]):
        a = np.concatenate([c.audio_float_array for c in voice.synthesize(s["say"], cfg)])
        a = a * (0.89 / max(np.abs(a).max(), 1e-9))  # piper peak-normalises to 1.0; leave headroom
        pcm = (np.clip(a, -1, 1) * 32767).astype(np.int16)
        with wave.open(str(A / f"voice_{i:02d}_{s['id']}.wav"), "wb") as w:
            w.setnchannels(1)
            w.setsampwidth(2)
            w.setframerate(voice.config.sample_rate)
            w.writeframes(pcm.tobytes())
        out[s["id"]] = len(a) / voice.config.sample_rate
        print(f"{s['id']}: {out[s['id']]:.2f}s peak={np.abs(a).max():.3f} rms={np.sqrt((a**2).mean()):.3f}")
    json.dump(out, open(A / "durations.json", "w"), indent=1)


if __name__ == "__main__":
    main()
