"""Mix narration + music with ffmpeg and mux with the re-timed silent reel.

Output: "adaptive-TDS reports/adaptive_TDS_reel_voice.mp4" (AAC 192k, 44.1 kHz stereo, ~-14 LUFS).
Run after make_voice.py, make_reel.py, make_music.py.
"""
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
A = ROOT / "experiments/adaptive/reel_audio"
VIDEO = ROOT / "adaptive-TDS reports/adaptive_TDS_reel.mp4"
OUT = ROOT / "adaptive-TDS reports/adaptive_TDS_reel_voice.mp4"


def main():
    ids = [s["id"] for s in json.load(open(ROOT / "experiments/adaptive/reel_voiceover.json"))["scenes"]]
    durs = json.load(open(A / "scene_durations.json"))
    cmd = ["ffmpeg", "-y", "-loglevel", "error", "-i", str(A / "music.wav")]
    parts, labels, t0 = [], [], 0.0
    for i, (sid, d) in enumerate(zip(ids, durs)):
        cmd += ["-i", str(next(A.glob(f"voice_{i:02d}_*.wav")))]
        ms = int(round((t0 + 0.3) * 1000))
        parts.append(f"[{i + 1}:a]aresample=44100,aformat=channel_layouts=stereo,adelay={ms}|{ms}[v{i}]")
        labels.append(f"[v{i}]")
        t0 += round(d * 30) / 30
    cmd += ["-i", str(VIDEO)]
    vid = len(ids) + 1
    f = ";".join(parts)
    f += f";{''.join(labels)}amix=inputs={len(ids)}:normalize=0:duration=longest,apad=whole_dur={t0:.3f},atrim=0:{t0:.3f},asplit=2[voice][key]"
    # music ~ -20 dB under voice; sidechain ducks it further while speaking
    f += (";[0:a]volume=-20dB[m];[m][key]sidechaincompress=threshold=0.02:ratio=6:attack=40:release=500:makeup=1[md]"
          ";[voice][md]amix=inputs=2:normalize=0,loudnorm=I=-14:TP=-1.5:LRA=11[a]")
    cmd += ["-filter_complex", f, "-map", f"{vid}:v", "-map", "[a]", "-c:v", "copy", "-c:a", "aac", "-b:a", "192k",
            "-ar", "44100", "-ac", "2", "-shortest", "-movflags", "+faststart", str(OUT)]
    subprocess.run(cmd, check=True)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
