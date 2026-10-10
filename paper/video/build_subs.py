"""Subtitles for the aTDS video from the per-sentence cues logged during rendering (cues/<Scene>.json).

Writes an .srt (sidecar) and an .ass (burned in below the animation, which is scaled into the top of the frame).
usage: python build_subs.py <scene video dir> <out base path, without extension> SCENE...
"""
import json
import os
import re
import subprocess
import sys

DISPLAY = [(r"\badaptive T D S\b", "adaptive TDS"), (r"\bT D S\b", "TDS"), (r"tau zero", "τ₀"), (r"plus or minus one second", "±1 second"), (r"plus or minus ", "±"),
           (r"\bcross correlation\b", "cross-correlation"), (r"\bone minute\b", "one-minute"), (r"(\d+) second delay", r"\1-second delay"),
           (r"sixty second window", "60-second window"), (r"thirty second (step|search)", r"30-second \1"), (r"five minute window", "five-minute window"),
           (r"two minute stretch", "two-minute stretch")]
MAXC = 78


def disp(t):
    for a, b in DISPLAY:
        t = re.sub(a, b, t)
    return t


def chunks(text, start, end):
    """Split a long sentence at commas/colons (else words) into pieces <= MAXC, timed by character share."""
    if len(text) <= MAXC:
        return [(start, end, text)]
    parts, cur = [], ""
    for tok in re.split(r"(?<=[,:;])\s+", text):
        if cur and len(cur) + 1 + len(tok) > MAXC:
            parts.append(cur)
            cur = tok
        else:
            cur = (cur + " " + tok).strip()
    parts.append(cur)
    out = []
    for p in parts:
        while len(p) > MAXC:
            cut = p.rfind(" ", 0, MAXC)
            out.append(p[:cut])
            p = p[cut + 1:]
        out.append(p)
    total = sum(len(p) for p in out)
    res, t = [], start
    for p in out:
        d = (end - start) * len(p) / total
        res.append((t, t + d, p))
        t += d
    return res


def ts(t, ass=False):
    h, m = int(t // 3600), int(t % 3600 // 60)
    s = t % 60
    return f"{h}:{m:02d}:{s:05.2f}" if ass else f"{h:02d}:{m:02d}:{int(s):02d},{int(round((s % 1) * 1000)) % 1000:03d}"


def main():
    vdir, base, scenes = sys.argv[1], sys.argv[2], sys.argv[3:]
    off, rows = 0.0, []
    for s in scenes:
        c = json.load(open(os.path.join("cues", s + ".json")))
        dur = float(subprocess.check_output(["ffprobe", "-v", "error", "-show_entries", "format=duration", "-of", "csv=p=0",
                                             os.path.join(vdir, s + ".mp4")]).decode())
        for a, b, text in c["cues"]:
            rows += chunks(disp(text), off + a, off + b)
        off += dur
    with open(base + ".srt", "w") as f:
        for i, (a, b, t) in enumerate(rows, 1):
            f.write(f"{i}\n{ts(a)} --> {ts(b)}\n{t}\n\n")
    head = """[Script Info]
ScriptType: v4.00+
PlayResX: 1280
PlayResY: 720

[V4+ Styles]
Format: Name, Fontname, Fontsize, PrimaryColour, SecondaryColour, OutlineColour, BackColour, Bold, Italic, Underline, StrikeOut, ScaleX, ScaleY, Spacing, Angle, BorderStyle, Outline, Shadow, Alignment, MarginL, MarginR, MarginV, Encoding
Style: Default,DejaVu Sans,27,&H00E4EAEC,&H000000FF,&H00171514,&H64171514,0,0,0,0,100,100,0,0,1,2,0,2,60,60,22,1

[Events]
Format: Layer, Start, End, Style, Name, MarginL, MarginR, MarginV, Effect, Text
"""
    with open(base + ".ass", "w") as f:
        f.write(head)
        for a, b, t in rows:
            f.write(f"Dialogue: 0,{ts(a, True)},{ts(b, True)},Default,,0,0,0,,{t}\n")
    print(len(rows), "subtitle lines,", round(off, 1), "s")


if __name__ == "__main__":
    main()
