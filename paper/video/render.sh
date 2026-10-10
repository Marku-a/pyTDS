#!/usr/bin/env bash
# Render the 3Blue1Brown-style aTDS explainer, join the chapters, add subtitles.
# usage: VDATA=video_data.json TTS=kokoro KOKORO_DIR=<dir with kokoro-v1.0.onnx, voices-v1.0.bin> KOKORO_VOICE=af_heart VCACHE=vo MANIM=<manim> ./render.sh <l|m|h> <out base>
#   (or TTS=piper VOICE=<piper .onnx>)
set -euo pipefail
Q=${1:-m}; OUT=${2:-atds_3b1b}; MANIM=${MANIM:-manim}
SCENES="S01Hook S02XCorr S03Stability S04Assumption S05Memory S06Window S07Lag S08Tolerance S09Example S10Real S11Recap"
for s in $SCENES; do "$MANIM" -q"$Q" --disable_caching atds_3b1b.py "$s" >/dev/null; done
DIR=$(ls -d media/videos/atds_3b1b/*/ | head -1)
: > list.txt; for s in $SCENES; do echo "file '$DIR$s.mp4'" >> list.txt; done
ffmpeg -y -loglevel error -f concat -safe 0 -i list.txt -c:v copy -c:a aac -b:a 160k "${OUT}_nosubs.mp4"
python3 build_subs.py "$DIR" "$OUT" $SCENES
# animation scaled into the top 640 px; subtitles burned into the band below
ffmpeg -y -loglevel error -i "${OUT}_nosubs.mp4" -vf "scale=1138:640,pad=1280:720:71:0:color=0x141517,ass=${OUT}.ass" \
  -c:v libx264 -preset medium -crf 20 -pix_fmt yuv420p -c:a copy "${OUT}.mp4"
echo "${OUT}.mp4 ${OUT}.srt"
