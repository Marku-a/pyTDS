#!/usr/bin/env bash
# Render the 3Blue1Brown-style aTDS explainer and join the chapters.
# usage: VDATA=video_data.json VOICE=en-us-lessac-medium.onnx VCACHE=vo MANIM=<manim binary> ./render.sh <quality: l|m|h> <out.mp4>
set -euo pipefail
Q=${1:-m}; OUT=${2:-atds_3b1b.mp4}; MANIM=${MANIM:-manim}
SCENES="S01Hook S02XCorr S03Stability S04Assumption S05Memory S06Window S07Lag S08Tolerance S09Example S10Real S11Recap"
for s in $SCENES; do "$MANIM" -q"$Q" --disable_caching atds_3b1b.py "$s" >/dev/null; done
DIR=$(ls -d media/videos/atds_3b1b/*/ | head -1)
: > list.txt; for s in $SCENES; do echo "file '$DIR$s.mp4'" >> list.txt; done
ffmpeg -y -loglevel error -f concat -safe 0 -i list.txt -c:v copy -c:a aac -b:a 160k "$OUT"
echo "$OUT"
