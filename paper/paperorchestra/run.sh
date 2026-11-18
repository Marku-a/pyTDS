#!/usr/bin/env bash
# Run PaperOrchestra (google-research/paper-orchestra) on the aTDS raw materials.
# Usage: bash paper/paperorchestra/run.sh [work_dir]
#   GEMINI_API_KEY   real key, or any placeholder when a network secret for
#                    generativelanguage.googleapis.com (header x-goog-api-key)
#                    is injected by the session proxy.
#   WRITER_MODEL     default gemini-3.1-pro-preview (use a Flash model on free tier)
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
WORK="${1:-$HERE/work}"
PO_COMMIT=d9bae997e793a015ed38bb1605be9eb58c00284f
WRITER_MODEL="${WRITER_MODEL:-gemini-3.1-pro-preview}"
REFLECT_MODEL="${REFLECT_MODEL:-$WRITER_MODEL}"
export GEMINI_API_KEY="${GEMINI_API_KEY:-injected-by-proxy}"

mkdir -p "$WORK"
if [ ! -d "$WORK/paper-orchestra" ]; then
  git clone -q https://github.com/google-research/paper-orchestra.git "$WORK/paper-orchestra"
  git -C "$WORK/paper-orchestra" checkout -q "$PO_COMMIT"
  git -C "$WORK/paper-orchestra" apply "$HERE/s2_retry.patch"
fi
if [ ! -x "$WORK/venv/bin/python" ]; then
  python3 -m venv "$WORK/venv"
  "$WORK/venv/bin/pip" install -q -r "$WORK/paper-orchestra/requirements.txt"
fi
cd "$WORK/paper-orchestra"
"$WORK/venv/bin/python" paper_writing_cli.py \
  --raw_materials_dir "$HERE/raw_materials" \
  --latex_template_dir "$HERE/templates/netphys" \
  --output_dir "$HERE/output" \
  --idea_filename idea.md \
  --experimental_log_filename experimental_log.md \
  --writer_model_name "$WRITER_MODEL" \
  --reflection_model_name "$REFLECT_MODEL" \
  --research_cutoff "$(date +%Y-%m)" \
  --use_plotting false
