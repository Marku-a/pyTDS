# aTDS paper

Manuscript: **Adaptive Time Delay Stability: letting the signals choose the window,
lag range and tolerance of TDS** (`atds_paper.tex` → `atds_paper.pdf`).

It is built from the approved aTDS brief in `../atds_brief/` (technical report,
experiment scripts, raw results, referee record). No experiment was re-run and
no number was added; every figure is drawn from `../atds_brief/results/*.json`.

## How it was written

The workflow follows [PaperOrchestra](https://github.com/google-research/paper-orchestra)
(Song et al. 2026), with its agents played by Claude Code instead of Gemini/OpenAI:

| PaperOrchestra step | Here |
|---|---|
| raw materials (`idea_sparse.md`, `experimental_log.md`) | `raw_materials/idea.md`, `raw_materials/experimental_log.md` — extracted verbatim from the approved report |
| plotting agent | `make_figures.py` → `figures/*.pdf` (from the JSON results only) |
| literature review agent (Semantic Scholar) | `references.bib` — every entry checked by web search, DOIs from index pages |
| section writing agent | `atds_paper.tex` |
| content refinement / reviewer | number-by-number fact check against the JSON results and code |

To run the real PaperOrchestra on the same inputs, pass
`--raw_materials_dir paper/raw_materials --idea_filename idea.md` plus a template
and your own Gemini/OpenAI and Semantic Scholar keys.

## Build

```bash
python paper/make_figures.py          # needs matplotlib
cd paper && pdflatex atds_paper && bibtex atds_paper && pdflatex atds_paper && pdflatex atds_paper
```

## Status

Draft. Evidence is synthetic plus one real EEG pair; reviewed by an AI referee
only, not yet by a human expert. Author affiliation and target journal still to
be filled in.
