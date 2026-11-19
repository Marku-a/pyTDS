# PaperOrchestra review record

Run of [PaperOrchestra](https://github.com/google-research/paper-orchestra) (commit `d9bae99`) on the aTDS raw materials, with Claude answering every model call through `manual_backend.py` (no Gemini/OpenAI key). The reviewer, meta-reviewer, refinement and formatting agents are PaperOrchestra's own prompts; Claude wrote the responses. Because the same model wrote and reviewed the paper, these scores are a self-assessment under PaperOrchestra's rubric, not an independent review.

## Scores per round (meta-review that PaperOrchestra saved; rubric 1-10, sub-axes 1-4)

| Round | Overall | Decision | Soundness | Presentation | Contribution | Clarity | Significance |
|---|---|---|---|---|---|---|---|
| Initial draft | 6 | Reject | 3 | 2 | 3 | 2 | 3 |
| Revision 1 | 6 | Accept | 3 | 3 | 3 | 3 | 3 |
| Revision 2 | 6 | Accept | 3 | 3 | 3 | 3 | 3 |
| Revision 3 (final) | 6 | Accept | 3 | 3 | 3 | 3 | 3 |

PaperOrchestra's loop log reports "Score for initial: 6/10" and accepted every revision ("Overall same, internal metrics stable/improved"); final score 6.0. Each round used a 3-reviewer ensemble (network physiology, statistics/signal processing, applied physiology) plus a meta-review.

## What changed in each round

- **Initial draft:** Draft from the section-writing agent.
- **Revision 1:** Justified the max-B window vs joint B_xy; added calibration cost; added "Relation to the oracle" and "Biases of the null model"; added Practical guidance and a summary-of-evidence table; framed the 8.3 % FDP as inconclusive; stated which design constants were not varied.
- **Revision 2:** Stated that even a 300-s window spans ten 30-s sleep epochs; added an automatic out-of-range check (share of stable windows at the lag edge); stated that per-link power was measured only in networks; added clustered SE and per-link power to Next steps.
- **Revision 3 (final):** Added a notation and sign-convention table (negative tau0 = x leads; positive reported delay = x leads).

## Weaknesses still open in the final meta-review

- No staged or multi-subject real data.
- Constants other than 30 not varied; null biases qualitative.

## Formatting check

The formatting agent inspected page images: all 6 figures and 4 tables fit the margins. Found: legends overlapping data in Figures 2b and 6a, and BibTeX casing/venue errors in references (TDSpy, Pearson, ROC, Hurst, Davies & Harte, Rulkov). These are outside the LaTeX source, so they were fixed after the run in `make_figures.py` and `final/references.bib`.

## Changes made to PaperOrchestra for this run (`paperorchestra.patch`)

- File-queue backend hook in `utils/gemini_utils.py` (Claude answers each call).
- Semantic Scholar: retry on 429, then fall back to Crossref (title/author/year/venue verification) because the public S2 tier was rate-limited.
- Literature agent: keep verified records that have no abstract in the index (marked as such) instead of dropping them; this kept 34 papers instead of 12.
- Google Search grounding (Gemini-only) replaced by Claude's own web search plus a reference list verified earlier by web search.
- Plotting agent (Gemini image model) not used; the six figures are drawn from the raw JSON results by `paper/make_figures.py`.

## Independent fact check after the run

A separate reviewer agent checked the final paper against the raw JSON results, the experiment log, the code and the bibliography. All numbers matched. Six wording problems in the new sections were fixed in `atds_paperorchestra.tex` after the run:

1. The 60 s / 30 s window setting is now cited to Rizzo 2020 only.
2. Removed an unsupported claim that other TDS studies "carried over" the original parameters.
3. The edge-of-lag-range check is now marked untested, and the text notes it would miss the at-chance fast-signal case.
4. The joint B_xy recommendation now states it rests on one configuration (both rules at AUC 1.00).
5. Rotation keeps autocorrelations exactly only in the circular sense.
6. "Often with the wrong sign" is now "in 16 of 50 cases".

## Real-data update (after the PaperOrchestra run)

Added in response to the reviewers' main weakness (no real staged data), with all TDS computed by pyTDS only (`experiments/atds_real/`):

- **Section 4.8, Table 4, Figure 7:** 35-subject polysomnography cohort from the author's MSc project (10 channels, 1 Hz, scored stages), per-stage networks with the cross-subject stage-matched null.
- **Section 4.9, Table 5:** three public records from hydrology, energy and fuel prices with known delays.
- **Benchmark:** pyTDS reproduces the MATLAB 43.6 % benchmark (43.53 %) on `data/data-*.txt`; the M001602.dat columns are a different series.
- **Revised:** abstract, Introduction scope, Limitations, Practical guidance, Not tested, Next steps, summary table, Conclusion, data availability.

A reviewer agent checked the code against `experiments/atds_real/SPEC.md` (one labelling fix to the benchmark output). A second reviewer fact-checked the new text against the results files: all table values matched; one interval typo and five wording overstatements were fixed. The largest was that the "REM below deep sleep" inversion holds for mean TDS only, not for significant links.
