# PaperOrchestra run 3: review record

Third full run of PaperOrchestra (commit `d9bae99`, patched as in `../paperorchestra.patch`, now with **no-revert mode**: every revision is kept). Inputs: raw materials with Parts 1-4 of the experimental log (simulations, real data, run-3 experiments E1-E5, reviewer follow-ups X6 and the held-out validation X7). Claude (session model) answered the writing and revision calls; **every reviewer and meta-reviewer call was answered by a fresh, independent Claude agent** that saw only PaperOrchestra's request files.

| Version | Review | Overall | Decision | Soundness | Presentation | Contribution | Originality | Significance | Clarity |
|---|---|---|---|---|---|---|---|---|---|
| v0 (draft) | reviewer 1 | 5 | Reject | 3 | 3 | 2 | 2 | 2 | 3 |
| v0 (draft) | reviewer 2 | 5 | Reject | 3 | 3 | 2 | 2 | 2 | 3 |
| v0 (draft) | reviewer 3 | 4 | Reject | 3 | 3 | 2 | 2 | 2 | 3 |
| v0 (draft) | meta | 5 | Reject | 3 | 3 | 2 | 2 | 2 | 3 |
| v1 | reviewer 1 | 5 | Reject | 3 | 2 | 2 | 2 | 2 | 2 |
| v1 | reviewer 2 | 5 | Reject | 3 | 2 | 2 | 2 | 2 | 2 |
| v1 | reviewer 3 | 5 | Reject | 3 | 3 | 2 | 2 | 2 | 3 |
| v1 | meta | 5 | Reject | 3 | 2 | 2 | 2 | 2 | 2 |
| v2 | reviewer 1 | 4 | Reject | 2 | 3 | 2 | 2 | 2 | 3 |
| v2 | reviewer 2 | 5 | Reject | 3 | 3 | 2 | 2 | 2 | 3 |
| v2 | reviewer 3 | 5 | Reject | 3 | 3 | 2 | 2 | 2 | 3 |
| v2 | meta | 4 | Reject | 2 | 3 | 2 | 2 | 2 | 3 |
| v3 (final) | reviewer 1 | 5 | Reject | 3 | 2 | 2 | 2 | 2 | 2 |
| v3 (final) | reviewer 2 | 5 | Reject | 3 | 2 | 2 | 2 | 2 | 2 |
| v3 (final) | reviewer 3 | 5 | Reject | 3 | 2 | 2 | 2 | 2 | 2 |
| v3 (final) | meta | 5 | Reject | 3 | 2 | 2 | 2 | 2 | 2 |

## What changed between versions

- **v0:** structure idea / simulated / real; new sections on full-record baselines, sensitivity, regime-switching/nonlinear/drift, stage-aware cohort, three catchments.
- **v1:** after the v0 reviews, new experiments X6 (estimated rescaling everywhere, boundary margin, failing-catchment diagnosis, window estimators, per-link power) and X7 (held-out validation of the Bayley-Hammersley window); aTDS-BH introduced and validated.
- **v2:** contribution repositioned as a family of data-driven rules; Cauchy-Schwarz justification of the window rule; honest labels for the BH factor; generalisation table; statement that the BH gain equals a longer window and that default aTDS never beats the estimated rescaling with statistical support; three tables moved to an appendix.
- **v3:** abstract and conclusion rewritten to report default aTDS first and aTDS-BH as a validated variant.

## Why the score stays at 4-5

All twelve independent reviews name the same ceiling: every ingredient is standard, and the simplest data-driven rule (the published setting rescaled by an estimated correlation time) captures most of the benefit. Clarity and presentation moved between 2 and 3 across versions (3 after v2, back to 2 for the longer v3); the added experiments could not change originality, contribution or significance (2 throughout). Raising the score further needs a substantively new result (e.g. a real dataset with known, time-varying delays where time resolution matters), not more writing.
