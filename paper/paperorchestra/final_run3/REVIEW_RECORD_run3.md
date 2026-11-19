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

## Final version: aTDS-BH removed (author decision)

At the author's request the Bayley-Hammersley-type window variant (aTDS-BH) was removed from the paper after the PaperOrchestra loop. Its section was replaced by a held-out replication of the main benchmark and the other-generator comparison (published, estimated rescaling and aTDS only). A reviewer agent then fact-checked the edited paper against the results files; six fixes were applied (above all, "aTDS is equivalent to the estimated rescaling" was softened to "performs similarly; equivalence not established", and the boundary-margin claim was qualified because controls stayed at 5.5 %). The edited version was scored by a fresh independent review round with PaperOrchestra's reviewer prompts:

| Review | Overall | Decision | Soundness | Presentation | Contribution | Originality | Significance | Clarity |
|---|---|---|---|---|---|---|---|---|
| review1 | 4 | Reject | 3 | 2 | 2 | 2 | 2 | 2 |
| review2 | 5 | Reject | 3 | 2 | 2 | 2 | 2 | 2 |
| review3 | 5 | Reject | 3 | 2 | 2 | 2 | 2 | 2 |
| meta | 5 | Reject | 3 | 2 | 2 | 2 | 2 | 2 |

The main weaknesses named:

- The paper does not show a benefit over a trivial baseline. By the authors' own analysis, no synthetic test shows aTDS beating the estimated-r rescaling with statistical support. The rescaling has higher per-link power at r = 4 (0.82 vs 0.66) and lower false stability under regime switching (2.0% vs 8.1%). The large gains are measured against the published setting used outside its design regime.
- The equivalence claim in the abstract and conclusion is overstated. The per-link power results argue against strict equivalence, and one reviewer questions whether equivalence was formally tested on the held-out set and on every comparison table row.
- The advertised additions (explicit false-alarm budget, rhythm cap, single mixed-speed rule) are not shown to add value. On the d = 30 rhythm, the cap makes aTDS report the wrong (anti-phase) delay where the rescaling succeeds.
- Calibrated tolerances grow to about 5-10% of the lag range (up to +/-56 samples) and can exceed the delay itself. This coarsens the meaning of 'stable', makes aTDS scores incomparable with the TDS literature, and some of the detection gain may come from the looser criterion. The oracle prefers tolerance 0 with longer windows.
- The false-stability budget controls only the mean chance score. Single recordings reach a 95th percentile of 20%, regime boundaries give 8-9%, the SE treats correlated rotations as independent, and the rotation null ignores shared slow co-modulation. These biases favour wider tolerances and smaller p-values on non-stationary physiological data. Block, IAAFT or multivariate nulls were not evaluated.
- The window rule is fragile. Truncating the Bartlett sum at the first zero crossing ignores seasonal and negative lobes, which causes the L0123002 failure, and makes the window record-length dependent for long-memory signals. Per-signal truncation also weakens the Cauchy-Schwarz justification. The constant 30 is a convention, and 45-60 detect better.
- There is no improvement in the target domain. On 35 polysomnograms aTDS does not improve sleep-stage contrasts, the uncapped window cannot resolve stages, and the wake contrast weakens. There is no comparison with time-resolved alternatives such as sliding-window GCC-PHAT.
- Real-data evidence of an advantage is thin. It rests on one hourly rainfall-runoff record whose reference delay is partly defined by full-record cross-correlation, and several real analyses use very few windows (12-19).
- All main results use the non-causal stability anchor, which changes the benchmark score from 43.5% to 58.1%, and robustness to the causal anchor is untested. Long delays under circular correlation rely on wrap-around overlap, and a zero-padded alternative is not discussed.
- Originality is modest: the ingredients (effective sample size, surrogate calibration, BH-FDR) are standard, and the novelty lies in combining them for TDS.
- Presentation is dense and repetitive (duplicated guidance bullets, numbers restated across the text and tables), which dilutes the core message. The sleep cohort's source, inclusion criteria and ethics approval or consent are not documented.
