# PaperOrchestra run 2: review record

Second full run of [PaperOrchestra](https://github.com/google-research/paper-orchestra) (commit `d9bae99`, patched as in `../paperorchestra.patch`) on raw materials that now include the real-data results, with the structure: idea and algorithm, results on simulated data, results on real data, conclusion.

**Who answered which call.** Claude (the session model) answered the writing calls: outline, literature, Introduction and Related Work, sections, revision and formatting. **Every reviewer and meta-reviewer call was answered by a fresh, independent Claude agent** that saw only PaperOrchestra's request files (its reviewer instructions plus the paper text), with no access to this conversation, the repository or earlier reviews.

## Scores

| Round | Review | Overall | Decision | Soundness | Presentation | Contribution | Originality | Significance | Confidence |
|---|---|---|---|---|---|---|---|---|---|
| Draft (v0) | reviewer1 | 6 | Accept | 3 | 3 | 2 | 2 | 2 | 4 |
| Draft (v0) | reviewer2 | 5 | Reject | 3 | 3 | 2 | 2 | 2 | 4 |
| Draft (v0) | reviewer3 | 5 | Reject | 3 | 2 | 2 | 2 | 2 | 4 |
| Draft (v0) | meta | 6 | Accept | 3 | 3 | 2 | 2 | 2 | 3 |
| Revision 1 | reviewer1 | 5 | Reject | 3 | 3 | 2 | 2 | 2 | 4 |
| Revision 1 | reviewer2 | 5 | Reject | 3 | 3 | 2 | 2 | 2 | 4 |
| Revision 1 | reviewer3 | 5 | Reject | 3 | 3 | 2 | 2 | 2 | 4 |
| Revision 1 | meta | 5 | Reject | 3 | 3 | 2 | 2 | 2 | 4 |

PaperOrchestra's loop: v0 meta 6 → v1 meta 5, so it logged "STOP: Overall Score decreased. Reverting to previous iteration" and kept v0. The final paper (`atds_paper_run2.pdf`) is v0 plus formatting fixes (version labels removed, Bartlett 1946 reference corrected). **Final score: 6.0.**

## Revision 1 (kept as an alternative)

`alternative_revision1/` holds the author revision answering the v0 meta-review: tolerance as a fraction of the lag range (4.5–10.4 %), the 300-s cap stated as the user's choice, why widen the tolerance rather than lengthen the window, the REM change no longer counted as an improvement (no ground truth), and untested constants and missing GCC/coherence baselines listed as limitations. It is more cautious and more accurate on the REM point, but scored 5 (all three reviewers 5, Reject).

## Weaknesses in the last meta-review (what would raise the score)

- Limited methodological novelty: every ingredient (Bartlett effective sample size, circular-correlation aliasing bound, circular-shift surrogates, BH-FDR) is standard; the contribution is mainly their combination for one coupling measure.
- Modest gain over the fair baseline: against the setting rescaled by the true r, aTDS is equivalent or undecided except at r = 8, and the oracle (longer windows, tol = 0) beats it at every r >= 2. The benefit is mainly automation rather than better detection.
- The main synthetic benchmark scales the delay together with the correlation time (linear AR(1), constant delay, permanent coupling). This matches the assumptions built into the window and lag rules, so success there is partly by construction. Long-delay/short-memory, nonlinear, intermittent and regime-switching coupling are barely tested, and the long-delay case fails silently without a user-supplied dmax.
- The motivating physiological application is essentially negative. The uncapped 48-min window cannot resolve sleep stages. The manual 300-s cap undermines the parameter-free claim, and with it there is no significant gain and the wake significant-link contrast becomes significantly weaker.
- Calibrated tolerances widen to 5–10% of the lag range (up to ±56 samples, sometimes larger than the delay). This changes what 'stable delay' means, makes aTDS scores incomparable with the published TDS literature (undercutting the reproducibility motivation), and brings no significant AUC benefit.
- The rotation null is biased toward leniency: it removes shared slow co-modulation, and the SE treats correlated rotations as independent. Single recordings reach about 20% chance score, and the sleep-cohort tolerance is not calibrated with the stage-matched null.
- Missing simple baselines: estimated-time-scale rescaling or downsampling of the published setting, full-record cross-correlation with a surrogate test, and standard delay estimators (GCC-PHAT, coherence-based).
- Evidence from other fields rests on single records with very long windows (27–80 days), and the 'known' delays are partly defined by the full-record cross-correlation peak itself.
- Design constants (5% budget, floor of 61, 2-SE margin, 20 rotations, 0.3 rhythm threshold, first-zero truncation of the Bartlett sum) are never varied. Train/test separation in the synthetic tests is weak. The network FDR experiment is small (10 networks per design), and its 8.3% FDP is inconclusive.
- The default non-causal stability anchor departs from the original method (58.1% vs 43.6% on the benchmark pair) without sufficient justification, and there are minor inconsistencies in how the delay is defined (median vs mode). The presentation is dense and long relative to the core idea.

Most of these need new experiments, not rewriting: comparisons with general delay estimators (GCC-PHAT, coherence) and with simple rescaling of the published setting by an estimated time scale; per-stage tolerance calibration with the stage-matched null; sensitivity of the design constants; regime-switching and nonlinear benchmarks; and more real records per field.
