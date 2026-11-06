# aTDS brief — referee record

Five review rounds by one referee agent playing a senior reviewer from the
network-physiology school that introduced TDS (it signs as "Referee 1"; it does
not speak for any named person). Each report is copied word for word. After
each report there is a short note of what the author changed in response; the
full detail is in the brief itself (`atds_brief.html`).

Final decision (round 5): **APPROVED to show to a professor as an idea.**

---

## Round 1

**1. Summary**
The brief proposes deriving the three TDS parameters from the data. The window is L = 30·B, where B is the Bartlett factor of the slowest signal. The lag range is (L−1)/2, capped at P/2 for rhythmic signals. The tolerance is the widest value whose mean TDS on circularly rotated copies stays ≤ 5 %. The 4-of-5 rule is unchanged. On AR(1) pairs with linear coupling, aTDS keeps detection and delay recovery as the system is slowed by a factor r = 1–16, while the 60/±30/±1 setting fails.

**2. Overall assessment**
The idea is sound and worth pursuing. The useful insights are:
- Measuring the window in independent samples rather than in raw samples.
- Making the lag range scale with the window, which keeps the number of independent lags searched roughly constant at about L/B ≈ 30.
- Showing that a fixed lag range locks onto its edge and produces a confident but wrong "stable" delay (Fig. 4).

The individual ingredients are standard: effective sample size for autocorrelated series (Bartlett), and time-shift surrogates. What is new is calibrating the TDS tolerance against a null. However, the evidence so far shows mainly that parameters should scale with the time scale, which nobody disputes. It does not yet show that data-driven estimation beats simple rescaling, or that it works on physiological data.

**3. Major issues**

1. **The baseline is a straw man.** In Figs 6 and 7 the system is the r=1 system with time stretched by r, and r is known. Any practitioner would rescale to 60r/±30r/±r. At r=1, aTDS itself chose 65/±32/±3. So the large gaps (Fig. 7, r=8: AUC 0.343 vs 0.825) mostly measure unit conversion. *Fix:* add two baselines, (a) the rescaled fixed setting and (b) an oracle grid search on training pairs. Claim a benefit only where aTDS beats (a).
2. **The verdict overstates the results.** The brief says "aTDS beat the fixed parameters in every test, and the fixed parameters won none". The appendix shows otherwise:
   - Fig. 6, r=2: ±1 gives AUC 0.993 vs 0.975 calibrated.
   - Fig. 7, r=1: both 1.0. The sentence "aTDS is better at every slow-down factor" is contradicted in the same paragraph.
   - Fig. 4, d ≤ 30: the ±30 range scores higher (96.5–100 % vs 82.3–96.4 %).
   - None of the results has a confidence interval. With 20 vs 20 pairs, the standard error of an AUC is about 0.06–0.09, so the Fig. 6 gains of +0.04 to +0.10 lie within about one standard error.

   *Fix:* add bootstrap CIs, raise n to at least 100 per class, and rewrite the verdict.
3. **The synthetic tests are too easy and miss the actual regime.** AR(1) has exponential memory, and the tests use Gaussian, linear, stationary coupling with a constant delay. Cross-correlation is close to optimal in exactly that setting. Physiological signals differ:
   - Heart rate and EEG band powers carry long-range correlations (DFA α ≈ 0.5–1.2 depending on sleep stage).
   - For long-range-correlated signals, B truncated at the first zero crossing depends on record length.
   - For oscillatory signals, the first zero crossing falls at P/4, so B badly underestimates the variance inflation.
   - The period cap is admitted to be untested.
   - The repo's own Rulkov benchmark and the data-delta/data-sigma pair (reference TDS 43.6 %) were not used.

   *Fix:* add fGn/ARFIMA tests (α = 0.6–1.0), noisy oscillators, Rulkov, and the EEG benchmark.
4. **The window–resolution trade-off that justifies 30 is never tested.** Fig. 3 shows AUC still rising past the aTDS choice: 1.0 at L=240 for r=2, and 0.999 at L=1920 for r=8 vs 0.755 at 480. With stationary, always-on coupling, longer is always better. *Fix:* use the block-coupled generator already in generators.py (`block_len`) and score time-localised detection against window length. Otherwise present 30 as a convention rather than an optimum. The "±0.18" argument also ignores that τ₀ is the arg-max over about 30 independent lags. The expected noise maximum is about 0.18·√(2 ln 30) ≈ 0.47, which exceeds the c = 0.3 used in Figs 3, 6 and 7. The claim "0.4–0.6 stands clearly above the noise" is therefore optimistic.
5. **Nonstationarity and sleep stages.** At 1 Hz, 30·B for slow signals implies windows of many minutes, comparable to stage episodes. Windows will straddle stage transitions, and a whole-night rotation mixes stages. Parameters chosen per subject or per stage would change the "ruler" exactly where network-physiology conclusions are drawn: differences between stages. *Fix:* estimate B and the tolerance within stage-homogeneous segments, freeze one parameter set per cohort across stages, and report how stage contrasts change.
6. **The circular-shift null is not the original surrogate test.** Cross-subject, stage-matched surrogates (Bashan et al. 2012) keep the generic stage-driven co-modulation in the null. Within-recording rotation removes it, so the null score is lower and the tolerance comes out wider, which is anticonservative on real data. *Fix:* calibrate with both nulls on real data and report the tolerances and false-link rates.
7. **"5 %" is not a significance level.** A mean null score ≤ 5 % bounds the average fraction of false-stable windows, not the probability of a false link. For networks, the link threshold still needs per-pair surrogate percentiles plus a multiple-comparison correction (FDR) across N(N−1)/2 links. One tolerance pooled over up to 50 pairs with different B gives heterogeneous per-pair false-alarm rates. tests.py already computes per-τc chance means and 95th percentiles (testC) but never reports them. *Fix:* report them, and present the 5 % as a design budget.
8. **Mixed-speed systems.** Bartlett's variance for a cross-correlation estimate is 1 + 2Σ r_x(k) r_y(k), not max_i(1 + 2Σ r_i²). For a fast–slow pair this cross term is close to 1, yet aTDS imposes the slow signal's window (for example 1939 vs 39 samples in the Fig. 2 extremes). All tests used τc_x = τc_y. *Fix:* use pairwise B_xy, and test pairs with mismatched speeds.
9. **Some claims are not supported by the code.**
   - Methods says: "Pairs used to set parameters are never among the pairs scored." This is false for Figs 1, 4 and 5: theory.py and tests2.lagB call `adaptive_params` on the scored pair. It is true only for Figs 3, 6 and 7.
   - The tolerance exceeds the delay itself: ±22–39 for d = 5 (Fig. 4), and ±66 at r = 16. "Stable within ±30" cannot then resolve lead vs lag. The direction claim needs a separate test.

**4. Minor issues**
1. Fig. 2 is tautological, since L/B = 30 holds by construction. It is also built from 5 realisations per τc (median B), not "seven test signals". The phrase "delay found in a window is close to noise" is not tested there.
2. Fig. 3: AUC at the aTDS windows themselves (121 and 515) was not computed. Also state that the "fixed 60" line uses a calibrated tolerance (±3 / ±2), not ±1. The summary table's "Fixed TDS" column mislabels this.
3. Fig. 4: "reported as exactly ±30" hides that 6/30 pairs report the wrong direction: τ = +30 at d = 80 (3/5), 120 (1/5) and 160 (2/5). Say so. "97–100 %" should read 96.5–100 %.
4. Fig. 7: AUC 0.362 and 0.343 are below chance, not "guessing". This inversion likely reflects edge-locking in slow uncoupled pairs; report the fixed setting's chance scores (testC) to confirm. "0.82–0.98" should read 0.82–1.0.
5. Fig. 5: the real score is already 86.5 % at tolerance 0, so the Fig. 1 gain (5.2 % to 88.2 %) comes from the window, not the tolerance. State this.
6. Fig. 6: report the coupled mean scores (10.8–12.7 % vs chance 3.4–4.6 %). The effect is modest.
7. theory.py comment "c=0.3" contradicts the code (c = 0.6). Fix it.
8. Claiming the rules "cannot be tuned toward a hoped-for result": y's autocorrelation, and so B, depends on the coupling in general. Soften the wording.

**5. Next steps (prioritised)**
1. Rescaled-fixed and oracle baselines, with bootstrap CIs and n ≥ 100 pairs per class.
2. Real data: the EEG benchmark pair and one sleep-apnea cohort. Compare per-stage networks under aTDS vs the published setting, using cross-subject stage-matched surrogates.
3. Long-range-correlated, oscillatory (testing the period cap) and Rulkov test signals.
4. A block-coupling test to settle the choice of 30 (time resolution vs reliability).
5. Mixed-speed systems with pairwise B_xy.
6. A per-link significance procedure with FDR control, and a false-link rate measured on uncoupled networks.

**6. Recommendation**
A promising idea with honest mechanics, but not yet a result. Rework it against fair baselines and real sleep data before writing a paper.

Referee 1 (network physiology)

**Author's changes after round 1** (experiments `exp3.py`, tasks T1–T9): rescaled-by-hand and oracle baselines; 100 + 100 pairs per class with bootstrap and Wilson intervals; new Table 2 with long-memory (fGn), rhythmic (AR(2)), Rulkov, mixed-speed and EEG tests; Fig. 3b (on/off coupling); separate training pairs for Figs 1, 4, 5; verdict rewritten around "automation, not optimality"; sleep-stage and cross-subject-surrogate tests listed as not done (no staged data available).

---

## Round 2

**1. Status of the original points**

- **M1 (baseline): resolved in substance, overclaimed in wording.** The rescaled and oracle baselines were added. However, a confidence interval that includes 0 does not establish equivalence. At r=4 the interval for aTDS − rescaled is [−0.084, +0.001]. aTDS's AUC point estimate is below the rescaled setting's at r = 2, 4 and 8 (0.985/0.993, 0.912/0.955, 0.792/0.806). Say "not distinguishable", or run an equivalence test with a stated margin.
- **M2 (overstatement): mostly resolved.** Figs 6 and 7 now have n=100 per class with confidence intervals, and the verdict lists losses. Fig. 3 is still 20+20 per class and Fig. 3b uses 10 pairs.
- **M3 (easy tests): partly resolved.** The fGn, AR(2), Rulkov and EEG rows are welcome. Two problems remain:
  - The fGn rows cannot tell the methods apart: with c=0.5 the published setting already scores AUC 1.00 at every H.
  - The Rulkov generator couples the two neurons symmetrically in both directions (both use x[o][i−1−d]), so "true delay positive" is not well defined.
  - Nonstationary signals are still missing.
- **M4 (why 30): partly resolved.** Fig. 3b's "right" criterion is |τ+d| ≤ max(tol,1). Here tol is the calibrated tolerance, which grows with L (13 → 63), so longer windows face a laxer criterion. The test also uses a single block length (2000 samples, about 250 B), so "best at 45–60" holds only for that block length. At c=0.3 every window gives 0.48–0.59, which says nothing.
- **M5 (sleep stages): not resolved.** It is acknowledged.
- **M6 (cross-subject surrogates): not resolved.** It is acknowledged.
- **M7 (5 %): partly resolved.** The wording is now correct, but no per-link test or FDR procedure exists yet. See point N1.
- **M8 (mixed speeds): partly resolved.** Per-pair windows were compared, but no coupled fast–slow pair was tested, which was the case where B_xy matters. Only the chance score of the uncoupled cross pair is reported (2.85 %).
- **M9 (separate training pairs; direction): partly resolved.** Figs 1, 4 and 5 now use separate training pairs. The EEG row still calibrates on the scored pair (T7 calls `adaptive_params([a, b])`); see N2. The direction claim holds for Fig. 4, but not for Rulkov d=32 (N4).
- **Minor 1–6 and 8: resolved.** In minor 4, n=100 removed the inversion. The published setting's chance score is only 0.6–0.9 %, so my edge-locking hypothesis is not supported, and I withdraw it.
- **Minor 7: resolved.** The theory.py comment now reads c=0.6 and is marked superseded.

**2. New problems**

- **N1. The false-alarm level is not controlled.** The Verdict says: "Its value is automation and a controlled false-alarm level". The appendix shows the 5 % budget exceeded on new data in seven cases:
  - r=8: 5.48 %
  - r=16: 6.66 % (95th percentile 23.7 %)
  - fGn H=0.9: 6.18 %
  - EEG: 10.1 %
  - mixed system, system-wide window, fast pair: 6.25 %
  - mixed system, B_xy rule, slow pair: 5.56 %
  - aTDS at r=1: 95th percentile 6.61 %

  The likely cause is that T1 calibrates on a single training pair. Either calibrate on several training pairs or several rotations per pair, or drop "controlled".
- **N2. The Methods claim of separate training pairs is false for the EEG row.** Methods says: "In every figure and table, aTDS sets its parameters on separate training pairs." That is not true for EEG: the parameters were set on the scored pair, and the chance score still reached 10.1 %, double the budget.
- **N3. The EEG aTDS result is uninformative and its parameters are not reported.** Table 2 omits that the tolerance is ±116 and the lag range ±879 (±14 and ±149 in the version capped at 300 s). The reported 100 % over 33 windows, with a ±2-minute tolerance, says almost nothing about delay stability. Report these parameters, and do not count EEG under "better or equal": the brief itself says there is no ground truth beyond detection.
- **N4. Rulkov d=32 is counted as a win but rests on a failed calibration.**
  - The appendix gives aTDS `most_reported = −36`. That is opposite in sign to the 75 % "found", and sits exactly at the lag cap (P/2 = 36).
  - `calibration_failed` is true in both Rulkov rows.

  A flagged calibration failure should not count as "aTDS better".
- **N5. The rhythm explanation does not match the sign.** The table says: "delay wrong by one period (41 − 30 = 11)". Under the code's sign convention, a one-period alias of τ=−30 is τ=+11, which would be reported as −11. The reported +11 has the true sign. That fits the arg-max of |C| picking the anti-phase lobe near d − P/2 ≈ 10 better than a period alias. Please check, then correct both the table and the Limits sentence.
- **N6. Small factual errors.**
  - "using longer windows and tolerance 0" (Verdict table): at r=16 the oracle chose tolerance ±135.
  - Fig. 3b: "0.49–0.59" should read 0.48–0.59.
  - "Tolerance only … true delay found more often from r = 4 on": at r=8 the Wilson intervals overlap ([0.58, 0.76] vs [0.41, 0.61]).
- **N7. The tolerance comparison is partly circular.** With ±1, the coupled score is only about 4 %, so the "found" criterion (at least 3 stable windows) fails for lack of stable windows rather than for a wrong delay. Also report delay accuracy conditional on at least 3 stable windows.
- **N8. Reporting conventions.**
  - The EEG benchmark reproduces the 43.6 % reference only with the causal anchor; the default anchor gives 58.1 %. State which anchor all other results use.
  - T5 reports −mode while T6 reports the raw mode, so the signs in Table 2 are not comparable. Unify them.

**3. Still required before showing the idea**

1. Make the calibration hold its budget out of sample. Use more training pairs and rotations, then verify that mean chance stays ≤ 5 % on new pairs across every row of Figs 6–7 and Table 2.
2. Fix the long-delay failure (H=0.6) with a minimum window or decoupled lag range, and rerun Table 2.
3. Correct N2–N6, and replace "equivalent" with a stated equivalence margin.
4. Redo Fig. 3b with a correctness criterion that does not depend on the tolerance, and with two or three block lengths.
5. Test a coupled fast–slow pair to justify the per-pair B_xy rule.
6. Run a staged multi-subject dataset (any public PSG cohort) with cross-subject, stage-matched surrogates. This remains the decisive test for the target field.

**4. Recommendation**
A substantially improved and now honestly framed brief. It is fit to show as an idea once N1–N5 are corrected, but it is not yet evidence that aTDS works on sleep networks.

Referee 1 (network physiology), second round

**Author's changes after round 2** (experiments `exp4.py`, tasks U1–U12): rule v2.1 — minimum window 61 (was 20) and tolerance from the upper bound (mean + 2 SE) of the chance score with rotations pooled over all training recordings; equivalence test (TOST, ±0.05 AUC); Fig. 3 with 100 + 100 pairs; Fig. 3b with a fixed ±2 criterion and three block lengths; fGn with weaker coupling; one-way Rulkov; nonstationary test; coupled fast–slow pair; EEG parameters from the first half, scored on the second; per-link rotation test with FDR (Table 3); rhythm sign checked (anti-phase lobe); conventions unified.

---

## Round 3

**1. Status of each item**

- **M1 (baseline / equivalence): resolved.** The equivalence test (TOST, margin ±0.05) is applied correctly:
  - r=1: 90 % interval [−0.009, +0.006].
  - r=2: [−0.029, −0.002]. This is equivalent, but the rescaled setting is marginally better, so say so.
  - r=8: aTDS better, 95 % interval [+0.007, +0.142].
  - r=4 and r=16: undecided.

  One cost of v2.1 should be stated: at r=16 the aTDS AUC fell from 0.722 (v2) to 0.645, now below the rescaled setting's 0.688.
- **M2 (sample sizes): resolved.** Fig. 3 now uses 100+100 pairs and Fig. 3b 20 pairs per point, both with intervals.
- **M3 (test signals): partly resolved.**
  - At c=0.2 the fGn rows H=0.6 and 0.75 sit near chance for both methods (AUC 0.59–0.70), so they are again not informative. Only H=0.9 separates the methods.
  - The new nonstationary test is welcome.
  - There is still no test in which the coupling itself switches with the regime, which is the sleep-stage case.
- **M4 (Fig. 3b): resolved, and informative.** With a fixed ±2 criterion, the best L/B grows with episode length: 45 at 1000-sample blocks, 60 at 2000, and still rising at 120 for 4000 (0.946). So 30 suits episodes of about 4 to 8 windows. State that this assumption is built into the constant 30.
- **M7 (per-link test): partly resolved.** The rotation test plus Benjamini–Hochberg exists now, but its result is partly built in (see N-new 3).
- **M8 (mixed speeds): mostly resolved.** For the coupled fast–slow pair, B_xy = 11.3 gives a window of 339 and tolerance ±19, against 939 and ±51 with max-B, at the same AUC. However, both methods reach AUC 1.00, so detection is not tested; a weaker-coupling run is needed. The default rule has not changed.
- **M9 / N2 (separate training data): resolved for EEG** (parameters set on the first half, scored on the second). The Methods wording is now self-contradictory; see N-new 4.
- **N1 (false-alarm budget): resolved.** The mean chance score on new data is ≤ 5 % everywhere except the fast–slow pair with the max-B window (5.24 %). The values I checked range from 0.0 to 4.94 %.
- **N3 (EEG parameters): resolved.**
- **N4 (Rulkov): partly resolved.** Calibration no longer fails, but see N-new 2.
- **N5 (rhythm sign): resolved in substance.** The wording is wrong; see N-new 5.
- **N6 (factual errors): resolved.**
- **N7 (conditional delay accuracy): resolved.** The new numbers show the cost of calibration. Among pairs with at least 3 stable windows, ±1 is right 100 % of the time, against 90/91/79 % calibrated at r=4/8/16. The gain in "delay found" comes from more pairs having stable windows, not from more accurate delays. The text says this correctly.
- **N8 (conventions): resolved.**
- **M5 and M6 (sleep stages; cross-subject surrogates): not resolved.** This is acknowledged.

**2. New problems**

1. **The minimum window of 61 reintroduces a sample-count constant.** The brief argues that sample-count constants "mean something different for every system and every sampling rate". Yet 61 exists only to reproduce the published ±30, and it does nothing at H=0.6: delay found 23 % vs 28 % for the published setting. *Fix:* take the minimum lag from the longest plausible delay in seconds, given by the user, and present it as a required input.
2. **The Rulkov "win" and "miss" are the same result.** aTDS reports 19 for d=16 and 35 for d=32, the same +3 offset each time, which is probably the map's response latency. d=32 counts as "found" only because 3 ≤ 10 % of 32, while d=16 misses because 3 > 1.6. U7 also uses a mode-based hit rule, not the Methods definition ("at least 3 stable windows and their median τ₀"). Report the offset, and do not count d=32 as "aTDS better" and d=16 as "wrong".
3. **"0 false links" in Table 3 is partly built in.** With 99 rotations the smallest possible p is 0.01, and Benjamini–Hochberg over 15 links can then reject only if at least 3 p-values sit at that floor. Any network with fewer than 3 true links has zero power by construction. Two consequences:
   - Report the raw per-link false-positive rate instead: 5/150 = 3.3 % at p ≤ 0.05, which shows the per-link test is calibrated.
   - Rerun with 999 rotations. It is listed as a next step, but it is cheap and should be done now.

   The links are also disjoint, so indirect links (a common driver x→y, x→z producing a spurious y–z link) are untested. That is the central false-link mechanism in physiological networks.
4. **The Methods sentence on training data contradicts itself.** It says: "separate training data that are not scored: … the network's own signals (no labels) in Table 3". In Table 3 the network's own signals are the scored data. Reword it as "label-free calibration on the scored network".
5. **Rhythm row: the 93 % is a share of pairs, not of windows.** The table says "93 % of stable windows have negative C". U6 computes the sign only at each pair's first stable window, so the correct statement is "93 % of pairs (first stable window)".
6. **The algorithm section has a stale constant.** "Constants: … window between 20 samples and 10 % of T" should read 61.
7. **The SE in the tolerance rule ignores between-recording variance.** Rotations from 5 recordings are pooled as if they were independent, which explains the 5.24 % breach. Use an SE clustered by recording.
8. **Smaller points.**
   - Fig. 6: at r=1 the calibrated tolerance is significantly worse than ±1, with 95 % interval [−0.012, −0.001]. The difference is tiny, but "AUC gain not significant" omits it.
   - Table 2, nonstationary row: per-segment parameters use the known change point (oracle segmentation). State this.

**3. Still required before showing the idea**

1. Fix N-new 2, 4, 5 and 6 (text-only edits), and rerun Table 3 with 999 rotations, reporting the raw per-link false-positive rate (N-new 3).
2. Add a network with a common driver, and report its false-link rate.
3. Add a regime-switching test where the coupling is on in one regime only, scored per regime.
4. Replace the fixed minimum of 61 with a user-stated maximum delay in seconds.
5. Test a staged multi-subject cohort with cross-subject, stage-matched surrogates (still outstanding from M5 and M6).

**4. Recommendation**
Ready to show as an idea once items 1–2 of section 3 are done. It is not yet evidence for sleep-stage networks.

Referee 1 (network physiology), third round

**Author's changes after round 3** (experiments `exp5.py` V1–V2, `exp6.py` W1): Rulkov rerun with the median rule (constant +3 offset reported); Table 3 rerun with 999 rotations, raw per-link false-positive rate, and common-driver networks; text fixes N-new 2, 4, 5, 6 and the smaller points. Design decision by the author on item 4: the only mandatory input stays the two signals; the minimum window stays 61; the longest plausible delay d_max is an optional input (window ≥ max(61, 2·d_max + 1)), tested on a fast signal with delay 50.

---

## Round 4

**1. Items 1 and 2**

**Item 1 (text fixes and Table 3): resolved.**
- **Rulkov.** Under the Methods median rule, aTDS reports 19 for d=16 and 35 for d=32 (offset +3.0 at both), with AUC 1.00 vs 0.73 and 0.79 for the published setting. Both delays are now in one row, and the brief says that the 10 % rule splitting them into a hit and a miss is a criterion artefact.
- **Table 3 at 999 rotations.** I recomputed the figures from the appendix:
  - True links found: 27/30 (disjoint) and 26/30 (driver).
  - Mean false-discovery proportion: 8.3 % and 2.5 %.
  - False links in the 20 uncoupled networks: 0.
  - Raw rate of p ≤ 0.05 on uncoupled links: 4/150 = 2.7 % in uncoupled networks; 7/120 = 5.8 % and 6/110 = 5.5 % in coupled ones.

  This matches the stated "2.7–5.8 %". The per-link test is calibrated.
- **Other text fixes.** The Methods wording for Table 3, the rhythm-sign sentence ("in 93 % of pairs"), the constants line (61), and the r=1, r=2 and r=16 sentences all match the data: −0.005 [−0.012, −0.001]; 0.993 vs 0.979; 0.645 / 0.722 / 0.688.

**Item 2 (common driver): partly resolved, adequately for an idea brief.**
- The design and the bookkeeping are correct: the indirect 1–2 link is tracked separately, and its raw p ≤ 0.05 in 1/10 networks (p = 0.026 in that one).
- However, with c=0.4 the indirect correlation is only 0.16. The test therefore does not stress the method, and "not reported" must not be read as robustness to indirect links.
- The text says this correctly: "with stronger drivers indirect links are expected to appear". Keep that sentence.

**2. New problems**

1. **A wrong number in Limits.** The text says: "in recordings shorter than about 1 200 samples it leaves fewer than about 20 windows".
   - With L = 61 and step 30, the number of windows is ⌊(T − 61)/30⌋ + 1, which is 38 at T = 1200. It falls to about 20 only near T ≈ 630.
   - Whenever the 10 % cap binds, there are always 19 windows.
   - In the code (`hi = max(min_window, 0.1·T)`), the 61 floor overrides the 10 % cap only when T < 610.

   *Fix:* change the sentence to "shorter than about 600 samples (where 61 exceeds 10 % of the record)", and say in the pseudocode that L_min takes precedence over the cap.
2. **The optional d_max was tested only with a correct value.** It was tested at d_max = 60 for a true delay of 50. The cost of an overstated d_max (a longer window, so coarser timing) and the effect of an understated one are untested. One sentence in Limits is enough.
3. **The default fails silently.** For the fast signal with delay 50, default aTDS gives a coupled score of 4.81 % against a chance score of 4.37 %, with no flag. State plainly that nothing in the output warns the user, so d_max is the only safeguard for long delays in fast signals.

**3. Ready to show?**

Yes, after fixing the Limits sentence (point 1) and adding one sentence on point 2. The brief now:
- states its baselines fairly (TOST against the rescaled setting, an oracle upper bound);
- reports its losses (r=16, H=0.6/0.75, rhythm beyond P/2, the Rulkov offset);
- separates training data from scored data;
- includes a calibrated per-link test with FDR control.

These can stay as stated future work, provided they appear in the "Next steps" and "Not tested here" sections:
- regime-switching coupling;
- a staged multi-subject cohort with cross-subject, stage-matched surrogates;
- a standard error clustered by recording;
- a weaker-coupling fast–slow pair and stronger fGn coupling;
- stronger common drivers and chains.

The author should state up front that the central claim for sleep networks (that aTDS changes or sharpens stage contrasts) is untested.

**4. Recommendation**
Ready to show to a professor as a well-tested idea once the Limits sentence is corrected; the claims for sleep-stage networks remain future work.

Referee 1 (network physiology), fourth round

**Author's changes after round 4** (text only): the Limits sentence corrected (about 600 samples, fewer than 19 windows) and the pseudocode note that L_min wins over the 10 % cap; sentences on d_max tested only with a correct value and on the silent failure; future-work items listed; the untested sleep-network claim stated in the header; new section "How to read this brief" defining every term, including the exact computation of the oracle.

---

## Round 5

**1. Round-4 points**
All three are resolved.
- **Limits sentence.** It now gives about 600 samples and fewer than 19 windows, which matches ⌊(T − 61)/30⌋ + 1. The pseudocode now notes that L_min wins over the 10 % cap.
- **d_max.** Limits now states that d_max was tested only with a correct value, and what an overstated or understated d_max does.
- **Silent failure.** The brief now states it (4.8 % against a chance score of 4.4 %, which reads as "no coupling") and names d_max as the only safeguard.

The future-work list and the up-front statement that effects on sleep-stage contrasts are untested are both adequate.

**2. Remaining errors in the new "How to read this brief" section**
These are all text-level.

1. **TDS score.** "High score = stable delay = coupling." This contradicts the brief's own Fig. 4, where scores of 72–96 % came with a wrong delay. Replace with: "A score well above the chance score of the same data indicates a stable delay. The score alone says nothing about whether the delay is right, and aTDS scores are not comparable with published ones."
2. **Per-link p-value.** "the share of its 999 rotated copies that score at least as high." The code computes p = (1 + count)/1000. Use that formula.
3. **Benjamini–Hochberg.** "keeps the expected share of false links … at 5 %." Write "at most 5 %, assuming the links are independent or positively dependent".
4. **Intervals.** "If two methods' intervals do not overlap, the difference is real." Write "significant at roughly the 5 % level". Add that, for methods scored on the same pairs, the interval of the paired difference is the proper test, and overlapping intervals do not mean "no difference".
5. **Oracle.** "The best parameters that exist." It is the best of 91 grid settings, chosen on 40 training pairs. Call it "a practical upper bound", since a setting outside the grid could do better (at r=16 the oracle reaches only 0.84).
6. **Training vs scored pairs.** This entry omits Table 3. Add: "In Table 3, calibration is label-free on the scored network itself."
7. **aTDS.** "only the two signals go in." Add "(optionally d_max; for a network, all of its signals)".
8. **Balanced accuracy.** A window counts as correctly "coupled" only if it is stable and its τ₀ is within ±2 samples of the true delay. Say so.
9. **EEG pair.** Add that the 43.6 % reference is reproduced only with the causal stability anchor (43.5 %); the default anchor gives 58.1 %.
10. **Published setting.** Note "lag ±30 (±29 effective)", as in the Methods section.

**3. Decision**
**APPROVED** to show to a professor as an idea. Edit 1 should be made before sending, because as written it misstates what the TDS score means. Edits 2–10 are one-line clarifications. None of them needs another review round.

Referee 1 (network physiology), fifth round

**Author's changes after round 5:** all ten edits made as written. This is the version saved in `atds_brief.html` and published.
