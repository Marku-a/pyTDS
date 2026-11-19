# Experimental log: aTDS

All numbers below come from fixed-seed runs; raw JSON in the pyTDS repository, atds_brief/results/. Do not report any number that is not in this log.


# PART 1: SIMULATED DATA (synthetic signals with known coupling and delay)

## How to read this brief

Everything in the tests below is measured the same way. This section explains each term once, so no figure needs guessing.

    Known answerMost tests use synthetic signals that we generate ourselves. We decide whether a pair is coupled, and with what delay, so we can check whether a method gets it right. One test uses a real EEG pair, where only detection can be checked.
    Coupled / uncoupled pairA coupled pair: y copies part of x with a delay d, y(t) = c · x(t − d) + noise, with coupling strength c between 0 and 1. An uncoupled pair is generated the same way with c = 0, so y has no link to x. Each test scores many coupled and many uncoupled pairs.
    TDS scoreThe percentage of windows labelled stable (4 of 5 neighbouring delays agree within the tolerance). A score well above the chance score of the same data indicates a stable delay, and so coupling. The score alone says nothing about whether the delay is right (Figure 4), and aTDS scores are not comparable with published TDS scores.
    Chance scoreThe TDS score a pair gets when the coupling is destroyed by rotating y in time (its end wraps to its start). It tells how high a score luck alone gives. "Chance" in the tables is the mean over many uncoupled pairs or rotated copies; "95th pct" is the value only 5 % of them exceed.
    AUCHow well the TDS score separates coupled from uncoupled pairs: the probability that a randomly chosen coupled pair scores higher than a randomly chosen uncoupled one. 1 = always (perfect), 0.5 = no better than a coin flip. Example: AUC 0.9 means a coupled pair outscores an uncoupled one 9 times out of 10.
    Delay foundThe share of coupled pairs for which the method reports the right delay: at least 3 stable windows, and the median τ₀ of those windows within 10 % of the true delay (at least ±1 sample).
    95 % intervalThe range the true value very likely lies in, given the finite number of pairs. For AUC: 1 000 bootstrap resamples (draw the pairs again with replacement, recompute, take the middle 95 %). For shares such as "delay found": the Wilson interval. If two methods' intervals do not overlap, the difference is significant at roughly the 5 % level. For methods scored on the same pairs the proper test is the interval of the paired difference (next entry); overlapping intervals do not mean "no difference".
    "Better", "equivalent", "undecided"For two methods scored on the same pairs, we compute the interval of the difference in AUC. "Better": the 95 % interval of the difference is entirely above 0. "Equivalent": the 90 % interval lies inside ±0.05, i.e. any difference is smaller than 0.05 AUC (an equivalence test, TOST). "Undecided": neither.
    Training vs scored pairsaTDS chooses its parameters on separate training pairs (or, for the EEG, the first half of the record) and is then scored on new pairs it never saw, so the parameters cannot fit the answer by accident. In Table 3, calibration is label-free on the scored network itself (rotations only), as in a real analysis.
    Slow-down factor rOne test system made r times slower: memory 2r samples and delay 5r samples. This is what happens when the same physical system is recorded r times faster. At r = 1 the published values fit; the larger r, the worse they fit.

## The four settings compared

    PublishedThe fixed TDS setting of Bashan et al. (2012): window 60, step 30, lag ±30 (±29 effective in pyTDS's circular cross-correlation), tolerance ±1, the same for every system.
    Rescaled by handThe published setting multiplied by the slow-down factor r: window 60r, step 30r, lag ±30r, tolerance ±r (at r = 4: 240 / ±120 / ±4). It is what a careful analyst would choose if told how much slower the system is. In real data nobody knows r, so this is a fair but generous baseline.
    aTDSThe method of this brief: only the two signals go in (optionally the longest plausible delay dmax; for a network, all of its signals), and it sets window, lag range and tolerance itself.
    OracleThe best of 91 grid settings, chosen by checking the answers on labelled training pairs. It is a practical upper bound, not a method anyone can use, because it needs to know which pairs are coupled; a setting outside the grid could do even better (at r = 16 the oracle reaches only AUC 0.84). How it is computed, for each slow-down factor r:

Generate 20 coupled and 20 uncoupled training pairs, separate from the pairs that are scored later.

Try every combination of 7 windows (30, 60, 120, 240, 480, 960, 1 920 samples; lag range (L − 1)/2 each) and 13 tolerances (0, 1, 2, 3, 5, 8, 12, 18, 27, 40, 60, 90, 135; tolerances larger than the lag range are skipped): up to 91 settings.

For each setting, score the 40 training pairs and compute the AUC.

Keep the setting with the highest training AUC (the first one if several tie).

Score the same 100 + 100 test pairs as the other three settings with it.

      How it is shown: grey line in Figure 7; the setting it chose at each r is in the Figure 7 data table. If aTDS is close to the oracle, little is left to gain by tuning; if it is far below, better parameters exist.

## Signals used

    AR(1)Smooth random noise whose memory fades exponentially; τc is how many samples it takes to forget its past. The main test signal.
    fGn (long memory)Fractional Gaussian noise: memory that fades slowly, like heart rate or EEG power. Its Hurst exponent H (0.5 = no memory, closer to 1 = long memory) equals the DFA exponent α used in physiology.
    AR(2) (rhythm)A noisy oscillator with a period of about 40 samples, like breathing.
    Rulkov neuronsA standard model of two bursting neurons; neuron 0 drives neuron 1 with a delay.
    EEG pairReal sleep EEG δ- and σ-band power of one night at 1 Hz, from the pyTDS repository. The original MATLAB result on it, 43.6 %, is reproduced by pyTDS only with the causal stability anchor (43.5 %); the default anchor, used everywhere else in this brief, gives 58.1 %.

## Other terms

    B, independent samplesB (Bartlett factor) is how many consecutive samples count as one independent sample. A window of L samples holds L / B independent samples. Bxy is the same for a pair of signals.
    Balanced accuracyUsed in Figure 3b: the average of two hit rates. While coupling is on, a window counts as correct only if it is labelled stable and its τ₀ is within ±2 samples of the true delay; while coupling is off, a window counts as correct if it is not labelled stable. 0.5 = guessing, 1 = perfect.
    Per-link p-value, FDRUsed in Table 3. A link's p-value is p = (1 + number of its 999 rotated copies that score at least as high as the real pair) / 1 000; the smallest possible value is 0.001. When many links are tested at once, some pass by luck; the Benjamini–Hochberg procedure keeps the expected share of false links among the reported ones (the false-discovery rate) at most 5 %, assuming the links are independent or positively dependent.
    v2.1The current version of the aTDS rules, described in the algorithm section (minimum window 61; tolerance from the upper bound of the chance score).

## Test 1 · Window

Figure 2 · A fixed window holds less information as the signal gets slower

The x-axis says how slow a signal is. Its correlation time τc is the number of samples it takes for the signal to "forget" its past (for its autocorrelation to fall to 1/e ≈ 0.37). τc = 1 is a fast, noisy signal; τc = 64 is a slow, smooth one, for example the same signal recorded 64 times faster. For each τc, B is the median over 5 test signals of 30 000 samples. The y-axis gives the number of independent samples a window holds, L / B.

What it shows: with a fixed 60-sample window, slower systems are measured with less information, down to 0.93 independent samples at τc = 64. The orange line is flat at 30 by construction, since aTDS sets L = 30 · B; the figure illustrates the rule, it does not test it. Figures 3 and 7 test what the window does to detection.

Figure 3 · Detection against window length, with the aTDS choice marked (100 pairs per class)

Two systems with known weak coupling (strength 0.3): a faster one (τc = 4, delay 10) and one 4× slower (τc = 16, delay 40). For each window length, 100 coupled and 100 uncoupled pairs were scored. The lag range is (L − 1)/2 and the tolerance is calibrated for every window length, including the 60-sample point, so only the window varies. AUC measures how well the score separates coupled from uncoupled pairs (1 = perfect, 0.5 = guessing); whiskers are bootstrap 95 % intervals. Vertical lines: the fixed window (dashed, 60) and the window aTDS chose for each system separately (121 for the faster, 515 for the slower system), in that system's colour; the orange dots are the AUC measured at exactly those windows.

What it shows: the window aTDS picks moves with the time scale (4× longer for the 4× slower system). At its own windows aTDS reaches AUC 0.97 [0.95, 0.99] for the faster and 0.79 [0.73, 0.85] for the slower system, against 0.83 and 0.51 at 60 samples. Detection keeps rising past the aTDS choice: with coupling that is always on, a longer window is always better (0.99 at 1 920 samples for the slower system). The aTDS choice is a compromise with time resolution, not the best window for detection.

Figure 3b · When coupling switches on and off: is 30 independent samples the right window?

Coupling (strength 0.6, τc = 8, delay 20) switched on and off in blocks of 1 000, 2 000 or 4 000 samples (about 125, 250 and 500 independent samples per block), 40 000 samples, 20 pairs per point. A window is scored as right when, during coupling, it is labelled stable and its τ₀ is within ±2 samples of the true delay (the same strict criterion for every window length), or, when coupling is off, it is not labelled stable. The y-axis is the balanced accuracy of these labels (0.5 = guessing). The x-axis is the number of independent samples per window, L / B. The vertical line marks the aTDS choice of 30.

What it shows: the trade-off behind the window rule is real, and the best window depends on how long the coupling episodes last. For the shortest blocks the best is 45 independent samples (0.87; 30 gives 0.85), for 2 000-sample blocks 60 (0.92; 30 gives 0.88), and for the longest blocks a plateau of 0.94–0.95 from 45 on (30 gives 0.91). 30 is within 0.015–0.05 of the best in all three cases and on the short side. It is a convention that favours time resolution, not an optimum.

## Test 2 · Lag range

Figure 4 · Only the lag range changes: ±30 vs the aTDS lag range

Slow coupled signals (τc = 16, coupling strength 0.6) with true delays from 5 to 160 samples, 10 pairs per delay. For each delay, aTDS set its parameters on 5 separate training pairs (window 454–514, tolerance ±23 to ±29). Both variants use that window and tolerance; only the lag range differs: the published ±30, or the aTDS range (L − 1)/2, here ±226 to ±256. Left: the delay each variant reports (the most frequent τ₀ among its stable windows; positive = x leads). The diagonal is a perfect answer. Right: the TDS score of the same runs. Filled dots reported the right delay; hollow dots a wrong one.

What it shows: with the aTDS lag range, the reported delay is exact for all 90 pairs, from 5 to 160 samples, sign included. With ±30, all 50 pairs with delays beyond 30 report a wrong delay at the edge of the search (±30, once 18), 16 of them with the wrong sign, and their TDS scores stay high (72–96 %). With a fixed lag range, a high TDS score can come with a wrong delay: the TDS score measures how stable the delay is, not whether it is right.

Why aTDS scores fall at long delays. The longer the delay, the smaller the part of each window in which x and the delayed y overlap, so the correlation peak gets weaker and fewer windows pass the 4-of-5 rule (scores 27–62 % at delay 160). The delay reported from those windows is still right. Within ±30 both variants are right, and the narrow range scores somewhat higher (94–100 % vs 82–99 %) because it has fewer lags to pick wrongly by chance. The tolerance (±23 to ±29) is wider than the shortest delays; it only decides which windows count as stable. Direction comes from the reported delay itself, which was right in every aTDS pair.

## Test 3 · Tolerance

Figure 5 · How the tolerance is chosen

The pair from Figure 1, with the aTDS window. Left: the chance score C(tol), the mean TDS of 20 circularly rotated copies of this pair (blue), and its upper bound C + 2 SE (dashed blue), at each tolerance. Right: the real, unshifted pair at the same tolerances. aTDS chose ±13 from rotations of its 5 separate training pairs; on this pair the upper bound crosses 5 % between ±12 and ±13.

What it shows: the tolerance comes from the coupling-free copies alone. Here the real score barely depends on it (90.0 % at tolerance 0, 91.3 % at ±13), so the jump from 5.2 % to 91.3 % in Figure 1 comes from the window and lag range, not from the tolerance. The tolerance's job is to cap false stability while allowing the delay to wobble.

Figure 6 · Only the tolerance changes: ±1 vs calibrated, 100 pairs per class, 95 % intervals

The slow-down factor r. One test system is made slower step by step: at r = 1 its memory is τc = 2 samples and its delay 5 samples, a case the fixed parameters fit well. At factor r both are multiplied by r (memory 2r, delay 5r), so r = 16 means memory 32 and delay 80 samples. This is what happens when the same physical system is recorded r times faster (for example at 16 Hz instead of 1 Hz). The coupling is weak (strength 0.3). Both variants use the same aTDS window and lag range, set on 5 separate training pairs; only the tolerance differs. 100 coupled and 100 uncoupled pairs per r. Left: AUC with bootstrap 95 % intervals. Right: share of coupled pairs where the stable windows report the true delay (within 10 %), with Wilson 95 % intervals.

What it shows: calibrating the tolerance helps slow systems find the right delay: 85 % vs 66 % of pairs at r = 8 (intervals do not overlap), 88 % vs 77 % at r = 4 and 42 % vs 26 % at r = 16 (intervals overlap). The AUC gain (+0.04 at r = 4 and 8) is not significant at 95 %; at r = 1 the calibrated tolerance is very slightly worse (−0.005, interval −0.012 to −0.001). The gain comes from having enough stable windows: with ±1, only 77, 66 and 26 of 100 coupled pairs at r = 4, 8, 16 have 3 or more stable windows, and those few are almost always right (100 %); calibrated, 98, 93 and 53 pairs have them, of which 90 %, 91 % and 79 % are right. The mean chance score on new pairs is 3.3–4.4 % at every r, inside the 5 % budget.

## All three rules together, against fair baselines

Figure 7 · aTDS against the published setting, the setting rescaled by hand, and an oracle

The same slow-down test as Figure 6 (r = 1 to 16, 100 coupled and 100 uncoupled pairs per r). Four settings:

Published: 60 / ±30 / ±1, the same at every r.

Rescaled by hand: the published setting multiplied by the slow-down factor r: window 60r, step 30r, lag ±30r, tolerance ±r (at r = 4: 240 / ±120 / ±4). Nothing is tuned; the formula is applied once per r. It is what a careful analyst would choose if told that the system is exactly r times slower than the 1-Hz case the published values were made for. In real data nobody knows r, so this baseline is fair but generous. aTDS is not told r.

aTDS: all three rules (v2.1), set on 5 separate training pairs.

Oracle: the best window and tolerance from a grid of 7 windows × 13 tolerances, chosen by AUC on 20 + 20 separate, labelled training pairs. It uses knowledge no real analysis has: an upper bound.

What it shows: aTDS is far better than the published setting once the system is slower than that setting assumes (AUC 0.92 vs 0.59 at r = 4; 95 % interval of the difference +0.24 to +0.42). Against the setting rescaled by hand, with an equivalence margin of ±0.05 AUC (90 % interval of the difference inside ±0.05): equivalent at r = 1 and 2 (at r = 2 the rescaled setting is marginally higher, 0.99 vs 0.98), better at r = 8 (0.88 vs 0.81; 95 % interval +0.007 to +0.142; true delay found in 85 % vs 60 % of pairs), and undecided at r = 4 and 16 (90 % intervals −0.064 to +0.003 and −0.105 to +0.018). The oracle is better than aTDS at every r ≥ 2 (0.99 vs 0.88 at r = 8), mostly by choosing much longer windows. The v2.1 rule (stricter tolerance) cost some detection at r = 16: AUC 0.65, against 0.72 with the earlier rule and 0.69 for the rescaled setting. The published setting is close to chance from r = 4 on (AUC 0.51–0.59), with chance scores of 0.6–0.9 %.

## Test 4 · Other kinds of signals, and one real recording

The tests above use one kind of signal: smooth noise with short memory (AR(1)) and a constant delay. Physiological signals often have long memory (heart rate, EEG band power), rhythms (breathing, heartbeat), bursts (neurons), mixed speeds, or change their dynamics over time (sleep stages). Table 2 runs aTDS against the published setting on each of these, and on the real sleep-EEG benchmark pair in the repository. In every row aTDS set its parameters on data that were not scored.

       | 
 | Signals (known answer) | Published 60 / ±30 / ±1 | aTDS v2.1 | Reading

 | Long memory, H = 0.6 (fGn), delay 20, weak coupling 0.2 | AUC 0.60 [0.47, 0.72], delay found 28 % | AUC 0.59 [0.45, 0.71], 23 % (window 61, lag ±30) | Equal. The old minimum window of 20 gave window 32, lag ±15 and AUC 0.47; the new minimum of 61 fixes that failure.

 | Long memory, H = 0.75, delay 20, coupling 0.2 | AUC 0.70 [0.59, 0.81], 33 % | AUC 0.65 [0.53, 0.77], 23 % (window 73) | Not distinguishable (intervals overlap).

 | Long memory, H = 0.9, delay 20, coupling 0.2 | AUC 0.46 [0.33, 0.58], 15 % | AUC 0.82 [0.73, 0.91], 85 % (window 614, tol ±34) | aTDS better. B depends on record length here (13 → 23 → 84 for 5 000 → 20 000 → 80 000 samples).

 | Rhythm, period ≈ 40, delay 8 (below half a period) | AUC 0.97, delay found 43 % | AUC 1.00, 100 % (period found 41, lag capped at ±20) | aTDS better.

 | Rhythm, period ≈ 40, delay 30 (beyond half a period) | AUC 0.46, 0 % | AUC 0.995, delay 0 %: reports 11 | Coupling detected, delay wrong. 30 is outside the capped lag range (±20); the largest |C| inside it is the anti-phase lobe (in 93 % of pairs the first stable window has negative C), about half a period from the true delay, with the right sign.

 | Rulkov neurons, one-way coupling, delays 16 and 32 | AUC 0.73 and 0.79; median delay reported 0 and 15.5 | AUC 1.00 and 1.00; median delay reported 19 and 35 | aTDS detects better and reports the true delay plus a constant 3 steps at both delays, probably the map's response latency. (By the 10 % rule this counts as found at 32 and missed at 16; it is the same result.) Calibration did not fail with one-way coupling.

 | Fast signal (τc = 1) with a long delay, 50, coupling 0.5, 100 + 100 pairs | AUC 0.56, delay found 0 % | default (window 61, lag ±30): AUC 0.56, 0 %; with the optional dmax = 60 (window 121, lag ±60): AUC 1.00, 100 % | A delay beyond ±30 in a fast signal needs the optional input. At delay 20 all three settings are perfect.

 | Changing dynamics: memory τc = 4 in the first half, 16 in the second, delay 10, coupling 0.4 | first half AUC 0.98; second half 0.71, 45 % | one window for the whole record (252): 1.00 and 0.96, 90 %; a window per half, using the known change point (124 / 484): 1.00 and 0.98, 100 % | aTDS better; setting it per segment is slightly better still.

 | Mixed system: slow pair (τc 32, delay 40) + fast pair (τc 2, delay 5), coupling 0.4 | slow AUC 0.44, fast 1.00 | one window (972): slow 0.99, fast 1.00, but fast pair tol ±56; per-pair windows: fast 64, tol ±3 | A system-wide window costs fast pairs delay precision.

 | Coupled fast–slow pair: slow x (τc 32) drives fast y (τc 2), delay 10, coupling 0.5 | AUC 0.83, delay found 57 % | window from max B (939, ±51): AUC 1.00, 100 %, chance 5.2 %; window from the pair's joint Bxy (339, ±19): AUC 1.00, 100 %, chance 3.9 % | The joint Bxy = 1 + 2Σ rxry (here 11.3, between Bx 31.3 and By 5.0) gives the same detection with a window 2.8× shorter and a tolerance 2.7× tighter.

 | Real: sleep EEG δ and σ power, 1 Hz; parameters from the first half (14 970 s), scored on the second half | 50.8 % vs chance 1.2 % (95th pct 2.6); 37.3 % with the causal anchor | window 652 (11 min), lag ±325, tol ±22, 44 windows: 95.5 % vs chance 2.1 % (95th pct 11.4); capped at 300 s (lag ±149, tol ±9): 82.7 % vs 2.3 % | Both detect the known strong coupling (p ≤ 0.01, the smallest possible with 100 rotations), delay 0 in both. No ground truth beyond detection, so this row is not counted for or against aTDS.

Long-memory signals: fractional Gaussian noise with Hurst exponent H (DFA α = H), 40 + 40 pairs. Rhythms: noisy AR(2) oscillators, 30 + 30 pairs. Rulkov: two Rulkov map neurons, neuron 1 driven by neuron 0 with the given delay (coupling 0.12 vs 0), 20 + 20 pairs. Changing dynamics, mixed system and fast–slow pair: 40 + 40. Intervals in brackets are bootstrap 95 %. The EEG reference value 43.6 % from the original MATLAB code is reproduced by pyTDS only with the causal stability anchor (43.5 % on the whole record); every other result on this page uses the default anchor.

Table 3 · Links in a network: a per-link test with false-discovery-rate control

Networks of 6 signals (τc = 8, 10 000 samples), 15 possible links. Two designs, 10 coupled and 10 uncoupled networks each: disjoint (3 separate true links: 0→1, 2→3, 4→5, delays 10, 20, 15) and common driver (signal 0 drives 1 with delay 10 and 2 with delay 20, plus 3→4 with delay 15). A common driver creates an indirect link between 1 and 2 (delay 10) that is not a direct coupling. Coupling 0.4. aTDS sets one parameter set from the network's own signals by rotation only (no labels). Each link gets a p-value from its score against 999 rotated copies of the same pair; the 15 links are tested together with the Benjamini–Hochberg procedure at q = 0.05.

What it shows: the per-link test is calibrated: on links that are not coupled, 2.7–5.8 % reach p ≤ 0.05, close to the nominal 5 %. After FDR control there is no false link in any of the 20 uncoupled networks. In the coupled networks aTDS finds 27 of 30 true links (disjoint) and 26 of 30 (common driver); the false-discovery proportion averages 8.3 % and 2.5 %, against the 5 % target (10 networks each, so these are noisy). The indirect link of the common driver was not reported in any network (raw p ≤ 0.05 in 1 of 10): at coupling 0.4 the indirect correlation is only about 0.4 × 0.4 = 0.16. Pairwise TDS cannot tell an indirect link from a direct one, so with stronger drivers indirect links are expected to appear.

## Verdict

aTDS does automatically what a careful analyst would do by hand: it scales the TDS parameters to the system, without being told the system's time scale.

Against the published setting it is much better whenever the system is slower than that setting assumes. Against the setting rescaled by hand it is equivalent or better where the tests can decide (r = 1, 2, 8), and undecided elsewhere. It does not find the best possible parameters: an oracle that looks at labelled data does better. Its value is automation with a false-alarm budget that held on new data in the main test.

       | 
 | Comparison | Result

 | aTDS vs published 60 / ±30 / ±1 (Fig. 7) | aTDS better from r = 2 on (AUC difference, 95 %: +0.05 to +0.14 at r = 2, +0.24 to +0.42 at r = 4, +0.03 to +0.25 at r = 16); equal at r = 1.

 | aTDS vs the setting rescaled by hand (Fig. 7) | Equivalent (±0.05 AUC) at r = 1 and 2; aTDS better at r = 8; undecided at r = 4 and 16. aTDS is not told r.

 | aTDS vs oracle grid (Fig. 7) | Oracle better at r ≥ 2, mostly through longer windows.

 | Lag range only (Fig. 4) | ±30 gives a wrong delay for all 50 pairs with delays beyond 30 (16 with the wrong sign) at scores of 72–96 %; aTDS exact in all 90 pairs.

 | Tolerance only (Fig. 6) | Delay found more often at r = 8 (85 vs 66 %, significant), also higher at r = 4 and 16 (not significant); AUC gain not significant; chance 3.3–4.4 % on new pairs.

 | Window and on/off coupling (Fig. 3b) | 30 independent samples is within 0.015–0.05 of the best balanced accuracy for all three episode lengths, on the short side.

 | Other signals (Table 2) | Better on long memory H = 0.9, rhythms with delay below half a period, changing dynamics, and Rulkov detection (delay reported with a constant +3 offset). Equal or not distinguishable on H = 0.6 and 0.75. Wrong delay (coupling still detected) on rhythms with delay beyond half a period. A fast signal with delay 50 needs the optional dmax; with it, AUC 1.00.

 | Network links (Table 3) | Per-link rotation test (999 rotations) + FDR: 0 false links in 20 uncoupled networks; 87–90 % of true links found; false-discovery proportion 2.5–8.3 %; the common driver's indirect link not reported.

## Limits

Not optimal, only scaled. With coupling that is always on, longer windows detect better than aTDS (Figures 3 and 7, oracle). 30 independent samples trades some detection for time resolution; for longer coupling episodes, 45–60 did slightly better (Figure 3b).

The minimum window sets a minimum lag range. In v2.1 the window is never shorter than 61 samples, so the lag range is never narrower than ±30, the published range. Delays longer than half the window still cannot be found: a fast signal with delay 50 was missed by default (AUC 0.56). The optional longest plausible delay dmax fixes this (AUC 1.00 with dmax = 60). Nothing in the output warns the user when a delay is out of reach: for the fast signal with delay 50 the default score was 4.8 % against a chance score of 4.4 %, which reads as "no coupling". dmax is the only safeguard for long delays in fast signals. dmax was tested only with a correct value (60 for a true delay of 50): an overstated dmax makes the window longer than needed (coarser timing), an understated one leaves the delay out of reach. The minimum of 61 also has a cost: in recordings shorter than about 600 samples, 61 is more than 10 % of the record and takes precedence over the 10 % cap, so there are fewer than 19 windows; and for very fast signals the windows are longer than the 30 · B rule asks, so there are fewer of them.

One answer per window, so short events can be missed. aTDS uses longer windows for slow signals (259 samples instead of 60 in Figure 1; 11 minutes on the EEG pair). A maximum window can be set, trading reliability for finer timing. For records whose dynamics change, setting the parameters per segment helped (Table 2).

The delay is pinned down only within the tolerance. The calibrated tolerance grows with the time scale (±13 in Figure 1, up to ±56 at r = 16) and can exceed the delay itself. The reported delay was right in Figure 4, but the stability judgement is coarse.

5 % is a budget for the mean, not a significance level. With rotations pooled over several recordings and the upper-bound rule, the mean chance score on new data stayed ≤ 5 % in all rows except one (5.2 %, fast–slow pair with the max-B window). Single recordings still exceed it by luck (95th percentile up to 20 %), so each link needs its own test (Table 3).

B depends on the record length for long-memory signals. For H = 0.9, B grew from 13 to 84 as the record grew from 5 000 to 80 000 samples, so the window depends on how long one records.

Rhythms and bursts. Delays longer than half a cycle are reported wrongly (as the anti-phase lobe); for Rulkov neurons the delay is reported 3 steps too long. Such signals should first be turned into cycle-by-cycle series.

One window for all pairs. Setting the window from the slowest signal keeps all scores on one ruler, but fast pairs then get a wide tolerance (±56 instead of ±3). For pairs whose two signals differ in speed, the pair's joint Bxy gave the same detection with a much shorter window (Table 2).

Scores are not comparable with published TDS numbers. Compare an aTDS score with the aTDS chance level of the same data.

## Not tested here

Coupling that switches with the regime. No test switches the coupling itself on in one regime and off in another (the sleep-stage case); the changing-dynamics test changes only the signals' memory.

Sleep stages in real cohorts. Whether parameters should be set per stage or frozen per cohort, and how that changes stage contrasts, needs staged multi-subject recordings, which are not in this repository. The synthetic changing-dynamics test (Table 2) is only a first step.

The original surrogate test. Bashan et al. (2012) build the null by mixing signals across subjects within the same sleep stage. That null keeps stage-driven co-modulation, which rotating one signal within a recording removes, so the rotation null can be too lenient on real data. Comparing both nulls needs multi-subject data.

## Next steps

A public staged sleep cohort: per-stage networks under aTDS and the published setting, with cross-subject, stage-matched surrogates next to the rotation null.

Make the pair's joint Bxy the window rule for pairs of different speed, and retest Figures 3–7 with it.

Networks with stronger common drivers and chains, to measure how often indirect links appear, and whether partial (conditional) TDS is needed.

A regime-switching test in which the coupling itself is on in one regime and off in another, scored per regime (the sleep-stage case).

A standard error clustered by recording in the tolerance rule, since rotations of one recording are not independent.

A weaker-coupling fast–slow pair (to compare the max-B and Bxy windows on detection, not only on window length), and long-memory signals with stronger coupling.

## Methods

Signals. AR(1): x has correlation time τc (φ = e−1/τc), z-scored; the coupled partner is y(t) = c · x(t − d) + √(1 − c²) · η(t), where η is an independent AR(1) process with the same τc; uncoupled pairs use c = 0. Long memory: fractional Gaussian noise (exact circulant-embedding method) with Hurst exponent H. Rhythms: AR(2) with period 40 and pole radius 0.97. Rulkov: two Rulkov maps, neuron 1 driven by neuron 0 through delayed diffusive coupling. Changing dynamics: two AR(1) halves with τc 4 and 16. Real data: the δ and σ power series in pyTDS/data. Records are 20 000 samples long unless stated (30 000 in Figures 1, 2 and 5; 40 000 in Figure 3b; 6 000 for Rulkov; 10 000 in Table 3). Every pair has its own fixed random seed.

aTDS v2.1. As in the algorithm section, with two changes made in this revision: the minimum window is 61 samples (was 20), optionally max(61, 2·dmax + 1), and the tolerance is the largest value for which the mean chance score plus 2 standard errors, over rotations pooled across all available training recordings, stays ≤ 5 % (was: the mean over 20 rotations of one recording). In every figure and table, aTDS sets its parameters on separate training data that are not scored: 5 training pairs in the synthetic tests, the first half of the record for the EEG pair. In Table 3 the calibration is label-free on the scored network itself (rotations only, no knowledge of which links are true), as it would be in a real analysis.

Baselines. Published: window 60, step 30, lag ±30 (±29 effective), tolerance ±1. Rescaled by hand: window 60r, step 30r, lag ±30r, tolerance ±r, using the true slow-down factor r. Oracle: windows 30–1 920 × tolerances 0–135, chosen by AUC on 20 + 20 labelled training pairs.

TDS. pyTDS, circular cross-correlation, z-scored windows, τ₀ = lag of the largest |C(τ)|, stability rule 4 of 5 windows, default (non-causal) stability anchor.

Statistics. AUC: probability that a coupled pair outscores an uncoupled one. Intervals: 1 000 bootstrap resamples of the pairs (paired across methods for differences); Wilson intervals for shares. Equivalence: two one-sided tests, i.e. the 90 % interval of the AUC difference inside ±0.05. Delay found: at least 3 stable windows and their median τ₀ within 10 % (at least ±1 sample) of the true delay. Network links: p = (1 + number of rotations scoring at least the real score) / 1 000, Benjamini–Hochberg at q = 0.05. Rulkov delay: median τ₀ of the stable windows, as for all other rows.

Scope. Synthetic signals with known answers and one real benchmark pair. No staged sleep cohort and no cross-subject surrogates.

Bashan A, Bartsch RP, Kantelhardt JW, Havlin S, Ivanov PCh. Network physiology reveals relations between network topology and physiological function. Nature Communications 3:702 (2012). · Bartlett MS. On the theoretical specification and sampling properties of autocorrelated time-series. J. R. Stat. Soc. Suppl. 8:27–41 (1946).


## Raw data: slow-down test (exp4.json U1), all four settings

| r | setting | AUC | AUC 95% CI | delay found | found 95% CI | mean chance score % | window / lag / tol |
|---|---|---|---|---|---|---|---|
| 1 | fixed | 0.997 | [0.989, 1.0] | 1.0 | [0.963, 1.0] | 0.78 | 60 / 30 / 1 |
| 1 | rescaled | 0.997 | [0.989, 1.0] | 1.0 | [0.963, 1.0] | 0.78 | 60 / 30 / 1 |
| 1 | atds | 0.995 | [0.988, 1.0] | 1.0 | [0.963, 1.0] | 4.36 | 66 / 32 / 3 |
| 1 | atds_tol1 | 1.0 | [1.0, 1.0] | 1.0 | [0.963, 1.0] | 0.66 | 66 / 32 / 1 |
| 1 | oracle | 1.0 | [0.999, 1.0] | 1.0 | [0.963, 1.0] | 0.03 | 60 / 29 / 0 |
| 1 | AUC difference aTDS - rescaled | 90% CI [-0.009, 0.006] | 95% CI [-0.011, 0.007] | | | | |
| 2 | fixed | 0.888 | [0.839, 0.93] | 0.84 | [0.756, 0.899] | 0.64 | 60 / 30 / 1 |
| 2 | rescaled | 0.993 | [0.984, 0.999] | 1.0 | [0.963, 1.0] | 0.52 | 120 / 60 / 2 |
| 2 | atds | 0.979 | [0.961, 0.992] | 0.99 | [0.946, 0.998] | 3.63 | 126 / 62 / 6 |
| 2 | atds_tol1 | 0.988 | [0.969, 0.999] | 0.97 | [0.915, 0.99] | 0.12 | 126 / 62 / 1 |
| 2 | oracle | 1.0 | [1.0, 1.0] | 1.0 | [0.963, 1.0] | 0.0 | 240 / 119 / 0 |
| 2 | AUC difference aTDS - rescaled | 90% CI [-0.029, -0.002] | 95% CI [-0.032, 0.001] | | | | |
| 4 | fixed | 0.59 | [0.513, 0.658] | 0.18 | [0.117, 0.267] | 0.72 | 60 / 30 / 1 |
| 4 | rescaled | 0.955 | [0.926, 0.978] | 0.91 | [0.838, 0.952] | 0.39 | 240 / 120 / 4 |
| 4 | atds | 0.923 | [0.885, 0.959] | 0.88 | [0.802, 0.93] | 3.65 | 237 / 118 / 12 |
| 4 | atds_tol1 | 0.885 | [0.84, 0.925] | 0.77 | [0.678, 0.842] | 0.0 | 237 / 118 / 1 |
| 4 | oracle | 1.0 | [1.0, 1.0] | 1.0 | [0.963, 1.0] | 0.0 | 960 / 479 / 0 |
| 4 | AUC difference aTDS - rescaled | 90% CI [-0.064, 0.003] | 95% CI [-0.073, 0.008] | | | | |
| 8 | fixed | 0.541 | [0.464, 0.618] | 0.0 | [0.0, 0.037] | 0.59 | 60 / 30 / 1 |
| 8 | rescaled | 0.806 | [0.752, 0.857] | 0.6 | [0.502, 0.691] | 0.34 | 480 / 240 / 8 |
| 8 | atds | 0.877 | [0.829, 0.923] | 0.85 | [0.767, 0.907] | 3.31 | 563 / 281 / 27 |
| 8 | atds_tol1 | 0.83 | [0.785, 0.875] | 0.66 | [0.563, 0.745] | 0.0 | 563 / 281 / 1 |
| 8 | oracle | 0.985 | [0.965, 1.0] | 0.97 | [0.915, 0.99] | 0.0 | 1920 / 959 / 0 |
| 8 | AUC difference aTDS - rescaled | 90% CI [0.018, 0.131] | 95% CI [0.007, 0.142] | | | | |
| 16 | fixed | 0.509 | [0.428, 0.585] | 0.0 | [0.0, 0.037] | 0.91 | 60 / 30 / 1 |
| 16 | rescaled | 0.688 | [0.631, 0.737] | 0.41 | [0.319, 0.508] | 0.42 | 960 / 480 / 16 |
| 16 | atds | 0.645 | [0.58, 0.716] | 0.42 | [0.328, 0.518] | 4.11 | 1083 / 541 / 56 |
| 16 | atds_tol1 | 0.63 | [0.585, 0.675] | 0.26 | [0.184, 0.354] | 0.0 | 1083 / 541 / 1 |
| 16 | oracle | 0.843 | [0.788, 0.891] | 0.73 | [0.636, 0.807] | 4.95 | 1920 / 959 / 135 |
| 16 | AUC difference aTDS - rescaled | 90% CI [-0.105, 0.018] | 95% CI [-0.117, 0.029] | | | | |

(atds_tol1 = aTDS window and lag range with tolerance fixed at +-1; oracle = best of 7 windows x 13 tolerances by AUC on 20+20 labelled training pairs.)

## Raw data: window length test (exp4.json U2; r=2 is tau_c=4, delay 10; r=8 is tau_c=16, delay 40; coupling 0.3; 100+100 pairs)

| system r | window L | calibrated tol | AUC | 95% CI |
|---|---|---|---|---|
| 2 | 30 | 1 | 0.589 | [0.507, 0.666] |
| 2 | 60 | 2 | 0.829 | [0.768, 0.88] |
| 2 | 120 | 5 | 0.974 | [0.956, 0.99] |
| 2 | 240 | 13 | 1.0 | [0.999, 1.0] |
| 2 | 480 | 27 | 1.0 | [1.0, 1.0] |
| 2 | 960 | 51 | 1.0 | [1.0, 1.0] |
| 2 | 1920 | 116 | 1.0 | [1.0, 1.0] |
| 8 | 30 | 1 | 0.469 | [0.395, 0.551] |
| 8 | 60 | 2 | 0.513 | [0.432, 0.597] |
| 8 | 120 | 6 | 0.518 | [0.438, 0.602] |
| 8 | 240 | 13 | 0.659 | [0.586, 0.737] |
| 8 | 480 | 26 | 0.755 | [0.689, 0.818] |
| 8 | 960 | 47 | 0.915 | [0.874, 0.952] |
| 8 | 1920 | 101 | 0.99 | [0.976, 0.999] |
| 2 | 121 | 5 | 0.97 | [0.947, 0.988] |
| 8 | 515 | 25 | 0.79 | [0.728, 0.85] |

(The last two rows are the windows aTDS chose: 121 and 515.)

---

# PART 2: REAL DATA (all TDS computed with the author's pyTDS; scripts and raw results in experiments/atds_real/)

## R1. Benchmark reproduction
pyTDS with the published setting and the causal stability anchor gives 43.53 % on the benchmark delta-sigma sleep-EEG pair (data/data-delta.txt, data/data-sigma.txt, 29 941 s at 1 Hz); the original MATLAB implementation gives 43.6052 %. With the default (non-causal) anchor pyTDS gives 58.1 %. Every other analysis uses the default anchor.

## R2. Sleep cohort: 35 overnight polysomnograms (author's MSc project data, prepared to reproduce the network analysis of Bashan et al. 2012)
Data: 35 recordings, 23 150-31 752 s (median 28 680 s), 10 channels at 1 Hz: heart rate, respiratory rate, chin EMG std, leg EMG std, EOG std, delta, theta, alpha, sigma, beta EEG power. Per-second sleep stages: wake (W), light sleep (LS), deep sleep (DS), REM. 45 links per subject.

Settings (one parameter set per setting for the whole cohort):
- published: window 60 s, step 30 s, lag +-30 s, tolerance +-1 s.
- aTDS (uncapped): window from the median over subjects of max_i B_i = 228 -> 30 x 228 = 6 846 s, clipped to 10 % of the median length: L = 2 868 s (48 min), step 1 434, lag +-1 433 s, tolerance +-136 s (calibrated on rotations pooled over 50 random subject-pair combinations, 20 rotations each; calibration did not fail).
- aTDS capped at 300 s: L = 300 s, step 150, lag +-149 s, tolerance +-14 s.

Stage assignment: delay series and stability labels computed over the whole night; a window is assigned to a stage if >= 80 % of its seconds carry that stage; per-stage TDS = % stable among assigned windows; at least 3 windows required.
Significance: cross-subject, stage-matched null as in Bashan et al. 2012: per link, 100 surrogate pairs (channel i of subject a, channel j of subject b, a != b, truncated to the shorter length, stages of subject a); link significant if TDS > 95th percentile of its null (per link, stage, setting).
Effects: paired Cohen dz across subjects, 95 % bootstrap intervals (5 000 resamples of subjects); gain = paired difference of dz between settings.

Results:
- Uncapped aTDS: no subject had >= 3 windows that are >= 80 % wake, REM or deep sleep; only 24 subjects had >= 3 in light sleep. The uncapped window cannot resolve sleep stages.
- Mean over subjects (published / aTDS-300): mean link TDS W 16.90 / 40.15, LS 19.76 / 46.91, REM 11.44 / 36.64, DS 13.73 / 35.82 %; significant links (of 45) W 25.23 / 26.63, LS 34.20 / 38.20, REM 23.66 / 28.51, DS 18.31 / 23.03. Subjects with enough windows under aTDS-300: LS 35, REM 35, DS 34, W 27 (8 subjects lack three pure 300-s wake windows). Raw scores are not comparable across settings (aTDS tolerance is wider); significant links against each setting's own null are the fair comparison.

| Measure | Contrast | Published dz [95 % CI] (% subjects positive) | aTDS-300 dz [95 % CI] (% positive) (n) | Gain aTDS-300 minus published [95 % CI] |
|---|---|---|---|---|
| Mean link TDS | LS-DS | 1.37 [1.02, 1.91] (91 %) | 1.65 [1.27, 2.25] (94 %) (34) | +0.31 [-0.06, 0.75] |
| | W-DS | 0.57 [0.25, 0.97] (69 %) | 0.29 [-0.09, 0.76] (65 %) (26) | -0.28 [-0.70, 0.13] |
| | REM-DS | -0.50 [-0.88, -0.17] (40 %) | 0.11 [-0.25, 0.42] (47 %) (34) | +0.67 [0.40, 0.98] |
| | LS-REM | 1.90 [1.53, 2.53] (100 %) | 1.37 [0.97, 1.97] (91 %) (35) | -0.53 [-1.18, 0.05] |
| Significant links | LS-DS | 2.58 [2.11, 3.41] (100 %) | 2.47 [1.95, 3.34] (100 %) (34) | -0.14 [-0.74, 0.47] |
| | W-DS | 1.16 [0.69, 1.88] (89 %) | 0.34 [-0.04, 0.95] (69 %) (26) | -1.17 [-1.65, -0.77] |
| | REM-DS | 0.75 [0.45, 1.13] (74 %) | 0.88 [0.58, 1.25] (74 %) (34) | +0.15 [-0.20, 0.50] |
| | LS-REM | 1.67 [1.26, 2.40] (94 %) | 1.75 [1.35, 2.39] (94 %) (35) | +0.08 [-0.48, 0.54] |

Reading: both settings show the strong LS > DS contrast; aTDS-300 gives no significant gain. Under the published setting REM has LOWER mean link TDS than DS (dz -0.50) but MORE significant links (dz 0.75); aTDS-300 removes the mean-TDS inversion (gain +0.67) and leaves the significant-link contrast unchanged. The W-DS contrast is weaker under aTDS-300 (significant links gain -1.17) and rests on 26 subjects. On real staged data aTDS does not uniformly sharpen sleep-stage contrasts; it needs a window cap, and per-stage calibration was not tested (parameters frozen per cohort). Tolerance calibration still used within-recording rotations.

## R3. Three public records from other fields with known delays
Delays known from the domain and from the peak of the full-record cross-correlation (sources and preprocessing in experiments/atds_real/data_fields/MANIFEST.json). aTDS calibrated on rotations of each record itself (20 rotations); significance from 200 rotations, p = (1+k)/(1+N). Delay = most frequent tau0 among stable windows, signed so the driver leads.

| Record | Samples | Known delay | Setting | L / lag / tol | Score (null95) | p | Delay found |
|---|---|---|---|---|---|---|---|
| Hourly catchment rainfall -> river flow (airGR catchment L0123003, hydrology) | 43 848 h | about 8 h | published | 60 h / 30 h / 1 h | 1.2 % (1.1) | 0.02 | 13 h (wrong) |
| | | | aTDS | 79.8 d / 39.8 d / 1.79 d | 79.5 % (9.1) | 0.005 | 8 h (correct) |
| Melbourne temperature -> Victoria electricity demand, summer, mean daily profile removed (energy) | 13 008 half-hours | about 1 h | published | 30 h / 15 h / 30 min | 11.8 % (3.0) | 0.005 | 0 (wrong) |
| | | | aTDS | 27.1 d / 13.5 d / 1.02 d | 100 % (21.1) | 0.005 | 1 h (correct) |
| Weekly changes, WTI crude oil -> US retail gasoline (economics) | 1 746 weeks | 1-2 weeks | published | 60 / 30 / 1 weeks | 86.0 % (7.0) | 0.005 | 0 weeks |
| | | | aTDS | 64 / 31 / 3 weeks | 92.5 % (15.1) | 0.005 | 0 weeks |

Reading: both settings detect coupling in all three (published only marginally for rainfall-runoff). aTDS reports the known delay for rainfall-runoff and temperature-demand; the published setting does not. Oil-gasoline does not discriminate: in weekly changes the cross-correlation already peaks at lag 0 (r = 0.55, partly because WTI and retail prices are recorded on different weekdays; lag 1: r = 0.39). aTDS windows are long (80 and 27 days), acceptable because these delays are constant.

---

# PART 3: ADDITIONAL EXPERIMENTS (baselines, sensitivity, harder benchmarks, stage-aware cohort analysis, more records; scripts and raw results in experiments/atds_v3/; all TDS with pyTDS)

## X1. Baselines beyond TDS settings (slow-down benchmark of Part 1: tc = 2r, d = 5r, c = 0.3, 100 + 100 pairs per r, same seeds)
New methods: xcorr_full = full-record normalised cross-correlation, |lag| <= 1000, score = max |C|, delay = argmax; gccphat_full = full-record GCC-PHAT, same lag range; rescaled_est = published setting multiplied by r_hat = max(1, round(tc_hat/2)), tc_hat = median over the 5 training pairs of the first lag where the autocorrelation falls below 1/e (r_hat equalled the true r for r = 1-8 and was 15 at r = 16). Recomputed published / rescaled / aTDS AUCs reproduce Part 1 exactly.

AUC by r (1 / 2 / 4 / 8 / 16):
- published .997 / .888 / .590 / .541 / .509
- rescaled (true r) .997 / .993 / .955 / .806 / .688
- aTDS .995 / .979 / .923 / .877 / .645
- rescaled_est .997 / .993 / .955 / .806 / .691 [.64, .74 at r=16]
- xcorr_full 1.0 at every r; gccphat_full 1.0 at every r
Delay found: aTDS 1.0 / .99 / .88 / .85 / .42; published 1.0 / .84 / .18 / 0 / 0; xcorr_full and gccphat_full 1.0 at every r.
aTDS minus baseline AUC [95 % CI]: vs published +.091 [.05, .14] (r=2), +.333 [.24, .42] (r=4), +.336 [.25, .44] (r=8), +.136 [.03, .25] (r=16); vs rescaled -.014 [-.03, .00], -.032 [-.07, .01], +.071 [.007, .142], -.043 [-.12, .03]; vs xcorr_full -.021 [-.04, -.01], -.077 [-.12, -.04], -.123 [-.17, -.08], -.355 [-.43, -.29]. At r = 1 all differences include 0 except vs xcorr_full (-.005 [-.013, -.001]).
Reading: on this stationary benchmark (one constant delay, coupling always on) full-record estimators are perfect and beat every windowed method, and a published setting rescaled by an estimated time scale is as good as aTDS. TDS (and aTDS) give a time-resolved answer, which full-record estimators cannot; see X3.

## X2. Sensitivity of the design constants (slow-down benchmark, r = 1 / 4 / 8, 100 + 100 pairs; one constant varied at a time; default: n_eff 30, budget 5 %, z 2, floor 61, 20 rotations; default AUC .995 / .923 / .877, chance 3.3-4.4 %)
- n_eff 15: AUC .997 / .717 / .684; n_eff 45: 1.0 / .946 / .899; n_eff 60: 1.0 / .991 / .954 (chance rises to 5.0-5.4 % at n_eff 45-60 for r = 4, 8; windows longer).
- budget 1 %: 1.0 / .920 / .884 (chance 0.2-0.7 %); budget 10 %: .981 / .874 / .851.
- z 0: .995 / .911 / .861; z 3: .995 / .923 / .881.
- floor 31: identical to 61; floor 121: 1.0 / .923 / .877.
- rotations 10: .995 / .921 / .888; rotations 40: .995 / .923 / .862.
- Rhythm threshold (AR(2) rhythm, period about 40, delays 8 and 30, 30 + 30 pairs, window 316-318): threshold 0.2 and 0.3 -> rhythm detected, lag capped at +-20, tolerance 1: AUC 1.0 and .995, delay found 1.0 and 0 (delay 30 out of range); threshold 0.4 -> rhythm not detected, lag +-158, tolerance 16: AUC 1.0 and 1.0, delay found 1.0 and 1.0.
Reading: results are robust to the budget, z, floor and number of rotations; the number of independent samples per window matters most (more = better detection, coarser time resolution, slightly higher chance); the rhythm cap trades delay range for alias protection.

## X3. Harder synthetic benchmarks (20 000 samples)
1. Regime switching (the sleep-stage case): blocks of 2 000 samples alternating regime A (tc 8; y = 0.6 x(t-20) + noise) and regime B (tc 32; uncoupled); controls uncoupled throughout; 40 + 40 pairs; scores per regime with the 80 % purity rule.
   - published (60/30/30/1): regime-A score coupled 5.62 % vs control 0.71 %; AUC .944 [.89, .99]; delay found .95; regime-B false stability 0.89 % (control 0.81 %).
   - aTDS whole record (window 516, lag 257, tol 29; the 1000 cap was not binding): regime-A score 98.94 % vs control 4.39 %; AUC 1.0; delay found 1.0; regime-B false stability 8.12 % (control 5.78 %).
   - aTDS per regime (A: window 251, tol 13; B: window 1000, tol 53): regime-A score 84.33 % vs 3.60 %; AUC 1.0; delay found 1.0; regime-B false stability 9.00 % (control 4.50 %).
   Reading: aTDS detects regime-limited coupling far more strongly (99 % vs 5.6 % of regime-A windows stable) but its false stability in the uncoupled regime (8-9 %) exceeds the 5 % budget, because windows near regime boundaries inherit stability.
2. Nonlinear coupling (tc 8, d 20, 100 + 100 pairs): monotone (tanh): published AUC .997, aTDS 1.0, xcorr_full 1.0, delay found 1.0 for all. Non-monotone (quadratic; zero linear correlation): published AUC .803 (delay found .74), aTDS .747 (.62), xcorr_full .341 (.01). Cross-correlation-based TDS is not designed for non-monotone coupling; aTDS does not help there.
3. Drifting delay (tc 8, c 0.6, delay drifting linearly 10 -> 30 samples, 40 pairs): published score 8.5 %, median absolute error of tau0 in stable windows 0.30 samples; aTDS score 86.5 %, error 0.25 samples (capped at 1000: identical). aTDS tracks a slowly drifting delay with nearly all windows stable.

## X4. Stage-aware aTDS on the 35-subject sleep cohort (same data, stage assignment, 100-draw cross-subject stage-matched null and bootstrap as Part 2)
Stage-episode lengths, median [IQR] in s: W 90 [60-210] (331 episodes), LS 420 [120-1118] (708), DS 480 [150-1050] (229), REM 480 [180-900] (305); all stages 300 [90-870] (1 573).
Setting atds_stage: window = clip(0.5 x median episode length (all stages), 61, 300) = 150 s, lag +-74 s; tolerance per stage calibrated on 50 cross-subject stage-matched surrogates (mean + 2 SE <= 5 %): W 6, LS 6, DS 6, REM 4 s; no calibration failed.
Per-stage mean link TDS (W / LS / DS / REM): published 16.9 / 19.8 / 13.7 / 11.4; atds_cap300 40.2 / 46.9 / 35.8 / 36.6; atds_stage 31.8 / 36.4 / 27.0 / 22.6.
Effects, atds_stage dz [95 % CI] (n): mean link TDS LS-DS 1.62 [1.29, 2.13] (34); W-DS 0.53 [0.20, 0.97] (33); REM-DS -0.70 [-1.17, -0.36] (34); LS-REM 2.16 [1.63, 3.02] (35). Significant links LS-DS 2.55 [2.07, 3.37]; W-DS 0.47 [0.13, 0.92]; REM-DS 0.80 [0.48, 1.24]; LS-REM 1.43 [0.98, 2.24].
Gains atds_stage minus published: mean link TDS LS-DS +0.28 [-0.01, 0.57], W-DS -0.05 [-0.33, 0.22], REM-DS -0.13 [-0.35, 0.03], LS-REM +0.26 [-0.21, 0.88]; significant links LS-DS -0.05 [-0.72, 0.57], W-DS -1.02 [-1.48, -0.68], REM-DS +0.07 [-0.31, 0.51], LS-REM -0.25 [-0.68, 0.25].
Reading: the stage-aware setting keeps wake usable (33 subjects instead of 26 with the 300-s cap) and reproduces all published contrasts with the same signs, but gives no significant gain over the published setting and a significant loss on the wake contrast in significant links. The published 60-s setting suits this 1-Hz sleep data, for which it was designed.

## X5. More records in one field: three airGR catchments (rainfall -> river flow; 200 rotations)
| Catchment | Expected delay | Setting | Score (null95) | p | Delay found |
|---|---|---|---|---|---|
| L0123001, daily (10 593 days) | 0-3 days (domain); xcorr peak +2 d | published | 59.1 % (2.6) | 0.005 | 1 day (in range) |
| | | aTDS (window 501 d) | 97.6 % (12.3) | 0.005 | 2 days (in range) |
| L0123002, daily | no domain expectation; xcorr peak at the edge of the search (unreliable) | published | 32.4 % (2.3) | 0.005 | 1 day |
| | | aTDS | 21.1 % (21.1) | 0.104 | 162 days (implausible) |
| L0123003, hourly (Part 2) | about 8 h | published | 1.2 % (1.1) | 0.02 | 13 h (wrong) |
| | | aTDS | 79.5 % (9.1) | 0.005 | 8 h (correct) |
Reading: aTDS recovers the delay where the published window is far too short for the time scale (hourly data), matches the published setting on daily data, and fails on L0123002 (not significant, implausible delay; the full-record cross-correlation is also unreliable there).
