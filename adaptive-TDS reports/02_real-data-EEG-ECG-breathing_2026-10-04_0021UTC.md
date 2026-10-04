# Report 2: Real data (sleep EEG, ECG→PPG pulse, heart rate↔breathing)

*Published 2026-10-04 00:21 UTC. Code: `experiments/adaptive/exp_real.py`. Raw numbers: `experiments/adaptive/results/real_R*.json/.csv`. Figures: `figures/R*.png`.*

## How to read this report

Each signal pair is scored three ways:
- **default:** window 60, step 30, max lag 30, tolerance ±1.
- **adaptive-fixedtol:** window, step and max lag from the data; tolerance stays ±1.
- **adaptive-calibrated:** tolerance also chosen from the data, to keep the chance score at or below 5%.

Report 1 explains how these are derived.

For every real pair I also built a **null**: the second signal is rotated in time 200 times, which destroys any real coupling, and TDS is re-run each time. "null95" is the score that only 5% of these coupling-free versions exceed. A real score far above null95 means the coupling is real, not chance.

**Limitation:** real data has no exact ground truth for "how much" coupling there is. So this report can show **consistency** (does the answer stay the same when only the sampling rate changes?) and **real vs. chance**. It cannot show which score is "more correct". The synthetic tests (Report 3) answer that.

## R1: Sleep EEG, delta vs. sigma band power (`data/`, 1 Hz, 8.3 h)

| variant | window / step / max lag / tolerance (s) | score | null95 |
|---|---|---|---|
| default | 60 / 30 / 30 / ±1 | 58.1% | 2.5% |
| adaptive-fixedtol | 337 / 168 / 168 / ±1 | 80.8% | 0.0% |
| adaptive-calibrated | 337 / 168 / 168 / ±19 | 89.8% | 12.5% |

- **The two signals live on very different time scales.** The delta band "remembers" about 34 s of its past; the sigma band about 3 s. The adaptive rule takes the slower one, so the window grows to 5.6 minutes.
- **Every variant says the coupling is real** (all p = 0.005, the smallest possible with 200 shifts), and the stable delay is 0 s in all variants.
- **Scores rise with the adaptive window.** That is partly real, because longer windows estimate delays more reliably. It is also partly "smoothing": the null95 of the calibrated variant rises to 12.5% too. A higher score is therefore *not* automatically "better".
- **Side finding (not about the adaptive idea):** `CLAUDE.md` says the default reproduces the MATLAB benchmark of 43.6%. With the current defaults it gives **58.1%**. The 43.5% value comes only with `stability_anchor="end"`. The repo's test only checks that the score is between 0 and 100, so nobody noticed. I checked this myself.

## R4: The same EEG, resampled (the cleanest real-data test)

I artificially changed the sampling rate (×0.5, ×2, ×4). The physiology is identical, so a good method should give the same answer.

| sampling factor | 0.5 | 1 | 2 | 4 |
|---|---|---|---|---|
| default | 74.7 | 58.1 | 41.3 | 26.8 |
| adaptive-fixedtol | 93.7 | 80.8 | 76.4 | 70.7 |
| adaptive-calibrated | 94.7 | 89.8 | 90.2 | 90.2 |

- **The default's score depends on how the data was recorded:** 75% → 27% for the same sleep night.
- **Adaptive-calibrated is almost invariant**, which is what you want if TDS is to be used across labs and fields with different sampling rates.
- **Adaptive-fixedtol is only partly invariant.** ±1 sample means a different amount of time at each rate, so the tolerance must adapt too.

## R2: ECG → finger pulse (PPG), 100 Hz, 5 min

The "tag" here is the **pulse arrival time**: the delay from the heart's electrical beat (R peak in the ECG) to the pulse peak at the finger. I measured it directly on 430 beats: **430 ms** (IQR 430–440).

**Result: no variant recovers 430 ms.** All of them find a stable delay of about −80 ms.

The reason is a fundamental limit, not a bug. The heart beats every ~690 ms, so the waveforms repeat. A delay of 430 ms looks exactly like a delay of −260 ms (one beat earlier). For repeating signals, cross-correlation can only see delays up to **half a period** (~345 ms here). The adaptive rule detected the beat period and correctly capped the max lag at 350 ms, but that cap also makes the true value unreachable. (Why exactly −80 ms: probably the correlation lining up other parts of the waveform, such as the ECG T wave, with the pulse. I did not investigate this further.)

**Lesson:** for rhythmic signals TDS measures "phase", not absolute delay. No window choice fixes that. You would need beat-to-beat series, such as R-R intervals vs. pulse amplitudes, instead of raw waveforms.

Two further observations:
- The **calibrated** variant chose tolerance 0, and its score fell to 3% (null95 1.6%). Rhythmic signals look "stable" even by chance (default null95 10–29%), so the calibration tightened the rule to keep chance stability low. It worked as intended, but at this rate it left almost no stable windows.
- Across 100/50/25 Hz, the delay in milliseconds was stable for all variants, while the scores changed with rate for the default and the fixed-tolerance variant.

## R3: Heart rate vs. breathing, 3 recordings, resampled to 4 Hz

Breathing is known to modulate heart rate (respiratory sinus arrhythmia).

| variant | rec 1 | rec 2 | rec 3 | wrong pairings (heart rate of one person + breathing of another) |
|---|---|---|---|---|
| default | 31.6 | 65.2 | 50.0 | 0 (all 6 pairings) |
| adaptive (both versions) | 32.5 | 75.0 | 52.6 | 0 (all 6 pairings) |

- All variants separate real from wrong pairings perfectly, with a stable delay of 0–0.25 s.
- **Adaptive ≈ default here.** At 4 Hz the data happens to sit near the scale the defaults were designed for: the adaptive window was 59 samples vs. the default 60. This is a good sanity check: when the defaults already fit the data, the adaptive rule reproduces them.

## Bottom line from real data

- **Clear win: invariance to sampling rate (R4).** The calibrated adaptive version gives the same answer however the data was sampled; the default does not.
- **Neutral:** when the data is already near the default's scale (R3), adaptive changes almost nothing.
- **Limit found:** rhythmic waveforms (R2) need a different input, not a different window.
- Real data cannot tell us whether the higher adaptive scores are *more accurate*. The synthetic ground-truth tests do that.
