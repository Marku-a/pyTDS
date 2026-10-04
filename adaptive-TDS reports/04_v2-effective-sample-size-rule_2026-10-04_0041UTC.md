# Report 4: The improved rule (v2), based on effective sample size

*Published 2026-10-04 00:41 UTC. Same experiments as Reports 2 and 3, re-run with two new variants. Numbers: `experiments/adaptive/results/`. Figures: `figures/`.*

## What changed from v1 and why

- **v1:** window = 10 × the slowest signal's decorrelation time. Report 3 showed this is too short for fast and spiky signals.
- **v2:** window = 30 × the **Bartlett factor**. The Bartlett factor measures how many samples count as "one independent sample", because neighbouring samples of a smooth signal repeat each other. v2 therefore sizes the window to hold **about 30 independent samples**, using the slowest signal in the system. Unlike the v1 measure, it sums the *whole* autocorrelation, so slow structure counts. Step, max lag and the two tolerance options are as before.

The two new variants are **v2-fixedtol** (tolerance ±1) and **v2-calibrated** (tolerance chosen so the chance score is ≤ 5%).

## Ground-truth results

Each cell below is AUC / hit rate.
- **AUC:** how well the TDS score separates coupled from uncoupled pairs (1 = perfect).
- **Hit rate:** how often the true delay is found.

**E1: growing time scale, weak coupling** (the hardest case):

| scale r | default | v1 | v2-fixedtol | v2-calibrated | oracle (peeks at answers) |
|---|---|---|---|---|---|
| 1 | 1 / 1 | 0.90 / 0.85 | **1 / 1** | **1 / 1** | 1 / 1 |
| 2 | 0.89 / 0.9 | 0.79 / 0.6 | **1 / 1** | 0.99 / 1 | 1 / 1 |
| 4 | 0.60 / 0.1 | 0.78 / 0.6 | **1 / 1** | 0.95 / 1 | 1 / 1 |
| 8 | 0.35 / 0 | 0.51 / 0.25 | **0.92 / 0.85** | 0.86 / 0.9 | 1 / 1 |
| 16 | 0.44 / 0 | 0.53 / 0.4 | 0.68 / 0.35 | **0.74 / 0.65** | 0.92 / 0.85 |

With strong coupling, both v2 variants score 1 / 1 at every scale. The default drops to about 0.55 / 0 at r = 8–16.

**Other ground-truth tests:**

| test | default | v1 | v2-fixedtol | v2-calibrated |
|---|---|---|---|---|
| E1c: signals of different speeds (worst of 4 combinations, AUC) | 0.76 | 0.78 | **1.00** | **1.00** |
| E1d: on/off coupling, fast signals (balanced accuracy) | 0.89 | 0.56 | **0.91** | 0.90 |
| E1d: on/off coupling, slower signals | 0.53 | 0.61 | **0.90** | 0.88 |
| E1b: 8-signal network, edge AUC | 0.87 | **0.98** | 0.97 (v2-calibrated) | – |
| E2: Rulkov neurons, delay 32 (AUC / hit) | 0 / 0 | 0 / 0 | **1 / 1** | **1 / 1** |
| E2: Rulkov, delay 16 (AUC / hit) | 1 / **0.9** | 0.64 / 0 | 1 / 0 | 1 / 0 |

- **E1d** asks whether TDS marks the right moments as coupled; 0.5 means guessing.
- **In E1b** all adaptive variants found every delay (hit 1.0), against 0.6 for the default.
- **Rulkov:**
  - **Delay 16:** v2 separates coupled from uncoupled perfectly but reports a delay of 19–20 instead of 16, which counts as a miss. The bursts repeat about every 41 steps, so the oscillation guard caps the max lag at 20, and the true delay sits near that edge.
  - **Delay 32:** only v2 detects it at all.

## Real data

**R4: sampling-rate invariance**, the same EEG resampled ×0.5 to ×4. Scores in %:

| variant | ×0.5 | ×1 | ×2 | ×4 |
|---|---|---|---|---|
| default | 75 | 58 | 41 | 27 |
| v2-calibrated | 100 | 100 | 100 | 100 |

This invariance is **partly misleading**. On this EEG, v2 chose a 29-minute window, and with only 33 windows in 8 hours the score saturates at 100%. The next test shows what happens when you cap the window.

**User window cap on the EEG** (v2-calibrated):

| max window | windows in the night | score | chance score (mean / 95th pct) |
|---|---|---|---|
| 2 min | 498 | 79% | 5% / 8.5% |
| 5 min | 198 | 90% | 5% / 10% |
| 10 min | 98 | 94% | 5% / 12% |
| no cap (29 min) | 33 | 100% | 6% / 24% |

- The calibration keeps the chance score near 5% whatever the window. So **every row is statistically valid**; they differ only in time resolution.
- Sleep stages change over minutes, so for sleep EEG you would cap the window at a few minutes. **That cap is the one parameter that should stay with the user**, because it encodes what you want to resolve, not something the data can tell you.

**R2: ECG → PPG, 100 Hz.**
- v2 scores 89% against a 7% null95. The default gets 42% against 10.5%, and v1 only 3%.
- The delay found is still −80 ms rather than the physical 430 ms, the rhythm limit explained in Report 2.
- Across 100/50/25 Hz, **v2-calibrated keeps chance stability low** (null95 3–9%), whereas v2-fixedtol's null95 rises to 42% at 25 Hz. This is where calibrating the tolerance earns its place.

**R3: heart rate vs. breathing at 4 Hz.** v2 chose window 59 with tolerance 1, essentially the default. It scores the same as v1 and separates real from mismatched pairs perfectly.

## Verdict on v2

- **Better than or equal to the default in every ground-truth test except one delay-accuracy case** (Rulkov at delay 16). It is often dramatically better: on slow signals, mixed speeds, on/off coupling and long delays.
- **Close to the oracle** (the best possible parameters, found by peeking) everywhere except the slowest, weakest case.
- **Matches the defaults automatically** when the data is already at the scale the defaults were made for (R3).
- **Calibrated vs. fixed tolerance:** they are similar on clean synthetic data. Calibration matters on rhythmic real data, where it keeps false stability under control, so it is the safer default.
- **One user input remains: the longest window you accept**, which is the time resolution you need.
