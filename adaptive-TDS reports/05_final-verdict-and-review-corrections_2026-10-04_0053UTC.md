# Report 5: Final verdict (with corrections after an independent code review)

*Published 2026-10-04 00:53 UTC. Read this one first if you only read one.*

## Short answer

**Yes, an adaptive TDS is real and makes the algorithm better.** It must be done carefully, and one parameter should stay with the user.

The best rule found (**v2-calibrated**) works like this:
1. It sizes the window so it holds about 30 *independent* samples of the slowest signal in the system (the Bartlett effective-sample-size formula).
2. It derives the step and the maximum delay from the window, and caps the delay for rhythmic signals.
3. It picks the stability tolerance so that the chance of "fake stability" stays at 5%. That chance is measured on time-shifted copies of the data, which contain no coupling.

On ground-truth data it **matched or beat the fixed defaults in every test except one delay-accuracy case** (Rulkov at delay 16). It was often far better, and it came close to the "oracle" (the best parameters found by peeking at the answers).

## What the review found and how it was fixed

A separate reviewer agent audited the experiment code. It confirmed:
- **No cheating:** the adaptive rules never see which pairs are coupled.
- **Correct signs and metrics.**

It also found four problems. All are fixed, and the numbers below are after the fix.

1. **Network test:** one pair (signals 0 and 5) was counted as "not linked", although it *is* linked indirectly (0 → 1 → 5, total delay 15). It is now left out of the scoring and reported on its own.
2. **Calibration could fail silently.** If even the strictest tolerance (0) gives more than 5% chance stability, the code now flags `calibration_failed`.
   - Synthetic data: it fired only on the Rulkov neurons (always for v1, once for v2).
   - Real data: only v1 on the 25 Hz heart signals.
3. **The Rulkov "trials" were nearly identical copies.** The generator is deterministic and only 1% noise was added. Each trial now starts from random initial conditions.
4. **The signals used to choose the parameters were also among the evaluated trials.** They are now taken from separate random seeds.

The corrections changed some numbers in Reports 3 and 4 (mainly the network test and Rulkov), but **no conclusion changed**. Where they moved results, they made the case for v2 stronger.

## Final numbers (after fixes)

| test (what is known) | default | v1 rule | **v2-calibrated** | oracle |
|---|---|---|---|---|
| Weak coupling, fast signals (AUC) | 1.00 | 0.90 | **1.00** | 1.00 |
| Weak coupling, 8× slower (AUC / delay found) | 0.35 / 0% | 0.67 / 30% | **0.91 / 95%** | 1.00 / 100% |
| Weak coupling, 16× slower (AUC / delay found) | 0.44 / 0% | 0.62 / 50% | **0.79 / 65%** | 0.92 / 85% |
| Signals of different speeds (worst case, AUC) | 0.76 | 0.50 | **1.00** | 1.00 |
| Coupling switching on/off (accuracy, slower signals) | 0.53 | 0.59 | **0.87** | – |
| 8-signal network (link AUC / delays found) | 0.88 / 60% | 1.00 / 100% | **1.00 / 100%** (v2-fixedtol) | – |
| Rulkov neurons, delay 32 (AUC / delay found) | 0.66 / 0% | 0 / 0% | **1.00 / 100%** | 1.00 / 100% |
| Rulkov neurons, delay 16 (AUC / delay found) | 1.00 / **20%** | 1.00 / 0% | 0.94 / 0% | 1.00 / 0% |

How to read these:
- **AUC:** 1 = coupled and uncoupled pairs perfectly separated; 0.5 = guessing.
- **Delay found:** the share of coupled pairs where the right delay was reported.
- **Rulkov at delay 16:** every adaptive rule reports 19–20 instead of 16. Even the oracle misses. The neuron bursts repeat every ~41 steps, so delays near half that period are hard to pin down.

## What carries over to "many science fields"

- **Sampling rate stops mattering.** The same sleep EEG at 0.5×, 1×, 2× and 4× the rate gives 75 / 58 / 41 / 27% with the defaults, but a stable answer with the adaptive rule. A lab recording at 256 Hz and one at 1 Hz would no longer get different conclusions from the same physiology.
- **Slow systems** (long memory, long delays) are where the defaults fail outright: the delay is longer than the fixed 30-sample search range, or the window is shorter than the signal's memory. Climate, ecology and slow physiology are typical examples. The adaptive rule handles them.
- **When the defaults already fit** (heart rate vs. breathing at 4 Hz), the adaptive rule picks almost the default values by itself (window 59 vs. 60). It does no harm.
- **Mixed systems** (fast and slow signals together): one system-wide parameter set based on the slowest signal worked best, better than per-pair parameters, and it keeps all scores on the same ruler.

## Limits found (honest list)

1. **Time resolution vs. reliability can't be decided by data.** v2 chose a 29-minute window on the 8-hour sleep EEG, so the score saturated at 100%. Capping the window at 2/5/10 minutes gave 79/90/94%, always with chance around 5%. **The longest window you accept should stay a user setting**, because it expresses what you want to resolve.
2. **Rhythmic waveforms** (raw ECG/PPG): delays longer than half a cycle can't be measured by any correlation method. The heart-to-finger delay (430 ms) wasn't recovered by any variant. Use beat-to-beat series instead of raw waveforms.
3. **Calibration can fail** on strongly rhythmic signals. It now says so instead of silently returning tolerance 0.
4. **Data:** the sandbox blocked most public data hosts (PhysioNet, Zenodo, NOAA, UCI), so the real-data tests cover sleep EEG and cardio-respiratory physiology only. Other fields were tested only with synthetic signals at different time scales.

## Side finding (separate from this study)

The docs claim the default reproduces the MATLAB 43.6% benchmark. Today it gives 58.1%; 43.5% needs `stability_anchor="end"`. I filed this as a separate suggested task.

## If you want to take it further

1. Move `adaptive_params(window_rule="bartlett", tolerance="calibrated")` into the library as an **opt-in** option, e.g. `TDSParams.auto(signals, max_window=...)`, keeping the current defaults unchanged.
2. Run the cross-subject surrogate test (`surrogate.py`) with the adaptive rule applied **inside each surrogate**, so network thresholds stay valid.
3. Test on more real fields once data access allows (PhysioNet sleep with hypnograms; climate indices such as ENSO vs. temperature).

## What was done, by whom

- **I wrote:** the adaptive module and all reports.
- **Coder agents:** built and ran the experiments (two in parallel, then one re-run, then the review fixes).
- **Reviewer agent:** audited the code.
- **Everything is on branch `claude/practical-mccarthy-n2mec6`, not merged to master.** The library code in `pyTDS/` is unchanged, and all 67 existing tests plus 3 new ones pass.
