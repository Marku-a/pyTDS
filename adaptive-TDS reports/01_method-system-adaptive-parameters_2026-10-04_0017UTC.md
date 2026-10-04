# Report 1: The method, system-adaptive TDS parameters

*Published 2026-10-04 00:17 UTC · branch `claude/practical-mccarthy-n2mec6` (not merged to master)*

## The problem in one paragraph

TDS has four parameters, and all of them are counted in **samples**:

- **window** = 60: how much data each cross-correlation uses.
- **step** = 30: how far the window moves each time. In the code this is called `overlap`.
- **max_lag** = 30: the largest delay it looks for.
- **tolerance** = ±1: how much the delay may wobble and still count as "the same delay".

These numbers were chosen for 1 Hz EEG band-power data (Bashan et al. 2012). For other data the same numbers mean something completely different:

- **100 Hz heart signals:** 60 samples is 0.6 s, which is shorter than one heartbeat.
- **Slowly changing signals:** ±1 sample is far stricter than the precision a delay can actually be measured with.

So the question is: **can the parameters be derived from the data automatically, so TDS works in any field without hand-tuning?**

## What "system-adaptive" means here (your choice)

The algorithm computes **one parameter set for the whole system** (all the signals you analyse together), not a separate set per signal pair. This keeps every score in the network measured with the same ruler, so scores can be compared with each other.

## How the parameters are derived

Code: `experiments/adaptive/adaptive.py`. The library in `pyTDS/` is **not changed**.

1. **Measure each signal's time scale.** The *decorrelation time* τc is how many samples it takes until the signal's correlation with a shifted copy of itself falls below 1/e (≈ 0.37). Roughly, "how long the signal remembers its own past".
2. **Choose the system time scale.** By default this is the slowest signal's τc, so every pair gets enough independent data per window. A variant uses the median signal instead.
3. **window** = 10 × the system time scale, clipped to between 20 samples and 10% of the recording. **step** = window / 2, the same 50% overlap as the default. **max_lag** = (window − 1) / 2, the same ratio as the default.
4. **Oscillation guard.** If a signal oscillates with period P (it has a clear repeat in its autocorrelation), max_lag is capped at P/2. Beyond that, the delay is ambiguous: a delay of d and a delay of d − P look identical.
5. **Tolerance**, in two versions:
   - **fixed:** keep ±1, as in the original.
   - **calibrated:** build a "no coupling" version of the data by rotating one signal of each pair in time by a large random amount (a *circular shift*). This destroys any real coupling but keeps each signal's own rhythm. Run TDS on it many times. Then pick the **largest tolerance whose chance score stays at or below 5%**.

    So the tolerance is set by a false-positive budget, not by hand. Slow signals automatically get a wider tolerance; fast signals get a narrow one.

**Important safeguard:** none of these steps looks at the real coupling. Steps 1–4 use each signal alone, and step 5 uses only coupling-free shifted data. The parameters therefore cannot be "tuned to find coupling", which would inflate the scores.

## How it is tested (details in the following reports)

| Data | Field | What is known ("tagged") |
|---|---|---|
| Synthetic linear signals (AR(1) processes) | generic | exact delay, coupled or not, when coupling is on or off |
| Rulkov chaotic neurons (from `examples/`) | neuroscience model | exact delay, coupled or not |
| `data/` EEG delta/sigma | sleep EEG | MATLAB reference score 43.6% |
| NeuroKit2 ECG + PPG, 100 Hz | cardiology | the heart-to-finger pulse delay, measured independently from the beats |
| NeuroKit2 heart rate + breathing, 3 recordings | cardio-respiratory physiology | breathing is known to modulate heart rate (respiratory sinus arrhythmia) |

Each method is compared with: **fixed defaults**, **adaptive (tolerance fixed)**, **adaptive (tolerance calibrated)**, and in the synthetic tests an **oracle**. The oracle is the best parameter set found by trying a grid while looking at the true answers. It can't be used in practice and only shows the best that is achievable.

Measures used throughout:
- **AUC:** the chance that a randomly chosen coupled pair scores higher than an uncoupled one. 1.0 is perfect separation; 0.5 is a coin flip.
- **Hit rate:** the fraction of coupled pairs where the median stable delay equals the true delay (within 10%).

**Data access note:** most public data hosts (PhysioNet, Zenodo, NOAA, UCI) are blocked by this sandbox's network policy. Only GitHub-hosted data could be downloaded, which is why the external data is the NeuroKit2 sample recordings.
