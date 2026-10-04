# Instructions: how adaptive TDS works and how to use it

*Published 2026-10-04 01:14 UTC. Code: `experiments/adaptive/adaptive.py` (experimental, not in the `pyTDS` library yet).*

## 1. What TDS does, in one minute

TDS (Time Delay Stability) asks: **do two signals keep a fixed time delay between them?** If they do, they are probably coupled, meaning one drives the other.

1. Cut both signals into short pieces called **windows**. Each window starts a **step** later than the previous one.
2. In each window, slide one signal against the other, up to a **max lag**, and find the shift where they match best. That shift is the window's **delay τ₀**.
3. Look at 5 consecutive windows. If at least 4 delays agree within a **tolerance**, those windows are **stable**.
4. **TDS score** = % of windows that are stable. A high score means a stable delay, which means coupling.

## 2. Why the parameters must adapt

The original parameters are fixed sample counts: window 60, step 30, max lag 30, tolerance ±1. They fit 1-sample-per-second sleep data, but in other settings they break:

| your data | what goes wrong with the fixed values |
|---|---|
| sampled faster (100 Hz heart signals) | the window is 0.6 s, shorter than a heartbeat |
| slow signals (climate, slow physiology) | the window is shorter than the signal's own memory, so delays look jittery |
| long delays (> 30 samples) | the delay is outside the search range and can never be found |
| the same data at a different sampling rate | a different score for the same physics (58% → 27% in our EEG test) |

## 3. How the adaptive rule sets each parameter

It computes **one parameter set for the whole system** (all the signals you analyse together), so all scores stay comparable.

**Window: "make every window hold ~30 independent samples".**
Neighbouring samples of a smooth signal are near-copies, so they don't count as new information. The **Bartlett factor** B says how many samples make one independent sample:

- B = 1 + 2·Σ r(k)², where r(k) is the signal's correlation with itself shifted by k. The sum runs until r first drops to 0.
- A white-noise signal has B ≈ 1; a slow, smooth signal has a large B.
- **window = 30 × B of the slowest signal in the system.**

**Step** = window / 2, the same 50% overlap as the original.

**Max lag** = (window − 1) / 2, the same ratio as the original. If a signal is **rhythmic** (period P), the max lag is capped at P/2. Beyond half a cycle, a delay of d and d − P look identical.

**Tolerance: "allow at most 5% fake stability".**
1. Make coupling-free copies of the data: rotate one signal of each pair in time by a large random amount. Its own rhythm is kept; the coupling is gone.
2. Run TDS on these copies with tolerance 0, 1, 2, …
3. Keep the **largest tolerance whose chance score is ≤ 5%**.

Slow signals get a wider tolerance; fast signals a narrower one.

**Safety rule:** nothing in steps 1–4 looks at the real coupling. The rule cannot be "tuned to find what you hope to find".

## 4. Fully automatic: you only give it the signals

You don't have to choose anything. `adaptive_params([s1, s2])` sets the window, the step, the max delay and the tolerance by itself. A built-in safety limit keeps the window at most 10% of the recording, so there are always enough windows. All nine field results (Report 7) were produced this way, except that the sleep EEG there was capped at 5 min; uncapped, it chooses 29 min and gives 100% against a chance level of 24%.

**Optional override, `max_window`:** use it only if you need to know *when* the coupling happens at a finer time scale than the automatic window. For example, on 8 h of sleep EEG the automatic window is 29 min. That is fine for "are these signals coupled?", but too coarse for "in which sleep stage?". Capping it at 2/5/10 min gave 79/90/94%, with chance still about 5%. If you don't have that question, ignore this option.

## 5. How to run it

```python
import sys; sys.path.insert(0, "experiments/adaptive")
from adaptive import adaptive_params
from pyTDS.core import tds

# all signals of your system (list of 1-D arrays, or a T x N array)
params, info = adaptive_params(
    [s1, s2, s3],
    window_rule="bartlett",     # the v2 rule (recommended)
    tolerance="calibrated",     # 5% fake-stability budget
    # max_window=300,           # OPTIONAL: only for finer timing (samples)
)
result = tds(s1, s2, params)   # use the SAME params for every pair
print(params)                   # window, overlap(=step), max_lag, tolerance
print(result["score"])
```

`info` explains every choice:

| key | meaning |
|---|---|
| `bartlett` | B per signal (samples per independent sample) |
| `tau_c` | decorrelation time per signal: the lag where self-correlation drops below 0.37 |
| `periods` | detected rhythm period per signal (None = not rhythmic) |
| `window`, `step`, `max_lag`, `tolerance` | the chosen parameters |
| `window_clipped` | True if your `max_window` (or the 10%-of-recording limit) cut the window |
| `null_score_by_tolerance` | chance score for each tolerance tried |
| `calibration_failed` | **True = warning.** Even tolerance 0 gives > 5% fake stability, so the result is not trustworthy as is. Raise `max_window`, or use a different input (see below). |

## 6. When it will not help

- **Raw rhythmic waveforms** (ECG, PPG, any oscillation). Correlation can only see delays shorter than half a cycle. Convert to beat-by-beat or cycle-by-cycle series first.
- **Very short recordings.** With fewer than ~20 windows the score is coarse. The window is never allowed above 10% of the recording unless you set `max_window`.
- **Comparing with published TDS values.** Adaptive scores are on a different (fairer) scale than the original fixed-parameter scores, so compare like with like.
