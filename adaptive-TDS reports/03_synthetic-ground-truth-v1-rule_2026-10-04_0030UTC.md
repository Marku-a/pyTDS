# Report 3: Ground-truth tests of the first adaptive rule (v1)

*Published 2026-10-04 00:30 UTC. Code: `experiments/adaptive/exp_synthetic.py`, `generators.py`. Numbers: `experiments/adaptive/results/synthetic_*.csv`. Figures: `figures/E*.png`.*

> **Note added 00:53 UTC:** a later code review corrected some numbers in this report (mainly the network test and the Rulkov tests). The conclusions did not change. The corrected numbers are in Report 5.


## Why synthetic data

In synthetic data I **build the coupling myself**, so I know:
- whether two signals are coupled,
- the exact delay,
- when the coupling is switched on or off.

That makes it possible to say which method is *more accurate*, which real data (Report 2) cannot.

## The test signals

- **AR(1) signals:** random signals with a controllable "memory", or time scale τc. The signal at each step keeps a fraction of its previous value and adds new noise. Small τc gives fast, jittery signals; large τc gives slow, smooth ones.
- **Coupling:** signal y = c × (signal x delayed by d samples) + independent noise of its own. The coupling strength c is 0 (none), 0.3 (weak) or 0.6 (strong).
- **Rulkov neurons:** a standard model of a bursting neuron. Two of them are coupled with a known delay; this reuses `examples/rulkov_validation.py`.

Every condition is repeated 20 times with different random seeds: 20 coupled and 20 uncoupled pairs. The adaptive parameters are computed **once per condition** from a pool of coupled and uncoupled signals, so they cannot know which pairs are coupled.

## Scores used

- **AUC:** how well the TDS score separates coupled from uncoupled pairs (1 = perfect, 0.5 = guessing).
- **Hit rate:** the fraction of coupled pairs where TDS reports the right delay (within 10%).
- **Oracle:** the best parameter set from a grid, chosen **while looking at the answers**. It can't be used in practice; it only shows what is achievable.

## Results

**E1: same coupling, growing time scale.** Delay = 5r and τc = 2r, where r is the scale factor (1, 2, 4, 8, 16). Each cell is AUC / hit rate.

| r (scale) | default | adaptive v1 (calibrated) | oracle |
|---|---|---|---|
| 1, weak coupling | **1.00 / 1.00** | 0.90 / 0.85 | 1.00 / 1.00 |
| 4, weak | 0.60 / 0.10 | **0.78 / 0.60** | 1.00 |
| 16, weak | 0.44 / 0.00 | **0.53 / 0.40** | 0.93 / 0.85 |
| 8–16, strong | 0.52–0.60 / 0.00 | **1.00** | 1.00 |

- **Slow signals:** the default breaks down. The delay grows past its fixed 30-sample maximum lag, so it can never be found (hit = 0). Adaptive keeps working with strong coupling.
- **Fast signals:** v1 is *worse* than the default. Its window was only 20 samples, against the default's 60.

**E1c: signals with different speeds** (AUC):

| (τc of x, τc of y) | default | adaptive v1 |
|---|---|---|
| fast, fast | **1.00** | 0.78 |
| fast, slow | 1.00 | 1.00 |
| slow, fast | 0.76 | **1.00** |
| slow, slow | 0.96 | **1.00** |

**E1d: coupling switched on and off in blocks.** This measures whether TDS marks the right moments as coupled (balanced accuracy, 0.5 = guessing):

| scale r | default | adaptive v1 |
|---|---|---|
| 1 | **0.89** | 0.56 |
| 4 | 0.53 | 0.61 |

**E1b: a network of 8 signals** (4 fast, 4 slow, 5 true links). This is your "whole system" setting. Edge AUC measures how well true links are ranked above non-links.

| variant | edge AUC | delay hit |
|---|---|---|
| default | 0.87 | 0.60 |
| **one adaptive set for the whole system** (window from slowest signal) | **0.98** | **1.00** |
| same, window from the median signal | 0.92 | 0.82 |
| separate set per pair | 0.96 | 1.00 |
| one set per time-scale group | 0.96 | 1.00 |

**E2: Rulkov neurons, delay 16.** The default finds the delay (hit 0.9). Adaptive v1 fails: AUC 0.64, hit 0. Uncoupled pairs scored 15–37% because the window (~40 samples) was too short for the neurons' slow bursting.

## What this tells us

1. **The whole-system idea works.** On the mixed network, one adaptive parameter set beat the defaults (AUC 0.98 vs. 0.87) and found every delay. It was also as good as or better than per-pair parameters, while keeping scores comparable.
2. **The v1 window rule is too short.** "10 × decorrelation time" gives tiny windows for fast or spiky signals. It loses to the default there and fails on the Rulkov neurons.
3. **The oracle always preferred long windows with tolerance 0.** That suggests longer windows are better when coupling is steady. When coupling turns on and off (E1d), shorter windows are needed to see when it happens. This is the one trade-off data cannot settle for you: **how short a coupling episode you care about**.

## What I changed next (v2)

The new window rule is based on **effective sample size** (Bartlett's formula). It sums a signal's whole autocorrelation, so it sees slow structure such as the bursts in Rulkov neurons, and it sizes the window to hold about 30 independent samples. For simple signals this triples the v1 window: fast signals get ≈ 64 samples, close to the default 60.

A quick check on Rulkov:
- **v1:** uncoupled 31%, coupled 18%. The rule failed: the pair without coupling scored higher.
- **v2:** uncoupled 4%, coupled 38%, delay 19 vs. true 16.

Report 4 covers the full re-run.
