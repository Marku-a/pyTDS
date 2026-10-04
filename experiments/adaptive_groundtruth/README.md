# Does adaptive TDS work? A ground-truth test on public multi-subject data

Question: on Garmin runs the adaptive rule (v2-calibrated, from branch
`claude/practical-mccarthy-n2mec6`) found more speed → HR links than the fixed
setting. Earlier tests could only show that adaptive *scores* are higher, and a
higher score is not automatically better. This study asks the sharper question:
**does adaptive TDS reproduce known effects better, across many people, and what part
of it does the work?**

Design and hypotheses were committed before any full run (`SPEC.md`, commit f533372).
All data come from PhysioNet's open-data mirror on AWS (`physionet-open.s3.amazonaws.com`);
raw files are not committed.

## Data and ground truth

| dataset | people | signals | known effect (source) |
|---|---|---|---|
| Fantasia | 20 young, 20 elderly, 110 min supine rest | breathing → heart rate, 4 Hz | cardio-respiratory coupling weakens with age |
| Sleep-EDF Expanded | 20 nights (18–19 usable per contrast) | 10 EEG band powers (2 channels × δ θ α σ β), 1 Hz | TDS links weakest in deep sleep (DS); stronger in light sleep (LS) and wake (W); REM intermediate (Bashan et al. 2012 Nat Commun; Lin et al. 2020 Commun Biol) |
| Treadmill CPET | 197 graded exercise tests | treadmill speed → HR, VO2, 0.5 Hz | speed is imposed, so HR must follow it |

## Variants (an ablation of what v2-calibrated changes)

| variant | window | max lag | tolerance |
|---|---|---|---|
| published | 60 s | 30 s | ±1 s (Bashan 2012 setting) |
| default_samples | pyTDS defaults in samples (60 / 30 / ±1) | | |
| win_only | v2 window (Bartlett rule, capped) | window / 2 | ±1 s |
| tol_only | 60 s | 30 s | calibrated on a circular-shift null (5 % chance) |
| **v2_calibrated** | v2 window | window / 2 | calibrated |
| v2_cal_lag3 | v2 window | window / 3 | calibrated |

Caps fixed in advance: 300 s (Fantasia, sleep), 240 s (treadmill, as in GarminTDS).
The v2 rule hit the cap on almost every subject in every dataset.

## Results

Effect sizes with bootstrap 95 % CIs (5,000 resamples of subjects; `bootstrap_effects.py`,
`results/effects_bootstrap.csv`). "Gain" = difference from the published setting, paired.

| known effect | published | win_only | tol_only | **v2_calibrated** | v2_cal_lag3 |
|---|---|---|---|---|---|
| Sleep LS > DS (dz) | 2.67 | 3.64 | 2.78 | **3.97** (gain +1.30, CI +0.52…+2.50) | 4.04 |
| Sleep REM > DS (dz) | 0.41 (CI −0.03…0.92: not detected) | 1.00 | 0.49 | **1.13** (gain +0.72, CI +0.31…+1.49) | 1.15 |
| Sleep W > DS (dz) | **1.40** | 1.14 | 1.47 | 1.18 (gain −0.23, CI −1.07…+0.15) | 1.19 |
| Fantasia young > elderly (d) | 1.25 | 0.73 | 1.28 | **1.54** (gain +0.29, CI −0.37…+1.07) | 0.86 |

Other checks:

- **Fantasia:** all variants detect breathing → HR coupling in 85–100 % of people (lag 0.5 s), so
  detection is at ceiling. False positives on cross-subject pairs (breathing of one
  person vs HR of another): published 4/40, v2_calibrated 1/40.
- **Treadmill: uninformative.** The graded protocol is flat → ramp → drop; the 120-s
  high-pass removes the flat parts and the ramp, leaving only the end-of-ramp drop. Every
  variant detects speed → HR in 1.5–7 % of tests, about the false-positive level. This was
  a flaw in the pre-registered design (found by review), so no conclusion is drawn from it.
  One weak observation: among the few significant tests, the published setting put HR
  *before* speed (impossible) in 11/16, v2_calibrated in 2/13.
- **Selection control** (`sleep_selection_control.py`): 300-s windows only count inside
  long stage episodes. Scoring the 60-s setting only inside those same 300-s pure stretches
  raises LS > DS from 2.56 to 2.78 and REM > DS from 0.40 to 0.50, far short of v2's 3.97 and
  1.13. The gain is not window selection.

## What the adaptive rule is actually doing

**It gives each window enough independent data.** Slow physiological signals have long
memory. Measured with the Bartlett factor of the slowest signal (`independent_samples.py`):

| dataset | independent samples in a 60-s window | in the v2 window |
|---|---|---|
| Fantasia (HR) | 1.8 | 8.9 |
| Sleep EEG band power | 0.2 | 0.8 |
| Treadmill (high-passed) | 3.4 | 13.5 |

With 1–2 independent samples per window, the delay found in each window is mostly noise,
so "stable" windows are partly luck and real state differences are blurred. A longer
window fixes that. The ablation agrees:

1. **The window is the main ingredient.** win_only gets most of the sleep gain (LS 3.64,
   REM 1.00). Calibrating the tolerance alone at 60 s changes almost nothing (2.78, 0.49).
2. **The calibrated tolerance is what makes the long window safe.** On Fantasia the long
   window with ±1 s tolerance is the *worst* variant (d 0.73): a delay that wobbles by a
   sample or two in a 300-s window breaks stability. Calibration widens the tolerance to
   what chance allows (±15 s here) and recovers the best result (1.54).
3. **The wide lag search matters on some data.** Limiting max lag to window / 3 costs
   nothing in sleep but drops Fantasia to d 0.86, as it did on the Garmin runs
   (8/30 → 4/30 runs). With a loose tolerance, a narrow lag range makes chance
   stability more likely and raises the bar for real coupling.
4. **The price is time resolution.** v2 windows are 5× longer, so there are 5× fewer
   windows. Short states suffer: wake bouts are brief and fragmented, and W > DS is the
   one contrast where v2 is weaker (not significant). Short recordings suffer too
   (treadmill tests give about 8 windows).

## Verdict

- **Adaptive TDS works for what it is designed for: slow signals measured in long recordings.**
  It reproduces the published sleep-stage ordering of TDS networks better than the published
  fixed setting. It detects the REM vs deep-sleep difference, which the fixed setting misses,
  and it gives the largest young vs elderly difference.
- **It is not magic.** The gain comes from one thing: windows long enough to hold real
  information. The calibrated tolerance and wide lag search are what keep those long
  windows from failing.
- **It does not help for short records or short events.** Rapid state changes and
  recordings with fewer than about 20 windows are its weak spot, which applies to
  25-minute runs.

## Reproduce

```
pip install wfdb mne pandas scipy
python experiments/adaptive_groundtruth/fantasia.py --jobs 4      # ~10 min, 306 MB download
python experiments/adaptive_groundtruth/sleep_edf.py --jobs 4     # ~20 min, ~1 GB download (deleted after caching)
python experiments/adaptive_groundtruth/treadmill.py --jobs 4     # ~5 min, 23 MB
python experiments/adaptive_groundtruth/bootstrap_effects.py
python experiments/adaptive_groundtruth/sleep_selection_control.py
python experiments/adaptive_groundtruth/independent_samples.py
```
