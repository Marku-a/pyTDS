# Adaptive Time Delay Stability (aTDS)

## Problem Statement
Time Delay Stability (TDS; Bashan et al., Nature Communications 2012) is the standard coupling measure of network physiology: two signals are linked when the delay tau0 that maximises their windowed circular cross-correlation stays nearly constant over consecutive windows (a window is stable if at least 4 of 5 neighbouring delays agree within a tolerance; the TDS score is the percentage of stable windows). TDS depends on three parameters counted in samples - window L, lag range tau_max, tolerance tol - whose published values (60 / +-30 / +-1) were chosen for 1-Hz sleep recordings. For any other time scale or sampling rate they are a guess, and a fixed choice fails silently (e.g. a lag range narrower than the true delay still gives high TDS scores with a wrong delay). The authors of the open-source TDSpy implementation (Steinbrinker et al., SoftwareX 2024) note that differing parameter choices across implementations make TDS results hard to reproduce.

## Core Hypothesis
All three TDS parameters can be set from the signals themselves - from each signal's own autocorrelation and from coupling-free circularly rotated copies - without any coupling labels or known delays, so that TDS behaves consistently across time scales and sampling rates while keeping a fixed false-stability budget.

## Contributions
1. A window rule: L = 30 x max_i B_i, with B the Bartlett variance-inflation factor of each signal's autocorrelation, so every window holds ~30 independent samples of the slowest signal (floor 61 samples, cap 10 % of the record; optional longest plausible delay d_max).
2. A lag-range rule: tau_max = (L-1)/2 (the published ratio and the circular-correlation aliasing limit), capped at half the shortest detected rhythm period.
3. A tolerance rule: the widest tolerance whose chance score on circularly rotated copies, mean + 2 SE, stays <= 5 %.
4. A controlled evaluation against the published setting, a hand-rescaled setting told the true slow-down factor, and a labelled-data oracle, plus a per-link rotation test with Benjamini-Hochberg FDR control for networks.
5. Open-source implementation on top of pyTDS with all scripts, raw JSON results and fixed seeds.

## Scope (must be respected by the paper)
- Evidence is synthetic signals with known answers plus ONE real sleep-EEG pair (delta/sigma power, 1 Hz). No staged multi-subject sleep cohort; no cross-subject stage-matched surrogates; coupling that switches with physiological regime is untested.
- aTDS automates a time-scale-consistent parameter choice with a false-alarm budget; it is NOT optimal (an oracle tuned on labels does better) and the 5 % is a budget for the mean chance score, NOT a significance level.
- aTDS scores are not comparable with published TDS scores.
- Limitations listed in the experimental log must be reported, not hidden.

---

# Full method description (from the approved technical report)

## Adaptive TDS: Time Delay Stability that chooses its own parameters

TDS needs a window, a lag range and a tolerance, and nobody knows the right values for a new system. aTDS measures them from the signals. This brief tests each of the three rules against signals whose coupling and delay are known. It does not yet show whether aTDS changes or sharpens sleep-stage contrasts in real sleep networks; that is untested.

## The idea

Time Delay Stability (Bashan et al. 2012) links two signals when the delay τ₀ that best aligns them stays nearly constant over consecutive windows. A window counts as stable when at least 4 of 5 neighbouring delays agree within a tolerance. The TDS score is the percentage of stable windows.

The result depends on three numbers: the window, the lag range searched and the tolerance. The published values (60 / ±30 / ±1 samples) were chosen for 1‑Hz sleep recordings. They are counted in samples, so they mean something different for every system and every sampling rate. For a new system they are a guess.

aTDS: you give it the signals, and it sets all three numbers itself.

The 4-of-5 stability rule stays as published. The three rules use only each signal's own structure and coupling-free copies; they never use which pairs are coupled or what the delay is. (A coupled signal's own autocorrelation partly reflects its partner, so B is not entirely blind to coupling.)

Figure 1 · Same coupled pair, two parameter sets

Signal y follows x with a true delay of 20 samples (memory τc = 8 samples, coupling strength 0.6). Each dot is τ₀ in one window, over the first 12 000 samples. Filled dots are stable windows and hollow dots unstable ones. The dashed line is the true delay (τ = −20; negative means x leads). Top: fixed parameters (window 60, lag ±30, tolerance ±1). Bottom: aTDS, which chose window 259, lag ±129 and tolerance ±13 by itself, on a separate training pair of the same kind.

What it shows: the coupling is identical in both panels; only the parameters differ. With the fixed window, τ₀ scatters and rarely holds for 4 of 5 windows: score 5.2 %, with τ₀ within ±2 samples of the true delay in 28 % of windows. With the aTDS parameters, τ₀ sits on the true delay: score 91.3 %, within ±2 samples in 90 % of windows. Most of this gain comes from the window and lag range (Figure 5).

## How aTDS chooses the three parameters

## Window

L = 30 · B

B = 1 + 2 Σ r(k)², summed over the autocorrelation r(k) up to its first zero crossing (the Bartlett factor). For a smooth signal, B neighbouring samples carry about one sample's worth of new information. So L holds about 30 independent samples of the slowest signal in the system. Limits: at least 61 samples (so the lag range is never narrower than the published ±30), at most 10 % of the recording. Optional: if the user knows the longest plausible delay dmax, the window is at least 2·dmax + 1.

Why 30: with n independent samples each correlation value has a random error of about 1/√n (±0.18 at n = 30). 30 is a convention that balances reliability against time resolution; Figure 3b tests it.

Tested in Figures 2, 3 and 3b.

## Lag range and step

τmax = (L − 1) / 2

The same ratio as the published setting (±30 in a 60-sample window). Step = L / 2, the published 50 % overlap. If a signal has a rhythm with period P, τmax ≤ P / 2, because beyond half a cycle a delay d and d − P look the same.

Tested in Figure 4.

## Tolerance

largest tol with chance + 2 SE ≤ 5 %

One signal is circularly shifted by a large random amount. That keeps each signal's own structure but destroys the coupling. TDS is scored on these shifted copies at tolerance 0, 1, 2, …, and the largest tolerance whose mean chance score, plus two standard errors, stays ≤ 5 % is kept. Rotations are pooled over all recordings available (for a cohort, all subjects). If even 0 exceeds 5 %, the result is flagged calibration_failed.

Tested in Figures 5 and 6.

## The algorithm

Input: the signals of one system, s1 … sn, each T samples long. Output: one parameter set used for every pair in the system. The stability rule (4 of 5 windows) is not changed. Constants: 30 independent samples per window, 5 % chance budget, window between 61 samples and 10 % of T, 20 random shifts per pair.

## Step 0 · Autocorrelation of each signal

ri(k) = Σt xi(t) xi(t+k) / Σt xi(t)²,   xi = si − mean(si),   k = 0 … T/4

## Step 1 · Window

Bi = 1 + 2 Σk=1Ki−1 ri(k)²,   Ki = first k with ri(k) ≤ 0

L = clip( round(30 · maxi Bi), Lmin, 0.1 T ),   step = ⌊L / 2⌋

Lmin = 61,   or max(61, 2·dmax + 1) if the longest plausible delay dmax is given (optional)

The only input is the signals. dmax is optional. Without it, aTDS searches at least ±30 samples, the published range; with it, the lag range is guaranteed to reach dmax (Table 2, fast signal with delay 50).

What B means. Neighbouring samples of a smooth signal nearly repeat each other. Bi (the Bartlett factor) says how many consecutive samples of signal i carry as much new information as one independent sample: B = 1 for white noise, B = 8 means every 8 samples count as one. A window of L samples therefore holds L / Bi independent samples of signal i.

Both signals, one window. B is computed for every signal in the system (for a pair, for x and for y). The largest B, i.e. the slowest signal, sets the window, so every signal, and every pair built from them, has at least 30 independent samples in each window. One window is used for all pairs so that all scores are measured on the same ruler. When signals differ strongly in speed, this gives fast pairs a longer window and a wider tolerance than they need; per-pair windows are the alternative (Table 2).

Why 30. Each value of C(τ) in a window carries a random error of about 1/√n, where n = L / B is the number of independent samples (Bartlett's formula): ±0.18 at n = 30, ±1 at n = 1. But τ₀ is the largest of roughly 30 independent lags, and the largest of 30 pure-noise values is higher, about 0.18 · √(2 ln 30) ≈ 0.47. So a weak coupling (0.3) does not stand out in every single window. What separates it from noise is the 4-of-5 stability rule: noise peaks land at random lags, while the true peak keeps returning to the same lag. More independent samples make τ₀ more reliable but leave fewer windows, and so a coarser picture over time. 30 is a fixed design constant chosen as a balance, not an optimum: with on/off coupling, 45–60 did slightly better (Figure 3b).


## Step 2 · Lag range

τmax = ⌊(L − 1) / 2⌋

if some signal has a rhythm with period Pi:   τmax = min(τmax, ⌊mini Pi / 2⌋)

Rhythm test: after the first zero crossing of ri, take the highest autocorrelation peak; if it is at least 0.3 and not at the edge of the range, its lag is Pi. Otherwise the signal is treated as not rhythmic.


## Step 3 · Tolerance

1 · Make coupling-free copies. Rotate y in time by a large random amount h (the end wraps around to the start). Each signal keeps its own shape and rhythm, but y no longer follows x, so any stability found between them is pure chance. Repeat with 20 different h. When several recordings of the same kind are available (a cohort), pool the rotated copies of all of them.

2 · Score each copy at tolerance tol. For copy j, find τ₀ in every window (the same window and lag range as above), apply the 4-of-5 rule with this tolerance, and count the stable windows:

scorej(tol) = 100 % × (stable windows of copy j) / (all windows of copy j)

3 · Average over the copies.

C(tol) = meanj scorej(tol)

Because the copies contain no coupling, every stable window in them is a false alarm. C(tol) is therefore the false-alarm rate of TDS at that tolerance: the score you would get from uncoupled signals just by luck.

4 · Choose the tolerance. Try tol = 0, 1, 2, … and keep the largest tol whose upper bound, C(tol) plus two standard errors of that mean, is ≤ 5 %:

tol* = max { tol : C(tol) + 2 · SE(tol) ≤ 5 % },   SE = SDj(scorej) / √(number of copies)

The upper bound is used, not the mean alone, so that the budget also holds on new recordings, not only on the copies used to choose the tolerance. With the mean alone, the chance score on new pairs reached 6.7 %; with the upper bound and pooled rotations it stayed at 3.3–4.4 % (Figure 6).

Why a wider tolerance is wanted at all. A real delay is never perfectly constant; it wobbles from window to window. A wider tolerance lets TDS accept that wobble and see real coupling, but it also lets random τ₀ values agree by luck. The tolerance has to be as wide as possible without letting luck take over.

What the 5 % means. It is a design budget for false stability: on coupling-free copies, at most 5 % of windows on average may be labelled stable. It is not the probability that a link is false. Whether a particular link is significant still needs its own test (its score against its own rotated copies) and, in a network, a correction for testing many links at once. The 5 % is a fixed design constant; a stricter user could use 1 %.

Why stopping at the first tol above 5 % is right. A wider tolerance can only add stable windows, never remove them, so C(tol) only goes up. Once it passes 5 %, every larger tolerance is also above 5 %, and the previous one is the widest tolerance inside the budget. (The upper bound behaves the same way in practice.)


aTDS_parameters(s1 … sn):
    # Step 1: window
    for each signal i:
        ri = autocorrelation(si)
        Ki = first k where ri(k) ≤ 0
        Bi = 1 + 2 · sum(ri(k)² for k = 1 … Ki−1)
    Lmin = 61   (or max(61, 2·d_max + 1) if the optional d_max is given)
    L    = clip(round(30 · max(B)), Lmin, 0.1·T)   # Lmin wins if 0.1·T < Lmin
    step = L // 2

    # Step 2: lag range
    lag = (L − 1) // 2
    P   = [dominant_period(ri) for each signal, if rhythmic]
    if P not empty: lag = max(1, min(lag, min(P) // 2))

    # Step 3: tolerance
    null = []
    for each pair (i, j) in the system, in every recording available (at most 50 pairs):
        repeat 20 times:
            h = random integer in [a, T − a),  a = max(L, T // 10)
            null.append( τ₀ series of (si, rotate(sj, h)) with L, step, lag )
    tol = 0
    for t = 0, 1, 2, … up to max(1, lag // 4):
        s  = [ TDS score of each null series, 4-of-5 rule, tolerance t ]
        C  = mean(s);  SE = sd(s) / √len(s)
        if C + 2·SE ≤ 5 %: tol = t
        else: break
    calibration_failed = (C + 2·SE at t = 0) > 5 %

    return window = L, step, max_lag = lag, tolerance = tol
