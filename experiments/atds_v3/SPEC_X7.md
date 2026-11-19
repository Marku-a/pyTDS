# X7 — validate the Bayley-Hammersley window estimator outside the benchmark it was found on

X6 item 4 found that replacing the first-zero Bartlett factor by a Bayley-Hammersley-type factor
`B_BH = max(1, 1 + 2 * sum_{k=1}^{K} (1 - k/n) r(k))`, K = first lag with r(k) <= 0 (estimator (c) in `x6_followups.py`),
raises aTDS's AUC on the slow-down benchmark. Because it was chosen on that benchmark, test it **unchanged** elsewhere.
Define `atds_bh`: identical to aTDS (30 independent samples, floor 61, cap 10 %, lag (L-1)//2 with rhythm cap, tolerance calibration unchanged) except the window uses B_BH (max over signals).

Same rules as SPEC.md (pyTDS only, pyTDS unchanged, fixed seeds, 4 workers, aggregates only). Reuse `x6_followups.py`.

1. **Held-out slow-down set**: the U1 design (tc = 2r, d = 5r, c = 0.3, r in {1, 2, 4, 8, 16}) with NEW seeds (offset all seeds by +5 000 000), 100 + 100 pairs: published, rescaled_est, atds, atds_bh — AUC (CI), delay found, chance.
2. **Other generators** (seeds as in their original experiments): long memory fGn H = 0.6 / 0.75 / 0.9 (exp4 U5), rhythms d = 8 / 30 (U6), mixed system and fast-slow pair (U9), changing dynamics (U10), Rulkov (exp5 V1): atds vs atds_bh vs rescaled_est with the original metrics.
3. **Time-resolved cases**: E3 regime switching (with and without the X6 boundary margin), drifting delay: atds_bh next to the X6 numbers.
4. **Per-link power** as X6 item 5 (r in {2, 4, 8}, 50 + 50 pairs, 99 rotations) for atds_bh.
5. **Real records**: the five records of X6 item 1 (200 rotations) and the sleep cohort `atds_cap300` analysis repeated with the BH window (cap 300 still applies, so report whether the window changes at all).

Report a single comparison table per item with atds, atds_bh and rescaled_est (and published), and state plainly where atds_bh is better, equal or worse than atds and than rescaled_est.
