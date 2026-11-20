# Treadmill exercise (PhysioNet treadmill-exercise-cardioresp 1.0.1) — read SPEC_COMMON.md first
Data: `test_measure.csv` (breath-by-breath rows per test: inspect columns — expect time, speed, HR, VO2, VCO2, RR, VE, ID) and
`subject-info.csv`. ~992 tests. Use tests with ≥ 8 minutes of data and complete HR/VO2/VCO2; take a random (seed 0) sample of 200 tests
to keep runtime bounded (report the selection).
Preprocessing (fixed): per test, interpolate each variable to an even 1-Hz grid; remove the slow exercise ramp by subtracting a
60-s moving average from every signal (so only breath-to-breath and minute-scale fluctuations remain) — report this choice and also
run a variant WITHOUT detrending for aTDS and fixed TDS only.
Pairs (run each): (1) x = VO2, y = VCO2 (physiology: VCO2 lags VO2 by ~10-20 s because of CO2 storage);
(2) x = speed, y = HR (HR lags workload changes).
Reference delays (per test, independent of TDS): for pair (1) and (2), the lag (s) of the maximum of the cross-correlation computed on
first differences of the 1-Hz series within ±60 s (state clearly that this is a cross-correlation reference, not an annotation); for
pair (2) also, if speed has steps (changes > 0.3 km/h within 5 s), the step-response delay: time from each speed step to the HR reaching
63 % of its subsequent change (median over steps), else report n/a.
Extra metrics: distribution of delays found vs references; whether fixed TDS (±30 s) hits the edge of its lag range.
