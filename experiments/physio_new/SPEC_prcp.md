# Posture change and baroreflex (PhysioNet prcp 1.0.0) — read SPEC_COMMON.md first
Data: 10 subjects, records NNNNN (.hea/.dat WFDB; .wabp blood-pressure beat annotations, .wqrs QRS annotations, .anI posture/
event annotations). Write a small WFDB reader (or pip install wfdb from PyPI). Inspect signals (ECG, ABP).
Preprocessing (fixed): beat series: RR intervals from QRS; systolic BP per beat (max ABP in each beat from .wabp onsets); interpolate
both to 4 Hz; subtract a 60-s moving average. Pair: x = SBP, y = RR (baroreflex: RR lengthens after SBP rises; known delay 0-2 beats,
~0-2 s).
Reference: beat-domain cross-correlation (SBP_k vs RR_{k+j}, j = 0..5 beats) → lag of max correlation in beats × median RR (s).
Extra: if .anI marks posture changes, report TDS stable-window fraction before/after standing for each setting.
Note n = 10 → controls k = 3; AUC CI will be wide; say so.
