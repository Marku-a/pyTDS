# Fantasia (PhysioNet fantasia 1.0.0) — read SPEC_COMMON.md first
Data: 40 records f1y01..f1y10, f2y01..f2y10 (young), f1o01..f1o10, f2o01..f2o10 (elderly); WFDB format (.hea/.dat, 250 Hz; read
the header to find the respiration and ECG channels and the storage format; write a small reader with numpy — no wfdb package
installed; pip install of wfdb from PyPI is allowed if simpler). Beat annotations (.ecg) exist; you may use them or detect R peaks
(scipy find_peaks on band-passed ECG).
Preprocessing (fixed): instantaneous heart rate from R-R intervals, evenly resampled at 4 Hz (linear interpolation of RR at beat
times), respiration band-passed 0.05-1 Hz and resampled to 4 Hz; both z-scored. Use the first 60 minutes of each record. Pair:
x = respiration, y = heart rate.
Reference delay (independent of TDS, per subject): breath-locked phase analysis: detect inspiration onsets (troughs → rising) of
the respiration signal; average the heart-rate response over breaths aligned at inspiration onset (window 0..1.5 × median breath
period); reference = latency of the HR maximum after inspiration onset (s) and the median breath period. Report both.
Extra metrics: young vs elderly comparison of TDS score and n above controls (RSA known weaker in elderly); whether the rhythm cap was
applied by aTDS (lag ≈ half breath period) and what rescaled fixed TDS's lag range is relative to the breath period; delay found
compared modulo one breath period.
