# Apnea (PhysioNet ucddb 1.0.0) — read SPEC_COMMON.md first
Data: 25 subjects `ucddbNNN.rec` (EDF; header start time at bytes 176-184 'hh.mm.ss'; signals 'Flow' and 'SpO2' at 8 Hz) and
`ucddbNNN_respevt.txt` (events: time hh:mm:ss, type APNEA-*/HYP-*, duration s; parse lines matching ^\s*hh:mm:ss\s+TYPE ... duration).
Preprocessing (fixed): airflow band-pass 0.05-2 Hz (2nd-order Butterworth, zero phase); breathing envelope = sqrt(moving mean of
flow² over 4 s); SpO2 values < 50 or > 100 → NaN → linear interpolation; then subtract a 300-s moving average from BOTH signals
(remove slow drifts); keep 8 Hz. Pair: x = envelope, y = SpO2.
References (report both, per subject median):
 (a) event-locked trough lag: per event, SpO2-minimum time within [onset, end+60 s] minus envelope-minimum time within the event;
 (b) event-period cross-correlation lag: cross-correlation of the two signals restricted to event periods (events ±30 s), lag of max |C|
     within ±60 s.
Extra metric (coupling switches): event mask = events ±30 s; for each setting the % of stable windows (window centre) inside vs outside
event periods (coupled pairs), and for controls the same split. Note sign: low ventilation → low SpO2 later (positive correlation of
envelope and SpO2 at a lag where SpO2 follows).
Known probe on subject 002 (300-s detrend): reference (a) 24.2 s; aTDS L 359 s, lag ±180 s, tol 12.8 s, delay found ~5.8 s, stable
61.9 % in events vs 34.1 % outside; fixed TDS 8.4 % vs 1.7 %. Investigate and report why the TDS delay differs from reference (a)
(e.g. compare with reference (b)).
