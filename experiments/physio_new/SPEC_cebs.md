# ECG -> seismocardiogram (PhysioNet cebsdb 1.0.0) — read SPEC_COMMON.md first
Data: 20 subjects × b (basal, 5 min), m (music, 50 min), p (post-music); WFDB, 5 kHz; channels ECG I, ECG II, respiration, SCG;
.atr beat annotations. Use b and p records (5 min each) of all subjects (first 5 min of m optional).
Preprocessing (fixed): ECG II band-pass 0.5-30 Hz, SCG band-pass 1-30 Hz (or as in the README); downsample to 250 Hz.
Pair: x = ECG, y = SCG (electromechanical delay: the aortic-valve-opening wave of the SCG follows the R peak by ~50-120 ms).
Reference: per beat, time from R peak (annotation or detected) to the maximum of |SCG| (or the AO peak) within 20-200 ms after R;
median per record; compare delays modulo one beat.
Extra: basal vs post-music records; whether the rhythm cap applies (lag ≈ half RR).
