# Pulse transit time during exercise (PhysioNet pulse-transit-time-ppg 1.1.0) — read SPEC_COMMON.md first
Data: 22 subjects × activities sit / walk / run (sNN_sit, sNN_walk, sNN_run; WFDB .hea/.dat, 500 Hz; also csv/ folder). Inspect the
header/README for channels (ECG, several PPG channels, accelerometers); use ECG and one finger PPG channel (the one closest to the
fingertip, same for all).
Preprocessing (fixed): ECG band-pass 0.5-30 Hz, PPG 0.5-8 Hz; downsample both to 125 Hz.
Pair: x = ECG, y = PPG.
Reference: beat-by-beat pulse arrival time = R peak to steepest PPG upstroke in the following beat window (median per record) and the
median RR; compare delays modulo one beat (as experiments/atds_real/bidmc.py does). Also the full-record cross-correlation within half
a beat.
Extra (time-varying across activity): for each setting, the delay found per activity vs the reference per activity; does the found
delay shorten from sit to run as the reference does? Controls: ECG of subject i with PPG of subject j, same activity.
Run aTDS common set separately per activity and also one common set over all activities.
