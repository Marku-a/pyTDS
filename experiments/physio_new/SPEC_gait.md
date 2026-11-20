# Gait (PhysioNet gaitpdb 1.0.0) — read SPEC_COMMON.md first
Data: text files GaCoXX_YY.txt (controls), GaPtXX_YY.txt / JuPt / SiPt (Parkinson's); columns: time, 8 left-foot sensors, 8 right-foot
sensors, total left force, total right force; 100 Hz. Use the first trial (_01) of every subject (~165). Exclude dual-task trials
(_10) if present.
Preprocessing: none beyond z-scoring (keep 100 Hz).
Pairs: (1) x = total left force, y = total right force (known: the right foot lands about half a stride after the left);
(2) x = left heel sensor (L1 or the most posterior sensor per the database description), y = left toe sensor (most anterior)
(within-foot roll-over delay).
Reference (exact, from the data): heel strikes = upward crossings of 20 % of max total force per foot; pair (1): median time from each
left heel strike to the next right heel strike, and median stride time; pair (2): median time from heel-sensor rise to toe-sensor rise
within each stance. Delays compared modulo the stride time.
Extra: controls vs Parkinson's comparison of delays found; whether aTDS's rhythm cap applies (lag ≈ half stride).
