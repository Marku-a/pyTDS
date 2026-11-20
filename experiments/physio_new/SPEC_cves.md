# Cerebral vasoregulation (PhysioNet cves 1.0.0) — read SPEC_COMMON.md first
Data: older adults with/without stroke; inspect RECORDS, files-and-variables.xlsx (read with pandas/openpyxl if available, else
unzip XML) and data/ to find per-subject recordings with arterial blood pressure (ABP), cerebral blood flow velocity (CBFV, MCA,
transcranial Doppler) and CO2 (end-tidal / capnography). Use whatever protocol segment is available for most subjects (prefer one
containing CO2 changes such as hyperventilation/rebreathing; else rest/sit-to-stand). Cap at 40 subjects.
Preprocessing (fixed): low-pass (0.5 Hz) or beat-average ABP and CBFV to remove the cardiac pulse, resample all to 4 Hz, subtract a
120-s moving average.
Pairs: (1) x = CO2, y = CBFV (known: CBFV follows CO2 with a delay of several seconds, ~3-10 s); (2) x = ABP (mean), y = CBFV
(autoregulation: CBFV follows ABP within ~0-3 s).
Reference delay: if protocol annotations give the onset of a CO2 manoeuvre, the step-response delay (time from CO2 change onset to
63 % of the CBFV change); otherwise the full-record cross-correlation lag on first differences, stated as such. Report n subjects used.
