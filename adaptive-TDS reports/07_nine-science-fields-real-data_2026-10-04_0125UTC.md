# Report 7: Adaptive TDS on real data from nine science fields

*Published 2026-10-04 01:25 UTC. Code: `experiments/adaptive/exp_fields.py`. Numbers: `experiments/adaptive/results/fields.json` / `fields_summary.csv`. Data sources: `experiments/adaptive/data_fields/MANIFEST.json` (the raw files are not committed).*

## What I did

1. **Found the data.** An agent searched the internet for public time series in which one signal is *known* to drive another. Most data portals are blocked in this sandbox (NOAA, PhysioNet, USGS, Kaggle, EIA…), so it used GitHub-hosted copies. Two of them are personal mirrors of official data (the El Niño index from NOAA, and petrol prices from EIA); the rest come from NIST, CRAN package data, TensorFlow's dataset store and the "datasets" org.
2. **Ran both versions on each pair:**
   - **Fixed:** window 60, step 30, max delay 30, tolerance ±1 sample.
   - **Adaptive:** the v2 rule from Reports 4 and 6. No time-resolution cap, except for sleep EEG (capped at 5 min).

   The preprocessing was the same for both versions.
3. **Built a "chance" baseline for every result.** The second signal is rotated in time 200 times, which destroys the coupling, and TDS is re-run each time. **chance** below is the score that only 5% of those coupling-free runs exceed. A result counts only if it is well above chance.

## Results

| field: leader → follower | n, interval | fixed: score / chance / delay | adaptive: score / chance / delay | verdict |
|---|---|---|---|---|
| Hydrology: rain → river flow | 43,848 × 1 h | 1.2% / 1.1% / −9 h | **82% / 14% / −10 h ✓** | ✓ **better**: fixed can't tell it from chance |
| Hydrology: same catchment, daily | 10,593 × 1 day | 59% / 2.6% / −1 d ✓ | 98% / 20% / −2 d ✓ | ≈ both fine |
| Energy: heat → electricity demand (summer) | 13,008 × 30 min | 12% / 3% / 0 ✗ | **100% / 26% / −1 h ✓** | ✓ **better**: right delay |
| Air quality: wind → PM2.5 | 43,824 × 1 h | 5.4% / 1.6% / 0 | **57% / 12% / 0** | ✓ stronger detection |
| Climate: El Niño → global temperature | 1,595 × 1 month | 37% / 8% / **−3 mo ✓** | 68% / 26% / −1 mo ✗ | ≈ **mixed** |
| Engineering: gas feed → CO₂ (furnace) | **296** × 9 s | **100% / 0% / −40 s ✓** | 55% / 25% / −45 s | ✗ **worse**: too short |
| Economics: crude oil → petrol (weekly changes) | 1,746 × 1 week | 86% / 7% / 0 | 94% / 19% / 0 | ≈ tie: neither sees the 1–2 week lag |
| Weather: temperature ↔ humidity (zero-lag control) | 70,091 × 1 h | 69% / 19% / 0 ✓ | 100% / 42% / 0 ✓ | ≈ tie: adaptive raises a calibration warning |
| Neuroscience: sleep EEG delta ↔ sigma (5-min cap) | 29,941 × 1 s | 58% / 2.5% / 0 | 90% / 10% / 0 | ✓ same answer at any sampling rate (Report 2) |
| Physiology: heart rate ↔ breathing | 2,019 × 0.25 s | 65% / 6% / 0 | 79% / 21% / 0 | ≈ tie |

How to read the table:
- **Delay:** negative means the leader comes first; ✓ means it matches the literature value.
- **"better":** the adaptive score is much further above its chance level, and its delay is no worse.

**Overall: 4 better, 4 tie or mixed, 1 worse** (counting the two hydrology rows as one field).

## What the adaptive rule chose, and why (its own explanations)

- **River flow:** "1 independent sample per 2.7 days → window = 30 × 2.7 days = 80 days; tolerance ±3 days keeps fake stability ≤ 5%."
- **Electricity demand:** "1 independent sample per 1.6 days → window 48 days, clipped to 27 days (10% of the record)."
- **El Niño:** "1 independent sample per 7 months → window 18 years, clipped to 13 years; tolerance ±10 months."
- **Heart rate:** "1 independent sample per 2.5 s → window 50 s; breathing rhythm detected (7.5 s) → max delay capped at 3.75 s."
- **Jena temperature:** "1 independent sample per 39 days (because of the yearly cycle) → window 292 days; even tolerance 0 gives > 5% fake stability" → calibration warning.

## What I learned

1. **Big wins where the fixed window is far too short for the system's memory.** River flow keeps a "memory" of about 2.7 days. A 60-hour window holds only about 23 independent samples, too few for a reliable delay, so the fixed version couldn't detect a coupling every hydrologist knows exists. The same happened with electricity demand.
2. **The price of long windows is coarse timing.** With tolerance ±10 months, the El Niño delay came out at 1 month instead of the known ~3. The chance level also rises (12–42%), because there are fewer, longer windows. The adaptive version still detected coupling far above chance everywhere except the furnace, but **for delay precision, set a time-resolution cap** (`max_window`), as I did for EEG.
3. **Seasonal cycles and trends inflate the window.** Temperature's yearly cycle made "one independent sample per 39 days". Removing known cycles first (deseasonalising) is good practice before adaptive TDS. It's a candidate for an automatic step in a v3.
4. **Tiny datasets:** with 296 samples the adaptive window is clipped to 29 samples, so keep the original settings.
5. **Lags shorter than the sampling interval** (petrol at weekly resolution) can't be resolved by any setting.

## Caveats

- **Expected delays:** these come from the literature or the data agent's own cross-correlation check; they are approximate.
- **Energy:** the summers were joined end to end, which creates a few artificial seams.
- **Beijing:** "Iws" is cumulative wind speed, a rough proxy.
- **Weather control:** Jena used 100 null shifts instead of 200 because the series is long.
