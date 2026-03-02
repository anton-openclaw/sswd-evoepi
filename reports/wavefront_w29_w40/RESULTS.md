# Wavefront Calibration W29-W40 Results

**Date:** 2026-03-02
**Design:** activation_threshold [50, 100] × K_half [200K, 400K, 600K] × T_vbnc [12, 14]
**Fixed:** k_vbnc=2.0, s0=0.002, P_env_max=2000, D_P=50, D_P_max_range=500km, 5 Channel Islands origins

## Critical Finding: Wavefront Cannot Reach Alaska

Across all 12 runs (and all 12 from W17-W28), disease **never reaches Alaska**.
The wavefront propagates from Channel Islands to BC-N (~24-34 months) but stalls there.
AK-PWS remains at 92-93% in every scenario — completely untouched by disease.

### Why the Wavefront Stalls
The BC coast north of BC-N has a long stretch (~1000-1500km) of remote coastline
with widely spaced sites. Even at D_P_max_range=500km, pathogen concentrations
dilute below any activation threshold before reaching Alaska. The dispersal
kernel is exponential decay with distance — beyond ~200km, concentrations
approach zero regardless of activation_threshold.

## What DID Work: Mid-Coast Recovery Gradient

### T_vbnc=14 is a Powerful Knob
Shifting the VBNC midpoint from 12°C to 14°C dramatically increases recovery:
- BC-N: 3.4% → 11.5% (K½=200K), 21.1% → 40.6% (K½=400K), 42.3% → 62.2% (K½=600K)
- The VBNC sigmoid now gives meaningful winter disease suppression at mid-latitudes

### K_half Controls Recovery Magnitude
- K½=200K: BC-N 3-14%, SS-S 5-28%
- K½=400K: BC-N 21-48%, SS-S 22-66%
- K½=600K: BC-N 42-65%, SS-S 39-73%

### Best Target Matches (Ignoring Alaska)
- **W31** (act=50, K½=400K, T_vbnc=12): BC-N=21.1% ✅, SS-S=22.1% (high), JDF=11.8% (high)
- **W37** (act=100, K½=400K, T_vbnc=12): BC-N=29.3% (close), SS-S=22.7% (high)
- **W30** (act=50, K½=200K, T_vbnc=14): BC-N=11.5% (low), SS-S=27.4% (high), JDF=21.3% (high)

## Recommendation: Abandon Pure Wavefront for Alaska

The waterborne dispersal mechanism cannot bridge the BC→AK gap with any
reasonable parameters. Options:
1. **Region-timed seeding**: Seed disease at observed arrival times per region
   (Gravem 2021: CA-S Jun 2013, OR Jan 2014, WA/BC Aug 2015, AK Jan 2017)
2. **Hybrid**: Keep wavefront for CA→BC, add timed seeding for AK
3. **Extended dispersal kernel**: Heavy-tailed (Lévy) instead of exponential

Option 1 is most pragmatic and biologically defensible (multiple introduction
pathways likely contributed to the actual northward spread).

## Summary Table

| Round | act | K_half | T_vbnc | AK-PWS | BC-N | SS-S | JDF | OR | CA-N |
|-------|-----|--------|--------|--------|------|------|-----|-----|------|
| W29 | 50 | 200K | 12 | 92.5 | 93.2 | 92.5 | 3.4 | 5.5 | 6.2 | 0.5 | 0.0 |
| W30 | 50 | 200K | 14 | 92.6 | 93.2 | 92.9 | 11.5 | 27.4 | 33.0 | 6.4 | 4.3 |
| W31 | 50 | 400K | 12 | 92.8 | 93.4 | 92.7 | 21.1 | 22.1 | 24.7 | 3.5 | 1.4 |
| W32 | 50 | 400K | 14 | 92.2 | 92.9 | 92.5 | 40.6 | 62.8 | 69.5 | 25.0 | 19.0 |
| W33 | 50 | 600K | 12 | 92.6 | 93.2 | 92.7 | 42.3 | 38.7 | 41.9 | 9.6 | 5.1 |
| W34 | 50 | 600K | 14 | 92.3 | 93.2 | 92.9 | 62.2 | 68.9 | 79.8 | 45.9 | 36.5 |
| W35 | 100 | 200K | 12 | 92.6 | 93.1 | 93.0 | 3.6 | 5.4 | 5.8 | 0.4 | 0.0 |
| W36 | 100 | 200K | 14 | 92.6 | 93.4 | 92.8 | 13.5 | 28.1 | 35.0 | 6.8 | 4.5 |
| W37 | 100 | 400K | 12 | 92.5 | 93.2 | 92.7 | 29.3 | 22.7 | 24.6 | 3.4 | 1.3 |
| W38 | 100 | 400K | 14 | 92.5 | 93.3 | 93.2 | 47.6 | 66.1 | 76.3 | 26.1 | 19.0 |
| W39 | 100 | 600K | 12 | 92.8 | 93.2 | 92.5 | 50.8 | 39.0 | 42.4 | 9.7 | 5.4 |
| W40 | 100 | 600K | 14 | 92.1 | 93.2 | 93.0 | 64.7 | 72.9 | 87.0 | 46.2 | 35.9 |

**Targets:** AK-PWS=50%, AK-FN=50%, AK-FS=20%, BC-N=20%, SS-S=5%, JDF=2%, OR=0.25%, CA-N=0.1%