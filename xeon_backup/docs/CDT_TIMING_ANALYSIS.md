# CDT Timing Analysis — March 15, 2026

## Executive Summary

The 25-month disease onset delay at origin nodes in W249-W268 is **not caused by CDT=1000**. Origin nodes are immediately marked `disease_reached=True` and exempt from CDT checks. The delay was caused by **missing `seed_vibrio`** — without environmental Vibrio seeding, the 5 initial E individuals die within ~25 days, pathogen decays to near-zero, and only the slowly-building VBNC reservoir eventually restarts infections ~25 months later.

**Fix already in code**: `seed_vibrio` default changed from `None` → `2000.0` (committed to config.py).

## Root Cause

### What CDT actually controls
CDT (cumulative_dose_threshold) only governs wavefront activation of **non-origin nodes**. When a node hasn't been reached yet, it accumulates dispersal input daily. When cumulative dose > CDT, the node activates. Origin nodes skip this entirely — they're set to `disease_reached=True` at disease_year.

### Why origin nodes had 25-month delay
Without `seed_vibrio`:
1. Disease seeded at year 1 (day 365): 5 individuals set to E state
2. E→I1→I2→death takes ~25 days at 16°C (mean durations: E=16d, I1=6d, I2=4d)
3. 5 I1/I2 individuals shed ~4-40 bact/mL/day each → peak P_k ~200-300 bact/mL
4. With K_half=1,200,000 and P_k=300: λ ≈ 0.00015/day → ~0.08 new infections/day
5. Too slow for chain transmission. All 5 die, P_k decays (half-life 1.3d at 16°C)
6. Environmental VBNC reservoir slowly rebuilds over ~25 months until P_k sustains transmission

### With seed_vibrio=2000
1. Origin nodes start with P_k=2000 bact/mL (typical P_env_max)
2. Day 0: λ ≈ 0.001/day, ~5 new infections/day across K=5000 population
3. First secondary infection: day 1
4. 10% infected: ~day 135 (4.5 months)
5. P_k stabilizes around 1200-1900 bact/mL from environmental reservoir + shedding

## Wavefront Propagation Analysis

### Current parameters
- `wavefront_D_P = 300 km` (long-range dispersal kernel for CDT)
- `wavefront_D_P_max_range = 3000 km`
- `CDT = 1000` (no decay)
- `D_P = 15 km` (regular pathogen dispersal)

### Observed wavefront timing (W249/W257 data)
Disease appears at origin nodes (month 25 due to seed_vibrio bug), then:

| Region | W249 arrival | W257 arrival | Time after origin |
|--------|-------------|-------------|-------------------|
| CA-S (non-origin) | month 25-26 | month 25-26 | 0-1 months |
| CA-C | month 26 | month 26-27 | 1-2 months |
| CA-N | month 26-29 | month 27-29 | 1-4 months |
| OR | month 27-29 | month 27-29 | 2-4 months |
| WA-O | month 28-29 | month 28-29 | 3-4 months |
| SS-S/SS-N | month 28-32 | month 29-32 | 3-7 months |
| BC-C | month 28-30 | month 29-31 | 3-6 months |
| BC-N | month 29-31 | month 29-32 | 4-7 months |
| JDF | month 29-31 | month 29-31 | 4-6 months |
| AK-FS | month 30-31 | month 30-33 | 5-8 months |
| AK-FN | month 31 | month 32-34 | 6-9 months |
| AK-PWS | month 32 | month 32-34 | 7-9 months |
| AK-WG | month 32-33 | month 33-44 | 7-19 months |
| AK-AL | month 33-94 | month 44-59 | 8-34 months |

### Comparison to observed SSWD spread
- Channel Islands → PNW (OR/WA): 3-7 months (observed: 3-7 months) ✓
- Channel Islands → BC: 7-11 months (observed: 7-12 months) ✓
- Channel Islands → SE Alaska: 12-18 months (model: 5-8 months) ⚡ faster
- Channel Islands → PWS: 18-24 months (model: 7-9 months) ⚡ ~2× faster

### CDT sensitivity (analytical)
With daily dispersal dose of 20 bact/mL/day (typical for nearby node, no decay):

| CDT | Days to activation | Months |
|-----|-------------------|--------|
| 100 | 5 | 0.2 |
| 250 | 12.5 | 0.4 |
| 500 | 25 | 0.8 |
| 1000 | 50 | 1.7 |
| 1500 | 75 | 2.5 |
| 2000 | 100 | 3.3 |

## Force of Infection Analysis

At different Vibrio concentrations (K_half=1,200,000, r=0.15, K=5000):

| P_k (bact/mL) | λ (d⁻¹) | Daily prob | Expected new E/day |
|----------------|---------|-----------|-------------------|
| 100 | 0.000051 | 0.005% | 0.3 |
| 500 | 0.000257 | 0.026% | 1.3 |
| 1,000 | 0.000514 | 0.051% | 2.6 |
| 2,000 | 0.001028 | 0.103% | 5.1 |
| 10,000 | 0.005105 | 0.509% | 25.4 |
| 100,000 | 0.047518 | 4.64% | 231.8 |
| 1,000,000 | 0.280787 | 24.5% | 1,222.8 |

Key insight: With K_half=1.2M, P_k needs to be >10,000 for rapid epidemic spread. P_env reservoir (~1,700 at 16°C) drives slow, steady infections consistent with observed multi-month outbreak dynamics.

## Recommendations

1. **Keep CDT=1000** — wavefront timing is reasonable (slightly fast for Alaska, but secondary to other parameters)
2. **Next sweep MUST use seed_vibrio=2000** — the default is already set, just ensure configs don't override it to None
3. **Optional future tuning**: If AK arrival is too fast, increase CDT to 1500-2000 or decrease wavefront_D_P to 200km
4. **K_half > 1.2M remains the dominant lever** for AK recovery — higher K_half means lower force of infection, more recovery

## Disease Stage Durations (reference)

| Stage | Rate param | T=10°C | T=13°C | T=16°C | T=20°C |
|-------|-----------|--------|--------|--------|--------|
| E→I1 | mu_EI1 | 20.8d | 18.0d | 15.6d | 12.9d |
| I1→I2 | mu_I1I2 | 8.4d | 7.0d | 5.8d | 4.6d |
| I2→D | mu_I2D | 4.5d | 4.2d | 3.9d | 3.6d |
| **Total** | | **33.8d** | **29.2d** | **25.3d** | **21.0d** |

## Files
- Analysis scripts: `scripts/cdt_timing_analysis.py`, `scripts/cdt_wavefront_analysis.py`
- Data: W249-W268 NPZ files in `results/calibration/W249-W268/`
