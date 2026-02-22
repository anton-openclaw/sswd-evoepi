# Phase 7: R→S Validation — COMPLETE ✅

**Date**: February 21, 2026, 6:03 PM PT  
**Duration**: ~3.2 min total (97.6s sinusoidal + 96.1s satellite)

## What Was Done

1. **Verified R→S fix** (commit 119497b): Recovery transitions set `DiseaseState.S` instead of `.R`
2. **Verified satellite SST** (commit d875490): Climatology data for all 5 validation nodes
3. **Test suite**: 659/659 tests pass (193.7s)
4. **Ran two validation simulations**:
   - R→S with sinusoidal SST (direct comparison to baseline)
   - R→S with satellite SST (new capability)
5. **Generated 6 comparison figures**
6. **Wrote detailed comparison analysis**

## Key Results

### R→S vs Permanent Immunity (sinusoidal SST)

| | Permanent Immunity | R→S |
|---|:---:|:---:|
| Overall crash | 98.5% | **99.7%** |
| Final pop | 365 | **122** |
| Extinctions | 0 nodes | **2 nodes** (SJI, Monterey) |
| Recovery events | 365 | 276 |
| Recovery Δc | +0.029 to +0.154 | **−0.021 to +0.030** |

### Critical Finding
Recovery trait no longer evolves upward with R→S. Selection shifts from recovery → resistance. Two of five nodes go locally extinct. The model now correctly represents echinoderm biology (no adaptive immunity) but predicts much worse conservation outcomes.

### Sinusoidal vs Satellite SST
Marginal differences in overall crash rate (99.7% vs 99.9%). Satellite SST shifts which specific nodes persist but doesn't change the qualitative outcome.

## Output Files

- `results/validation_rs_fix/results.json` — Full data
- `results/validation_rs_fix/comparison.md` — Detailed analysis
- `results/validation_rs_fix/run.log` — Raw simulation output
- 6 PNG figures in `results/validation_rs_fix/`

## Next Steps

1. Fix paper sections describing immunity/recovery behavior
2. Reconsider calibration targets given R→S dynamics
3. K=100K validation with R→S (does stochastic rescue persist?)
4. Satellite SST should become default for all future runs
