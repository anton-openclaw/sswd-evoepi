# Performance Summary — Continuous Mortality + Optimizations

**Date:** 2026-02-19

## Runtime Progression

| Stage | Time (s) | vs Baseline | Notes |
|-------|----------|-------------|-------|
| Baseline (pre-daily-mort) | 51s | — | Continuous settlement + spawning overhaul |
| After daily mortality | 107s | +109% | Correct behavior but 2× slower |
| **After optimization** | **76.7s** | **+50%** | Batched movement, grid cascade, shared alive mask |

## What Was Optimized

### 1. Movement: Batched substeps (saved ~24s)
- Was: 24 separate `update_movement()` calls per node per day, each recomputing alive mask, speeds, and RNG
- Now: Alive mask + speed computed once, all 24 substeps of RNG generated in single call, tight inner loop
- Gravity path unchanged (falls back to per-substep calls)

### 2. Spawning: Grid-based cascade induction (saved ~5s)
- Was: O(T×I) pairwise distance computation for cascade induction, chunked at 1M
- Now: Grid-cell binning with cascade_radius cells, only 9 neighbor cells checked per target
- Small problems (T×I < 500) use brute-force (lower overhead)

### 3. Shared alive mask + disease precompute (saved ~2s)
- daily_natural_mortality and daily_growth_and_aging accept optional `alive_idx` parameter
- Computed once per node per day, shared across both calls
- Recomputed between mortality→growth only if deaths occurred

## Profile Breakdown (post-optimization estimate)

| Component | Time (s) | % |
|-----------|----------|---|
| Movement | ~32s | 42% |
| Spawning | ~17s | 22% |
| Disease | ~12s | 16% |
| Daily mortality | ~2s | 3% |
| Daily growth | ~2s | 3% |
| Other | ~12s | 15% |

## SA Round 3 Timing Projections

At 77s/run, 8 cores:

| Method | Runs | Wall time |
|--------|------|-----------|
| Morris (20 traj, 39 params) | 840 | **2.2 hours** |
| Sobol N=128 | 10,023 | **27 hours** |
| Sobol N=256 | 20,046 | **53 hours** |

## Test Results

582/582 passing, zero regressions.

## Commits

- `ef25756` — daily mortality + growth functions (Phase 1)
- `4fe62d0` — wire into spatial sim (Phase 3)
- `863ae7e` — tests (Phase 4)
- `1eaffa0` — optimizations (batched movement + grid cascade + shared alive mask)
