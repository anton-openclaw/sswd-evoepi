# Parallelization Report — Node-Level Threading

**Date:** 2026-02-26  
**Author:** Automated (subagent)  
**Scope:** `sswd_evoepi/model.py` — daily loop in `run_spatial_simulation()`

## Summary

Added `parallel_workers` config parameter and `ThreadPoolExecutor`-based parallel
processing for all independent per-node loops in the daily simulation step.

## Changes Made

### 1. `sswd_evoepi/config.py`
- Added `parallel_workers: int = 1` to `SimulationSection`
- Default is 1 (serial) for full backward compatibility

### 2. `sswd_evoepi/model.py`
Parallelized 6 independent per-node loops in the daily step:

| Loop | Description | RNG-dependent | Parallelized |
|------|-------------|:---:|:---:|
| 2a | CRW Movement | ✅ | ✅ |
| 2b | Spatial grid build | ❌ | ✅ |
| 3  | Disease update | ✅ | ✅ |
| 4  | Continuous settlement | ✅ | ✅ |
| 5  | Daily spawning | ✅ | ✅ |
| 6  | Mortality + growth | ✅ | ✅ |

**Kept serial:**
- Environment update (step 1) — trivial cost
- Pathogen dispersal (D.T @ P) — couples all nodes
- Vibrio tracking / disease prevalence tracking — reads from shared state
- Annual demographic step — follows daily loop

### 3. Thread-Safe RNG Design
When `parallel_workers > 1`:
- Creates N independent `np.random.Generator` objects, one per node
- Each seeded from the master `rng` via `rng.integers(2**63)`
- Per-node RNGs are deterministic for a given master seed
- Serial path (`parallel_workers=1`) uses the original single `rng` unchanged

### 4. `tests/test_parallel.py`
6 new tests verifying:
- Config default is serial (workers=1)
- Serial path produces valid results
- Parallel path (workers=4) completes without errors
- Parallel and serial produce statistically comparable results
- workers > nodes (workers=8, nodes=3) works correctly

## Benchmark Results

### Small test: 5 nodes × K=100, 3 years

| Workers | Time (s) | Speedup | Final Pop | Disease Deaths |
|---------|----------|---------|-----------|----------------|
| 1       | 4.21     | 1.00×   | 145       | 319            |
| 2       | 4.98     | 0.85×   | 152       | 328            |
| 4       | 5.33     | 0.79×   | 152       | 328            |
| 8       | 5.65     | 0.75×   | 152       | 328            |

### Medium test: 10 nodes × K=500, 2 years

| Workers | Time (s) | Final Pop | Disease Deaths |
|---------|----------|-----------|----------------|
| 1       | 11.19    | 4973      | 0              |
| 2       | 11.05    | 4973      | 0              |
| 4       | 11.54    | 4973      | 0              |

## GIL Analysis

**Threading provides no speedup and slight overhead due to CPython's GIL.**

Why:
1. The per-node functions (`daily_movement`, `daily_disease_update`, etc.) are
   Python-heavy with many small NumPy array operations
2. While NumPy releases the GIL for large BLAS ops (matrix multiply, etc.), the
   per-node workload consists of many small vectorized ops (100–500 element arrays)
   where Python bookkeeping dominates
3. ThreadPoolExecutor adds ~1ms overhead per task submission, which accumulates
   over 896 nodes × 365 days × 6 phases = ~2M task submissions per year
4. On this hardware (2-core sandbox VM), even perfect threading would only give 2× max

## Recommendation: Switch to ProcessPoolExecutor

**The threading infrastructure is a dead end for this workload.** To get real
speedups, the next step should be:

### Option A: `ProcessPoolExecutor` with shared memory
- Use `multiprocessing.shared_memory` for agent arrays
- Each worker process operates on its slice of nodes
- No GIL contention → linear speedup up to core count
- Complexity: medium (need to manage shared memory lifecycle)

### Option B: Batch node processing (no parallelism needed)
- Instead of processing nodes one at a time, batch-vectorize across all nodes
- Example: `all_agents = np.concatenate([node.agents for node in nodes])`
- Run a single vectorized `daily_movement()` on the combined array
- This avoids Python loop overhead entirely and is likely 10-50× faster
- Complexity: high (need to track node boundaries within combined arrays)

### Option C: Numba/Cython JIT compilation
- JIT-compile the inner loops of `daily_movement`, `daily_disease_update`
- Releases GIL automatically → threading would then work
- Complexity: medium (need to annotate hot functions)

**Recommended path:** Option B (batch vectorization) > Option C (JIT) > Option A (multiprocessing)

Batch vectorization eliminates the parallelization problem entirely — instead of
896 serial function calls, you'd have 1 call operating on all agents at once.
NumPy's vectorized operations on arrays of 500K+ elements are extremely efficient.

## Test Results

- **776 tests collected, 775 passed, 1 pre-existing failure**
  - The single failure (`test_reintroduction_experiment::test_load_node_definitions`)
    is a pre-existing data mismatch (expects 907 nodes, data has 896) — unrelated
    to parallelization changes
- All 6 new parallel tests pass
- Serial path behavior is completely unchanged (parallel_workers=1 default)

## Config Usage

```yaml
simulation:
  parallel_workers: 4   # or any positive integer
```

Or in Python:
```python
config.simulation.parallel_workers = 4
```

Setting `parallel_workers=1` (the default) uses the original serial code path
with zero overhead.
