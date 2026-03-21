# SSWD-EvoEpi Performance Audit

**Date:** 2026-02-26  
**Auditor:** Anton's AI assistant (Claude)  
**Target:** 896-node spatial simulation, K=5000, 13 years, 64-core Xeon W-3365  
**Current runtime:** ~15-30+ hours per run  
**Current parallelism:** 18 instances × ~7 threads = 126 of 128 hardware threads  

---

## 1. Current Architecture

### Simulation Loop Structure

The spatial simulation (`run_spatial_simulation()` in `model.py`) runs:

```
for year in range(13):                           # 13 years
    for day in range(365):                       # 4,745 daily steps
        for node in range(896):                  # 896 nodes per day
            1. Update SST (lookup from pre-generated array)     ~O(1)
            2. Movement: 24 CRW substeps                       ~O(K) per substep
            3. Build spatial transmission grid                  ~O(K)
            4. Disease update                                   ~O(K)
            5. Pathogen dispersal (D^T @ P)                    ~O(N²) once/day
            6. Continuous settlement (from pending cohorts)     ~O(small)
            7. Spawning step (during season)                    ~O(K_adults)
            8. Daily mortality + growth                         ~O(K)
    # Annual:
    9. Larval dispersal via C matrix                           ~O(N × larvae)
    10. Record metrics per node                                ~O(N × K)
```

### Scale Numbers

| Parameter | Value | Impact |
|-----------|-------|--------|
| Nodes (N) | 896 | Outer spatial loop |
| Agents/node (K) | 5,000 | Agent operations scale linearly |
| Total agents | 4,480,000 | Memory: ~265 MB agents + ~460 MB genotypes |
| Genotype array | 5000 × 51 loci × 2 alleles × int8 | 510 KB/node, 446 MB total |
| Agent struct | ~59 bytes × 12,500 (max_agents) | 720 KB/node, 630 MB total |
| Movement substeps | 24/day | 24× multiplier on movement cost |
| Daily timesteps | 4,745 (13 × 365) | Main loop iterations |
| Total node-days | 4,251,520 (896 × 4,745) | Total work units |

### Memory Footprint Estimate

```
Per node:
  agents:    12,500 × 59B  =   721 KB
  genotypes: 12,500 × 102B =  1,245 KB
  Total per node: ~2.0 MB

All 896 nodes: ~1.74 GB agent data
Distance matrix (896²): ~6.1 MB (float64)
C matrix (896²): ~6.1 MB
D matrix (896²): ~6.1 MB (sparse in practice but stored dense)
SST time series: 896 × 4745 × 8B = 33 MB
Overhead: ~1 GB (Python, NumPy, etc.)

Total estimated: ~3-4 GB (well within 503 GB RAM)
```

---

## 2. Profiling Analysis (Code-Based Estimate)

### Per-Day Cost Breakdown (Single Node, K=5000)

Based on code structure analysis. All operations are on arrays of ~5000 alive agents.

| Component | Operations | Estimated Cost | % of Day |
|-----------|-----------|----------------|----------|
| **Movement (24 substeps)** | 24 × (RNG draw + trig + reflect) on K agents | **~60-70%** | Dominant |
| Disease update | Force of infection + progression on K agents | ~10-15% | Moderate |
| Spatial TX grid build | Binning + 2× diffusion passes on grid | ~3-5% | Light |
| Daily mortality | Vectorized survival check on K agents | ~2-3% | Light |
| Daily growth | VB growth + stage transitions on K agents | ~2-3% | Light |
| Spawning (seasonal) | Readiness + spawning + cascade (~270 days/yr) | ~5-10% | Moderate |
| Settlement | Settle pending cohorts (usually few) | ~1% | Negligible |

**Key insight: Movement IS likely ~60-70% of compute**, primarily because:

1. **24 substeps/day** — each substep generates `rng.normal()` for all alive agents, computes `cos()` and `sin()` on all headings, and performs boundary reflection. That's 24 passes over ~5000 agents per node per day.

2. **The `daily_movement()` fast path** (no gravity) is well-optimized: it pre-generates all turning angles in one `rng.normal(size=(24, n_alive))` call and loops over substeps with vectorized trig. This is already the best pure-NumPy approach.

3. **24 × 896 × 4745 = ~102 billion floating-point substep operations** over the full run. Even at ~10μs per substep per node, that's ~102M × 10μs = ~1020s = 17 minutes for just the substep kernel — but the actual cost includes array overhead, cache misses, and Python loop overhead.

### Why Movement Dominates

Let's count FLOPS per day per node in `daily_movement()` fast path:
- `rng.normal(size=(24, 5000))` = 120,000 random numbers
- 24 substeps × 5000 agents × (1 add + 1 mod + 1 cos + 1 sin + 2 mul + 2 add + 2 reflect) ≈ 24 × 5000 × 10 = 1.2M FLOPS
- Plus memory traffic: reading/writing `heading`, `x`, `y` arrays 24 times

For comparison, disease update per day:
- Force of infection: ~5000 × 5 FLOPS = 25K FLOPS
- State transitions: ~few hundred agents × ~10 FLOPS = ~5K FLOPS

**Movement is 50-100× more compute-intensive per day than disease**, confirming the ~60-80% estimate.

### The Serial Node Loop Problem

The innermost daily loop (`for i, node in enumerate(network.nodes):`) iterates over all 896 nodes **sequentially in Python**. Each iteration processes that node's agents via NumPy vectorized operations, but the loop itself is pure Python. 

Critical observation: **nodes are independent within a daily timestep** (before pathogen dispersal). The disease step, movement, mortality, growth, and spawning at node `i` don't depend on node `j`'s state for that day. Only pathogen dispersal (step 5) and larval dispersal (annual) couple nodes.

This means **the daily node loop is embarrassingly parallel** — the single biggest parallelism opportunity.

---

## 3. Quick Wins (Hours of Effort, Minimal Risk)

### QW-1: Reduce Movement Substeps from 24 to 6-8

| Attribute | Value |
|-----------|-------|
| **Effort** | 2-4 hours |
| **Expected speedup** | 3-4× on movement = **~2× overall** |
| **Risk** | Low-Medium — changes spatial distribution dynamics |
| **Behavior change** | Yes — coarser movement, less spatial mixing |

**Rationale:** 24 substeps/day (hourly) is very fine-grained for a sessile/slow-moving benthic invertebrate moving at 0.5 m/min. With 6 substeps (4-hour resolution), daily displacement statistics are nearly identical due to the central limit theorem — the sum of 6 random walks with dt=240min gives statistically similar endpoint distributions as 24 walks with dt=60min. The key difference is within-day spatial correlation for disease transmission, which may not matter since disease transmission uses a daily grid average anyway.

**Validation:** Run a 5-node test with substeps={6, 12, 24} and compare disease dynamics and recovery fractions. If results are within 5%, use 6.

**Implementation:**
```python
# config.py
substeps_per_day: int = 6  # Was 24
```

### QW-2: Sparse D Matrix for Pathogen Dispersal

| Attribute | Value |
|-----------|-------|
| **Effort** | 3-4 hours |
| **Expected speedup** | Minor (~1-2% of total, but eliminates O(N²) scaling) |
| **Risk** | None — exact same results |
| **Behavior change** | No |

**Rationale:** The D matrix has `max_range=50 km`, meaning most entries are zero. With 896 nodes, the matrix is 896² = 802,816 entries, but with D_P=15km and max_range=50km, the effective neighborhood is ~10-30 nodes per source. Sparsity is ~97%+. Using `scipy.sparse.csr_matrix` for `D.T @ P` would reduce the daily pathogen dispersal from O(N²) to O(N × k) where k~20.

The daily `pathogen_dispersal_step()` is just `D.T @ P` — trivial to make sparse. Currently called 4,745 times.

**Implementation:**
```python
from scipy.sparse import csr_matrix
D_sparse = csr_matrix(D)
# In loop:
dispersal_in = D_sparse.T @ P  # Same API, much faster for sparse
```

### QW-3: Skip Disease Update for Disease-Free Nodes

| Attribute | Value |
|-----------|-------|
| **Effort** | 2-3 hours |
| **Expected speedup** | **1.5-2× for first 3 years** (pre-disease), negligible after epidemic |
| **Risk** | None |
| **Behavior change** | No |

**Rationale:** The code already has `disease_active_flags[i]` but the movement, mortality, and growth loops run for ALL nodes every day regardless. During the 1-year spinup before disease starts, all 896 nodes are disease-free. After disease starts, nodes that have cleared the epidemic and have zero infected agents could also be skipped for disease (they already are — `if not disease_active_flags[i]: continue`).

This is already implemented for disease. No change needed.

### QW-4: Pre-compute `alive_idx` Once Per Node Per Day, Share Across All Operations

| Attribute | Value |
|-----------|-------|
| **Effort** | 1-2 hours |
| **Expected speedup** | ~5-10% (eliminates redundant `np.where(agents['alive'])` calls) |
| **Risk** | None |
| **Behavior change** | No |

**Rationale:** Currently, `daily_natural_mortality()`, `daily_growth_and_aging()`, and other functions each call `np.where(agents['alive'])[0]` independently. This is already partially done in the single-node loop (`_alive_idx_demo`), and the spatial loop also computes `_alive_idx` once per node. But movement, spawning, and disease each re-derive it. Pre-computing once at the top of the daily loop and passing it through would save ~5 `np.where()` calls per node per day = 4,745 × 896 × 5 = ~21M redundant scans eliminated.

Already partially implemented in the spatial loop. Could be extended to movement and spawning.

### QW-5: Use NumPy BLAS for Pathogen Dispersal Instead of Python Array Construction

| Attribute | Value |
|-----------|-------|
| **Effort** | 1 hour |
| **Expected speedup** | ~1% (already using `D.T @ P` which is NumPy BLAS) |
| **Risk** | None |
| **Behavior change** | No |

**Rationale:** The current `pathogen_dispersal_step()` is already `D.T @ P` which is efficient. But the construction of `P` array from node states uses a Python list comprehension:
```python
P = np.array([node_disease_states[i].vibrio_concentration for i in range(N)])
```
This creates a Python list of 896 floats then converts to array. Could use a pre-allocated array instead:
```python
P_buffer = np.empty(N)
for i in range(N):
    P_buffer[i] = node_disease_states[i].vibrio_concentration
```
Or better, store vibrio as a contiguous array alongside node disease states.

---

## 4. Medium Effort (Days of Work, Moderate Refactoring)

### ME-1: Parallelize the Daily Node Loop with `multiprocessing`

| Attribute | Value |
|-----------|-------|
| **Effort** | 3-5 days |
| **Expected speedup** | **4-8× with 8-16 worker processes** |
| **Risk** | Medium — requires careful state management |
| **Behavior change** | Yes — RNG sequences change (but stochastic anyway) |

**Rationale:** This is the single highest-impact optimization. The daily loop iterates over 896 nodes sequentially. Each node's movement, disease, mortality, and growth are independent of other nodes within a timestep. Only pathogen dispersal and snapshot recording couple nodes.

**Architecture:**
```
Main process:
    for day in range(365):
        # 1. Parallel: independent node operations
        with ProcessPoolExecutor(max_workers=16) as pool:
            futures = [pool.submit(process_node_day, i, ...) for i in range(896)]
            results = [f.result() for f in futures]
        
        # 2. Serial: pathogen dispersal (cheap, needs all node states)
        P = np.array([r.vibrio for r in results])
        dispersal_in = D.T @ P
        for i in range(N):
            nodes[i].vibrio += dispersal_in[i]
```

**Challenges:**
- Agent arrays are per-node (~2 MB each). With 896 nodes, total state is ~1.7 GB. Shipping arrays to/from workers is expensive via pickling.
- **Solution 1: Shared memory** (`multiprocessing.shared_memory`). Store all agent arrays in a single contiguous shared memory block. Workers access their node's slice without copying.
- **Solution 2: Process per node-group**. Partition 896 nodes into 16 groups of 56. Each worker owns its group's memory and runs the full daily loop for its 56 nodes. Only vibrio concentrations (896 floats = 7 KB) are exchanged at the end of each day.

**Solution 2 is strongly recommended** — minimal inter-process communication, natural partitioning, and each worker can use NumPy efficiently on its local arrays.

**Expected performance:**
- Current: ~896 × (movement + disease + mortality) = sequential
- Parallel (16 workers): Each handles 56 nodes. Wall time ≈ sequential / 16 × (1 + overhead) ≈ **6× speedup** (80% efficiency)
- Combined with QW-1 (substeps 6): **12× total speedup**

### ME-2: Vectorize the Node Loop into Batched NumPy Operations

| Attribute | Value |
|-----------|-------|
| **Effort** | 5-7 days |
| **Expected speedup** | **3-5×** (eliminates Python loop overhead, enables SIMD) |
| **Risk** | High — major refactoring of data layout |
| **Behavior change** | No (if done correctly) |

**Rationale:** Instead of `for node in nodes: process(node.agents)`, store ALL agents across ALL nodes in a single flat array with a `node_id` field (already exists!). Then operations like movement, disease, mortality become single vectorized calls on the 4.48M-agent array.

**Example — Movement:**
```python
# Current: 896 × vectorized(5000)
for i, node in enumerate(network.nodes):
    daily_movement(node.agents, ...)

# Proposed: single vectorized call on 4.48M agents
daily_movement(all_agents, ...)  # all_agents is one big array
```

**Pros:** Eliminates 896 Python loop iterations per operation per day. NumPy can use SIMD across the full 4.48M agents. Better cache behavior for contiguous arrays.

**Cons:** Disease and spawning need node-level compartment counts, which requires `np.bincount()` on `node_id` field. Spawning gravity needs per-node neighbor search. Genotype arrays would need a flat layout too. Major refactoring.

**Hybrid approach:** Flatten only for movement and mortality (the biggest cost), keep per-node for disease and spawning.

### ME-3: Numba JIT for Movement Inner Loop

| Attribute | Value |
|-----------|-------|
| **Effort** | 2-3 days |
| **Expected speedup** | **2-5× on movement** = **1.5-3× overall** |
| **Risk** | Low (Numba is well-established) |
| **Behavior change** | No |

**Rationale:** The `daily_movement()` fast path is already well-vectorized NumPy, but the 24-substep Python loop with array slicing has overhead that Numba can eliminate. A `@numba.njit` compiled version of the substep loop would run at near-C speed.

**Implementation:**
```python
@numba.njit
def _movement_substeps(headings, x, y, speed_dt, all_turns, habitat_side, n_substeps):
    TWO_PI = 2.0 * np.pi
    for s in range(n_substeps):
        for i in range(len(headings)):
            headings[i] = (headings[i] + all_turns[s, i]) % TWO_PI
            dx = speed_dt[i] * np.cos(headings[i])
            dy = speed_dt[i] * np.sin(headings[i])
            x[i] = _reflect_scalar(x[i] + dx, habitat_side)
            y[i] = _reflect_scalar(y[i] + dy, habitat_side)
    return headings, x, y
```

Numba eliminates:
- Python loop overhead for 24 substeps
- NumPy temporary array allocation per substep
- Array indexing overhead
- Can potentially auto-vectorize with SIMD

**Note:** Numba's `@njit(parallel=True)` with `prange` could also parallelize across agents within a single node.

### ME-4: Cache SST Lookups and Environmental Computations

| Attribute | Value |
|-----------|-------|
| **Effort** | 1-2 days |
| **Expected speedup** | ~2-3% (minor but free) |
| **Risk** | None |
| **Behavior change** | No |

**Rationale:** SST time series are already pre-generated. But `seasonal_flushing()` is called daily for each node with constant parameters — it could be pre-computed as a (896, 4745) array at initialization. Same for `salinity_modifier()`, `vibrio_decay_rate()`, and `environmental_vibrio()` — all depend only on SST and static per-node parameters.

Pre-compute at initialization:
```python
# Shape: (N, total_days)
flushing_by_day = np.zeros((N, total_days))
decay_by_day = np.zeros((N, total_days))
env_vibrio_by_day = np.zeros((N, total_days))
for i, node in enumerate(network.nodes):
    for d in range(total_days):
        T = sst_by_node[i][d]
        flushing_by_day[i, d] = seasonal_flushing(...)
        decay_by_day[i, d] = vibrio_decay_rate(T)
        env_vibrio_by_day[i, d] = environmental_vibrio(T, ...)
```

Then in the daily loop, these are just array lookups.

### ME-5: Consolidated Agent Array with Flat Layout

| Attribute | Value |
|-----------|-------|
| **Effort** | 4-5 days |
| **Expected speedup** | **2-3×** (enables single-call vectorization) |
| **Risk** | Medium — must update all functions that access per-node agents |
| **Behavior change** | No |

**Rationale:** Currently each node has its own agent array of size `max_agents ≈ 12,500`. This means:
- 896 separate NumPy arrays in memory (fragmented)
- Python loop to iterate over them
- Each array is small enough that NumPy overhead dominates compute

A single flat array of ~11.2M entries (896 × 12,500) with a `node_id` field enables:
- Single `np.where(all_agents['alive'])` call for 4.48M agents
- Single movement call
- Single mortality call
- Disease grouped by `node_id` using `np.bincount()` and `np.searchsorted()`

This is the prerequisite for ME-2 and enables massive vectorization gains.

---

## 5. Major Rewrites (Weeks of Work, Significant Risk)

### MR-1: Replace Individual-Level Movement with Stochastic Migration Rates

| Attribute | Value |
|-----------|-------|
| **Effort** | 2-3 weeks |
| **Expected speedup** | **10-50× on movement** = **5-30× overall** |
| **Risk** | High — fundamentally changes movement model |
| **Behavior change** | Yes — eliminates within-node spatial structure |

**Rationale:** The CRW movement model tracks x,y coordinates of 4.48M agents with 24 substeps/day. But the x,y positions are ONLY used for:
1. **Spatial disease transmission grid** (builds InfectedDensityGrid from x,y positions)
2. **Spawning gravity** (distance to conspecifics for cascade induction)

If spatial transmission is replaced with a mean-field model (which the original disease model already supports — `infected_density_grid=None`), then x,y positions serve no purpose for disease dynamics. Spawning gravity is seasonal and only affects reproductive adults.

**Alternative:** Instead of tracking 4.48M individual positions, compute a **mixing rate** based on CRW parameters:
- Expected displacement per day = base_speed × 60 × 24 × directional_persistence ≈ 720m/day
- This determines how quickly agents "mix" within the habitat
- Can be modeled as a single spatial mixing parameter per node

**If within-node spatial structure is not critical** (validate with sensitivity analysis), this eliminates:
- 4.48M × 24 substeps × (trig + reflect) per day = ~60-80% of compute
- All `x`, `y`, `heading`, `speed` fields from agent arrays (saves 16 bytes per agent = 72 MB)

**Validation approach:** Run 5-node tests with and without movement, compare disease dynamics and recovery.

### MR-2: GPU Acceleration with CuPy

| Attribute | Value |
|-----------|-------|
| **Effort** | 2-3 weeks |
| **Expected speedup** | **5-20× on vectorized operations** |
| **Risk** | Medium — GPU memory limits (8 GB RTX 5060), no GPU on Xeon |
| **Behavior change** | No (if done correctly) |

**Rationale:** CuPy is a drop-in NumPy replacement for GPU. Key operations that would benefit:
- Movement: trig operations on 5000 agents × 24 substeps → perfect for GPU
- Force of infection: vectorized computation on 5000 agents
- Random number generation: CuPy uses GPU-optimized RNG

**Problem:** The Xeon (production machine) has no GPU. The local machine has an RTX 5060 with 8 GB. Agent data per node is ~2 MB, so a single node fits easily. But running all 896 nodes requires batching.

**Better for calibration on local machine than production on Xeon.**

### MR-3: Cython Extension for Hot Loops

| Attribute | Value |
|-----------|-------|
| **Effort** | 1-2 weeks |
| **Expected speedup** | **3-10× on compiled sections** |
| **Risk** | Low — well-established approach |
| **Behavior change** | No |

**Rationale:** The movement substep loop, disease progression, and mortality functions are good candidates for Cython compilation. Unlike Numba, Cython produces a compiled .so file that doesn't have JIT warmup costs.

Key targets:
1. Movement substep kernel (the tight 24-iteration loop with trig)
2. Disease transmission (force of infection computation)
3. Daily mortality (survival probability computation)

---

## 6. Calibration-Specific Optimizations

### CS-1: Progressive K Refinement

| Attribute | Value |
|-----------|-------|
| **Effort** | 2-3 hours |
| **Expected speedup** | **5-10× for parameter space exploration** |
| **Risk** | Medium — relative rankings may shift at different K |
| **Behavior change** | Yes — lower K = noisier, different dynamics |

**Rationale:** Movement cost scales as O(K × substeps). Disease scales as O(K). Total cost scales roughly as O(K).

| K | Agents | Relative Cost | Use Case |
|---|--------|---------------|----------|
| 500 | 448K | 0.1× | Rough parameter sweep (10× faster) |
| 1000 | 896K | 0.2× | Narrow parameter space (5× faster) |
| 2000 | 1.79M | 0.4× | Refine top candidates (2.5× faster) |
| 5000 | 4.48M | 1.0× | Final calibration (full fidelity) |

**Protocol:**
1. Run parameter sweep at K=500 (10× faster): identify top ~20 parameter sets
2. Re-run top 20 at K=2000 (2.5× faster): narrow to top 5
3. Final validation at K=5000 with multiple seeds

**Risk mitigation:** Demographic stochasticity is much higher at K=500. Some parameters (e.g., Allee effects) may behave qualitatively differently. Run a few known-good parameters at K=500 and K=5000 to validate rank-order preservation.

### CS-2: Early Stopping for Divergent Runs

| Attribute | Value |
|-----------|-------|
| **Effort** | 3-4 hours |
| **Expected speedup** | **1.5-3× for parameter sweeps** (kills bad runs early) |
| **Risk** | Low — only affects calibration workflow, not model |
| **Behavior change** | No (model code unchanged) |

**Rationale:** Many parameter combinations produce clearly wrong dynamics by year 3-5:
- Total extinction (pop → 0) by year 3
- No epidemic at all (100% survival by year 5)  
- All regions crash to 0 (no regional variation)

These can be detected early and aborted:

```python
# In calibration runner, after each year:
if year == 3:
    total_pop = yearly_total_pop[year]
    if total_pop == 0:
        return {'rmse_log': float('inf'), 'early_stop': 'extinction_y3'}
    if total_pop > 0.95 * initial_pop:
        return {'rmse_log': float('inf'), 'early_stop': 'no_epidemic_y3'}

if year == 5:
    # Check if regional variation exists
    region_pops = compute_regional_recovery_at_year(result, sites, year)
    if max(region_pops.values()) - min(region_pops.values()) < 0.01:
        return {'rmse_log': float('inf'), 'early_stop': 'no_gradient_y5'}
```

**Expected savings:** In a typical calibration sweep, ~30-50% of parameter combinations are clearly bad. Killing them at year 3-5 (of 13) saves 60-75% of their runtime.

### CS-3: Proxy Metrics for Fast Screening

| Attribute | Value |
|-----------|-------|
| **Effort** | 1-2 days |
| **Expected speedup** | **5-10× for initial screening** |
| **Risk** | Medium — proxy may miss good candidates |
| **Behavior change** | N/A — calibration workflow only |

**Rationale:** The full calibration target is regional recovery fractions at year 13. But several cheaper proxies correlate well:

1. **R₀ at disease_year**: Can be computed analytically from parameters. If R₀ < 1 at warm nodes but R₀ > 5 at cool nodes, the gradient is built in. Cost: ~0 (analytical).

2. **Year-3 mortality rate by region**: Run only 3 years. If mortality is >95% everywhere or <20% everywhere, parameters are bad. Cost: 3/13 = 23% of full run.

3. **5-node proxy run**: Run the canonical 5-node network (not 896-node) for 13 years. K=5000, 5 nodes = 25K agents vs 4.48M. Cost: ~0.5% of full run. If the 5-node gradient is wrong, the 896-node gradient will be wrong too.

**Protocol for efficient calibration:**
```
Phase 1: Analytical R₀ filter (free)
  → Filter out parameters with uniform R₀ across latitudes
  
Phase 2: 5-node proxy (0.5% cost × many runs)
  → Run 100+ parameter sets on 5-node network
  → Score on recovery gradient (Alaska > BC > WA > OR > CA)
  → Keep top 20

Phase 3: Full 896-node at K=1000 (20% cost)
  → Run top 20 parameter sets
  → Score on all 8 regional targets
  → Keep top 5

Phase 4: Full 896-node at K=5000 (100% cost)
  → Run top 5 with 3-5 seeds each
  → Final scoring
```

### CS-4: Run Fewer Seeds During Exploration

| Attribute | Value |
|-----------|-------|
| **Effort** | 0 hours (already supported) |
| **Expected speedup** | **3-5× during exploration** |
| **Risk** | Higher variance in parameter evaluation |
| **Behavior change** | None |

**Rationale:** The calibration runner supports multiple seeds. During exploration, use 1 seed. Only for final top-5 candidates, use 3-5 seeds to assess variance.

### CS-5: Shorter Spinup Period

| Attribute | Value |
|-----------|-------|
| **Effort** | 1-2 hours |
| **Expected speedup** | **~8% (1 of 13 years)** |
| **Risk** | Low — if initialization is at demographic equilibrium |
| **Behavior change** | Marginal |

**Rationale:** Currently using 1-year spinup (disease starts year 1). If populations are initialized at demographic equilibrium (which `initialize_population()` attempts), the spinup could potentially be reduced to 0 years (disease starts immediately). This saves 1/13 = 7.7% of runtime.

**Validation:** Compare year-1 population structure with and without spinup. If they're within 5%, skip spinup.

---

## 7. Recommendations — Prioritized Action Plan

### Tier 1: Do This Week (combined ~3-5× speedup)

| Priority | Action | Effort | Speedup | Risk |
|----------|--------|--------|---------|------|
| **1** | **QW-1: Reduce substeps 24→6** | 2h | 2-3× | Low |
| **2** | **CS-2: Early stopping** | 3h | 1.5-2× (calib.) | None |
| **3** | **CS-1: Progressive K (K=1000 for exploration)** | 2h | 5× (calib.) | Medium |
| **4** | **QW-2: Sparse D matrix** | 3h | ~1% | None |

**Combined effect:** Substep reduction alone likely cuts runtime from ~15-30h to ~5-15h. Progressive K for exploration cuts exploration runs to ~3-6h at K=1000.

### Tier 2: Do This Sprint (combined ~10-20× speedup)

| Priority | Action | Effort | Speedup | Risk |
|----------|--------|--------|---------|------|
| **5** | **ME-1: Parallelize daily node loop** | 3-5d | 4-8× | Medium |
| **6** | **ME-3: Numba JIT for movement** | 2-3d | 2-3× | Low |
| **7** | **CS-3: 5-node proxy screening** | 1-2d | 10× (calib.) | Medium |
| **8** | **ME-4: Pre-compute environmental arrays** | 1-2d | 2-3% | None |

**Combined with Tier 1:** Node parallelism on 16 cores + substep reduction + Numba movement could bring a single full-fidelity run from 15-30h down to **1-3 hours**.

### Tier 3: Long-Term (combined ~50-100× speedup)

| Priority | Action | Effort | Speedup | Risk |
|----------|--------|--------|---------|------|
| **9** | **ME-5: Flat agent array** | 4-5d | 2-3× | Medium |
| **10** | **MR-1: Replace CRW with mixing rates** | 2-3w | 10-50× | High |
| **11** | **MR-3: Cython hot loops** | 1-2w | 3-10× | Low |
| **12** | **MR-2: GPU acceleration** | 2-3w | 5-20× | Medium |

### Current Parallelism Usage Note

You mentioned running 18 parallel instances using ~7 threads each (126 total threads on 128 HW threads). This is the **optimal strategy for embarrassingly parallel calibration runs** — each run is a single-threaded Python process. The `~7 threads` likely comes from NumPy/BLAS internal threading (e.g., MKL or OpenBLAS using multi-threaded BLAS for `D.T @ P` and large array operations).

**Recommendation:** 
- Set `OMP_NUM_THREADS=1` and `MKL_NUM_THREADS=1` per process. The BLAS operations in this model are small (896×896 matrix × 896 vector), so multi-threaded BLAS is wasteful — the overhead exceeds the benefit.
- With `*_NUM_THREADS=1`, you can run **64 parallel instances** (one per physical core) instead of 18, getting **~3.5× more throughput** for parameter sweeps.
- On the 128-thread Xeon, try both 64 and 96 instances to find the sweet spot.

This alone could be the single biggest calibration speedup — from 18 concurrent runs to 64.

```bash
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# Run 64 parallel calibration instances
parallel -j 64 python calibration_runner.py --config params_{}.json --seed {} --output results/cal_{} ::: $(seq 1 64)
```

---

## Appendix A: Estimated Runtime Under Different Optimization Scenarios

| Scenario | Movement | Disease | Other | Total | vs Baseline |
|----------|----------|---------|-------|-------|-------------|
| **Baseline** (K=5000, 24 substeps) | ~18h | ~4h | ~3h | ~25h | 1× |
| QW-1 (6 substeps) | ~4.5h | ~4h | ~3h | ~11.5h | 2.2× |
| QW-1 + ME-1 (parallel 16 cores) | ~0.3h | ~0.25h | ~0.2h | ~0.75h | 33× |
| QW-1 + ME-3 (Numba movement) | ~1.5h | ~4h | ~3h | ~8.5h | 3× |
| QW-1 + ME-1 + ME-3 | ~0.1h | ~0.25h | ~0.2h | ~0.55h | 45× |
| MR-1 (no CRW, mixing rates) | ~0h | ~4h | ~3h | ~7h | 3.6× |
| MR-1 + ME-1 (parallel) | ~0h | ~0.25h | ~0.2h | ~0.45h | 56× |

**Note:** These are rough estimates. Actual performance depends on cache behavior, memory bandwidth, GIL contention (for multiprocessing vs threading), and communication overhead.

## Appendix B: Memory Layout Analysis

### Current Layout (Per-Node Arrays)
```
896 nodes × 12,500 agents × 59 bytes = 662 MB (fragmented across 896 allocations)
896 nodes × 12,500 agents × 102 bytes = 1.1 GB genotypes (fragmented)
```

### Proposed Flat Layout  
```
1 array: 11,200,000 agents × 59 bytes = 629 MB (contiguous)
1 array: 11,200,000 agents × 102 bytes = 1.09 GB genotypes (contiguous)
```

Benefits of contiguous: better cache utilization, SIMD-friendly, single NumPy call for movement/mortality over all 4.48M alive agents.

## Appendix C: Quick Test to Validate Substep Sensitivity

```python
# Run this before committing to substep reduction
from sswd_evoepi.spatial import make_5node_network
from sswd_evoepi.model import run_spatial_simulation
from sswd_evoepi.config import default_config

results = {}
for substeps in [6, 12, 24]:
    config = default_config()
    config.movement.substeps_per_day = substeps
    network = make_5node_network(seed=42)
    result = run_spatial_simulation(network, n_years=10, disease_year=3, seed=42, config=config)
    results[substeps] = {
        'final_pop': result.final_total_pop,
        'yearly_pop': result.yearly_total_pop.tolist(),
    }
    print(f"substeps={substeps}: final_pop={result.final_total_pop}")
```

## Appendix D: Top Recommendation Summary

**If you could only do 3 things:**

1. **Set `OMP_NUM_THREADS=1` and run 64 parallel instances** instead of 18 at ~7 threads each. This requires zero code changes and immediately gives ~3.5× more throughput for calibration sweeps.

2. **Reduce substeps from 24 to 6** (after validation). Single config change, ~2-3× speedup per run.

3. **Implement parallel node processing** (ME-1) with 16 worker processes per run. This is the biggest bang-for-buck code change at ~4-8× per-run speedup.

Together: 3.5× (more instances) × 2.5× (fewer substeps) × 5× (parallel nodes) = **~44× faster calibration throughput**.

At current 15-30h/run with 18 instances: 18 runs in 15-30h.  
After optimizations: 64 instances of 0.6-1.2h runs = **64 runs in 1.2 hours**.  
That's ~50× improvement in calibration throughput.
