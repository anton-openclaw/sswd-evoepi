# OPTIMIZATION_PLAN.md — Profiling Analysis & Optimization Targets

**Date:** 2026-02-19  
**Config:** 3-node (Sitka/Howe Sound/Monterey), K=5000/node, 20 years, pathogen evolution ON, seed=42  
**Wall clock:** 97.67s (instrumented run), 110.51s (cProfile run w/ overhead)

---

## 1. Component Breakdown (Monkey-Patched Timers)

| Component | Total (s) | Per-day (ms) | % of Wall | Calls | Per-call (µs) |
|---|---|---|---|---|---|
| **movement** | 55.522 | 7.606 | 56.8% | 21,900 | 2,535 |
| daily_growth | 1.793 | 0.246 | 1.8% | 21,900 | 82 |
| daily_mortality | 1.616 | 0.221 | 1.7% | 21,900 | 74 |
| settlement | 0.183 | 0.025 | 0.2% | 1,472 | 124 |
| environment | 0.050 | 0.007 | 0.1% | 43,800 | 1.1 |
| pathogen_dispersal | 0.041 | 0.006 | <0.1% | 7,300 | 5.6 |
| larval_dispersal | 0.004 | 0.001 | <0.1% | 20 | 200 |
| **TRACKED** | **59.208** | 8.111 | **60.6%** | | |
| **UNTRACKED** | **38.462** | 5.269 | **39.4%** | | |
| **WALL CLOCK** | **97.670** | | **100%** | | |

### Resolving the Untracked 39.4%

The monkey patches failed to intercept disease and spawning because
`run_spatial_simulation` binds references at import time. cProfile reveals
what's hiding in that 38.5s:

| Hidden Component | Cumtime (s) | Self-time (s) | Calls | Source |
|---|---|---|---|---|
| **spawning_step** | 22.77 | 2.90 | 15,420 | spawning.py:125 |
| **daily_disease_update** | 11.90 | 7.70 | 17,520 | disease.py:596 |
| **run_spatial_simulation (loop overhead)** | — | 6.10 | 1 | model.py:1503 |

**Sum: 22.8 + 11.9 + 6.1 = 40.8s** — accounts for the untracked time (with cProfile overhead).

### Corrected Full Breakdown

| Component | Time (s) | % of Runtime |
|---|---|---|
| **Movement** | 55.5 | **56.8%** |
| **Spawning** | 22.8 | **23.3%** |
| **Disease** | 11.9 | **12.2%** |
| Main loop overhead | 6.1 | 6.2% |
| Growth + mortality | 3.4 | 3.5% |
| Everything else | <0.3 | <0.3% |

---

## 2. cProfile Deep Dive — Function-Level Hotspots

### Movement (55.5s = 56.8%)

| Function | Self-time (s) | Calls | Per-call (µs) | Notes |
|---|---|---|---|---|
| `update_movement` | 31.39 | 525,600 | 59.7 | 24 substeps × 3 nodes × 7,300 days |
| `_reflect` | 18.73 | 1,051,200 | 17.8 | 2 calls per substep (x, y) |
| `movement.build` | 0.93 | 21,900 | 42.5 | Per-node-day grid rebuild |
| `_diffuse_2d` | 0.33 | 13,436 | 24.3 | Grid diffusion for spatial disease |
| `np.pad` | 0.14 | 13,436 | 10.4 | Grid boundary padding |

**Root cause:** 525,600 `update_movement` calls, each doing:
- `alive.astype(bool)` — array copy every call
- `np.where(alive_mask)` — scan all agents
- `rng.normal()` — RNG draw
- 2× `_reflect()` — modular arithmetic + `np.where`
- `SPEED_MODIFIER[disease_states]` — fancy index

The per-call cost (60µs) is already lean. The problem is pure **call volume**: half a million invocations. This is O(N × substeps × K × T).

### Spawning (22.8s = 23.3%)

| Function | Self-time (s) | Calls | Per-call (ms) | Notes |
|---|---|---|---|---|
| `_generate_larval_cohort` | 7.87 | 2,311 | 3.41 | Genotype assembly |
| `_check_cascade_induction` | 2.97 | 6,068 | 0.49 | O(T×I) pairwise distances |
| `_cascade_induction_step` | 2.94 | 15,420 | 0.19 | Orchestration overhead |
| `spawning_step` (self) | 2.90 | 15,420 | 0.19 | Mask computation, conditionals |
| `_readiness_induction_step` | 0.82 | 2,311 | 0.36 | Readiness checks |

**Root cause for `_generate_larval_cohort`:** At 3.4ms/call × 2,311 calls = 7.9s. This is the genotype assembly for offspring — previously optimized from 22s→0.17s via batch RNG (Feb 16), but now called in a different context with possibly larger cohort sizes or different code path.

**Root cause for `_check_cascade_induction`:** Computes O(T×I) pairwise distance matrix between targets and inducers. With the memory guard at T×I > 1M, large spawning events trigger chunked computation. 6,068 calls suggests cascade events are frequent.

### Disease (11.9s = 12.2%)

| Function | Self-time (s) | Calls | Per-call (µs) | Notes |
|---|---|---|---|---|
| `daily_disease_update` (self) | 7.70 | 17,520 | 440 | Core engine |
| `sample_stage_duration` | 0.23 | 83,689 | 2.8 | Gamma RNG per transition |
| `arrhenius` | 0.18 | 220,009 | 0.8 | exp() per rate computation |
| `sigma_1_strain` | 0.07 | 11,702 | 5.8 | Pathogen evo shedding |
| `sigma_2_strain` | 0.06 | 12,589 | 5.1 | Pathogen evo shedding |

**Root cause:** 7.7s self-time in `daily_disease_update` is spent on:
- `alive.astype(bool)` — redundant array copy (also done in movement, growth, mortality)
- `np.sum(alive_mask & (ds == X))` — 3+ compartment counts per call
- Per-susceptible dose-response loop with individual resistance values
- Pathogen evolution: per-individual strain-weighted shedding via `sigma_1_strain`/`sigma_2_strain`

### Cross-Cutting: Redundant Alive Mask

`agents['alive'].astype(bool)` and `np.where(agents['alive'])` are called in:
- `update_movement` (525,600×)
- `daily_disease_update` (17,520×)
- `daily_growth_and_aging` (21,900×)
- `daily_natural_mortality` (21,900×)
- `spawning_step` (15,420×)
- `settle_daily_cohorts` (1,472×)

**Total: ~604,000 alive-mask recomputations per simulation.** Each is a full-array scan + allocation.

### Cross-Cutting: Type Conversions

| Operation | Self-time (s) | Calls |
|---|---|---|
| `.astype()` | 6.47 | 1,216,196 |
| `np.sum` (reduce) | 9.21 | 1,166,444 |
| `list.append` | 1.90 | 26,103,956 |

The `astype` and `reduce` overhead is distributed across all modules.
The 26M `list.append` calls are from recording — suggests per-individual per-day recording.

---

## 3. Top 3 Optimization Targets

### Target 1: Movement — Batch Substeps (56.8% of runtime)

**Problem:** 525,600 calls to `update_movement`, each with per-call overhead
(alive mask, RNG, reflect, speed lookup). The CRW is embarrassingly parallelizable
across agents and substeps are independent given no inter-agent interaction
(gravity is rare).

**Proposed optimization: Fuse N substeps into one vectorized call**

```python
def daily_movement_batched(agents, habitat_side, sigma_turn, base_speed, substeps, rng):
    alive_mask = agents['alive'].astype(bool)
    alive_idx = np.where(alive_mask)[0]
    if len(alive_idx) == 0:
        return

    n = len(alive_idx)
    dt = (24.0 * 60.0) / substeps

    # Disease speed modifiers (constant within a day)
    speed_mods = SPEED_MODIFIER[agents['disease_state'][alive_idx]]
    eff_speed = base_speed * speed_mods  # (n,)

    # Generate ALL turning angles at once: (substeps, n)
    all_turns = rng.normal(0.0, sigma_turn, size=(substeps, n)).astype(np.float32)

    # Accumulate headings and displacements
    headings = agents['heading'][alive_idx].copy()
    x = agents['x'][alive_idx].copy()
    y = agents['y'][alive_idx].copy()

    for s in range(substeps):
        headings = (headings + all_turns[s]) % TWO_PI
        dx = eff_speed * np.cos(headings) * dt
        dy = eff_speed * np.sin(headings) * dt
        x = _reflect(x + dx, habitat_side)
        y = _reflect(y + dy, habitat_side)

    agents['heading'][alive_idx] = headings
    agents['x'][alive_idx] = x
    agents['y'][alive_idx] = y
    agents['speed'][alive_idx] = eff_speed
```

**What this eliminates per day per node:**
- 23 redundant `alive.astype(bool)` + `np.where` (keep 1, drop 23)
- 23 redundant `SPEED_MODIFIER` lookups
- 1 bulk RNG call replaces 24 individual calls
- Python loop overhead reduced from 24 function calls to 24 tight iterations

**Expected speedup on movement:** 2–3× (alive mask + dispatch overhead dominates at small N).
With 55.5s → ~20–28s. **Save 28–36s (29–37% of total).**

**If spawning gravity is disabled** (which it typically is), the entire inner loop
can be further fused with cumulative heading sums. If gravity IS active, it only
applies to a subset of days — detect and fall back to per-substep for those days.

**Implementation difficulty:** Moderate — the logic is straightforward but needs
careful handling of the spawning gravity branch and grid-based spatial disease
(`_diffuse_2d` / `InfectedDensityGrid`).

### Target 2: Spawning Cascade — Spatial Index (23.3% of runtime)

**Problem:** `_check_cascade_induction` computes O(T×I) pairwise distances
between cascade targets and recent spawners. Called 6,068 times with variable
T and I sizes. At large populations, this becomes quadratic.

**Also:** `_generate_larval_cohort` at 7.9s (2,311 calls × 3.4ms) — likely
dominated by per-offspring genotype assembly with NumPy structured array writes.

**Proposed optimizations:**

**(a) Replace pairwise distance with grid-based proximity check:**

Instead of computing all T×I distances, bin agents into grid cells of
`cascade_radius` side length. Only check targets in the same or adjacent
cells as inducers. For cascade_radius = ~50m in a 700m habitat, this is
a 14×14 grid — most cells empty.

```python
# Bin agents into grid cells
cell_size = cascade_radius
n_cells = int(np.ceil(habitat_side / cell_size))
target_cells = (target_positions / cell_size).astype(int)
inducer_cells = (inducer_positions / cell_size).astype(int)

# Build inducer lookup: cell → list of inducer indices
# Only check targets in adjacent cells (9 cells max per inducer)
```

This reduces O(T×I) → O(T × avg_inducers_per_neighborhood), typically O(T × ~5).

**(b) Skip cascade on zero-inducer days:**

15,420 calls to `_cascade_induction_step` but only 6,068 calls to
`_check_cascade_induction` — meaning 60% of cascade-step calls are
short-circuited early. The early exit is already there but the
function-call + mask overhead still costs 0.19ms × 9,352 = 1.8s.

**(c) Batch `_generate_larval_cohort` genotype writes:**

Profile suggests 3.4ms/call with structured array field-by-field writes.
Pre-allocate a contiguous buffer and write all offspring genotypes in one
vectorized operation.

**Expected speedup on spawning:** 2–3× (cascade distance is the main target).
22.8s → ~8–11s. **Save 12–15s (12–15% of total).**

**Implementation difficulty:** Moderate — grid binning is straightforward;
the tricky part is handling edge cases (agents near cell boundaries need
to check neighbor cells).

### Target 3: Shared Alive Mask + Disease Precomputation (12.2% + cross-cutting)

**Problem:** The alive mask (`agents['alive'].astype(bool)`) is recomputed
~604,000 times per simulation. Each call allocates a new boolean array and
scans all agents. Additionally, `daily_disease_update` recomputes
temperature-dependent rates (Arrhenius) for the same temperature 3× per day
(once per node), and counts compartments with separate `np.sum` calls.

**Proposed optimizations:**

**(a) Compute alive mask once per node per day, pass it down:**

```python
# In the daily loop of run_spatial_simulation:
for node_idx in range(K):
    node_agents = agents_per_node[node_idx]
    alive_mask = node_agents['alive'].astype(bool)
    alive_idx = np.where(alive_mask)[0]

    # Pass to ALL daily functions:
    daily_movement(node_agents, ..., alive_mask=alive_mask, alive_idx=alive_idx)
    daily_disease_update(node_agents, ..., alive_mask=alive_mask)
    daily_growth_and_aging(node_agents, ..., alive_mask=alive_mask)
    daily_natural_mortality(node_agents, ..., alive_mask=alive_mask)
```

This eliminates ~600K redundant array allocations + scans.

**(b) Cache Arrhenius rates per node per day:**

Disease calls `arrhenius()` 220,009 times. Many of these are for the same
(rate, Ea, T) tuple within a single day at a single node. Precompute all
5 temperature-dependent rates once per node-day:

```python
rates = {
    'mu_EI1': arrhenius(cfg.mu_EI1_ref, cfg.Ea_EI1, T),
    'mu_I1I2': arrhenius(cfg.mu_I1I2_ref, cfg.Ea_I1I2, T),
    'mu_I2D': arrhenius(cfg.mu_I2D_ref, cfg.Ea_I2D, T),
    'rho_rec': arrhenius(cfg.rho_rec_ref, cfg.Ea_rec, T),
    'growth_rate': arrhenius(cfg.growth_ref, cfg.Ea_growth, T),
}
```

5 calls replace 220,009. Saves ~0.18s directly + eliminates repeated exp() overhead deeper in the call chain.

**(c) Vectorize compartment counts:**

Replace separate `np.sum(alive & (ds == X))` for S, E, I1, I2, D, R with:

```python
counts = np.bincount(ds[alive_mask], minlength=6)
n_S, n_E, n_I1, n_I2, n_D, n_R = counts
```

One pass instead of 6.

**Expected speedup on disease:** 1.5–2× (11.9s → ~6–8s). **Save 4–6s.**
**Cross-cutting alive mask savings:** ~1–2s additional from growth/mortality/spawning.
**Total save: 5–8s (5–8% of total).**

**Implementation difficulty:** Trivial (alive mask passthrough) to Moderate
(Arrhenius caching requires refactoring function signatures).

---

## 4. Priority Summary

| Priority | Target | Current (s) | Expected After (s) | Saved (s) | % Speedup | Difficulty |
|---|---|---|---|---|---|---|
| 🥇 1 | Movement: batch substeps | 55.5 | 20–28 | 28–36 | 29–37% | Moderate |
| 🥈 2 | Spawning: grid-based cascade | 22.8 | 8–11 | 12–15 | 12–15% | Moderate |
| 🥉 3 | Shared alive mask + disease precompute | 11.9 + 3.4 | 8–10 | 5–8 | 5–8% | Trivial–Moderate |

### Combined Expected Result

**Current runtime:** ~98s (3 nodes, 15K agents, 20 years)

**After all 3 optimizations:** ~36–51s → **1.9–2.7× overall speedup**

**Best case (aggressive estimates):** 36s (2.7×)  
**Conservative case:** 51s (1.9×)

---

## 5. What NOT to Optimize

| Component | Time | Why Skip |
|---|---|---|
| daily_growth | 1.8s (1.8%) | Already fast, batch-every-7-days would save <1s and risk drift |
| daily_mortality | 1.6s (1.7%) | Already fast, simple operations |
| environment | 0.05s (<0.1%) | Negligible — precomputing SST saves nothing meaningful |
| pathogen_dispersal | 0.04s (<0.1%) | Negligible |
| settlement | 0.18s (0.2%) | Negligible |

---

## 6. Scaling Implications

At full coastline scale (150 nodes, 75K agents), these optimizations become
critical:

- **Movement** scales as O(N × substeps × K), so 150 nodes × 10× agents = 1500×
  current movement time without optimization. Batching substeps is essential.
- **Cascade induction** scales as O(T × I) per node. With larger populations
  per node, quadratic pairwise distance becomes the bottleneck. Grid indexing
  is mandatory at scale.
- **Alive mask** overhead scales linearly with total agents. At 75K agents,
  600K redundant scans × 75K array → significant memory bandwidth waste.

---

## 7. Implementation Order

1. **Movement batching** (biggest win, no API changes to other modules)
2. **Shared alive mask** (trivial, touches all modules but simple signature change)
3. **Spawning grid index** (most complex, but second-biggest win)
4. **Disease Arrhenius cache** (small win, pairs naturally with alive mask work)

Estimated implementation time: 2–3 hours for all four changes.
