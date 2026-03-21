# FIX: Daily Larval Dispersal Pipeline

## Bug Summary

The daily spawning system produces larval cohorts and stores them in
`pending_cohorts[SOURCE_NODE]`. When PLD elapses, they are settled at
the **source node** via `settle_daily_cohorts()` — bypassing both the
**C connectivity matrix** (no spatial dispersal) and **Beverton-Holt
density-dependent recruitment**. This makes every site an independent
population with unchecked self-recruitment. The C matrix is dead code
for ~80% of larvae.

## Root Cause

In `model.py` (spatial simulation loop):
1. **Line ~2490**: `_insort_cohort(pending_cohorts[i], c)` stores cohorts at source `i`
2. **Line ~2439**: Daily settlement pops from `pending_cohorts[i]` and settles at node `i`
3. **`settle_daily_cohorts()`**: Docstring claims "BH applied upstream" but no caller applies BH
4. The annual C-matrix step (line ~2555) only catches leftovers whose PLD extends past year-end

## Required Fix

### 1. Distribute larvae through C at PRODUCTION time (not settlement time)

When `spawning_step()` produces a cohort at source node `i`:

```python
# BEFORE (broken):
for c in cohorts_today:
    c.source_node = nd.node_id
    _insort_cohort(pending_cohorts[i], c)

# AFTER (fixed):
for c in cohorts_today:
    c.source_node = nd.node_id
    # Distribute through C matrix immediately
    # Use distribute_larvae_counts for efficiency (count-based, no genotype arrays)
    dest_map = distribute_larvae_counts(
        [i], [c.n_competent], network.C, rng
    )
    for dest_id, arrivals in dest_map.items():
        for (n_arriving, src_id) in arrivals:
            if n_arriving <= 0:
                continue
            dest_cohort = LarvalCohort(
                source_node=src_id,
                n_competent=n_arriving,
                genotypes=None,  # genotypes created at settlement via SRS from source parents
                parent_pairs=None,
                pld_days=c.pld_days,
                spawn_day=c.spawn_day,
                sst_at_spawn=c.sst_at_spawn,
            )
            _insort_cohort(pending_cohorts[dest_id], dest_cohort)
```

This means:
- Larvae are immediately dispersed through C (binomial survival + multinomial destination)
- ~98% of larvae are lost (C row sum ≈ 0.02) — matching biological reality
- Surviving larvae are stored at their DESTINATION node
- PLD timer starts from spawn day (they drift during PLD, arrive at destination)

### 2. Apply Beverton-Holt at daily settlement

When PLD elapses and cohorts are settled:

```python
# BEFORE (broken):
ready = _pop_ready_cohorts(pending_cohorts[i], sim_day)
if ready:
    n_settled = settle_daily_cohorts(ready, node.agents, ...)

# AFTER (fixed):
ready = _pop_ready_cohorts(pending_cohorts[i], sim_day)
if ready:
    # Apply BH to each cohort before settlement
    bh_filtered = []
    for c in ready:
        current_alive = int(np.sum(node.agents['alive']))
        available_slots = max(0, nd.carrying_capacity - current_alive)
        if available_slots == 0:
            break
        n_adults = int(np.sum(
            node.agents['alive'] & (node.agents['stage'] == Stage.ADULT)
        ))
        cue_mod = settlement_cue_modifier(n_adults)
        effective_settlers = max(0, int(c.n_competent * cue_mod))
        if effective_settlers <= 0:
            continue
        n_recruits = beverton_holt_recruitment(
            effective_settlers, nd.carrying_capacity, pop_cfg.settler_survival
        )
        n_recruits = min(n_recruits, available_slots)
        if n_recruits <= 0:
            continue
        # Create genotypes via SRS from source parents
        src_node = network.nodes[c.source_node]
        # ... SRS lottery from source node's current adults ...
        # Create new cohort with n_recruits and proper genotypes
        bh_filtered.append(filtered_cohort)
    if bh_filtered:
        n_settled = settle_daily_cohorts(bh_filtered, node.agents, ...)
```

### 3. Handle genotypes correctly

Since cohorts are now distributed as COUNTS (no genotype arrays), genotypes must
be created at settlement time via SRS lottery from the source node's adults. This
matches the annual C-matrix pathway pattern (lines 2647-2680 in model.py).

The pattern is:
1. At dispersal: only counts travel through C (efficient)
2. At settlement: SRS from source node's current breeding adults produces genotypes
3. Only BH-filtered recruits get genotypes (no wasted computation)

### 4. Remove or simplify annual C-matrix step

The annual step (line ~2555) becomes redundant since all larvae now go through C
at production time. Options:
- **Remove entirely** (clean, since daily system handles everything)
- **Keep as safety net** for any edge cases (cohorts spanning year boundary)

Recommendation: **Remove it.** The daily system with C dispersal handles all cases.
Any cohort still in pending_cohorts at year-end will settle in the next year's daily
loop when its PLD elapses. No need for annual redistribution.

### 5. Update `settle_daily_cohorts` docstring and behavior

- Remove the misleading "BH applied upstream" comment
- `settle_daily_cohorts` should ONLY do slot placement (no filtering)
- All filtering (C, BH, cue) happens in the caller
- OR: add BH inside `settle_daily_cohorts` itself (simpler, self-contained)

**Recommendation**: Add BH inside `settle_daily_cohorts` since it already has
`pop_cfg` (which contains `settler_survival`) and `carrying_capacity`. This makes
the function self-contained and correct regardless of caller.

## Implementation Notes

### Performance considerations
- `distribute_larvae_counts` does binomial + multinomial per source node per day
- With 896 nodes spawning daily for ~200 days = ~180K calls per year
- Each call iterates N=896 destinations — but most get 0 larvae
- The C matrix is sparse in practice; consider sparse representation
- Profile before optimizing; may be fine as-is

### Optimization: batch daily dispersal
Instead of calling `distribute_larvae_counts` per-cohort, batch all same-day cohorts:
```python
# Collect all day's spawning across nodes
day_sources = []
day_counts = []
for i, node in enumerate(network.nodes):
    if cohorts_today_for_node_i:
        day_sources.append(i)
        day_counts.append(sum(c.n_competent for c in cohorts_today_for_node_i))
if day_sources:
    dest_map = distribute_larvae_counts(day_sources, day_counts, network.C, rng)
    # Store at destinations
    for dest_id, arrivals in dest_map.items():
        for (n_arriving, src_id) in arrivals:
            ...
```
This does ONE multinomial draw per source per day instead of one per cohort.

### Test impact
- Existing recruitment tests may break (they likely assumed self-recruitment)
- Need to verify C matrix is loaded and accessible during daily loop
- Settlement numbers will drop dramatically (98% larval loss via C)
- Recalibration will be needed

## Files to modify
1. **`sswd_evoepi/model.py`** — main fix: daily dispersal + BH in settlement
2. **`sswd_evoepi/model.py`** — remove/simplify annual C-matrix step
3. **`sswd_evoepi/model.py`** — update `settle_daily_cohorts` (add BH or update docstring)
4. **Tests** — update/add tests for daily C-matrix dispersal

## Verification
After fix, verify:
1. Larvae from node j appear at node k (k ≠ j) — spatial connectivity works
2. Total settlement is ~2% of production (matching C row sums)
3. BH limits recruitment at high settler density
4. Self-recruitment fraction matches C[j,j] diagonal
5. PLD timing still works (larvae settle after PLD days at destination)
6. Year-boundary cohorts still work (spawned Dec, settle Feb next year)
