# Spawning System Overhaul — Specification

**Date:** 2026-02-15
**Authors:** Willem Weertman & Anton 🔬
**Status:** APPROVED — ready for implementation

---

## Summary

Replace the single-day pulse spawning model with a biologically realistic system incorporating:
1. Extended spawning season (November–July, ~270 days)
2. Cascading spawning bouts triggered by social chemical cues
3. Sex-asymmetric induction (females strongly trigger males, males moderately trigger females)
4. Male multi-bout spawning (2–3 per season) vs female single spawning
5. Pre-spawning aggregation movement ("spawning gravity")
6. Post-spawning immunosuppression (~2× disease susceptibility for 3–4 weeks)

---

## New Parameters

### Spawning Season (config: `population.spawning`)

| Parameter | Symbol | Value | Units | Confidence | Description |
|-----------|--------|-------|-------|------------|-------------|
| `season_start_doy` | — | 305 | day of year | ★★★☆☆ | Season start (~Nov 1) |
| `season_end_doy` | — | 196 | day of year | ★★★☆☆ | Season end (~Jul 15) |
| `peak_doy` | — | 105 | day of year | ★★☆☆☆ | Peak spawning activity (~Apr 15) |
| `peak_width_days` | σ_peak | 45 | days | ★★☆☆☆ | Std dev of seasonal peak (Normal) |
| `lat_shift_per_deg` | — | 3.0 | days/°N | ★★☆☆☆ | Latitude shift of peak (retained) |

**Notes:** Season wraps across year boundary. `peak_doy` and `peak_width_days` define the probability distribution for spontaneous spawning readiness. Latitude shifts the peak but not the season boundaries.

### Spontaneous Spawning Rates (config: `population.spawning`)

| Parameter | Symbol | Value | Units | Confidence | Description |
|-----------|--------|-------|-------|------------|-------------|
| `p_spontaneous_female` | p_sf | 0.005 | d⁻¹ | ★★☆☆☆ | Daily probability a ready female spawns spontaneously |
| `p_spontaneous_male` | p_sm | 0.008 | d⁻¹ | ★★☆☆☆ | Daily probability a ready male initiates a bout spontaneously |

**Notes:** These are base rates modulated by seasonal readiness (see Spawning Readiness below). At peak season, effective rates are p_s × seasonal_modifier(doy). Off-peak, rates drop toward zero.

### Cascade Induction (config: `population.spawning`)

| Parameter | Symbol | Value | Units | Confidence | Description |
|-----------|--------|-------|-------|------------|-------------|
| `induction_female_to_male` | κ_fm | 0.80 | — | ★★★☆☆ | Probability a male spawns when a female nearby has spawned in last `cascade_window` days |
| `induction_male_to_female` | κ_mf | 0.30 | — | ★★☆☆☆ | Probability a female spawns when a male nearby has spawned in last `cascade_window` days |
| `cascade_window` | τ_cascade | 3 | days | ★★☆☆☆ | Duration of chemical spawning cue persistence |
| `cascade_radius` | r_cascade | 50.0 | m | ★★☆☆☆ | Effective range of chemical spawning cue |

**Notes:** "Nearby" = within `cascade_radius` meters (Euclidean distance in node-local x,y space). κ_fm is high because females very strongly induce males. κ_mf is moderate because males moderately trigger females. Both sexes CAN spawn spontaneously regardless.

### Male Multi-Bout (config: `population.spawning`)

| Parameter | Symbol | Value | Units | Confidence | Description |
|-----------|--------|-------|-------|------------|-------------|
| `male_max_bouts` | B_max | 3 | — | ★★★☆☆ | Maximum spawning bouts per male per season |
| `male_refractory_days` | τ_refrac | 21 | days | ★★☆☆☆ | Minimum days between male spawning bouts |

**Notes:** After each bout, male enters refractory period. Males contribute sperm to the SRS lottery each time they spawn. Females spawn exactly once per season (their full egg mass).

### Spawning Aggregation / Gravity (config: `movement.spawning_gravity`)

| Parameter | Symbol | Value | Units | Confidence | Description |
|-----------|--------|-------|-------|------------|-------------|
| `gravity_enabled` | — | true | bool | — | Enable pre-spawning aggregation movement |
| `gravity_strength` | G_s | 0.3 | m/min | ★★☆☆☆ | Maximum speed bias toward nearest conspecifics |
| `gravity_range` | R_grav | 100.0 | m | ★★☆☆☆ | Sensory detection range for conspecifics |
| `pre_spawn_gravity_days` | τ_pre | 14 | days | ★★☆☆☆ | Days before spawning readiness that gravity activates |
| `post_spawn_gravity_days` | τ_post | 14 | days | ★★☆☆☆ | Days after spawning that gravity persists |

**Notes:** Gravity modifies the CRW heading with a bias vector pointing toward the center-of-mass of conspecifics within `gravity_range`. Strength ramps linearly from 0 at τ_pre to G_s at spawning, then decays linearly back to 0 over τ_post. This naturally creates aggregation clumps without overriding the CRW entirely.

### Post-Spawning Immunosuppression (config: `disease.spawning_immunity`)

| Parameter | Symbol | Value | Units | Confidence | Description |
|-----------|--------|-------|-------|------------|-------------|
| `immunosuppression_enabled` | — | true | bool | — | Enable post-spawning susceptibility increase |
| `susceptibility_multiplier` | ψ_spawn | 2.0 | — | ★★★☆☆ | Force-of-infection multiplier during immunosuppression |
| `immunosuppression_duration` | τ_immuno | 28 | days | ★★☆☆☆ | Duration of post-spawning immunosuppression |

**Notes:** After each spawning event, the individual's effective resistance is reduced: r_eff = r_i / ψ_spawn for τ_immuno days. For males spawning multiple bouts, the timer resets each time — they spend more total time immunosuppressed. This stacks with the density-dependent transmission increase from physical aggregation.

---

## New Agent Fields

Add to `AGENT_DTYPE` in `types.py`:

| Field | dtype | Bytes | Description |
|-------|-------|-------|-------------|
| `spawning_ready` | int8 | 1 | 0=not ready, 1=ready to spawn this season |
| `has_spawned` | int8 | 1 | Females: 0/1. Males: count of bouts this season (0–3) |
| `spawn_refractory` | int16 | 2 | Days remaining in male refractory period (countdown) |
| `spawn_gravity_timer` | int16 | 2 | Days remaining in gravity phase (countdown, pre+post) |
| `immunosuppression_timer` | int16 | 2 | Days remaining in post-spawning immunosuppression |
| `last_spawn_day` | int16 | 2 | Day-of-year of last spawning event (for cascade tracking) |

**Total added:** 10 bytes per agent (AGENT_DTYPE grows from ~41 to ~51 bytes).

---

## Spawning Readiness Model

Each adult enters "spawning ready" state based on a probability distribution across the season:

```
readiness_prob(doy) = seasonal_envelope(doy) × size_maturity(L)
```

Where:
- `seasonal_envelope(doy)` = Normal PDF centered on `peak_doy` (shifted by latitude) with σ = `peak_width_days`, evaluated daily and normalized so max = 1.0
- `size_maturity(L)` = 1.0 if L ≥ L_min_repro, else 0.0

At the start of each day during the season, each non-ready adult rolls against readiness_prob. Once ready, they remain ready for the rest of the season (or until they spawn, for females).

---

## Daily Spawning Step (pseudocode)

```python
def spawning_step(node, day_of_year, agents, ...):
    if not in_spawning_season(day_of_year):
        return []
    
    # 1. Update readiness
    for each adult not yet ready:
        if random() < readiness_prob(day_of_year, lat):
            agent.spawning_ready = 1
    
    # 2. Tick down refractory timers
    for males with spawn_refractory > 0:
        spawn_refractory -= 1
    
    # 3. Spontaneous spawning attempts
    new_spawners_today = []
    for each ready female (has_spawned == 0):
        if random() < p_spontaneous_female:
            spawn(female)
            new_spawners_today.append(female)
    
    for each ready male (has_spawned < B_max, spawn_refractory == 0):
        if random() < p_spontaneous_male:
            spawn(male)
            new_spawners_today.append(male)
    
    # 4. Cascade induction (iterate for cascade_window recent spawners)
    recent_female_spawners = females who spawned in last cascade_window days
    recent_male_spawners = males who spawned in last cascade_window days
    
    for each ready male (has_spawned < B_max, refractory == 0):
        if any recent_female_spawner within cascade_radius:
            if random() < κ_fm:
                spawn(male)
    
    for each ready female (has_spawned == 0):
        if any recent_male_spawner within cascade_radius:
            if random() < κ_mf:
                spawn(female)
    
    # 5. Collect spawned gametes for SRS lottery
    # ... produce larval cohort from today's spawners
```

---

## Movement Modification (gravity)

In `movement.py`, the CRW heading update becomes:

```python
# Existing CRW
new_heading = heading + Normal(0, sigma_turn)

# Spawning gravity bias (if agent has gravity_timer > 0)
if agent.spawn_gravity_timer > 0:
    # Find center of mass of conspecifics within R_grav
    neighbors = agents_within_radius(agent, R_grav)
    if len(neighbors) > 0:
        com = center_of_mass(neighbors)
        angle_to_com = atan2(com.y - agent.y, com.x - agent.x)
        
        # Blend: gravity_fraction ramps up then down
        gravity_frac = gravity_strength_at_time(agent) / base_speed
        new_heading = circular_blend(new_heading, angle_to_com, gravity_frac)
```

---

## Disease Module Modification

In the force-of-infection calculation, modify individual susceptibility:

```python
# Current: lambda_i = a × P/(K_half+P) × (1 - r_i) × S_sal × f_size
# New:     lambda_i = a × P/(K_half+P) × (1 - r_eff_i) × S_sal × f_size
# where r_eff_i = r_i / psi_spawn if immunosuppression_timer > 0, else r_i
```

---

## Implementation Phases

### Phase 1: Agent fields + season window + spontaneous spawning
- Add new fields to AGENT_DTYPE
- Add SpawningSection to config
- Implement spawning season check (wrapping Nov–Jul)
- Implement readiness probability
- Implement spontaneous spawning (no cascade yet)
- Integrate daily spawning step into model loop
- Tests: season boundaries, readiness distribution, total fecundity normalization

### Phase 2: Cascade induction + male multi-bout
- Implement cascade triggering with spatial proximity
- Implement male refractory period and bout counting
- Implement sex-asymmetric induction strengths
- Tests: cascade propagation at different densities, male bout counts, female single-spawn enforcement

### Phase 3: Spawning gravity (movement modification)
- Implement gravity bias in CRW
- Implement pre/post spawning gravity timers
- Tests: clump formation, density increase during aggregation, dispersal after spawning

### Phase 4: Post-spawning immunosuppression
- Implement immunosuppression timer
- Modify force-of-infection to use r_eff
- Tests: susceptibility increase, timer reset for males, disease rate comparison

### Phase 5: Integration test + validation
- Full 5-node 20-year simulation with new spawning system
- Compare against old pulse model results
- Verify: total annual recruitment similar, but timing and composition differ
- Verify: post-spawning disease spike emerges naturally
- Verify: low-density cascade failure (additional Allee effect)

---

## Validation Targets

| Metric | Expected Behavior |
|--------|-------------------|
| Total annual recruitment | Within 20% of pulse model at equilibrium |
| Spawning bouts per node per season | 2–5 major bouts at healthy density |
| Male bout count | Mean ~2.2 per male per season |
| Ne/N | Similar to current (~10⁻³) — SRS dominates |
| Post-spawning disease incidence | ~2× higher than pre-spawning baseline |
| Low-density cascade failure | < 50 adults → sporadic individual spawning, no bouts |
| Aggregation clump size | 5–20 individuals within 50m during bouts |

---

*Specification by Anton 🔬, based on Willem Weertman's biological observations*
*University of Washington, Sea Star Lab, Friday Harbor Laboratories*
