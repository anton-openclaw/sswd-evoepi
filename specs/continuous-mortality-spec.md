# Continuous Mortality & Daily Growth — Specification

**Date:** 2026-02-19
**Authors:** Willem Weertman & Anton 🔬
**Status:** APPROVED — implementation in progress

---

## Problem

Natural mortality and growth are applied as year-boundary lump sums in Phase B. This creates:
1. **Artificial population drops** at year boundaries (e.g., 237→55 in one timestep)
2. **Sawtooth artifacts** in population trajectories — same problem continuous settlement fixed for epidemic curves
3. **Stage transition discontinuities** — all aging/growth happens in one pulse
4. **Interaction artifacts with daily disease** — disease kills continuously but natural mortality competes only once/year

## Solution

Move all demographic processes into the daily loop:
- **Natural mortality**: daily probability derived from annual survival rates
- **Growth**: daily VB increment (dt=1/365)
- **Aging**: daily (1/365 year per day)
- **Stage transitions**: checked daily after growth
- **Senescence**: daily probability check

Eliminate Phase B mortality/growth entirely. Phase B retains only reproduction/settlement cleanup.

---

## Mathematical Conversion

### Daily Mortality Probability

Annual survival rate `S_annual` → daily survival → daily death probability:

```
S_daily = S_annual^(1/365)
p_death_daily = 1 - S_daily = 1 - S_annual^(1/365)
```

For senescence (age > senescence_age):
```
annual_mort = base_mort + senescence_mortality * (age - senescence_age) / 20
daily_mort = 1 - (1 - annual_mort)^(1/365)
```

**Verification**: Over 365 days, compound daily survival = S_daily^365 = S_annual ✓

### Daily Growth

Von Bertalanffy with dt=1/365:
```
L(t+dt) = L_inf - (L_inf - L(t)) * exp(-k_growth * dt)
```
where `dt = 1/365` years.

Growth noise: scale σ proportionally (σ_daily = σ_annual / √365).

### Daily Aging

```
age += 1/365  (each day)
```

No annual increment anymore — age advances continuously.

---

## Implementation Plan

### New Functions

#### `daily_natural_mortality(agents, pop_cfg, rng) -> (n_killed, n_senescence)`
- Vectorized: compute daily death probability for ALL alive agents at once
- `annual_surv = pop_cfg.annual_survival[agents['stage']]`  (vectorized stage lookup)
- `daily_surv = annual_surv ** (1/365)`
- Senescence overlay for age > senescence_age
- Single vectorized random draw: `rng.random(n_alive) < p_death_daily`
- Stamp cause_of_death (NATURAL=2 or SENESCENCE=3)

#### `daily_growth_and_aging(agents, pop_cfg, rng) -> None`
- Increment age: `agents['age'][alive] += 1.0/365.0`
- VB growth with dt=1/365: `new_size = L_inf - (L_inf - old_size) * exp(-k * 1/365)`
- Add noise: `new_size += rng.normal(0, sigma_daily, n_alive)` where σ_daily = 2.0/√365 ≈ 0.105mm
- Stage assignment: vectorized via size thresholds
- All vectorized — no Python loops

### Changes to Coupled Sim (`run_coupled_simulation`)

**In daily loop (Phase A)**, add after disease step, before settlement:
```python
# Daily demographics
with perf.track("daily_mortality"):
    n_mort, n_senes = daily_natural_mortality(agents, pop_cfg, rng)
    daily_nat_deaths_accum += n_mort
    daily_senes_deaths_accum += n_senes

with perf.track("daily_growth"):
    daily_growth_and_aging(agents, pop_cfg, rng)
```

**Remove from Phase B:**
- Delete `annual_natural_mortality(agents, pop_cfg, rng)` call
- Delete `annual_growth_and_aging(agents, pop_cfg, rng)` call
- Keep yearly accumulators: `yearly_natural_deaths[year] = daily_nat_deaths_accum`

### Changes to Spatial Sim (`run_spatial_simulation`)

Same pattern — move mortality/growth into the per-node daily loop.

### Backward Compatibility

- `annual_natural_mortality()` and `annual_growth_and_aging()` functions kept (not deleted)
  — they may be used by tests or other code paths
- New daily functions are the default path in both sims

---

## Vectorization Strategy

The old `annual_natural_mortality` uses a Python `for idx in alive_idx:` loop. The new daily version MUST be vectorized or it will be 365× slower.

```python
def daily_natural_mortality(agents, pop_cfg, rng):
    alive = np.where(agents['alive'])[0]
    if len(alive) == 0:
        return 0, 0
    
    stages = agents['stage'][alive]
    ages = agents['age'][alive]
    
    # Vectorized annual survival lookup
    annual_surv = np.array(pop_cfg.annual_survival, dtype=np.float64)
    # Clip stages to valid range
    s_clipped = np.clip(stages, 0, len(annual_surv) - 1)
    base_annual_mort = 1.0 - annual_surv[s_clipped]
    
    # Senescence overlay
    total_annual_mort = base_annual_mort.copy()
    senes_mask = ages > pop_cfg.senescence_age
    if np.any(senes_mask):
        extra = pop_cfg.senescence_mortality * (ages[senes_mask] - pop_cfg.senescence_age) / 20.0
        total_annual_mort[senes_mask] = np.minimum(1.0, base_annual_mort[senes_mask] + extra)
    
    # Convert to daily
    daily_mort = 1.0 - (1.0 - total_annual_mort) ** (1.0/365.0)
    
    # Roll
    rolls = rng.random(len(alive))
    dies = rolls < daily_mort
    
    # Apply
    dead_idx = alive[dies]
    agents['alive'][dead_idx] = False
    
    # Stamp cause of death
    dead_ages = agents['age'][dead_idx]
    senes_dead = dead_ages > pop_cfg.senescence_age
    agents['cause_of_death'][dead_idx[senes_dead]] = 3   # SENESCENCE
    agents['cause_of_death'][dead_idx[~senes_dead]] = 2  # NATURAL
    
    return int(np.sum(dies)), int(np.sum(senes_dead))
```

Growth is similarly vectorized — no loops.

---

## Metrics Impact

- `yearly_natural_deaths` and `yearly_senescence_deaths`: now accumulated daily (same semantics)
- `daily_pop` timeseries: will be MUCH smoother (no year-boundary drops)
- **Expected**: population curves become gradual declines/increases rather than staircases

---

## Testing

1. **Statistical equivalence**: Run 100 seeds with daily vs annual mortality. Mean population at year 20 should match within 5% (same expected deaths, different timing).
2. **No year-boundary discontinuity**: max(|pop[day] - pop[day-1]|) should never exceed ~5% of pop (disease days excepted).
3. **Stage transitions work daily**: agents progress through stages smoothly.
4. **Senescence works daily**: old agents die at correct elevated rate.
5. **Cause-of-death accounting**: all deaths stamped correctly.
6. **Performance**: daily vectorized mortality should be ≤ 2ms/day overhead (vs old ~5ms/year loop).

---

## Performance Notes

Old annual mortality: Python loop over all alive agents, once per year.
New daily mortality: NumPy vectorized over all alive agents, 365× per year.

The vectorized version processes ~15,000 agents in <0.5ms per call.
365 × 0.5ms = 182ms/year — acceptable (old loop was ~5ms but only 1×/year = 5ms).
Net overhead: ~177ms/year × 20 years = ~3.5s. Acceptable for 20-50s runtime.

Growth overhead is similar. Total expected: +5-8s per 20-year run.
