# Dynamic P_env Decay Spec

## Problem

W42 nails Alaska (AK-PWS=50.4%) but southern regions recover too much:
- SS-S = 22.7% (target 5%)
- JDF = 24.4% (target 2%)

After the crash, populations rebuild through evolutionary rescue (resistant alleles).
Currently P_env is static — the environmental vibrio reservoir stays the same
regardless of how many infected hosts remain. This means even after a crash,
disease pressure stays high... but evolved resistance lets populations bounce back.

## Solution: Dynamic P_env with Host-Amplification Feedback

Replace the static `P_env_max` with a dynamic per-node `P_env_pool` that:

1. **Builds up** from local shedding (infected animals amplify the reservoir)
2. **Decays** at rate `delta_env` (protist grazing, UV, dilution)
3. **Has a floor** from community maintenance: `P_env_floor × VBNC(SST)`

This creates differential recovery:
- **Warm south**: High community floor → P_env stays suppressive → populations stay crushed
- **Cool north**: Low community floor → P_env decays in winter → populations recover

## Current Mechanics (what to change)

In `update_vibrio_concentration()`:

```python
# Current: static P_env
env = environmental_vibrio(T_celsius, salinity, cfg)  # = P_env_max * VBNC * thermal * sal
dP = shed - decay - flush + env + dispersal_input
```

**New approach**: Split `environmental_vibrio()` into two components:

```python
# Community floor (always present, SST-modulated)
floor = P_env_floor * vbnc_activation * g_peak * sal_mod

# Dynamic host-amplified pool (tracked per node)
# Pool gets contributions from shedding, decays over time
pool_input = alpha_env * shed  # fraction of shedding that enters environmental pool
pool_decay = delta_env * P_env_pool  # natural decay of environmental pool

# Update pool
P_env_pool_new = max(0.0, P_env_pool + (pool_input - pool_decay) * dt)

# Total environmental vibrio contribution
env = floor + P_env_pool
```

## Config Parameters (in DiseaseSection)

```python
# Dynamic P_env parameters
P_env_dynamic: bool = False          # Enable dynamic P_env (False = backward compat)
P_env_floor: float = 50.0            # Community-maintained vibrio floor (bact/mL/d)
alpha_env: float = 0.1               # Fraction of shedding entering environmental pool
delta_env: float = 0.05              # Environmental pool decay rate (d⁻¹) → ~20 day half-life
```

When `P_env_dynamic=False`: behavior identical to current (uses static `P_env_max`).
When `P_env_dynamic=True`: uses floor + dynamic pool instead of P_env_max.

## State Changes

### NodeDiseaseState
Add field:
```python
P_env_pool: float = 0.0  # Dynamic environmental vibrio pool (bact/mL)
```

### update_vibrio_concentration()
Add parameter `P_env_pool: float = 0.0` and return tuple `(P_new, P_env_pool_new)` when dynamic is enabled.

**IMPORTANT**: To minimize disruption, keep the return type as `float` when `P_env_dynamic=False`.
When `P_env_dynamic=True`, return a tuple `(P_new, P_env_pool_new)`.

Actually, cleaner approach: always return `float` for P_new. Update `P_env_pool` on the NodeDiseaseState directly inside `daily_disease_update()`, not in `update_vibrio_concentration()`. This keeps the vibrio update function's signature unchanged.

### Revised approach for minimal disruption:

1. In `daily_disease_update()`, BEFORE calling `update_vibrio_concentration()`:
   - Compute `pool_input = alpha_env * total_shedding`
   - Compute `pool_decay = delta_env * node_state.P_env_pool`
   - Update `node_state.P_env_pool = max(0.0, node_state.P_env_pool + (pool_input - pool_decay) * dt)`

2. In `environmental_vibrio()`, add optional `P_env_pool` parameter:
   - If `cfg.P_env_dynamic` and `P_env_pool is not None`:
     - Return `cfg.P_env_floor * vbnc_activation * g_peak * sal_mod + P_env_pool`
   - Else: return current behavior (`cfg.P_env_max * ...`)

3. Pass `node_state.P_env_pool` to the `environmental_vibrio()` call inside `update_vibrio_concentration()`.

Wait, that's still messy. Let me think about the cleanest way...

### Cleanest approach:

**In `update_vibrio_concentration()`**, add `P_env_pool` parameter (default 0.0):

```python
def update_vibrio_concentration(
    P_k, n_I1, n_I2, n_D_fresh, T_celsius, salinity, phi_k,
    dispersal_input, cfg, dt=1.0, override_shedding=None,
    disease_reached=True, P_env_pool=0.0,
) -> float:
    ...
    if disease_reached:
        if cfg.P_env_dynamic:
            # Floor component (community maintenance, SST-modulated)
            k = getattr(cfg, 'k_vbnc', K_VBNC)
            vbnc_activation = 1.0 / (1.0 + np.exp(-k * (T_celsius - cfg.T_vbnc)))
            g_peak = thermal_performance(3000.0, T_celsius, rate_ref=1.0)
            sal_mod = salinity_modifier(salinity, cfg.s_min, cfg.s_full)
            floor = cfg.P_env_floor * vbnc_activation * g_peak * sal_mod
            env = floor + P_env_pool
        else:
            env = environmental_vibrio(T_celsius, salinity, cfg)
    ...
```

**In `daily_disease_update()`**, after computing shedding and BEFORE calling `update_vibrio_concentration()`:

```python
if getattr(cfg, 'P_env_dynamic', False):
    # Compute shedding for pool update
    total_shed = (shedding_rate_I1(T_celsius, cfg) * n_I1
                  + shedding_rate_I2(T_celsius, cfg) * n_I2
                  + cfg.sigma_D * n_D_fresh)
    pool_input = cfg.alpha_env * total_shed
    pool_decay = cfg.delta_env * node_state.P_env_pool
    node_state.P_env_pool = max(0.0, node_state.P_env_pool + (pool_input - pool_decay))
```

Then pass `P_env_pool=node_state.P_env_pool` to `update_vibrio_concentration()`.

## Files to Change

1. **`sswd_evoepi/config.py`** — Add 4 fields to `DiseaseSection`
2. **`sswd_evoepi/disease.py`** — Modify `update_vibrio_concentration()`, `daily_disease_update()`, `NodeDiseaseState`
3. **`tests/test_dynamic_penv.py`** — New test file

## Tests

1. `test_backward_compat` — P_env_dynamic=False gives identical results to current
2. `test_pool_builds_from_shedding` — With infected hosts, P_env_pool increases
3. `test_pool_decays_without_hosts` — Without infected hosts, P_env_pool decays exponentially
4. `test_floor_provides_minimum` — Even with pool=0, floor still provides disease pressure
5. `test_floor_modulated_by_sst` — Floor is higher at warm SST, lower at cold SST
6. `test_warm_vs_cold_differential` — Warm node retains more disease pressure than cold node after crash
7. `test_dynamic_penv_config_validation` — P_env_floor >= 0, alpha_env >= 0, delta_env > 0

## Backward Compatibility

- `P_env_dynamic` defaults to `False` → zero behavior change for all existing configs
- `P_env_max` still used when `P_env_dynamic=False`
- NodeDiseaseState.P_env_pool defaults to 0.0 → no effect when not enabled
- All existing tests should pass unchanged
