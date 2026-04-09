# Design Spec: Load-Dependent Disease Model (Tier 2)

## Motivation

Recommendation from Alyssa Gehman: the model should have LD50-like dynamics
where pathogen dose/load determines disease outcome, not just binary
infection status. Currently, once infected, all agents follow the same
timer-based trajectory regardless of exposure intensity.

## Key Principle: Backward Compatibility

The load-dependent model is activated by a config flag:
  `disease.load_dependent: true`  (default: false)

When false, the existing timer-based SEIPD system runs unchanged.
All existing tests, configs, and results remain valid.

## Current System (Timer-Based)

```
S → E → I1 → I2 → D  (default at timer expiry)
                  ↘ S  (daily recovery check: p = rho_rec × c_i)
```

- Infection: dose-response on P_k (Michaelis-Menten)
- Progression: Erlang-distributed timers, temperature-scaled
- Mortality: deterministic at I2 timer expiry
- Recovery: daily Bernoulli trial during I2 (and early I1)
- Shedding: fixed rates per compartment (sigma_1, sigma_2)
- No within-host pathogen dynamics

## New System (Load-Dependent)

```
S → Infected (with pathogen_load L_i)
    dL/dt = r_growth(T) × L × (1 - L/L_max) - delta_clear(r_i, c_i) × L
    P(death per day) = Hill(L_i, LD50_i, n_hill)
    shedding_i = sigma_base × L_i / L_ref
    If L_i < L_clear → recovered (S)
```

### Within-Host ODE (daily Euler step per infected agent)

```python
# Pathogen growth (logistic, temperature-dependent)
growth = r_growth * arrhenius(1.0, Ea_growth, T) * L * (1 - L / L_max)

# Immune clearance (resistance + recovery traits)
clearance = delta_clear * (0.3 + 0.7 * r_i) * (0.5 + 0.5 * c_i) * L

# Reinfection from environment
reinfection = alpha_reinfect * P_k / (K_half + P_k) * (1 - r_i)

# Update
dL = growth - clearance + reinfection
L_new = max(0, L + dL)
```

### Load-Dependent Mortality (daily Hill function)

```python
# Per-agent LD50 shifted by tolerance
LD50_i = LD50_base * (1 + t_i * tau_max)  # tolerant hosts need higher load to die

# Daily mortality probability
p_death = p_death_max * L^n_hill / (LD50_i^n_hill + L^n_hill)
```

### Load-Proportional Shedding

```python
# Shedding proportional to load (replaces fixed sigma_1, sigma_2)
shedding_i = sigma_load * L_i / L_ref
# Sum over all infected agents at node → total shedding into P_k ODE
```

### Clearance Threshold

```python
# If load drops below threshold → recovered (S)
if L_i < L_clear:
    disease_state → S
    pathogen_load → 0
```

### Initial Dose at Infection

When S → Infected:
- Draw initial load from: L_0 = L_init_base × P_k / (K_half + P_k)
- Higher environmental pathogen → higher initial dose → harder to clear
- This is where "dose" enters the system

## Compartmental Mapping

The load model REPLACES the E/I1/I2 timer system but we keep disease_state
for behavioral effects and recording:

| Load Level        | Mapped State | Speed | Feeding | Can Spawn |
|-------------------|-------------|-------|---------|-----------|
| L = 0             | S (0)       | 1.0   | 1.0     | Yes       |
| 0 < L < L_symp   | I1 (2)      | 0.5   | 0.5     | No        |
| L >= L_symp       | I2 (3)      | 0.1   | 0.0     | No        |
| Dead              | D (4)       | -     | -       | -         |

- E state is eliminated — infection is immediate with initial load
- disease_timer is repurposed or unused (set to 0)
- The I1/I2 distinction is now emergent from load level, not timer-driven

## New Parameters (LoadDependentSection in config.py)

```python
@dataclass
class LoadDependentSection:
    enabled: bool = False

    # Within-host dynamics
    r_growth: float = 0.5        # Pathogen growth rate (d^-1) at T_ref
    Ea_growth: float = 5000.0    # Growth activation energy (K)
    L_max: float = 1e6           # Within-host carrying capacity
    delta_clear: float = 0.3     # Base immune clearance rate (d^-1)
    alpha_reinfect: float = 0.01 # Environmental reinfection strength

    # Initial infection
    L_init_base: float = 1000.0  # Initial load scale at max dose-response

    # Load thresholds
    L_clear: float = 10.0        # Below this → cleared (recovered)
    L_symp: float = 1e4          # Above this → symptomatic (I2 behavior)

    # Mortality
    LD50_base: float = 5e5       # Load for 50% daily mortality (no tolerance)
    n_hill: float = 3.0          # Hill coefficient (steepness)
    p_death_max: float = 0.15    # Max daily mortality probability at saturating load

    # Shedding
    sigma_load: float = 50.0     # Shedding rate at L_ref load
    L_ref: float = 1e5           # Reference load for shedding normalization
```

## Trait Remapping

Current traits get new mechanistic roles:

| Trait     | Current Role              | New Role (Load Model)                    |
|-----------|--------------------------|------------------------------------------|
| r_i       | Reduces FOI: ×(1-r_i)   | Reduces FOI AND increases clearance rate |
| t_i       | Extends I2 timer         | Increases LD50 (need more load to die)   |
| c_i       | Daily recovery prob      | Increases clearance rate                 |

- r_i=1.0 → clearance = delta_clear × 1.0 → clears any load. TRULY immune.
- t_i=1.0 → LD50 = LD50_base × (1 + tau_max) → very hard to kill
- c_i=1.0 → clearance boosted → faster load reduction

This naturally fixes the "R5 immune but not immune" problem.

## Agent Dtype Changes

```python
# ADD to AGENT_DTYPE:
('pathogen_load', np.float64),  # 8 B — within-host pathogen load
```

Note: pathogen_load already NOT in dtype. The existing pathogen_virulence
field (float32) tracks strain identity, not load. Both coexist.

## Files to Modify

1. **types.py**: Add pathogen_load field to AGENT_DTYPE
2. **config.py**: Add LoadDependentSection dataclass + wire into SimConfig
3. **disease.py**: New function `load_dependent_disease_update()` parallel
   to existing `daily_disease_update()`. Keep old function intact.
4. **model.py**: Branch on config.disease.load_dependent to call correct
   disease function. Update monthly/annual recording.
5. **monthly_recorder.py**: Record pathogen_load stats (mean, max per node)
6. **tests/**: New test file test_load_dependent.py. Existing tests unchanged
   (they use default config where load_dependent=false).

## Files NOT Modified

- genetics.py (traits still computed the same way)
- movement.py (uses disease_state which we still set)
- spawning.py (uses can_spawn which checks disease_state)
- reproduction.py (no disease logic)
- conservation.py (no disease logic)
- viz/ (reads disease_state, still valid)

## Migration Strategy

1. Add pathogen_load to dtype → all agents init with 0.0
2. Add LoadDependentSection to config
3. Implement load_dependent_disease_update() as NEW function
4. Modify model.py to dispatch: if load_dependent → new path, else → old path
5. Old daily_disease_update() untouched
6. Tests: new file for load-dependent, old tests pass unchanged

## Calibration Notes

The load-dependent model has ~12 new parameters. However, many map to
existing concepts:
- r_growth + L_max → determines epidemic speed (analogous to mu_EI1/mu_I1I2)
- delta_clear → determines recovery rate (analogous to rho_rec)
- LD50_base + n_hill → determines case fatality rate (analogous to timer I2→D)
- sigma_load → determines transmission (analogous to sigma_1/sigma_2)

Initial calibration: set parameters so that an average agent (r=0.08, t=0,
c=0) at average temperature (13°C) has similar time-to-death and case
fatality as the timer model. Then sweep from there.

## Sentinel Agent Compatibility

Sentinels in the load model:
- Get infected normally (pathogen_load increases)
- But p_death is set to 0 (immortal, same as current timer restart)
- Still shed proportional to load (with shedding_fraction modifier)
- Load still governed by within-host ODE (can clear, get reinfected)
