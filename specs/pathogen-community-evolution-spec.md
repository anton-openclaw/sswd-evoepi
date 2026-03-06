# Pathogen Community Multi-Trait Evolution — Unified Spec

**Date:** 2026-03-05  
**Status:** DRAFT  
**Builds on:** pathogen-community-adaptation-spec.md (Chunk 1 ✅), pathogen-evolution-spec.md (Feb 17)

---

## Motivation

We currently evolve one pathogen trait per site: **T_vbnc_local** (temperature tolerance).
But real pathogen communities evolve along multiple axes simultaneously, constrained by
the **virulence-transmission trade-off** (Anderson & May 1982):

- A pathogen that kills too fast doesn't spread (short infectious period)
- A pathogen that's too mild doesn't shed enough to transmit
- Cold adaptation may come at a metabolic cost to virulence

Willem's key insight: the site's bacterial community is the evolutionary unit. Individual
bacteria are short-lived, but the community composition persists. The community evolves
along ALL trait axes simultaneously, with trade-offs constraining the evolutionary space.

### Why This Matters for the North-South Gradient

**Southern sites** (warm, dense pre-crash):
- Pathogen already temperature-adapted (T_vbnc ≈ SST)
- Can afford high virulence — many hosts, fast transmission
- Acute outbreaks, massive crashes, little recovery

**Northern sites** (cold, sparse post-crash):
- Pathogen needs cold tolerance (T_vbnc must drop)
- High virulence is *maladaptive* — kills host before finding next one
- Natural selection favors LOWER virulence: chronic, sublethal infection
- Populations suppressed but not annihilated

This gradient **self-organizes** from the trade-off without us hard-coding it.

---

## Design: Site-Level Community Traits

### Per-Site Pathogen State

Each site's pathogen community is characterized by THREE evolving traits stored
in `NodeDiseaseState`:

```python
T_vbnc_local: float = 12.0     # Temperature tolerance (EXISTING, Chunk 1 ✅)
v_local: float = 0.5           # Community mean virulence (NEW)
a_local: float = 1.0           # Community mean infectivity modifier (NEW)
```

**v_local** (virulence, 0–1): Scales disease progression rates and shedding.
- v=0: minimally virulent — slow progression, low shedding, long infectious period
- v=0.5: ancestral/baseline — current default disease parameters
- v=1: maximally virulent — rapid killing, high shedding, short infectious period

**a_local** (infectivity modifier, 0–2): Scales the exposure coefficient a_exposure.
- a=1.0: baseline infectivity (default a_exposure)
- a<1.0: reduced infectivity (cost of cold adaptation or low virulence)
- a>1.0: enhanced infectivity (specialized to local hosts)

### Virulence Modulates Disease Parameters

Following the existing pathogen-evolution-spec.md §2, virulence scales rates:

```
Δv = v_local - v_anchor    (deviation from ancestral 0.5)

mu_I2D_eff    = mu_I2D_ref    × exp(α_kill × Δv)     # Kill rate
mu_I1I2_eff   = mu_I1I2_ref   × exp(α_prog × Δv)     # Progression rate
sigma_2_eff   = sigma_2_ref   × exp(α_shed × Δv)     # Late shedding
sigma_1_eff   = sigma_1_ref   × exp(α_shed × Δv × γ_early)  # Early shedding
```

**Key constraint**: α_kill > α_shed ensures intermediate virulence optimum.
At defaults (α_kill=2.0, α_shed=1.5):
- High v: kills fast, sheds a lot per day, but host dies quickly → less total output
- Low v: sheds less per day, but host lives longer → MORE total lifetime output
- Optimal v balances instantaneous shedding rate × infectious duration

### Infectivity Modulates Exposure

```
a_effective = a_exposure × a_local
```

Simple multiplicative modifier on the base exposure coefficient.

### The Trade-Off Constraint

Adapting to cold temperatures has a metabolic cost. The community can't simultaneously
maximize cold tolerance, virulence, AND infectivity. We model this as:

**Option A — Soft trade-off (recommended):**
Each trait evolves independently toward its local optimum, but cold adaptation
naturally reduces virulence because:
1. Cold-adapted bacteria have slower metabolism → can't sustain high replication
2. Slower replication → less tissue destruction → lower shedding
3. This is captured by making v_local drift toward a **temperature-dependent
   equilibrium**: colder sites favor lower v

```python
# Virulence equilibrium depends on temperature
v_eq(T) = v_max_warm × sigmoid((T - T_v_mid) / T_v_width)
# Warm (T >> T_v_mid): v_eq → v_max_warm (high virulence viable)
# Cold (T << T_v_mid): v_eq → 0 (low virulence favored)
```

**Option B — Hard budget (simpler):**
Total adaptation "distance" from ancestral state is bounded:

```
|ΔT_vbnc/T_scale|² + |Δv/v_scale|² + |Δa/a_scale|² ≤ budget²
```

Pushing one trait further constrains the others.

**Recommendation:** Start with Option A (soft/emergent). The temperature-virulence
coupling is more biologically honest and avoids arbitrary budget parameters.

---

## Daily Evolution Update

### Temperature Adaptation (EXISTING — Chunk 1 ✅)

```python
# Already implemented: pool-driven, Michaelis-Menten saturation
if P_env_pool > 0 and T < T_vbnc_local:
    pool_factor = min(P_env_pool / P_adapt_half, 1.0)
    delta_T = adapt_rate * pool_factor * (T_vbnc_local - T)
    T_vbnc_local -= delta_T
```

### Virulence Adaptation (NEW)

The community's mean virulence drifts toward the local optimum, which depends on:
1. **Host density** — more hosts = higher optimal v (can afford fast transmission)
2. **Temperature** — colder = lower optimal v (metabolic constraint)

```python
def adapt_virulence(node_state, T_celsius, n_hosts, K, cfg):
    """Evolve site virulence toward local optimum."""
    if node_state.P_env_pool <= 0:
        return
    
    # Density-dependent optimal virulence
    # High density → high v optimal; low density → low v optimal
    density_ratio = n_hosts / max(K, 1)
    
    # Temperature constraint: cold → lower max virulence
    temp_factor = sigmoid((T_celsius - cfg.T_v_mid) / cfg.T_v_width)
    
    # Local optimum: product of density and temperature effects
    v_optimal = cfg.v_max_warm * density_ratio * temp_factor
    v_optimal = clip(v_optimal, cfg.v_min, cfg.v_max)
    
    # Drift toward optimum
    pool_factor = min(node_state.P_env_pool / max(cfg.P_adapt_half, 1e-9), 1.0)
    delta_v = cfg.v_adapt_rate * pool_factor * (v_optimal - node_state.v_local)
    node_state.v_local += delta_v
```

The **density dependence** is key: 
- Pre-crash (high density): v_optimal high → virulent strains favored
- Post-crash (low density): v_optimal low → attenuated strains favored
- This creates the temporal dynamics: virulence DECREASES after the crash

### Infectivity Adaptation (NEW, optional)

Infectivity drifts toward local optimum based on host resistance:
- If hosts have high resistance (r̄ high): pathogen evolves higher infectivity to compensate
- If hosts have low resistance: baseline infectivity sufficient

```python
def adapt_infectivity(node_state, mean_resistance, cfg):
    """Evolve site infectivity in response to host resistance."""
    if node_state.P_env_pool <= 0:
        return
    
    # Higher host resistance → selection for higher infectivity
    a_optimal = 1.0 + cfg.a_resistance_response * mean_resistance
    a_optimal = clip(a_optimal, cfg.a_min, cfg.a_max)
    
    pool_factor = min(node_state.P_env_pool / max(cfg.P_adapt_half, 1e-9), 1.0)
    delta_a = cfg.a_adapt_rate * pool_factor * (a_optimal - node_state.a_local)
    node_state.a_local += delta_a
```

**Note:** Infectivity adaptation could be deferred (Phase 2) — it's the arms race
component. Temperature + virulence are the priority for the gradient problem.

---

## Dispersal Carries All Traits

Extending the Chunk 2 dispersal mixing to carry virulence alongside temperature:

```python
# Vectorized: mix all traits proportionally to pathogen mass transfer
T_vec = np.array([nds.T_vbnc_local for nds in node_disease_states])
V_vec = np.array([nds.v_local for nds in node_disease_states])

PT_vec = P * T_vec   # concentration-weighted T_vbnc
PV_vec = P * V_vec   # concentration-weighted virulence

incoming_PT = D_T_sparse @ PT_vec
incoming_PV = D_T_sparse @ PV_vec
incoming_P = dispersal_in  # already computed

P_before = P.copy()
total = P_before + incoming_P
mask = total > 0

# Mix both traits
T_vec[mask] = (P_before[mask] * T_vec[mask] + incoming_PT[mask]) / total[mask]
V_vec[mask] = (P_before[mask] * V_vec[mask] + incoming_PV[mask]) / total[mask]

# Write back
for i in range(N):
    if total[i] > 0:
        node_disease_states[i].T_vbnc_local = T_vec[i]
        node_disease_states[i].v_local = V_vec[i]
```

This means cold-adapted, low-virulence strains from northern sites will slowly
spread south via dispersal, and vice versa. But local selection will maintain the
gradient — warm southern sites will select back toward high virulence.

---

## Wavefront Inheritance

When a new site is activated by the wavefront, it inherits ALL community traits
from source nodes (extending existing `_inherit_T_vbnc`):

```python
def _inherit_pathogen_traits(target_idx, node_disease_states, disease_reached,
                              D_T_sparse, P, cfg):
    """Inherit weighted-average pathogen community traits from source nodes."""
    # Same weighted-average logic as _inherit_T_vbnc, but for all traits
    weighted_T = weighted_V = weighted_A = total_weight = 0.0
    
    for each source with disease_reached and dispersal weight:
        w = dispersal_weight * P[source]
        weighted_T += w * node_disease_states[source].T_vbnc_local
        weighted_V += w * node_disease_states[source].v_local
        weighted_A += w * node_disease_states[source].a_local
        total_weight += w
    
    if total_weight > 0:
        node_disease_states[target_idx].T_vbnc_local = weighted_T / total_weight
        node_disease_states[target_idx].v_local = weighted_V / total_weight
        node_disease_states[target_idx].a_local = weighted_A / total_weight
```

---

## How Virulence Enters the Disease Engine

Currently, disease parameters are global (same at every node). With community
virulence, they become **per-node**:

```python
# In daily_disease_update(), BEFORE computing rates:
v = node_state.v_local
dv = v - cfg_pe.v_anchor  # deviation from ancestral

# Adjust rates for this node's pathogen community
mu_I2D_local = mu_I2D_base * exp(cfg_pe.alpha_kill * dv)
mu_I1I2_local = mu_I1I2_base * exp(cfg_pe.alpha_prog * dv)
sigma_2_local = sigma_2_base * exp(cfg_pe.alpha_shed * dv)
sigma_1_local = sigma_1_base * exp(cfg_pe.alpha_shed * dv * cfg_pe.gamma_early)
a_local = a_exposure * node_state.a_local
```

This is clean: instead of per-agent strain tracking, we just modify the effective
rates at each node based on the site's community virulence. Much simpler than the
per-agent approach in the original pathogen-evolution-spec.md, and consistent with
the community-level paradigm.

**Key simplification vs original spec:** No per-agent `pathogen_virulence` field,
no strain inheritance at transmission, no weighted shedder sampling. The community
IS the strain. This is appropriate because bacteria turn over fast — the community
mean is what matters, not individual lineages.

---

## Config Changes

```python
# In DiseaseSection (or PathogenEvolutionSection):

# Virulence evolution
virulence_evolution: bool = False       # Enable community virulence evolution
v_adapt_rate: float = 0.001            # Rate of virulence drift toward optimum
v_max_warm: float = 0.7               # Max optimal virulence at warm sites
T_v_mid: float = 12.0                 # Temperature midpoint for virulence sigmoid
T_v_width: float = 3.0                # Width of virulence-temperature sigmoid

# Infectivity evolution (Phase 2)  
infectivity_evolution: bool = False     # Enable infectivity co-evolution
a_adapt_rate: float = 0.0005           # Rate of infectivity drift
a_resistance_response: float = 1.0     # How much infectivity responds to resistance
a_min: float = 0.5                     # Minimum infectivity modifier
a_max: float = 2.0                     # Maximum infectivity modifier

# Existing PathogenEvolutionSection fields used for trade-off curve:
# alpha_kill, alpha_prog, alpha_shed, gamma_early, v_anchor
# These are already in config.py — no changes needed
```

---

## Implementation Chunks

### Chunk 2: Dispersal trait mixing (from community-adaptation-spec)
- Carry T_vbnc AND v_local through dispersal
- Vectorized weighted-average mixing
- Extend `_inherit_T_vbnc` → `_inherit_pathogen_traits`

### Chunk 3: Community virulence evolution
- Add `v_local` to `NodeDiseaseState`
- Add `adapt_virulence()` function
- Wire v_local into disease rate calculations (per-node rate adjustment)
- Config fields for virulence evolution
- Tests

### Chunk 4: Infectivity evolution (optional, Phase 2)
- Add `a_local` to `NodeDiseaseState`
- Add `adapt_infectivity()` function  
- Wire a_local into exposure calculation
- Tests

---

## Expected Emergent Dynamics

### During Initial Outbreak (Year 1-3)
- Disease arrives from south with ancestral v=0.5
- High host density → v_optimal high → virulence stays high or increases
- Acute crashes everywhere disease reaches

### Post-Crash Recovery Phase (Year 3-13)
- **South (warm, crashed)**: 
  - Few hosts → v_optimal drops
  - But warm temps allow high virulence metabolically
  - Net: v stays moderate-high → continued suppression
  - P_env floor keeps pathogen present

- **North (cold, recovering)**:
  - Few hosts → v_optimal drops  
  - Cold temps further constrain virulence
  - Net: v drops substantially → chronic sublethal disease
  - Populations can recover WHILE infected
  - T_vbnc also drops → longer disease season BUT lower virulence

- **JDF/SS-S (intermediate)**:
  - Intermediate density, intermediate temperature
  - v moderate → moderate suppression
  - Semi-enclosed → pathogen retention → sustained pressure

### The Gradient Emerges
- **South**: high v × warm × floor → acute, populations stay crashed
- **North**: low v × cold-adapted × sparse → chronic, populations recover
- **Intermediate**: moderate everything → moderate suppression

This is the gradient we've been trying to build! And it emerges from first principles
of evolutionary epidemiology, not from parameter tuning.

---

## Recording

Per-region annual output (extend existing calibration recording):
- `final_mean_v_local`: mean community virulence per region
- `yearly_v_local`: yearly snapshots like yearly_T_vbnc
- `final_mean_a_local`: mean community infectivity per region (if enabled)
