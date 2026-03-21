# Pathogen Co-Evolution — Specification

**Date:** 2026-02-17
**Authors:** Willem Weertman & Anton 🔬
**Status:** DRAFT — awaiting review

---

## Summary

Extend the disease module with evolvable pathogen virulence, enabling co-evolutionary
dynamics between *V. pectenicida* and *Pycnopodia helianthoides*. The pathogen evolves
along a mechanistically constrained virulence-transmission tradeoff: more aggressive
strains kill faster and shed more pathogen per day, but their hosts die sooner — limiting
total lifetime transmission. This prevents the model from reaching biologically implausible
parameter space while allowing emergent evolutionary dynamics.

### Motivation

The current model treats *V. pectenicida* as a static entity with fixed disease parameters.
In reality, host-pathogen systems co-evolve. The host is already evolving (resistance loci
under selection), but the pathogen is frozen. This asymmetry biases results toward
inevitably increasing host resistance without the evolutionary counter-pressure that
drives real arms races.

**Key biological observation (Willem):** "When stars begin to waste they put out more
particles, but when not as sick still contagious but less virulence — it's a counterbalance."
This is the virulence-transmission tradeoff in action.

### Key Design Principles

1. **Single evolvable trait** — pathogen virulence phenotype *v* ∈ [0, 1]
2. **Mechanistic tradeoff** — *v* modulates disease parameters through a constrained curve,
   not independently
3. **Small mutation steps** — adaptive dynamics; no evolutionary jumps
4. **Density-dependent selection** — optimal virulence depends on host population density
5. **Backward compatible** — disabled by default; existing behavior unchanged

### Biological Analogies

| System | Observation | Relevance |
|--------|-------------|-----------|
| Myxomatosis–rabbit | CFR dropped 99.8% → 70–90% and stabilized | Classic intermediate virulence optimum |
| DFTD–Tasmanian devil | Coevolution predicted host-pathogen coexistence (Clement 2024) | Eco-evo IBM closest to our architecture |
| HIV | Virulence-transmission tradeoff shapes set-point viral load | R₀ maximized at intermediate virulence |
| *V. aestuarianus*–oyster | Temperature modulates virulence expression | Analogous thermal performance |

---

## §1 — Virulence Phenotype

### 1.1 Definition

Each pathogen strain is characterized by a single scalar virulence phenotype:

**v ∈ [v_min, v_max]** (default [0, 1])

- **v = 0**: minimally virulent — slow progression, low shedding, long infectious period
- **v = 0.5**: ancestral/baseline — matches current default disease parameters
- **v = 1**: maximally virulent — rapid killing, high shedding, short infectious period

The phenotype is NOT a fitness measure — it's a strategy. Fitness (R₀) emerges from the
interaction of *v* with host density, temperature, and resistance.

### 1.2 Per-Agent Tracking

Each infected agent carries the virulence phenotype of the strain that infected it:

| Field | Dtype | Default | Description |
|-------|-------|---------|-------------|
| `pathogen_virulence` | float32 | 0.0 | Virulence of infecting strain (0 when S or R) |

This is a new field in `AGENT_DTYPE` (+4 bytes per agent).

When an agent is in state **S** or **R**, `pathogen_virulence = 0.0` (sentinel).
When transitioning **S → E**, the field is set to the inherited strain virulence (§3).
The field persists through **E → I₁ → I₂ → D** (the same strain throughout infection).
On **recovery (→ R)**, the field resets to 0.0.

---

## §2 — Mechanistic Tradeoff Curve

### 2.1 Rationale

*V. pectenicida* destroys host tissue to replicate. More aggressive tissue destruction:
- Produces more pathogen per unit time (higher shedding)
- Kills the host faster (shorter infectious period)
- Total lifetime pathogen output is mechanistically constrained

This is not an assumed tradeoff — it emerges from the biology of tissue-destructive pathogens.

### 2.2 Parameter Mapping

Virulence *v* modulates disease progression and shedding through exponential scaling
centered on the ancestral phenotype (v = 0.5):

```
Δv = v - v_anchor    (deviation from ancestral, where v_anchor = 0.5)

mu_I2D(v)    = mu_I2D_ref    × exp(α_kill × Δv)     # I₂→D death rate
mu_I1I2(v)   = mu_I1I2_ref   × exp(α_prog × Δv)     # I₁→I₂ progression rate
sigma_2(v)   = sigma_2_ref   × exp(α_shed × Δv)     # Late-stage shedding rate
sigma_1(v)   = sigma_1_ref   × exp(α_shed × Δv × γ_early)  # Early shedding (attenuated)
```

Where:
- **α_kill** — sensitivity of death rate to virulence (curvature of killing)
- **α_prog** — sensitivity of I₁→I₂ progression to virulence
- **α_shed** — sensitivity of shedding to virulence
- **γ_early** — attenuation factor for early-stage shedding (< 1; early disease less affected)

### 2.3 Tradeoff Constraint

The **total lifetime pathogen output** (TLO) of a strain:

```
TLO(v) ≈ σ₁(v)/μ_I1I2(v) + σ₂(v)/μ_I2D(v)
        ≈ σ₁_ref × exp(α_shed×Δv×γ_early − α_prog×Δv) / μ_I1I2_ref
        + σ₂_ref × exp(α_shed×Δv − α_kill×Δv) / μ_I2D_ref
```

For the late-stage dominant term:

```
TLO_late(v) ∝ exp((α_shed − α_kill) × Δv)
```

**If α_kill > α_shed**: TLO decreases with virulence → diminishing returns on aggression.
This creates an **intermediate optimum** because:
- Low *v*: low shedding rate limits instantaneous transmission
- High *v*: fast death limits cumulative transmission
- Intermediate *v*: optimizes shedding rate × duration

### 2.4 R₀ and Density Dependence

The effective reproduction number for a strain at a node:

```
R₀(v) ∝ [a × S × (1−r̄) × S_sal] / [K_half + P] × [σ₁(v)/μ_I1I2(v) + σ₂(v)/μ_I2D(v) + σ_D×τ_D]
```

The dose-response saturation term `P/(K_half + P)` creates **density-dependent selection**:
- **High host density** (pre-epidemic): P can be large → saturation → higher shedding (higher *v*) gives diminishing returns → selection favors moderate *v*
- **Low host density** (post-crash): P is low → linear regime → higher shedding helps more → selection may temporarily favor higher *v*
- **Very low density**: few susceptibles → absolute R₀ low regardless → drift dominates

This means the optimal virulence **shifts during an epidemic**, creating evolutionary dynamics
that can't be captured by static-parameter models.

### 2.5 Default Parameters

| Parameter | Symbol | Default | Range | Confidence | Description |
|-----------|--------|---------|-------|------------|-------------|
| `alpha_kill` | α_kill | 2.0 | [1.0, 4.0] | ★☆☆☆☆ | Death rate scaling exponent |
| `alpha_prog` | α_prog | 1.0 | [0.5, 2.0] | ★☆☆☆☆ | I₁→I₂ progression scaling |
| `alpha_shed` | α_shed | 1.5 | [0.5, 3.0] | ★☆☆☆☆ | Shedding rate scaling exponent |
| `gamma_early` | γ_early | 0.3 | [0.0, 1.0] | ★★☆☆☆ | Early shedding attenuation |
| `v_anchor` | v_anchor | 0.5 | — | — | Ancestral virulence (identity point) |

**Key constraint:** α_kill > α_shed ensures intermediate optimum. If α_kill ≤ α_shed,
runaway virulence or complete attenuation results (may be biologically valid in some
scenarios but not the default).

**Sensitivity:** At defaults (α_kill=2.0, α_shed=1.5), for v=0→1:
- mu_I2D: ×0.37 to ×2.72 (7.3× range)
- sigma_2: ×0.47 to ×2.12 (4.5× range)
- TLO_late: ×1.28 to ×0.78 (modest, ~40% range — intentionally constrained)

---

## §3 — Strain Inheritance and Mutation

### 3.1 Transmission Source Attribution

When an agent transitions **S → E**, we must determine which strain infected it.
Two sources of pathogen:

1. **Local shedders** — I₁ and I₂ individuals at the same node, each with their own *v*
2. **Environmental reservoir** — background VBNC Vibrio with ancestral virulence

Attribution is probabilistic, weighted by contribution to local Vibrio concentration:

```
P_shed = σ₁(T) × n_I1 + σ₂(T) × n_I2   # Total shedding from local hosts
P_env  = environmental_vibrio(T, S, cfg)   # Background environmental input

p_from_shedder = P_shed / (P_shed + P_env)
p_from_env     = P_env  / (P_shed + P_env)
```

With probability `p_from_shedder`: inherit strain from a random shedder (weighted by
individual shedding rate — I₂ individuals shed 10× more, so are 10× more likely to be
the source strain).

With probability `p_from_env`: inherit ancestral strain (v = v_env, default 0.5).

### 3.2 Shedder Selection (Weighted Sampling)

When inheriting from local shedders:

```python
# Weight by individual shedding contribution
weights = []
for each alive agent at node:
    if disease_state == I1:
        w = sigma_1(v_i, T)   # shedding rate of THIS strain
    elif disease_state == I2:
        w = sigma_2(v_i, T)
    else:
        w = 0
weights = normalize(weights)
source = rng.choice(shedders, p=weights)
v_parent = source.pathogen_virulence
```

This naturally selects for high-shedding strains as transmission sources — an evolutionary
force favoring higher shedding, counterbalanced by the tradeoff curve.

### 3.3 Mutation

At each transmission event, the inherited virulence undergoes small mutation:

```
v_offspring = clip(v_parent + N(0, σ_v_mut), v_min, v_max)
```

| Parameter | Symbol | Default | Range | Confidence | Description |
|-----------|--------|---------|-------|------------|-------------|
| `sigma_v_mutation` | σ_v_mut | 0.02 | [0.005, 0.10] | ★☆☆☆☆ | Mutation step size (std dev) |
| `v_init` | v_init | 0.5 | [0.0, 1.0] | ★★☆☆☆ | Initial virulence of all strains |
| `v_min` | v_min | 0.0 | — | — | Minimum virulence bound |
| `v_max` | v_max | 1.0 | — | — | Maximum virulence bound |
| `v_env` | v_env | 0.5 | [0.0, 1.0] | ★★☆☆☆ | Virulence of environmental reservoir |

**Note on σ_v_mut:** This controls evolutionary speed. At 0.02, reaching v=0 from v=0.5
requires ~25 successive low-mutations — many pathogen "generations" (transmission events).
This prevents unrealistic evolutionary jumps while allowing adaptation over epidemic
timescales.

### 3.4 Invasion Scenario

In the invasion scenario (`disease.scenario = "invasion"`):
- Initial pathogen introduced with `v_init` virulence
- No environmental reservoir → all transmission is shedder-to-susceptible
- Mutation proceeds normally

---

## §4 — Integration with Disease Module

### 4.1 Modified Functions

**`daily_disease_update()`** — Main changes:

1. **S → E transition (Step 2):** After determining which susceptibles get infected:
   - For each newly infected agent: determine source strain (§3.1), apply mutation (§3.3)
   - Set `agents['pathogen_virulence'][idx] = v_new`
   - Stage durations sampled using virulence-modified rates: `sample_stage_duration(mu_EI1(v), ...)`

2. **Disease progression (Step 3):** Virulence-dependent progression:
   - E → I₁ timer: sampled using `mu_EI1_ref` (incubation not affected by virulence — pathogen replicating silently)
   - I₁ → I₂ timer: sampled using `mu_I1I2(v)` (§2.2)
   - I₂ → D timer: sampled using `mu_I2D(v)` (§2.2)

3. **Shedding (Step 1):** Vibrio concentration update uses strain-weighted shedding:
   - Can't use single σ₁, σ₂ values anymore — each individual has different effective rates
   - Total shedding = Σ_I1[σ₁(v_i, T)] + Σ_I2[σ₂(v_i, T)]

4. **Recovery:** Recovery probability unchanged — it's a host trait (r_i), not pathogen trait.
   (Pathogen could evolve to evade immunity, but that's a different axis — future extension.)

### 4.2 Modified Shedding Calculation

Currently:
```python
shed = shedding_rate_I1(T, cfg) * n_I1 + shedding_rate_I2(T, cfg) * n_I2 + cfg.sigma_D * n_D_fresh
```

With pathogen evolution:
```python
# Per-individual shedding weighted by strain virulence
v_I1 = agents['pathogen_virulence'][alive & (ds == I1)]
v_I2 = agents['pathogen_virulence'][alive & (ds == I2)]

shed_I1 = np.sum(sigma_1_strain(v_I1, T, cfg))  # vectorized
shed_I2 = np.sum(sigma_2_strain(v_I2, T, cfg))  # vectorized

# Carcass shedding: use virulence of strain that killed the host
# Approximation: use mean v of recent deaths (tracked in CarcassTracker)
shed_D = cfg.sigma_D * n_D_fresh * exp(alpha_shed * (v_mean_carcass - v_anchor))

shed = shed_I1 + shed_I2 + shed_D
```

### 4.3 R₀ Computation

`compute_R0()` becomes strain-dependent. The existing function computes a population-level R₀;
with pathogen evolution, it should accept an optional `v` parameter:

```python
def compute_R0(T, S_0, phi_k, cfg, v=None, ...):
    if v is None:
        v = cfg.pathogen_evolution.v_init  # or mean circulating v
    # Use virulence-adjusted rates
    sigma1 = sigma_1_strain(v, T, cfg)
    sigma2 = sigma_2_strain(v, T, cfg)
    mu_I1I2 = mu_I1I2_strain(v, cfg, T)
    mu_I2D = mu_I2D_strain(v, cfg, T)
    ...
```

### 4.4 New Helper Functions

```python
def sigma_1_strain(v: float|np.ndarray, T: float, cfg) -> float|np.ndarray:
    """Early shedding rate for strain with virulence v at temperature T."""
    dv = v - cfg.pathogen_evolution.v_anchor
    base = arrhenius(cfg.sigma_1_eff, cfg.Ea_sigma, T)
    return base * np.exp(cfg.pathogen_evolution.alpha_shed * dv * cfg.pathogen_evolution.gamma_early)

def sigma_2_strain(v: float|np.ndarray, T: float, cfg) -> float|np.ndarray:
    """Late shedding rate for strain with virulence v at temperature T."""
    dv = v - cfg.pathogen_evolution.v_anchor
    base = arrhenius(cfg.sigma_2_eff, cfg.Ea_sigma, T)
    return base * np.exp(cfg.pathogen_evolution.alpha_shed * dv)

def mu_I2D_strain(v: float|np.ndarray, cfg, T: float) -> float|np.ndarray:
    """I2→D rate for strain with virulence v at temperature T."""
    dv = v - cfg.pathogen_evolution.v_anchor
    base = arrhenius(cfg.mu_I2D_ref, cfg.Ea_I2D, T)
    return base * np.exp(cfg.pathogen_evolution.alpha_kill * dv)

def mu_I1I2_strain(v: float|np.ndarray, cfg, T: float) -> float|np.ndarray:
    """I1→I2 rate for strain with virulence v at temperature T."""
    dv = v - cfg.pathogen_evolution.v_anchor
    base = arrhenius(cfg.mu_I1I2_ref, cfg.Ea_I1I2, T)
    return base * np.exp(cfg.pathogen_evolution.alpha_prog * dv)
```

All functions accept scalar or array *v* for vectorized per-individual computation.

---

## §5 — Configuration

### 5.1 New Config Section

```python
@dataclass
class PathogenEvolutionSection:
    """Pathogen co-evolution parameters.

    When enabled, V. pectenicida evolves virulence along a mechanistic
    tradeoff curve. Host-pathogen coevolution emerges from the interaction
    of host resistance evolution (genetics module) with pathogen virulence
    evolution (this module).
    """
    enabled: bool = False          # Off by default (backward compatible)

    # Initial conditions
    v_init: float = 0.5            # Initial virulence of all strains
    v_env: float = 0.5             # Virulence of environmental reservoir (ubiquitous)

    # Bounds
    v_min: float = 0.0             # Minimum virulence
    v_max: float = 1.0             # Maximum virulence

    # Mutation
    sigma_v_mutation: float = 0.02 # Mutation step size (normal std dev)

    # Tradeoff curve
    v_anchor: float = 0.5          # Ancestral virulence (identity point for curve)
    alpha_kill: float = 2.0        # Death rate scaling exponent
    alpha_prog: float = 1.0        # I1→I2 progression scaling exponent
    alpha_shed: float = 1.5        # Shedding rate scaling exponent
    gamma_early: float = 0.3       # Early shedding attenuation (0=no change, 1=full)

    # Output
    track_strain_history: bool = False  # Record per-transmission strain lineages
```

### 5.2 Config Integration

Add `pathogen_evolution` as a sub-section of `DiseaseSection` (conceptually part of disease),
OR as a top-level section in `SimulationConfig` (parallel to `genetics`).

**Recommended:** Top-level section, because pathogen evolution is the mirror of host
evolution — they are peers, not parent-child.

```python
@dataclass
class SimulationConfig:
    ...
    genetics: GeneticsSection
    pathogen_evolution: PathogenEvolutionSection  # NEW
    ...
```

---

## §6 — Output and Recording

### 6.1 Per-Node Annual Metrics

Add to the annual recording output:

| Metric | Description |
|--------|-------------|
| `mean_virulence` | Mean *v* of all currently infected (E+I₁+I₂) individuals |
| `std_virulence` | Std dev of *v* across infected individuals |
| `mean_virulence_new` | Mean *v* of newly infected individuals this year |
| `virulence_of_deaths` | Mean *v* of individuals who died of disease this year |
| `n_strains_unique` | Number of distinct strain clusters (v binned to 0.01) |

### 6.2 CoupledSimResult Extension

```python
@dataclass
class CoupledSimResult:
    ...
    # Pathogen evolution (annual timeseries)
    yearly_mean_virulence: Optional[np.ndarray] = None       # (n_years,)
    yearly_virulence_new_infections: Optional[np.ndarray] = None
    yearly_virulence_of_deaths: Optional[np.ndarray] = None
```

### 6.3 Spatial Sim Extension

```python
# Per-node per-year virulence tracking
yearly_mean_virulence: np.ndarray   # (n_nodes, n_years)
```

### 6.4 Optional Lineage Tracking

When `track_strain_history = True`, record each transmission event:

```python
@dataclass
class TransmissionEvent:
    day: int
    node_id: int
    source_agent: int          # -1 if environmental
    dest_agent: int
    v_parent: float
    v_offspring: float
    source_disease_state: int  # I1 or I2
```

Stored as a list, written to file at end. Enables phylogenetic visualization of strain
evolution. **WARNING:** Memory-intensive at scale — use only for targeted analysis runs.

---

## §7 — Validation Targets

### 7.1 Analytical

| Test | Expected | Tolerance | Description |
|------|----------|-----------|-------------|
| **Static baseline** | v=0.5 gives identical results to pathogen_evolution.enabled=False | Exact | Ancestral virulence reproduces current model |
| **TLO constraint** | TLO varies < 2× across v=[0,1] at default params | Check | Tradeoff limits total pathogen output range |
| **Intermediate optimum** | R₀ maximized at v ∈ (0.2, 0.8) for default α_kill, α_shed | Calculate | Confirms non-trivial evolutionary dynamics |

### 7.2 Evolutionary Dynamics

| Test | Expected | Notes |
|------|----------|-------|
| **No mutation (σ=0)** | Mean v stays at v_init forever | Sanity check |
| **Drift only** | Mean v diffuses randomly with σ ∝ √(gen) | Neutral evolution |
| **Selection with tradeoff** | v converges toward intermediate optimum | Core prediction |
| **Density dependence** | Optimal v shifts between pre/post epidemic | Key novel prediction |
| **Myxomatosis pattern** | Starting from v=1, virulence decreases toward intermediate | Classic analog |

### 7.3 Comparison with Literature

| Comparison | Expected | Reference |
|------------|----------|-----------|
| DFTD co-evolution | Host-pathogen coexistence (vs extinction in host-only) | Clement et al. 2024 |
| Myxomatosis | Virulence decrease from high initial | Fenner & Fantini 1999 |
| Intermediate virulence | ESS at v* where dR₀/dv = 0 | Anderson & May 1982 |

---

## §8 — Implementation Phases

### Phase 1: Data Structures (1 hour)
- Add `pathogen_virulence` to AGENT_DTYPE
- Add `PathogenEvolutionSection` to config.py
- Add to SimulationConfig
- Add config validation (alpha_kill > 0, sigma_v_mutation ≥ 0, v_min < v_max)
- Tests: field exists, config loads, backward compatibility

### Phase 2: Tradeoff Functions (1 hour)
- Implement `sigma_1_strain()`, `sigma_2_strain()`, `mu_I2D_strain()`, `mu_I1I2_strain()`
- Both scalar and vectorized (np.ndarray) inputs
- At v = v_anchor, output equals base rate (identity test)
- Numerical tests: monotonicity, range, TLO constraint
- Update `compute_R0()` to accept optional `v` parameter

### Phase 3: Strain Inheritance (2 hours)
- Source attribution: shedder vs environmental (§3.1)
- Shedder-weighted sampling (§3.2)
- Mutation at transmission (§3.3)
- Modify `daily_disease_update()` Step 2 (S→E)
- Tests: inheritance works, mutation bounded, environmental fallback

### Phase 4: Virulence-Dependent Dynamics (2 hours)
- Modify shedding calculation in `update_vibrio_concentration()` or caller
- Per-individual stage duration sampling using virulence-adjusted rates
- Modify `daily_disease_update()` Steps 1 and 3
- Tests: high-v kills faster, low-v kills slower, TLO approximately constant

### Phase 5: Output Recording (1 hour)
- Annual virulence metrics in CoupledSimResult
- Per-node virulence tracking in spatial sim
- Optional lineage tracking
- Tests: metrics recorded correctly

### Phase 6: Integration & Validation (2 hours)
- Wire into `run_coupled_simulation()` and `run_spatial_simulation()`
- Backward compatibility test: disabled = identical to current
- Multi-seed evolutionary dynamics validation
- Performance benchmark (expected: <10% overhead — only extra work is per-individual
  exp() calls and weighted sampling at transmission)

**Estimated total: ~9 hours implementation + testing**

---

## §9 — Performance Considerations

### 9.1 Overhead Analysis

The main new costs per daily timestep:

1. **Strain-weighted shedding** (Step 1): Replace `n_I1 × σ₁` with `Σ σ₁(v_i)`.
   Cost: one `exp()` per infected agent per day. At 1000 infected agents: ~1000 exp() calls.
   NumPy vectorized → negligible.

2. **Source attribution at transmission** (Step 2): Weighted random choice among shedders.
   Cost: O(n_shedders) per new infection. Typically <100 shedders → negligible.

3. **Virulence-adjusted timers** (Step 3): One extra `exp()` per stage transition.
   Cost: O(n_transitions/day) → negligible.

**Expected total overhead: < 5%** of current disease module runtime.
Disease module is currently <1% of total runtime (spawning dominates), so net impact
on simulation speed: **< 0.05%**.

### 9.2 Memory

One float32 per agent: +4 bytes × max_agents.
At 10,000 agents: +40 KB. Negligible.

Optional lineage tracking: ~32 bytes per transmission event.
At 10,000 infections: ~320 KB. Manageable.

---

## §10 — Future Extensions

### 10.1 Multi-Trait Evolution (Not This Spec)

The current single-trait (virulence) model could be extended to:
- **Immune evasion**: pathogen evolving to reduce effectiveness of host resistance genes
- **Temperature adaptation**: pathogen T_opt evolving (relevant under climate change)
- **Transmissibility vs virulence**: separating dose-response from tissue destruction

Each would add one evolvable trait dimension. Start with virulence alone;
add dimensions only if single-trait predictions are insufficient.

### 10.2 Explicit Phylogenetics

With lineage tracking enabled, strain phylogenies can be reconstructed post-hoc.
This enables:
- Molecular clock calibration
- Identification of selective sweeps in pathogen population
- Comparison with real V. pectenicida genomic diversity (when available)

### 10.3 Superinfection

Currently not modeled: a host infected with strain A cannot be re-infected with strain B.
This is a simplification. Superinfection (within-host competition between strains) would
add complexity but could affect evolutionary dynamics if strains compete for the same
host tissues.

**Recommendation:** Defer unless single-strain results are qualitatively wrong.

---

## §11 — Key References

1. **Anderson & May 1982** — "Coevolution of hosts and parasites" — original
   virulence-transmission tradeoff theory
2. **Alizon & van Baalen 2005** — "Emergence of pathogen virulence" — tradeoff curve
   shapes and evolutionary stable strategies
3. **Boots & Haraguchi 1999** — "Adaptive dynamics and the evolution of virulence" —
   small-step adaptive dynamics for host-pathogen systems
4. **Clement et al. 2024** (Evolution) — DFTD eco-evo IBM; coevolution predicted
   coexistence, host-only predicted extinction. CRITICAL REFERENCE (awaiting PDF).
5. **Clement et al. 2025** (PLOS Pathogens) — Expanded DFTD coevolution model
6. **Fenner & Fantini 1999** — "Biological control of vertebrate pests" — myxomatosis
   empirical system, virulence attenuation from CFR 99.8% → 70-90%
7. **Prentice 2025** (Nature E&E) — Koch's postulates for V. pectenicida FHCF-3;
   confirms single-pathogen etiology for Pycnopodia
8. **Royal Society B 2022** — "Population size affects host-pathogen coevolution" —
   small populations constrain host adaptation rate

---

## Appendix A — Worked Example

### Scenario: Ancestral Strain (v=0.5)

At v=0.5, Δv=0: all rates equal their reference values.
- mu_I2D(0.5) = mu_I2D_ref × exp(0) = 0.173 d⁻¹
- sigma_2(0.5) = sigma_2_ref × exp(0) = 50.0 bact/mL/d/host
- Mean I₂ duration = 1/0.173 ≈ 5.8 days
- TLO_late ≈ 50/0.173 ≈ 289 bact·d/mL/host

### Scenario: Hypervirulent Strain (v=0.9)

Δv = +0.4:
- mu_I2D(0.9) = 0.173 × exp(2.0 × 0.4) = 0.173 × 2.23 = 0.386 d⁻¹
- sigma_2(0.9) = 50.0 × exp(1.5 × 0.4) = 50.0 × 1.82 = 91.1 bact/mL/d/host
- Mean I₂ duration = 1/0.386 ≈ 2.6 days
- TLO_late ≈ 91.1/0.386 ≈ 236 bact·d/mL/host (18% LESS than ancestral)

### Scenario: Attenuated Strain (v=0.1)

Δv = −0.4:
- mu_I2D(0.1) = 0.173 × exp(2.0 × −0.4) = 0.173 × 0.449 = 0.078 d⁻¹
- sigma_2(0.1) = 50.0 × exp(1.5 × −0.4) = 50.0 × 0.549 = 27.4 bact/mL/d/host
- Mean I₂ duration = 1/0.078 ≈ 12.9 days
- TLO_late ≈ 27.4/0.078 ≈ 352 bact·d/mL/host (22% MORE than ancestral)

**→ Attenuated strain produces MORE total pathogen despite lower shedding rate.**
**→ But at lower instantaneous concentration, dose-response saturation may limit its advantage.**
**→ The equilibrium depends on host density — this is the emergent dynamics we want.**

---

## Appendix B — Sensitivity Analysis Implications

When pathogen evolution is enabled, 5 new parameters enter the model:
- α_kill, α_prog, α_shed, γ_early, σ_v_mut

Plus initial conditions: v_init, v_env.

**Recommendation:** Run structural comparison first (co-evolution ON vs OFF) before
sweeping co-evolution parameters. If co-evolution qualitatively changes outcomes
(e.g., prevents extinction, as in Clement 2024 DFTD), then sweep the tradeoff curve
parameters in a focused SA.
