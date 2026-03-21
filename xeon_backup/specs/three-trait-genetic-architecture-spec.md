# Three-Trait Genetic Architecture Specification

**Version:** 1.0  
**Date:** 2026-02-19  
**Authors:** Anton 🔬 & Willem Weertman  
**Status:** APPROVED — all open questions resolved (Feb 19)

## 1. Summary

Replace the current 52-locus architecture (51 additive resistance + 1 overdominant EF1A) with a **51-locus, three-trait architecture** partitioning loci equally into:

| Trait | Loci | Indices | Symbol | Mechanistic role |
|-------|------|---------|--------|-----------------|
| **Resistance** | 17 | 0–16 | r_i | Reduce probability of infection |
| **Tolerance** | 17 | 17–33 | t_i | Reduce mortality while infected |
| **Recovery** | 17 | 34–50 | c_i | Clear pathogen and return to S/R |

**EF1A (locus 51) is removed.** It was based on *Pisaster ochraceus* data (Wares 2016) with no evidence of overdominance in *Pycnopodia helianthoides*.

Total loci: **51** (down from 52). Genotype array: `(max_n, 51, 2)` int8.

## 2. Biological Rationale

These three traits represent distinct immune strategies with different evolutionary dynamics:

- **Resistance** = immune exclusion. Receptor polymorphisms, barrier defenses, innate recognition. Prevents the pathogen from establishing. Reduces force of infection on the individual AND reduces pathogen pressure on the population (fewer shedding hosts).

- **Tolerance** = damage limitation. Tissue repair, anti-inflammatory regulation, metabolic compensation. The host gets infected and sheds pathogen normally but survives longer. **Epidemiological consequence: tolerant hosts are silent spreaders.** They maintain pathogen pressure on the population while saving themselves.

- **Recovery** = pathogen clearance. Coelomocyte phagocytosis, immune effector mobilization. The host actively clears the infection and transitions to recovered (immune) state. Intermediate epidemiological effect: removes shedding host from the pool, but only after a period of infection.

**Selection pressure differences:**
- Resistance is selected *before* infection (avoided entirely)
- Tolerance is selected *during* infection (survive the damage)
- Recovery is selected *during late infection* (clear before death)

All three traits are under different selection intensities depending on epidemic severity, pathogen virulence, and population density.

**Empirical basis:** Schiebelhut et al. (2018, 2024) identified ~51 loci under selection in SSWD survivors. We have no GWAS data distinguishing resistance from tolerance from recovery loci. The equal 17/17/17 partition is a simplifying assumption; the partition ratio itself can be explored as an SA parameter in future work.

## 3. Changes to `types.py`

### 3.1 Constants

```python
# REMOVE
N_LOCI = 52
N_ADDITIVE = 51
IDX_EF1A = 51

# REPLACE WITH
N_LOCI = 51            # Total diploid loci (no EF1A). Fixed constant.

# Default partition (configurable via GeneticsSection):
N_RESISTANCE_DEFAULT = 17
N_TOLERANCE_DEFAULT = 17
N_RECOVERY_DEFAULT = 17

# Index slices computed dynamically from config:
def trait_slices(n_r: int, n_t: int, n_c: int):
    """Compute locus index slices from partition sizes.
    
    Validates n_r + n_t + n_c == N_LOCI.
    """
    assert n_r + n_t + n_c == N_LOCI, f"Partition must sum to {N_LOCI}"
    return (
        slice(0, n_r),                    # resistance
        slice(n_r, n_r + n_t),            # tolerance
        slice(n_r + n_t, n_r + n_t + n_c) # recovery
    )
```

**Decision: Variable partition.** The 17/17/17 split is the default but configurable via `GeneticsSection`. Total is always constrained to N_LOCI=51. For SA, parameterize as two free variables (n_resistance, n_tolerance); n_recovery = 51 − n_resistance − n_tolerance.

### 3.2 AGENT_DTYPE changes

```python
# REMOVE
('fecundity_mod',  np.float32),   # Always 1.0 (CE-1), dead weight

# ADD
('tolerance',      np.float32),   # Tolerance score t_i ∈ [0, 1]
('recovery_ability', np.float32), # Recovery/clearance score c_i ∈ [0, 1]
```

The `resistance` field (float32) is retained unchanged.

**Net dtype change:** remove `fecundity_mod` (4B), add `tolerance` (4B) + `recovery_ability` (4B) = **+4 bytes** per agent (55→59 bytes, offset by removing EF1A column from genotype array).

### 3.3 Genotype array

```python
def allocate_genotypes(max_n: int) -> np.ndarray:
    return np.zeros((max_n, N_LOCI, 2), dtype=np.int8)
    # Shape: (max_n, 51, 2) — was (max_n, 52, 2)
```

## 4. Changes to `config.py`

### 4.1 `GeneticsSection` updates

```python
@dataclass
class GeneticsSection:
    # REMOVE
    # n_additive: int = 51
    # n_loci: int = 52
    # s_het: float = 0.19
    # q_ef1a_init: float = 0.24

    # REPLACE WITH
    n_resistance: int = 17        # Loci 0 .. n_resistance-1
    n_tolerance: int = 17         # Loci n_resistance .. n_resistance+n_tolerance-1
    n_recovery: int = 17          # Loci n_resistance+n_tolerance .. sum-1
    # Constraint: n_resistance + n_tolerance + n_recovery = N_LOCI (validated at init)

    # Per-trait initialization
    target_mean_r: float = 0.15       # Population-mean resistance at t=0
    target_mean_t: float = 0.10       # Population-mean tolerance at t=0
    target_mean_c: float = 0.08       # Population-mean recovery at t=0

    # Tolerance mechanics
    tau_max: float = 0.85             # Max I₂→D mortality reduction at t_i=1

    # Shared genetics parameters (unchanged)
    mu_per_locus: float = 1.0e-8
    n_bank: int = 100
    effect_size_seed: int = 12345
    q_init_mode: str = "beta"
    q_init_beta_a: float = 2.0
    q_init_beta_b: float = 8.0
```

**Design note on initial trait means:**
- `target_mean_r = 0.15` — pre-outbreak resistance (established in genetics overhaul)
- `target_mean_t = 0.10` — low initial tolerance; most of the reduction in μ_I2D comes from high-t_i tail
- `target_mean_c = 0.08` — low initial clearance; recovery is rare pre-outbreak (consistent with field data)

### 4.2 `DiseaseSection` updates

```python
# REMOVE
recovery_exponent: int = 2    # No longer used; recovery is c_i-based

# No other DiseaseSection changes needed.
# rho_rec stays as the base clearance rate.
```

## 5. Changes to `genetics.py`

### 5.1 Remove all EF1A code

Delete:
- `S_HET`, `W_OD`, `W_ADD` constants
- EF1A heterozygote bonus in `compute_resistance_single()` and `compute_resistance_batch()`
- `eliminate_ef1a_lethals()` function
- EF1A initialization in `initialize_genotypes()` and `initialize_genotypes_beta()`
- EF1A statistics in `GeneticDiagnostics`

### 5.2 Effect size initialization

Three independent effect size arrays, one per trait:

```python
def initialize_trait_effect_sizes(
    rng: np.random.Generator,
    n_loci: int,
    total_weight: float = 1.0,
) -> np.ndarray:
    """Draw and normalize effect sizes for one trait.

    Exp(λ) distribution, sorted descending, normalized to sum = total_weight.
    """
    raw = rng.exponential(scale=1.0, size=n_loci)
    normalized = raw / raw.sum() * total_weight
    normalized.sort()
    return normalized[::-1].copy()
```

Called three times with independent RNG draws:
```python
effects_resistance = initialize_trait_effect_sizes(rng, N_RESISTANCE, total_weight=1.0)
effects_tolerance  = initialize_trait_effect_sizes(rng, N_TOLERANCE,  total_weight=1.0)
effects_recovery   = initialize_trait_effect_sizes(rng, N_RECOVERY,   total_weight=1.0)
```

Each trait score ranges [0, 1] independently. The `total_weight=1.0` ensures that a fully homozygous-derived individual at all 17 loci has trait score = 1.0.

### 5.3 Trait score computation

Generalized batch function:

```python
def compute_trait_batch(
    genotypes: np.ndarray,
    effects: np.ndarray,
    alive_mask: np.ndarray,
    locus_slice: slice,
) -> np.ndarray:
    """Vectorized trait score computation for one trait.

    score_i = Σ e_l × (g_l0 + g_l1) / 2, summed over trait's loci.

    Args:
        genotypes: (max_agents, N_LOCI, 2) int8.
        effects: (n_trait_loci,) float64 effect sizes.
        alive_mask: (max_agents,) bool.
        locus_slice: slice object for this trait's loci (e.g., RESISTANCE_SLICE).

    Returns:
        (max_agents,) float32 — trait score for alive agents, 0.0 for dead.
    """
    max_agents = genotypes.shape[0]
    scores = np.zeros(max_agents, dtype=np.float32)
    alive_idx = np.where(alive_mask)[0]
    if len(alive_idx) == 0:
        return scores
    allele_sums = genotypes[alive_idx, locus_slice, :].sum(axis=2)
    allele_means = allele_sums.astype(np.float64) * 0.5
    scores[alive_idx] = (allele_means @ effects).astype(np.float32)
    return scores
```

Called as (slices computed from config at simulation init):
```python
res_slice, tol_slice, rec_slice = trait_slices(cfg.n_resistance, cfg.n_tolerance, cfg.n_recovery)
r_vals = compute_trait_batch(genotypes, effects_resistance, alive, res_slice)
t_vals = compute_trait_batch(genotypes, effects_tolerance,  alive, tol_slice)
c_vals = compute_trait_batch(genotypes, effects_recovery,   alive, rec_slice)
```

### 5.4 Update all trait scores

```python
def update_all_trait_scores(agents, genotypes, effects_r, effects_t, effects_c):
    """Recompute and write r_i, t_i, c_i for all alive agents."""
    alive = agents['alive']
    agents['resistance']       = compute_trait_batch(genotypes, effects_r, alive, RESISTANCE_SLICE)
    agents['tolerance']        = compute_trait_batch(genotypes, effects_t, alive, TOLERANCE_SLICE)
    agents['recovery_ability'] = compute_trait_batch(genotypes, effects_c, alive, RECOVERY_SLICE)
```

### 5.5 Genotype initialization

```python
def initialize_genotypes_three_trait(
    n_agents: int,
    effects_r: np.ndarray,
    effects_t: np.ndarray,
    effects_c: np.ndarray,
    rng: np.random.Generator,
    target_mean_r: float = 0.15,
    target_mean_t: float = 0.10,
    target_mean_c: float = 0.08,
    beta_a: float = 2.0,
    beta_b: float = 8.0,
) -> np.ndarray:
    """Initialize 51-locus genotypes with per-trait target means.

    Each trait block gets Beta-distributed per-locus allele frequencies,
    scaled so population-mean trait score ≈ target.
    """
    geno = np.zeros((n_agents, N_LOCI, 2), dtype=np.int8)

    for locus_slice, effects, target in [
        (RESISTANCE_SLICE, effects_r, target_mean_r),
        (TOLERANCE_SLICE,  effects_t, target_mean_t),
        (RECOVERY_SLICE,   effects_c, target_mean_c),
    ]:
        n_trait = effects.shape[0]
        raw_q = rng.beta(beta_a, beta_b, size=n_trait)
        current = np.dot(effects, raw_q)
        if current > 0:
            scale = target / current
            q_vals = np.clip(raw_q * scale, 0.001, 0.5)
        else:
            q_vals = np.full(n_trait, 0.01)

        start = locus_slice.start
        for i in range(n_trait):
            l_idx = start + i
            geno[:, l_idx, 0] = (rng.random(n_agents) < q_vals[i]).astype(np.int8)
            geno[:, l_idx, 1] = (rng.random(n_agents) < q_vals[i]).astype(np.int8)

    return geno
```

### 5.6 GeneticDiagnostics

```python
@dataclass
class GeneticDiagnostics:
    allele_freq: np.ndarray  # (N_LOCI=51,)

    # Per-trait distribution stats
    mean_resistance: float = 0.0
    var_resistance: float = 0.0
    va_resistance: float = 0.0

    mean_tolerance: float = 0.0
    var_tolerance: float = 0.0
    va_tolerance: float = 0.0

    mean_recovery: float = 0.0
    var_recovery: float = 0.0
    va_recovery: float = 0.0

    # Diversity (across all 51 loci)
    heterozygosity_obs: float = 0.0
    heterozygosity_exp: float = 0.0

    # Allele frequency change
    delta_q_top3: np.ndarray  # Top 3 absolute Δq across all loci
    mean_abs_delta_q: float = 0.0

    n_alive: int = 0

    # REMOVED: ef1a_allele_freq, ef1a_het_freq
```

### 5.7 Additive variance

Compute V_A per trait:

```python
def compute_additive_variance(allele_freq, effects, locus_slice):
    q = allele_freq[locus_slice]
    return float(2.0 * np.sum(effects**2 * q * (1.0 - q)))
```

### 5.8 Genotype bank

`GenotypeBankState` adds `bank_tolerance` and `bank_recovery` arrays (float32, shape `(N_BANK,)`). `update_summary()` computes mean/var for all three traits. `compress_to_genotype_bank()` and `expand_genotype_bank()` pass all three effect arrays.

### 5.9 Mutation

`apply_mutations()` is unchanged — it operates on the full genotype array and doesn't distinguish trait blocks. Mutations at any locus are equally likely. This is correct: all loci share the same per-locus mutation rate.

### 5.10 F_ST

`compute_fst()` is unchanged — operates on all N_LOCI. Could optionally compute per-trait F_ST (resistance loci only, etc.) as a diagnostic.

## 6. Changes to `disease.py`

### 6.1 Force of infection — UNCHANGED

```python
λ_i = a × P/(K_half + P) × (1 − r_i) × S_sal × f_size(L_i)
```

Only resistance `r_i` appears. Tolerance and recovery do not affect infection probability. This is the defining feature of resistance vs tolerance.

### 6.2 I₂→D mortality — tolerance modulates

Current (timer-based): when `disease_timer` reaches 0 in I₂, the agent dies.

**New mechanism:** Daily survival probability during I₂, modulated by tolerance:

```python
def daily_I2_death_probability(t_i: float, mu_I2D: float, tau_max: float) -> float:
    """Daily probability of disease death in I₂ state.

    p_death = μ_I2D × (1 − t_i × τ_max)

    At t_i=0: p_death = μ_I2D (no tolerance benefit)
    At t_i=1: p_death = μ_I2D × (1 − τ_max) (maximum reduction)

    Args:
        t_i: Individual tolerance score [0, 1].
        mu_I2D: Base I₂→D rate at current temperature (d⁻¹).
        tau_max: Maximum mortality reduction (default 0.85).
    """
    return mu_I2D * (1.0 - t_i * tau_max)
```

**Implementation in `daily_disease_update()`:**

In the I₂ timer-expired block, instead of automatic death:

```python
elif state == DiseaseState.I2:
    # Timer expired — but tolerance may prevent death
    t_i = agents['tolerance'][idx]
    effective_rate = mu_I2D * (1.0 - t_i * tau_max)
    p_death = 1.0 - np.exp(-effective_rate)  # convert rate to probability
    if rng.random() < p_death:
        # Death
        ds[idx] = DiseaseState.D
        agents['alive'][idx] = False
        agents['cause_of_death'][idx] = 1
        new_deaths += 1
    else:
        # Survived this round — reset timer for another I₂ period
        dt_rem[idx] = sample_stage_duration(mu_I2D, K_SHAPE_I2, rng)
```

**Alternative (simpler, preserves timer semantics):** Scale the I₂ timer duration by tolerance. Tolerant individuals get longer timers → more chances for recovery before death:

```python
# When entering I₂:
effective_rate = mu_I2D * (1.0 - t_i * tau_max)
dt_rem[idx] = sample_stage_duration(effective_rate, K_SHAPE_I2, rng)
# At timer=0: death (deterministic, as before)
```

**Decision: Timer-scaling.** Tolerant individuals get longer I₂ timers → more recovery chances before deterministic death. Simpler, preserves existing control flow, caps I₂ extension predictably (~7× at max tolerance). Avoids geometric superspreader tails from daily-probability approach.

### 6.3 Recovery — recovery ability (`c_i`) modulates

**Decision: Linear in c_i.** Gives more even selection pressure across the c_i distribution than the old quadratic (r_i²) approach.

```python
def recovery_probability_I2(c_i: float, rho_rec: float = 0.05) -> float:
    """Daily probability of recovery from I₂.

    p_rec = ρ_rec × c_i

    Linear in clearance ability. At c_i=0: no recovery possible.
    At c_i=1: p_rec = ρ_rec (base rate, ~5%/day).
    """
    return rho_rec * c_i


def recovery_probability_I1(c_i: float, rho_rec: float = 0.05) -> float:
    """Daily probability of early recovery from I₁.

    Only meaningful for high-clearance individuals (c_i > 0.5).
    Linear above threshold.
    p_early = ρ_rec × max(0, 2(c_i − 0.5))
    """
    if c_i <= 0.5:
        return 0.0
    return rho_rec * 2.0 * (c_i - 0.5)
```

**Key change:** Recovery now uses `c_i` (recovery ability) instead of `r_i` (resistance). Resistance no longer influences recovery at all — each trait does exactly one thing.

In `daily_disease_update()`, all `recovery_probability_*(r_i, ...)` calls become `recovery_probability_*(c_i, ...)`:

```python
# I₂ recovery check
c_i = agents['recovery_ability'][idx]
p_rec = recovery_probability_I2(c_i, cfg.rho_rec)

# I₁ early recovery check
c_i = agents['recovery_ability'][idx]
p_early = recovery_probability_I1(c_i, cfg.rho_rec)
```

### 6.4 Disease progression E→I₁→I₂ — UNCHANGED

Incubation and early progression are NOT affected by any host trait. This is biologically correct: the pathogen replicates regardless of host genotype. Only the *consequences* differ (death rate, clearance probability).

### 6.5 Shedding — UNCHANGED

Shedding rates σ₁, σ₂ depend on disease state and temperature only. Tolerant hosts shed exactly as much as non-tolerant hosts. This is the key epidemiological feature of tolerance.

### 6.6 R₀ computation

Update `compute_R0()` to accept `mean_tolerance` and `mean_recovery`:

```python
def compute_R0(..., mean_tolerance=0.0, mean_recovery=0.0, tau_max=0.85):
    # Effective I₂→D rate accounting for population-mean tolerance
    mu_I2D_eff = mu_I2D * (1.0 - mean_tolerance * tau_max)
    # Effective I₂ duration accounting for recovery
    # Mean time in I₂ ≈ 1 / (μ_I2D_eff + ρ_rec × mean_recovery)
    effective_I2_exit = mu_I2D_eff + cfg.rho_rec * mean_recovery
    ...
```

### 6.7 Pathogen evolution interactions

When pathogen evolution is enabled, strain-specific rates interact with host tolerance:

```python
# Per-individual effective I₂→D rate with both PE and tolerance:
rate_I2D = mu_I2D_strain(v_i, T, cfg, pe_cfg) * (1.0 - t_i * tau_max)
```

This creates a three-way coevolutionary dynamic:
- Pathogen evolves higher virulence → higher μ_I2D
- Host evolves tolerance → lower effective μ_I2D
- Host evolves recovery → higher clearance before death

## 7. Changes to other modules

### 7.1 `reproduction.py`

- Remove any reference to `fecundity_mod` (CE-1: already always 1.0)
- Mendelian inheritance operates on full 51-locus genotype — no changes needed if using `N_LOCI`
- After offspring genotypes are produced, call `update_all_trait_scores()` for settlers

### 7.2 Recorder / spatial simulation

- Track per-node per-year: `mean_resistance`, `mean_tolerance`, `mean_recovery`
- Track `delta_r`, `delta_t`, `delta_c` (allele frequency shifts in each trait block)
- `CoupledSimResult` and `SpatialSimResult` add tolerance/recovery time series

### 7.3 Visualization

- Existing resistance evolution plots generalize to three-panel (r/t/c) versions
- Allele frequency heatmaps: color-code by trait block
- New plot: trait correlation scatter (r_i vs t_i vs c_i) — do they evolve independently?

### 7.4 Sensitivity analysis

New/modified parameters for SA:
- `n_resistance` (default 17, range 5–30)
- `n_tolerance` (default 17, range 5–30, constrained: sum = 51)
- `n_recovery` (default 17, range 5–30, constrained: sum = 51)
- `target_mean_t` (range 0.02–0.30)
- `target_mean_c` (range 0.02–0.25)
- `tau_max` (range 0.3–0.95)

Removed from SA:
- `s_het` (EF1A, removed)
- `q_ef1a_init` (EF1A, removed)
- `n_additive` (replaced by n_resistance; total fixed at 51)

**Note:** The locus partition (n_resistance/n_tolerance/n_recovery) is constrained to sum to 51. For SA, parameterize as two free variables (e.g., n_resistance and n_tolerance; n_recovery = 51 − n_resistance − n_tolerance).

## 8. Implementation phases

### Phase 1: Types + Config (est. 30 min)
- Update constants in `types.py`
- Remove `fecundity_mod`, add `tolerance` + `recovery_ability` to AGENT_DTYPE
- Update `GeneticsSection` in `config.py`
- Update `SimulationConfig.from_dict()` parsing

### Phase 2: Genetics core (est. 1.5 hr)
- Remove all EF1A code
- Implement `initialize_trait_effect_sizes()`
- Implement `compute_trait_batch()` (generalized)
- Implement `initialize_genotypes_three_trait()`
- Update `update_all_trait_scores()` (was `update_resistance_scores()`)
- Update `GeneticDiagnostics` and `compute_genetic_diagnostics()`
- Update genotype bank (compress/expand/update_summary)

### Phase 3: Disease wiring (est. 1 hr)
- Modify `recovery_probability_I2()` and `recovery_probability_I1()` to use `c_i`
- Add tolerance modulation of I₂→D mortality (timer-scaling approach)
- Update `daily_disease_update()` to read `tolerance` and `recovery_ability` fields
- Update `compute_R0()` for population-mean tolerance/recovery

### Phase 4: Integration (est. 1 hr)
- Update reproduction module offspring trait computation
- Update spatial simulation to track three traits
- Update recorder/result dataclasses
- Wire through `SimulationConfig` → `SpatialSim` → `CoupledSim`

### Phase 5: Tests (est. 1.5 hr)
- **Invariant tests** (CRITICAL):
  - Tolerant individuals get infected at the same rate as non-tolerant (resistance held equal)
  - Tolerant individuals shed at the same rate as non-tolerant
  - High-c_i individuals recover more often than low-c_i (resistance and tolerance held equal)
  - High-t_i individuals survive I₂ longer than low-t_i
  - Resistance still reduces infection probability
- Update all existing genetics tests (remove EF1A expectations)
- Update all existing disease tests (recovery uses c_i not r_i)
- Trait independence: verify three trait blocks don't cross-contaminate
- Genotype initialization: verify per-trait target means

### Phase 6: Validation run (est. 30 min)
- Run 5-node 20-year simulation
- Verify all three traits evolve under selection
- Compare population dynamics to pre-refactor baseline
- Check for emergent dynamics (does tolerance create silent spreaders?)

## 9. Backward compatibility

**Breaking change.** The genotype array shape changes from `(N, 52, 2)` → `(N, 51, 2)` and the interpretation of locus indices changes entirely. Old checkpoints and SA results are incompatible.

- SA Rounds 1–3: completed, results archived. R3 Sobol is still running but uses `n_additive` as a parameter — those results remain interpretable as "total resistance loci" even under the new architecture.
- Old saved simulations: cannot be loaded. Not a concern (no production runs yet).
- Tests: all will need updating. No partial migration path.

## 10. Design decisions (resolved Feb 19)

1. **Timer-scaling** for tolerance. Tolerant individuals get longer I₂ timers → more recovery chances. Simpler, predictable I₂ duration extension (~7× at max tolerance). Avoids geometric superspreader tails.

2. **Linear recovery** in c_i. `p_rec = ρ_rec × c_i`. Even selection pressure across the c_i distribution. Old quadratic (r_i²) created threshold effects that concentrated selection on the high tail.

3. **Variable partition** in config. Default 17/17/17 but configurable. Total constrained to 51. SA: two free variables (n_resistance, n_tolerance); n_recovery = 51 − sum.

4. **Initial means:** `target_mean_r = 0.15`, `target_mean_t = 0.10`, `target_mean_c = 0.08`. Low tolerance/recovery = weak initial effect → strong selection signal when epidemic hits.

5. **Tolerance does NOT slow E→I₁→I₂ progression.** Each trait does exactly one thing. Resistance blocks infection, tolerance extends survival in I₂, recovery clears pathogen. Clean separation, no blurred boundaries.
