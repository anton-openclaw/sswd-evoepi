"""Coupled disease ↔ population ↔ genetics simulation.

Single-node simulation (Phases 5+7) and multi-node spatial simulation
(Phase 9):
  - Daily loop: environment update → disease step → pathogen dispersal
  - Annual loop: natural mortality → growth → stage transitions →
                 reproduction → larval dispersal → recruitment
  - Disease kills individuals (sets alive=False), reducing population
  - Reduced population → reduced reproduction → Allee effects
  - Post-epidemic recovery dynamics
  - Genetics ↔ Disease coupling:
      r_i feeds into force of infection λ_i (higher r_i → lower infection prob)
      Disease kills low-r_i individuals → survivors enriched for resistance alleles
      Post-epidemic reproduction passes resistance alleles to offspring via SRS
      Allele frequencies tracked pre/post-epidemic for calibration
  - Spatial coupling (Phase 9):
      Daily: pathogen exchange between nodes via D matrix
      Annual: larvae dispersed between nodes via C matrix

References:
  - integration-architecture-spec.md §2.1 (master loop pseudocode)
  - spatial-connectivity-spec.md §8 (node-level simulation loop)
  - population-dynamics-spec.md §2 (lifecycle)
  - disease-module-spec.md (SEIPD+R)
  - genetics-evolution-spec.md §6 (selection mechanisms), §9 (calibration)
  - CODE_ERRATA CE-1: no cost of resistance
  - CODE_ERRATA CE-4: corrected Beverton-Holt formula
  - CODE_ERRATA CE-5: demographic Allee for high-fecundity species

Build target: Phase 5 (single-node) + Phase 7 (genetics) + Phase 9 (spatial).
"""

from __future__ import annotations

import bisect
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

from sswd_evoepi.config import (
    DiseaseSection,
    PopulationSection,
    SimulationConfig,
    SpawningSection,
    default_config,
)
from sswd_evoepi.disease import (
    CarcassTracker,
    NodeDiseaseState,
    daily_disease_update,
    environmental_vibrio,
    vibrio_decay_rate,
    arrhenius,
    compute_R0,
)
from sswd_evoepi.genetics import (
    compute_allele_frequencies,
    compute_additive_variance,
    compute_genetic_diagnostics,
    initialize_trait_effect_sizes,
    update_all_trait_scores,
    update_resistance_scores,
    RESISTANCE_SLICE,
    TOLERANCE_SLICE,
    RECOVERY_SLICE,
)
from sswd_evoepi.reproduction import (
    beverton_holt_recruitment,
    fecundity,
    fertilization_success,
    get_spawning_day,
    is_spawning_season,
    larval_survival,
    mendelian_inherit_batch,
    pelagic_larval_duration,
    settle_recruits,
    settlement_cue_modifier,
    srs_reproductive_lottery,
    _compute_resistance,
)
from sswd_evoepi.perf import PerfMonitor
from sswd_evoepi.spawning import (
    spawning_step,
    reset_spawning_season,
    in_spawning_season,
)
from sswd_evoepi.types import (
    AGENT_DTYPE,
    ANNUAL_SURVIVAL,
    LarvalCohort,
    N_LOCI,
    N_RESISTANCE_DEFAULT,
    STAGE_SIZE_THRESHOLDS,
    DiseaseState,
    Stage,
    allocate_agents,
    allocate_genotypes,
    trait_slices,
)


# ═══════════════════════════════════════════════════════════════════════
# CONSTANTS
# ═══════════════════════════════════════════════════════════════════════

DAYS_PER_YEAR = 365


# ═══════════════════════════════════════════════════════════════════════
# EFFECT SIZES (deterministic from seed)
# ═══════════════════════════════════════════════════════════════════════

def make_effect_sizes(seed: int = 12345, n_loci: int = N_RESISTANCE_DEFAULT, total_weight: float = 1.0) -> np.ndarray:
    """Create canonical effect-size vector (exponential, normalized).

    Sorted descending so first loci have largest effects.
    Total weight sums to total_weight (default 1.0).

    CE-3: Exponential distribution per Lotterhos & Whitlock 2016.
    """
    rng = np.random.default_rng(seed)
    raw = rng.exponential(1.0, size=n_loci)
    normalized = raw / raw.sum() * total_weight
    normalized.sort()
    return normalized[::-1].copy()


# ═══════════════════════════════════════════════════════════════════════
# VON BERTALANFFY GROWTH
# ═══════════════════════════════════════════════════════════════════════

def von_bertalanffy(age: float, L_inf: float = 1000.0,
                    k: float = 0.08, t0: float = -0.5) -> float:
    """Expected diameter (mm) from age (years)."""
    return L_inf * (1.0 - np.exp(-k * (age - t0)))


def grow_individual(current_size: float, L_inf: float = 1000.0,
                    k: float = 0.08, dt_years: float = 1.0,
                    rng: Optional[np.random.Generator] = None) -> float:
    """Advance size by dt_years using VB differential form with noise.

    dL/dt = k × (L_inf − L)
    Individual variation: multiplicative log-normal noise (CV ≈ 0.15).
    """
    expected = k * (L_inf - current_size) * dt_years
    if expected <= 0:
        return current_size
    if rng is not None:
        noise = np.exp(rng.normal(0.0, 0.15))
        expected *= noise
    return max(current_size, current_size + expected)


# ═══════════════════════════════════════════════════════════════════════
# STAGE TRANSITIONS
# ═══════════════════════════════════════════════════════════════════════

def assign_stage(size_mm: float, current_stage: int) -> int:
    """Assign life stage based on size thresholds. One-directional only."""
    if current_stage == Stage.EGG_LARVA:
        return Stage.EGG_LARVA  # handled externally
    if size_mm >= 400.0:
        return Stage.ADULT
    if size_mm >= 150.0:
        return max(current_stage, Stage.SUBADULT)
    if size_mm >= 10.0:
        return max(current_stage, Stage.JUVENILE)
    return max(current_stage, Stage.SETTLER)


# ═══════════════════════════════════════════════════════════════════════
# NATURAL MORTALITY
# ═══════════════════════════════════════════════════════════════════════

def natural_mortality_prob(stage: int, age: float,
                           annual_survival: np.ndarray = ANNUAL_SURVIVAL,
                           senescence_age: float = 50.0,
                           senescence_mortality: float = 0.10) -> float:
    """Annual probability of natural death for one individual."""
    base_mort = 1.0 - annual_survival[stage]
    if age > senescence_age:
        extra = senescence_mortality * (age - senescence_age) / 20.0
        return min(1.0, base_mort + extra)
    return base_mort


# ═══════════════════════════════════════════════════════════════════════
# POPULATION INITIALIZATION
# ═══════════════════════════════════════════════════════════════════════

def _make_per_locus_q(
    n_loci: int,
    effect_sizes: np.ndarray,
    target_mean: float,
    mode: str,
    rng: np.random.Generator,
    beta_a: float = 2.0,
    beta_b: float = 8.0,
    **kwargs,
) -> np.ndarray:
    """Compute per-locus allele frequencies for a trait block.

    Args:
        n_loci: Number of loci for this trait.
        effect_sizes: Effect size vector.
        target_mean: Target population-mean trait score.
        mode: "uniform" or "beta".
        rng: Random generator.
        beta_a, beta_b: Beta distribution shape parameters.
        **kwargs: Absorbs legacy arguments (q_ef1a, w_od, n_additive, etc.).

    Returns:
        (n_loci,) float64 per-locus allele frequencies.
    """
    if mode == "beta":
        raw_q = rng.beta(beta_a, beta_b, size=n_loci)
        current = np.dot(effect_sizes, raw_q)
        if current > 0:
            scale = target_mean / current
            q_vals = np.clip(raw_q * scale, 0.001, 0.999)
        else:
            q_vals = np.full(n_loci, 0.01)
    else:
        q_uniform = target_mean / effect_sizes.sum() if effect_sizes.sum() > 0 else 0.01
        q_uniform = np.clip(q_uniform, 0.001, 0.999)
        q_vals = np.full(n_loci, q_uniform)

    return q_vals


def initialize_population(
    n_individuals: int,
    max_agents: int,
    habitat_area: float,
    effect_sizes: np.ndarray,
    pop_cfg: PopulationSection,
    rng: np.random.Generator,
    q_init: float = 0.05,
    genetics_cfg: Optional['GeneticsSection'] = None,
    effects_t: Optional[np.ndarray] = None,
    effects_c: Optional[np.ndarray] = None,
    **kwargs,
) -> Tuple[np.ndarray, np.ndarray]:
    """Initialize a population at demographic equilibrium.

    Creates individuals with age/size drawn from an approximate
    stable age distribution, random genotypes with given allele
    frequencies, and computed resistance/tolerance/recovery scores.

    Three-trait initialization: if ``genetics_cfg`` and all three
    effect arrays are provided, uses ``initialize_genotypes_three_trait``
    to set up all 51 loci with per-trait target means. Otherwise
    falls back to resistance-only initialization.

    Args:
        n_individuals: Target number of live individuals.
        max_agents: Array capacity (must be >= n_individuals).
        habitat_area: Habitat area (m²) for spatial placement.
        effect_sizes: Resistance effect sizes.
        pop_cfg: Population configuration.
        rng: Random generator.
        q_init: Initial resistant allele frequency (legacy fallback).
        genetics_cfg: GeneticsSection with initialization config.
        effects_t: Tolerance effect sizes (optional).
        effects_c: Recovery effect sizes (optional).
        **kwargs: Absorbs legacy arguments (q_ef1a, w_od, etc.).

    Returns:
        (agents, genotypes) tuple.
    """
    from sswd_evoepi.genetics import initialize_genotypes_three_trait

    agents = allocate_agents(max_agents)
    genotypes = allocate_genotypes(max_agents)

    n_res = len(effect_sizes)

    hab_side = np.sqrt(max(habitat_area, 1.0))

    # --- Vectorized demographic initialization ---
    u_ages = rng.random(n_individuals)
    ages = np.empty(n_individuals)
    ages[u_ages < 0.05] = rng.uniform(0.0, 1.0, size=(u_ages < 0.05).sum())
    mask_juv = (u_ages >= 0.05) & (u_ages < 0.20)
    ages[mask_juv] = rng.uniform(1.0, 3.0, size=mask_juv.sum())
    mask_sub = (u_ages >= 0.20) & (u_ages < 0.35)
    ages[mask_sub] = rng.uniform(3.0, 7.0, size=mask_sub.sum())
    mask_adult = u_ages >= 0.35
    ages[mask_adult] = rng.uniform(7.0, 30.0, size=mask_adult.sum())

    sizes = von_bertalanffy(ages, pop_cfg.L_inf, pop_cfg.k_growth,
                            pop_cfg.t0_growth)
    sizes *= rng.lognormal(0.0, 0.10, size=n_individuals)
    sizes = np.maximum(sizes, 0.5)

    stages = np.full(n_individuals, Stage.ADULT, dtype=np.int8)
    for stage, threshold_size in sorted(
        STAGE_SIZE_THRESHOLDS.items(), key=lambda x: -x[1]
    ):
        stages[sizes < threshold_size] = stage

    # Genotype initialization: three-trait if all effect arrays provided
    if genetics_cfg is not None and effects_t is not None and effects_c is not None:
        # Full three-trait initialization
        init_geno = initialize_genotypes_three_trait(
            n_agents=n_individuals,
            effects_r=effect_sizes,
            effects_t=effects_t,
            effects_c=effects_c,
            rng=rng,
            target_mean_r=genetics_cfg.target_mean_r,
            target_mean_t=genetics_cfg.target_mean_t,
            target_mean_c=genetics_cfg.target_mean_c,
            beta_a=genetics_cfg.q_init_beta_a,
            beta_b=genetics_cfg.q_init_beta_b,
            n_resistance=genetics_cfg.n_resistance,
            n_tolerance=genetics_cfg.n_tolerance,
            n_recovery=genetics_cfg.n_recovery,
        )
        genotypes[:n_individuals] = init_geno
    else:
        # Legacy: resistance-only initialization
        if genetics_cfg is not None:
            q_per_locus = _make_per_locus_q(
                n_loci=n_res,
                effect_sizes=effect_sizes,
                target_mean=genetics_cfg.target_mean_r,
                mode=genetics_cfg.q_init_mode,
                rng=rng,
                beta_a=genetics_cfg.q_init_beta_a,
                beta_b=genetics_cfg.q_init_beta_b,
            )
        else:
            q_per_locus = np.full(n_res, q_init)

        for l_idx in range(n_res):
            genotypes[:n_individuals, l_idx, 0] = (rng.random(n_individuals) < q_per_locus[l_idx]).astype(np.int8)
            genotypes[:n_individuals, l_idx, 1] = (rng.random(n_individuals) < q_per_locus[l_idx]).astype(np.int8)

    # Fill agent fields
    sl = slice(0, n_individuals)
    agents['alive'][sl] = True
    agents['x'][sl] = rng.uniform(0, hab_side, size=n_individuals)
    agents['y'][sl] = rng.uniform(0, hab_side, size=n_individuals)
    agents['heading'][sl] = rng.uniform(0, 2 * np.pi, size=n_individuals)
    agents['speed'][sl] = 0.1
    agents['size'][sl] = sizes
    agents['age'][sl] = ages
    agents['stage'][sl] = stages
    agents['sex'][sl] = rng.integers(0, 2, size=n_individuals)
    agents['disease_state'][sl] = DiseaseState.S
    agents['disease_timer'][sl] = 0
    agents['node_id'][sl] = 0
    agents['origin'][sl] = 0  # WILD
    agents['settlement_day'][sl] = 0

    # Compute trait scores for initial population
    if effects_t is not None and effects_c is not None:
        # Full three-trait score computation
        res_s, tol_s, rec_s = trait_slices(
            genetics_cfg.n_resistance if genetics_cfg else N_RESISTANCE_DEFAULT,
            genetics_cfg.n_tolerance if genetics_cfg else N_TOLERANCE_DEFAULT,
            genetics_cfg.n_recovery if genetics_cfg else N_RECOVERY_DEFAULT,
        )
        update_all_trait_scores(
            agents, genotypes, effect_sizes, effects_t, effects_c,
            res_slice=res_s, tol_slice=tol_s, rec_slice=rec_s,
        )
    else:
        # Resistance-only (legacy)
        for i in range(n_individuals):
            agents[i]['resistance'] = _compute_resistance(
                genotypes[i], effect_sizes
            )

    return agents, genotypes


# ═══════════════════════════════════════════════════════════════════════
# ANNUAL DEMOGRAPHIC UPDATE
# ═══════════════════════════════════════════════════════════════════════

def annual_natural_mortality(
    agents: np.ndarray,
    pop_cfg: PopulationSection,
    rng: np.random.Generator,
) -> Tuple[int, int]:
    """Apply annual natural mortality to all alive individuals.

    Returns (n_killed, n_senescence) — total natural deaths and
    the subset that died from senescence (age > senescence_age).
    """
    alive_mask = agents['alive'].astype(bool)
    alive_idx = np.where(alive_mask)[0]
    n_killed = 0
    n_senescence = 0

    annual_surv = np.array(pop_cfg.annual_survival, dtype=np.float64)

    for idx in alive_idx:
        stage = int(agents['stage'][idx])
        age = float(agents['age'][idx])
        p_death = natural_mortality_prob(
            stage, age, annual_surv,
            pop_cfg.senescence_age, pop_cfg.senescence_mortality,
        )
        if rng.random() < p_death:
            agents['alive'][idx] = False
            # Distinguish senescence from stage-specific natural mortality
            if age > pop_cfg.senescence_age:
                agents['cause_of_death'][idx] = 3  # DeathCause.SENESCENCE
                n_senescence += 1
            else:
                agents['cause_of_death'][idx] = 2  # DeathCause.NATURAL
            n_killed += 1

    return n_killed, n_senescence


def annual_growth_and_aging(
    agents: np.ndarray,
    pop_cfg: PopulationSection,
    rng: np.random.Generator,
) -> None:
    """Advance age by 1 year, grow via VB, update stages."""
    alive_mask = agents['alive'].astype(bool)
    alive_idx = np.where(alive_mask)[0]

    for idx in alive_idx:
        agents['age'][idx] += 1.0
        new_size = grow_individual(
            float(agents['size'][idx]),
            pop_cfg.L_inf,
            pop_cfg.k_growth,
            dt_years=1.0,
            rng=rng,
        )
        agents['size'][idx] = new_size
        agents['stage'][idx] = assign_stage(
            new_size, int(agents['stage'][idx])
        )


# ═══════════════════════════════════════════════════════════════════════
# CONTINUOUS (DAILY) DEMOGRAPHIC FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════
# These replace the annual lump-sum mortality/growth in the simulation
# loop. The old annual_* functions above are retained for backward
# compatibility and testing.
#
# Reference: specs/continuous-mortality-spec.md

_INV_365 = 1.0 / 365.0  # pre-computed for daily timestep


def daily_natural_mortality(
    agents: np.ndarray,
    pop_cfg: PopulationSection,
    rng: np.random.Generator,
    alive_idx: np.ndarray = None,
) -> Tuple[int, int]:
    """Apply daily natural mortality to all alive individuals (vectorized).

    Converts annual stage-specific survival rates to daily death
    probabilities via ``p_death_daily = 1 - S_annual^(1/365)``.
    Senescence overlay for agents older than ``senescence_age``.

    Args:
        agents: Structured agent array (modified in place).
        pop_cfg: Population configuration with annual_survival,
            senescence_age, senescence_mortality.
        rng: NumPy random generator.
        alive_idx: Pre-computed alive indices (optional; computed if None).

    Returns:
        (n_killed, n_senescence) — total natural deaths this day and
        the subset that died from senescence.
    """
    alive = alive_idx if alive_idx is not None else np.where(agents['alive'])[0]
    if len(alive) == 0:
        return 0, 0

    stages = agents['stage'][alive]
    ages = agents['age'][alive]

    # Vectorized annual survival lookup by stage
    annual_surv = np.array(pop_cfg.annual_survival, dtype=np.float64)
    s_clipped = np.clip(stages, 0, len(annual_surv) - 1)
    base_annual_mort = 1.0 - annual_surv[s_clipped]

    # Senescence overlay: extra mortality for age > senescence_age
    total_annual_mort = base_annual_mort.copy()
    senes_mask = ages > pop_cfg.senescence_age
    if np.any(senes_mask):
        extra = (pop_cfg.senescence_mortality
                 * (ages[senes_mask] - pop_cfg.senescence_age) / 20.0)
        total_annual_mort[senes_mask] = np.minimum(
            1.0, base_annual_mort[senes_mask] + extra
        )

    # Convert annual mortality → daily mortality
    daily_mort = 1.0 - (1.0 - total_annual_mort) ** _INV_365

    # Single vectorized random draw
    rolls = rng.random(len(alive))
    dies = rolls < daily_mort

    # Apply deaths
    dead_idx = alive[dies]
    agents['alive'][dead_idx] = False

    # Stamp cause of death: SENESCENCE (3) or NATURAL (2)
    dead_ages = agents['age'][dead_idx]
    senes_dead = dead_ages > pop_cfg.senescence_age
    # Use field[slot] = value pattern (CE-critical for structured arrays)
    agents['cause_of_death'][dead_idx[senes_dead]] = 3   # SENESCENCE
    agents['cause_of_death'][dead_idx[~senes_dead]] = 2  # NATURAL

    n_killed = int(np.sum(dies))
    n_senescence = int(np.sum(senes_dead))
    return n_killed, n_senescence


def daily_growth_and_aging(
    agents: np.ndarray,
    pop_cfg: PopulationSection,
    rng: np.random.Generator,
    alive_idx: np.ndarray = None,
) -> None:
    """Advance age by 1/365 year, grow via VB, update stages (vectorized).

    Von Bertalanffy differential form with dt=1/365:
        ``new_size = L_inf - (L_inf - old_size) * exp(-k * 1/365)``

    Growth noise: multiplicative log-normal with daily-scaled sigma
    (σ_daily = 2.0 / √365 ≈ 0.105 mm).

    Stage transitions checked via size thresholds (one-directional).
    Size is clamped to never decrease (no shrinking).

    Args:
        agents: Structured agent array (modified in place).
        pop_cfg: Population configuration with L_inf, k_growth.
        rng: NumPy random generator.
        alive_idx: Pre-computed alive indices (optional; computed if None).
    """
    alive = alive_idx if alive_idx is not None else np.where(agents['alive'])[0]
    if len(alive) == 0:
        return

    n = len(alive)

    # Daily aging
    agents['age'][alive] += _INV_365

    # VB growth with dt = 1/365
    old_size = agents['size'][alive].astype(np.float64)
    L_inf = pop_cfg.L_inf
    k = pop_cfg.k_growth
    decay = np.exp(-k * _INV_365)  # scalar, same for all
    deterministic = L_inf - (L_inf - old_size) * decay

    # Growth noise: multiplicative log-normal (daily-scaled σ)
    sigma_daily = 2.0 / np.sqrt(365.0)
    noise = np.exp(rng.normal(0.0, sigma_daily * _INV_365, size=n))
    # Apply noise to the increment only (not to existing size)
    increment = (deterministic - old_size) * noise
    new_size = old_size + np.maximum(increment, 0.0)

    # Clamp: no shrinking
    new_size = np.maximum(new_size, old_size)
    agents['size'][alive] = new_size

    # Vectorized stage assignment (one-directional only)
    # Thresholds: SETTLER (<10) → JUVENILE (<150) → SUBADULT (<400) → ADULT
    current_stage = agents['stage'][alive].copy()

    # Only promote, never demote (use np.maximum against current stage)
    # EGG_LARVA (0) handled externally — skip those
    not_egg = current_stage != Stage.EGG_LARVA

    new_stage = current_stage.copy()
    # Start from smallest → largest threshold
    new_stage[not_egg & (new_size >= 10.0)] = np.maximum(
        new_stage[not_egg & (new_size >= 10.0)], Stage.JUVENILE
    )
    new_stage[not_egg & (new_size >= 150.0)] = np.maximum(
        new_stage[not_egg & (new_size >= 150.0)], Stage.SUBADULT
    )
    new_stage[not_egg & (new_size >= 400.0)] = Stage.ADULT

    agents['stage'][alive] = new_stage


def annual_reproduction(
    agents: np.ndarray,
    genotypes: np.ndarray,
    effect_sizes: np.ndarray,
    habitat_area: float,
    sst: float,
    carrying_capacity: int,
    pop_cfg: PopulationSection,
    rng: np.random.Generator,
    sim_day: int = 0,
    w_od: float = 0.0,  # deprecated, ignored
) -> dict:
    """Full annual reproduction: spawning → fertilization → SRS → settlement.

    Returns diagnostics dict.
    """
    alive = agents['alive']
    stage = agents['stage']
    sex = agents['sex']
    disease = agents['disease_state']

    # Spawning requires: alive, ADULT, healthy (S only — R→S means
    # recovered agents are already S; R check kept for snapshot compat)
    can_spawn = alive & (stage == Stage.ADULT) & (
        (disease == DiseaseState.S) | (disease == DiseaseState.R)
    )
    female_mask = can_spawn & (sex == 0)
    male_mask = can_spawn & (sex == 1)

    females = np.where(female_mask)[0]
    males = np.where(male_mask)[0]

    diag = {
        'n_spawning_females': len(females),
        'n_spawning_males': len(males),
        'fertilization_success': 0.0,
        'n_competent': 0,
        'n_recruits': 0,
    }

    if len(females) == 0 or len(males) == 0:
        return diag

    # Fertilization (Allee effect)
    n_adults = len(females) + len(males)
    male_density = len(males) / habitat_area if habitat_area > 0 else 0.0
    fert = fertilization_success(male_density, pop_cfg.gamma_fert)
    diag['fertilization_success'] = fert

    if fert < 1e-10:
        return diag

    # Total eggs
    sizes = agents['size'][females].astype(np.float64)
    eggs = np.where(
        sizes >= pop_cfg.L_min_repro,
        pop_cfg.F0 * (sizes / pop_cfg.L_ref) ** pop_cfg.fecundity_exp,
        0.0,
    ).sum()

    n_zygotes = int(eggs * fert)
    if n_zygotes <= 0:
        return diag

    # Pelagic larval phase
    pld = pelagic_larval_duration(sst)
    surv = larval_survival(pld)
    n_competent = max(1, int(n_zygotes * surv))
    n_competent = min(n_competent, 100_000)  # memory cap
    diag['n_competent'] = n_competent

    # SRS lottery
    offspring_geno, parent_pairs = srs_reproductive_lottery(
        females=females,
        males=males,
        agents=agents,
        genotypes=genotypes,
        n_offspring_target=n_competent,
        alpha_srs=pop_cfg.alpha_srs,
        F0=pop_cfg.F0,
        L_ref=pop_cfg.L_ref,
        fecundity_exp=pop_cfg.fecundity_exp,
        L_min_repro=pop_cfg.L_min_repro,
        rng=rng,
    )

    if len(offspring_geno) == 0:
        return diag

    # Settlement cue
    n_adults_present = n_adults
    cue_mod = settlement_cue_modifier(n_adults_present)
    effective_settlers = max(0, int(len(offspring_geno) * cue_mod))

    # Beverton-Holt density-dependent recruitment (use full K)
    # Density regulation: BH asymptotes to K; slot availability further limits
    current_alive = int(np.sum(agents['alive']))
    n_recruits = beverton_holt_recruitment(
        effective_settlers, carrying_capacity, pop_cfg.settler_survival,
    )
    # Cap by available slots (dead agent array entries)
    available_slots = max(0, carrying_capacity - current_alive)
    n_recruits = min(n_recruits, available_slots, effective_settlers, len(offspring_geno))

    if n_recruits <= 0:
        return diag

    # Select which settlers survive
    if n_recruits < len(offspring_geno):
        keep_idx = rng.choice(len(offspring_geno), size=n_recruits, replace=False)
        settler_geno = offspring_geno[keep_idx]
    else:
        settler_geno = offspring_geno[:n_recruits]

    # Find dead/empty slots
    dead_slots = np.where(~agents['alive'])[0]
    n_slots = min(n_recruits, len(dead_slots))

    hab_side = np.sqrt(max(habitat_area, 1.0))

    for j in range(n_slots):
        slot = dead_slots[j]
        agents[slot]['alive'] = True
        agents[slot]['x'] = rng.uniform(0, hab_side)
        agents[slot]['y'] = rng.uniform(0, hab_side)
        agents[slot]['heading'] = rng.uniform(0, 2 * np.pi)
        agents[slot]['speed'] = 0.1
        agents[slot]['size'] = 0.5  # mm at settlement
        agents[slot]['age'] = 0.0
        agents[slot]['stage'] = Stage.SETTLER
        agents[slot]['sex'] = rng.integers(0, 2)
        agents[slot]['disease_state'] = DiseaseState.S
        agents[slot]['disease_timer'] = 0
        agents[slot]['node_id'] = 0
        agents[slot]['origin'] = 0
        agents[slot]['cause_of_death'] = 0  # DeathCause.ALIVE
        agents[slot]['settlement_day'] = sim_day  # Phase 11: juvenile immunity

        genotypes[slot] = settler_geno[j]
        agents[slot]['resistance'] = _compute_resistance(
            settler_geno[j], effect_sizes
        )

    diag['n_recruits'] = n_slots
    return diag


# ═══════════════════════════════════════════════════════════════════════
# CONTINUOUS SETTLEMENT (daily cohort settlement)
# ═══════════════════════════════════════════════════════════════════════


def _cohort_settlement_day(c: LarvalCohort) -> int:
    """Key function for sorted insertion by settlement day."""
    return c.spawn_day + int(c.pld_days)


def _insort_cohort(pending: List[LarvalCohort], cohort: LarvalCohort) -> None:
    """Insert cohort into sorted pending list (by settlement_day, ascending)."""
    key = _cohort_settlement_day(cohort)
    # bisect on the key — find insertion point
    lo, hi = 0, len(pending)
    while lo < hi:
        mid = (lo + hi) // 2
        if _cohort_settlement_day(pending[mid]) < key:
            lo = mid + 1
        else:
            hi = mid
    pending.insert(lo, cohort)


def _pop_ready_cohorts(
    pending: List[LarvalCohort], sim_day: int
) -> List[LarvalCohort]:
    """Pop all cohorts with settlement_day <= sim_day from sorted list.

    O(k) where k = number of ready cohorts (typically 0-3 per day).
    The pending list remains sorted after popping.
    """
    if not pending or _cohort_settlement_day(pending[0]) > sim_day:
        return []
    # Find cutoff index using bisect
    # All cohorts with settlement_day <= sim_day are at the front
    cut = 0
    while cut < len(pending) and _cohort_settlement_day(pending[cut]) <= sim_day:
        cut += 1
    ready = pending[:cut]
    del pending[:cut]
    return ready


def settle_daily_cohorts(
    cohorts: List[LarvalCohort],
    agents: np.ndarray,
    genotypes: np.ndarray,
    carrying_capacity: int,
    pop_cfg: PopulationSection,
    rng: np.random.Generator,
    effect_sizes: np.ndarray = None,
    node_id: int = 0,
    habitat_area: float = 1e6,
    sim_day: int = 0,
    w_od: float = 0.0,  # deprecated, ignored
    effects_t: np.ndarray = None,
    effects_c: np.ndarray = None,
    res_slice: slice = None,
    tol_slice: slice = None,
    rec_slice: slice = None,
) -> int:
    """Settle competent larval cohorts into the population.

    Called daily in the simulation loop for cohorts whose PLD has elapsed.
    Applies Beverton-Holt density-dependent recruitment per cohort,
    modulated by settlement cue (adult biofilm Allee effect).

    Args:
        cohorts: List of LarvalCohort objects ready to settle (PLD elapsed).
        agents: Agent structured array (modified in place).
        genotypes: Genotype array (modified in place).
        carrying_capacity: K for this node.
        pop_cfg: Population dynamics configuration (for settler_survival).
        rng: Random number generator.
        effect_sizes: Resistance effect sizes.
            If None, uses ``make_effect_sizes()`` default.
        node_id: Node index for placed agents.
        habitat_area: Habitat area in m² (for spatial placement of settlers).
        sim_day: Current simulation day (stamped as settlement_day for
            juvenile immunity tracking).

    Returns:
        Total number of recruits successfully settled across all cohorts.
    """
    if not cohorts:
        return 0

    if effect_sizes is None:
        effect_sizes = make_effect_sizes()

    total_settled = 0

    for cohort in cohorts:
        if cohort.n_competent <= 0:
            continue

        current_alive = int(np.sum(agents['alive']))
        available_slots = max(0, carrying_capacity - current_alive)
        if available_slots == 0:
            break

        # Count adults for settlement cue modifier (Allee effect)
        n_adults = int(np.sum(
            agents['alive'] & (agents['stage'] == Stage.ADULT)
        ))
        cue_mod = settlement_cue_modifier(n_adults)
        effective_settlers = max(0, int(cohort.n_competent * cue_mod))

        if effective_settlers <= 0:
            continue

        # Beverton-Holt density-dependent recruitment (full K)
        n_recruits = beverton_holt_recruitment(
            effective_settlers, carrying_capacity, pop_cfg.settler_survival,
        )
        n_recruits = min(n_recruits, available_slots, effective_settlers, cohort.n_competent)

        if n_recruits <= 0:
            continue

        # Select which settlers survive (random subsample if needed)
        if n_recruits < cohort.n_competent:
            keep_idx = rng.choice(cohort.n_competent, size=n_recruits, replace=False)
            settler_geno = cohort.genotypes[keep_idx]
        else:
            settler_geno = cohort.genotypes[:n_recruits]

        # Find dead/empty slots in agent array
        dead_slots = np.where(~agents['alive'])[0]
        n_slots = min(n_recruits, len(dead_slots))
        if n_slots <= 0:
            continue

        hab_side = np.sqrt(max(habitat_area, 1.0))
        slots = dead_slots[:n_slots]

        # Batch field assignments (vectorized — same pattern as settle_recruits)
        agents['alive'][slots] = True
        agents['x'][slots] = rng.uniform(0, hab_side, n_slots)
        agents['y'][slots] = rng.uniform(0, hab_side, n_slots)
        agents['heading'][slots] = rng.uniform(0, 2 * np.pi, n_slots)
        agents['speed'][slots] = 0.1
        agents['size'][slots] = 0.5  # mm at settlement
        agents['age'][slots] = 0.0
        agents['stage'][slots] = Stage.SETTLER
        agents['sex'][slots] = rng.integers(0, 2, n_slots)
        agents['disease_state'][slots] = DiseaseState.S
        agents['disease_timer'][slots] = 0
        agents['node_id'][slots] = node_id
        agents['origin'][slots] = 0  # WILD
        agents['cause_of_death'][slots] = 0  # DeathCause.ALIVE
        agents['settlement_day'][slots] = sim_day  # Phase 11: juvenile immunity

        # Batch genotype assignment
        genotypes[slots] = settler_geno[:n_slots]

        # Batch trait score computation for settlers
        n_res = len(effect_sizes)
        allele_means_r = (
            settler_geno[:n_slots, :n_res, :].sum(axis=2).astype(np.float64) * 0.5
        )
        agents['resistance'][slots] = (allele_means_r @ effect_sizes).astype(np.float32)

        # Tolerance and recovery if effect arrays provided
        if effects_t is not None and tol_slice is not None:
            allele_means_t = (
                settler_geno[:n_slots, tol_slice, :].sum(axis=2).astype(np.float64) * 0.5
            )
            agents['tolerance'][slots] = (allele_means_t @ effects_t).astype(np.float32)
        if effects_c is not None and rec_slice is not None:
            allele_means_c = (
                settler_geno[:n_slots, rec_slice, :].sum(axis=2).astype(np.float64) * 0.5
            )
            agents['recovery_ability'][slots] = (allele_means_c @ effects_c).astype(np.float32)

        total_settled += n_slots

    return total_settled


# ═══════════════════════════════════════════════════════════════════════
# COUPLED SIMULATION RESULT
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class CoupledSimResult:
    """Results from a coupled disease + population + genetics simulation."""
    n_years: int = 0
    # Annual timeseries (length = n_years)
    yearly_pop: Optional[np.ndarray] = None
    yearly_adults: Optional[np.ndarray] = None
    yearly_recruits: Optional[np.ndarray] = None
    yearly_natural_deaths: Optional[np.ndarray] = None
    yearly_disease_deaths: Optional[np.ndarray] = None
    yearly_mean_resistance: Optional[np.ndarray] = None
    yearly_mean_tolerance: Optional[np.ndarray] = None
    yearly_mean_recovery: Optional[np.ndarray] = None
    yearly_fert_success: Optional[np.ndarray] = None

    # Genetics tracking (Phase 7)
    yearly_allele_freq_top3: Optional[np.ndarray] = None   # (n_years, 3)
    yearly_va: Optional[np.ndarray] = None                  # (n_years,) additive variance
    yearly_va_tolerance: Optional[np.ndarray] = None        # (n_years,)
    yearly_va_recovery: Optional[np.ndarray] = None         # (n_years,)
    pre_epidemic_allele_freq: Optional[np.ndarray] = None   # (N_LOCI,) snapshot
    post_epidemic_allele_freq: Optional[np.ndarray] = None  # (N_LOCI,) snapshot

    # Daily timeseries (length = n_years * 365)
    daily_pop: Optional[np.ndarray] = None
    daily_infected: Optional[np.ndarray] = None
    daily_vibrio: Optional[np.ndarray] = None

    # Death accounting by cause (annual timeseries)
    yearly_senescence_deaths: Optional[np.ndarray] = None

    # Spawning event tracking (daily)
    daily_spawning_counts: Optional[np.ndarray] = None  # (total_days,) spawners per day

    # Pathogen evolution (annual timeseries)
    yearly_mean_virulence: Optional[np.ndarray] = None          # (n_years,)
    yearly_virulence_new_infections: Optional[np.ndarray] = None  # (n_years,)
    yearly_virulence_of_deaths: Optional[np.ndarray] = None     # (n_years,)

    # Summary
    initial_pop: int = 0
    final_pop: int = 0
    min_pop: int = 0
    min_pop_year: int = 0
    total_disease_deaths: int = 0
    total_natural_deaths: int = 0
    total_senescence_deaths: int = 0
    peak_mortality_fraction: float = 0.0
    recovered: bool = False


# ═══════════════════════════════════════════════════════════════════════
# COUPLED SINGLE-NODE SIMULATION
# ═══════════════════════════════════════════════════════════════════════

def run_coupled_simulation(
    n_individuals: int = 500,
    carrying_capacity: int = 500,
    habitat_area: float = 10000.0,
    T_celsius: float = 15.0,
    salinity: float = 30.0,
    phi_k: float = 0.02,
    n_years: int = 10,
    disease_year: int = 3,
    initial_infected: int = 5,
    seed: int = 42,
    record_daily: bool = False,
    config: Optional[SimulationConfig] = None,
    sst_by_day: Optional[np.ndarray] = None,
    perf: Optional['PerfMonitor'] = None,
) -> CoupledSimResult:
    """Run a coupled disease ↔ population simulation on a single node.

    The simulation runs for n_years:
      - Years 0 to disease_year-1: disease-free spinup (reach equilibrium)
      - Year disease_year: pathogen introduced (initial infections seeded)
      - Remaining years: epidemic + recovery

    Daily loop (365 days/year):
      1. Update SST (constant or from sst_by_day)
      2. Disease step: Vibrio dynamics, transmission, progression
      3. Track daily metrics

    Annual step (end of each year):
      1. Natural mortality (stage-specific, senescence)
      2. Growth & aging (VB, stage transitions)
      3. Reproduction (spawning, fertilization, SRS, settlement)
      4. Disease year: seed initial infections
      5. Record annual metrics

    Args:
        n_individuals: Starting population size.
        carrying_capacity: K for the node.
        habitat_area: Habitat area (m²).
        T_celsius: Constant SST (°C), unless sst_by_day provided.
        salinity: Constant salinity (psu).
        phi_k: Flushing rate (d⁻¹).
        n_years: Total years to simulate.
        disease_year: Year to introduce disease (0-indexed).
        initial_infected: Number of initially infected individuals.
        seed: RNG seed.
        record_daily: If True, record daily pop/infected/vibrio.
        config: Optional SimulationConfig; uses default if None.
        sst_by_day: Optional SST array (n_years*365,) for variable temperature.

    Returns:
        CoupledSimResult with annual and optionally daily timeseries.
    """
    if config is None:
        config = default_config()
    if perf is None:
        perf = PerfMonitor(enabled=False)
    perf.start()

    pop_cfg = config.population
    dis_cfg = config.disease
    pe_cfg = (config.pathogen_evolution
              if hasattr(config, 'pathogen_evolution')
              and config.pathogen_evolution.enabled
              else None)

    rng = np.random.default_rng(seed)

    # Initialize effect sizes for all three traits
    gen_cfg = config.genetics
    effect_sizes = make_effect_sizes(gen_cfg.effect_size_seed, n_loci=gen_cfg.n_resistance)
    effects_t = initialize_trait_effect_sizes(
        np.random.default_rng(gen_cfg.effect_size_seed + 1),
        gen_cfg.n_tolerance, total_weight=1.0,
    )
    effects_c = initialize_trait_effect_sizes(
        np.random.default_rng(gen_cfg.effect_size_seed + 2),
        gen_cfg.n_recovery, total_weight=1.0,
    )
    res_slice, tol_slice, rec_slice = trait_slices(
        gen_cfg.n_resistance, gen_cfg.n_tolerance, gen_cfg.n_recovery,
    )

    # Allocate arrays — extra capacity for population fluctuation
    max_agents = max(int(carrying_capacity * 2.5), n_individuals * 3)
    agents, genotypes = initialize_population(
        n_individuals=n_individuals,
        max_agents=max_agents,
        habitat_area=habitat_area,
        effect_sizes=effect_sizes,
        pop_cfg=pop_cfg,
        rng=rng,
        q_init=0.05,
        genetics_cfg=gen_cfg,
        effects_t=effects_t,
        effects_c=effects_c,
    )

    # Disease state
    node_disease = NodeDiseaseState(node_id=0)

    # Vibrio: steady-state background if ubiquitous
    if dis_cfg.scenario == "ubiquitous":
        env = environmental_vibrio(T_celsius, salinity, dis_cfg)
        xi = vibrio_decay_rate(T_celsius)
        node_disease.vibrio_concentration = env / (xi + phi_k) if (xi + phi_k) > 0 else 0.0
    else:
        node_disease.vibrio_concentration = 0.0

    # Track whether disease is active
    disease_active = False
    cumulative_disease_deaths = 0

    # ── Spawning system variables ────────────────────────────────────
    spawning_enabled = config.spawning is not None
    accumulated_cohorts = []  # List of LarvalCohort objects collected during season
    previous_in_season = False  # Track season transitions
    spawning_diagnostics = {
        'n_spawning_events': 0,
        'n_cohorts': 0,
        'total_larvae': 0,
    }

    # ── Allocate result arrays ────────────────────────────────────────
    total_days = n_years * DAYS_PER_YEAR

    # Spawning event tracking (daily)
    daily_spawning_counts = np.zeros(total_days, dtype=np.int32)

    yearly_pop = np.zeros(n_years, dtype=np.int32)
    yearly_adults = np.zeros(n_years, dtype=np.int32)
    yearly_recruits = np.zeros(n_years, dtype=np.int32)
    yearly_natural_deaths = np.zeros(n_years, dtype=np.int32)
    yearly_disease_deaths = np.zeros(n_years, dtype=np.int32)
    yearly_senescence_deaths = np.zeros(n_years, dtype=np.int32)
    yearly_mean_resistance = np.zeros(n_years, dtype=np.float64)
    yearly_mean_tolerance = np.zeros(n_years, dtype=np.float64)
    yearly_mean_recovery = np.zeros(n_years, dtype=np.float64)
    yearly_fert_success = np.zeros(n_years, dtype=np.float64)

    # Genetics tracking (Phase 7)
    yearly_allele_freq_top3 = np.zeros((n_years, 3), dtype=np.float64)
    yearly_va = np.zeros(n_years, dtype=np.float64)
    yearly_va_tolerance = np.zeros(n_years, dtype=np.float64)
    yearly_va_recovery = np.zeros(n_years, dtype=np.float64)
    pre_epidemic_allele_freq = None
    post_epidemic_allele_freq = None
    epidemic_started_year = None  # track when epidemic begins for snapshot timing

    # Pathogen evolution virulence tracking
    if pe_cfg is not None:
        yearly_mean_v = np.zeros(n_years, dtype=np.float64)
        yearly_v_new_inf = np.zeros(n_years, dtype=np.float64)
        yearly_v_deaths = np.zeros(n_years, dtype=np.float64)
        year_new_inf_v_accum = []   # list of v values for new infections this year
        year_death_v_accum = []     # list of v values for disease deaths this year
    else:
        yearly_mean_v = None
        yearly_v_new_inf = None
        yearly_v_deaths = None

    if record_daily:
        daily_pop = np.zeros(total_days, dtype=np.int32)
        daily_infected = np.zeros(total_days, dtype=np.int32)
        daily_vibrio = np.zeros(total_days, dtype=np.float64)
    else:
        daily_pop = daily_infected = daily_vibrio = None

    initial_pop = int(np.sum(agents['alive']))
    min_pop = initial_pop
    min_pop_year = 0
    peak_mort_frac = 0.0

    # ── Daily demographic accumulators (continuous mortality) ────────
    daily_nat_deaths_accum = 0
    daily_senes_deaths_accum = 0

    # ── Main simulation loop ─────────────────────────────────────────
    for year in range(n_years):
        year_disease_deaths_before = cumulative_disease_deaths
        yearly_recruits_accum = 0  # Daily settlement accumulator for this year

        # Reset daily demographic accumulators for this year
        daily_nat_deaths_accum = 0
        daily_senes_deaths_accum = 0

        # ── PHASE A: DAILY DISEASE LOOP (365 days) ───────────────────
        for day in range(DAYS_PER_YEAR):
            sim_day = year * DAYS_PER_YEAR + day
            global_day = sim_day

            # Get today's temperature
            if sst_by_day is not None and sim_day < len(sst_by_day):
                today_T = float(sst_by_day[sim_day])
            else:
                today_T = T_celsius

            # Disease step (only when active)
            if disease_active:
                with perf.track("disease"):
                    deaths_before = node_disease.cumulative_deaths
                    node_disease = daily_disease_update(
                        agents=agents,
                        node_state=node_disease,
                        T_celsius=today_T,
                        salinity=salinity,
                        phi_k=phi_k,
                        dispersal_input=0.0,
                        day=global_day,
                        cfg=dis_cfg,
                        rng=rng,
                        pe_cfg=pe_cfg,
                    )
                    new_deaths = node_disease.cumulative_deaths - deaths_before
                    cumulative_disease_deaths += new_deaths

            # ── Daily settlement (continuous — settle cohorts that reached PLD) ─
            if spawning_enabled and accumulated_cohorts:
                ready = [c for c in accumulated_cohorts
                         if (sim_day - c.spawn_day) >= c.pld_days]
                accumulated_cohorts = [c for c in accumulated_cohorts
                                       if (sim_day - c.spawn_day) < c.pld_days]
                if ready:
                    with perf.track("settlement"):
                        n_settled = settle_daily_cohorts(
                            ready, agents, genotypes, carrying_capacity,
                            pop_cfg, rng, effect_sizes,
                            habitat_area=habitat_area,
                            sim_day=sim_day,
                            effects_t=effects_t,
                            effects_c=effects_c,
                            res_slice=res_slice,
                            tol_slice=tol_slice,
                            rec_slice=rec_slice,
                        )
                        yearly_recruits_accum += n_settled

            # ── Spawning step (if enabled) ────────────────────────────
            if spawning_enabled:
                day_of_year = (day + 1)  # Convert 0-based to 1-based day of year
                currently_in_season = in_spawning_season(day_of_year, 
                    config.spawning.season_start_doy, 
                    config.spawning.season_end_doy)
                
                # Reset spawning season when transitioning out of season
                if previous_in_season and not currently_in_season:
                    reset_spawning_season(agents)
                
                # Daily spawning step during season
                if currently_in_season:
                    with perf.track("spawning"):
                        cohorts_today = spawning_step(
                            agents=agents,
                            genotypes=genotypes,
                            day_of_year=day_of_year,
                            node_latitude=48.0,  # Default latitude for single-node
                            spawning_config=config.spawning,
                            disease_config=config.disease,
                            rng=rng,
                            current_sim_day=sim_day,
                            current_sst=today_T,
                        )
                    if cohorts_today:
                        accumulated_cohorts.extend(cohorts_today)
                        spawning_diagnostics['n_spawning_events'] += len(cohorts_today)
                        spawning_diagnostics['n_cohorts'] += len(cohorts_today)
                        spawning_diagnostics['total_larvae'] += sum(c.n_competent for c in cohorts_today)
                    
                    # Count spawners today: agents whose last_spawn_day == today's DOY
                    alive_spawned = (
                        agents['alive'] &
                        (agents['last_spawn_day'] == day_of_year)
                    )
                    daily_spawning_counts[sim_day] = int(np.sum(alive_spawned))
                
                # Tick down immunosuppression timers
                with perf.track("immunosuppression_tick"):
                    immuno_mask = agents['immunosuppression_timer'] > 0
                    agents['immunosuppression_timer'][immuno_mask] -= 1
                
                previous_in_season = currently_in_season

            # Daily demographics (continuous mortality + growth)
            # Compute alive_idx once, share across both calls
            _alive_idx_demo = np.where(agents['alive'])[0]
            with perf.track('daily_mortality'):
                n_mort, n_senes = daily_natural_mortality(agents, pop_cfg, rng, alive_idx=_alive_idx_demo)
                daily_nat_deaths_accum += n_mort
                daily_senes_deaths_accum += n_senes

            with perf.track('daily_growth'):
                if n_mort > 0:
                    _alive_idx_demo = np.where(agents['alive'])[0]
                daily_growth_and_aging(agents, pop_cfg, rng, alive_idx=_alive_idx_demo)

            # Daily recording
            if record_daily:
                alive_mask = agents['alive'].astype(bool)
                ds = agents['disease_state']
                daily_pop[sim_day] = int(np.sum(alive_mask))
                daily_infected[sim_day] = int(np.sum(
                    alive_mask & ((ds == DiseaseState.I1) | (ds == DiseaseState.I2))
                ))
                daily_vibrio[sim_day] = node_disease.vibrio_concentration

        # ── PHASE B: ANNUAL DEMOGRAPHIC UPDATE ───────────────────────

        # B1. Natural mortality now runs daily (Phase A) — just record accumulators
        n_nat_dead = daily_nat_deaths_accum
        n_senes_dead = daily_senes_deaths_accum
        # B2. Growth now runs daily (Phase A)

        # B3. Record pre-reproduction population
        alive_mask = agents['alive'].astype(bool)
        pop_now = int(np.sum(alive_mask))

        # B4. Reproduction
        if spawning_enabled:
            # Continuous settlement already happened daily in Phase A.
            # Settle any remaining cohorts that matured on the last day.
            if accumulated_cohorts:
                end_of_year_ready = [c for c in accumulated_cohorts
                                     if ((year + 1) * DAYS_PER_YEAR - c.spawn_day) >= c.pld_days]
                accumulated_cohorts = [c for c in accumulated_cohorts
                                       if ((year + 1) * DAYS_PER_YEAR - c.spawn_day) < c.pld_days]
                if end_of_year_ready:
                    end_of_year_sim_day = (year + 1) * DAYS_PER_YEAR - 1
                    n_eoy = settle_daily_cohorts(
                        end_of_year_ready, agents, genotypes,
                        carrying_capacity, pop_cfg, rng, effect_sizes,
                        habitat_area=habitat_area,
                        sim_day=end_of_year_sim_day,
                        effects_t=effects_t,
                        effects_c=effects_c,
                        res_slice=res_slice,
                        tol_slice=tol_slice,
                        rec_slice=rec_slice,
                    )
                    yearly_recruits_accum += n_eoy

            total_competent_larvae = spawning_diagnostics['total_larvae']
            repro_diag = {
                'n_spawning_females': spawning_diagnostics['n_spawning_events'],
                'n_spawning_males': spawning_diagnostics['n_spawning_events'],
                'fertilization_success': 1.0 if total_competent_larvae > 0 else 0.0,
                'n_competent': total_competent_larvae,
                'n_recruits': yearly_recruits_accum,
            }

            # Reset diagnostics for next year (keep pending cohorts across years)
            spawning_diagnostics = {'n_spawning_events': 0, 'n_cohorts': 0, 'total_larvae': 0}

        elif pop_now > 0:
            # Fall back to annual reproduction (backward compatibility)
            with perf.track("recruitment"):
                repro_diag = annual_reproduction(
                    agents=agents,
                    genotypes=genotypes,
                    effect_sizes=effect_sizes,
                    habitat_area=habitat_area,
                    sst=T_celsius,
                    carrying_capacity=carrying_capacity,
                    pop_cfg=pop_cfg,
                    rng=rng,
                    sim_day=(year + 1) * DAYS_PER_YEAR - 1,
                )
        else:
            repro_diag = {
                'n_spawning_females': 0,
                'n_spawning_males': 0,
                'fertilization_success': 0.0,
                'n_competent': 0,
                'n_recruits': 0,
            }

        # B5. Seed disease at disease_year
        if year == disease_year and not disease_active:
            # ── Pre-epidemic allele frequency snapshot (Phase 7) ──────
            alive_mask_pre = agents['alive'].astype(bool)
            if int(np.sum(alive_mask_pre)) > 0:
                pre_epidemic_allele_freq = compute_allele_frequencies(
                    genotypes, alive_mask_pre
                )

            disease_active = True
            epidemic_started_year = year
            alive_idx = np.where(agents['alive'])[0]
            n_to_infect = min(initial_infected, len(alive_idx))
            if n_to_infect > 0:
                infect_idx = rng.choice(alive_idx, size=n_to_infect,
                                        replace=False)
                mu_EI1 = arrhenius(dis_cfg.mu_EI1_ref, dis_cfg.Ea_EI1,
                                   T_celsius)
                from sswd_evoepi.disease import sample_stage_duration, K_SHAPE_E
                for idx in infect_idx:
                    agents['disease_state'][idx] = DiseaseState.E
                    agents['disease_timer'][idx] = sample_stage_duration(
                        mu_EI1, K_SHAPE_E, rng
                    )
                    if pe_cfg is not None:
                        agents['pathogen_virulence'][idx] = pe_cfg.v_init

        # ── Post-epidemic allele frequency snapshot (Phase 7) ────────
        # Take snapshot 2 years after epidemic start (when initial wave is over)
        if (epidemic_started_year is not None
                and year == epidemic_started_year + 2
                and post_epidemic_allele_freq is None):
            alive_mask_post = agents['alive'].astype(bool)
            if int(np.sum(alive_mask_post)) > 1:
                post_epidemic_allele_freq = compute_allele_frequencies(
                    genotypes, alive_mask_post
                )

        # ── Record annual metrics ────────────────────────────────────
        alive_mask = agents['alive'].astype(bool)
        pop_after = int(np.sum(alive_mask))
        adult_mask = alive_mask & (agents['stage'] == Stage.ADULT)
        n_adults = int(np.sum(adult_mask))

        yearly_pop[year] = pop_after
        yearly_adults[year] = n_adults
        yearly_recruits[year] = repro_diag.get('n_recruits', 0)
        yearly_natural_deaths[year] = n_nat_dead
        yearly_senescence_deaths[year] = n_senes_dead
        yearly_disease_deaths[year] = (
            cumulative_disease_deaths - year_disease_deaths_before
        )
        yearly_fert_success[year] = repro_diag.get('fertilization_success', 0.0)

        # Mean trait scores of living individuals
        if pop_after > 0:
            yearly_mean_resistance[year] = float(
                np.mean(agents['resistance'][alive_mask])
            )
            yearly_mean_tolerance[year] = float(
                np.mean(agents['tolerance'][alive_mask])
            )
            yearly_mean_recovery[year] = float(
                np.mean(agents['recovery_ability'][alive_mask])
            )
        else:
            yearly_mean_resistance[year] = 0.0
            yearly_mean_tolerance[year] = 0.0
            yearly_mean_recovery[year] = 0.0

        # ── Pathogen evolution virulence tracking ─────────────────────
        if pe_cfg is not None:
            # Time-weighted mean virulence across all infected agent-days
            if node_disease.virulence_count_daily > 0:
                yearly_mean_v[year] = (
                    node_disease.virulence_sum_daily
                    / node_disease.virulence_count_daily
                )

            # Mean virulence of new infections this year (from accumulator)
            if node_disease.virulence_count_new_infections > 0:
                yearly_v_new_inf[year] = (
                    node_disease.virulence_sum_new_infections
                    / node_disease.virulence_count_new_infections
                )

            # Mean virulence of disease deaths this year (from accumulator)
            if node_disease.virulence_count_deaths > 0:
                yearly_v_deaths[year] = (
                    node_disease.virulence_sum_deaths
                    / node_disease.virulence_count_deaths
                )

            # Reset accumulators for next year
            node_disease.virulence_sum_new_infections = 0.0
            node_disease.virulence_count_new_infections = 0
            node_disease.virulence_sum_deaths = 0.0
            node_disease.virulence_count_deaths = 0
            node_disease.virulence_sum_daily = 0.0
            node_disease.virulence_count_daily = 0

        # ── Genetics tracking (Phase 7) ──────────────────────────────
        if pop_after > 0:
            allele_freq = compute_allele_frequencies(genotypes, alive_mask)
            yearly_allele_freq_top3[year] = allele_freq[:3]
            yearly_va[year] = compute_additive_variance(
                allele_freq, effect_sizes, res_slice
            )
            yearly_va_tolerance[year] = compute_additive_variance(
                allele_freq, effects_t, tol_slice
            )
            yearly_va_recovery[year] = compute_additive_variance(
                allele_freq, effects_c, rec_slice
            )

        # Track minimum population
        if pop_after < min_pop:
            min_pop = pop_after
            min_pop_year = year

        # Peak disease mortality fraction (this year)
        year_dd = yearly_disease_deaths[year]
        year_start_pop = yearly_pop[year - 1] if year > 0 else initial_pop
        if year_start_pop > 0:
            mort_frac = year_dd / year_start_pop
            peak_mort_frac = max(peak_mort_frac, mort_frac)

    # ── Compile result ────────────────────────────────────────────────
    final_pop = int(np.sum(agents['alive']))
    recovered = final_pop > n_individuals * 0.5

    perf.stop()

    result = CoupledSimResult(
        n_years=n_years,
        yearly_pop=yearly_pop,
        yearly_adults=yearly_adults,
        yearly_recruits=yearly_recruits,
        yearly_natural_deaths=yearly_natural_deaths,
        yearly_disease_deaths=yearly_disease_deaths,
        yearly_senescence_deaths=yearly_senescence_deaths,
        yearly_mean_resistance=yearly_mean_resistance,
        yearly_mean_tolerance=yearly_mean_tolerance,
        yearly_mean_recovery=yearly_mean_recovery,
        yearly_fert_success=yearly_fert_success,
        # Genetics (Phase 7)
        yearly_allele_freq_top3=yearly_allele_freq_top3,
        yearly_va=yearly_va,
        yearly_va_tolerance=yearly_va_tolerance,
        yearly_va_recovery=yearly_va_recovery,
        pre_epidemic_allele_freq=pre_epidemic_allele_freq,
        post_epidemic_allele_freq=post_epidemic_allele_freq,
        # Pathogen evolution
        yearly_mean_virulence=yearly_mean_v,
        yearly_virulence_new_infections=yearly_v_new_inf,
        yearly_virulence_of_deaths=yearly_v_deaths,
        # Daily
        daily_pop=daily_pop,
        daily_infected=daily_infected,
        daily_vibrio=daily_vibrio,
        # Spawning event tracking
        daily_spawning_counts=daily_spawning_counts,
        # Summary
        initial_pop=initial_pop,
        final_pop=final_pop,
        min_pop=min_pop,
        min_pop_year=min_pop_year,
        total_disease_deaths=cumulative_disease_deaths,
        total_natural_deaths=int(np.sum(yearly_natural_deaths)),
        total_senescence_deaths=int(np.sum(yearly_senescence_deaths)),
        peak_mortality_fraction=peak_mort_frac,
        recovered=recovered,
    )
    return result


# ═══════════════════════════════════════════════════════════════════════
# SPATIAL MULTI-NODE SIMULATION (Phase 9)
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class SpatialSimResult:
    """Results from a multi-node spatial simulation."""
    n_years: int = 0
    n_nodes: int = 0
    node_names: Optional[List[str]] = None
    node_K: Optional[np.ndarray] = None
    snapshot_recorder: Optional[object] = None  # SnapshotRecorder if enabled
    # Per-node, per-year: shape (n_nodes, n_years)
    yearly_pop: Optional[np.ndarray] = None
    yearly_adults: Optional[np.ndarray] = None
    yearly_recruits: Optional[np.ndarray] = None
    yearly_natural_deaths: Optional[np.ndarray] = None
    yearly_disease_deaths: Optional[np.ndarray] = None
    yearly_recoveries: Optional[np.ndarray] = None
    yearly_senescence_deaths: Optional[np.ndarray] = None
    yearly_mean_resistance: Optional[np.ndarray] = None
    yearly_mean_tolerance: Optional[np.ndarray] = None
    yearly_mean_recovery: Optional[np.ndarray] = None
    yearly_vibrio_max: Optional[np.ndarray] = None
    # Per-node genetics: shape (n_nodes, n_years)
    yearly_va: Optional[np.ndarray] = None
    yearly_va_tolerance: Optional[np.ndarray] = None
    yearly_va_recovery: Optional[np.ndarray] = None
    yearly_allele_freq_top3: Optional[np.ndarray] = None  # (n_nodes, n_years, 3)
    yearly_ne_ratio: Optional[np.ndarray] = None
    # Pre/post-epidemic snapshots per node: shape (n_nodes, N_LOCI)
    pre_epidemic_allele_freq: Optional[np.ndarray] = None
    post_epidemic_allele_freq: Optional[np.ndarray] = None
    # Per-year totals
    yearly_total_pop: Optional[np.ndarray] = None
    yearly_total_larvae_dispersed: Optional[np.ndarray] = None
    # Peak disease prevalence per node (fraction infected at peak)
    peak_disease_prevalence: Optional[np.ndarray] = None
    # Pathogen evolution per-node virulence tracking: (n_nodes, n_years)
    yearly_mean_virulence: Optional[np.ndarray] = None
    # Spawning event tracking (daily, per-node)
    daily_spawning_counts: Optional[np.ndarray] = None  # (n_nodes, total_days)
    # Summary
    initial_total_pop: int = 0
    final_total_pop: int = 0
    disease_year: Optional[int] = None
    seed: int = 0


def run_spatial_simulation(
    network: 'MetapopulationNetwork',
    n_years: int = 10,
    disease_year: int = 3,
    initial_infected_per_node: int = 5,
    seed: int = 42,
    config: Optional[SimulationConfig] = None,
    progress_callback=None,
    snapshot_recorder=None,
) -> SpatialSimResult:
    """Run multi-node spatial simulation with full genetics tracking.

    Each node runs independently within a daily timestep. At the end
    of each daily step, pathogen is exchanged between nodes via the D
    matrix. At the end of each annual step, larvae are dispersed via
    the C matrix.

    Full integration (Phase 10):
      - 5 nodes, each with: population + disease + genetics + environmental forcing
      - Daily: update_environment → disease_step → pathogen_dispersal
      - Annual: natural_mortality → growth_aging → reproduction_with_srs →
                larval_dispersal → genetics_recording
      - Genetics tracking: allele frequencies, EF1A, V_A, Ne/N per node per year
      - Pre/post-epidemic allele frequency snapshots

    Args:
        network: MetapopulationNetwork (from spatial.build_network).
        n_years: Number of years to simulate.
        disease_year: Year to introduce disease (default 3; None = no disease).
        initial_infected_per_node: Infections per node at disease_year.
        seed: RNG seed.
        config: SimulationConfig; uses default if None.
        progress_callback: Optional callable(year, n_years) for progress.

    Returns:
        SpatialSimResult with per-node annual timeseries + genetics.
    """
    from sswd_evoepi.spatial import (
        SpatialNode,
        MetapopulationNetwork,
        distribute_larvae,
        pathogen_dispersal_step,
    )
    from sswd_evoepi.environment import (
        sst_with_trend, seasonal_flushing, generate_satellite_sst_series,
        generate_yearly_sst_series,
    )
    from sswd_evoepi.movement import daily_movement, InfectedDensityGrid

    if config is None:
        config = default_config()

    pop_cfg = config.population
    dis_cfg = config.disease
    pe_cfg = (config.pathogen_evolution
              if hasattr(config, 'pathogen_evolution')
              and config.pathogen_evolution.enabled
              else None)
    mov_cfg = config.movement
    mov_enabled = mov_cfg.enabled
    spatial_tx = mov_enabled and mov_cfg.spatial_transmission

    rng = np.random.default_rng(seed)
    N = network.n_nodes

    # ── Initialize populations at each node ──────────────────────────
    gen_cfg = config.genetics
    effect_sizes = make_effect_sizes(gen_cfg.effect_size_seed, n_loci=gen_cfg.n_resistance)
    effects_t = initialize_trait_effect_sizes(
        np.random.default_rng(gen_cfg.effect_size_seed + 1),
        gen_cfg.n_tolerance, total_weight=1.0,
    )
    effects_c = initialize_trait_effect_sizes(
        np.random.default_rng(gen_cfg.effect_size_seed + 2),
        gen_cfg.n_recovery, total_weight=1.0,
    )
    res_slice, tol_slice, rec_slice = trait_slices(
        gen_cfg.n_resistance, gen_cfg.n_tolerance, gen_cfg.n_recovery,
    )

    for node in network.nodes:
        nd = node.definition
        K = nd.carrying_capacity
        max_agents = max(int(K * 2.5), K + 500)
        agents, geno = initialize_population(
            n_individuals=K,
            max_agents=max_agents,
            habitat_area=nd.habitat_area,
            effect_sizes=effect_sizes,
            pop_cfg=pop_cfg,
            rng=rng,
            q_init=0.05,
            genetics_cfg=gen_cfg,
            effects_t=effects_t,
            effects_c=effects_c,
        )
        # Stamp node_id on all agents
        agents['node_id'][:K] = nd.node_id
        node.agents = agents
        node.genotypes = geno

        # Init vibrio from environmental background
        if dis_cfg.scenario == "ubiquitous":
            env = environmental_vibrio(nd.mean_sst, nd.salinity, dis_cfg)
            xi = vibrio_decay_rate(nd.mean_sst)
            phi = nd.flushing_rate
            node.vibrio_concentration = (
                env / (xi + phi) if (xi + phi) > 0 else 0.0
            )
        else:
            node.vibrio_concentration = 0.0

    # ── Allocate result arrays ───────────────────────────────────────
    yearly_pop = np.zeros((N, n_years), dtype=np.int32)
    yearly_adults = np.zeros((N, n_years), dtype=np.int32)
    yearly_recruits = np.zeros((N, n_years), dtype=np.int32)
    yearly_nat_deaths = np.zeros((N, n_years), dtype=np.int32)
    yearly_dis_deaths = np.zeros((N, n_years), dtype=np.int32)
    yearly_recoveries = np.zeros((N, n_years), dtype=np.int32)
    yearly_senes_deaths = np.zeros((N, n_years), dtype=np.int32)
    yearly_mean_r = np.zeros((N, n_years), dtype=np.float64)
    yearly_mean_t = np.zeros((N, n_years), dtype=np.float64)
    yearly_mean_c = np.zeros((N, n_years), dtype=np.float64)
    yearly_vibrio_max = np.zeros((N, n_years), dtype=np.float64)
    yearly_total_pop = np.zeros(n_years, dtype=np.int32)
    yearly_total_larvae = np.zeros(n_years, dtype=np.int32)

    # Genetics tracking per node per year
    yearly_va = np.zeros((N, n_years), dtype=np.float64)
    yearly_va_t = np.zeros((N, n_years), dtype=np.float64)
    yearly_va_c = np.zeros((N, n_years), dtype=np.float64)
    yearly_allele_freq_top3 = np.zeros((N, n_years, 3), dtype=np.float64)
    yearly_ne_ratio = np.zeros((N, n_years), dtype=np.float64)

    # Pre/post-epidemic allele frequency snapshots (n_nodes, N_LOCI)
    from sswd_evoepi.types import N_LOCI as _N_LOCI
    pre_epidemic_af = np.zeros((N, _N_LOCI), dtype=np.float64)
    post_epidemic_af = np.zeros((N, _N_LOCI), dtype=np.float64)
    pre_snapshot_taken = False
    post_snapshot_taken = False
    epidemic_started_year = None

    # Peak disease prevalence tracker per node
    peak_disease_prev = np.zeros(N, dtype=np.float64)

    # Pathogen evolution per-node virulence tracking
    if pe_cfg is not None:
        yearly_mean_v_spatial = np.zeros((N, n_years), dtype=np.float64)
    else:
        yearly_mean_v_spatial = None

    initial_total = network.total_population()
    node_names = [n.definition.name for n in network.nodes]
    node_K = np.array([n.definition.carrying_capacity for n in network.nodes])

    # Per-node disease tracking
    node_disease_states = []
    disease_active_flags = [False] * N
    cumulative_dis_deaths = [0] * N
    cumulative_recoveries = [0] * N

    for i in range(N):
        node_disease_states.append(NodeDiseaseState(node_id=i))
        node_disease_states[i].vibrio_concentration = (
            network.nodes[i].vibrio_concentration
        )

    # ── Spatial transmission grids (one per node) ────────────────────
    density_grids = [None] * N
    if spatial_tx:
        for i, node in enumerate(network.nodes):
            hab_side = np.sqrt(max(node.definition.habitat_area, 1.0))
            density_grids[i] = InfectedDensityGrid(
                habitat_side=hab_side,
                cell_size=mov_cfg.cell_size,
            )

    # ── Spawning event tracking (daily, per-node) ─────────────────
    total_sim_days = n_years * DAYS_PER_YEAR
    daily_spawning_counts_spatial = np.zeros((N, total_sim_days), dtype=np.int32)

    # ── Continuous settlement: per-node pending cohort lists ────────
    # Each node accumulates LarvalCohort objects from daily spawning.
    # At year-end, unsettled cohorts are collected for C matrix dispersal.
    # Dispersed cohorts land in pending_cohorts[dest] and settle daily
    # when PLD elapses.
    #
    # OPTIMIZATION (Phase 7): Lists are kept sorted by settlement_day
    # (ascending). Daily check pops from front while settlement_day <=
    # sim_day — O(1) amortized instead of O(n) filter per day.
    pending_cohorts: List[List[LarvalCohort]] = [[] for _ in range(N)]
    spawning_enabled = config.spawning is not None

    # Per-node spawning season tracking
    previous_in_season = [False] * N
    yearly_recruits_accum = [0] * N  # track daily settlement per node per year

    # ── Per-node daily mortality accumulators ────────────────────────
    node_daily_nat_deaths = np.zeros(N, dtype=np.int32)
    node_daily_senes_deaths = np.zeros(N, dtype=np.int32)

    # ── Pre-generate SST timeseries if configured (satellite or monthly) ──
    use_precomputed_sst = config.simulation.sst_source in ('satellite', 'monthly')
    satellite_sst_by_node: List[Optional[np.ndarray]] = [None] * N
    if use_precomputed_sst:
        for i, node in enumerate(network.nodes):
            nd = node.definition
            if config.simulation.sst_source == 'monthly':
                satellite_sst_by_node[i] = generate_yearly_sst_series(
                    node_name=nd.name,
                    start_year=config.simulation.sst_start_year,
                    n_years=n_years,
                    data_dir=config.simulation.sst_data_dir,
                    climatology_dir=config.simulation.sst_data_dir,
                    trend_per_year=nd.sst_trend,
                    reference_year=2015,
                    scenario=config.simulation.sst_scenario,
                    projection_dir=config.simulation.sst_projection_dir,
                    obs_end_year=config.simulation.sst_obs_end_year,
                )
            else:
                satellite_sst_by_node[i] = generate_satellite_sst_series(
                    n_years=n_years,
                    start_year=2000,  # matches start_year below
                    node_name=nd.name,
                    trend_per_year=nd.sst_trend,
                    reference_year=2015,
                    data_dir=config.simulation.sst_data_dir,
                )

    # ── Main simulation loop ─────────────────────────────────────────
    start_year = 2000  # reference year for SST

    for year in range(n_years):
        cal_year = start_year + year
        year_dd_before = list(cumulative_dis_deaths)
        year_rec_before = list(cumulative_recoveries)

        # Reset per-year accumulators
        yearly_recruits_accum = [0] * N
        node_daily_nat_deaths[:] = 0
        node_daily_senes_deaths[:] = 0

        # Track max vibrio per node this year
        max_vibrio_year = np.zeros(N, dtype=np.float64)

        # ── DAILY LOOP (365 days) ────────────────────────────────────
        for day in range(DAYS_PER_YEAR):
            month = min(day // 30, 11)

            # 1. Update environment at each node
            sim_day = year * DAYS_PER_YEAR + day
            for i, node in enumerate(network.nodes):
                nd = node.definition
                if use_precomputed_sst and satellite_sst_by_node[i] is not None:
                    node.current_sst = float(satellite_sst_by_node[i][sim_day])
                else:
                    node.current_sst = sst_with_trend(
                        day, cal_year, nd.mean_sst, nd.sst_amplitude,
                        nd.sst_trend, reference_year=start_year,
                    )
                node.current_flushing = seasonal_flushing(
                    nd.flushing_rate, month, nd.is_fjord,
                )
                node.current_salinity = nd.salinity

            # 2a. Movement (CRW sub-steps within the day)
            if mov_enabled:
                for i, node in enumerate(network.nodes):
                    hab_side = np.sqrt(max(node.definition.habitat_area, 1.0))
                    daily_movement(
                        agents=node.agents,
                        habitat_side=hab_side,
                        sigma_turn=mov_cfg.sigma_turn,
                        base_speed=mov_cfg.base_speed,
                        substeps=mov_cfg.substeps_per_day,
                        rng=rng,
                    )

            # 2b. Build spatial transmission grids (if enabled)
            if spatial_tx:
                for i, node in enumerate(network.nodes):
                    if density_grids[i] is not None:
                        density_grids[i].build(
                            node.agents,
                            diffusion_passes=mov_cfg.diffusion_passes,
                        )

            # 3. Per-node disease step (if disease is active at that node)
            for i, node in enumerate(network.nodes):
                if not disease_active_flags[i]:
                    continue
                deaths_before = node_disease_states[i].cumulative_deaths
                node_disease_states[i] = daily_disease_update(
                    agents=node.agents,
                    node_state=node_disease_states[i],
                    T_celsius=node.current_sst,
                    salinity=node.current_salinity,
                    phi_k=node.current_flushing,
                    dispersal_input=0.0,  # handled below
                    day=year * DAYS_PER_YEAR + day,
                    cfg=dis_cfg,
                    rng=rng,
                    infected_density_grid=(
                        density_grids[i] if spatial_tx else None
                    ),
                    pe_cfg=pe_cfg,
                )
                new_deaths = (
                    node_disease_states[i].cumulative_deaths - deaths_before
                )
                cumulative_dis_deaths[i] += new_deaths
                cumulative_recoveries[i] = (
                    node_disease_states[i].cumulative_recoveries
                )

            # 3. Pathogen dispersal between nodes (D matrix)
            P = np.array([
                node_disease_states[i].vibrio_concentration
                for i in range(N)
            ])
            dispersal_in = pathogen_dispersal_step(P, network.D)
            for i in range(N):
                node_disease_states[i].vibrio_concentration += dispersal_in[i]
                network.nodes[i].vibrio_concentration = (
                    node_disease_states[i].vibrio_concentration
                )

            # 4. Continuous settlement: settle cohorts whose PLD has elapsed
            sim_day = year * DAYS_PER_YEAR + day
            if spawning_enabled:
                for i, node in enumerate(network.nodes):
                    if not pending_cohorts[i]:
                        continue
                    # Phase 7 optimization: pop from sorted front instead
                    # of filtering the entire list each day
                    ready = _pop_ready_cohorts(pending_cohorts[i], sim_day)
                    if ready:
                        nd = node.definition
                        n_settled = settle_daily_cohorts(
                            ready, node.agents, node.genotypes,
                            nd.carrying_capacity, pop_cfg, rng,
                            effect_sizes, node_id=nd.node_id,
                            habitat_area=nd.habitat_area,
                            sim_day=sim_day,
                            effects_t=effects_t,
                            effects_c=effects_c,
                            res_slice=res_slice,
                            tol_slice=tol_slice,
                            rec_slice=rec_slice,
                        )
                        yearly_recruits_accum[i] += n_settled

            # 5. Daily spawning step (if enabled)
            if spawning_enabled:
                day_of_year = (day + 1)  # Convert 0-based to 1-based
                for i, node in enumerate(network.nodes):
                    nd = node.definition
                    currently_in_season = in_spawning_season(
                        day_of_year,
                        config.spawning.season_start_doy,
                        config.spawning.season_end_doy,
                    )
                    # Reset spawning season when transitioning out
                    if previous_in_season[i] and not currently_in_season:
                        reset_spawning_season(node.agents)
                    # Daily spawning during season
                    if currently_in_season:
                        cohorts_today = spawning_step(
                            agents=node.agents,
                            genotypes=node.genotypes,
                            day_of_year=day_of_year,
                            node_latitude=nd.lat,
                            spawning_config=config.spawning,
                            disease_config=config.disease,
                            rng=rng,
                            current_sim_day=sim_day,
                            current_sst=node.current_sst,
                        )
                        if cohorts_today:
                            # Tag cohorts with source node and insert sorted
                            for c in cohorts_today:
                                c.source_node = nd.node_id
                                _insort_cohort(pending_cohorts[i], c)
                        # Count spawners today at this node
                        alive_spawned = (
                            node.agents['alive'] &
                            (node.agents['last_spawn_day'] == day_of_year)
                        )
                        daily_spawning_counts_spatial[i, sim_day] = int(np.sum(alive_spawned))
                    # Tick down immunosuppression timers
                    immuno_mask = node.agents['immunosuppression_timer'] > 0
                    node.agents['immunosuppression_timer'][immuno_mask] -= 1
                    previous_in_season[i] = currently_in_season

            # 6. Daily demographics (continuous mortality + growth)
            #    Compute alive_idx ONCE per node, share across both calls
            for i, node in enumerate(network.nodes):
                nd = node.agents
                _alive_idx = np.where(nd['alive'])[0]
                n_mort, n_senes = daily_natural_mortality(nd, pop_cfg, rng, alive_idx=_alive_idx)
                node_daily_nat_deaths[i] += n_mort
                node_daily_senes_deaths[i] += n_senes
                # Recompute alive_idx after mortality (some agents may have died)
                if n_mort > 0:
                    _alive_idx = np.where(nd['alive'])[0]
                daily_growth_and_aging(nd, pop_cfg, rng, alive_idx=_alive_idx)

            # Optional individual-level snapshot recording
            if snapshot_recorder is not None:
                snapshot_recorder.capture_all_nodes(
                    sim_day, year, network.nodes
                )

            # Track max vibrio and disease prevalence
            for i in range(N):
                v = node_disease_states[i].vibrio_concentration
                if v > max_vibrio_year[i]:
                    max_vibrio_year[i] = v
                # Track peak disease prevalence
                if disease_active_flags[i]:
                    node = network.nodes[i]
                    alive_m = node.agents['alive'].astype(bool)
                    n_alive = int(np.sum(alive_m))
                    if n_alive > 0:
                        ds = node.agents['disease_state']
                        n_inf = int(np.sum(
                            alive_m & (
                                (ds == DiseaseState.E) |
                                (ds == DiseaseState.I1) |
                                (ds == DiseaseState.I2)
                            )
                        ))
                        prev = n_inf / n_alive
                        if prev > peak_disease_prev[i]:
                            peak_disease_prev[i] = prev

        # ── ANNUAL DEMOGRAPHIC STEP ──────────────────────────────────

        # 1. Natural mortality + growth now handled daily (Phase A above)
        node_nat_deaths = node_daily_nat_deaths.copy()
        node_senes_deaths = node_daily_senes_deaths.copy()

        # 2. Larval dispersal via C matrix (annual)
        #    Collect unsettled cohorts from each source node, disperse via C,
        #    place into dest pending_cohorts. Cohorts with PLD elapsed settle
        #    immediately; others settle in next year's daily loop.
        source_ids = []
        source_n = []
        source_geno = []
        source_cohort_meta = []  # (spawn_day, pld_days, sst_at_spawn) per source batch

        for i in range(N):
            # Merge all pending cohorts at this node into one batch for C matrix
            node_cohorts = pending_cohorts[i]
            if not node_cohorts:
                continue
            all_geno = []
            all_spawn_days = []
            all_pld_days = []
            all_sst = []
            for c in node_cohorts:
                if c.n_competent > 0 and c.genotypes is not None:
                    for j in range(min(c.n_competent, len(c.genotypes))):
                        all_geno.append(c.genotypes[j])
                        all_spawn_days.append(c.spawn_day)
                        all_pld_days.append(c.pld_days)
                        all_sst.append(c.sst_at_spawn)
            pending_cohorts[i] = []  # clear source pending

            if not all_geno:
                continue
            n_larvae = len(all_geno)
            # Cap per node to avoid memory issues
            if n_larvae > 5000:
                idx = rng.choice(n_larvae, size=5000, replace=False)
                all_geno = [all_geno[j] for j in idx]
                all_spawn_days = [all_spawn_days[j] for j in idx]
                all_pld_days = [all_pld_days[j] for j in idx]
                all_sst = [all_sst[j] for j in idx]
                n_larvae = 5000
            source_ids.append(i)
            source_n.append(n_larvae)
            source_geno.append(np.array(all_geno))
            source_cohort_meta.append(
                (all_spawn_days, all_pld_days, all_sst)
            )

        total_larvae = sum(source_n)
        if total_larvae > 0:
            settler_map = distribute_larvae(
                source_ids, source_n, source_geno,
                network.C, rng,
            )

            # Build a lookup: source_idx -> meta arrays for reassignment
            source_meta_by_id = {}
            for idx, sid in enumerate(source_ids):
                source_meta_by_id[sid] = source_cohort_meta[idx]

            # Place dispersed cohorts into dest pending_cohorts
            end_of_year_day = (year + 1) * DAYS_PER_YEAR
            for k in range(N):
                if not settler_map[k]:
                    continue
                for geno_batch, src in settler_map[k]:
                    n_batch = len(geno_batch)
                    if n_batch == 0:
                        continue
                    # Get meta from source (spawn_day, pld_days, sst)
                    meta = source_meta_by_id.get(src)
                    if meta is None:
                        # Fallback: use default PLD from dest node SST
                        fallback_pld = pelagic_larval_duration(
                            network.nodes[k].current_sst
                        )
                        spawn_days_batch = [end_of_year_day - int(fallback_pld)] * n_batch
                        pld_days_batch = [fallback_pld] * n_batch
                        sst_batch = [network.nodes[k].current_sst] * n_batch
                    else:
                        all_sd, all_pld, all_sst_m = meta
                        n_meta = len(all_sd)
                        # Sample meta for each settler (dispersed larvae
                        # are random draws from source genotypes)
                        meta_idx = rng.choice(n_meta, size=n_batch, replace=True)
                        spawn_days_batch = [all_sd[j] for j in meta_idx]
                        pld_days_batch = [all_pld[j] for j in meta_idx]
                        sst_batch = [all_sst_m[j] for j in meta_idx]

                    # Create individual LarvalCohort objects per settler
                    # Group by (spawn_day, pld_days) for efficiency
                    from collections import defaultdict
                    groups = defaultdict(list)
                    for j in range(n_batch):
                        key = (spawn_days_batch[j], pld_days_batch[j])
                        groups[key].append(j)

                    for (sd, pld), indices in groups.items():
                        cohort_geno = geno_batch[indices]
                        cohort = LarvalCohort(
                            source_node=src,
                            n_competent=len(indices),
                            genotypes=cohort_geno,
                            parent_pairs=np.zeros((len(indices), 2), dtype=np.int32),
                            pld_days=pld,
                            spawn_day=sd,
                            sst_at_spawn=sst_batch[indices[0]],
                        )
                        # If PLD already elapsed, settle immediately
                        if (end_of_year_day - sd) >= pld:
                            nd = network.nodes[k].definition
                            n_settled = settle_daily_cohorts(
                                [cohort], network.nodes[k].agents,
                                network.nodes[k].genotypes,
                                nd.carrying_capacity, pop_cfg, rng,
                                effect_sizes, node_id=nd.node_id,
                                habitat_area=nd.habitat_area,
                                sim_day=end_of_year_day,
                                effects_t=effects_t,
                                effects_c=effects_c,
                                res_slice=res_slice,
                                tol_slice=tol_slice,
                                rec_slice=rec_slice,
                            )
                            yearly_recruits_accum[k] += n_settled
                        else:
                            _insort_cohort(pending_cohorts[k], cohort)

        # 5. Pre-epidemic allele frequency snapshot
        if (disease_year is not None and year == disease_year
                and not pre_snapshot_taken):
            for i, node in enumerate(network.nodes):
                alive_mask = node.agents['alive'].astype(bool)
                if int(np.sum(alive_mask)) > 0:
                    pre_epidemic_af[i] = compute_allele_frequencies(
                        node.genotypes, alive_mask
                    )
            pre_snapshot_taken = True

        # 6. Seed disease if it's the epidemic year
        if disease_year is not None and year == disease_year:
            epidemic_started_year = year
            from sswd_evoepi.disease import sample_stage_duration, K_SHAPE_E
            for i, node in enumerate(network.nodes):
                disease_active_flags[i] = True
                alive_idx = np.where(node.agents['alive'])[0]
                n_to_infect = min(initial_infected_per_node, len(alive_idx))
                if n_to_infect > 0:
                    infect_idx = rng.choice(alive_idx, size=n_to_infect,
                                            replace=False)
                    mu_EI1 = arrhenius(dis_cfg.mu_EI1_ref, dis_cfg.Ea_EI1,
                                       node.current_sst)
                    for idx in infect_idx:
                        node.agents['disease_state'][idx] = DiseaseState.E
                        node.agents['disease_timer'][idx] = (
                            sample_stage_duration(mu_EI1, K_SHAPE_E, rng)
                        )
                        if pe_cfg is not None:
                            node.agents['pathogen_virulence'][idx] = pe_cfg.v_init

        # 7. Post-epidemic allele frequency snapshot (2 years after)
        if (epidemic_started_year is not None
                and year == epidemic_started_year + 2
                and not post_snapshot_taken):
            for i, node in enumerate(network.nodes):
                alive_mask = node.agents['alive'].astype(bool)
                if int(np.sum(alive_mask)) > 1:
                    post_epidemic_af[i] = compute_allele_frequencies(
                        node.genotypes, alive_mask
                    )
            post_snapshot_taken = True

        # ── Record annual metrics ────────────────────────────────────
        for i, node in enumerate(network.nodes):
            alive_mask = node.agents['alive'].astype(bool)
            pop_now = int(np.sum(alive_mask))
            adult_mask = alive_mask & (node.agents['stage'] == Stage.ADULT)
            yearly_pop[i, year] = pop_now
            yearly_adults[i, year] = int(np.sum(adult_mask))
            yearly_recruits[i, year] = yearly_recruits_accum[i]
            yearly_nat_deaths[i, year] = node_nat_deaths[i]
            yearly_senes_deaths[i, year] = node_senes_deaths[i]
            yearly_dis_deaths[i, year] = (
                cumulative_dis_deaths[i] - year_dd_before[i]
            )
            yearly_recoveries[i, year] = (
                cumulative_recoveries[i] - year_rec_before[i]
            )
            yearly_vibrio_max[i, year] = max_vibrio_year[i]
            if pop_now > 0:
                yearly_mean_r[i, year] = float(
                    np.mean(node.agents['resistance'][alive_mask])
                )
                yearly_mean_t[i, year] = float(
                    np.mean(node.agents['tolerance'][alive_mask])
                )
                yearly_mean_c[i, year] = float(
                    np.mean(node.agents['recovery_ability'][alive_mask])
                )
                # Genetics tracking
                af = compute_allele_frequencies(node.genotypes, alive_mask)
                yearly_allele_freq_top3[i, year] = af[:3]
                yearly_va[i, year] = compute_additive_variance(af, effect_sizes, res_slice)
                yearly_va_t[i, year] = compute_additive_variance(af, effects_t, tol_slice)
                yearly_va_c[i, year] = compute_additive_variance(af, effects_c, rec_slice)

            # Per-node virulence tracking (time-weighted daily mean)
            if pe_cfg is not None:
                nds = node_disease_states[i]
                if nds.virulence_count_daily > 0:
                    yearly_mean_v_spatial[i, year] = (
                        nds.virulence_sum_daily / nds.virulence_count_daily
                    )
                # Reset accumulators for next year
                nds.virulence_sum_new_infections = 0.0
                nds.virulence_count_new_infections = 0
                nds.virulence_sum_deaths = 0.0
                nds.virulence_count_deaths = 0
                nds.virulence_sum_daily = 0.0
                nds.virulence_count_daily = 0

        yearly_total_pop[year] = int(yearly_pop[:, year].sum())
        yearly_total_larvae[year] = total_larvae if total_larvae > 0 else 0

        if progress_callback is not None:
            progress_callback(year, n_years)

    # ── Compile result ────────────────────────────────────────────────
    return SpatialSimResult(
        n_years=n_years,
        n_nodes=N,
        node_names=node_names,
        node_K=node_K,
        yearly_pop=yearly_pop,
        yearly_adults=yearly_adults,
        yearly_recruits=yearly_recruits,
        yearly_natural_deaths=yearly_nat_deaths,
        yearly_disease_deaths=yearly_dis_deaths,
        yearly_recoveries=yearly_recoveries,
        yearly_senescence_deaths=yearly_senes_deaths,
        yearly_mean_resistance=yearly_mean_r,
        yearly_mean_tolerance=yearly_mean_t,
        yearly_mean_recovery=yearly_mean_c,
        yearly_vibrio_max=yearly_vibrio_max,
        yearly_va=yearly_va,
        yearly_va_tolerance=yearly_va_t,
        yearly_va_recovery=yearly_va_c,
        yearly_allele_freq_top3=yearly_allele_freq_top3,
        yearly_ne_ratio=yearly_ne_ratio,
        pre_epidemic_allele_freq=pre_epidemic_af if pre_snapshot_taken else None,
        post_epidemic_allele_freq=post_epidemic_af if post_snapshot_taken else None,
        yearly_total_pop=yearly_total_pop,
        yearly_total_larvae_dispersed=yearly_total_larvae,
        peak_disease_prevalence=peak_disease_prev,
        yearly_mean_virulence=yearly_mean_v_spatial,
        daily_spawning_counts=daily_spawning_counts_spatial,
        initial_total_pop=initial_total,
        final_total_pop=int(yearly_total_pop[-1]),
        disease_year=disease_year,
        seed=seed,
        snapshot_recorder=snapshot_recorder,
    )
