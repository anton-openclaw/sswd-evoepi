"""Broadcast spawning, SRS lottery, fertilization kinetics, and settlement.

Implements Phase 3 of the population dynamics module:
  - Spawning phenology (spring/summer, latitude-dependent)
  - Size-dependent fecundity (allometric)
  - Quadratic Allee fertilization: F(D) per Lundquist & Botsford 2004
  - Sweepstakes reproductive success: Pareto(α) offspring distribution
  - Mendelian inheritance with independent assortment
  - EF1A lethal homozygote elimination
  - Ne/N effective population size tracking
  - Larval production and pelagic survival
  - Settlement with density-dependent recruitment (Beverton-Holt)
  - Settlement cue Allee effect (adult biofilm)

References:
  - population-dynamics-spec.md §3 (reproduction), §4 (larval), §5–6 (density/Allee)
  - genetics-evolution-spec.md §4 (SRS), §5 (Mendelian inheritance)
  - CODE_ERRATA CE-1: cost_resistance removed; fecundity_mod always 1.0

Authors: Anton 🔬 & Willem Weertman
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Tuple

import numpy as np

from sswd_evoepi.types import (
    AGENT_DTYPE,
    IDX_EF1A,
    N_ADDITIVE,
    N_LOCI,
    DiseaseState,
    LarvalCohort,
    Stage,
    allocate_agents,
    allocate_genotypes,
)


# ═══════════════════════════════════════════════════════════════════════
# SPAWNING PHENOLOGY
# ═══════════════════════════════════════════════════════════════════════

# Base spawning day-of-year at reference latitude (48°N, San Juan Islands)
_SPAWNING_DAY_BASE = 150     # ~May 30
_SPAWNING_LAT_REF = 48.0     # °N
_SPAWNING_LAT_SHIFT = 3.0    # days per degree latitude
_SPAWNING_DAY_MIN = 120      # earliest (late April)
_SPAWNING_DAY_MAX = 200      # latest (mid-July)


def get_spawning_day(lat: float) -> int:
    """Day-of-year for spawning at a given latitude.

    Base: day 150 (May 30) at 48°N (San Juan Islands).
    Shift: +3 days per degree north, −3 days per degree south.

    Args:
        lat: Latitude in degrees north.

    Returns:
        Day of year (1–365) for spawning pulse.
    """
    shift = int(_SPAWNING_LAT_SHIFT * (lat - _SPAWNING_LAT_REF))
    return max(_SPAWNING_DAY_MIN, min(_SPAWNING_DAY_MAX, _SPAWNING_DAY_BASE + shift))


def is_spawning_season(day_of_year: int, spawning_day: int, window: int = 0) -> bool:
    """Check if current day falls within the spawning window.

    For the annual pulse model, spawning occurs on exactly one day.
    The window parameter allows a ±window-day tolerance for testing.

    Args:
        day_of_year: Current day of year (1–365).
        spawning_day: Target spawning day of year.
        window: Half-width of spawning window in days (default 0 = pulse).

    Returns:
        True if day_of_year is within [spawning_day − window, spawning_day + window].
    """
    return abs(day_of_year - spawning_day) <= window


# ═══════════════════════════════════════════════════════════════════════
# FECUNDITY
# ═══════════════════════════════════════════════════════════════════════


def fecundity(
    size_mm: float,
    F0: float = 1.0e7,
    L_ref: float = 500.0,
    fecundity_exp: float = 2.5,
    L_min_repro: float = 400.0,
) -> float:
    """Eggs per female per annual spawning event.

    Allometric scaling: F(L) = F0 × (L / L_ref)^b

    CE-1: cost_of_resistance removed — fecundity depends only on size.

    Args:
        size_mm: Individual arm-tip-to-arm-tip diameter (mm).
        F0: Reference fecundity at L_ref (eggs).
        L_ref: Reference body size (mm).
        fecundity_exp: Allometric exponent b.
        L_min_repro: Minimum reproductive size (mm).

    Returns:
        Number of eggs (float, ≥ 0).
    """
    if size_mm < L_min_repro:
        return 0.0
    return F0 * (size_mm / L_ref) ** fecundity_exp


def total_eggs(
    agents: np.ndarray,
    female_idx: np.ndarray,
    F0: float = 1.0e7,
    L_ref: float = 500.0,
    fecundity_exp: float = 2.5,
    L_min_repro: float = 400.0,
) -> float:
    """Compute total eggs from all spawning females (vectorized).

    Args:
        agents: Structured array with AGENT_DTYPE.
        female_idx: Indices of spawning females in agents.
        F0, L_ref, fecundity_exp, L_min_repro: Fecundity parameters.

    Returns:
        Total number of eggs (float).
    """
    if len(female_idx) == 0:
        return 0.0
    sizes = agents['size'][female_idx]
    fec = np.where(
        sizes >= L_min_repro,
        F0 * (sizes / L_ref) ** fecundity_exp,
        0.0,
    )
    return float(fec.sum())


# ═══════════════════════════════════════════════════════════════════════
# FERTILIZATION KINETICS — ALLEE EFFECT
# ═══════════════════════════════════════════════════════════════════════


def fertilization_success(
    male_density: float,
    gamma_fert: float = 4.5,
    aggregation_factor: float = 1.0,
) -> float:
    """Fraction of eggs fertilized at given male adult density.

    F(ρ) = 1 − exp(−γ × ρ_eff)

    This is the mean-field approximation of the Lundquist & Botsford (2004)
    broadcast spawner fertilization model. Produces a quadratic relationship
    between zygote production and density at low density (the Allee effect).

    Args:
        male_density: Male adult density (ind/m²).
        gamma_fert: Sperm contact parameter (m²). Default 4.5.
        aggregation_factor: Multiplier ≥ 1.0 for spawning aggregation.

    Returns:
        Fertilization success F ∈ [0, 1].
    """
    rho_eff = male_density * aggregation_factor
    return 1.0 - np.exp(-gamma_fert * rho_eff)


def compute_aggregation_factor(
    n_adults: int,
    habitat_area: float,
    n_agg_threshold: int = 20,
    rho_max_agg: float = 0.5,
) -> float:
    """Compute factor by which aggregation increases effective spawning density.

    If adults aggregate into clumps during spawning season, effective local
    density within the clump is higher than the spatially uniform average.

    Args:
        n_adults: Total number of adults at the node.
        habitat_area: Habitat area in m².
        n_agg_threshold: Adults needed for full aggregation behavior.
        rho_max_agg: Effective density within aggregation clump (ind/m²).

    Returns:
        Multiplier ≥ 1.0.
    """
    if n_adults < 2 or habitat_area <= 0:
        return 1.0

    ambient_density = n_adults / habitat_area
    if ambient_density >= rho_max_agg:
        return 1.0  # already dense; no benefit

    # Fraction of population that successfully aggregates
    p_agg = min(1.0, n_adults / n_agg_threshold)

    # Blended effective density
    rho_eff = p_agg * rho_max_agg + (1.0 - p_agg) * ambient_density
    return rho_eff / max(ambient_density, 1e-10)


# ═══════════════════════════════════════════════════════════════════════
# SWEEPSTAKES REPRODUCTIVE SUCCESS (SRS)
# ═══════════════════════════════════════════════════════════════════════


def srs_reproductive_lottery(
    females: np.ndarray,
    males: np.ndarray,
    agents: np.ndarray,
    genotypes: np.ndarray,
    n_offspring_target: int,
    alpha_srs: float = 1.35,
    F0: float = 1.0e7,
    L_ref: float = 500.0,
    fecundity_exp: float = 2.5,
    L_min_repro: float = 400.0,
    rng: Optional[np.random.Generator] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Sweepstakes reproductive lottery with Mendelian inheritance.

    Each spawning parent receives a random Pareto(α) weight; females are
    further weighted by fecundity (size-dependent). Parent pairs are sampled
    with replacement, offspring receive Mendelian-inherited genotypes, and
    EF1A lethal homozygotes are eliminated.

    CE-1: No cost_of_resistance in quality weighting.

    Args:
        females: Indices of spawning females in agents/genotypes arrays.
        males: Indices of spawning males.
        agents: Agent structured array (AGENT_DTYPE).
        genotypes: Genotype array (max_agents, N_LOCI, 2) int8.
        n_offspring_target: Number of viable offspring to produce.
        alpha_srs: Pareto shape parameter (lower = more skew).
        F0, L_ref, fecundity_exp, L_min_repro: Fecundity parameters.
        rng: NumPy random Generator.

    Returns:
        offspring_genotypes: (n_valid, N_LOCI, 2) int8 array.
        parent_pairs: (n_valid, 2) int32 array — (mother_idx, father_idx).
    """
    if rng is None:
        rng = np.random.default_rng()

    n_f = len(females)
    n_m = len(males)

    empty_geno = np.empty((0, N_LOCI, 2), dtype=np.int8)
    empty_pairs = np.empty((0, 2), dtype=np.int32)

    if n_f == 0 or n_m == 0 or n_offspring_target <= 0:
        return empty_geno, empty_pairs

    # --- Step 1: Pareto weights ---
    raw_f = rng.pareto(alpha_srs, size=n_f) + 1.0
    raw_m = rng.pareto(alpha_srs, size=n_m) + 1.0

    # --- Step 2: Female quality = fecundity(size) ---
    # CE-1: no cost-of-resistance penalty
    sizes = agents['size'][females].astype(np.float64)
    quality_f = np.where(
        sizes >= L_min_repro,
        (sizes / L_ref) ** fecundity_exp,
        1e-6,  # near-zero for undersized (shouldn't happen if filtered)
    )
    weighted_f = raw_f * quality_f

    # --- Step 3: Normalize to probability distributions ---
    sum_f = weighted_f.sum()
    sum_m = raw_m.sum()
    if sum_f < 1e-12 or sum_m < 1e-12:
        return empty_geno, empty_pairs

    prob_f = weighted_f / sum_f
    prob_m = raw_m / sum_m

    # --- Step 4: Oversample to account for EF1A lethal elimination (~6%) ---
    n_sample = int(n_offspring_target * 1.10) + 20
    mother_idx = rng.choice(females, size=n_sample, p=prob_f)
    father_idx = rng.choice(males, size=n_sample, p=prob_m)

    # --- Step 5: Mendelian inheritance ---
    offspring_geno = mendelian_inherit_batch(
        genotypes, mother_idx, father_idx, n_sample, rng
    )

    # --- Step 6: EF1A lethal elimination ---
    ef1a_sum = offspring_geno[:, IDX_EF1A, :].sum(axis=1)
    viable_mask = ef1a_sum < 2  # eliminate ins/ins (sum == 2)
    offspring_geno = offspring_geno[viable_mask]
    parent_pairs = np.column_stack(
        [mother_idx[viable_mask], father_idx[viable_mask]]
    ).astype(np.int32)

    # Trim to target
    if len(offspring_geno) > n_offspring_target:
        offspring_geno = offspring_geno[:n_offspring_target]
        parent_pairs = parent_pairs[:n_offspring_target]

    return offspring_geno, parent_pairs


# ═══════════════════════════════════════════════════════════════════════
# MENDELIAN INHERITANCE
# ═══════════════════════════════════════════════════════════════════════


def mendelian_inherit_batch(
    genotypes: np.ndarray,
    mother_idx: np.ndarray,
    father_idx: np.ndarray,
    n_offspring: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Vectorized Mendelian inheritance for a batch of offspring.

    At each locus, offspring receives one randomly chosen allele from each
    parent. Independent assortment (no linkage) is assumed.

    Args:
        genotypes: (max_agents, N_LOCI, 2) int8 parental genotypes.
        mother_idx: (n_offspring,) indices of mothers.
        father_idx: (n_offspring,) indices of fathers.
        n_offspring: Number of offspring.
        rng: Random generator.

    Returns:
        offspring_geno: (n_offspring, N_LOCI, 2) int8.
    """
    # Pre-draw all allele choices: 0 or 1 for maternal and paternal
    choices = rng.integers(0, 2, size=(n_offspring, N_LOCI, 2), dtype=np.int8)

    offspring = np.empty((n_offspring, N_LOCI, 2), dtype=np.int8)
    idx_range = np.arange(n_offspring)

    for l in range(N_LOCI):
        # Maternal alleles at locus l
        mat_geno = genotypes[mother_idx, l, :]   # (n_offspring, 2)
        mat_choice = choices[:, l, 0]             # (n_offspring,)
        offspring[:, l, 0] = mat_geno[idx_range, mat_choice]

        # Paternal alleles at locus l
        pat_geno = genotypes[father_idx, l, :]
        pat_choice = choices[:, l, 1]
        offspring[:, l, 1] = pat_geno[idx_range, pat_choice]

    return offspring


# ═══════════════════════════════════════════════════════════════════════
# Ne/N EFFECTIVE POPULATION SIZE
# ═══════════════════════════════════════════════════════════════════════


def compute_ne(offspring_counts: np.ndarray) -> float:
    """Compute effective population size from offspring count distribution.

    Ne = (4N − 2) / (Vk + 2)

    Where N = number of parents, Vk = variance in offspring number.
    Hedgecock & Pudovkin 2011.

    Args:
        offspring_counts: (N_parents,) array of offspring per parent.

    Returns:
        Effective population size Ne (float).
    """
    N = len(offspring_counts)
    if N < 2:
        return float(N)
    Vk = float(np.var(offspring_counts, ddof=1))
    return (4.0 * N - 2.0) / (Vk + 2.0)


def compute_ne_ratio(parent_pairs: np.ndarray, n_parents: int) -> float:
    """Compute Ne/N from parent-pair records.

    Args:
        parent_pairs: (n_offspring, 2) array of (mother_idx, father_idx).
        n_parents: Total number of potential parents (females + males).

    Returns:
        Ne/N ratio. Returns 1.0 if n_parents < 2.
    """
    if n_parents < 2 or len(parent_pairs) == 0:
        return 1.0

    # Count offspring per unique parent
    all_parent_ids = parent_pairs.ravel()
    unique_parents = np.unique(all_parent_ids)
    offspring_per_parent = np.zeros(n_parents, dtype=np.int64)
    for pid in all_parent_ids:
        offspring_per_parent[pid % n_parents] += 1

    Ne = compute_ne(offspring_per_parent)
    return Ne / n_parents


def compute_ne_from_pairs(
    parent_pairs: np.ndarray,
    females: np.ndarray,
    males: np.ndarray,
) -> Tuple[float, float]:
    """Compute Ne separately for females and males, return harmonic mean Ne and Ne/N.

    Args:
        parent_pairs: (n_offspring, 2) — col 0 = mother idx, col 1 = father idx.
        females: Array of female parent indices.
        males: Array of male parent indices.

    Returns:
        (Ne, Ne_over_N) tuple.
    """
    n_f = len(females)
    n_m = len(males)
    N_total = n_f + n_m

    if N_total < 2 or len(parent_pairs) == 0:
        return float(N_total), 1.0

    # Offspring counts per female
    female_counts = np.zeros(n_f, dtype=np.int64)
    for i, f in enumerate(females):
        female_counts[i] = np.sum(parent_pairs[:, 0] == f)

    # Offspring counts per male
    male_counts = np.zeros(n_m, dtype=np.int64)
    for i, m in enumerate(males):
        male_counts[i] = np.sum(parent_pairs[:, 1] == m)

    Ne_f = compute_ne(female_counts) if n_f >= 2 else float(n_f)
    Ne_m = compute_ne(male_counts) if n_m >= 2 else float(n_m)

    # Harmonic mean: Ne = 4 * Ne_f * Ne_m / (Ne_f + Ne_m)
    if Ne_f + Ne_m > 0:
        Ne = 4.0 * Ne_f * Ne_m / (Ne_f + Ne_m)
    else:
        Ne = 0.0

    Ne_over_N = Ne / N_total if N_total > 0 else 0.0
    return Ne, Ne_over_N


# ═══════════════════════════════════════════════════════════════════════
# LARVAL PRODUCTION & PELAGIC SURVIVAL
# ═══════════════════════════════════════════════════════════════════════

# Temperature-dependent PLD parameters
_PLD_REF = 63.0        # days at reference temperature
_T_REF_DEV = 10.5      # °C (Hodin 2021 lab conditions)
_Q_DEV = 0.05          # temperature coefficient (°C⁻¹)
_MU_LARVA_DAILY = 0.05 # daily larval mortality rate


def pelagic_larval_duration(sst: float) -> float:
    """Temperature-dependent pelagic larval duration (days).

    PLD(T) = PLD_ref × exp(−Q_dev × (T − T_ref))

    Warmer → shorter PLD. Clamped to [30, 150] days.

    Args:
        sst: Sea surface temperature (°C).

    Returns:
        PLD in days.
    """
    pld = _PLD_REF * np.exp(-_Q_DEV * (sst - _T_REF_DEV))
    return float(np.clip(pld, 30.0, 150.0))


def larval_survival(pld_days: float, mu_daily: float = _MU_LARVA_DAILY) -> float:
    """Fraction of larvae surviving the pelagic phase.

    survival = exp(−μ × PLD)

    Args:
        pld_days: Pelagic larval duration (days).
        mu_daily: Daily mortality rate (d⁻¹).

    Returns:
        Survival fraction ∈ (0, 1).
    """
    return float(np.exp(-mu_daily * pld_days))


def produce_larval_cohort(
    node_id: int,
    agents: np.ndarray,
    genotypes: np.ndarray,
    habitat_area: float,
    sst: float,
    lat: float,
    alpha_srs: float = 1.35,
    gamma_fert: float = 4.5,
    F0: float = 1.0e7,
    L_ref: float = 500.0,
    fecundity_exp: float = 2.5,
    L_min_repro: float = 400.0,
    rng: Optional[np.random.Generator] = None,
) -> Tuple[Optional[LarvalCohort], dict]:
    """Full spawning pipeline: identify spawners → eggs → fertilization → larvae → SRS.

    Returns a LarvalCohort for spatial dispersal and a diagnostics dict.

    Args:
        node_id: Node index.
        agents: Agent structured array.
        genotypes: Genotype array.
        habitat_area: Node habitat area (m²).
        sst: SST at spawning time (°C).
        lat: Node latitude (°N).
        alpha_srs: Pareto SRS shape.
        gamma_fert: Fertilization kinetics parameter.
        F0, L_ref, fecundity_exp, L_min_repro: Fecundity params.
        rng: Random generator.

    Returns:
        (cohort, diagnostics): LarvalCohort or None, plus dict of stats.
    """
    if rng is None:
        rng = np.random.default_rng()

    alive = agents['alive']
    stage = agents['stage']
    sex = agents['sex']
    disease = agents['disease_state']

    # Spawning requires: alive, ADULT, healthy (S or R)
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
        'total_eggs': 0.0,
        'fertilization_success': 0.0,
        'n_zygotes': 0,
        'pld_days': 0.0,
        'pelagic_survival': 0.0,
        'n_competent': 0,
        'ne': 0.0,
        'ne_over_n': 0.0,
    }

    if len(females) == 0 or len(males) == 0:
        return None, diag

    # Total eggs
    eggs = total_eggs(agents, females, F0, L_ref, fecundity_exp, L_min_repro)
    diag['total_eggs'] = eggs

    if eggs <= 0:
        return None, diag

    # Fertilization (Allee effect)
    n_adults = len(females) + len(males)
    male_density = len(males) / habitat_area if habitat_area > 0 else 0.0
    agg_factor = compute_aggregation_factor(n_adults, habitat_area)
    fert = fertilization_success(male_density, gamma_fert, agg_factor)
    diag['fertilization_success'] = fert

    n_zygotes = int(eggs * fert)
    diag['n_zygotes'] = n_zygotes

    if n_zygotes <= 0:
        return None, diag

    # Pelagic larval phase
    pld = pelagic_larval_duration(sst)
    surv = larval_survival(pld)
    diag['pld_days'] = pld
    diag['pelagic_survival'] = surv

    n_competent = max(1, int(n_zygotes * surv))
    diag['n_competent'] = n_competent

    # Cap competent larvae to avoid memory explosion — reasonable biological max
    n_competent = min(n_competent, 100_000)

    # SRS lottery: sample parents and produce offspring genotypes
    offspring_geno, parent_pairs = srs_reproductive_lottery(
        females=females,
        males=males,
        agents=agents,
        genotypes=genotypes,
        n_offspring_target=n_competent,
        alpha_srs=alpha_srs,
        F0=F0,
        L_ref=L_ref,
        fecundity_exp=fecundity_exp,
        L_min_repro=L_min_repro,
        rng=rng,
    )

    n_actual = len(offspring_geno)
    diag['n_competent'] = n_actual

    if n_actual == 0:
        return None, diag

    # Ne/N diagnostics
    Ne, Ne_over_N = compute_ne_from_pairs(parent_pairs, females, males)
    diag['ne'] = Ne
    diag['ne_over_n'] = Ne_over_N

    cohort = LarvalCohort(
        source_node=node_id,
        n_competent=n_actual,
        genotypes=offspring_geno,
        parent_pairs=parent_pairs,
        pld_days=pld,
    )

    return cohort, diag


# ═══════════════════════════════════════════════════════════════════════
# SETTLEMENT & RECRUITMENT
# ═══════════════════════════════════════════════════════════════════════


def settlement_cue_modifier(n_adults: int, halfsat: int = 5) -> float:
    """Modifier on settlement success based on adult presence (Allee effect).

    Michaelis-Menten form: baseline 0.2 (coralline algae) to 1.0 (full adult biofilm).

    Args:
        n_adults: Number of adults present at receiving node.
        halfsat: Half-saturation adult count.

    Returns:
        Modifier ∈ [0.2, 1.0].
    """
    baseline = 0.2
    additional = 0.8 * n_adults / (halfsat + n_adults) if (halfsat + n_adults) > 0 else 0.0
    return baseline + additional


def beverton_holt_recruitment(
    n_settlers: int,
    K: int,
    s0: float = 0.03,
) -> int:
    """Density-dependent settler-to-juvenile survival (Beverton-Holt).

    R = S × s0 / (1 + S / K)

    At low S: R ≈ S × s0 (supply-limited; each settler has s0 survival).
    At high S: R → K × s0 (habitat-limited; asymptote).

    The asymptote K × s0 represents the maximum annual recruitment rate
    at a node at full settlement pressure. For K=500, s0=0.03: max ~15
    recruits/year, which balances the ~10 annual adult deaths (0.02 × 500).

    Args:
        n_settlers: Number of settlers arriving.
        K: Available carrying capacity (slots).
        s0: Maximum per-settler survival (density-independent).

    Returns:
        Number of juveniles that survive first year.
    """
    if n_settlers <= 0 or K <= 0:
        return 0
    S = float(n_settlers)
    denominator = 1.0 + S / K
    R = S * s0 / denominator
    return max(0, int(np.round(R)))


def settle_recruits(
    agents: np.ndarray,
    genotypes: np.ndarray,
    settler_genotypes: np.ndarray,
    node_id: int,
    carrying_capacity: int,
    n_adults_present: int,
    habitat_area: float,
    effect_sizes: np.ndarray,
    w_od: float = 0.160,
    settler_survival: float = 0.03,
    rng: Optional[np.random.Generator] = None,
) -> int:
    """Add settled larvae to the node's agent arrays.

    Pipeline:
      1. Apply settlement cue modifier
      2. Beverton-Holt density-dependent recruitment
      3. Initialize new SETTLER-stage agents

    Args:
        agents: Agent structured array (modified in place).
        genotypes: Genotype array (modified in place).
        settler_genotypes: (n_settlers, N_LOCI, 2) int8 genotypes of arriving settlers.
        node_id: Node index.
        carrying_capacity: K for this node.
        n_adults_present: Adults present (for settlement cue).
        habitat_area: Habitat area in m² (for spatial placement).
        effect_sizes: (N_ADDITIVE,) float64 effect size vector.
        w_od: Overdominant weight for EF1A.
        settler_survival: Beverton-Holt s0.
        rng: Random generator.

    Returns:
        Number of recruits actually added.
    """
    if rng is None:
        rng = np.random.default_rng()

    n_arriving = len(settler_genotypes)
    if n_arriving == 0:
        return 0

    # Settlement cue Allee effect
    cue_mod = settlement_cue_modifier(n_adults_present)
    effective_settlers = max(0, int(n_arriving * cue_mod))

    # Density-dependent recruitment
    current_alive = int(np.sum(agents['alive']))
    available_K = max(0, carrying_capacity - current_alive)
    n_recruits = beverton_holt_recruitment(effective_settlers, available_K, settler_survival)
    n_recruits = min(n_recruits, effective_settlers, n_arriving)

    if n_recruits <= 0:
        return 0

    # Select which settlers survive (random subsample)
    if n_recruits < n_arriving:
        keep_idx = rng.choice(n_arriving, size=n_recruits, replace=False)
        settler_genotypes = settler_genotypes[keep_idx]

    # Find dead/empty slots in agent array
    dead_slots = np.where(~agents['alive'])[0]
    n_slots = min(n_recruits, len(dead_slots))
    if n_slots <= 0:
        return 0

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
        agents[slot]['fecundity_mod'] = 1.0
        agents[slot]['node_id'] = node_id
        agents[slot]['origin'] = 0  # WILD

        genotypes[slot] = settler_genotypes[j]

        # Compute resistance score
        agents[slot]['resistance'] = _compute_resistance(
            settler_genotypes[j], effect_sizes, w_od
        )

    return n_slots


def _compute_resistance(
    genotype: np.ndarray,
    effects: np.ndarray,
    w_od: float = 0.160,
) -> float:
    """Compute individual resistance score r_i from genotype.

    r_i = Σ ẽ_l × (g_l0 + g_l1)/2  +  δ_het × w_od

    Args:
        genotype: (N_LOCI, 2) int8 array.
        effects: (N_ADDITIVE,) float64 effect sizes.
        w_od: Overdominant weight.

    Returns:
        r_i ∈ [0, 1].
    """
    # Additive component
    allele_means = genotype[:N_ADDITIVE, :].sum(axis=1).astype(np.float64) * 0.5
    additive = float(np.dot(allele_means, effects))

    # Overdominant component
    ef1a_sum = int(genotype[IDX_EF1A, 0]) + int(genotype[IDX_EF1A, 1])
    od_bonus = w_od if ef1a_sum == 1 else 0.0

    return additive + od_bonus


# ═══════════════════════════════════════════════════════════════════════
# PER-CAPITA GROWTH RATE (for Allee effect analysis)
# ═══════════════════════════════════════════════════════════════════════


def per_capita_growth_rate(
    density: float,
    gamma_fert: float = 4.5,
    F0: float = 1.0e7,
    L_ref: float = 500.0,
    mean_size: float = 600.0,
    fecundity_exp: float = 2.5,
    mu_larva_daily: float = 0.05,
    pld_days: float = 63.0,
    settler_survival: float = 0.03,
    adult_mortality: float = 0.02,
    sex_ratio: float = 0.5,
) -> float:
    """Compute per-capita population growth rate at a given density.

    r(D) = F(D) × (eggs/2) × larval_surv × s0 − m

    Each individual (averaging over sexes) produces (1/sex_ratio) × F(ρ_male) ×
    fecundity(L) × larval_survival × settler_survival offspring per year.
    The Allee effect manifests as r(D) increasing with density: at low density,
    F(D) → 0 and r(D) drops dramatically, making the population vulnerable to
    stochastic extinction even though the deterministic rate may remain positive.

    For high-fecundity broadcast spawners, the deterministic Allee threshold
    (where r = 0) is extremely low — near zero density. The practical Allee
    effect operates through stochastic processes at low N.

    Args:
        density: Adult density (ind/m²).
        gamma_fert: Fertilization kinetics parameter.
        F0: Reference fecundity (eggs).
        L_ref: Reference size (mm).
        mean_size: Mean adult size (mm).
        fecundity_exp: Fecundity allometric exponent.
        mu_larva_daily: Daily larval mortality.
        pld_days: Pelagic larval duration.
        settler_survival: Annual settler survival.
        adult_mortality: Annual adult natural mortality.
        sex_ratio: Fraction male.

    Returns:
        Per-capita growth rate (can be negative at extremely low density).
    """
    male_density = density * sex_ratio
    F = fertilization_success(male_density, gamma_fert)
    eggs_per_female = fecundity(mean_size, F0, L_ref, fecundity_exp)
    larval_surv = larval_survival(pld_days, mu_larva_daily)

    # Per-capita births: half the population is female, each producing eggs
    per_capita_births = (1.0 - sex_ratio) * eggs_per_female * F * larval_surv * settler_survival
    return per_capita_births - adult_mortality
