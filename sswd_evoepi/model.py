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

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

from sswd_evoepi.config import (
    DiseaseSection,
    PopulationSection,
    SimulationConfig,
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
    update_resistance_scores,
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
    settlement_cue_modifier,
    srs_reproductive_lottery,
    _compute_resistance,
)
from sswd_evoepi.types import (
    AGENT_DTYPE,
    ANNUAL_SURVIVAL,
    IDX_EF1A,
    N_ADDITIVE,
    N_LOCI,
    STAGE_SIZE_THRESHOLDS,
    DiseaseState,
    Stage,
    allocate_agents,
    allocate_genotypes,
)


# ═══════════════════════════════════════════════════════════════════════
# CONSTANTS
# ═══════════════════════════════════════════════════════════════════════

DAYS_PER_YEAR = 365


# ═══════════════════════════════════════════════════════════════════════
# EFFECT SIZES (deterministic from seed)
# ═══════════════════════════════════════════════════════════════════════

def make_effect_sizes(seed: int = 12345) -> np.ndarray:
    """Create canonical effect-size vector (exponential, normalized).

    Sorted descending so first loci have largest effects.
    Total additive weight sums to W_add ≈ 0.840.

    CE-3: Exponential distribution per Lotterhos & Whitlock 2016.
    """
    rng = np.random.default_rng(seed)
    raw = rng.exponential(1.0, size=N_ADDITIVE)
    W_add = 0.840  # total additive contribution to r_i
    normalized = raw / raw.sum() * W_add
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

def initialize_population(
    n_individuals: int,
    max_agents: int,
    habitat_area: float,
    effect_sizes: np.ndarray,
    pop_cfg: PopulationSection,
    rng: np.random.Generator,
    q_init: float = 0.05,
    q_ef1a: float = 0.24,
    w_od: float = 0.160,
) -> Tuple[np.ndarray, np.ndarray]:
    """Initialize a population at demographic equilibrium.

    Creates individuals with age/size drawn from an approximate
    stable age distribution, random genotypes with given allele
    frequencies, and computed resistance scores.

    Args:
        n_individuals: Target number of live individuals.
        max_agents: Array capacity (must be >= n_individuals).
        habitat_area: Habitat area (m²) for spatial placement.
        effect_sizes: (N_ADDITIVE,) effect size vector.
        pop_cfg: Population configuration.
        rng: Random generator.
        q_init: Initial resistant allele frequency at additive loci.
        q_ef1a: Initial EF1A allele frequency.
        w_od: Overdominant weight for EF1A.

    Returns:
        (agents, genotypes) tuple.
    """
    agents = allocate_agents(max_agents)
    genotypes = allocate_genotypes(max_agents)

    hab_side = np.sqrt(max(habitat_area, 1.0))

    for i in range(n_individuals):
        # Sample age from approximate stable age distribution
        # Most individuals are young (high juvenile mortality)
        # but adults are long-lived
        u = rng.random()
        if u < 0.05:
            age = rng.uniform(0.0, 1.0)     # settlers
        elif u < 0.20:
            age = rng.uniform(1.0, 3.0)     # juveniles
        elif u < 0.35:
            age = rng.uniform(3.0, 7.0)     # subadults
        else:
            age = rng.uniform(7.0, 30.0)    # adults
        
        size = von_bertalanffy(age, pop_cfg.L_inf, pop_cfg.k_growth,
                               pop_cfg.t0_growth)
        # Add individual variation (±15%)
        size *= rng.lognormal(0.0, 0.10)
        size = max(0.5, size)

        stage = assign_stage(size, Stage.SETTLER)

        agents[i]['alive'] = True
        agents[i]['x'] = rng.uniform(0, hab_side)
        agents[i]['y'] = rng.uniform(0, hab_side)
        agents[i]['heading'] = rng.uniform(0, 2 * np.pi)
        agents[i]['speed'] = 0.1
        agents[i]['size'] = size
        agents[i]['age'] = age
        agents[i]['stage'] = stage
        agents[i]['sex'] = rng.integers(0, 2)
        agents[i]['disease_state'] = DiseaseState.S
        agents[i]['disease_timer'] = 0
        agents[i]['fecundity_mod'] = 1.0
        agents[i]['node_id'] = 0
        agents[i]['origin'] = 0  # WILD

        # Random genotypes
        for l in range(N_ADDITIVE):
            genotypes[i, l, 0] = 1 if rng.random() < q_init else 0
            genotypes[i, l, 1] = 1 if rng.random() < q_init else 0
        # EF1A
        genotypes[i, IDX_EF1A, 0] = 1 if rng.random() < q_ef1a else 0
        genotypes[i, IDX_EF1A, 1] = 1 if rng.random() < q_ef1a else 0
        # Fix lethal homozygotes
        if genotypes[i, IDX_EF1A, :].sum() == 2:
            genotypes[i, IDX_EF1A, 1] = 0

        agents[i]['resistance'] = _compute_resistance(
            genotypes[i], effect_sizes, w_od
        )

    return agents, genotypes


# ═══════════════════════════════════════════════════════════════════════
# ANNUAL DEMOGRAPHIC UPDATE
# ═══════════════════════════════════════════════════════════════════════

def annual_natural_mortality(
    agents: np.ndarray,
    pop_cfg: PopulationSection,
    rng: np.random.Generator,
) -> int:
    """Apply annual natural mortality to all alive individuals.

    Returns number killed.
    """
    alive_mask = agents['alive'].astype(bool)
    alive_idx = np.where(alive_mask)[0]
    n_killed = 0

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
            n_killed += 1

    return n_killed


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


def annual_reproduction(
    agents: np.ndarray,
    genotypes: np.ndarray,
    effect_sizes: np.ndarray,
    habitat_area: float,
    sst: float,
    carrying_capacity: int,
    pop_cfg: PopulationSection,
    rng: np.random.Generator,
    w_od: float = 0.160,
) -> dict:
    """Full annual reproduction: spawning → fertilization → SRS → settlement.

    Returns diagnostics dict.
    """
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
        agents[slot]['fecundity_mod'] = 1.0
        agents[slot]['node_id'] = 0
        agents[slot]['origin'] = 0

        genotypes[slot] = settler_geno[j]
        agents[slot]['resistance'] = _compute_resistance(
            settler_geno[j], effect_sizes, w_od
        )

    diag['n_recruits'] = n_slots
    return diag


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
    yearly_fert_success: Optional[np.ndarray] = None

    # Genetics tracking (Phase 7)
    yearly_allele_freq_top3: Optional[np.ndarray] = None   # (n_years, 3)
    yearly_ef1a_freq: Optional[np.ndarray] = None           # (n_years,)
    yearly_va: Optional[np.ndarray] = None                  # (n_years,) additive variance
    pre_epidemic_allele_freq: Optional[np.ndarray] = None   # (N_LOCI,) snapshot
    post_epidemic_allele_freq: Optional[np.ndarray] = None  # (N_LOCI,) snapshot

    # Daily timeseries (length = n_years * 365)
    daily_pop: Optional[np.ndarray] = None
    daily_infected: Optional[np.ndarray] = None
    daily_vibrio: Optional[np.ndarray] = None

    # Summary
    initial_pop: int = 0
    final_pop: int = 0
    min_pop: int = 0
    min_pop_year: int = 0
    total_disease_deaths: int = 0
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

    pop_cfg = config.population
    dis_cfg = config.disease

    rng = np.random.default_rng(seed)

    # Initialize effect sizes
    effect_sizes = make_effect_sizes(config.genetics.effect_size_seed)

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
        q_ef1a=config.genetics.q_ef1a_init,
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

    # ── Allocate result arrays ────────────────────────────────────────
    total_days = n_years * DAYS_PER_YEAR

    yearly_pop = np.zeros(n_years, dtype=np.int32)
    yearly_adults = np.zeros(n_years, dtype=np.int32)
    yearly_recruits = np.zeros(n_years, dtype=np.int32)
    yearly_natural_deaths = np.zeros(n_years, dtype=np.int32)
    yearly_disease_deaths = np.zeros(n_years, dtype=np.int32)
    yearly_mean_resistance = np.zeros(n_years, dtype=np.float64)
    yearly_fert_success = np.zeros(n_years, dtype=np.float64)

    # Genetics tracking (Phase 7)
    yearly_allele_freq_top3 = np.zeros((n_years, 3), dtype=np.float64)
    yearly_ef1a_freq = np.zeros(n_years, dtype=np.float64)
    yearly_va = np.zeros(n_years, dtype=np.float64)
    pre_epidemic_allele_freq = None
    post_epidemic_allele_freq = None
    epidemic_started_year = None  # track when epidemic begins for snapshot timing

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

    # ── Main simulation loop ─────────────────────────────────────────
    for year in range(n_years):
        year_disease_deaths_before = cumulative_disease_deaths

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
                )
                new_deaths = node_disease.cumulative_deaths - deaths_before
                cumulative_disease_deaths += new_deaths

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

        # B1. Natural mortality
        n_nat_dead = annual_natural_mortality(agents, pop_cfg, rng)

        # B2. Growth, aging, stage transitions
        annual_growth_and_aging(agents, pop_cfg, rng)

        # B3. Record pre-reproduction population
        alive_mask = agents['alive'].astype(bool)
        pop_now = int(np.sum(alive_mask))

        # B4. Reproduction (only if population > 0)
        if pop_now > 0:
            repro_diag = annual_reproduction(
                agents=agents,
                genotypes=genotypes,
                effect_sizes=effect_sizes,
                habitat_area=habitat_area,
                sst=T_celsius,
                carrying_capacity=carrying_capacity,
                pop_cfg=pop_cfg,
                rng=rng,
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
        yearly_disease_deaths[year] = (
            cumulative_disease_deaths - year_disease_deaths_before
        )
        yearly_fert_success[year] = repro_diag.get('fertilization_success', 0.0)

        # Mean resistance of living individuals
        if pop_after > 0:
            yearly_mean_resistance[year] = float(
                np.mean(agents['resistance'][alive_mask])
            )
        else:
            yearly_mean_resistance[year] = 0.0

        # ── Genetics tracking (Phase 7) ──────────────────────────────
        if pop_after > 0:
            allele_freq = compute_allele_frequencies(genotypes, alive_mask)
            yearly_allele_freq_top3[year] = allele_freq[:3]
            yearly_ef1a_freq[year] = float(allele_freq[IDX_EF1A])
            yearly_va[year] = compute_additive_variance(
                allele_freq, effect_sizes
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

    result = CoupledSimResult(
        n_years=n_years,
        yearly_pop=yearly_pop,
        yearly_adults=yearly_adults,
        yearly_recruits=yearly_recruits,
        yearly_natural_deaths=yearly_natural_deaths,
        yearly_disease_deaths=yearly_disease_deaths,
        yearly_mean_resistance=yearly_mean_resistance,
        yearly_fert_success=yearly_fert_success,
        # Genetics (Phase 7)
        yearly_allele_freq_top3=yearly_allele_freq_top3,
        yearly_ef1a_freq=yearly_ef1a_freq,
        yearly_va=yearly_va,
        pre_epidemic_allele_freq=pre_epidemic_allele_freq,
        post_epidemic_allele_freq=post_epidemic_allele_freq,
        # Daily
        daily_pop=daily_pop,
        daily_infected=daily_infected,
        daily_vibrio=daily_vibrio,
        # Summary
        initial_pop=initial_pop,
        final_pop=final_pop,
        min_pop=min_pop,
        min_pop_year=min_pop_year,
        total_disease_deaths=cumulative_disease_deaths,
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
    yearly_mean_resistance: Optional[np.ndarray] = None
    yearly_vibrio_max: Optional[np.ndarray] = None
    # Per-node genetics: shape (n_nodes, n_years)
    yearly_ef1a_freq: Optional[np.ndarray] = None
    yearly_va: Optional[np.ndarray] = None
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
    # Summary
    initial_total_pop: int = 0
    final_total_pop: int = 0
    disease_year: Optional[int] = None
    seed: int = 0


def run_spatial_simulation(
    network: 'MetapopulationNetwork',
    n_years: int = 10,
    disease_year: Optional[int] = None,
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
        disease_year: Year to introduce disease (None = no disease).
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
    from sswd_evoepi.environment import sst_with_trend, seasonal_flushing
    from sswd_evoepi.movement import daily_movement, InfectedDensityGrid

    if config is None:
        config = default_config()

    pop_cfg = config.population
    dis_cfg = config.disease
    mov_cfg = config.movement
    mov_enabled = mov_cfg.enabled
    spatial_tx = mov_enabled and mov_cfg.spatial_transmission

    rng = np.random.default_rng(seed)
    N = network.n_nodes

    # ── Initialize populations at each node ──────────────────────────
    effect_sizes = make_effect_sizes(config.genetics.effect_size_seed)

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
            q_ef1a=config.genetics.q_ef1a_init,
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
    yearly_mean_r = np.zeros((N, n_years), dtype=np.float64)
    yearly_vibrio_max = np.zeros((N, n_years), dtype=np.float64)
    yearly_total_pop = np.zeros(n_years, dtype=np.int32)
    yearly_total_larvae = np.zeros(n_years, dtype=np.int32)

    # Genetics tracking per node per year
    yearly_ef1a_freq = np.zeros((N, n_years), dtype=np.float64)
    yearly_va = np.zeros((N, n_years), dtype=np.float64)
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

    initial_total = network.total_population()
    node_names = [n.definition.name for n in network.nodes]
    node_K = np.array([n.definition.carrying_capacity for n in network.nodes])

    # Per-node disease tracking
    node_disease_states = []
    disease_active_flags = [False] * N
    cumulative_dis_deaths = [0] * N

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

    # ── Main simulation loop ─────────────────────────────────────────
    start_year = 2000  # reference year for SST

    for year in range(n_years):
        cal_year = start_year + year
        year_dd_before = list(cumulative_dis_deaths)

        # Track max vibrio per node this year
        max_vibrio_year = np.zeros(N, dtype=np.float64)

        # ── DAILY LOOP (365 days) ────────────────────────────────────
        for day in range(DAYS_PER_YEAR):
            month = min(day // 30, 11)

            # 1. Update environment at each node
            for i, node in enumerate(network.nodes):
                nd = node.definition
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
                )
                new_deaths = (
                    node_disease_states[i].cumulative_deaths - deaths_before
                )
                cumulative_dis_deaths[i] += new_deaths

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

            # Optional individual-level snapshot recording
            if snapshot_recorder is not None:
                sim_day = year * DAYS_PER_YEAR + day
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

        # 1. Natural mortality + growth at each node
        node_nat_deaths = np.zeros(N, dtype=np.int32)
        for i, node in enumerate(network.nodes):
            nd = node.agents
            n_killed = annual_natural_mortality(nd, pop_cfg, rng)
            node_nat_deaths[i] = n_killed
            annual_growth_and_aging(nd, pop_cfg, rng)

        # 2. Reproduction at each node → produce larval cohorts
        source_ids = []
        source_n = []
        source_geno = []
        node_repro_diags = []

        for i, node in enumerate(network.nodes):
            pop_now = node.n_alive
            if pop_now == 0:
                node_repro_diags.append({'n_recruits': 0})
                continue
            repro_diag = annual_reproduction(
                agents=node.agents,
                genotypes=node.genotypes,
                effect_sizes=effect_sizes,
                habitat_area=node.definition.habitat_area,
                sst=node.current_sst,
                carrying_capacity=node.definition.carrying_capacity,
                pop_cfg=pop_cfg,
                rng=rng,
            )
            node_repro_diags.append(repro_diag)

            # Collect larvae for dispersal
            n_comp = repro_diag.get('n_competent', 0)
            if n_comp > 0:
                # Get genotypes of alive individuals for sampling
                alive_mask = node.agents['alive'].astype(bool)
                alive_geno = node.genotypes[alive_mask]
                # Sample n_comp genotypes (representing competent larvae)
                if len(alive_geno) > 0:
                    idx = rng.choice(len(alive_geno), size=min(n_comp, 1000),
                                     replace=True)
                    source_ids.append(i)
                    source_n.append(min(n_comp, 1000))
                    source_geno.append(alive_geno[idx])

        # 3. Larval dispersal via C matrix
        total_larvae = sum(source_n)
        if total_larvae > 0:
            settler_map = distribute_larvae(
                source_ids, source_n, source_geno,
                network.C, rng,
            )

            # 4. Settle larvae at destination nodes
            for k in range(N):
                if not settler_map[k]:
                    continue
                node = network.nodes[k]
                nd = node.definition
                for geno_batch, src in settler_map[k]:
                    n_settlers = len(geno_batch)
                    # Find empty slots
                    dead_slots = np.where(~node.agents['alive'])[0]
                    n_to_place = min(n_settlers, len(dead_slots))
                    if n_to_place == 0:
                        continue
                    # B-H density check
                    current_alive = node.n_alive
                    available = max(0, nd.carrying_capacity - current_alive)
                    n_to_place = min(n_to_place, available)
                    if n_to_place <= 0:
                        continue

                    hab_side = np.sqrt(max(nd.habitat_area, 1.0))
                    for j in range(n_to_place):
                        slot = dead_slots[j]
                        node.agents[slot]['alive'] = True
                        node.agents[slot]['x'] = rng.uniform(0, hab_side)
                        node.agents[slot]['y'] = rng.uniform(0, hab_side)
                        node.agents[slot]['heading'] = rng.uniform(0, 2*np.pi)
                        node.agents[slot]['speed'] = 0.1
                        node.agents[slot]['size'] = 0.5
                        node.agents[slot]['age'] = 0.0
                        node.agents[slot]['stage'] = Stage.SETTLER
                        node.agents[slot]['sex'] = rng.integers(0, 2)
                        node.agents[slot]['disease_state'] = DiseaseState.S
                        node.agents[slot]['disease_timer'] = 0
                        node.agents[slot]['fecundity_mod'] = 1.0
                        node.agents[slot]['node_id'] = k
                        node.agents[slot]['origin'] = 0
                        node.genotypes[slot] = geno_batch[j]
                        node.agents[slot]['resistance'] = _compute_resistance(
                            geno_batch[j], effect_sizes, 0.160,
                        )

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
            yearly_recruits[i, year] = node_repro_diags[i].get('n_recruits', 0)
            yearly_nat_deaths[i, year] = node_nat_deaths[i]
            yearly_dis_deaths[i, year] = (
                cumulative_dis_deaths[i] - year_dd_before[i]
            )
            yearly_vibrio_max[i, year] = max_vibrio_year[i]
            if pop_now > 0:
                yearly_mean_r[i, year] = float(
                    np.mean(node.agents['resistance'][alive_mask])
                )
                # Genetics tracking
                af = compute_allele_frequencies(node.genotypes, alive_mask)
                yearly_allele_freq_top3[i, year] = af[:3]
                yearly_ef1a_freq[i, year] = float(af[IDX_EF1A])
                yearly_va[i, year] = compute_additive_variance(af, effect_sizes)

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
        yearly_mean_resistance=yearly_mean_r,
        yearly_vibrio_max=yearly_vibrio_max,
        yearly_ef1a_freq=yearly_ef1a_freq,
        yearly_va=yearly_va,
        yearly_allele_freq_top3=yearly_allele_freq_top3,
        yearly_ne_ratio=yearly_ne_ratio,
        pre_epidemic_allele_freq=pre_epidemic_af if pre_snapshot_taken else None,
        post_epidemic_allele_freq=post_epidemic_af if post_snapshot_taken else None,
        yearly_total_pop=yearly_total_pop,
        yearly_total_larvae_dispersed=yearly_total_larvae,
        peak_disease_prevalence=peak_disease_prev,
        initial_total_pop=initial_total,
        final_total_pop=int(yearly_total_pop[-1]),
        disease_year=disease_year,
        seed=seed,
        snapshot_recorder=snapshot_recorder,
    )
