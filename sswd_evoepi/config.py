"""Configuration system for SSWD-EvoEpi.

Hierarchical YAML configuration with deep-merge support:
  base.yaml → scenario override → climate override → sweep overrides

References:
  - integration-architecture-spec.md §5
  - data-parameterization-plan.md §1 (authoritative parameter inventory)

Design decisions (from Willem):
  - NO cost_of_resistance parameter (CE-1)
  - Both etiological scenarios supported: "ubiquitous" | "invasion" (CE-2)
  - Exponential decay for effect sizes (CE-3)
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import yaml


# ═══════════════════════════════════════════════════════════════════════
# CONFIGURATION DATACLASSES
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class SimulationSection:
    """Top-level simulation timing and control."""
    start_year: int = 1910
    epidemic_year: int = 2013
    end_year: int = 2100
    spinup_years: int = 100
    seed: int = 42
    checkpoint_interval: int = 10
    sst_source: str = 'sinusoidal'  # 'sinusoidal', 'satellite', or 'monthly'
    sst_data_dir: str = 'data/sst'  # directory with *_climatology.csv / *_monthly.csv files
    sst_start_year: int = 2002     # start calendar year for 'monthly' SST source
    sst_scenario: str = 'observed_only'  # 'observed_only', 'ssp245', 'ssp585', etc.
    sst_projection_dir: str = 'data/sst/projections'  # directory with *_{scenario}_monthly.csv
    sst_obs_end_year: int = 2025   # last year of OISST observations


@dataclass
class SpatialSection:
    """Spatial network and connectivity parameters."""
    node_file: str = "data/nodes/node_definitions.csv"
    forcing_dir: str = "data/forcing/"
    C_matrix_file: str = "data/connectivity/C_matrix.npz"
    D_matrix_file: str = "data/connectivity/D_matrix.npz"
    D_L: float = 400.0           # Larval dispersal scale (km)
    D_P: float = 15.0            # Pathogen dispersal scale (km)
    D_P_max_range: float = 50.0  # Maximum pathogen dispersal range (km); cutoff for D matrix
    r_total: float = 0.02        # Total settlement success fraction
    alpha_self_fjord: float = 0.30   # Larval self-recruitment fraction for fjord nodes
    alpha_self_open: float = 0.10    # Larval self-recruitment fraction for open coast


@dataclass
class PopulationSection:
    """Population dynamics parameters."""
    L_inf: float = 1000.0        # VB asymptotic size (mm)
    k_growth: float = 0.08       # VB growth rate (yr⁻¹)
    t0_growth: float = -0.5      # VB age offset (yr)
    F0: float = 1.0e7            # Reference fecundity (eggs)
    fecundity_exp: float = 2.5   # Allometric exponent
    L_ref: float = 500.0         # Reference size for fecundity (mm)
    L_min_repro: float = 400.0   # Minimum reproductive size (mm)
    gamma_fert: float = 4.5      # Fertilization kinetics parameter (m²)
    alpha_srs: float = 1.35      # Pareto SRS shape
    sigma_alpha: float = 0.10    # Annual SRS shape variation
    settler_survival: float = 0.03  # Annual settler survival (B-H s₀)
    annual_survival: List[float] = field(
        default_factory=lambda: [0.001, 0.03, 0.90, 0.95, 0.98]
    )
    senescence_age: float = 50.0
    senescence_mortality: float = 0.10
    urchin_submodel: bool = False


@dataclass
class SpawningSection:
    """Spawning system parameters.
    
    Extended spawning season with cascading spawning bouts, spawning gravity,
    and post-spawning immunosuppression. Replaces single-day pulse model.
    
    References:
      - spawning-overhaul-spec.md (all parameters with confidence ratings)
    """
    # Spawning season parameters
    season_start_doy: int = 305      # Season start (~Nov 1)
    season_end_doy: int = 196        # Season end (~Jul 15) - wraps year boundary
    peak_doy: int = 105              # Peak spawning activity (~Apr 15)
    peak_width_days: float = 60.0    # Std dev of seasonal peak (Normal) - calibrated in Phase 1A2
    lat_shift_per_deg: float = 3.0   # Latitude shift of peak (days/°N)
    
    # Spontaneous spawning rates (base daily probabilities, modulated by season)
    p_spontaneous_female: float = 0.012   # Daily probability ready female spawns spontaneously - calibrated for 80%+ participation  
    p_spontaneous_male: float = 0.0125    # Daily probability ready male initiates bout spontaneously - calibrated for 2.2 mean bouts
    
    # Cascade induction parameters (Phase 2 - not used in Phase 1)
    induction_female_to_male: float = 0.80    # Probability male spawns when female nearby has spawned
    induction_male_to_female: float = 0.60    # Probability female spawns when male nearby has spawned  
    cascade_window: int = 3                   # Duration of chemical spawning cue persistence (days)
    cascade_radius: float = 200.0  # Effective range of chemical spawning cue (m)
    
    # Multi-bout parameters
    female_max_bouts: int = 2                 # Maximum spawning bouts per female per season
    male_max_bouts: int = 3                   # Maximum spawning bouts per male per season
    male_refractory_days: int = 0             # Minimum days between male spawning bouts (0 = next day OK)
    
    # Spawning gravity parameters (Phase 3 - not used in Phase 1)
    gravity_enabled: bool = True              # Enable pre-spawning aggregation movement
    gravity_strength: float = 0.3             # Maximum speed bias toward conspecifics (m/min)
    gravity_range: float = 100.0              # Sensory detection range for conspecifics (m)
    pre_spawn_gravity_days: int = 14          # Days before readiness that gravity activates
    post_spawn_gravity_days: int = 14         # Days after spawning that gravity persists
    
    # Readiness induction parameters (Phase 12)
    readiness_induction_radius: float = 300.0   # Meters — chemical detection range for spawning cues
    readiness_induction_prob: float = 0.5        # Probability of induced readiness per day if within radius


@dataclass
class DiseaseSection:
    """Disease dynamics parameters.

    scenario: "ubiquitous" — pathogen always present, SST-triggered epidemics
              "invasion"   — pathogen absent until explicit introduction
    """
    scenario: str = "ubiquitous"

    # Transmission
    a_exposure: float = 0.75      # Exposure rate (d⁻¹)
    K_half: float = 87000.0       # Half-infective dose (bact/mL)

    # Shedding (field-effective; ERRATA E2)
    sigma_1_eff: float = 5.0      # I₁ shedding (bact/mL/d/host)
    sigma_2_eff: float = 50.0     # I₂ shedding
    sigma_D: float = 15.0         # Saprophytic burst, field-effective (CE-6)
    Ea_sigma: float = 5000.0      # Shedding E_a/R (K)

    # Disease progression rates at T_ref=20°C
    # Calibrated to Prentice et al. 2025 (Nature E&E): mean 11.6d exposure→death,
    # 5.6d symptoms→death at ~13°C, Arrhenius-corrected to 20°C reference.
    mu_EI1_ref: float = 0.233     # E→I₁ rate (d⁻¹); 4.3d at 20°C, 6.0d at 13°C
    mu_I1I2_ref: float = 0.434    # I₁→I₂ rate (d⁻¹); 2.3d at 20°C, 3.5d at 13°C
    mu_I2D_ref: float = 0.563     # I₂→D rate (d⁻¹); 1.8d at 20°C, 2.1d at 13°C

    # Activation energies (E_a/R in K)
    Ea_EI1: float = 4000.0
    Ea_I1I2: float = 5000.0
    Ea_I2D: float = 2000.0        # ERRATA E1: was 6000

    # Recovery
    rho_rec: float = 0.05         # Base recovery rate (d⁻¹) — no empirical basis

    # Environmental pathogen
    P_env_max: float = 500.0      # Background Vibrio input (bact/mL/d)

    # Temperature
    T_vbnc: float = 12.0          # VBNC midpoint (°C)
    T_ref: float = 20.0           # V. pectenicida T_opt (°C)

    # Salinity
    s_min: float = 10.0           # Salinity minimum for Vibrio (psu)
    s_full: float = 28.0          # Full-marine salinity (psu)

    # Post-spawning immunosuppression (Phase 4)
    immunosuppression_enabled: bool = True    # Enable post-spawning susceptibility increase
    susceptibility_multiplier: float = 2.0   # Force-of-infection multiplier during immunosuppression
    immunosuppression_duration: int = 28      # Duration of post-spawning immunosuppression (days)

    # Tolerance (mirrored from GeneticsSection for disease-module convenience)
    tau_max: float = 0.85             # Max I₂→D mortality reduction at t_i=1

    # Juvenile immunity (Phase 11)
    min_susceptible_age_days: int = 0  # Days post-settlement before susceptible
                                       # 0 = immediate (backward compatible)
                                       # SA range: [0, 180]

    # Invasion scenario extras
    invasion_year: Optional[int] = None
    invasion_nodes: Optional[List[int]] = None

    # Wavefront disease spread (spatial spread from origin)
    wavefront_enabled: bool = False           # If True, P_env gated per-node; disease spreads as wavefront
    disease_origin_nodes: Optional[List[int]] = None  # Node IDs where disease starts. None = all nodes (backward compat)
    activation_threshold: float = 1.0         # Vibrio concentration (bact/mL) that triggers node activation


@dataclass
class GeneticsSection:
    """Genetics and evolution parameters — three-trait architecture.

    51 loci partitioned into resistance (R), tolerance (T), recovery (C).
    Default: 17/17/17. Constraint: n_resistance + n_tolerance + n_recovery = 51.

    Allele frequency initialization:
      - q_init_mode="uniform": all loci start at q_init (old behavior)
      - q_init_mode="beta": per-locus frequencies drawn from Beta(q_init_beta_a,
        q_init_beta_b), then scaled so population-mean trait ≈ target_mean.
        This produces realistic among-locus variation (some loci common, some rare),
        consistent with independent evolutionary histories of immune genes.

    References:
      - three-trait-genetic-architecture-spec.md §4
    """
    # Trait partition (must sum to 51)
    n_resistance: int = 17        # Loci 0 .. n_resistance-1
    n_tolerance: int = 17         # Loci n_resistance .. n_resistance+n_tolerance-1
    n_recovery: int = 17          # Loci n_resistance+n_tolerance .. sum-1

    # Per-trait initialization targets
    target_mean_r: float = 0.15   # Population-mean resistance at t=0
    target_mean_t: float = 0.10   # Population-mean tolerance at t=0
    target_mean_c: float = 0.02   # Population-mean recovery at t=0 — rare standing variation

    # Tolerance mechanics
    tau_max: float = 0.85         # Max I₂→D mortality reduction at t_i=1

    # Shared genetics parameters
    mu_per_locus: float = 1.0e-8  # Mutation rate (per allele per generation)
    n_bank: int = 100             # Tier 2 genotype bank size
    effect_size_seed: int = 12345 # Seed for reproducible effect sizes
    q_init_mode: str = "beta"     # "uniform" or "beta"
    q_init_beta_a: float = 2.0    # Beta shape a (only used if mode="beta")
    q_init_beta_b: float = 8.0    # Beta shape b (only used if mode="beta")


@dataclass
class MovementSection:
    """Agent movement and spatial transmission parameters.

    Movement: correlated random walk (CRW).
    Spatial transmission: grid-based local Vibrio exposure — susceptibles
    near infected agents get proportionally higher exposure.

    References:
      - modeling-approach.md §2.4 (hourly substeps)
      - disease-module-spec §5.5 (speed modifiers by disease state)
    """
    enabled: bool = True           # Enable CRW movement
    base_speed: float = 0.5        # m/min (Pycnopodia undisturbed; Kay & Emlet 2002)
    sigma_turn: float = 0.6        # Turning angle std dev (radians)
    substeps_per_day: int = 24     # Sub-steps per day (24 = hourly)
    spatial_transmission: bool = True   # Use local Vibrio exposure grid
    cell_size: float = 20.0        # Grid cell size for spatial transmission (m)
    diffusion_passes: int = 2      # Smoothing passes (3×3 averaging)
    speed_sigma: float = 0.0       # Log-normal std dev for step-length variability (0 = fixed steps)


@dataclass
class ReleaseEvent:
    """A captive-bred reintroduction event.

    Specifies when, where, and how many individuals to release,
    along with their genetic composition and demographics.

    genetics_mode controls genotype assignment:
      'allele_freqs': Per-locus protective allele frequencies (allele_freqs).
      'trait_targets': Target mean trait values, internally converted to
                       allele frequencies (trait_targets dict).
      'genotypes': Explicit genotype array provided directly.
    """
    time_step: int                          # Simulation day (0-indexed from start)
    node_id: int = 0                        # Target node index
    n_individuals: int = 100                # Number of individuals to release
    genetics_mode: str = 'trait_targets'    # 'allele_freqs', 'trait_targets', 'genotypes'
    allele_freqs: Optional[np.ndarray] = None    # (N_LOCI,) per-locus frequencies
    trait_targets: Optional[Dict[str, float]] = None  # e.g. {'resistance': 0.3, ...}
    genotypes: Optional[np.ndarray] = None       # (n_individuals, N_LOCI, 2) int8
    age_range: Tuple[int, int] = (365, 730)      # (min_days, max_days) for released ages
    mark_released: bool = True                   # Tag with Origin.CAPTIVE_BRED


@dataclass
class ConservationSection:
    """Conservation intervention parameters."""
    enabled: bool = False
    schedule_file: Optional[str] = None
    breeding_protocol: str = "standard"
    n_broodstock: int = 130
    production_scaling: bool = True
    releases: List[Dict[str, Any]] = field(default_factory=list)


@dataclass
class PathogenEvolutionSection:
    """Pathogen co-evolution parameters.

    When enabled, V. pectenicida evolves virulence along a mechanistic
    tradeoff curve. Host-pathogen coevolution emerges from the interaction
    of host resistance evolution (genetics module) with pathogen virulence
    evolution (this module).

    References:
      - pathogen-evolution-spec.md §5.1
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


@dataclass
class OutputSection:
    """Output control."""
    directory: str = "results/"
    weekly_state: bool = True
    daily_disease: bool = True
    annual_genetics: bool = True
    decadal_genotype_snapshots: bool = True
    netcdf_compression: int = 4


@dataclass
class SimulationConfig:
    """Complete simulation configuration.

    Load from YAML via `load_config()`. Sections map 1:1 to YAML top-level keys.
    """
    simulation: SimulationSection = field(default_factory=SimulationSection)
    spatial: SpatialSection = field(default_factory=SpatialSection)
    population: PopulationSection = field(default_factory=PopulationSection)
    spawning: SpawningSection = field(default_factory=SpawningSection)
    disease: DiseaseSection = field(default_factory=DiseaseSection)
    genetics: GeneticsSection = field(default_factory=GeneticsSection)
    pathogen_evolution: PathogenEvolutionSection = field(default_factory=PathogenEvolutionSection)
    movement: MovementSection = field(default_factory=MovementSection)
    conservation: ConservationSection = field(default_factory=ConservationSection)
    output: OutputSection = field(default_factory=OutputSection)
    release_events: List[ReleaseEvent] = field(default_factory=list)


# ═══════════════════════════════════════════════════════════════════════
# YAML LOADING & MERGING
# ═══════════════════════════════════════════════════════════════════════

def deep_merge(base: Dict, override: Dict) -> Dict:
    """Recursively merge override into base. Modifies base in place.

    - Dict values are merged recursively
    - Non-dict values are replaced
    - Keys in override but not base are added

    Args:
        base: Base dictionary (modified in place).
        override: Override dictionary.

    Returns:
        The merged base dictionary.
    """
    for key, value in override.items():
        if (
            key in base
            and isinstance(base[key], dict)
            and isinstance(value, dict)
        ):
            deep_merge(base[key], value)
        else:
            base[key] = value
    return base


def _dict_to_section(section_cls, data: Dict) -> Any:
    """Convert a dict to a dataclass, ignoring unknown keys."""
    import dataclasses
    valid_fields = {f.name for f in dataclasses.fields(section_cls)}
    filtered = {k: v for k, v in data.items() if k in valid_fields}
    return section_cls(**filtered)


def _yaml_to_config(data: Dict) -> SimulationConfig:
    """Convert a merged YAML dict to a SimulationConfig."""
    sections = {}
    section_map = {
        'simulation': SimulationSection,
        'spatial': SpatialSection,
        'population': PopulationSection,
        'spawning': SpawningSection,
        'disease': DiseaseSection,
        'genetics': GeneticsSection,
        'pathogen_evolution': PathogenEvolutionSection,
        'movement': MovementSection,
        'conservation': ConservationSection,
        'output': OutputSection,
    }
    for key, cls in section_map.items():
        if key in data and isinstance(data[key], dict):
            sections[key] = _dict_to_section(cls, data[key])
        else:
            sections[key] = cls()

    # Parse release events (top-level list, not a section)
    release_events = []
    if 'release_events' in data and isinstance(data['release_events'], list):
        for re_dict in data['release_events']:
            if isinstance(re_dict, dict):
                re_dict = dict(re_dict)  # don't mutate original
                if 'age_range' in re_dict and isinstance(re_dict['age_range'], list):
                    re_dict['age_range'] = tuple(re_dict['age_range'])
                release_events.append(_dict_to_section(ReleaseEvent, re_dict))
    sections['release_events'] = release_events

    return SimulationConfig(**sections)


def validate_config(config: SimulationConfig) -> None:
    """Validate configuration constraints. Raises ValueError on failure.

    Checks:
      - Disease scenario is valid
      - Year ordering is consistent
      - Genetic parameters are consistent
      - No cost_of_resistance sneaking in
    """
    # SST source
    valid_sst_sources = {"sinusoidal", "satellite", "monthly"}
    if config.simulation.sst_source not in valid_sst_sources:
        raise ValueError(
            f"simulation.sst_source must be one of {valid_sst_sources}, "
            f"got '{config.simulation.sst_source}'"
        )
    if config.simulation.sst_source in ("satellite", "monthly"):
        import os
        import warnings
        if not os.path.isdir(config.simulation.sst_data_dir):
            warnings.warn(
                f"simulation.sst_data_dir '{config.simulation.sst_data_dir}' "
                f"does not exist. SST loading will fail at runtime.",
                UserWarning,
                stacklevel=2,
            )

    # SST scenario validation
    valid_sst_scenarios = {"observed_only", "ssp126", "ssp245", "ssp370", "ssp585"}
    if config.simulation.sst_scenario not in valid_sst_scenarios:
        raise ValueError(
            f"simulation.sst_scenario must be one of {valid_sst_scenarios}, "
            f"got '{config.simulation.sst_scenario}'"
        )
    # If scenario != observed_only and sim extends beyond obs, projection dir must exist
    if config.simulation.sst_scenario != 'observed_only':
        import os
        import warnings
        if not os.path.isdir(config.simulation.sst_projection_dir):
            warnings.warn(
                f"simulation.sst_projection_dir "
                f"'{config.simulation.sst_projection_dir}' does not exist. "
                f"SST projection loading will fail at runtime.",
                UserWarning,
                stacklevel=2,
            )

    # Disease scenario
    valid_scenarios = {"ubiquitous", "invasion"}
    if config.disease.scenario not in valid_scenarios:
        raise ValueError(
            f"disease.scenario must be one of {valid_scenarios}, "
            f"got '{config.disease.scenario}'"
        )

    # Invasion scenario requires year and nodes
    if config.disease.scenario == "invasion":
        if config.disease.invasion_year is None:
            raise ValueError(
                "disease.invasion_year required for 'invasion' scenario"
            )

    # Wavefront validation
    if config.disease.wavefront_enabled:
        if config.disease.disease_origin_nodes is None or len(config.disease.disease_origin_nodes) == 0:
            raise ValueError(
                "disease.disease_origin_nodes required when wavefront_enabled=True"
            )
        if config.disease.activation_threshold <= 0:
            raise ValueError(
                "disease.activation_threshold must be positive"
            )

    # Year ordering
    if config.simulation.start_year >= config.simulation.end_year:
        raise ValueError(
            f"start_year ({config.simulation.start_year}) must be < "
            f"end_year ({config.simulation.end_year})"
        )
    if config.simulation.epidemic_year < config.simulation.start_year:
        raise ValueError(
            f"epidemic_year ({config.simulation.epidemic_year}) must be >= "
            f"start_year ({config.simulation.start_year})"
        )

    # Genetics — three-trait partition must sum to 51
    g = config.genetics
    partition_sum = g.n_resistance + g.n_tolerance + g.n_recovery
    if partition_sum != 51:
        raise ValueError(
            f"genetics partition must sum to 51, got "
            f"{g.n_resistance}+{g.n_tolerance}+{g.n_recovery}={partition_sum}"
        )
    if g.n_resistance < 1 or g.n_tolerance < 1 or g.n_recovery < 1:
        raise ValueError(
            f"Each trait must have ≥1 locus, got "
            f"R={g.n_resistance}, T={g.n_tolerance}, C={g.n_recovery}"
        )

    # Annual survival array length
    if len(config.population.annual_survival) != 5:
        raise ValueError(
            f"population.annual_survival must have 5 elements (one per stage), "
            f"got {len(config.population.annual_survival)}"
        )

    # Pathogen evolution
    pe = config.pathogen_evolution
    if pe.alpha_kill <= 0:
        raise ValueError(
            f"pathogen_evolution.alpha_kill must be > 0, got {pe.alpha_kill}"
        )
    if pe.sigma_v_mutation < 0:
        raise ValueError(
            f"pathogen_evolution.sigma_v_mutation must be >= 0, "
            f"got {pe.sigma_v_mutation}"
        )
    if pe.v_min >= pe.v_max:
        raise ValueError(
            f"pathogen_evolution.v_min ({pe.v_min}) must be < "
            f"v_max ({pe.v_max})"
        )
    if not (pe.v_min <= pe.v_init <= pe.v_max):
        raise ValueError(
            f"pathogen_evolution.v_init ({pe.v_init}) must be in "
            f"[v_min={pe.v_min}, v_max={pe.v_max}]"
        )

    # Positive parameters
    if config.simulation.seed < 0:
        raise ValueError("simulation.seed must be non-negative")
    if config.population.L_inf <= 0:
        raise ValueError("population.L_inf must be positive")
    if config.disease.K_half <= 0:
        raise ValueError("disease.K_half must be positive")

    # Release events
    valid_genetics_modes = {'allele_freqs', 'trait_targets', 'genotypes'}
    for i, re in enumerate(config.release_events):
        if re.time_step < 0:
            raise ValueError(
                f"release_events[{i}].time_step must be >= 0, got {re.time_step}"
            )
        if re.n_individuals < 1:
            raise ValueError(
                f"release_events[{i}].n_individuals must be >= 1, "
                f"got {re.n_individuals}"
            )
        if re.genetics_mode not in valid_genetics_modes:
            raise ValueError(
                f"release_events[{i}].genetics_mode must be one of "
                f"{valid_genetics_modes}, got '{re.genetics_mode}'"
            )
        if re.genetics_mode == 'allele_freqs' and re.allele_freqs is None:
            raise ValueError(
                f"release_events[{i}]: allele_freqs required when "
                f"genetics_mode='allele_freqs'"
            )
        if re.genetics_mode == 'genotypes' and re.genotypes is None:
            raise ValueError(
                f"release_events[{i}]: genotypes required when "
                f"genetics_mode='genotypes'"
            )
        if (len(re.age_range) != 2
                or re.age_range[0] < 0
                or re.age_range[0] > re.age_range[1]):
            raise ValueError(
                f"release_events[{i}].age_range must be (min, max) with "
                f"0 <= min <= max, got {re.age_range}"
            )


def load_config(
    base_path: Union[str, Path],
    scenario_path: Optional[Union[str, Path]] = None,
    climate_path: Optional[Union[str, Path]] = None,
    sweep_overrides: Optional[Dict] = None,
) -> SimulationConfig:
    """Load and merge hierarchical YAML configuration.

    Merge order: base → scenario → climate → sweep overrides.
    Each layer overrides only the fields it specifies.

    Args:
        base_path: Path to base configuration YAML.
        scenario_path: Optional scenario override YAML.
        climate_path: Optional climate override YAML.
        sweep_overrides: Optional dict of parameter sweep overrides.

    Returns:
        Validated SimulationConfig.

    Raises:
        FileNotFoundError: If base_path doesn't exist.
        ValueError: If validation fails.
    """
    base_path = Path(base_path)
    if not base_path.exists():
        raise FileNotFoundError(f"Config file not found: {base_path}")

    with open(base_path) as f:
        config_dict = yaml.safe_load(f) or {}

    if scenario_path is not None:
        scenario_path = Path(scenario_path)
        if scenario_path.exists():
            with open(scenario_path) as f:
                scenario = yaml.safe_load(f) or {}
            deep_merge(config_dict, scenario)

    if climate_path is not None:
        climate_path = Path(climate_path)
        if climate_path.exists():
            with open(climate_path) as f:
                climate = yaml.safe_load(f) or {}
            deep_merge(config_dict, climate)

    if sweep_overrides is not None:
        deep_merge(config_dict, sweep_overrides)

    config = _yaml_to_config(config_dict)
    validate_config(config)
    return config


def default_config() -> SimulationConfig:
    """Return a SimulationConfig with all default values."""
    config = SimulationConfig()
    validate_config(config)
    return config
