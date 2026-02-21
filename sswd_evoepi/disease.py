"""Disease dynamics module — SEIPD+R compartmental model for SSWD.

Implements:
  - SEIPD+R compartments: S→E→I₁→I₂→D (+ recovery from I₁/I₂ → R)
  - Force of infection: λ_i = a × P/(K_half+P) × (1−r_i) × S_sal × f_size(L_i)
  - Three-trait host genetics:
      * Resistance (r_i): reduces force of infection
      * Tolerance (t_i): timer-scales I₂ duration (longer survival while infected)
      * Recovery (c_i): linear clearance probability from I₁/I₂ → R
  - Temperature-dependent progression rates (Arrhenius, T_ref=20°C)
  - Environmental pathogen dynamics: shedding − decay − flushing + reservoir
  - VBNC sigmoidal reservoir: P_env(T)
  - Salinity modifier for fjord protection
  - Size-dependent susceptibility (Eisenlord OR=1.23 per 10mm)
  - Both etiological scenarios via config: ubiquitous (default) vs invasion
  - Stochastic S→E transitions: prob = 1 − exp(−λ_i × Δt)
  - Erlang-distributed stage durations (sampled countdown timers)
  - Carcass tracking for saprophytic shedding

References:
  - disease-module-spec.md (authoritative specification)
  - three-trait-genetic-architecture-spec.md §6
  - ERRATA E1: E_a/R for I₂→D = 2,000 K (not 6,000)
  - ERRATA E2: Field-effective shedding σ₁=5, σ₂=50 at 20°C
  - ERRATA E14: σ_D = 150 bact/mL/d/carcass
  - CODE_ERRATA CE-2: ubiquitous vs invasion via config

Build target: Phase 3 (three-trait disease wiring).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Tuple

import numpy as np

from sswd_evoepi.config import DiseaseSection, PathogenEvolutionSection, SimulationConfig
from sswd_evoepi.types import DiseaseState, AGENT_DTYPE


# ═══════════════════════════════════════════════════════════════════════
# CONSTANTS
# ═══════════════════════════════════════════════════════════════════════

T_REF_K = 293.15       # 20°C in Kelvin — V. pectenicida T_opt (Lambert 1998)
T_OPT = 20.0           # °C — growth optimum
T_MAX = 30.0           # °C — upper lethal limit

# Erlang shape parameters for stage durations (disease-module-spec §5.2)
K_SHAPE_E = 3          # Incubation fairly regular (CV = 0.58)
K_SHAPE_I1 = 2         # More variable (CV = 0.71)
K_SHAPE_I2 = 2         # Variable end-stage (CV = 0.71)

# Size-dependent susceptibility (Eisenlord 2016: OR = 1.23 per 10mm)
BETA_L = 0.021         # per mm — ln(1.23)/10
L_BAR = 300.0          # mm — reference size
SIGMA_L = 100.0        # mm — normalization

# Salinity modifier exponent
ETA_SAL = 2.0

# VBNC dynamics
K_VBNC = 1.0           # °C⁻¹ — transition steepness

# Carcass shedding duration (days)
CARCASS_SHED_DAYS = 3

# Recovery thresholds
C_EARLY_THRESH = 0.5   # Minimum c_i for early recovery from I₁

# Behavioral modifiers indexed by DiseaseState value
SPEED_MODIFIER = np.array([1.0, 1.0, 0.5, 0.1, 0.0, 1.0], dtype=np.float64)
FEEDING_MODIFIER = np.array([1.0, 1.0, 0.5, 0.0, 0.0, 1.0], dtype=np.float64)
CAN_SPAWN = np.array([True, True, False, False, False, True])

# Vibrio decay parameters (log-linear interpolation between known points)
XI_10C = 1.0           # d⁻¹ at 10°C (half-life ~0.7 d)
XI_20C = 0.33          # d⁻¹ at 20°C (half-life ~2.1 d)


# ═══════════════════════════════════════════════════════════════════════
# TEMPERATURE FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════

def arrhenius(rate_ref: float, Ea_over_R: float, T_celsius: float) -> float:
    """Arrhenius temperature scaling.

    k(T) = k_ref × exp(E_a/R × (1/T_ref − 1/T))

    Args:
        rate_ref: Rate at T_ref (20°C).
        Ea_over_R: Activation energy / gas constant (K).
        T_celsius: Current temperature (°C).

    Returns:
        Rate at T_celsius.
    """
    T_K = T_celsius + 273.15
    return rate_ref * np.exp(Ea_over_R * (1.0 / T_REF_K - 1.0 / T_K))


def arrhenius_vec(rate_ref: float, Ea_over_R: float,
                  T_celsius: np.ndarray) -> np.ndarray:
    """Vectorized Arrhenius for arrays of temperatures."""
    T_K = T_celsius + 273.15
    return rate_ref * np.exp(Ea_over_R * (1.0 / T_REF_K - 1.0 / T_K))


def thermal_performance(
    Ea_over_R: float,
    T_celsius: float,
    rate_ref: float = 1.0,
) -> float:
    """Thermal performance curve with decline above T_opt.

    Below T_opt: standard Arrhenius increase.
    Above T_opt: quadratic decline to zero at T_max.
    """
    base = arrhenius(rate_ref, Ea_over_R, T_celsius)
    if T_celsius <= T_OPT:
        return base
    elif T_celsius >= T_MAX:
        return 0.0
    else:
        decline = 1.0 - ((T_celsius - T_OPT) / (T_MAX - T_OPT)) ** 2
        return base * max(0.0, decline)


# ═══════════════════════════════════════════════════════════════════════
# FORCE OF INFECTION COMPONENTS
# ═══════════════════════════════════════════════════════════════════════

def salinity_modifier(salinity: float, s_min: float = 10.0,
                      s_full: float = 28.0) -> float:
    """Salinity suppression of Vibrio viability. Returns [0, 1].

    At salinity > s_full: S_sal = 1 (full marine, no suppression).
    At salinity < s_min:  S_sal = 0 (Vibrio unviable in freshwater).
    """
    if salinity >= s_full:
        return 1.0
    if salinity <= s_min:
        return 0.0
    frac = (salinity - s_min) / (s_full - s_min)
    return frac ** ETA_SAL


def size_susceptibility(size_mm: float) -> float:
    """Size-dependent susceptibility multiplier. Larger = more susceptible.

    Eisenlord 2016: OR = 1.23 per 10mm radius.
    f_size(L) = exp(β_L × (L − L̄) / σ_L)
    """
    return np.exp(BETA_L * (size_mm - L_BAR) / SIGMA_L)


def force_of_infection(
    P_k: float,
    r_i: float,
    salinity: float,
    size_mm: float,
    cfg: DiseaseSection,
) -> float:
    """Per-individual instantaneous hazard rate of infection (d⁻¹).

    λ_i = a × P/(K_half+P) × (1−r_i) × S_sal × f_size(L_i)

    Args:
        P_k: Vibrio concentration at node (bacteria/mL).
        r_i: Individual resistance score [0, 1].
        salinity: Local salinity (psu).
        size_mm: Individual diameter (mm).
        cfg: Disease configuration section.

    Returns:
        Instantaneous hazard rate (d⁻¹).
    """
    if P_k <= 0.0:
        return 0.0
    dose_response = P_k / (cfg.K_half + P_k)
    resistance = 1.0 - r_i
    sal_mod = salinity_modifier(salinity, cfg.s_min, cfg.s_full)
    size_mod = size_susceptibility(size_mm)
    return cfg.a_exposure * dose_response * resistance * sal_mod * size_mod


def infection_probability(lambda_i: float, dt: float = 1.0) -> float:
    """Convert hazard rate to discrete probability over dt days.

    p = 1 − exp(−λ × Δt)
    """
    return 1.0 - np.exp(-lambda_i * dt)


# ═══════════════════════════════════════════════════════════════════════
# ENVIRONMENTAL PATHOGEN DYNAMICS
# ═══════════════════════════════════════════════════════════════════════

def vibrio_decay_rate(T_celsius: float) -> float:
    """Vibrio natural decay rate in seawater (d⁻¹).

    Vibrio survives LONGER at warm temperatures (Lupo 2020):
      10°C: ξ ≈ 1.0 d⁻¹ (half-life ~0.7 d)
      20°C: ξ ≈ 0.33 d⁻¹ (half-life ~2.1 d)

    Uses log-linear interpolation between known points,
    clamped outside the 10–20°C range.
    """
    if T_celsius <= 10.0:
        return XI_10C
    elif T_celsius >= 20.0:
        return XI_20C
    else:
        frac = (T_celsius - 10.0) / 10.0
        return np.exp(np.log(XI_10C) * (1.0 - frac) + np.log(XI_20C) * frac)


def environmental_vibrio(
    T_celsius: float,
    salinity: float,
    cfg: DiseaseSection,
) -> float:
    """Background environmental V. pectenicida input rate (bact/mL/d).

    Sigmoidal VBNC dynamics: near-zero below 10°C,
    activating 10–15°C, peak at 20°C, declining above.

    In invasion scenario, returns 0 (pathogen absent until arrival).
    """
    if cfg.scenario == "invasion":
        return 0.0

    # VBNC resuscitation sigmoid
    vbnc_activation = 1.0 / (1.0 + np.exp(-K_VBNC * (T_celsius - cfg.T_vbnc)))

    # Thermal performance (peak at T_opt=20°C)
    g_peak = thermal_performance(3000.0, T_celsius, rate_ref=1.0)

    sal_mod = salinity_modifier(salinity, cfg.s_min, cfg.s_full)

    return cfg.P_env_max * vbnc_activation * g_peak * sal_mod


def shedding_rate_I1(T_celsius: float, cfg: DiseaseSection) -> float:
    """Pre-symptomatic shedding rate (bact/mL/d/host) at temperature T."""
    return arrhenius(cfg.sigma_1_eff, cfg.Ea_sigma, T_celsius)


def shedding_rate_I2(T_celsius: float, cfg: DiseaseSection) -> float:
    """Symptomatic shedding rate (bact/mL/d/host) at temperature T."""
    return arrhenius(cfg.sigma_2_eff, cfg.Ea_sigma, T_celsius)


def update_vibrio_concentration(
    P_k: float,
    n_I1: int,
    n_I2: int,
    n_D_fresh: int,
    T_celsius: float,
    salinity: float,
    phi_k: float,
    dispersal_input: float,
    cfg: DiseaseSection,
    dt: float = 1.0,
    override_shedding: "float | None" = None,
) -> float:
    """Euler-step update of Vibrio concentration at one node.

    dP/dt = shedding − decay − flushing + background + dispersal

    Args:
        P_k: Current concentration (bacteria/mL).
        n_I1: Count of I₁ individuals at node.
        n_I2: Count of I₂ individuals at node.
        n_D_fresh: Count of fresh carcasses (< 3 days dead).
        T_celsius: Current SST (°C).
        salinity: Current salinity (psu).
        phi_k: Flushing rate (d⁻¹).
        dispersal_input: Σ d_jk × P_j from neighbor nodes.
        cfg: Disease configuration.
        dt: Timestep (days).
        override_shedding: If provided, use this value as total shedding
            instead of computing from n_I1/n_I2/n_D_fresh counts.  Used
            by pathogen evolution to inject per-individual strain-weighted
            shedding.

    Returns:
        Updated P_k (non-negative).
    """
    # Shedding from infectious hosts
    if override_shedding is not None:
        shed = override_shedding
    else:
        shed = (shedding_rate_I1(T_celsius, cfg) * n_I1
                + shedding_rate_I2(T_celsius, cfg) * n_I2
                + cfg.sigma_D * n_D_fresh)

    # Decay (faster at cold T)
    decay = vibrio_decay_rate(T_celsius) * P_k

    # Flushing
    flush = phi_k * P_k

    # Environmental reservoir
    env = environmental_vibrio(T_celsius, salinity, cfg)

    # Euler update
    dP = shed - decay - flush + env + dispersal_input
    P_new = max(0.0, P_k + dP * dt)

    return P_new


# ═══════════════════════════════════════════════════════════════════════
# DISEASE PROGRESSION
# ═══════════════════════════════════════════════════════════════════════

def sample_stage_duration(
    mu_rate: float,
    k_shape: int,
    rng: np.random.Generator,
) -> int:
    """Sample duration (days, integer) for a disease stage.

    Erlang(k, k×μ) has mean = 1/μ, variance = 1/(k×μ²).

    Args:
        mu_rate: Rate parameter at current T (d⁻¹).
        k_shape: Erlang shape parameter.
        rng: NumPy random generator.

    Returns:
        Duration in days (≥ 1).
    """
    if mu_rate <= 0.0:
        return 1
    rate = k_shape * mu_rate
    total = 0.0
    for _ in range(k_shape):
        total += rng.exponential(1.0 / rate)
    return max(1, int(np.round(total)))


def batch_sample_stage_duration(
    mu_rates: np.ndarray,
    k_shape: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Vectorized Erlang duration sampling for multiple individuals.

    Equivalent to calling sample_stage_duration() for each element,
    but uses batch Gamma sampling for ~10-100× speedup.

    Gamma(k, 1/(k×μ)) has mean 1/μ — same as Erlang(k, k×μ).

    Args:
        mu_rates: Array of rate parameters (d⁻¹), shape (n,).
        k_shape: Erlang shape parameter (same for all).
        rng: NumPy random generator.

    Returns:
        Array of durations in days (≥ 1), shape (n,), dtype int32.
    """
    n = len(mu_rates)
    if n == 0:
        return np.array([], dtype=np.int32)
    rates = np.maximum(np.asarray(mu_rates, dtype=np.float64), 1e-10)
    scales = 1.0 / (k_shape * rates)
    totals = rng.gamma(k_shape, scales, size=n)
    return np.maximum(1, np.round(totals)).astype(np.int32)


def recovery_probability_I2(c_i: float, rho_rec: float = 0.05) -> float:
    """Daily probability of recovery from I₂ state.

    p_rec = ρ_rec × c_i

    Linear in clearance ability. At c_i=0: no recovery possible.
    At c_i=1: p_rec = ρ_rec (base rate, ~5%/day).

    Args:
        c_i: Individual recovery/clearance score [0, 1].
        rho_rec: Base recovery rate (d⁻¹).
    """
    return rho_rec * c_i


def recovery_probability_I1(c_i: float, rho_rec: float = 0.05) -> float:
    """Daily probability of early recovery from I₁ state.

    Only meaningful for high-clearance individuals (c_i > 0.5).
    Linear above threshold:
        p_early = ρ_rec × 2 × (c_i − 0.5)

    At c_i ≤ 0.5: 0. At c_i=1.0: ρ_rec.

    Args:
        c_i: Individual recovery/clearance score [0, 1].
        rho_rec: Base recovery rate (d⁻¹).
    """
    if c_i <= C_EARLY_THRESH:
        return 0.0
    return rho_rec * 2.0 * (c_i - C_EARLY_THRESH)


# ═══════════════════════════════════════════════════════════════════════
# R₀ COMPUTATION
# ═══════════════════════════════════════════════════════════════════════

def compute_R0(
    T_celsius: float,
    S_0: int,
    phi_k: float,
    cfg: DiseaseSection,
    salinity: float = 30.0,
    mean_resistance: float = 0.08,
    v: "float | None" = None,
    pe_cfg: "PathogenEvolutionSection | None" = None,
    mean_tolerance: float = 0.0,
    mean_recovery: float = 0.0,
    tau_max: "float | None" = None,
) -> float:
    """Compute basic reproduction number at a node.

    R₀ = (a × S₀ × (1−r̄) × S_sal) / (K_half × (ξ(T) + φ))
         × (σ₁(T)/μ_I₁I₂(T) + σ₂(T)/μ_I₂D_eff(T) + σ_D × τ_D)

    Includes carcass saprophytic shedding contribution (CE-6).
    Tolerance reduces effective I₂→D rate (longer infectious period).
    Recovery adds an exit pathway from I₂ (reduces effective I₂ duration).

    When *v* and *pe_cfg* are provided and pathogen evolution is enabled,
    uses virulence-adjusted rates for the given strain.

    Args:
        T_celsius: Temperature (°C).
        S_0: Number of susceptible individuals.
        phi_k: Flushing rate (d⁻¹).
        cfg: Disease configuration.
        salinity: Local salinity (psu).
        mean_resistance: Population mean resistance score.
        v: Optional pathogen virulence phenotype for strain-specific R₀.
        pe_cfg: Optional PathogenEvolutionSection config.
        mean_tolerance: Population mean tolerance score [0, 1].
        mean_recovery: Population mean recovery/clearance score [0, 1].
        tau_max: Max mortality reduction. Defaults to cfg.tau_max.

    Returns:
        R₀ estimate.
    """
    if tau_max is None:
        tau_max = cfg.tau_max

    sal_mod = salinity_modifier(salinity, cfg.s_min, cfg.s_full)
    suscept = (1.0 - mean_resistance) * sal_mod

    xi = vibrio_decay_rate(T_celsius)
    total_removal = xi + phi_k

    if v is not None and pe_cfg is not None and pe_cfg.enabled:
        sigma1 = sigma_1_strain(v, T_celsius, cfg, pe_cfg)
        sigma2 = sigma_2_strain(v, T_celsius, cfg, pe_cfg)
        mu_i1i2 = mu_I1I2_strain(v, T_celsius, cfg, pe_cfg)
        mu_i2d = mu_I2D_strain(v, T_celsius, cfg, pe_cfg)
    else:
        sigma1 = arrhenius(cfg.sigma_1_eff, cfg.Ea_sigma, T_celsius)
        sigma2 = arrhenius(cfg.sigma_2_eff, cfg.Ea_sigma, T_celsius)
        mu_i1i2 = arrhenius(cfg.mu_I1I2_ref, cfg.Ea_I1I2, T_celsius)
        mu_i2d = arrhenius(cfg.mu_I2D_ref, cfg.Ea_I2D, T_celsius)

    # Effective I₂→D rate accounting for population-mean tolerance
    mu_i2d_eff = mu_i2d * (1.0 - mean_tolerance * tau_max)
    mu_i2d_eff = max(mu_i2d_eff, mu_i2d * 0.05)  # floor at 5%

    # Effective I₂ exit rate accounting for recovery
    # Mean time in I₂ ≈ 1 / (μ_I2D_eff + ρ_rec × mean_recovery)
    effective_I2_exit = mu_i2d_eff + cfg.rho_rec * mean_recovery

    # Total pathogen shed per infection event: live stages + carcass burst
    shedding_integral = (sigma1 / mu_i1i2
                         + sigma2 / effective_I2_exit
                         + cfg.sigma_D * CARCASS_SHED_DAYS)

    R0 = (cfg.a_exposure * S_0 * suscept) / (cfg.K_half * total_removal) * shedding_integral
    return R0


# ═══════════════════════════════════════════════════════════════════════
# PATHOGEN VIRULENCE TRADEOFF FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════

def sigma_1_strain(
    v: "float | np.ndarray",
    T_celsius: float,
    disease_cfg: DiseaseSection,
    pe_cfg: "PathogenEvolutionSection | None" = None,
) -> "float | np.ndarray":
    """Early (I₁) shedding rate for pathogen strain with virulence *v*.

    At v = v_anchor the result equals the base Arrhenius rate exactly.
    When pathogen evolution is disabled (pe_cfg is None or not enabled),
    returns the base rate regardless of *v*.

    Accepts both scalar and ndarray *v* for vectorized per-individual
    computation.
    """
    base = arrhenius(disease_cfg.sigma_1_eff, disease_cfg.Ea_sigma, T_celsius)
    if pe_cfg is None or not pe_cfg.enabled:
        return base
    dv = v - pe_cfg.v_anchor
    return base * np.exp(pe_cfg.alpha_shed * dv * pe_cfg.gamma_early)


def sigma_2_strain(
    v: "float | np.ndarray",
    T_celsius: float,
    disease_cfg: DiseaseSection,
    pe_cfg: "PathogenEvolutionSection | None" = None,
) -> "float | np.ndarray":
    """Late (I₂) shedding rate for pathogen strain with virulence *v*.

    At v = v_anchor the result equals the base Arrhenius rate exactly.
    """
    base = arrhenius(disease_cfg.sigma_2_eff, disease_cfg.Ea_sigma, T_celsius)
    if pe_cfg is None or not pe_cfg.enabled:
        return base
    dv = v - pe_cfg.v_anchor
    return base * np.exp(pe_cfg.alpha_shed * dv)


def mu_I2D_strain(
    v: "float | np.ndarray",
    T_celsius: float,
    disease_cfg: DiseaseSection,
    pe_cfg: "PathogenEvolutionSection | None" = None,
) -> "float | np.ndarray":
    """I₂→D death rate for pathogen strain with virulence *v*.

    At v = v_anchor the result equals the base Arrhenius rate exactly.
    """
    base = arrhenius(disease_cfg.mu_I2D_ref, disease_cfg.Ea_I2D, T_celsius)
    if pe_cfg is None or not pe_cfg.enabled:
        return base
    dv = v - pe_cfg.v_anchor
    return base * np.exp(pe_cfg.alpha_kill * dv)


def mu_I1I2_strain(
    v: "float | np.ndarray",
    T_celsius: float,
    disease_cfg: DiseaseSection,
    pe_cfg: "PathogenEvolutionSection | None" = None,
) -> "float | np.ndarray":
    """I₁→I₂ progression rate for pathogen strain with virulence *v*.

    At v = v_anchor the result equals the base Arrhenius rate exactly.
    """
    base = arrhenius(disease_cfg.mu_I1I2_ref, disease_cfg.Ea_I1I2, T_celsius)
    if pe_cfg is None or not pe_cfg.enabled:
        return base
    dv = v - pe_cfg.v_anchor
    return base * np.exp(pe_cfg.alpha_prog * dv)


# ═══════════════════════════════════════════════════════════════════════
# CARCASS TRACKER
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class CarcassTracker:
    """Track fresh carcasses at a node for saprophytic shedding.

    Uses a ring buffer of daily death counts over the last
    CARCASS_SHED_DAYS days.
    """
    daily_deaths: np.ndarray = field(default=None)
    write_idx: int = 0

    def __post_init__(self):
        if self.daily_deaths is None:
            self.daily_deaths = np.zeros(CARCASS_SHED_DAYS, dtype=np.int32)

    def add_deaths(self, n: int) -> None:
        """Record n disease deaths today."""
        self.daily_deaths[self.write_idx] = n
        self.write_idx = (self.write_idx + 1) % CARCASS_SHED_DAYS

    @property
    def n_fresh(self) -> int:
        """Total fresh carcasses (last 3 days)."""
        return int(self.daily_deaths.sum())

    def reset(self) -> None:
        """Clear all carcass records."""
        self.daily_deaths[:] = 0
        self.write_idx = 0


# ═══════════════════════════════════════════════════════════════════════
# NODE DISEASE STATE
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class NodeDiseaseState:
    """Disease state tracked per node. Updated daily."""

    node_id: int = 0

    # Environmental pathogen concentration (bacteria/mL)
    vibrio_concentration: float = 0.0

    # Daily compartment counts
    n_S: int = 0
    n_E: int = 0
    n_I1: int = 0
    n_I2: int = 0
    n_R: int = 0
    n_D_fresh: int = 0

    # Cumulative epidemic statistics
    cumulative_infections: int = 0
    cumulative_deaths: int = 0
    cumulative_recoveries: int = 0
    peak_prevalence: float = 0.0
    peak_vibrio: float = 0.0
    epidemic_active: bool = False
    epidemic_start_day: int = -1

    # R₀ estimate
    R0_estimate: float = 0.0

    # Virulence tracking accumulators (for pathogen evolution output)
    virulence_sum_new_infections: float = 0.0
    virulence_count_new_infections: int = 0
    virulence_sum_deaths: float = 0.0
    virulence_count_deaths: int = 0

    # Daily running mean of all infected agents' virulence (time-weighted)
    virulence_sum_daily: float = 0.0
    virulence_count_daily: int = 0

    # Carcass tracker
    carcass_tracker: CarcassTracker = field(default_factory=CarcassTracker)

    def reset_epidemic_stats(self) -> None:
        """Reset per-epidemic counters."""
        self.cumulative_infections = 0
        self.cumulative_deaths = 0
        self.cumulative_recoveries = 0
        self.peak_prevalence = 0.0
        self.peak_vibrio = 0.0
        self.epidemic_active = False
        self.epidemic_start_day = -1


# ═══════════════════════════════════════════════════════════════════════
# DAILY DISEASE UPDATE — CORE ENGINE
# ═══════════════════════════════════════════════════════════════════════

def daily_disease_update(
    agents: np.ndarray,
    node_state: NodeDiseaseState,
    T_celsius: float,
    salinity: float,
    phi_k: float,
    dispersal_input: float,
    day: int,
    cfg: DiseaseSection,
    rng: np.random.Generator,
    infected_density_grid=None,
    pe_cfg: "PathogenEvolutionSection | None" = None,
) -> NodeDiseaseState:
    """One daily timestep of disease dynamics at a single Tier 1 node.

    Sequence:
      1. Update Vibrio concentration (pathogen ODE)
      2. Transmission (S → E for susceptible individuals)
      3. Progression (E → I₁, I₁ → I₂, I₂ → D/R)
      4. Carcass tracking
      5. Update compartment counts & diagnostics

    Args:
        agents: Structured array of agents at this node.
        node_state: Mutable NodeDiseaseState for this node.
        T_celsius: Current SST (°C).
        salinity: Current salinity (psu).
        phi_k: Flushing rate (d⁻¹).
        dispersal_input: Pathogen input from neighbor nodes (bact/mL/d).
        day: Current simulation day.
        cfg: Disease configuration section.
        rng: NumPy random generator for this node.
        infected_density_grid: Optional InfectedDensityGrid (from
            movement module).  When provided, per-susceptible dose-
            response is scaled by local infected density (spatial
            transmission).  When None, mean-field (well-mixed) is used.

    Returns:
        Updated NodeDiseaseState.
    """
    alive = agents['alive']
    ds = agents['disease_state']
    dt_rem = agents['disease_timer']
    resistance = agents['resistance']
    size = agents['size']

    new_infections = 0
    new_deaths = 0
    new_recoveries = 0

    # ── STEP 1: Update Vibrio concentration ──────────────────────────

    # Count current shedding individuals using bincount (single pass)
    alive_mask = alive.astype(bool)
    alive_ds = ds[alive_mask]
    if len(alive_ds) > 0:
        _counts_init = np.bincount(alive_ds, minlength=6)
        n_I1 = int(_counts_init[DiseaseState.I1])
        n_I2 = int(_counts_init[DiseaseState.I2])
    else:
        n_I1 = 0
        n_I2 = 0
    n_D_fresh = node_state.carcass_tracker.n_fresh

    # Per-individual strain-weighted shedding when pathogen evolution active
    pe_active = pe_cfg is not None and pe_cfg.enabled
    override_shed = None
    if pe_active:
        I1_mask = alive_mask & (ds == DiseaseState.I1)
        I2_mask = alive_mask & (ds == DiseaseState.I2)
        I1_indices = np.where(I1_mask)[0]
        I2_indices = np.where(I2_mask)[0]

        shed_I1 = 0.0
        if len(I1_indices) > 0:
            v_I1 = agents['pathogen_virulence'][I1_indices]
            shed_I1 = float(np.sum(sigma_1_strain(v_I1, T_celsius, cfg, pe_cfg)))

        shed_I2 = 0.0
        if len(I2_indices) > 0:
            v_I2 = agents['pathogen_virulence'][I2_indices]
            shed_I2 = float(np.sum(sigma_2_strain(v_I2, T_celsius, cfg, pe_cfg)))

        shed_D = cfg.sigma_D * n_D_fresh  # carcass shedding unchanged
        override_shed = shed_I1 + shed_I2 + shed_D

    P_k = update_vibrio_concentration(
        node_state.vibrio_concentration,
        n_I1, n_I2, n_D_fresh,
        T_celsius, salinity, phi_k, dispersal_input,
        cfg,
        override_shedding=override_shed,
    )
    node_state.vibrio_concentration = P_k

    # ── STEP 2: Transmission (S → E) ────────────────────────────────

    # Pre-compute shared values for efficiency
    if P_k > 0:
        sal_mod = salinity_modifier(salinity, cfg.s_min, cfg.s_full)
        mu_EI1 = arrhenius(cfg.mu_EI1_ref, cfg.Ea_EI1, T_celsius)

        susceptible = alive_mask & (ds == DiseaseState.S)

        # Phase 11: Juvenile immunity — exclude recently settled agents
        if cfg.min_susceptible_age_days > 0:
            age_days = day - agents['settlement_day']
            juvenile_mask = age_days < cfg.min_susceptible_age_days
            susceptible = susceptible & ~juvenile_mask

        susc_indices = np.where(susceptible)[0]

        if len(susc_indices) > 0:
            # Vectorized force of infection with Phase 4 immunosuppression
            # r_eff = r_i / psi_spawn if immunosuppression_timer > 0, else r_i
            r_raw = resistance[susc_indices]
            immuno_timers = agents['immunosuppression_timer'][susc_indices]
            
            # Apply immunosuppression: divide resistance by susceptibility multiplier
            if cfg.immunosuppression_enabled:
                immuno_mask = immuno_timers > 0
                r_eff = r_raw.copy()
                r_eff[immuno_mask] = r_raw[immuno_mask] / cfg.susceptibility_multiplier
                # Clamp effective resistance to [0, 1] 
                r_eff = np.clip(r_eff, 0.0, 1.0)
            else:
                r_eff = r_raw
            
            r_mod = 1.0 - r_eff
            s_mod = np.exp(BETA_L * (size[susc_indices] - L_BAR) / SIGMA_L)

            # Dose-response: spatial or mean-field
            if infected_density_grid is not None:
                # Spatial: scale P_k by local infected density factor
                local_factors = infected_density_grid.lookup(
                    agents['x'][susc_indices],
                    agents['y'][susc_indices],
                )
                P_local = P_k * local_factors
                dose_resp = P_local / (cfg.K_half + P_local)
            else:
                # Mean-field: uniform P_k for all susceptibles
                dose_resp = P_k / (cfg.K_half + P_k)

            lambda_arr = cfg.a_exposure * dose_resp * r_mod * sal_mod * s_mod
            p_inf = 1.0 - np.exp(-lambda_arr)

            # Stochastic infection draws
            draws = rng.random(len(susc_indices))
            infected_mask = draws < p_inf

            # ── Pre-compute strain inheritance data (once) ───────
            if pe_active:
                alive_I1_mask = alive_mask & (ds == DiseaseState.I1)
                alive_I2_mask = alive_mask & (ds == DiseaseState.I2)

                # Total shedding from local hosts vs environment
                P_shed_total = (shedding_rate_I1(T_celsius, cfg) * int(np.sum(alive_I1_mask))
                                + shedding_rate_I2(T_celsius, cfg) * int(np.sum(alive_I2_mask)))
                P_env_total = environmental_vibrio(T_celsius, salinity, cfg)
                P_total = P_shed_total + P_env_total

                # Pre-compute shedder sampling distribution
                shedder_mask = alive_I1_mask | alive_I2_mask
                shedder_indices = np.where(shedder_mask)[0]
                if len(shedder_indices) > 0 and P_total > 0:
                    p_from_shedder = P_shed_total / P_total
                    shedder_v = agents['pathogen_virulence'][shedder_indices]
                    weights = np.where(
                        ds[shedder_indices] == DiseaseState.I2,
                        sigma_2_strain(shedder_v, T_celsius, cfg, pe_cfg),
                        sigma_1_strain(shedder_v, T_celsius, cfg, pe_cfg),
                    )
                    w_sum = weights.sum()
                    if w_sum > 0:
                        weights = weights / w_sum
                    else:
                        weights = np.ones(len(shedder_indices)) / len(shedder_indices)
                else:
                    p_from_shedder = 0.0

            new_inf_idx = susc_indices[infected_mask]
            n_new = len(new_inf_idx)
            if n_new > 0:
                ds[new_inf_idx] = DiseaseState.E
                # Batch Erlang timer sampling
                new_timers = batch_sample_stage_duration(
                    np.full(n_new, mu_EI1), K_SHAPE_E, rng,
                )
                dt_rem[new_inf_idx] = new_timers
                new_infections += n_new

                # Strain inheritance (batch)
                if pe_active:
                    inherit_draws = rng.random(n_new)
                    from_shedder = (
                        (P_total > 0)
                        & (len(shedder_indices) > 0)
                        & (inherit_draws < p_from_shedder)
                    )
                    v_parents = np.full(n_new, pe_cfg.v_env)
                    n_from_shed = int(np.sum(from_shedder))
                    if n_from_shed > 0:
                        sources = rng.choice(
                            shedder_indices, size=n_from_shed, p=weights,
                        )
                        v_parents[from_shedder] = agents['pathogen_virulence'][sources]
                    v_new = v_parents + rng.normal(0, pe_cfg.sigma_v_mutation, size=n_new)
                    v_new = np.clip(v_new, pe_cfg.v_min, pe_cfg.v_max)
                    agents['pathogen_virulence'][new_inf_idx] = v_new
                    node_state.virulence_sum_new_infections += float(np.sum(v_new))
                    node_state.virulence_count_new_infections += n_new

    # ── STEP 3: Disease progression (vectorized) ──────────────────────

    # Pre-compute temperature-dependent base rates
    mu_I1I2 = arrhenius(cfg.mu_I1I2_ref, cfg.Ea_I1I2, T_celsius)
    mu_I2D = arrhenius(cfg.mu_I2D_ref, cfg.Ea_I2D, T_celsius)

    # All diseased individuals (E, I1, I2)
    diseased = alive_mask & (ds >= DiseaseState.E) & (ds <= DiseaseState.I2)
    diseased_idx = np.where(diseased)[0]

    if len(diseased_idx) > 0:
        # Batch timer decrement
        dt_rem[diseased_idx] -= 1

        states = ds[diseased_idx]
        timers = dt_rem[diseased_idx]

        # Split: expired (transition) vs not-expired (recovery check)
        expired_mask = timers <= 0
        not_expired_mask = ~expired_mask

        # ═══ EXPIRED TIMERS: state transitions ═══════════════════════

        exp_idx = diseased_idx[expired_mask]
        exp_states = states[expired_mask]

        # ── E → I1 ──────────────────────────────────────────────────
        e_mask = exp_states == DiseaseState.E
        e_idx = exp_idx[e_mask]
        if len(e_idx) > 0:
            ds[e_idx] = DiseaseState.I1
            if pe_active:
                v_arr = agents['pathogen_virulence'][e_idx]
                rates = mu_I1I2_strain(v_arr, T_celsius, cfg, pe_cfg)
            else:
                rates = np.full(len(e_idx), mu_I1I2)
            dt_rem[e_idx] = batch_sample_stage_duration(rates, K_SHAPE_I1, rng)

        # ── I1 → I2 ─────────────────────────────────────────────────
        i1_mask = exp_states == DiseaseState.I1
        i1_idx = exp_idx[i1_mask]
        if len(i1_idx) > 0:
            ds[i1_idx] = DiseaseState.I2
            if pe_active:
                v_arr = agents['pathogen_virulence'][i1_idx]
                rates = mu_I2D_strain(v_arr, T_celsius, cfg, pe_cfg)
            else:
                rates = np.full(len(i1_idx), mu_I2D)
            # Tolerance → longer I₂ timer
            t_arr = agents['tolerance'][i1_idx].astype(np.float64)
            effective_rates = rates * (1.0 - t_arr * cfg.tau_max)
            effective_rates = np.maximum(effective_rates, rates * 0.05)
            dt_rem[i1_idx] = batch_sample_stage_duration(
                effective_rates, K_SHAPE_I2, rng,
            )

        # ── I2 → D (death) ──────────────────────────────────────────
        i2d_mask = exp_states == DiseaseState.I2
        i2d_idx = exp_idx[i2d_mask]
        if len(i2d_idx) > 0:
            if pe_active:
                v_deaths = agents['pathogen_virulence'][i2d_idx]
                node_state.virulence_sum_deaths += float(np.sum(v_deaths))
                node_state.virulence_count_deaths += len(i2d_idx)
            ds[i2d_idx] = DiseaseState.D
            agents['alive'][i2d_idx] = False
            agents['cause_of_death'][i2d_idx] = 1  # DeathCause.DISEASE
            new_deaths += len(i2d_idx)

        # ═══ NOT-EXPIRED TIMERS: recovery checks ════════════════════

        ne_idx = diseased_idx[not_expired_mask]
        ne_states = states[not_expired_mask]

        # ── I2 recovery ──────────────────────────────────────────────
        i2_ne_mask = ne_states == DiseaseState.I2
        i2_ne_idx = ne_idx[i2_ne_mask]
        if len(i2_ne_idx) > 0:
            c_arr = agents['recovery_ability'][i2_ne_idx].astype(np.float64)
            p_rec = cfg.rho_rec * c_arr
            draws = rng.random(len(i2_ne_idx))
            rec_mask = draws < p_rec
            rec_idx = i2_ne_idx[rec_mask]
            if len(rec_idx) > 0:
                ds[rec_idx] = DiseaseState.R
                dt_rem[rec_idx] = 0
                agents['pathogen_virulence'][rec_idx] = 0.0
                new_recoveries += len(rec_idx)

        # ── I1 early recovery ────────────────────────────────────────
        i1_ne_mask = ne_states == DiseaseState.I1
        i1_ne_idx = ne_idx[i1_ne_mask]
        if len(i1_ne_idx) > 0:
            c_arr = agents['recovery_ability'][i1_ne_idx].astype(np.float64)
            above = c_arr > C_EARLY_THRESH
            i1_cand = i1_ne_idx[above]
            if len(i1_cand) > 0:
                c_cand = c_arr[above]
                p_early = cfg.rho_rec * 2.0 * (c_cand - C_EARLY_THRESH)
                draws = rng.random(len(i1_cand))
                rec_mask = draws < p_early
                rec_idx = i1_cand[rec_mask]
                if len(rec_idx) > 0:
                    ds[rec_idx] = DiseaseState.R
                    dt_rem[rec_idx] = 0
                    agents['pathogen_virulence'][rec_idx] = 0.0
                    new_recoveries += len(rec_idx)

    # ── STEP 4: Carcass tracking ─────────────────────────────────────

    node_state.carcass_tracker.add_deaths(new_deaths)

    # ── STEP 5: Update compartment counts & diagnostics ──────────────

    alive_mask = agents['alive'].astype(bool)  # refresh after deaths
    # Use bincount for single-pass compartment counting (replaces 5 separate scans)
    alive_ds = ds[alive_mask]
    if len(alive_ds) > 0:
        counts = np.bincount(alive_ds, minlength=6)
        node_state.n_S = int(counts[DiseaseState.S])
        node_state.n_E = int(counts[DiseaseState.E])
        node_state.n_I1 = int(counts[DiseaseState.I1])
        node_state.n_I2 = int(counts[DiseaseState.I2])
        node_state.n_R = int(counts[DiseaseState.R])
    else:
        node_state.n_S = 0
        node_state.n_E = 0
        node_state.n_I1 = 0
        node_state.n_I2 = 0
        node_state.n_R = 0
    node_state.n_D_fresh = node_state.carcass_tracker.n_fresh

    node_state.cumulative_infections += new_infections
    node_state.cumulative_deaths += new_deaths
    node_state.cumulative_recoveries += new_recoveries

    # Accumulate daily virulence of all infected agents (time-weighted mean)
    n_inf_total = node_state.n_E + node_state.n_I1 + node_state.n_I2
    if n_inf_total > 0 and pe_cfg is not None:
        inf_mask = alive_mask & (
            (ds == DiseaseState.E) | (ds == DiseaseState.I1) |
            (ds == DiseaseState.I2)
        )
        node_state.virulence_sum_daily += float(
            np.sum(agents['pathogen_virulence'][inf_mask])
        )
        node_state.virulence_count_daily += n_inf_total

    n_alive = int(np.sum(alive_mask))
    if n_alive > 0:
        prevalence = (node_state.n_I1 + node_state.n_I2) / n_alive
        node_state.peak_prevalence = max(node_state.peak_prevalence, prevalence)
    node_state.peak_vibrio = max(node_state.peak_vibrio, P_k)

    # Epidemic detection
    if new_infections > 0 and not node_state.epidemic_active:
        node_state.epidemic_active = True
        node_state.epidemic_start_day = day
    elif (node_state.epidemic_active
          and node_state.n_I1 + node_state.n_I2 + node_state.n_E == 0):
        node_state.epidemic_active = False

    return node_state


# ═══════════════════════════════════════════════════════════════════════
# EPIDEMIC SIMULATION (standalone, for testing / single-node runs)
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class EpidemicResult:
    """Results from a single-node epidemic simulation."""
    days: int = 0
    total_infected: int = 0
    total_deaths: int = 0
    total_recoveries: int = 0
    peak_prevalence: float = 0.0
    peak_vibrio: float = 0.0
    mortality_fraction: float = 0.0
    final_pop: int = 0
    initial_pop: int = 0
    daily_S: Optional[np.ndarray] = None
    daily_I: Optional[np.ndarray] = None
    daily_P: Optional[np.ndarray] = None
    daily_D_cumul: Optional[np.ndarray] = None


def run_single_node_epidemic(
    n_individuals: int,
    T_celsius: float,
    salinity: float,
    phi_k: float,
    cfg: DiseaseSection,
    n_days: int = 365,
    initial_infected: int = 5,
    mean_resistance: float = 0.08,
    resistance_std: float = 0.04,
    mean_size: float = 500.0,
    size_std: float = 100.0,
    seed: int = 42,
    record_daily: bool = False,
    initial_vibrio: Optional[float] = None,
    pe_cfg: "PathogenEvolutionSection | None" = None,
    mean_tolerance: float = 0.10,
    tolerance_std: float = 0.05,
    mean_recovery: float = 0.08,
    recovery_std: float = 0.04,
) -> EpidemicResult:
    """Run a standalone single-node epidemic simulation.

    Creates a population of susceptible individuals with given
    resistance/tolerance/recovery distributions, seeds initial
    infections, and runs daily disease updates.

    Args:
        n_individuals: Population size.
        T_celsius: Constant SST (°C).
        salinity: Constant salinity (psu).
        phi_k: Flushing rate (d⁻¹).
        cfg: Disease configuration.
        n_days: Number of days to simulate.
        initial_infected: Number of initially infected individuals.
        mean_resistance: Mean of resistance distribution.
        resistance_std: SD of resistance distribution.
        mean_size: Mean body size (mm).
        size_std: SD of body size (mm).
        seed: RNG seed.
        record_daily: If True, record daily compartment timeseries.
        initial_vibrio: Starting Vibrio concentration. If None, computed
            from steady-state background.
        pe_cfg: Optional pathogen evolution config.
        mean_tolerance: Mean of tolerance distribution.
        tolerance_std: SD of tolerance distribution.
        mean_recovery: Mean of recovery/clearance distribution.
        recovery_std: SD of recovery/clearance distribution.

    Returns:
        EpidemicResult with summary statistics.
    """
    rng = np.random.default_rng(seed)

    # Allocate agents
    agents = np.zeros(n_individuals, dtype=AGENT_DTYPE)
    agents['alive'] = True
    agents['disease_state'] = DiseaseState.S
    agents['disease_timer'] = 0

    # Assign resistance (truncated normal, [0, 1])
    r_values = rng.normal(mean_resistance, resistance_std, n_individuals)
    r_values = np.clip(r_values, 0.0, 1.0)
    agents['resistance'] = r_values.astype(np.float32)

    # Assign tolerance (truncated normal, [0, 1])
    if tolerance_std > 0:
        t_values = rng.normal(mean_tolerance, tolerance_std, n_individuals)
        t_values = np.clip(t_values, 0.0, 1.0)
        agents['tolerance'] = t_values.astype(np.float32)
    else:
        agents['tolerance'] = np.float32(mean_tolerance)

    # Assign recovery ability (truncated normal, [0, 1])
    if recovery_std > 0:
        c_values = rng.normal(mean_recovery, recovery_std, n_individuals)
        c_values = np.clip(c_values, 0.0, 1.0)
        agents['recovery_ability'] = c_values.astype(np.float32)
    else:
        agents['recovery_ability'] = np.float32(mean_recovery)

    # Assign sizes (truncated normal, min 50mm)
    sizes = rng.normal(mean_size, size_std, n_individuals)
    sizes = np.clip(sizes, 50.0, 1000.0)
    agents['size'] = sizes.astype(np.float32)

    # Seed initial infections
    if initial_infected > 0:
        infect_idx = rng.choice(n_individuals, size=min(initial_infected, n_individuals),
                                replace=False)
        mu_EI1 = arrhenius(cfg.mu_EI1_ref, cfg.Ea_EI1, T_celsius)
        pe_active = pe_cfg is not None and pe_cfg.enabled
        for idx in infect_idx:
            agents['disease_state'][idx] = DiseaseState.E
            agents['disease_timer'][idx] = sample_stage_duration(
                mu_EI1, K_SHAPE_E, rng
            )
            if pe_active:
                agents['pathogen_virulence'][idx] = pe_cfg.v_init

    # Initialize Vibrio
    if initial_vibrio is not None:
        P_init = initial_vibrio
    else:
        env = environmental_vibrio(T_celsius, salinity, cfg)
        xi = vibrio_decay_rate(T_celsius)
        P_init = env / (xi + phi_k) if (xi + phi_k) > 0 else 0.0

    node_state = NodeDiseaseState(node_id=0, vibrio_concentration=P_init)

    # Daily recording
    if record_daily:
        daily_S = np.zeros(n_days, dtype=np.int32)
        daily_I = np.zeros(n_days, dtype=np.int32)
        daily_P = np.zeros(n_days, dtype=np.float64)
        daily_D_cumul = np.zeros(n_days, dtype=np.int32)
    else:
        daily_S = daily_I = daily_P = daily_D_cumul = None

    # Run simulation
    for day in range(n_days):
        node_state = daily_disease_update(
            agents, node_state,
            T_celsius, salinity, phi_k,
            dispersal_input=0.0,
            day=day, cfg=cfg, rng=rng,
            pe_cfg=pe_cfg,
        )

        if record_daily:
            daily_S[day] = node_state.n_S
            daily_I[day] = node_state.n_I1 + node_state.n_I2
            daily_P[day] = node_state.vibrio_concentration
            daily_D_cumul[day] = node_state.cumulative_deaths

    # Compile results
    final_alive = int(np.sum(agents['alive']))
    result = EpidemicResult(
        days=n_days,
        total_infected=node_state.cumulative_infections,
        total_deaths=node_state.cumulative_deaths,
        total_recoveries=node_state.cumulative_recoveries,
        peak_prevalence=node_state.peak_prevalence,
        peak_vibrio=node_state.peak_vibrio,
        mortality_fraction=(node_state.cumulative_deaths / n_individuals
                            if n_individuals > 0 else 0.0),
        final_pop=final_alive,
        initial_pop=n_individuals,
        daily_S=daily_S,
        daily_I=daily_I,
        daily_P=daily_P,
        daily_D_cumul=daily_D_cumul,
    )
    return result


# ═══════════════════════════════════════════════════════════════════════
# INVASION SCENARIO HELPERS
# ═══════════════════════════════════════════════════════════════════════

def initialize_pathogen_ubiquitous(
    T_celsius: float,
    salinity: float,
    phi_k: float,
    cfg: DiseaseSection,
) -> float:
    """Compute steady-state Vibrio at a node (ubiquitous scenario).

    P* = P_env(T) / (ξ(T) + φ_k)
    """
    env = environmental_vibrio(T_celsius, salinity, cfg)
    xi = vibrio_decay_rate(T_celsius)
    denom = xi + phi_k
    return env / denom if denom > 0 else 0.0


def initialize_pathogen_invasion(
    node_id: int,
    cfg: DiseaseSection,
    invasion_dose: float = 10000.0,
) -> float:
    """Initialize Vibrio for invasion scenario.

    Returns invasion_dose if node_id is an invasion node, else 0.
    """
    if cfg.invasion_nodes and node_id in cfg.invasion_nodes:
        return invasion_dose
    return 0.0


# ═══════════════════════════════════════════════════════════════════════
# BEHAVIORAL MODIFIER QUERIES
# ═══════════════════════════════════════════════════════════════════════

def get_speed_modifier(disease_state: int) -> float:
    """Get movement speed modifier for a disease state."""
    return float(SPEED_MODIFIER[disease_state])


def get_feeding_modifier(disease_state: int) -> float:
    """Get feeding rate modifier for a disease state."""
    return float(FEEDING_MODIFIER[disease_state])


def can_spawn(disease_state: int) -> bool:
    """Check if individual in given disease state can spawn."""
    return bool(CAN_SPAWN[disease_state])
