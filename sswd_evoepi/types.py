"""Core data types for SSWD-EvoEpi.

This module is the SINGLE SOURCE OF TRUTH for:
  - AGENT_DTYPE: NumPy structured array dtype for individual agents
  - Stage, DiseaseState, Origin, Tier enumerations
  - Genotype constants (N_LOCI, N_ADDITIVE, IDX_EF1A)
  - Inter-module data transfer objects (LarvalCohort, SettlerPacket, NodeSnapshot)

All modules import these types from here. No other module defines agent fields.

References:
  - integration-architecture-spec.md §4 (canonical AGENT_DTYPE)
  - population-dynamics-spec.md §1 (Stage enum)
  - disease-module-spec.md (DiseaseState enum)
  - genetics-evolution-spec.md §1 (genotype constants)
  - ERRATA E4: disease_timer uses COUNTDOWN semantics
  - ERRATA E13: origin field for release tracking
  - CODE_ERRATA CE-1: cost_resistance removed; fecundity_mod always 1.0
"""

from dataclasses import dataclass, field
from enum import IntEnum
from typing import Optional

import numpy as np


# ═══════════════════════════════════════════════════════════════════════
# ENUMERATIONS
# ═══════════════════════════════════════════════════════════════════════

class Stage(IntEnum):
    """Life stages for Pycnopodia helianthoides.

    Transition criteria (size thresholds in mm diameter):
      EGG_LARVA  →  SETTLER:    settlement event (larval module)
      SETTLER    →  JUVENILE:   size ≥ 10 mm
      JUVENILE   →  SUBADULT:   size ≥ 150 mm
      SUBADULT   →  ADULT:      size ≥ 400 mm (reproductive maturity)
    """
    EGG_LARVA = 0   # Planktonic; 0–0.5 mm; 49–146 days PLD
    SETTLER   = 1   # Benthic recruit; 0.5–10 mm; ~1 year
    JUVENILE  = 2   # 10–150 mm; ~1–5 years
    SUBADULT  = 3   # 150–400 mm; ~5–10 years
    ADULT     = 4   # >400 mm; 10+ years; reproductive


class DiseaseState(IntEnum):
    """SEIPD+R compartments for SSWD.

    S  → E  (exposure via force of infection)
    E  → I1 (latent period, Erlang-distributed)
    I1 → I2 (early → late infectious)
    I2 → D  (death from disease) or I2 → R (recovery)
    """
    S  = 0   # Susceptible
    E  = 1   # Exposed (latent, not shedding)
    I1 = 2   # Early infectious (pre-symptomatic shedding)
    I2 = 3   # Late infectious (wasting, high shedding)
    D  = 4   # Dead from disease
    R  = 5   # Recovered (immune)


class DeathCause(IntEnum):
    """Cause of death tracking for demographic analysis."""
    ALIVE           = 0   # Not dead (default)
    DISEASE         = 1   # Killed by SSWD (I2 → D transition)
    NATURAL         = 2   # Stage-specific annual mortality
    SENESCENCE      = 3   # Age > senescence_age
    STARVATION      = 4   # Density-dependent (future use)


class Origin(IntEnum):
    """Origin tracking for conservation analysis (ERRATA E13)."""
    WILD            = 0   # Wild-born
    CAPTIVE_BRED    = 1   # Captive breeding program
    AGF_CROSS       = 2   # Assisted gene flow cross
    WILD_SOURCE_REL = 3   # Wild-sourced release


class Tier(IntEnum):
    """Two-tier fidelity system for spatial nodes."""
    TIER1 = 1   # Full ABM with individual agents
    TIER2 = 2   # Simplified demographics with genotype bank


# ═══════════════════════════════════════════════════════════════════════
# STAGE TRANSITION THRESHOLDS
# ═══════════════════════════════════════════════════════════════════════

# Size thresholds (mm) for stage transitions — one-directional only
STAGE_SIZE_THRESHOLDS = {
    Stage.SETTLER:  10.0,    # → JUVENILE
    Stage.JUVENILE: 150.0,   # → SUBADULT
    Stage.SUBADULT: 400.0,   # → ADULT
}


# ═══════════════════════════════════════════════════════════════════════
# GENOTYPE CONSTANTS
# ═══════════════════════════════════════════════════════════════════════

N_LOCI = 52          # 51 additive + 1 overdominant (EF1A analog)
N_ADDITIVE = 51      # Additive resistance loci: indices 0–50
IDX_EF1A = 51        # Overdominant locus index (EF1A analog)


# ═══════════════════════════════════════════════════════════════════════
# AGENT_DTYPE — Canonical structured array for individual agents
# ═══════════════════════════════════════════════════════════════════════

AGENT_DTYPE = np.dtype([
    # --- Spatial (POP writes) ---
    ('x',              np.float32),   #  4 B — position X within node (m)
    ('y',              np.float32),   #  4 B — position Y within node (m)
    ('heading',        np.float32),   #  4 B — movement heading (radians)
    ('speed',          np.float32),   #  4 B — instantaneous speed (m/min)

    # --- Life history (POP writes) ---
    ('size',           np.float32),   #  4 B — arm-tip diameter (mm)
    ('age',            np.float32),   #  4 B — age in years (fractional)
    ('stage',          np.int8),      #  1 B — Stage enum (0=EGG_LARVA..4=ADULT)
    ('sex',            np.int8),      #  1 B — 0=female, 1=male

    # --- Disease (DIS writes) ---
    ('disease_state',  np.int8),      #  1 B — DiseaseState enum (0=S..5=R)
    ('disease_timer',  np.int16),     #  2 B — days REMAINING in current disease state
                                      #         COUNTDOWN semantics (ERRATA E4)
                                      #         Managed exclusively by disease module.

    # --- Genetics (GEN writes) ---
    ('resistance',     np.float32),   #  4 B — polygenic resistance r_i ∈ [0, 1]
    ('fecundity_mod',  np.float32),   #  4 B — fecundity modifier; always 1.0
                                      #         (CE-1: cost of resistance removed)

    # --- Spawning (SPAWN writes) ---
    ('spawning_ready',          np.int8),   #  1 B — 0=not ready, 1=ready to spawn this season
    ('has_spawned',             np.int8),   #  1 B — females: 0/1, males: count of bouts (0-3)
    ('spawn_refractory',        np.int16),  #  2 B — days remaining in male refractory (countdown)
    ('spawn_gravity_timer',     np.int16),  #  2 B — days remaining in gravity phase (countdown)
    ('immunosuppression_timer', np.int16),  #  2 B — days remaining in post-spawning immunosuppression
    ('last_spawn_day',          np.int16),  #  2 B — day-of-year of last spawning event

    # --- Administrative ---
    ('node_id',        np.int16),     #  2 B — home node index
    ('alive',          np.bool_),     #  1 B — active flag (POP + DIS can set False)
    ('origin',         np.int8),      #  1 B — Origin enum (ERRATA E13)
                                      #         0=wild, 1=captive, 2=agf, 3=wild_source
    ('cause_of_death', np.int8),      #  1 B — DeathCause enum (0=alive, 1=disease,
                                      #         2=natural, 3=senescence)
    ('pathogen_virulence', np.float32),  #  4 B — virulence of infecting strain (0 when S or R)
    ('settlement_day', np.int32),     #  4 B — absolute sim day when settled/initialized
                                      #         0 for initial pop (always susceptible)
                                      #         Used for juvenile immunity (Phase 11)
])
# Total: ~55 bytes per agent (was ~51, +4 bytes for settlement_day)


def allocate_agents(max_n: int) -> np.ndarray:
    """Allocate a zeroed agent array.

    Args:
        max_n: Maximum number of agents (array capacity).

    Returns:
        Zeroed structured array of shape (max_n,) with AGENT_DTYPE.
    """
    return np.zeros(max_n, dtype=AGENT_DTYPE)


def allocate_genotypes(max_n: int) -> np.ndarray:
    """Allocate a zeroed genotype array.

    Stored separately from AGENT_DTYPE for cache-line efficiency during
    non-genetic operations (disease transmission, movement).

    Args:
        max_n: Maximum number of agents.

    Returns:
        Zeroed array of shape (max_n, N_LOCI, 2), dtype int8.
        Axis 0: agents, Axis 1: loci, Axis 2: allele copies (diploid).
    """
    return np.zeros((max_n, N_LOCI, 2), dtype=np.int8)


# ═══════════════════════════════════════════════════════════════════════
# INTER-MODULE DATA TRANSFER OBJECTS
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class LarvalCohort:
    """Produced by POP+GEN at source node; consumed by SPATIAL for dispersal."""
    source_node: int
    n_competent: int
    genotypes: np.ndarray          # (n_competent, N_LOCI, 2) int8
    parent_pairs: np.ndarray       # (n_competent, 2) int32 — mother, father indices
    pld_days: float
    spawn_day: int = 0             # absolute simulation day when spawned
    sst_at_spawn: float = 10.5    # SST at spawning (determines PLD)

    @property
    def settlement_day(self) -> int:
        """Absolute sim day when this cohort becomes competent to settle."""
        return self.spawn_day + int(self.pld_days)


@dataclass
class SettlerPacket:
    """Produced by SPATIAL dispersal; consumed by POP+GEN at receiving node."""
    dest_node: int
    genotypes: np.ndarray          # (n_settlers, N_LOCI, 2) int8
    source_nodes: np.ndarray       # (n_settlers,) int16 — origin node for each settler


@dataclass
class NodeSnapshot:
    """Read-only summary of a node's state for cross-module queries."""
    node_id: int
    tier: int
    n_alive: int
    n_adults: int
    mean_resistance: float
    var_resistance: float
    vibrio_concentration: float
    sst: float
    salinity: float
    flushing_rate: float
    R0_estimate: float
    epidemic_active: bool


# ═══════════════════════════════════════════════════════════════════════
# ANNUAL SURVIVAL (by stage)
# ═══════════════════════════════════════════════════════════════════════

# Annual survival probabilities indexed by Stage enum value
# Stage 0 (EGG_LARVA) is a placeholder — larvae are LarvalCohorts, not agents
ANNUAL_SURVIVAL = np.array([
    0.001,   # EGG_LARVA (placeholder; handled by larval module)
    0.03,    # SETTLER: 3% annual survival (Sewell & Watson / Miner 2018)
    0.90,    # JUVENILE
    0.95,    # SUBADULT
    0.98,    # ADULT
], dtype=np.float64)
