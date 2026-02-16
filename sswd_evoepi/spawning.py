"""Spawning system for SSWD-EvoEpi.

Phase 1: Extended spawning season with spontaneous spawning.
Future phases will add cascade induction, spawning gravity, and immunosuppression.

References:
  - spawning-overhaul-spec.md (complete specification)
  - reproduction-srs-spec.md (SRS lottery implementation)
  - population-dynamics-spec.md (fecundity calculations)

Author: Anton 🔬
"""

from typing import List, Tuple, Optional
import numpy as np
from scipy.stats import norm

from .types import AGENT_DTYPE, LarvalCohort, Stage, N_LOCI
from .config import SpawningSection, DiseaseSection


# ═══════════════════════════════════════════════════════════════════════
# SPAWNING SEASON UTILITIES
# ═══════════════════════════════════════════════════════════════════════

def in_spawning_season(doy: int, season_start: int = 305, season_end: int = 196) -> bool:
    """Check if a day-of-year is within the spawning season.
    
    Handles year-boundary wrapping (Nov-July season spans across calendar year).
    Default season: day 305 (Nov 1) through day 196 (Jul 15) = ~270 days.
    
    Args:
        doy: Day of year (1-365/366).
        season_start: Season start day (default 305 = ~Nov 1).
        season_end: Season end day (default 196 = ~Jul 15).
        
    Returns:
        True if doy is within spawning season.
        
    Examples:
        in_spawning_season(100, 305, 196)  # True - April (spring)
        in_spawning_season(320, 305, 196)  # True - November (early season)  
        in_spawning_season(250, 305, 196)  # False - September (out of season)
    """
    if season_start <= season_end:
        # Normal case: season doesn't cross year boundary
        return season_start <= doy <= season_end
    else:
        # Wrapped case: season crosses year boundary (Nov-Jul)
        return doy >= season_start or doy <= season_end


def seasonal_readiness_prob(doy: int, peak_doy: int, peak_width: float) -> float:
    """Calculate seasonal spawning readiness probability envelope.
    
    Uses Normal PDF centered on peak_doy with standard deviation peak_width.
    Normalized so maximum probability = 1.0 at peak_doy.
    
    For day-of-year wrapping (spawning season crosses year boundary), the
    Normal PDF handles the circular nature correctly by evaluating both
    the direct distance and the wrapped distance, taking the minimum.
    
    Args:
        doy: Day of year (1-365).
        peak_doy: Peak spawning day of year.
        peak_width: Standard deviation of spawning peak (days).
        
    Returns:
        Seasonal readiness probability [0.0, 1.0].
        
    Examples:
        seasonal_readiness_prob(105, 105, 45)  # 1.0 (peak day)
        seasonal_readiness_prob(150, 105, 45)  # ~0.61 (45 days from peak)
        seasonal_readiness_prob(305, 105, 45)  # ~0.13 (early season)
    """
    # Handle year wrapping by computing both direct and wrapped distances
    direct_dist = abs(doy - peak_doy)
    wrapped_dist = 365 - direct_dist  # Distance going the other way around the year
    
    # Use the shorter of the two distances
    min_dist = min(direct_dist, wrapped_dist)
    
    # Normal PDF normalized to max=1.0 at peak
    prob = np.exp(-0.5 * (min_dist / peak_width) ** 2)
    return prob


def latitude_adjusted_peak(base_peak_doy: int, latitude: float, lat_shift_per_deg: float = 3.0) -> int:
    """Adjust spawning peak day-of-year based on latitude.
    
    Higher latitudes have later spawning peaks. Default shift of 3 days per degree
    northward follows observed patterns in other marine invertebrates.
    
    Args:
        base_peak_doy: Base peak day (e.g., 105 for ~Apr 15).
        latitude: Node latitude (degrees North).
        lat_shift_per_deg: Days to shift peak per degree latitude (default 3.0).
        
    Returns:
        Adjusted peak day-of-year (wrapped to 1-365 range).
        
    Examples:
        latitude_adjusted_peak(105, 50.0)  # Later peak at higher latitude
        latitude_adjusted_peak(105, 35.0)  # Earlier peak at lower latitude
    """
    # Assume reference latitude is ~40°N (central California)
    reference_latitude = 40.0
    lat_offset = (latitude - reference_latitude) * lat_shift_per_deg
    
    adjusted_peak = base_peak_doy + int(round(lat_offset))
    
    # Wrap to valid day-of-year range
    if adjusted_peak < 1:
        adjusted_peak += 365
    elif adjusted_peak > 365:
        adjusted_peak -= 365
        
    return adjusted_peak


# ═══════════════════════════════════════════════════════════════════════
# PHASE 1 SPAWNING STEP (Spontaneous spawning only)
# ═══════════════════════════════════════════════════════════════════════

def spawning_step(
    agents: np.ndarray,
    genotypes: np.ndarray,
    day_of_year: int,
    node_latitude: float,
    spawning_config: SpawningSection,
    disease_config: DiseaseSection,
    rng: np.random.Generator,
) -> List[LarvalCohort]:
    """Execute daily spawning step for Phase 1-4 (complete spawning system).
    
    Updates agent spawning state fields and produces larval cohorts from
    spawning events. Includes:
    - Seasonal readiness updates
    - Spontaneous spawning 
    - Cascade induction from recent spawners (Phase 2)
    - Spawning gravity timers (Phase 3)
    - Post-spawning immunosuppression (Phase 4)
    - Female single-spawn enforcement
    - Male multi-bout tracking with refractory periods
    
    Args:
        agents: Agent array to modify.
        genotypes: Genotype array (parallel to agents).
        day_of_year: Current day of year (1-365).
        node_latitude: Node latitude for peak adjustment.
        spawning_config: Spawning configuration parameters.
        disease_config: Disease configuration parameters (for immunosuppression).
        rng: Random number generator.
        
    Returns:
        List of LarvalCohort objects from today's spawning events.
        
    Side effects:
        - Updates agents['spawning_ready'] (new ready individuals)
        - Updates agents['has_spawned'] (spawning counts/flags)
        - Updates agents['spawn_refractory'] (male refractory timers)
        - Updates agents['spawn_gravity_timer'] (gravity phase timers)
        - Updates agents['immunosuppression_timer'] (post-spawn immunosuppression)
        - Updates agents['last_spawn_day'] (tracking for future cascade)
    """
    cohorts = []
    
    # Only process during spawning season
    if not in_spawning_season(day_of_year, spawning_config.season_start_doy, spawning_config.season_end_doy):
        return cohorts
    
    # Get alive adult mask
    alive_mask = agents['alive']
    adult_mask = agents['stage'] == Stage.ADULT
    mature_mask = alive_mask & adult_mask
    
    if not np.any(mature_mask):
        return cohorts
    
    # Get mature adults
    mature_adults = agents[mature_mask]
    mature_indices = np.where(mature_mask)[0]
    
    # 1. Update spawning readiness for non-ready adults
    not_ready_mask = mature_adults['spawning_ready'] == 0
    if np.any(not_ready_mask):
        # Calculate latitude-adjusted peak and readiness probability
        adjusted_peak = latitude_adjusted_peak(
            spawning_config.peak_doy, 
            node_latitude, 
            spawning_config.lat_shift_per_deg
        )
        readiness_prob = seasonal_readiness_prob(
            day_of_year, 
            adjusted_peak, 
            spawning_config.peak_width_days
        )
        
        # Roll for readiness
        not_ready_indices = mature_indices[not_ready_mask]
        readiness_rolls = rng.random(len(not_ready_indices))
        newly_ready = readiness_rolls < readiness_prob
        
        if np.any(newly_ready):
            newly_ready_indices = not_ready_indices[newly_ready]
            agents['spawning_ready'][newly_ready_indices] = 1
            # Phase 3: Set spawning gravity timer when entering readiness
            if spawning_config.gravity_enabled:
                gravity_duration = spawning_config.pre_spawn_gravity_days + spawning_config.post_spawn_gravity_days
                agents['spawn_gravity_timer'][newly_ready_indices] = gravity_duration
    
    # 2. Decrement timers
    # 2a. Male refractory timers
    refractory_mask = agents['spawn_refractory'] > 0
    agents['spawn_refractory'][refractory_mask] -= 1
    
    # 2b. Spawning gravity timers (Phase 3)
    gravity_mask = agents['spawn_gravity_timer'] > 0
    agents['spawn_gravity_timer'][gravity_mask] -= 1
    
    # 2c. Post-spawning immunosuppression timers (Phase 4)
    immuno_mask = agents['immunosuppression_timer'] > 0
    agents['immunosuppression_timer'][immuno_mask] -= 1
    
    # 3. Spontaneous spawning attempts
    spawners_today = []
    
    # 3a. Ready females who haven't spawned yet
    female_mask = (mature_adults['sex'] == 0) & (mature_adults['spawning_ready'] == 1) & (mature_adults['has_spawned'] == 0)
    if np.any(female_mask):
        female_indices = mature_indices[female_mask]
        spawn_rolls = rng.random(len(female_indices))
        spawning_females = spawn_rolls < spawning_config.p_spontaneous_female
        
        if np.any(spawning_females):
            spawning_female_indices = female_indices[spawning_females]
            _execute_spawning_events(agents, spawning_female_indices, day_of_year, spawning_config, disease_config)
            spawners_today.extend(spawning_female_indices)
    
    # 3b. Ready males who can spawn (not in refractory, under bout limit)
    male_mask = (
        (mature_adults['sex'] == 1) & 
        (mature_adults['spawning_ready'] == 1) & 
        (mature_adults['has_spawned'] < spawning_config.male_max_bouts) &
        (mature_adults['spawn_refractory'] == 0)
    )
    if np.any(male_mask):
        male_indices = mature_indices[male_mask]
        spawn_rolls = rng.random(len(male_indices))
        spawning_males = spawn_rolls < spawning_config.p_spontaneous_male
        
        if np.any(spawning_males):
            spawning_male_indices = male_indices[spawning_males]
            _execute_spawning_events(agents, spawning_male_indices, day_of_year, spawning_config, disease_config, spawning_config.male_refractory_days)
            spawners_today.extend(spawning_male_indices)
    
    # 4. Phase 2: Cascade induction from recent spawners
    cascade_spawners = _cascade_induction_step(
        agents,
        mature_indices,
        day_of_year,
        spawning_config,
        disease_config,
        rng
    )
    spawners_today.extend(cascade_spawners)

    # 5. Generate larval cohorts from today's spawners
    if spawners_today:
        cohort = _generate_larval_cohort(
            agents, 
            genotypes, 
            spawners_today, 
            spawning_config, 
            rng
        )
        if cohort.n_competent > 0:
            cohorts.append(cohort)
    
    return cohorts


def _cascade_induction_step(
    agents: np.ndarray,
    mature_indices: np.ndarray,
    day_of_year: int,
    spawning_config: SpawningSection,
    disease_config: DiseaseSection,
    rng: np.random.Generator,
) -> List[int]:
    """Execute cascade induction from recent spawners (Phase 2).
    
    Finds agents who spawned within the cascade window and checks if they can
    trigger spawning in nearby ready individuals. Sex-asymmetric induction:
    - Females strongly induce males (κ_fm = 0.80)
    - Males moderately induce females (κ_mf = 0.30)
    
    Args:
        agents: Agent array to modify.
        mature_indices: Indices of mature adults to consider.
        day_of_year: Current day of year.
        config: Spawning configuration parameters.
        rng: Random number generator.
        
    Returns:
        List of agent indices that spawned due to cascade induction.
        
    Side effects:
        - Updates agents['has_spawned'] for triggered spawners
        - Updates agents['spawn_refractory'] for triggered males
        - Updates agents['last_spawn_day'] for triggered spawners
    """
    cascade_spawners = []
    
    if len(mature_indices) == 0:
        return cascade_spawners
    
    mature_adults = agents[mature_indices]
    
    # Find recent spawners within cascade window
    recent_spawners_mask = _get_recent_spawners_mask(
        mature_adults, day_of_year, spawning_config.cascade_window
    )
    
    if not np.any(recent_spawners_mask):
        return cascade_spawners  # No recent spawners to trigger cascade
    
    recent_spawner_indices = mature_indices[recent_spawners_mask]
    recent_spawners = agents[recent_spawner_indices]
    
    # Separate recent spawners by sex
    recent_female_mask = recent_spawners['sex'] == 0
    recent_male_mask = recent_spawners['sex'] == 1
    
    recent_female_indices = recent_spawner_indices[recent_female_mask]
    recent_male_indices = recent_spawner_indices[recent_male_mask]
    
    # Find potential cascade targets (ready individuals who haven't spawned today)
    target_mask = (
        (mature_adults['spawning_ready'] == 1) &
        (mature_adults['last_spawn_day'] != day_of_year)  # Haven't spawned today
    )
    
    if not np.any(target_mask):
        return cascade_spawners
    
    target_indices = mature_indices[target_mask]
    targets = agents[target_indices]
    
    # Process male targets (induced by recent female spawners)
    if len(recent_female_indices) > 0:
        male_target_mask = (
            (targets['sex'] == 1) &
            (targets['has_spawned'] < spawning_config.male_max_bouts) &
            (targets['spawn_refractory'] == 0)
        )
        
        if np.any(male_target_mask):
            male_target_indices = target_indices[male_target_mask]
            triggered_males = _check_cascade_induction(
                agents, male_target_indices, recent_female_indices,
                spawning_config.cascade_radius, spawning_config.induction_female_to_male, rng
            )
            
            if triggered_males:
                _execute_spawning_events(
                    agents, triggered_males, day_of_year, spawning_config, disease_config, spawning_config.male_refractory_days
                )
                cascade_spawners.extend(triggered_males)
    
    # Process female targets (induced by recent male spawners)
    if len(recent_male_indices) > 0:
        female_target_mask = (
            (targets['sex'] == 0) &
            (targets['has_spawned'] == 0)  # Females spawn only once
        )
        
        if np.any(female_target_mask):
            female_target_indices = target_indices[female_target_mask]
            triggered_females = _check_cascade_induction(
                agents, female_target_indices, recent_male_indices,
                spawning_config.cascade_radius, spawning_config.induction_male_to_female, rng
            )
            
            if triggered_females:
                _execute_spawning_events(
                    agents, triggered_females, day_of_year, spawning_config, disease_config
                )
                cascade_spawners.extend(triggered_females)
    
    return cascade_spawners


def _get_recent_spawners_mask(
    agents: np.ndarray, 
    current_doy: int, 
    cascade_window: int
) -> np.ndarray:
    """Get mask of agents that spawned within the cascade window.
    
    Handles day-of-year wrapping across year boundaries.
    
    Args:
        agents: Agent array to check.
        current_doy: Current day of year (1-365).
        cascade_window: Number of days to look back.
        
    Returns:
        Boolean mask indicating recent spawners.
    """
    last_spawn_days = agents['last_spawn_day']
    
    # Handle agents that never spawned (last_spawn_day = 0)
    never_spawned_mask = last_spawn_days == 0
    
    # Calculate days since last spawn, handling year wrapping
    days_since_spawn = np.zeros_like(last_spawn_days)
    
    for i, last_spawn in enumerate(last_spawn_days):
        if last_spawn == 0:
            days_since_spawn[i] = 999  # Large number for never-spawned
        elif last_spawn <= current_doy:
            # Same year
            days_since_spawn[i] = current_doy - last_spawn
        else:
            # Last spawn was later in previous year, current is early this year
            days_since_spawn[i] = (365 - last_spawn) + current_doy
    
    # Recent spawners are those within cascade window
    recent_mask = (days_since_spawn <= cascade_window) & ~never_spawned_mask
    
    return recent_mask


def _check_cascade_induction(
    agents: np.ndarray,
    target_indices: np.ndarray,
    inducer_indices: np.ndarray,
    cascade_radius: float,
    induction_probability: float,
    rng: np.random.Generator,
) -> List[int]:
    """Check which targets are induced by nearby recent spawners.
    
    Args:
        agents: Full agent array.
        target_indices: Indices of potential targets for induction.
        inducer_indices: Indices of recent spawners that could induce.
        cascade_radius: Maximum distance for induction (meters).
        induction_probability: Probability of induction if within range.
        rng: Random number generator.
        
    Returns:
        List of target indices that were induced to spawn.
    """
    induced_targets = []
    
    if len(target_indices) == 0 or len(inducer_indices) == 0:
        return induced_targets
    
    # Get positions
    target_positions = np.column_stack([
        agents['x'][target_indices], 
        agents['y'][target_indices]
    ])
    inducer_positions = np.column_stack([
        agents['x'][inducer_indices], 
        agents['y'][inducer_indices]
    ])
    
    # Check each target
    for i, target_idx in enumerate(target_indices):
        target_pos = target_positions[i]
        
        # Calculate distances to all recent spawners of opposite sex
        distances = np.sqrt(np.sum((inducer_positions - target_pos)**2, axis=1))
        
        # Check if any inducer is within cascade radius
        if np.any(distances <= cascade_radius):
            # Roll for induction
            if rng.random() < induction_probability:
                induced_targets.append(target_idx)
    
    return induced_targets


def _execute_spawning_events(
    agents: np.ndarray, 
    spawner_indices: np.ndarray, 
    day_of_year: int,
    spawning_config: SpawningSection,
    disease_config: DiseaseSection,
    refractory_days: Optional[int] = None
) -> None:
    """Execute spawning events for specified agents.
    
    Updates spawning status fields. For females, sets has_spawned=1 (single spawn).
    For males, increments has_spawned bout count and sets refractory period.
    Phase 4: Sets post-spawning immunosuppression timers.
    
    Args:
        agents: Agent array to modify.
        spawner_indices: Indices of spawning agents.
        day_of_year: Current day of year.
        spawning_config: Spawning configuration parameters.
        disease_config: Disease configuration parameters (for immunosuppression).
        refractory_days: Refractory period for males (None for females).
    """
    # Record spawning day for all spawners
    agents['last_spawn_day'][spawner_indices] = day_of_year
    
    # Update spawning counts
    if refractory_days is None:
        # Females: single spawn
        agents['has_spawned'][spawner_indices] = 1
    else:
        # Males: increment bout count and set refractory
        agents['has_spawned'][spawner_indices] += 1
        agents['spawn_refractory'][spawner_indices] = refractory_days
    
    # Phase 4: Set post-spawning immunosuppression timer
    if disease_config.immunosuppression_enabled:
        # For males spawning multiple bouts, timer resets each time
        agents['immunosuppression_timer'][spawner_indices] = disease_config.immunosuppression_duration


def _generate_larval_cohort(
    agents: np.ndarray,
    genotypes: np.ndarray,
    spawner_indices: List[int],
    config: SpawningSection,
    rng: np.random.Generator,
) -> LarvalCohort:
    """Generate larval cohort from spawning agents using SRS lottery.
    
    Implements sweepstakes reproductive success (SRS) where a small fraction
    of spawners contributes most offspring. Uses Pareto distribution for
    reproductive success as in reproduction-srs-spec.md.
    
    Args:
        agents: Full agent array.
        genotypes: Full genotype array.
        spawner_indices: Indices of agents that spawned today.
        config: Spawning configuration (contains SRS parameters from population config).
        rng: Random number generator.
        
    Returns:
        LarvalCohort with competent larvae and parent tracking.
        
    Note:
        This function will need to access population config parameters
        for fecundity and SRS calculations. For now using placeholder values.
    """
    if not spawner_indices:
        return LarvalCohort(
            source_node=agents['node_id'][spawner_indices[0]] if spawner_indices else 0,
            n_competent=0,
            genotypes=np.array([], dtype=np.int8).reshape(0, N_LOCI, 2),
            parent_pairs=np.array([], dtype=np.int32).reshape(0, 2),
            pld_days=120.0  # Default PLD
        )
    
    spawners = agents[spawner_indices]
    source_node = spawners['node_id'][0]  # Assume all same node
    
    # Separate by sex
    female_mask = spawners['sex'] == 0
    male_mask = spawners['sex'] == 1
    
    females = spawners[female_mask]
    males = spawners[male_mask]
    
    female_indices = np.array(spawner_indices)[female_mask]
    male_indices = np.array(spawner_indices)[male_mask]
    
    if len(females) == 0 or len(males) == 0:
        # Need both sexes for reproduction
        return LarvalCohort(
            source_node=source_node,
            n_competent=0,
            genotypes=np.array([], dtype=np.int8).reshape(0, N_LOCI, 2),
            parent_pairs=np.array([], dtype=np.int32).reshape(0, 2),
            pld_days=120.0
        )
    
    # Calculate total fecundity from females
    # Using simplified fecundity calculation for Phase 1
    # F_i = F0 * (L_i / L_ref)^fecundity_exp for L_i >= L_min_repro
    F0 = 1.0e7  # Will be from population config in full integration
    L_ref = 500.0
    fecundity_exp = 2.5
    L_min_repro = 400.0
    
    reproductive_females = females[females['size'] >= L_min_repro]
    if len(reproductive_females) == 0:
        return LarvalCohort(
            source_node=source_node,
            n_competent=0,
            genotypes=np.array([], dtype=np.int8).reshape(0, N_LOCI, 2),
            parent_pairs=np.array([], dtype=np.int32).reshape(0, 2),
            pld_days=120.0
        )
    
    female_fecundities = F0 * (reproductive_females['size'] / L_ref) ** fecundity_exp
    total_eggs = np.sum(female_fecundities)
    
    # Apply SRS lottery - simplified version for Phase 1
    # Pareto-distributed reproductive success with α=1.35
    alpha_srs = 1.35
    n_breeding_pairs = min(len(reproductive_females), len(males))
    
    # Sample breeding pairs with SRS weights
    if n_breeding_pairs == 1:
        # Single pair case
        selected_females = reproductive_females[:1]
        selected_males = males[:1]
        srs_weights = np.array([1.0])
        female_order = np.array([0])
        male_order = np.array([0])
    else:
        # Multiple pairs: use Pareto sampling for SRS
        pareto_samples = rng.pareto(alpha_srs, n_breeding_pairs)
        srs_weights = pareto_samples / np.sum(pareto_samples)
        
        # Select top pairs by size (simplified selection)
        female_order = np.argsort(reproductive_females['size'])[::-1]
        male_order = np.argsort(males['size'])[::-1]
        
        selected_females = reproductive_females[female_order[:n_breeding_pairs]]
        selected_males = males[male_order[:n_breeding_pairs]]
    
    # Calculate larvae per breeding pair
    eggs_per_pair = total_eggs * srs_weights
    
    # Convert to competent larvae (simplified - no detailed larval dynamics)
    competency_rate = 0.1  # 10% of eggs become competent larvae (placeholder)
    larvae_per_pair = (eggs_per_pair * competency_rate).astype(int)
    
    # Cap total larvae for computational efficiency (especially in tests)
    # This prevents performance issues when testing with high spawning probabilities
    max_larvae_per_cohort = 10000  # Reasonable limit for simulation performance
    if np.sum(larvae_per_pair) > max_larvae_per_cohort:
        # Scale down proportionally to stay under cap
        scale_factor = max_larvae_per_cohort / np.sum(larvae_per_pair)
        larvae_per_pair = (larvae_per_pair * scale_factor).astype(int)
        # Ensure at least one larva per pair if originally > 0
        larvae_per_pair = np.maximum(larvae_per_pair, (eggs_per_pair > 0).astype(int))
    
    total_competent = np.sum(larvae_per_pair)
    
    if total_competent == 0:
        return LarvalCohort(
            source_node=source_node,
            n_competent=0,
            genotypes=np.array([], dtype=np.int8).reshape(0, N_LOCI, 2),
            parent_pairs=np.array([], dtype=np.int32).reshape(0, 2),
            pld_days=120.0
        )
    
    # Generate offspring genotypes through sexual reproduction
    offspring_genotypes = np.zeros((total_competent, N_LOCI, 2), dtype=np.int8)
    parent_pairs = np.zeros((total_competent, 2), dtype=np.int32)
    
    # Map selected individuals back to original indices
    selected_female_indices = []
    selected_male_indices = []
    
    for i in range(n_breeding_pairs):
        # Find the original indices for selected breeding pairs
        # Use simple indexing since we selected by size order
        female_original_idx = female_indices[female_order[i]] if i < len(female_indices) else female_indices[0]
        male_original_idx = male_indices[male_order[i]] if i < len(male_indices) else male_indices[0]
        
        selected_female_indices.append(female_original_idx)
        selected_male_indices.append(male_original_idx)
    
    larva_idx = 0
    for i in range(n_breeding_pairs):
        n_offspring = larvae_per_pair[i]
        if n_offspring == 0:
            continue
            
        female_idx = selected_female_indices[i]
        male_idx = selected_male_indices[i]
        
        # Generate offspring through Mendelian inheritance
        for j in range(n_offspring):
            if larva_idx >= total_competent:
                break
                
            # Sexual reproduction: random allele from each parent at each locus
            for locus in range(N_LOCI):
                maternal_allele = genotypes[female_idx, locus, rng.integers(0, 2)]
                paternal_allele = genotypes[male_idx, locus, rng.integers(0, 2)]
                
                offspring_genotypes[larva_idx, locus, 0] = maternal_allele
                offspring_genotypes[larva_idx, locus, 1] = paternal_allele
            
            parent_pairs[larva_idx] = [female_idx, male_idx]
            larva_idx += 1
    
    return LarvalCohort(
        source_node=source_node,
        n_competent=int(total_competent),
        genotypes=offspring_genotypes,
        parent_pairs=parent_pairs,
        pld_days=120.0  # Placeholder PLD
    )


# ═══════════════════════════════════════════════════════════════════════
# SEASON RESET
# ═══════════════════════════════════════════════════════════════════════

def reset_spawning_season(agents: np.ndarray) -> None:
    """Reset spawning fields at the end of each spawning season.
    
    Called when transitioning from spawning season to non-spawning season.
    Resets readiness and spawning counts for the next season.
    
    Args:
        agents: Agent array to reset.
        
    Side effects:
        - Sets spawning_ready = 0 for all agents
        - Sets has_spawned = 0 for all agents  
        - Keeps other timers (refractory, gravity, immunosuppression) running
    """
    agents['spawning_ready'][:] = 0
    agents['has_spawned'][:] = 0
    # Note: Refractory, gravity, and immunosuppression timers continue
    # counting down as they may extend beyond the spawning season