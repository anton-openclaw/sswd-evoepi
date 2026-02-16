#!/usr/bin/env python3
"""
Simplified spawning diagnostic script.

Author: Anton 🔬
Date: 2026-02-15
"""

import numpy as np
import pandas as pd
from pathlib import Path

from sswd_evoepi.types import AGENT_DTYPE, Stage, N_LOCI
from sswd_evoepi.config import SimulationConfig, SpawningSection, DiseaseSection
from sswd_evoepi.spawning import spawning_step, seasonal_readiness_prob, latitude_adjusted_peak


def create_test_agents(n_agents: int, habitat_area_m2: float, node_id: int = 1) -> np.ndarray:
    """Create test population of spawning-ready adults."""
    
    # Distribute agents uniformly over square habitat
    side_length = np.sqrt(habitat_area_m2)
    
    agents = np.zeros(n_agents, dtype=AGENT_DTYPE)
    
    # Initialize basic agent properties
    for i in range(n_agents):
        agents[i]['node_id'] = node_id
        agents[i]['x'] = np.random.uniform(0, side_length)
        agents[i]['y'] = np.random.uniform(0, side_length)
        agents[i]['alive'] = True
        agents[i]['stage'] = Stage.ADULT
        agents[i]['age'] = 2.0  # Mature adults (2 years)
        agents[i]['size'] = 450.0  # Above reproduction threshold
        agents[i]['sex'] = np.random.choice([0, 1])  # Random sex
        agents[i]['disease_state'] = 0  # Susceptible
        
        # Initialize spawning fields
        agents[i]['spawning_ready'] = 0
        agents[i]['has_spawned'] = 0
        agents[i]['spawn_refractory'] = 0
        agents[i]['spawn_gravity_timer'] = 0
        agents[i]['immunosuppression_timer'] = 0
        agents[i]['last_spawn_day'] = 0
    
    return agents


def create_test_genotypes(n_agents: int) -> np.ndarray:
    """Create simple random genotypes."""
    genotypes = np.zeros((n_agents, N_LOCI, 2), dtype=np.int8)
    
    # Random alleles (0 or 1) with 50% frequency
    for i in range(n_agents):
        for locus in range(N_LOCI):
            genotypes[i, locus, 0] = np.random.choice([0, 1])
            genotypes[i, locus, 1] = np.random.choice([0, 1])
    
    return genotypes


def diagnose_spawning_parameters():
    """Run comprehensive spawning diagnosis."""
    
    print("=" * 60)
    print("SPAWNING SYSTEM DIAGNOSTIC ANALYSIS")
    print("=" * 60)
    
    # Test parameters
    n_agents = 500
    habitat_area = 10000  # m² (100m x 100m)
    node_latitude = 48.5  # Friday Harbor latitude
    test_day = 105  # Peak spawning day (Apr 15)
    n_test_days = 30
    
    print(f"Test setup:")
    print(f"  - {n_agents} agents in {habitat_area} m² habitat")
    print(f"  - Habitat side length: {np.sqrt(habitat_area):.1f} m")
    print(f"  - Expected nearest neighbor distance: {np.sqrt(habitat_area/n_agents):.1f} m")
    print(f"  - Node latitude: {node_latitude}°N")
    print(f"  - Testing days {test_day} to {test_day + n_test_days - 1}")
    print()
    
    # Load config
    config = SimulationConfig()
    spawning_config = config.spawning
    disease_config = config.disease
    
    print("Current spawning parameters:")
    print(f"  - p_spontaneous_female: {spawning_config.p_spontaneous_female}")
    print(f"  - p_spontaneous_male: {spawning_config.p_spontaneous_male}")
    print(f"  - cascade radius: {spawning_config.cascade_radius} m")
    print(f"  - induction_female_to_male: {spawning_config.induction_female_to_male}")
    print(f"  - induction_male_to_female: {spawning_config.induction_male_to_female}")
    print(f"  - cascade window: {spawning_config.cascade_window} days")
    print(f"  - peak_doy: {spawning_config.peak_doy}")
    print(f"  - peak_width_days: {spawning_config.peak_width_days}")
    print()
    
    # Set random seed for reproducibility
    np.random.seed(42)
    
    # Create test population
    agents = create_test_agents(n_agents, habitat_area)
    genotypes = create_test_genotypes(n_agents)
    
    rng = np.random.default_rng(42)
    
    # Run daily simulation
    daily_results = []
    total_spawning_events = 0
    total_larvae = 0
    
    for day in range(test_day, test_day + n_test_days):
        day_result = {"day": day}
        
        # Calculate seasonal readiness probability
        adjusted_peak = latitude_adjusted_peak(
            spawning_config.peak_doy, node_latitude, spawning_config.lat_shift_per_deg
        )
        readiness_prob = seasonal_readiness_prob(day, adjusted_peak, spawning_config.peak_width_days)
        day_result["readiness_prob"] = readiness_prob
        
        # Count agents by status before spawning step
        alive_adults = np.sum((agents['alive']) & (agents['stage'] == Stage.ADULT))
        spawning_ready = np.sum(agents['spawning_ready'])
        ready_females = np.sum((agents['spawning_ready']) & (agents['sex'] == 0) & (agents['has_spawned'] == 0))
        ready_males = np.sum((agents['spawning_ready']) & (agents['sex'] == 1) & 
                           (agents['has_spawned'] < spawning_config.male_max_bouts) &
                           (agents['spawn_refractory'] == 0))
        
        day_result.update({
            "alive_adults": alive_adults,
            "spawning_ready": spawning_ready,
            "ready_females": ready_females, 
            "ready_males": ready_males
        })
        
        # Calculate effective spawning probabilities
        eff_p_female = spawning_config.p_spontaneous_female * readiness_prob
        eff_p_male = spawning_config.p_spontaneous_male * readiness_prob
        
        day_result.update({
            "eff_p_female": eff_p_female,
            "eff_p_male": eff_p_male,
            "expected_female_spawns": ready_females * eff_p_female,
            "expected_male_spawns": ready_males * eff_p_male
        })
        
        # Run spawning step
        spawners_before = np.sum(agents['has_spawned'])
        
        cohorts = spawning_step(
            agents, genotypes, day, node_latitude,
            spawning_config, disease_config, rng
        )
        
        spawners_after = np.sum(agents['has_spawned'])
        
        # Count events today
        new_spawning_events = spawners_after - spawners_before
        
        # Count larvae from cohorts
        larvae_today = sum(cohort.n_competent for cohort in cohorts)
        
        day_result.update({
            "spawning_events_today": new_spawning_events,
            "larvae_today": larvae_today,
            "n_cohorts": len(cohorts)
        })
        
        # Accumulate totals
        total_spawning_events += new_spawning_events
        total_larvae += larvae_today
        
        daily_results.append(day_result)
        
        # Print detailed info for first few days
        if day <= test_day + 4:
            print(f"Day {day} (DOY):")
            print(f"  Readiness prob: {readiness_prob:.3f}")
            print(f"  Ready agents: {spawning_ready} ({ready_females} F, {ready_males} M)")
            print(f"  Effective spawn prob: F={eff_p_female:.4f}, M={eff_p_male:.4f}")
            print(f"  Expected spawns: F={ready_females * eff_p_female:.2f}, M={ready_males * eff_p_male:.2f}")
            print(f"  Actual spawning events: {new_spawning_events}")
            print(f"  Larvae produced: {larvae_today}")
            print()
    
    # Analysis summary
    print("=" * 60)
    print("SUMMARY ANALYSIS")
    print("=" * 60)
    
    df = pd.DataFrame(daily_results)
    
    print(f"Mean readiness probability: {df['readiness_prob'].mean():.3f}")
    print(f"Peak readiness probability: {df['readiness_prob'].max():.3f}")
    print(f"Mean effective female spawn prob: {df['eff_p_female'].mean():.4f}")
    print(f"Mean effective male spawn prob: {df['eff_p_male'].mean():.4f}")
    print()
    
    print(f"Total spawning events in {n_test_days} days: {total_spawning_events}")
    print(f"Mean spawning events per day: {total_spawning_events / n_test_days:.2f}")
    print(f"Total larvae produced: {total_larvae}")
    print(f"Mean larvae per day: {total_larvae / n_test_days:.0f}")
    print()
    
    # Calculate male bout statistics
    male_agents = agents[agents['sex'] == 1]
    male_bout_counts = male_agents['has_spawned']
    mean_male_bouts = np.mean(male_bout_counts)
    males_that_spawned = np.sum(male_bout_counts > 0)
    
    print(f"Male spawning statistics:")
    print(f"  Total males: {len(male_agents)}")
    print(f"  Males that spawned: {males_that_spawned}")
    print(f"  Mean bouts per male: {mean_male_bouts:.3f}")
    print(f"  Target mean bouts: ~2.2")
    print(f"  Current/Target ratio: {mean_male_bouts / 2.2:.3f}")
    print()
    
    # ROOT CAUSE ANALYSIS
    print("=" * 60)
    print("ROOT CAUSE ANALYSIS")
    print("=" * 60)
    
    print("1. SPAWNING FREQUENCY:")
    expected_daily_spawns = (n_agents/2) * df['eff_p_female'].mean()  # ~250 females
    actual_daily_spawns = total_spawning_events / n_test_days
    print(f"   Expected daily female spawns: {expected_daily_spawns:.2f}")
    print(f"   Actual daily spawning events: {actual_daily_spawns:.2f}")
    if expected_daily_spawns > 0:
        print(f"   Efficiency ratio: {actual_daily_spawns / expected_daily_spawns:.3f}")
    print()
    
    if actual_daily_spawns < expected_daily_spawns * 0.5:
        print("   ❌ ISSUE: Spawning probability too low")
        if expected_daily_spawns > 0:
            multiplier = expected_daily_spawns / actual_daily_spawns
            print(f"   RECOMMENDATION: Increase p_spontaneous by {multiplier:.1f}×")
    else:
        print("   ✅ Spawning frequency roughly as expected")
    print()
    
    print("2. SEASONAL READINESS ENVELOPE:")
    print(f"   Peak readiness prob: {df['readiness_prob'].max():.3f}")
    print(f"   Mean readiness prob: {df['readiness_prob'].mean():.3f}")
    if df['readiness_prob'].mean() < 0.5:
        print("   ❌ ISSUE: Seasonal envelope may be too narrow")
        print("   RECOMMENDATION: Consider wider peak_width_days or check peak timing")
    else:
        print("   ✅ Seasonal envelope seems reasonable")
    print()
    
    print("3. MALE BOUT TARGET ANALYSIS:")
    target_bouts = 2.2
    current_bouts = mean_male_bouts
    print(f"   Target: {target_bouts} bouts/male/season")
    print(f"   Current: {current_bouts:.3f} bouts/male")
    print(f"   Deficit: {target_bouts - current_bouts:.3f}")
    
    if current_bouts < target_bouts * 0.5:
        # Very low - major issue
        print("   ❌ CRITICAL ISSUE: Male bout count extremely low")
        multiplier = target_bouts / max(current_bouts, 0.01)
        suggested_p_male = spawning_config.p_spontaneous_male * multiplier
        print(f"   RECOMMENDATION: Increase p_spontaneous_male to ~{suggested_p_male:.3f}")
    elif current_bouts < target_bouts * 0.8:
        # Moderately low
        print("   ⚠️ ISSUE: Male bout count below target")
        multiplier = target_bouts / current_bouts
        suggested_p_male = spawning_config.p_spontaneous_male * multiplier
        print(f"   RECOMMENDATION: Increase p_spontaneous_male to ~{suggested_p_male:.3f}")
    else:
        print("   ✅ Male bout count acceptable")
    print()
    
    print("4. LARVAE PRODUCTION ANALYSIS:")
    if total_larvae > 0:
        larvae_per_spawner = total_larvae / total_spawning_events if total_spawning_events > 0 else 0
        print(f"   Larvae per spawning event: {larvae_per_spawner:.0f}")
        print(f"   Total larvae in {n_test_days} days: {total_larvae}")
        
        # Extrapolate to full season
        season_larvae = total_larvae * (270 / n_test_days)
        print(f"   Projected season larvae: {season_larvae:.0f}")
        
        if season_larvae > 100000:  # Very high threshold
            print("   ❌ ISSUE: Larvae production extremely high")
            print("   RECOMMENDATION: Check larval cohort generation parameters")
        elif season_larvae > 10000:
            print("   ⚠️ WARNING: Larvae production may be high")
            print("   RECOMMENDATION: Verify larval cohort scaling")
        else:
            print("   ✅ Larvae production seems reasonable")
    else:
        print("   ❌ ISSUE: No larvae produced")
        print("   RECOMMENDATION: Check larval cohort generation logic")
    
    print()
    
    # SPECIFIC RECOMMENDATIONS
    print("=" * 60)  
    print("SPECIFIC PARAMETER RECOMMENDATIONS")
    print("=" * 60)
    
    print("Based on diagnosis:")
    print()
    
    # Calculate specific recommendations
    if current_bouts < target_bouts:
        # Male spawning probability
        if current_bouts > 0:
            male_multiplier = target_bouts / current_bouts
        else:
            male_multiplier = 10  # If zero, big increase needed
        new_p_male = min(spawning_config.p_spontaneous_male * male_multiplier, 0.2)  # Cap at 20%
        print(f"1. p_spontaneous_male: {spawning_config.p_spontaneous_male:.3f} → {new_p_male:.3f}")
        
        # Female spawning probability (proportional increase)
        new_p_female = min(spawning_config.p_spontaneous_female * male_multiplier * 0.8, 0.15)  # Cap at 15%
        print(f"2. p_spontaneous_female: {spawning_config.p_spontaneous_female:.3f} → {new_p_female:.3f}")
    
    # Seasonal envelope
    if df['readiness_prob'].mean() < 0.6:
        new_peak_width = spawning_config.peak_width_days * 1.3
        print(f"3. peak_width_days: {spawning_config.peak_width_days} → {new_peak_width:.0f}")
    
    # Larvae production
    if total_larvae / total_spawning_events > 5000:
        print(f"4. WARNING: Consider reducing larval cohort size scaling")
    
    print()
    print("Test completed successfully!")
    
    # Save results
    results_dir = Path("results/spawning_diagnosis")
    results_dir.mkdir(parents=True, exist_ok=True)
    
    df.to_csv(results_dir / "daily_spawning_analysis.csv", index=False)
    print(f"Results saved to {results_dir}/")


if __name__ == "__main__":
    diagnose_spawning_parameters()