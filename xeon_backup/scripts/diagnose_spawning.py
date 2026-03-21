#!/usr/bin/env python3
"""
Diagnostic script to understand spawning parameter issues.

Author: Anton 🔬
Date: 2026-02-15
"""

import numpy as np
import pandas as pd
from pathlib import Path

from sswd_evoepi.types import AGENT_DTYPE, Stage, N_LOCI
from sswd_evoepi.config import SimulationConfig, SpawningSection, DiseaseSection
from sswd_evoepi.spawning import spawning_step, seasonal_readiness_prob, latitude_adjusted_peak
from sswd_evoepi.agents import create_agent
from sswd_evoepi.genetics import create_genotypes


def create_test_population(n_agents: int, habitat_area_m2: float, node_id: int = 1) -> tuple:
    """Create test population of spawning-ready adults."""
    
    # Distribute agents uniformly over square habitat
    side_length = np.sqrt(habitat_area_m2)
    
    agents = np.zeros(n_agents, dtype=AGENT_DTYPE)
    genotypes = create_genotypes(n_agents, N_LOCI, np.random.default_rng(42),
                                allele_freq=0.5)
    
    for i in range(n_agents):
        agents[i] = create_agent(
            node_id=node_id,
            position=(
                np.random.uniform(0, side_length),
                np.random.uniform(0, side_length)
            ),
            stage=Stage.ADULT,
            age_days=365 * 2,  # Mature adults
            size=450.0,  # Above reproduction threshold
            sex=np.random.choice([0, 1]),  # Random sex
            genotype_idx=i
        )
        # All start alive and adult
        agents[i]['alive'] = True
        agents[i]['stage'] = Stage.ADULT
    
    return agents, genotypes


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
    
    # Create test population
    agents, genotypes = create_test_population(n_agents, habitat_area)
    
    rng = np.random.default_rng(42)
    
    # Run daily simulation
    daily_results = []
    total_spawning_events = 0
    total_male_bouts = 0
    total_cascade_events = 0
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
        spawners_before = {
            "total_spawning_events": np.sum(agents['has_spawned']),
            "spawning_ready": np.sum(agents['spawning_ready'])
        }
        
        cohorts = spawning_step(
            agents, genotypes, day, node_latitude,
            spawning_config, disease_config, rng
        )
        
        spawners_after = {
            "total_spawning_events": np.sum(agents['has_spawned']),
            "spawning_ready": np.sum(agents['spawning_ready'])
        }
        
        # Count events today
        new_spawning_events = spawners_after["total_spawning_events"] - spawners_before["total_spawning_events"]
        new_ready = spawners_after["spawning_ready"] - spawners_before["spawning_ready"]
        
        # Count larvae from cohorts
        larvae_today = sum(cohort.n_competent for cohort in cohorts)
        
        day_result.update({
            "new_ready": new_ready,
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
            print(f"  New ready: {new_ready}")
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
    
    # Calculate distance analysis for cascade effectiveness
    positions = np.column_stack([agents['x'], agents['y']])
    
    # Sample a few agents to check distances
    sample_size = min(50, len(agents))
    sample_indices = np.random.choice(len(agents), sample_size, replace=False)
    
    distances = []
    for i in sample_indices:
        pos = positions[i]
        other_distances = np.sqrt(np.sum((positions - pos)**2, axis=1))
        # Find nearest neighbors (excluding self)
        other_distances[i] = np.inf
        nearest_k = np.partition(other_distances, min(10, len(agents)-1))[:10]
        distances.extend(nearest_k[:5])  # Take 5 nearest for each sample
    
    distances = np.array(distances)
    print(f"Spatial distribution analysis:")
    print(f"  Mean distance to 5 nearest neighbors: {np.mean(distances):.1f} m")
    print(f"  Cascade radius: {spawning_config.cascade_radius} m")
    print(f"  Fraction within cascade range: {np.mean(distances <= spawning_config.cascade_radius):.3f}")
    print()
    
    # ROOT CAUSE ANALYSIS
    print("=" * 60)
    print("ROOT CAUSE ANALYSIS")
    print("=" * 60)
    
    print("1. SPAWNING FREQUENCY ISSUES:")
    expected_daily_spawns = n_agents/2 * df['eff_p_female'].mean()  # ~250 females
    actual_daily_spawns = total_spawning_events / n_test_days
    print(f"   Expected daily female spawns: {expected_daily_spawns:.2f}")
    print(f"   Actual daily spawning events: {actual_daily_spawns:.2f}")
    print(f"   Ratio: {actual_daily_spawns / expected_daily_spawns:.3f}")
    print()
    
    if actual_daily_spawns < expected_daily_spawns * 0.5:
        print("   ❌ ISSUE: Spawning probability too low")
        print(f"   RECOMMENDATION: Increase p_spontaneous by {expected_daily_spawns / actual_daily_spawns:.1f}×")
    else:
        print("   ✅ Spawning frequency roughly as expected")
    print()
    
    print("2. SEASONAL READINESS ENVELOPE:")
    print(f"   Peak readiness prob: {df['readiness_prob'].max():.3f}")
    print(f"   Mean readiness prob: {df['readiness_prob'].mean():.3f}")
    if df['readiness_prob'].mean() < 0.5:
        print("   ❌ ISSUE: Seasonal envelope may be too narrow")
        print("   RECOMMENDATION: Consider wider peak_width_days or check peak timing")
    print()
    
    print("3. MALE BOUT TARGET:")
    target_bouts = 2.2
    current_bouts = mean_male_bouts
    print(f"   Target: {target_bouts} bouts/male/season")
    print(f"   Current: {current_bouts:.3f} bouts/male")
    print(f"   Deficit: {target_bouts - current_bouts:.3f}")
    
    # Assuming 270-day season, what daily prob needed?
    season_days = 270
    current_daily_male_prob = df['eff_p_male'].mean()
    needed_daily_prob = target_bouts / (season_days * (1 - spawning_config.male_refractory_days/season_days))
    
    print(f"   Current male daily prob: {current_daily_male_prob:.4f}")
    print(f"   Needed daily prob (crude est): {needed_daily_prob:.4f}")
    print(f"   Suggested p_spontaneous_male: {needed_daily_prob / df['readiness_prob'].mean():.4f}")
    print()
    
    print("4. LARVAE PRODUCTION:")
    if total_larvae > 0:
        larvae_per_spawner = total_larvae / total_spawning_events if total_spawning_events > 0 else 0
        print(f"   Larvae per spawning event: {larvae_per_spawner:.0f}")
        print(f"   Total larvae in {n_test_days} days: {total_larvae}")
        
        # Extrapolate to full season
        season_larvae = total_larvae * (270 / n_test_days)
        print(f"   Projected season larvae: {season_larvae:.0f}")
        
        if season_larvae > 5000:  # Rough threshold
            print("   ❌ ISSUE: Larvae production may be too high")
            print("   RECOMMENDATION: Check larval cohort generation parameters")
    
    print()
    print("5. CASCADE EFFECTIVENESS:")
    within_range_frac = np.mean(distances <= spawning_config.cascade_radius)
    print(f"   Agents within cascade range: {within_range_frac:.1%}")
    if within_range_frac > 0.8:
        print("   ✅ Cascade radius appropriate for density")
    else:
        print("   ⚠️ WARNING: Cascade radius may be too small for this density")
    print()
    
    # SPECIFIC RECOMMENDATIONS
    print("=" * 60)
    print("SPECIFIC PARAMETER RECOMMENDATIONS")
    print("=" * 60)
    
    print("Based on diagnosis:")
    print()
    
    # Female spawning probability
    if current_bouts < target_bouts * 0.5:
        new_p_female = spawning_config.p_spontaneous_female * 3
        print(f"1. Increase p_spontaneous_female: {spawning_config.p_spontaneous_female} → {new_p_female:.3f}")
    
    # Male spawning probability  
    if current_bouts < target_bouts:
        multiplier = target_bouts / max(current_bouts, 0.1)
        new_p_male = spawning_config.p_spontaneous_male * multiplier
        print(f"2. Increase p_spontaneous_male: {spawning_config.p_spontaneous_male} → {new_p_male:.3f}")
    
    # Seasonal envelope
    if df['readiness_prob'].mean() < 0.6:
        print(f"3. Consider wider seasonal envelope: peak_width_days {spawning_config.peak_width_days} → {spawning_config.peak_width_days * 1.5}")
    
    print()
    print("Test completed successfully!")
    
    # Save results
    results_dir = Path("results/spawning_diagnosis")
    results_dir.mkdir(parents=True, exist_ok=True)
    
    df.to_csv(results_dir / "daily_spawning_analysis.csv", index=False)
    print(f"Results saved to {results_dir}/")


if __name__ == "__main__":
    diagnose_spawning_parameters()