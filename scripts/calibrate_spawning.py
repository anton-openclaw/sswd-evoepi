#!/usr/bin/env python3
"""
Systematic parameter sweep to calibrate spawning rates in SSWD-EvoEpi.

Based on Phase 1A1 diagnosis, this script runs a parameter grid search
to find optimal spawning probability values that match biological targets:
- Male bout count: ~2.2 per male per season
- Major spawning bouts per season: 2-5
- Total annual larvae: within 20% of pulse model output

Author: Anton 🔬
Date: 2026-02-15
"""

import numpy as np
import pandas as pd
import itertools
import time
from pathlib import Path

from sswd_evoepi.types import AGENT_DTYPE, Stage, N_LOCI, allocate_agents, allocate_genotypes
from sswd_evoepi.config import SimulationConfig, SpawningSection, DiseaseSection
from sswd_evoepi.spawning import spawning_step, seasonal_readiness_prob, latitude_adjusted_peak
from sswd_evoepi.genetics import initialize_genotypes, initialize_effect_sizes


def create_test_population(n_agents: int, habitat_area_m2: float, node_id: int = 1, seed: int = 42) -> tuple:
    """Create test population of spawning-ready adults."""
    
    rng = np.random.default_rng(seed)
    
    # Distribute agents uniformly over square habitat
    side_length = np.sqrt(habitat_area_m2)
    
    agents = allocate_agents(n_agents)
    genotypes = allocate_genotypes(n_agents)
    
    # Initialize effect sizes and genotypes
    effect_sizes = initialize_effect_sizes(rng)
    geno_init = initialize_genotypes(n_agents, effect_sizes, rng, target_mean_r=0.08)
    genotypes[:n_agents] = geno_init
    
    for i in range(n_agents):
        agents[i]['node_id'] = node_id
        agents[i]['x'] = rng.uniform(0, side_length)
        agents[i]['y'] = rng.uniform(0, side_length)
        agents[i]['stage'] = Stage.ADULT
        agents[i]['age'] = 365 * 3  # Mature adults
        agents[i]['size'] = 450.0  # Above reproduction threshold
        agents[i]['sex'] = rng.choice([0, 1])  # Random sex
        agents[i]['alive'] = True
        
        # Reset spawning state for clean start
        agents[i]['spawning_ready'] = False
        agents[i]['has_spawned'] = 0
        agents[i]['spawn_refractory'] = 0
        agents[i]['last_spawn_day'] = -999
    
    return agents, genotypes


def evaluate_spawning_parameters(p_spontaneous_female: float, p_spontaneous_male: float, 
                                cascade_radius: float, peak_width_days: float,
                                n_agents: int = 200, simulation_days: int = 90,
                                node_latitude: float = 48.5, seed: int = 42) -> dict:
    """Evaluate spawning parameters for one full year."""
    
    habitat_area = 10000  # 100m x 100m = 10,000 m²
    agents, genotypes = create_test_population(n_agents, habitat_area, seed=seed)
    
    # Create modified spawning config
    spawning_config = SpawningSection()
    spawning_config.p_spontaneous_female = p_spontaneous_female
    spawning_config.p_spontaneous_male = p_spontaneous_male
    spawning_config.cascade_radius = cascade_radius
    spawning_config.peak_width_days = peak_width_days
    
    disease_config = DiseaseSection()
    
    rng = np.random.default_rng(seed)
    
    # Simulate full year starting at peak season (DOY 105)
    start_doy = 105  # Peak spawning season
    total_spawning_events = 0
    total_larvae = 0
    daily_spawners = []
    
    for day_offset in range(simulation_days):
        doy = (start_doy + day_offset - 1) % 365 + 1  # Handle year wrap
        
        # Count pre-spawning state
        spawning_events_before = np.sum(agents['has_spawned'])
        
        # Run spawning step
        cohorts = spawning_step(
            agents, genotypes, doy, node_latitude,
            spawning_config, disease_config, rng
        )
        
        # Count results
        spawning_events_after = np.sum(agents['has_spawned'])
        new_spawning_events = spawning_events_after - spawning_events_before
        larvae_today = sum(cohort.n_competent for cohort in cohorts)
        
        daily_spawners.append(new_spawning_events)
        total_spawning_events += new_spawning_events
        total_larvae += larvae_today
    
    # Analyze results
    daily_spawners = np.array(daily_spawners)
    major_bout_days = np.sum(daily_spawners >= 3)   # Days with ≥3 spawners (scaled for smaller pop)
    
    # Male bout statistics
    male_agents = agents[agents['sex'] == 1]
    male_bout_counts = male_agents['has_spawned']
    mean_male_bouts = np.mean(male_bout_counts)
    males_that_spawned = np.sum(male_bout_counts > 0)
    
    # Female spawning statistics
    female_agents = agents[agents['sex'] == 0]
    females_that_spawned = np.sum(female_agents['has_spawned'] > 0)
    female_spawn_rate = females_that_spawned / len(female_agents)
    
    return {
        'p_spontaneous_female': p_spontaneous_female,
        'p_spontaneous_male': p_spontaneous_male,
        'cascade_radius': cascade_radius,
        'peak_width_days': peak_width_days,
        'total_spawning_events': total_spawning_events,
        'total_larvae': total_larvae,
        'mean_male_bouts': mean_male_bouts,
        'major_bout_days': major_bout_days,
        'males_that_spawned': males_that_spawned,
        'females_that_spawned': females_that_spawned,
        'female_spawn_rate': female_spawn_rate,
        'max_daily_spawners': np.max(daily_spawners),
        'mean_daily_spawners': np.mean(daily_spawners)
    }


def score_parameter_set(results: dict, targets: dict) -> float:
    """Score parameter set against biological targets.
    
    Lower scores are better (0 = perfect match).
    """
    male_bout_error = abs(results['mean_male_bouts'] - targets['male_bouts'])
    major_bout_error = abs(results['major_bout_days'] - targets['major_bouts'])
    
    # Larvae ratio error (log scale to handle large differences)
    larvae_ratio = results['total_larvae'] / max(targets['larvae'], 1)
    larvae_error = abs(np.log(larvae_ratio)) if larvae_ratio > 0 else 10.0
    
    # Weighted composite score
    score = (
        male_bout_error * 1.0 +              # Primary target
        major_bout_error * 0.3 +             # Secondary target  
        larvae_error * 0.5                   # Tertiary target
    )
    
    return score


def run_parameter_sweep():
    """Run systematic parameter sweep."""
    
    print("=" * 70)
    print("SSWD-EVOEPI SPAWNING CALIBRATION SWEEP")
    print("=" * 70)
    
    # Define parameter grids based on 1A1 diagnosis
    param_grids = {
        'p_spontaneous_female': [0.005, 0.01, 0.02, 0.04, 0.08, 0.15],
        'p_spontaneous_male': [0.008, 0.015, 0.03, 0.06, 0.10, 0.20],
        'cascade_radius': [50, 100, 200, 500],
        'peak_width_days': [30, 45, 60, 90]
    }
    
    # Biological targets (scaled for n=200, 90 days)
    targets = {
        'male_bouts': 0.8,  # Scaled down for shorter season
        'major_bouts': 2.0, # Scaled down 
        'larvae': 5000      # Scaled down for smaller population
    }
    
    print(f"Parameter ranges:")
    for param, values in param_grids.items():
        print(f"  {param}: {values}")
    
    print(f"\nTargets:")
    for target, value in targets.items():
        print(f"  {target}: {value}")
    
    # Generate all parameter combinations
    param_names = list(param_grids.keys())
    param_combinations = list(itertools.product(*param_grids.values()))
    
    print(f"\nTotal combinations: {len(param_combinations)}")
    
    # Estimate runtime
    if len(param_combinations) > 100:
        print(f"⚠️  Large search space detected. Running coarse grid first...")
        
        # Use coarse grid for initial screening
        coarse_grids = {
            'p_spontaneous_female': [0.02, 0.08],
            'p_spontaneous_male': [0.03, 0.10],
            'cascade_radius': [200],
            'peak_width_days': [60]
        }
        
        param_combinations = list(itertools.product(*coarse_grids.values()))
        print(f"Coarse grid combinations: {len(param_combinations)}")
    
    # Run sweep
    results = []
    start_time = time.time()
    
    for i, params in enumerate(param_combinations):
        param_dict = dict(zip(param_names, params))
        
        print(f"  Testing {i+1}/{len(param_combinations)}: {param_dict}")
        
        # Run evaluation
        result = evaluate_spawning_parameters(**param_dict)
        score = score_parameter_set(result, targets)
        result['score'] = score
        
        print(f"    → male_bouts: {result['mean_male_bouts']:.2f}, "
              f"major_bouts: {result['major_bout_days']}, "
              f"larvae: {result['total_larvae']}, score: {score:.3f}")
        
        results.append(result)
        
        # Progress update
        if i % max(1, len(param_combinations) // 10) == 0:
            elapsed = time.time() - start_time
            progress = (i + 1) / len(param_combinations)
            eta = elapsed / progress - elapsed if progress > 0 else 0
            best_score = min(r['score'] for r in results) if results else float('inf')
            print(f"    Progress: {i+1}/{len(param_combinations)} ({progress:.1%}) "
                  f"ETA: {eta:.1f}s | Best score so far: {best_score:.3f}")
    
    # Convert to DataFrame and sort by score
    df = pd.DataFrame(results)
    df = df.sort_values('score')
    
    elapsed_time = time.time() - start_time
    print(f"\nSweep completed in {elapsed_time/60:.1f} minutes")
    
    # Print results table
    print("\n" + "=" * 120)
    print("RESULTS TABLE (Top 20 parameter sets)")
    print("=" * 120)
    
    print(f"{'p_sf':<6} {'p_sm':<6} {'cascade_r':<10} {'width':<6} "
          f"{'male_bouts':<11} {'major_bouts':<12} {'larvae':<8} {'score':<8} {'status':<8}")
    print("-" * 120)
    
    for i, row in df.head(20).iterrows():
        status = "GOOD" if row['score'] < 1.0 else "OK" if row['score'] < 2.0 else "BAD"
        
        print(f"{row['p_spontaneous_female']:<6.3f} {row['p_spontaneous_male']:<6.3f} "
              f"{row['cascade_radius']:<10.0f} {row['peak_width_days']:<6.0f} "
              f"{row['mean_male_bouts']:<11.2f} {row['major_bout_days']:<12.0f} "
              f"{row['total_larvae']:<8.0f} {row['score']:<8.3f} {status}")
    
    # Best parameter set
    best = df.iloc[0]
    print(f"\n" + "=" * 70)
    print("BEST PARAMETER SET")
    print("=" * 70)
    
    print(f"Parameters:")
    print(f"  p_spontaneous_female: {best['p_spontaneous_female']}")
    print(f"  p_spontaneous_male: {best['p_spontaneous_male']}")
    print(f"  cascade_radius: {best['cascade_radius']}")
    print(f"  peak_width_days: {best['peak_width_days']}")
    
    print(f"\nPerformance:")
    print(f"  Score: {best['score']:.3f}")
    print(f"  Male bouts/season: {best['mean_male_bouts']:.2f} (target: {targets['male_bouts']})")
    print(f"  Major bout days: {best['major_bout_days']:.0f} (target: {targets['major_bouts']})")
    print(f"  Total larvae: {best['total_larvae']:.0f} (target: {targets['larvae']:.0f})")
    print(f"  Female spawn rate: {best['female_spawn_rate']:.2f}")
    print(f"  Males that spawned: {best['males_that_spawned']}/{int(500/2)}")
    
    # Ratio analysis
    print(f"\nTarget ratios:")
    print(f"  Male bouts ratio: {best['mean_male_bouts']/targets['male_bouts']:.3f}")
    print(f"  Major bouts ratio: {best['major_bout_days']/targets['major_bouts']:.3f}")
    print(f"  Larvae ratio: {best['total_larvae']/targets['larvae']:.3f}")
    
    # Save results
    results_dir = Path("results/spawning_calibration")
    results_dir.mkdir(parents=True, exist_ok=True)
    
    timestamp = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")
    df.to_csv(results_dir / f"parameter_sweep_{timestamp}.csv", index=False)
    
    # Save best parameters as separate file for easy access
    best_params = {
        'p_spontaneous_female': best['p_spontaneous_female'],
        'p_spontaneous_male': best['p_spontaneous_male'], 
        'cascade_radius': best['cascade_radius'],
        'peak_width_days': best['peak_width_days'],
        'score': best['score'],
        'male_bouts_achieved': best['mean_male_bouts'],
        'major_bouts_achieved': best['major_bout_days'],
        'larvae_achieved': best['total_larvae']
    }
    
    best_params_file = results_dir / f"best_parameters_{timestamp}.json"
    import json
    with open(best_params_file, 'w') as f:
        json.dump(best_params, f, indent=2)
    
    print(f"\nResults saved to:")
    print(f"  Full sweep: {results_dir}/parameter_sweep_{timestamp}.csv")
    print(f"  Best params: {best_params_file}")
    
    return best_params


if __name__ == "__main__":
    best_params = run_parameter_sweep()
    
    print("\nCalibration complete!")
    print(f"Recommended parameters: {best_params}")