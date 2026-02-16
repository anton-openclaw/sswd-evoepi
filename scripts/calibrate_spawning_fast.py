#!/usr/bin/env python3
"""
Fast spawning calibration for testing and initial exploration.

This is a simplified version that runs quickly to verify the approach works.
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


def evaluate_spawning_parameters_fast(p_spontaneous_female: float, p_spontaneous_male: float, 
                                      cascade_radius: float, peak_width_days: float,
                                      n_agents: int = 100, simulation_days: int = 90,
                                      node_latitude: float = 48.5, seed: int = 42) -> dict:
    """Fast evaluation with fewer agents and shorter simulation."""
    
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
    
    # Simulate season centered on peak spawning (DOY 105)
    start_doy = 60  # March 1, near spawning season
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
    major_bout_days = np.sum(daily_spawners >= 5)  # Lower threshold for small population
    
    # Male bout statistics
    male_agents = agents[agents['sex'] == 1]
    male_bout_counts = male_agents['has_spawned']
    mean_male_bouts = np.mean(male_bout_counts)
    males_that_spawned = np.sum(male_bout_counts > 0)
    
    # Female spawning statistics
    female_agents = agents[agents['sex'] == 0]
    females_that_spawned = np.sum(female_agents['has_spawned'] > 0)
    female_spawn_rate = females_that_spawned / len(female_agents) if len(female_agents) > 0 else 0
    
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
        'max_daily_spawners': np.max(daily_spawners) if len(daily_spawners) > 0 else 0,
        'mean_daily_spawners': np.mean(daily_spawners) if len(daily_spawners) > 0 else 0
    }


def run_fast_sweep():
    """Run fast parameter sweep for testing."""
    
    print("=" * 60)
    print("FAST SPAWNING CALIBRATION (Testing)")
    print("=" * 60)
    
    # Reduced parameter grids for speed
    param_grids = {
        'p_spontaneous_female': [0.02, 0.08, 0.15],
        'p_spontaneous_male': [0.03, 0.10, 0.20],
        'cascade_radius': [100, 500],
        'peak_width_days': [45, 90]
    }
    
    # Adjusted targets for smaller population
    targets = {
        'male_bouts': 1.0,  # Scale down for shorter simulation
        'major_bouts': 2.0,
        'larvae': 2000     # Scale down for smaller pop
    }
    
    print(f"Parameter ranges:")
    for param, values in param_grids.items():
        print(f"  {param}: {values}")
    
    print(f"\nTargets (scaled for n=100, 90 days):")
    for target, value in targets.items():
        print(f"  {target}: {value}")
    
    # Generate all parameter combinations
    param_names = list(param_grids.keys())
    param_combinations = list(itertools.product(*param_grids.values()))
    
    print(f"\nTotal combinations: {len(param_combinations)}")
    
    # Run sweep
    results = []
    start_time = time.time()
    
    for i, params in enumerate(param_combinations):
        param_dict = dict(zip(param_names, params))
        
        print(f"  Testing {i+1}/{len(param_combinations)}: {param_dict}")
        
        # Run evaluation
        result = evaluate_spawning_parameters_fast(**param_dict)
        
        # Simple score
        male_bout_error = abs(result['mean_male_bouts'] - targets['male_bouts'])
        major_bout_error = abs(result['major_bout_days'] - targets['major_bouts'])
        larvae_error = abs(result['total_larvae'] - targets['larvae']) / targets['larvae']
        
        score = male_bout_error + major_bout_error * 0.5 + larvae_error
        result['score'] = score
        
        print(f"    → male_bouts: {result['mean_male_bouts']:.2f}, "
              f"major_bouts: {result['major_bout_days']}, "
              f"larvae: {result['total_larvae']}, score: {score:.3f}")
        
        results.append(result)
    
    elapsed_time = time.time() - start_time
    print(f"\nFast sweep completed in {elapsed_time:.1f} seconds")
    
    # Sort by score and print results
    df = pd.DataFrame(results)
    df = df.sort_values('score')
    
    print("\n" + "=" * 80)
    print("RESULTS (Best to Worst)")
    print("=" * 80)
    
    for i, (idx, row) in enumerate(df.iterrows()):
        status = "GOOD" if row['score'] < 1.0 else "OK" if row['score'] < 2.0 else "BAD"
        print(f"{i+1:2d}. p_f={row['p_spontaneous_female']:.2f}, p_m={row['p_spontaneous_male']:.2f}, "
              f"r={row['cascade_radius']:3.0f}, w={row['peak_width_days']:2.0f} → "
              f"male_bouts={row['mean_male_bouts']:.2f}, major={row['major_bout_days']:2.0f}, "
              f"larvae={row['total_larvae']:4.0f}, score={row['score']:.3f} {status}")
    
    best = df.iloc[0]
    print(f"\n" + "=" * 50)
    print("BEST PARAMETER SET")
    print("=" * 50)
    print(f"p_spontaneous_female: {best['p_spontaneous_female']}")
    print(f"p_spontaneous_male: {best['p_spontaneous_male']}")
    print(f"cascade_radius: {best['cascade_radius']}")
    print(f"peak_width_days: {best['peak_width_days']}")
    print(f"Score: {best['score']:.3f}")
    
    return df.iloc[0].to_dict()


if __name__ == "__main__":
    best_params = run_fast_sweep()
    print(f"\nFast calibration complete! Best params: {best_params}")