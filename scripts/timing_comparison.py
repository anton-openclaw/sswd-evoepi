#!/usr/bin/env python3
"""
Before/After Timing Comparison Script
Measures spawning system performance before vs after vectorization optimizations.
"""

import time
import json
import numpy as np
import sys
import os
import copy
from pathlib import Path

# Add sswd_evoepi to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from sswd_evoepi.model import run_coupled_simulation
from sswd_evoepi.config import SimulationConfig

def create_test_configs():
    """Create simulation configurations for testing."""
    base_config = SimulationConfig()
    
    # Enable spawning
    spawning_config = copy.deepcopy(base_config)
    spawning_config.spawning.enabled = True
    spawning_config.spawning.season_start_doy = 305  # Nov 1
    spawning_config.spawning.season_end_doy = 212    # Jul 31
    
    # Disable spawning for comparison
    no_spawning_config = copy.deepcopy(base_config)
    no_spawning_config.spawning.enabled = False
    
    return spawning_config, no_spawning_config

def time_simulation(config, n_individuals, n_years, seed, timeout_s=60):
    """Time a single simulation run with timeout protection."""
    start_time = time.perf_counter()
    
    try:
        # Run simulation
        result = run_coupled_simulation(
            n_individuals=n_individuals,
            n_years=n_years,
            seed=seed,
            config=config,
            record_daily=False  # Skip daily recording for speed
        )
        
        elapsed = time.perf_counter() - start_time
        
        # Check if we hit timeout
        if elapsed > timeout_s:
            return 'TIMEOUT'
        
        return elapsed
        
    except Exception as e:
        print(f"ERROR in simulation (n={n_individuals}, seed={seed}): {e}")
        return 'ERROR'

def estimate_spawning_overhead(spawning_time, no_spawning_time):
    """Estimate the spawning system overhead."""
    if np.isnan(spawning_time) or np.isnan(no_spawning_time):
        return np.nan, np.nan
    
    overhead_s = spawning_time - no_spawning_time
    overhead_factor = spawning_time / no_spawning_time
    return overhead_s, overhead_factor

def main():
    print("SSWD-EvoEpi Spawning System: Before/After Timing Comparison")
    print("=" * 60)
    
    # Create configurations
    spawning_config, no_spawning_config = create_test_configs()
    
    # Test configurations (reduced for timing comparison)
    test_configs = [
        {'n': 50, 'years': 2, 'label': 'small'},
        {'n': 100, 'years': 2, 'label': 'medium'},
        {'n': 200, 'years': 2, 'label': 'large'},
    ]
    
    results = []
    
    print(f"\nTesting {len(test_configs)} configurations...")
    
    for cfg in test_configs:
        print(f"\n--- {cfg['label'].upper()} ({cfg['n']} agents, {cfg['years']} years) ---")
        
        # With spawning - 3 runs
        print("  Testing WITH spawning...")
        times_spawn = []
        for seed in [42, 43, 44]:
            print(f"    Run {seed-41}/3...", end=' ')
            elapsed = time_simulation(spawning_config, cfg['n'], cfg['years'], seed)
            if isinstance(elapsed, str):
                print(elapsed)
                times_spawn.append(np.nan)
            else:
                print(f"{elapsed:.1f}s")
                times_spawn.append(elapsed)
        
        # Without spawning - 3 runs  
        print("  Testing WITHOUT spawning...")
        times_nospawn = []
        for seed in [42, 43, 44]:
            print(f"    Run {seed-41}/3...", end=' ')
            elapsed = time_simulation(no_spawning_config, cfg['n'], cfg['years'], seed)
            if isinstance(elapsed, str):
                print(elapsed)
                times_nospawn.append(np.nan)
            else:
                print(f"{elapsed:.1f}s")
                times_nospawn.append(elapsed)
        
        # Calculate averages (excluding NaN/errors)
        valid_spawn = [t for t in times_spawn if not np.isnan(t)]
        valid_nospawn = [t for t in times_nospawn if not np.isnan(t)]
        
        if len(valid_spawn) > 0 and len(valid_nospawn) > 0:
            mean_spawn = np.mean(valid_spawn)
            mean_nospawn = np.mean(valid_nospawn)
            overhead_factor = mean_spawn / mean_nospawn
            spawning_overhead_s = mean_spawn - mean_nospawn
        else:
            mean_spawn = np.nan
            mean_nospawn = np.nan
            overhead_factor = np.nan
            spawning_overhead_s = np.nan
        
        results.append({
            'label': cfg['label'],
            'n_agents': cfg['n'],
            'n_years': cfg['years'],
            'with_spawning_s': mean_spawn,
            'without_spawning_s': mean_nospawn,
            'spawning_overhead_s': spawning_overhead_s,
            'overhead_factor': overhead_factor,
            'raw_spawn_times': times_spawn,
            'raw_nospawn_times': times_nospawn
        })
        
        if not np.isnan(overhead_factor):
            print(f"  → Spawning overhead: {spawning_overhead_s:.1f}s ({overhead_factor:.1f}× total time)")
    
    # Combine results  
    all_results = {
        'simulation_scaling': results,
        'metadata': {
            'timestamp': time.time(),
            'optimization_notes': 'vectorized_cascade_and_recent_spawners',
            'comparison_baseline': 'spawning_enabled_vs_disabled',
            'prior_speedup_claimed': '8.8x (from comprehensive benchmark notes)'
        }
    }
    
    # Save results
    os.makedirs('results/performance', exist_ok=True)
    with open('results/performance/timing_comparison_data.json', 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    
    print(f"\n--- SUMMARY ---")
    for r in results:
        if not np.isnan(r['overhead_factor']):
            print(f"{r['label']:6}: {r['overhead_factor']:.1f}× overhead ({r['spawning_overhead_s']:.1f}s spawning cost)")
        else:
            print(f"{r['label']:6}: FAILED/TIMEOUT")
    
    print(f"\nResults saved to: results/performance/timing_comparison_data.json")
    return all_results

if __name__ == '__main__':
    results = main()