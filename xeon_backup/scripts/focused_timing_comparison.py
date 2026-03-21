#!/usr/bin/env python3
"""
Focused Before/After Timing Comparison
Tests spawning system performance at different scales.
"""

import time
import json
import numpy as np
import sys
import copy
from pathlib import Path

# Add sswd_evoepi to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from sswd_evoepi.model import run_coupled_simulation
from sswd_evoepi.config import SimulationConfig

def time_simulation_safely(config, n_individuals, n_years, seed, timeout_s=120):
    """Time a single simulation run with timeout and error handling."""
    start_time = time.perf_counter()
    
    try:
        result = run_coupled_simulation(
            n_individuals=n_individuals,
            n_years=n_years,
            seed=seed,
            config=config,
            record_daily=False
        )
        
        elapsed = time.perf_counter() - start_time
        
        if elapsed > timeout_s:
            return {'status': 'timeout', 'time': elapsed}
        
        return {
            'status': 'success', 
            'time': elapsed, 
            'final_pop': result.final_pop,
            'initial_pop': result.initial_pop,
            'total_disease_deaths': result.total_disease_deaths
        }
        
    except Exception as e:
        elapsed = time.perf_counter() - start_time
        return {'status': 'error', 'time': elapsed, 'error': str(e)}

def main():
    print("SSWD-EvoEpi: Focused Spawning System Timing Comparison")
    print("=" * 60)
    
    # Create configurations
    base_config = SimulationConfig()
    
    # With spawning enabled
    spawning_config = copy.deepcopy(base_config)
    spawning_config.spawning.enabled = True
    spawning_config.spawning.season_start_doy = 305  # Nov 1
    spawning_config.spawning.season_end_doy = 212    # Jul 31
    
    # Without spawning  
    no_spawning_config = copy.deepcopy(base_config)
    no_spawning_config.spawning.enabled = False
    
    # Test configurations - focusing on sizes where spawning matters
    test_configs = [
        {'n': 50, 'years': 1, 'label': 'small', 'seeds': [42, 43]},
        {'n': 100, 'years': 1, 'label': 'medium', 'seeds': [42, 43]},  
        {'n': 200, 'years': 1, 'label': 'large', 'seeds': [42]},  # Just 1 run for large
        {'n': 500, 'years': 1, 'label': 'xlarge', 'seeds': [42]}, # Just 1 run for xlarge
    ]
    
    results = []
    
    for cfg in test_configs:
        print(f"\n--- {cfg['label'].upper()} ({cfg['n']} agents, {cfg['years']} year) ---")
        
        spawn_times = []
        nospawn_times = []
        spawn_results = []
        nospawn_results = []
        
        # Test WITH spawning
        print("  With spawning:", end=' ')
        for seed in cfg['seeds']:
            print(f"{seed}...", end=' ')
            result = time_simulation_safely(spawning_config, cfg['n'], cfg['years'], seed)
            spawn_times.append(result['time'])
            spawn_results.append(result)
            if result['status'] != 'success':
                print(f"[{result['status']}]", end=' ')
        print()
        
        # Test WITHOUT spawning
        print("  Without spawning:", end=' ')
        for seed in cfg['seeds']:
            print(f"{seed}...", end=' ')
            result = time_simulation_safely(no_spawning_config, cfg['n'], cfg['years'], seed)
            nospawn_times.append(result['time'])
            nospawn_results.append(result)
            if result['status'] != 'success':
                print(f"[{result['status']}]", end=' ')
        print()
        
        # Calculate stats
        spawn_success = [r for r in spawn_results if r['status'] == 'success']
        nospawn_success = [r for r in nospawn_results if r['status'] == 'success']
        
        if len(spawn_success) > 0 and len(nospawn_success) > 0:
            spawn_mean = np.mean([r['time'] for r in spawn_success])
            nospawn_mean = np.mean([r['time'] for r in nospawn_success])
            overhead_s = spawn_mean - nospawn_mean
            overhead_factor = spawn_mean / nospawn_mean
            
            print(f"  → Spawning: {spawn_mean:.2f}s, No spawning: {nospawn_mean:.2f}s")
            print(f"  → Overhead: {overhead_s:+.2f}s ({overhead_factor:.2f}× total time)")
            
            # Population comparison
            spawn_pop = np.mean([r['final_pop'] for r in spawn_success])
            nospawn_pop = np.mean([r['final_pop'] for r in nospawn_success])
            print(f"  → Final pop: {spawn_pop:.0f} (spawn) vs {nospawn_pop:.0f} (no spawn)")
        else:
            spawn_mean = nospawn_mean = overhead_s = overhead_factor = np.nan
            spawn_pop = nospawn_pop = np.nan
            print(f"  → FAILED: {len(spawn_success)}/{len(cfg['seeds'])} spawning, {len(nospawn_success)}/{len(cfg['seeds'])} no-spawning")
        
        results.append({
            'label': cfg['label'],
            'n_agents': cfg['n'],
            'n_years': cfg['years'],
            'n_runs': len(cfg['seeds']),
            'spawning_time_s': spawn_mean,
            'no_spawning_time_s': nospawn_mean,
            'overhead_s': overhead_s,
            'overhead_factor': overhead_factor,
            'spawning_final_pop': spawn_pop,
            'no_spawning_final_pop': nospawn_pop,
            'raw_spawning_results': spawn_results,
            'raw_no_spawning_results': nospawn_results
        })
    
    # Historical comparison with pre-optimization data
    print(f"\n--- HISTORICAL COMPARISON ---")
    print("Pre-optimization baselines from MEMORY.md:")
    print("  - 50 agents, 2 years = 22.2s")
    print("  - 100 agents, 5 years = 78s")
    print()
    print("Current results (1 year equivalents):")
    for r in results:
        if not np.isnan(r['spawning_time_s']):
            # Scale to 2-year equivalent for 50 agents comparison
            if r['n_agents'] == 50:
                scaled_time = r['spawning_time_s'] * 2  # Scale to 2 years
                speedup = 22.2 / scaled_time
                print(f"  - {r['n_agents']} agents: {r['spawning_time_s']:.1f}s/yr → {scaled_time:.1f}s/2yr (vs 22.2s baseline = {speedup:.1f}× faster)")
    
    # Overall analysis
    print(f"\n--- SPAWNING SYSTEM ANALYSIS ---")
    valid_results = [r for r in results if not np.isnan(r['overhead_factor'])]
    
    if len(valid_results) > 0:
        mean_overhead = np.mean([r['overhead_s'] for r in valid_results])
        mean_factor = np.mean([r['overhead_factor'] for r in valid_results])
        
        print(f"Average spawning overhead: {mean_overhead:+.2f}s ({mean_factor:.2f}× factor)")
        print(f"Overhead scales with population size:")
        for r in valid_results:
            overhead_per_agent = r['overhead_s'] / r['n_agents'] * 1000  # ms per agent
            print(f"  - {r['n_agents']} agents: {overhead_per_agent:.1f} ms/agent")
        
        # Check if population dynamics differ
        pop_diffs = []
        for r in valid_results:
            if not np.isnan(r['spawning_final_pop']) and not np.isnan(r['no_spawning_final_pop']):
                diff_pct = 100 * (r['spawning_final_pop'] - r['no_spawning_final_pop']) / r['no_spawning_final_pop']
                pop_diffs.append(diff_pct)
                print(f"  - {r['label']}: {diff_pct:+.1f}% population difference")
        
        if len(pop_diffs) > 0:
            avg_pop_diff = np.mean(pop_diffs)
            print(f"Average population impact: {avg_pop_diff:+.1f}%")
        
    else:
        print("No successful comparisons completed.")
    
    # Save results
    import os
    os.makedirs('results/performance', exist_ok=True)
    
    output = {
        'comparison_results': results,
        'summary': {
            'mean_overhead_s': mean_overhead if len(valid_results) > 0 else None,
            'mean_overhead_factor': mean_factor if len(valid_results) > 0 else None,
            'successful_comparisons': len(valid_results),
            'total_comparisons': len(results)
        },
        'metadata': {
            'timestamp': time.time(),
            'optimization': 'vectorized_cascade_and_recent_spawners_checks',
            'test_type': 'spawning_enabled_vs_disabled',
            'notes': 'Post-optimization timing comparison focusing on spawning system overhead'
        }
    }
    
    with open('results/performance/focused_timing_comparison.json', 'w') as f:
        json.dump(output, f, indent=2, default=str)
    
    print(f"\nResults saved to: results/performance/focused_timing_comparison.json")
    return output

if __name__ == '__main__':
    results = main()