#!/usr/bin/env python3
"""
Simple benchmark of key SSWD-EvoEpi components.

This script measures the performance of critical model components at different scales
to identify bottlenecks and scaling behavior.
"""

import time
import numpy as np
import pandas as pd
import json
from pathlib import Path
from typing import Dict, List, Tuple
import matplotlib.pyplot as plt

# Add project root to path
import sys
sys.path.append(str(Path(__file__).parent.parent))

from sswd_evoepi.config import SimulationConfig, load_config, validate_config
from sswd_evoepi.model import make_effect_sizes, initialize_population, run_coupled_simulation
from sswd_evoepi.rng import create_rng_hierarchy


def time_function(func, *args, n_iterations=100, **kwargs):
    """Time a function call over multiple iterations"""
    start_time = time.perf_counter()
    
    for _ in range(n_iterations):
        try:
            func(*args, **kwargs)
        except Exception as e:
            print(f"  Error in {func.__name__}: {e}")
            return None
    
    end_time = time.perf_counter()
    return (end_time - start_time) / n_iterations


def benchmark_population_initialization(config: SimulationConfig, population_sizes: List[int]) -> Dict[str, List[float]]:
    """Benchmark population initialization at different scales"""
    print("🧪 Benchmarking population initialization...")
    
    results = {'n_agents': [], 'init_time': []}
    effect_sizes = make_effect_sizes(seed=42)
    
    for n_agents in population_sizes:
        print(f"  Testing {n_agents} agents...")
        
        rng_hierarchy = create_rng_hierarchy(master_seed=42, n_nodes=1)
        rng = rng_hierarchy['global']
        
        start_time = time.perf_counter()
        
        try:
            agents, genotypes = initialize_population(
                n_individuals=n_agents,
                max_agents=n_agents + 100,
                habitat_area=100_000_000,  # 100 km²
                effect_sizes=effect_sizes,
                pop_cfg=config.population,
                rng=rng
            )
            
            init_time = time.perf_counter() - start_time
            print(f"    {init_time:.4f}s")
            
        except Exception as e:
            print(f"    ERROR: {e}")
            init_time = float('nan')
        
        results['n_agents'].append(n_agents)
        results['init_time'].append(init_time)
    
    return results


def benchmark_full_simulation_scaling(config: SimulationConfig, population_sizes: List[int]) -> Dict[str, List[float]]:
    """Benchmark full simulations at different scales"""
    print("🔬 Benchmarking full simulation scaling...")
    
    results = {
        'n_agents': [],
        'time_2yr': [],
        'memory_mb': []
    }
    
    for n_agents in population_sizes:
        print(f"  Testing {n_agents} agents (2-year simulation)...")
        
        start_time = time.perf_counter()
        
        try:
            # Run a short 2-year simulation
            sim_results = run_coupled_simulation(
                n_individuals=n_agents,
                carrying_capacity=n_agents,
                n_years=2,
                config=config
            )
            
            simulation_time = time.perf_counter() - start_time
            
            # Estimate memory usage (rough)
            memory_mb = (n_agents * 250) / (1024 * 1024)  # ~250 bytes per agent
            
            print(f"    {simulation_time:.2f}s ({memory_mb:.2f} MB)")
            
        except Exception as e:
            print(f"    ERROR: {e}")
            simulation_time = float('nan')
            memory_mb = float('nan')
        
        results['n_agents'].append(n_agents)
        results['time_2yr'].append(simulation_time)
        results['memory_mb'].append(memory_mb)
    
    return results


def estimate_scaling_exponent(n_agents: List[int], times: List[float]) -> float:
    """Estimate scaling exponent from timing data"""
    valid_pairs = [(n, t) for n, t in zip(n_agents, times) if not np.isnan(t) and t > 0]
    
    if len(valid_pairs) < 2:
        return float('nan')
    
    log_n = np.log([n for n, t in valid_pairs])
    log_t = np.log([t for n, t in valid_pairs])
    
    # Linear regression in log space: log(t) = α + β*log(n)
    coeffs = np.polyfit(log_n, log_t, 1)
    return coeffs[0]  # β is the scaling exponent


def create_benchmark_report(init_results: Dict, sim_results: Dict, output_dir: Path):
    """Create benchmark report"""
    
    report_path = output_dir / 'SIMPLE_BENCHMARK_REPORT.md'
    
    with open(report_path, 'w') as f:
        f.write("# SSWD-EvoEpi Simple Benchmark Report\n\n")
        f.write(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        # Population initialization
        f.write("## Population Initialization Performance\n\n")
        f.write("| Population Size | Initialization Time (s) | Notes |\n")
        f.write("|-----------------|-------------------------|-------|\n")
        
        for n, t in zip(init_results['n_agents'], init_results['init_time']):
            if np.isnan(t):
                f.write(f"| {n:,} | ERROR | Failed |\n")
            else:
                f.write(f"| {n:,} | {t:.4f} | |\n")
        
        # Scaling exponent for initialization
        init_exp = estimate_scaling_exponent(init_results['n_agents'], init_results['init_time'])
        f.write(f"\n**Initialization scaling:** O(N^{init_exp:.2f})\n\n")
        
        # Full simulation performance
        f.write("## Full Simulation Performance (2 years)\n\n")
        f.write("| Population Size | Runtime (s) | Memory (MB) | Agents/second | Notes |\n")
        f.write("|-----------------|-------------|-------------|---------------|-------|\n")
        
        for i, n in enumerate(sim_results['n_agents']):
            t = sim_results['time_2yr'][i]
            mem = sim_results['memory_mb'][i]
            
            if np.isnan(t):
                f.write(f"| {n:,} | ERROR | {mem:.2f} | - | Failed |\n")
            else:
                agents_per_sec = n / t if t > 0 else 0
                f.write(f"| {n:,} | {t:.2f} | {mem:.2f} | {agents_per_sec:.0f} | |\n")
        
        # Scaling exponent for simulation
        sim_exp = estimate_scaling_exponent(sim_results['n_agents'], sim_results['time_2yr'])
        f.write(f"\n**Simulation scaling:** O(N^{sim_exp:.2f})\n\n")
        
        # Performance projections
        f.write("## Performance Projections\n\n")
        f.write("Based on scaling analysis, estimated runtimes for larger populations:\n\n")
        f.write("| Population Size | Estimated 20-year Runtime | Memory Usage |\n")
        f.write("|-----------------|---------------------------|-------------|\n")
        
        # Use the last successful measurement as a reference
        valid_sims = [(n, t) for n, t in zip(sim_results['n_agents'], sim_results['time_2yr']) 
                      if not np.isnan(t) and t > 0]
        
        if valid_sims and not np.isnan(sim_exp):
            ref_n, ref_t = valid_sims[-1]  # Use largest successful case
            
            for target_n in [1000, 5000, 10000, 50000]:
                # Scale time based on exponent
                projected_2yr = ref_t * (target_n / ref_n) ** sim_exp
                projected_20yr = projected_2yr * 10  # 10x for 20 years vs 2 years
                
                # Memory scaling is linear
                memory_mb = (target_n * 250) / (1024 * 1024)
                memory_gb = memory_mb / 1024
                
                if projected_20yr < 3600:  # Less than 1 hour
                    time_str = f"{projected_20yr:.0f}s"
                elif projected_20yr < 86400:  # Less than 1 day
                    time_str = f"{projected_20yr/3600:.1f}h"
                else:
                    time_str = f"{projected_20yr/86400:.1f}d"
                
                if memory_gb < 1:
                    mem_str = f"{memory_mb:.0f} MB"
                else:
                    mem_str = f"{memory_gb:.1f} GB"
                
                f.write(f"| {target_n:,} | {time_str} | {mem_str} |\n")
        else:
            f.write("| N/A | Insufficient data | N/A |\n")
        
        f.write("\n## Key Findings\n\n")
        
        if not np.isnan(sim_exp):
            if sim_exp < 1.5:
                complexity = "nearly linear (very good)"
            elif sim_exp < 2.0:
                complexity = "super-linear (acceptable)"
            elif sim_exp < 2.5:
                complexity = "quadratic-like (concerning)"
            else:
                complexity = "worse than quadratic (problematic)"
                
            f.write(f"- **Scaling behavior**: {complexity} - O(N^{sim_exp:.2f})\n")
        else:
            f.write("- **Scaling behavior**: Could not determine from available data\n")
        
        # Find most successful population size
        max_successful = 0
        for n, t in zip(sim_results['n_agents'], sim_results['time_2yr']):
            if not np.isnan(t):
                max_successful = max(max_successful, n)
        
        f.write(f"- **Maximum tested population**: {max_successful:,} agents\n")
        
        if max_successful >= 500:
            f.write("- **Assessment**: Model handles moderate population sizes well\n")
        elif max_successful >= 100:
            f.write("- **Assessment**: Model limited to small population sizes\n")
        else:
            f.write("- **Assessment**: Model has significant scalability issues\n")


def create_benchmark_visualization(init_results: Dict, sim_results: Dict, output_dir: Path):
    """Create benchmark visualization"""
    
    plt.style.use('dark_background')
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    
    # 1. Initialization performance
    n_agents = init_results['n_agents']
    init_times = init_results['init_time']
    
    # Filter out NaN values
    valid_init = [(n, t) for n, t in zip(n_agents, init_times) if not np.isnan(t)]
    
    if valid_init:
        n_init, t_init = zip(*valid_init)
        ax1.loglog(n_init, t_init, 'o-', color='#4ecdc4', linewidth=2, markersize=8)
        ax1.set_xlabel('Population Size')
        ax1.set_ylabel('Initialization Time (s)')
        ax1.set_title('Population Initialization Performance')
        ax1.grid(True, alpha=0.3)
    else:
        ax1.text(0.5, 0.5, 'No successful initialization data', ha='center', va='center', transform=ax1.transAxes)
    
    # 2. Simulation performance
    sim_n = sim_results['n_agents']
    sim_times = sim_results['time_2yr']
    
    valid_sim = [(n, t) for n, t in zip(sim_n, sim_times) if not np.isnan(t)]
    
    if valid_sim:
        n_sim, t_sim = zip(*valid_sim)
        ax2.loglog(n_sim, t_sim, 's-', color='#ff6b6b', linewidth=2, markersize=8)
        ax2.set_xlabel('Population Size')
        ax2.set_ylabel('2-Year Runtime (s)')
        ax2.set_title('Full Simulation Performance')
        ax2.grid(True, alpha=0.3)
    else:
        ax2.text(0.5, 0.5, 'No successful simulation data', ha='center', va='center', transform=ax2.transAxes)
    
    # 3. Memory usage
    memory_data = sim_results['memory_mb']
    valid_mem = [(n, m) for n, m in zip(sim_n, memory_data) if not np.isnan(m)]
    
    if valid_mem:
        n_mem, mem = zip(*valid_mem)
        ax3.loglog(n_mem, mem, 'd-', color='#feca57', linewidth=2, markersize=8)
        ax3.set_xlabel('Population Size')
        ax3.set_ylabel('Memory Usage (MB)')
        ax3.set_title('Memory Usage Scaling')
        ax3.grid(True, alpha=0.3)
        
        # Add reference lines
        ax3.axhline(y=1024, color='red', linestyle='--', alpha=0.7, label='1 GB')
        ax3.axhline(y=4096, color='orange', linestyle='--', alpha=0.7, label='4 GB')
        ax3.legend()
    else:
        ax3.text(0.5, 0.5, 'No memory data', ha='center', va='center', transform=ax3.transAxes)
    
    # 4. Throughput (agents per second)
    if valid_sim:
        throughput = [n / t for n, t in valid_sim]
        ax4.semilogx([n for n, t in valid_sim], throughput, '^-', 
                    color='#96ceb4', linewidth=2, markersize=8)
        ax4.set_xlabel('Population Size')
        ax4.set_ylabel('Throughput (agents/second)')
        ax4.set_title('Processing Throughput')
        ax4.grid(True, alpha=0.3)
    else:
        ax4.text(0.5, 0.5, 'No throughput data', ha='center', va='center', transform=ax4.transAxes)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'simple_benchmark.png', dpi=300, bbox_inches='tight',
                facecolor='black', edgecolor='none')
    plt.close()


def main():
    """Main benchmark execution"""
    print("🔬 Starting SSWD-EvoEpi Simple Benchmark...")
    
    # Load configuration
    config = load_config('configs/default.yaml')
    validate_config(config)
    
    # Create output directory
    output_dir = Path('results/performance')
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Define test population sizes
    population_sizes = [50, 100, 200, 500, 1000]
    
    print(f"\n📊 Testing population sizes: {population_sizes}")
    
    # Benchmark population initialization
    init_results = benchmark_population_initialization(config, population_sizes[:4])  # Smaller set for init
    
    # Benchmark full simulations
    sim_results = benchmark_full_simulation_scaling(config, population_sizes[:4])  # Limit to avoid timeouts
    
    # Save raw results
    raw_results = {
        'initialization': init_results,
        'simulation': sim_results,
        'timestamp': time.time()
    }
    
    with open(output_dir / 'simple_benchmark_data.json', 'w') as f:
        json.dump(raw_results, f, indent=2)
    
    # Create reports and visualizations
    print("\n📝 Generating reports...")
    create_benchmark_report(init_results, sim_results, output_dir)
    create_benchmark_visualization(init_results, sim_results, output_dir)
    
    print(f"\n✅ Simple benchmark complete!")
    print(f"   📄 Report: {output_dir}/SIMPLE_BENCHMARK_REPORT.md")
    print(f"   📊 Visualization: {output_dir}/simple_benchmark.png")
    print(f"   🔧 Raw data: {output_dir}/simple_benchmark_data.json")
    
    # Print quick summary
    print("\n🎯 Quick Summary:")
    
    # Find max successful population
    max_init = 0
    max_sim = 0
    
    for n, t in zip(init_results['n_agents'], init_results['init_time']):
        if not np.isnan(t):
            max_init = max(max_init, n)
    
    for n, t in zip(sim_results['n_agents'], sim_results['time_2yr']):
        if not np.isnan(t):
            max_sim = max(max_sim, n)
    
    print(f"   Max successful initialization: {max_init:,} agents")
    print(f"   Max successful simulation: {max_sim:,} agents")
    
    # Scaling analysis
    sim_exp = estimate_scaling_exponent(sim_results['n_agents'], sim_results['time_2yr'])
    if not np.isnan(sim_exp):
        print(f"   Simulation scaling: O(N^{sim_exp:.2f})")
        
        if sim_exp < 1.5:
            print("   → Nearly linear scaling (excellent)")
        elif sim_exp < 2.0:
            print("   → Super-linear but reasonable")
        elif sim_exp < 2.5:
            print("   → Quadratic-like (concerning for large N)")
        else:
            print("   → Worse than quadratic (major bottleneck)")
    else:
        print("   Scaling analysis: Insufficient data")


if __name__ == '__main__':
    main()