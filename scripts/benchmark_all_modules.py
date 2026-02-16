#!/usr/bin/env python3
"""
Comprehensive benchmark of ALL SSWD-EvoEpi model components.

This script measures the performance of each module independently and at different scales
to identify bottlenecks and scaling behavior across the entire model.
"""

import time
import numpy as np
import pandas as pd
import json
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from concurrent.futures import ThreadPoolExecutor, TimeoutError
import signal
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
# import seaborn as sns  # Not needed for this benchmark

# Add project root to path
import sys
sys.path.append(str(Path(__file__).parent.parent))

from sswd_evoepi.config import SimulationConfig, load_config, validate_config
from sswd_evoepi.disease import daily_disease_update, NodeDiseaseState
from sswd_evoepi.movement import daily_movement, update_movement
from sswd_evoepi.spawning import spawning_step, in_spawning_season
from sswd_evoepi.reproduction import srs_reproductive_lottery, mendelian_inherit_batch
from sswd_evoepi.spatial import compute_distance_matrix, construct_larval_connectivity, construct_pathogen_dispersal
from sswd_evoepi.environment import sinusoidal_sst, seasonal_flushing, compute_environmental_state
from sswd_evoepi.model import initialize_population, annual_natural_mortality, annual_growth_and_aging, annual_reproduction, run_coupled_simulation
from sswd_evoepi.types import AGENT_DTYPE
from sswd_evoepi.rng import create_rng_hierarchy


class BenchmarkTimeout(Exception):
    """Raised when a benchmark operation times out"""
    pass


def get_benchmark_rng():
    """Get a fresh RNG for benchmarking"""
    return np.random.default_rng(42)


def timeout_handler(signum, frame):
    raise BenchmarkTimeout("Benchmark operation timed out")


def with_timeout(timeout_seconds: int = 60):
    """Decorator to add timeout to benchmark functions"""
    def decorator(func):
        def wrapper(*args, **kwargs):
            # Set up timeout handler
            signal.signal(signal.SIGALRM, timeout_handler)
            signal.alarm(timeout_seconds)
            
            try:
                result = func(*args, **kwargs)
            except BenchmarkTimeout:
                print(f"  WARNING: {func.__name__} timed out after {timeout_seconds}s")
                return None
            finally:
                signal.alarm(0)  # Cancel timeout
            
            return result
        return wrapper
    return decorator


def create_test_agents(n_agents: int, config: SimulationConfig) -> Tuple[np.ndarray, np.ndarray]:
    """Create test population for benchmarking"""
    from sswd_evoepi.model import make_effect_sizes
    
    rng_hierarchy = create_rng_hierarchy(master_seed=42, n_nodes=1)
    rng = rng_hierarchy['global']
    
    # Create effect sizes
    effect_sizes = make_effect_sizes(seed=42)
    
    # Initialize population
    habitat_area = 100_000_000  # 100 km² in m²
    agents, genotypes = initialize_population(
        n_individuals=n_agents,
        max_agents=n_agents + 100,  # Some buffer
        habitat_area=habitat_area,
        effect_sizes=effect_sizes,
        pop_cfg=config.population,
        rng=rng
    )
    
    return agents, genotypes


@with_timeout(60)
def benchmark_disease_module(agents: np.ndarray, config: SimulationConfig, n_iterations: int = 1000) -> Optional[float]:
    """Benchmark the disease module"""
    print("  Benchmarking disease module...")
    
    # Initialize disease state
    disease_state = NodeDiseaseState(node_id=0)
    
    # Set up environmental conditions
    sst = 15.0  # Celsius
    salinity = 30.0  # psu
    vibrio_conc = np.array([1000.0])  # cells/mL
    
    # Time the disease updates
    start_time = time.perf_counter()
    
    for _ in range(n_iterations):
        daily_disease_update(
            agents=agents,
            disease_state=disease_state,
            sst=np.array([sst]),
            salinity=np.array([salinity]),
            vibrio_conc=vibrio_conc,
            doy=100,
            cfg=config,
            rng_state=get_benchmark_rng()
        )
    
    end_time = time.perf_counter()
    return (end_time - start_time) / n_iterations


@with_timeout(60)
def benchmark_movement_module(agents: np.ndarray, config: SimulationConfig, n_iterations: int = 1000) -> Optional[float]:
    """Benchmark the movement module"""
    print("  Benchmarking movement module...")
    
    start_time = time.perf_counter()
    
    for _ in range(n_iterations):
        daily_movement(
            agents=agents,
            cfg=config.movement,
            habitat_side=10000.0,  # 10 km
            rng_state=get_benchmark_rng()
        )
    
    end_time = time.perf_counter()
    return (end_time - start_time) / n_iterations


@with_timeout(60)
def benchmark_spawning_module(agents: np.ndarray, genotypes: np.ndarray, config: SimulationConfig, n_iterations: int = 270) -> Optional[float]:
    """Benchmark the spawning module (one spawning season)"""
    print("  Benchmarking spawning module...")
    
    start_time = time.perf_counter()
    
    for doy in range(305, 575):  # Nov-Jul spawning season
        if in_spawning_season(doy % 365):
            spawning_step(
                agents=agents,
                genotypes=genotypes,
                doy=doy % 365,
                cfg=config,
                rng_state=get_benchmark_rng()
            )
    
    end_time = time.perf_counter()
    return (end_time - start_time) / n_iterations


@with_timeout(60)
def benchmark_reproduction_module(agents: np.ndarray, genotypes: np.ndarray, config: SimulationConfig, n_iterations: int = 100) -> Optional[float]:
    """Benchmark the reproduction module (SRS + Mendelian inheritance)"""
    print("  Benchmarking reproduction module...")
    
    # Get breeding adults
    breeding_mask = (agents['stage'] >= 4) & (agents['health'] == 0)  # Healthy adults
    breeding_agents = agents[breeding_mask]
    breeding_genotypes = genotypes[breeding_mask]
    
    if len(breeding_agents) < 10:
        print("    WARNING: Too few breeding adults for meaningful benchmark")
        return None
    
    start_time = time.perf_counter()
    
    for _ in range(n_iterations):
        # SRS lottery
        n_offspring = 1000
        parent_pairs, offspring_counts = srs_reproductive_lottery(
            breeding_agents,
            n_offspring,
            alpha=1.35,
            rng_state=get_benchmark_rng()
        )
        
        # Mendelian inheritance
        if len(parent_pairs) > 0:
            offspring_genotypes = mendelian_inherit_batch(
                parent_pairs,
                breeding_genotypes,
                offspring_counts,
                rng_state=get_benchmark_rng()
            )
    
    end_time = time.perf_counter()
    return (end_time - start_time) / n_iterations


@with_timeout(60)
def benchmark_spatial_module(config: SimulationConfig, n_iterations: int = 10) -> Optional[float]:
    """Benchmark spatial operations (distance matrices, connectivity)"""
    print("  Benchmarking spatial module...")
    
    # Create test coordinates
    n_sites = 50
    lats = np.random.uniform(30, 60, n_sites)
    lons = np.random.uniform(-140, -120, n_sites)
    
    start_time = time.perf_counter()
    
    for _ in range(n_iterations):
        # Distance matrix computation
        distances = compute_distance_matrix(lats, lons)
        
        # Connectivity matrices
        larval_conn = construct_larval_connectivity(distances, config.spatial.larval_dispersal_km)
        pathogen_conn = construct_pathogen_dispersal(distances, config.spatial.pathogen_dispersal_km)
    
    end_time = time.perf_counter()
    return (end_time - start_time) / n_iterations


@with_timeout(60)
def benchmark_environment_module(config: SimulationConfig, n_iterations: int = 10000) -> Optional[float]:
    """Benchmark environment calculations"""
    print("  Benchmarking environment module...")
    
    start_time = time.perf_counter()
    
    for _ in range(n_iterations):
        # Daily environment calculations
        sst = sinusoidal_sst(100, 12.0, 5.0, phase_shift=0.0)
        flushing = seasonal_flushing(0.1, 4, phase_offset_months=2)
        env_state = compute_environmental_state(sst, 30.0, 0.1, 0.05)
    
    end_time = time.perf_counter()
    return (end_time - start_time) / n_iterations


@with_timeout(60)
def benchmark_growth_mortality_module(agents: np.ndarray, config: SimulationConfig, n_iterations: int = 100) -> Optional[float]:
    """Benchmark annual growth and mortality"""
    print("  Benchmarking growth/mortality module...")
    
    start_time = time.perf_counter()
    
    for _ in range(n_iterations):
        # Annual natural mortality
        annual_natural_mortality(agents, config, rng_state=get_benchmark_rng())
        
        # Annual growth and aging
        annual_growth_and_aging(agents, config)
    
    end_time = time.perf_counter()
    return (end_time - start_time) / n_iterations


def benchmark_full_simulation_scaling(config: SimulationConfig, n_agents_list: List[int]) -> Dict[str, List[float]]:
    """Benchmark full simulations at different scales"""
    print("  Benchmarking full simulation scaling...")
    
    results = {
        'n_agents': [],
        'time_with_spawning': [],
        'time_without_spawning': []
    }
    
    for n_agents in n_agents_list:
        print(f"    Testing {n_agents} agents...")
        
        # With spawning
        config.spawning.enabled = True
        start_time = time.perf_counter()
        try:
            run_coupled_simulation(
                config=config,
                n_agents=n_agents,
                n_years=2,
                output_dir=None,
                verbose=False
            )
            time_with_spawning = time.perf_counter() - start_time
        except Exception as e:
            print(f"      ERROR in spawning simulation: {e}")
            time_with_spawning = float('nan')
        
        # Without spawning
        config.spawning.enabled = False
        start_time = time.perf_counter()
        try:
            run_coupled_simulation(
                config=config,
                n_agents=n_agents,
                n_years=2,
                output_dir=None,
                verbose=False
            )
            time_without_spawning = time.perf_counter() - start_time
        except Exception as e:
            print(f"      ERROR in non-spawning simulation: {e}")
            time_without_spawning = float('nan')
        
        results['n_agents'].append(n_agents)
        results['time_with_spawning'].append(time_with_spawning)
        results['time_without_spawning'].append(time_without_spawning)
    
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


def create_benchmark_report(results: Dict[str, float], scaling_results: Dict[str, List[float]], output_dir: Path):
    """Create comprehensive benchmark report"""
    
    report_path = output_dir / 'BENCHMARK_REPORT.md'
    
    with open(report_path, 'w') as f:
        f.write("# SSWD-EvoEpi Comprehensive Benchmark Report\n\n")
        f.write(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        # Module performance table
        f.write("## Module Performance Breakdown\n\n")
        f.write("| Component | Time per call (ms) | Calls per year | Annual cost (s) | Notes |\n")
        f.write("|-----------|-------------------|----------------|-----------------|-------|\n")
        
        # Calculate annual costs (assuming daily calls for most modules)
        daily_modules = ['disease', 'movement', 'environment']
        annual_modules = ['growth_mortality']
        seasonal_modules = ['spawning', 'reproduction']  # 270 days, 1 time per year respectively
        
        total_annual_cost = 0.0
        
        for module, time_sec in results.items():
            if time_sec is None:
                f.write(f"| {module.replace('_', ' ').title()} | TIMEOUT | - | - | Timed out |\n")
                continue
            
            time_ms = time_sec * 1000
            
            if module in daily_modules:
                calls_per_year = 365
            elif module in annual_modules:
                calls_per_year = 1
            elif module == 'spawning':
                calls_per_year = 270  # Spawning season length
            elif module == 'reproduction':
                calls_per_year = 1
            elif module == 'spatial':
                calls_per_year = 1  # One-time setup
            else:
                calls_per_year = 365  # Default to daily
            
            annual_cost = time_sec * calls_per_year
            total_annual_cost += annual_cost
            
            f.write(f"| {module.replace('_', ' ').title()} | {time_ms:.3f} | {calls_per_year} | {annual_cost:.3f} | |\n")
        
        f.write(f"\n**Total annual simulation cost:** {total_annual_cost:.2f} seconds\n\n")
        
        # Scaling analysis
        f.write("## Scaling Analysis\n\n")
        f.write("| Population Size | With Spawning (s) | Without Spawning (s) | Spawning Overhead |\n")
        f.write("|-----------------|-------------------|---------------------|-------------------|\n")
        
        for i, n in enumerate(scaling_results['n_agents']):
            t_spawn = scaling_results['time_with_spawning'][i]
            t_no_spawn = scaling_results['time_without_spawning'][i]
            
            if np.isnan(t_spawn) or np.isnan(t_no_spawn):
                overhead = "N/A"
            else:
                overhead = f"{((t_spawn - t_no_spawn) / t_no_spawn) * 100:.1f}%"
            
            f.write(f"| {n} | {t_spawn:.2f} | {t_no_spawn:.2f} | {overhead} |\n")
        
        # Scaling exponents
        spawn_exp = estimate_scaling_exponent(scaling_results['n_agents'], scaling_results['time_with_spawning'])
        no_spawn_exp = estimate_scaling_exponent(scaling_results['n_agents'], scaling_results['time_without_spawning'])
        
        f.write(f"\n**Scaling exponents:**\n")
        f.write(f"- With spawning: O(N^{spawn_exp:.2f})\n")
        f.write(f"- Without spawning: O(N^{no_spawn_exp:.2f})\n\n")
        
        # Top bottlenecks
        f.write("## Top Performance Bottlenecks\n\n")
        
        # Sort modules by annual cost
        module_costs = []
        for module, time_sec in results.items():
            if time_sec is None:
                continue
            
            if module in daily_modules:
                calls_per_year = 365
            elif module == 'spawning':
                calls_per_year = 270
            elif module in ['reproduction', 'growth_mortality', 'spatial']:
                calls_per_year = 1
            else:
                calls_per_year = 365
                
            annual_cost = time_sec * calls_per_year
            module_costs.append((module, annual_cost, time_sec))
        
        module_costs.sort(key=lambda x: x[1], reverse=True)
        
        for i, (module, annual_cost, per_call) in enumerate(module_costs[:5], 1):
            pct_total = (annual_cost / total_annual_cost) * 100
            f.write(f"{i}. **{module.replace('_', ' ').title()}**: {annual_cost:.2f}s/year ({pct_total:.1f}% of total)\n")
        
        f.write("\n## Memory Usage Estimates\n\n")
        f.write("Approximate memory usage per agent:\n")
        f.write("- Agent struct: ~200 bytes\n")
        f.write("- Genotype data: ~50 bytes (50 loci)\n")
        f.write("- **Total per agent: ~250 bytes**\n\n")
        
        f.write("| Population Size | Memory Usage (MB) |\n")
        f.write("|-----------------|-------------------|\n")
        for n in [100, 500, 1000, 5000, 10000, 50000]:
            memory_mb = (n * 250) / (1024 * 1024)
            f.write(f"| {n:,} | {memory_mb:.2f} |\n")


def create_benchmark_visualization(results: Dict[str, float], scaling_results: Dict[str, List[float]], output_dir: Path):
    """Create benchmark visualization"""
    
    plt.style.use('dark_background')
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    
    # 1. Module performance breakdown (stacked bar)
    modules = []
    times_ms = []
    colors = []
    
    color_map = {
        'disease': '#ff6b6b',
        'movement': '#4ecdc4', 
        'spawning': '#45b7d1',
        'reproduction': '#96ceb4',
        'spatial': '#feca57',
        'environment': '#ff9ff3',
        'growth_mortality': '#54a0ff'
    }
    
    for module, time_sec in results.items():
        if time_sec is not None:
            modules.append(module.replace('_', ' ').title())
            times_ms.append(time_sec * 1000)
            colors.append(color_map.get(module, '#95a5a6'))
    
    bars = ax1.bar(modules, times_ms, color=colors, alpha=0.8)
    ax1.set_ylabel('Time per Call (ms)', fontsize=12)
    ax1.set_title('Module Performance Breakdown', fontsize=14, fontweight='bold')
    ax1.tick_params(axis='x', rotation=45)
    
    # Add value labels on bars
    for bar, time_ms in zip(bars, times_ms):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height + height*0.02,
                f'{time_ms:.2f}ms', ha='center', va='bottom', fontsize=10)
    
    # 2. Scaling behavior
    n_agents = scaling_results['n_agents']
    time_spawn = scaling_results['time_with_spawning']
    time_no_spawn = scaling_results['time_without_spawning']
    
    # Filter out NaN values
    valid_spawn = [(n, t) for n, t in zip(n_agents, time_spawn) if not np.isnan(t)]
    valid_no_spawn = [(n, t) for n, t in zip(n_agents, time_no_spawn) if not np.isnan(t)]
    
    if valid_spawn:
        n_spawn, t_spawn = zip(*valid_spawn)
        ax2.loglog(n_spawn, t_spawn, 'o-', color='#ff6b6b', linewidth=2, markersize=8, label='With Spawning')
    
    if valid_no_spawn:
        n_no_spawn, t_no_spawn = zip(*valid_no_spawn)
        ax2.loglog(n_no_spawn, t_no_spawn, 's-', color='#4ecdc4', linewidth=2, markersize=8, label='Without Spawning')
    
    ax2.set_xlabel('Population Size', fontsize=12)
    ax2.set_ylabel('Runtime (seconds)', fontsize=12)
    ax2.set_title('Scaling Behavior (2-year simulation)', fontsize=14, fontweight='bold')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # 3. Spawning overhead
    overhead_data = []
    n_valid = []
    
    for i, n in enumerate(n_agents):
        t_spawn = time_spawn[i]
        t_no_spawn = time_no_spawn[i]
        
        if not np.isnan(t_spawn) and not np.isnan(t_no_spawn) and t_no_spawn > 0:
            overhead_pct = ((t_spawn - t_no_spawn) / t_no_spawn) * 100
            overhead_data.append(overhead_pct)
            n_valid.append(n)
    
    if overhead_data:
        bars = ax3.bar(range(len(n_valid)), overhead_data, color='#feca57', alpha=0.8)
        ax3.set_xlabel('Population Size', fontsize=12)
        ax3.set_ylabel('Spawning Overhead (%)', fontsize=12)
        ax3.set_title('Spawning Module Overhead', fontsize=14, fontweight='bold')
        ax3.set_xticks(range(len(n_valid)))
        ax3.set_xticklabels([str(n) for n in n_valid])
        
        # Add value labels
        for i, (bar, overhead) in enumerate(zip(bars, overhead_data)):
            height = bar.get_height()
            ax3.text(bar.get_x() + bar.get_width()/2., height + height*0.02,
                    f'{overhead:.1f}%', ha='center', va='bottom', fontsize=10)
    else:
        ax3.text(0.5, 0.5, 'No valid overhead data', ha='center', va='center', transform=ax3.transAxes)
    
    # 4. Memory usage projection
    pop_sizes = np.array([100, 500, 1000, 5000, 10000, 50000, 100000])
    memory_mb = (pop_sizes * 250) / (1024 * 1024)  # 250 bytes per agent
    
    ax4.loglog(pop_sizes, memory_mb, 'o-', color='#96ceb4', linewidth=3, markersize=10)
    ax4.set_xlabel('Population Size', fontsize=12)
    ax4.set_ylabel('Memory Usage (MB)', fontsize=12)
    ax4.set_title('Memory Usage Projection', fontsize=14, fontweight='bold')
    ax4.grid(True, alpha=0.3)
    
    # Add reference lines
    ax4.axhline(y=1000, color='red', linestyle='--', alpha=0.7, label='1 GB')
    ax4.axhline(y=4000, color='orange', linestyle='--', alpha=0.7, label='4 GB')
    ax4.legend()
    
    plt.tight_layout()
    plt.savefig(output_dir / 'benchmark_breakdown.png', dpi=300, bbox_inches='tight', 
                facecolor='black', edgecolor='none')
    plt.close()


def main():
    """Main benchmark execution"""
    print("🔬 Starting comprehensive SSWD-EvoEpi benchmark...")
    
    # Load configuration
    config = load_config('configs/default.yaml')
    validate_config(config)
    
    # Create output directory
    output_dir = Path('results/performance')
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Test with 500 agents for module benchmarks
    print("\n📊 Creating test population (500 agents)...")
    agents, genotypes = create_test_agents(500, config)
    
    # Module benchmarks
    print("\n🧪 Running module benchmarks...")
    
    results = {}
    
    # Disease module
    results['disease'] = benchmark_disease_module(agents, config)
    
    # Movement module
    results['movement'] = benchmark_movement_module(agents, config)
    
    # Spawning module  
    results['spawning'] = benchmark_spawning_module(agents, genotypes, config)
    
    # Reproduction module
    results['reproduction'] = benchmark_reproduction_module(agents, genotypes, config)
    
    # Spatial module
    results['spatial'] = benchmark_spatial_module(config)
    
    # Environment module
    results['environment'] = benchmark_environment_module(config)
    
    # Growth/mortality module
    results['growth_mortality'] = benchmark_growth_mortality_module(agents, config)
    
    # Full simulation scaling
    print("\n📈 Running scaling benchmarks...")
    n_agents_list = [50, 100, 200, 500]
    scaling_results = benchmark_full_simulation_scaling(config, n_agents_list)
    
    # Save raw results
    raw_results = {
        'module_times': results,
        'scaling': scaling_results,
        'timestamp': time.time()
    }
    
    with open(output_dir / 'benchmark_raw_data.json', 'w') as f:
        json.dump(raw_results, f, indent=2)
    
    # Create reports and visualizations
    print("\n📝 Generating reports...")
    create_benchmark_report(results, scaling_results, output_dir)
    create_benchmark_visualization(results, scaling_results, output_dir)
    
    print(f"\n✅ Benchmark complete! Results saved to {output_dir}/")
    print(f"   📄 Report: {output_dir}/BENCHMARK_REPORT.md")
    print(f"   📊 Visualization: {output_dir}/benchmark_breakdown.png")
    print(f"   🔧 Raw data: {output_dir}/benchmark_raw_data.json")
    
    # Print quick summary
    print("\n🎯 Quick Summary:")
    print("   Top 3 most expensive modules (per year):")
    
    # Calculate annual costs for ranking
    module_costs = []
    for module, time_sec in results.items():
        if time_sec is None:
            continue
        
        if module in ['disease', 'movement', 'environment']:
            calls_per_year = 365
        elif module == 'spawning':
            calls_per_year = 270
        else:
            calls_per_year = 1
        
        annual_cost = time_sec * calls_per_year
        module_costs.append((module, annual_cost))
    
    module_costs.sort(key=lambda x: x[1], reverse=True)
    
    for i, (module, cost) in enumerate(module_costs[:3], 1):
        print(f"   {i}. {module.replace('_', ' ').title()}: {cost:.2f}s/year")


if __name__ == '__main__':
    main()