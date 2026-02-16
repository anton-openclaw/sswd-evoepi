#!/usr/bin/env python3
"""Quick epidemic comparison between Haversine and Overwater distances."""

import sys
import time
from pathlib import Path
import numpy as np

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from sswd_evoepi.config import default_config
from sswd_evoepi.model import run_spatial_simulation
from sswd_evoepi.spatial import make_5node_network, get_5node_definitions

def run_quick_comparison():
    print("=" * 60)
    print("QUICK EPIDEMIC COMPARISON")
    print("=" * 60)
    
    # Parameters
    n_years = 10
    disease_year = 3
    seed = 42
    
    print(f"Simulation: {n_years} years, disease at year {disease_year}")
    print(f"Seed: {seed}")
    print()
    
    scenarios = [
        ("Haversine×1.5", False),
        ("Overwater", True),
    ]
    
    results = {}
    
    for scenario_name, use_overwater in scenarios:
        print(f"Running {scenario_name}...")
        t0 = time.time()
        
        # Create network
        network = make_5node_network(seed=seed, use_overwater=use_overwater)
        
        # Progress callback
        progress_calls = []
        def progress(year, total):
            progress_calls.append((year, time.time() - t0))
            if year % 2 == 0:
                elapsed = time.time() - t0
                print(f"  Year {year}/{total} ({elapsed:.1f}s)")
        
        # Run simulation
        config = default_config()
        
        try:
            result = run_spatial_simulation(
                network=network,
                n_years=n_years,
                disease_year=disease_year,
                initial_infected_per_node=5,
                seed=seed,
                config=config,
                progress_callback=progress,
            )
            
            elapsed = time.time() - t0
            print(f"  Complete in {elapsed:.1f}s")
            print()
            
            results[scenario_name] = {
                'result': result,
                'runtime': elapsed,
                'network': network,
            }
            
        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    if len(results) >= 2:
        print_comparison(results, disease_year)
    
    return results

def print_comparison(results, disease_year):
    """Print epidemic comparison."""
    print("=" * 60)
    print("EPIDEMIC COMPARISON RESULTS")
    print("=" * 60)
    
    node_names = list(results.values())[0]['result'].node_names
    n_nodes = len(node_names)
    
    # Runtime comparison
    print("Runtime:")
    for scenario, data in results.items():
        print(f"  {scenario:<15s}: {data['runtime']:5.1f}s")
    print()
    
    # Final population comparison
    print("Final Population (Year 10):")
    print(f"{'Node':<25s} " + " ".join([f"{s:<12s}" for s in results.keys()]))
    print("─" * (25 + 13 * len(results)))
    
    for i, name in enumerate(node_names):
        row = [f"{name:<25s}"]
        for scenario, data in results.items():
            final_pop = int(data['result'].yearly_pop[i, -1])
            row.append(f"{final_pop:<12d}")
        print(" ".join(row))
    print()
    
    # Disease deaths comparison
    print(f"Total Disease Deaths (Years {disease_year}-10):")
    print(f"{'Node':<25s} " + " ".join([f"{s:<12s}" for s in results.keys()]))
    print("─" * (25 + 13 * len(results)))
    
    for i, name in enumerate(node_names):
        row = [f"{name:<25s}"]
        for scenario, data in results.items():
            total_dd = int(data['result'].yearly_disease_deaths[i, disease_year:].sum())
            row.append(f"{total_dd:<12d}")
        print(" ".join(row))
    print()
    
    # Population crash comparison
    print("Population Crash (% decline from pre-epidemic):")
    print(f"{'Node':<25s} " + " ".join([f"{s:<12s}" for s in results.keys()]))
    print("─" * (25 + 13 * len(results)))
    
    for i, name in enumerate(node_names):
        row = [f"{name:<25s}"]
        for scenario, data in results.items():
            pre_pop = data['result'].yearly_pop[i, disease_year]
            min_pop = np.min(data['result'].yearly_pop[i, disease_year:])
            crash_pct = (1 - min_pop / pre_pop) * 100 if pre_pop > 0 else 0
            row.append(f"{crash_pct:<12.1f}%")
        print(" ".join(row))
    print()
    
    # Epidemic arrival times (first disease death)
    print("Epidemic Arrival Time (first disease death):")
    print(f"{'Node':<25s} " + " ".join([f"{s:<12s}" for s in results.keys()]))
    print("─" * (25 + 13 * len(results)))
    
    for i, name in enumerate(node_names):
        row = [f"{name:<25s}"]
        for scenario, data in results.items():
            deaths = data['result'].yearly_disease_deaths[i, :]
            first_death_years = np.where(deaths > 0)[0]
            if len(first_death_years) > 0:
                arrival_year = first_death_years[0]
                row.append(f"Year {arrival_year:<7d}")
            else:
                row.append(f"{'Never':<12s}")
        print(" ".join(row))
    print()
    
    # Total metapopulation trajectory
    print("Total Metapopulation by Year:")
    print(f"{'Year':<6s} " + " ".join([f"{s:<12s}" for s in results.keys()]) + " Difference")
    print("─" * (6 + 13 * len(results) + 12))
    
    scenarios = list(results.keys())
    for year in range(10):
        row = [f"{year:<6d}"]
        pops = []
        for scenario in scenarios:
            pop = int(results[scenario]['result'].yearly_total_pop[year])
            pops.append(pop)
            row.append(f"{pop:<12d}")
        
        if len(pops) == 2:
            diff = pops[1] - pops[0]  # Overwater - Haversine
            row.append(f"{diff:+d}")
        
        row_str = " ".join(row)
        if year == disease_year:
            row_str += " ← disease"
        print(row_str)
    print()
    
    # Key findings
    print("KEY FINDINGS:")
    print()
    
    # Total mortality difference
    if len(results) == 2:
        scenarios = list(results.keys())
        haversine_data = results[scenarios[0]]
        overwater_data = results[scenarios[1]]
        
        h_total_deaths = int(haversine_data['result'].yearly_disease_deaths[:, disease_year:].sum())
        o_total_deaths = int(overwater_data['result'].yearly_disease_deaths[:, disease_year:].sum())
        
        print(f"Total disease deaths:")
        print(f"  {scenarios[0]:<15s}: {h_total_deaths:,}")
        print(f"  {scenarios[1]:<15s}: {o_total_deaths:,}")
        print(f"  Difference: {o_total_deaths - h_total_deaths:+,} ({(o_total_deaths/h_total_deaths - 1)*100:+.1f}%)")
        print()
        
        # Final population difference
        h_final = int(haversine_data['result'].yearly_total_pop[-1])
        o_final = int(overwater_data['result'].yearly_total_pop[-1])
        
        print(f"Final total population:")
        print(f"  {scenarios[0]:<15s}: {h_final:,}")
        print(f"  {scenarios[1]:<15s}: {o_final:,}")
        print(f"  Difference: {o_final - h_final:+,} ({(o_final/h_final - 1)*100:+.1f}%)")
        print()
        
        # Epidemic speed (Sitka to Monterey)
        h_deaths = haversine_data['result'].yearly_disease_deaths
        o_deaths = overwater_data['result'].yearly_disease_deaths
        
        # Monterey is node 4
        h_monterey_arrival = np.where(h_deaths[4, :] > 0)[0]
        o_monterey_arrival = np.where(o_deaths[4, :] > 0)[0]
        
        if len(h_monterey_arrival) > 0 and len(o_monterey_arrival) > 0:
            h_arrival = h_monterey_arrival[0]
            o_arrival = o_monterey_arrival[0]
            speed_diff = h_arrival - o_arrival
            
            print(f"Epidemic arrival at Monterey:")
            print(f"  {scenarios[0]:<15s}: Year {h_arrival}")
            print(f"  {scenarios[1]:<15s}: Year {o_arrival}")
            print(f"  Speed difference: {speed_diff:+d} years ({'faster' if speed_diff > 0 else 'slower'} with overwater)")


if __name__ == "__main__":
    run_quick_comparison()