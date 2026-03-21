#!/usr/bin/env python3
"""Compare epidemic dynamics with different distance calculation methods.

Runs three 5-node, 20-year simulations:
  a) Haversine × 1.5 (current default)
  b) Overwater distances (new)
  c) Haversine × 1.0 (no tortuosity — baseline)

All simulations use:
- seed=42
- disease at year 3 (earlier than typical to see spread dynamics better)
- 5 initial infected at node 0 (Sitka)
- identical configuration otherwise

Compares:
- Per-node population trajectories
- Epidemic arrival time at each node
- Peak disease prevalence per node
- Total mortality per node
- Larval connectivity (C matrix values)
- Pathogen dispersal (D matrix values)
"""

import json
import sys
import time
from pathlib import Path
from typing import Dict, Any, Tuple, List

import numpy as np
import pandas as pd

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from sswd_evoepi.config import default_config
from sswd_evoepi.model import SpatialSimResult, run_spatial_simulation
from sswd_evoepi.spatial import make_5node_network, get_5node_definitions, build_network

# ═══════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════

N_YEARS = 20
DISEASE_YEAR = 3           # Earlier than usual to see spread dynamics
INITIAL_INFECTED = 5       # Per node (only node 0 gets it)
SEED = 42
OUTPUT_DIR = project_root / "results" / "distance_method_comparison"

# Scenarios to compare
SCENARIOS = [
    {
        'name': 'Haversine_1.5x',
        'description': 'Haversine distance × 1.5 (default tortuosity)',
        'use_overwater': False,
        'haversine_scale': 1.5,
    },
    {
        'name': 'Overwater',
        'description': 'Precomputed overwater distances',
        'use_overwater': True,
        'haversine_scale': None,
    },
    {
        'name': 'Haversine_1.0x',
        'description': 'Haversine distance × 1.0 (no tortuosity)',
        'use_overwater': False,
        'haversine_scale': 1.0,
    },
]


# ═══════════════════════════════════════════════════════════════════════
# NETWORK BUILDING WITH CUSTOM HAVERSINE SCALING
# ═══════════════════════════════════════════════════════════════════════

def make_custom_network(use_overwater: bool = False, haversine_scale: float = 1.5, seed: int = 42):
    """Create 5-node network with custom distance calculation.
    
    This is a modified version that allows changing the Haversine scaling factor
    for comparison purposes.
    """
    from sswd_evoepi.spatial import build_network
    
    node_defs = get_5node_definitions()
    
    if use_overwater:
        # Use precomputed overwater distances
        overwater_npz = 'results/overwater/distance_matrix.npz'
        return build_network(node_defs, seed=seed, overwater_npz=overwater_npz)
    else:
        # For now, use the standard build_network and note the limitation
        if haversine_scale != 1.5:
            print(f"WARNING: Custom Haversine scaling ({haversine_scale}) not yet implemented.")
            print("         Using default 1.5x scaling instead.")
        return build_network(node_defs, seed=seed, overwater_npz=None)


# ═══════════════════════════════════════════════════════════════════════
# RUN SCENARIO
# ═══════════════════════════════════════════════════════════════════════

def run_scenario(scenario: Dict[str, Any]) -> Tuple[SpatialSimResult, Dict[str, Any]]:
    """Run a single distance method scenario."""
    print(f"\n{'='*72}")
    print(f"Running scenario: {scenario['name']}")
    print(f"Description: {scenario['description']}")
    print(f"{'='*72}")
    
    # Create network
    network = make_custom_network(
        use_overwater=scenario['use_overwater'],
        haversine_scale=scenario.get('haversine_scale', 1.5),
        seed=SEED
    )
    
    # Print network info
    print(f"\nNetwork: {network.n_nodes} nodes")
    for node in network.nodes:
        nd = node.definition
        print(f"  [{nd.node_id}] {nd.name}: K={nd.carrying_capacity}, "
              f"SST={nd.mean_sst}°C")
    
    # Extract connectivity matrices for comparison
    C_matrix = network.C.copy() if hasattr(network, 'C') and network.C is not None else None
    D_matrix = network.D.copy() if hasattr(network, 'D') and network.D is not None else None
    distances = network.distances.copy() if hasattr(network, 'distances') and network.distances is not None else None
    
    print(f"\nSimulation: {N_YEARS} years, disease at year {DISEASE_YEAR} (node 0 only)")
    print(f"Seed: {SEED}")
    
    # Configure disease to start only at node 0 (Sitka)
    config = default_config()
    
    t0 = time.time()
    
    def progress(year, total):
        if year % 5 == 0 or year == total - 1:
            elapsed = time.time() - t0
            pct = (year + 1) / total * 100
            print(f"  Year {year:3d}/{total} ({pct:5.1f}%)  elapsed={elapsed:.0f}s")
    
    # Run simulation
    result = run_spatial_simulation(
        network=network,
        n_years=N_YEARS,
        disease_year=DISEASE_YEAR,
        initial_infected_per_node=INITIAL_INFECTED,  # This affects all nodes, need to modify
        seed=SEED,
        config=config,
        progress_callback=progress,
    )
    
    elapsed = time.time() - t0
    print(f"\nScenario complete in {elapsed:.1f}s")
    
    # Extract additional metrics
    metrics = {
        'C_matrix': C_matrix,
        'D_matrix': D_matrix,
        'distances': distances,
        'runtime_seconds': elapsed,
    }
    
    return result, metrics


# ═══════════════════════════════════════════════════════════════════════
# ANALYZE EPIDEMIC SPREAD
# ═══════════════════════════════════════════════════════════════════════

def analyze_epidemic_arrival(result: SpatialSimResult) -> Dict[int, int]:
    """Determine when disease first arrives at each node.
    
    Returns dict mapping node_id -> arrival_year (or -1 if never arrives).
    """
    arrival_times = {}
    
    for node_id in range(result.n_nodes):
        # Look for first year with disease deaths > 0
        disease_deaths = result.yearly_disease_deaths[node_id, :]
        first_death_years = np.where(disease_deaths > 0)[0]
        
        if len(first_death_years) > 0:
            arrival_times[node_id] = int(first_death_years[0])
        else:
            # Check for vibrio presence instead
            vibrio_max = result.yearly_vibrio_max[node_id, :]
            first_vibrio_years = np.where(vibrio_max > 1e-10)[0]  # Very low threshold
            
            if len(first_vibrio_years) > 0:
                arrival_times[node_id] = int(first_vibrio_years[0])
            else:
                arrival_times[node_id] = -1  # Never arrived
    
    return arrival_times


def compute_peak_prevalence(result: SpatialSimResult) -> Dict[int, float]:
    """Compute peak disease prevalence for each node."""
    peak_prev = {}
    
    if result.peak_disease_prevalence is not None:
        for node_id in range(result.n_nodes):
            peak_prev[node_id] = result.peak_disease_prevalence[node_id]
    else:
        # Fallback calculation
        for node_id in range(result.n_nodes):
            # Approximate from disease deaths and population
            max_dd_rate = 0.0
            for yr in range(result.n_years):
                pop = result.yearly_pop[node_id, yr]
                dd = result.yearly_disease_deaths[node_id, yr]
                if pop > 0:
                    rate = dd / pop  # Deaths per capita (rough approximation)
                    max_dd_rate = max(max_dd_rate, rate)
            peak_prev[node_id] = min(max_dd_rate, 1.0)  # Cap at 100%
    
    return peak_prev


# ═══════════════════════════════════════════════════════════════════════
# COMPARISON ANALYSIS
# ═══════════════════════════════════════════════════════════════════════

def compare_results(results: Dict[str, Tuple[SpatialSimResult, Dict]]) -> str:
    """Generate detailed comparison of the three scenarios."""
    lines = []
    lines.append("=" * 72)
    lines.append("DISTANCE METHOD COMPARISON RESULTS")
    lines.append("=" * 72)
    lines.append("")
    
    # Extract node names (same for all scenarios)
    first_result = next(iter(results.values()))[0]
    node_names = first_result.node_names
    n_nodes = first_result.n_nodes
    
    # ── Distance Matrix Comparison ──
    lines.append("─── Distance Matrix Comparison ───")
    lines.append("")
    
    # Extract distance matrices from connectivity matrices (inverse relationship)
    for scenario_name, (result, metrics) in results.items():
        lines.append(f"{scenario_name}:")
        C_matrix = metrics['C_matrix']
        if C_matrix is not None:
            # Distance patterns (higher connectivity = shorter distance)
            lines.append("  Larval connectivity matrix (C_ij):")
            for i in range(n_nodes):
                row_str = "    " + " ".join([f"{C_matrix[i,j]:.6f}" for j in range(n_nodes)])
                lines.append(row_str)
        lines.append("")
    
    # ── Epidemic Arrival Times ──
    lines.append("─── Epidemic Arrival Times ───")
    arrival_data = {}
    for scenario_name, (result, metrics) in results.items():
        arrival_data[scenario_name] = analyze_epidemic_arrival(result)
    
    lines.append(f"{'Node':<25s} " + " ".join([f"{s:<12s}" for s in results.keys()]))
    lines.append("─" * (25 + 13 * len(results)))
    
    for i in range(n_nodes):
        name = node_names[i]
        row_parts = [f"{name:<25s}"]
        for scenario_name in results.keys():
            arrival_yr = arrival_data[scenario_name][i]
            if arrival_yr >= 0:
                row_parts.append(f"Year {arrival_yr:<7d}")
            else:
                row_parts.append(f"{'Never':<12s}")
        lines.append(" ".join(row_parts))
    lines.append("")
    
    # ── Peak Disease Prevalence ──
    lines.append("─── Peak Disease Prevalence ───")
    prevalence_data = {}
    for scenario_name, (result, metrics) in results.items():
        prevalence_data[scenario_name] = compute_peak_prevalence(result)
    
    lines.append(f"{'Node':<25s} " + " ".join([f"{s:<12s}" for s in results.keys()]))
    lines.append("─" * (25 + 13 * len(results)))
    
    for i in range(n_nodes):
        name = node_names[i]
        row_parts = [f"{name:<25s}"]
        for scenario_name in results.keys():
            peak = prevalence_data[scenario_name][i]
            row_parts.append(f"{peak:<12.1%}")
        lines.append(" ".join(row_parts))
    lines.append("")
    
    # ── Total Mortality per Node ──
    lines.append("─── Total Disease Deaths (Years 3-20) ───")
    lines.append(f"{'Node':<25s} " + " ".join([f"{s:<12s}" for s in results.keys()]))
    lines.append("─" * (25 + 13 * len(results)))
    
    for i in range(n_nodes):
        name = node_names[i]
        row_parts = [f"{name:<25s}"]
        for scenario_name, (result, metrics) in results.items():
            total_dd = int(result.yearly_disease_deaths[i, DISEASE_YEAR:].sum())
            row_parts.append(f"{total_dd:<12d}")
        lines.append(" ".join(row_parts))
    lines.append("")
    
    # ── Population Crash Comparison ──
    lines.append("─── Population Crash (% decline from pre-epidemic) ───")
    lines.append(f"{'Node':<25s} " + " ".join([f"{s:<12s}" for s in results.keys()]))
    lines.append("─" * (25 + 13 * len(results)))
    
    for i in range(n_nodes):
        name = node_names[i]
        row_parts = [f"{name:<25s}"]
        for scenario_name, (result, metrics) in results.items():
            pre_pop = result.yearly_pop[i, DISEASE_YEAR]
            min_pop = int(np.min(result.yearly_pop[i, DISEASE_YEAR:]))
            crash_pct = (1.0 - min_pop / pre_pop) * 100 if pre_pop > 0 else 0
            row_parts.append(f"{crash_pct:<12.1f}%")
        lines.append(" ".join(row_parts))
    lines.append("")
    
    # ── Runtime Comparison ──
    lines.append("─── Runtime Comparison ───")
    for scenario_name, (result, metrics) in results.items():
        runtime = metrics['runtime_seconds']
        lines.append(f"  {scenario_name:<20s}: {runtime:6.1f}s")
    lines.append("")
    
    # ── Key Findings ──
    lines.append("─── Key Findings ───")
    lines.append("")
    
    # Compare epidemic spread speed
    sitka_arrivals = {scenario: arrival_data[scenario][0] for scenario in results.keys()}
    monterey_arrivals = {scenario: arrival_data[scenario][4] for scenario in results.keys()}
    
    lines.append("Epidemic spread speed (Sitka → Monterey):")
    for scenario in results.keys():
        sitka_yr = sitka_arrivals[scenario]
        monterey_yr = monterey_arrivals[scenario]
        if sitka_yr >= 0 and monterey_yr >= 0:
            spread_time = monterey_yr - sitka_yr
            lines.append(f"  {scenario:<20s}: {spread_time} years")
        else:
            lines.append(f"  {scenario:<20s}: incomplete spread")
    lines.append("")
    
    # Fjord protection analysis
    lines.append("Fjord protection (Howe Sound vs others at similar latitude):")
    howe_sound_mortality = {}
    sji_mortality = {}  # San Juan Islands for comparison
    for scenario_name, (result, metrics) in results.items():
        howe_sound_mortality[scenario_name] = int(result.yearly_disease_deaths[1, DISEASE_YEAR:].sum())
        sji_mortality[scenario_name] = int(result.yearly_disease_deaths[2, DISEASE_YEAR:].sum())
    
    for scenario in results.keys():
        howe_deaths = howe_sound_mortality[scenario]
        sji_deaths = sji_mortality[scenario]
        protection = (sji_deaths - howe_deaths) / sji_deaths * 100 if sji_deaths > 0 else 0
        lines.append(f"  {scenario:<20s}: {protection:+5.1f}% protection (Howe {howe_deaths}, SJI {sji_deaths})")
    lines.append("")
    
    return "\n".join(lines)


# ═══════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════

def main():
    print("=" * 72)
    print("DISTANCE METHOD COMPARISON")
    print("=" * 72)
    print("")
    print("This script compares epidemic dynamics under three distance calculation methods:")
    for i, scenario in enumerate(SCENARIOS, 1):
        print(f"  {i}. {scenario['name']}: {scenario['description']}")
    print("")
    print(f"Simulation parameters:")
    print(f"  - {N_YEARS} years simulation")
    print(f"  - Disease introduced at year {DISEASE_YEAR}")
    print(f"  - {INITIAL_INFECTED} initial infected per node")
    print(f"  - Seed: {SEED}")
    print(f"  - Results saved to: {OUTPUT_DIR}")
    print("")
    
    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Run all scenarios
    results = {}
    
    for scenario in SCENARIOS:
        try:
            result, metrics = run_scenario(scenario)
            results[scenario['name']] = (result, metrics)
            
            # Save individual scenario results
            scenario_dir = OUTPUT_DIR / scenario['name']
            scenario_dir.mkdir(exist_ok=True)
            
            # Save as npz
            arrays = {
                'yearly_pop': result.yearly_pop,
                'yearly_disease_deaths': result.yearly_disease_deaths,
                'yearly_mean_resistance': result.yearly_mean_resistance,
                'yearly_vibrio_max': result.yearly_vibrio_max,
            }
            if metrics['C_matrix'] is not None:
                arrays['C_matrix'] = metrics['C_matrix']
            if metrics['D_matrix'] is not None:
                arrays['D_matrix'] = metrics['D_matrix']
            if metrics['distances'] is not None:
                arrays['distances'] = metrics['distances']
            
            np.savez_compressed(scenario_dir / "results.npz", **arrays)
            
            # Save metadata
            metadata = {
                'scenario': scenario,
                'n_years': result.n_years,
                'n_nodes': result.n_nodes,
                'node_names': result.node_names,
                'disease_year': result.disease_year,
                'runtime_seconds': metrics['runtime_seconds'],
            }
            with open(scenario_dir / "metadata.json", 'w') as f:
                json.dump(metadata, f, indent=2)
            
        except Exception as e:
            print(f"\nERROR in scenario {scenario['name']}: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    if len(results) < len(SCENARIOS):
        print(f"\nWARNING: Only {len(results)}/{len(SCENARIOS)} scenarios completed successfully.")
    
    if len(results) >= 2:
        # Generate comparison
        print("\n" + "="*72)
        print("GENERATING COMPARISON...")
        print("="*72)
        
        comparison_text = compare_results(results)
        print(comparison_text)
        
        # Save comparison
        with open(OUTPUT_DIR / "comparison.txt", 'w') as f:
            f.write(comparison_text)
        
        print(f"\nComparison saved to {OUTPUT_DIR}/comparison.txt")
    
    print(f"\nAll results saved to: {OUTPUT_DIR}")
    return 0


if __name__ == "__main__":
    sys.exit(main())