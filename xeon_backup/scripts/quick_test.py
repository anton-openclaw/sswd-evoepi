#!/usr/bin/env python3
"""Quick test to verify all systems are working."""

import time
from sswd_evoepi.config import default_config
from sswd_evoepi.model import run_spatial_simulation, run_coupled_simulation
from sswd_evoepi.spatial import make_5node_network

def test_single_node():
    """Test single-node simulation (fastest)."""
    print("=== Single-Node Test ===")
    start_time = time.time()
    
    result = run_coupled_simulation(
        n_individuals=100,  # Smaller population for speed
        carrying_capacity=100,
        habitat_area=5000,
        T_celsius=14.0,
        salinity=30.0,
        phi_k=0.02,
        n_years=5,  # Shorter simulation
        disease_year=2,
        initial_infected=3,
        seed=42
    )
    
    runtime = time.time() - start_time
    print(f"Single-node test completed in {runtime:.1f} seconds")
    print(f"  Initial: {result.initial_pop} → Final: {result.final_pop}")
    print(f"  Disease deaths: {result.total_disease_deaths}")
    print(f"  Final resistance: {result.yearly_mean_resistance[-1]:.4f}")
    
    return result

def test_spatial_short():
    """Test 5-node spatial simulation (short duration)."""
    print("\n=== 5-Node Spatial Test (5 years) ===")
    start_time = time.time()
    
    network = make_5node_network()
    config = default_config()
    
    results = run_spatial_simulation(
        network=network,
        n_years=5,  # Much shorter
        disease_year=2,
        initial_infected_per_node=3,
        seed=42,
        config=config
    )
    
    runtime = time.time() - start_time
    print(f"5-node spatial test completed in {runtime:.1f} seconds")
    
    # Show results for each node
    for i, node in enumerate(network.nodes):
        initial_pop = results.yearly_population[0][i]
        final_pop = results.yearly_population[-1][i]
        total_deaths = sum(results.yearly_disease_deaths[year][i] for year in range(5))
        print(f"  {node.definition.name}: {initial_pop} → {final_pop}, {total_deaths} deaths")
    
    return results

if __name__ == '__main__':
    # Test single node first (fastest)
    single_result = test_single_node()
    
    # Test spatial if single node works
    spatial_result = test_spatial_short()
    
    print("\n=== All Tests Passed ===")
    print("System is ready for full 20-year simulation.")