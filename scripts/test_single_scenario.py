#!/usr/bin/env python3
"""Test single scenario to debug comparison script."""

import sys
import time
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from sswd_evoepi.config import default_config
from sswd_evoepi.model import run_spatial_simulation
from sswd_evoepi.spatial import make_5node_network

def test_scenario():
    print("Testing single scenario...")
    
    # Create network
    print("Creating network (Haversine)...")
    network = make_5node_network(seed=42, use_overwater=False)
    print(f"Network created: {network.n_nodes} nodes")
    
    print("Starting simulation...")
    t0 = time.time()
    
    def progress(year, total):
        elapsed = time.time() - t0
        print(f"  Year {year}/{total} ({elapsed:.1f}s)")
    
    config = default_config()
    
    result = run_spatial_simulation(
        network=network,
        n_years=5,  # Short test
        disease_year=2,
        initial_infected_per_node=5,
        seed=42,
        config=config,
        progress_callback=progress,
    )
    
    elapsed = time.time() - t0
    print(f"Simulation complete in {elapsed:.1f}s")
    print(f"Final total population: {result.yearly_total_pop[-1]}")
    
    return result

if __name__ == "__main__":
    test_scenario()