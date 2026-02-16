#!/usr/bin/env python3
"""Simple timing test to verify basic functionality."""

import time
import sys
import copy
from pathlib import Path

# Add sswd_evoepi to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from sswd_evoepi.model import run_coupled_simulation
from sswd_evoepi.config import SimulationConfig

def main():
    print("Simple timing test...")
    
    # Create configs
    base_config = SimulationConfig()
    spawning_config = copy.deepcopy(base_config)
    spawning_config.spawning.enabled = True
    
    print("Testing 10 agents, 1 year with spawning...")
    start = time.perf_counter()
    result = run_coupled_simulation(
        n_individuals=10,
        n_years=1,
        seed=42,
        config=spawning_config,
        record_daily=False
    )
    elapsed = time.perf_counter() - start
    print(f"  Completed in {elapsed:.2f}s")
    print(f"  Final population: {result.final_pop}")
    
    print("Testing 10 agents, 1 year without spawning...")
    no_spawning_config = copy.deepcopy(base_config)
    no_spawning_config.spawning.enabled = False
    
    start = time.perf_counter()
    result2 = run_coupled_simulation(
        n_individuals=10,
        n_years=1,
        seed=42,
        config=no_spawning_config,
        record_daily=False
    )
    elapsed2 = time.perf_counter() - start
    print(f"  Completed in {elapsed2:.2f}s")
    print(f"  Final population: {result2.final_pop}")
    
    print("Success! Basic simulation works.")

if __name__ == '__main__':
    main()