#!/usr/bin/env python3
"""Simple profiler for spawning system."""

import cProfile
import sys
sys.path.insert(0, '.')

from sswd_evoepi.model import run_coupled_simulation
from sswd_evoepi.config import default_config

print("Running simple profile...")

cfg = default_config()
cfg.spawning.enabled = True

def simple_run():
    return run_coupled_simulation(
        n_individuals=50, 
        carrying_capacity=50, 
        habitat_area=2000, 
        T_celsius=14, 
        salinity=30, 
        phi_k=0.02, 
        n_years=2, 
        disease_year=999, 
        seed=42, 
        config=cfg
    )

# Profile it
print("Profiling 50-agent 2-year simulation...")
cProfile.run('simple_run()', 'simple_profile.prof')
print("Profiling complete.")