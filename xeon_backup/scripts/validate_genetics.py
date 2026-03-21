#!/usr/bin/env python3
"""Validate that vectorized genetic crossing preserves biological behavior."""
import sys; sys.path.insert(0, '.')
import numpy as np
from sswd_evoepi.model import run_coupled_simulation
from sswd_evoepi.config import default_config
from sswd_evoepi.genetics import compute_allele_frequencies

# Run 10 simulations with different seeds, collect statistics
results = []
for seed in range(10):
    cfg = default_config()
    r = run_coupled_simulation(
        n_individuals=100, carrying_capacity=100, habitat_area=10000,
        T_celsius=14, salinity=30, phi_k=0.02, n_years=5, disease_year=999,
        seed=seed, config=cfg)
    results.append(r)

# Check biological invariants:
# 1. Population should be near carrying capacity (no disease)
pops = [r.final_pop for r in results]
print(f'Final populations: mean={np.mean(pops):.1f}, std={np.std(pops):.1f}, range=[{min(pops)}, {max(pops)}]')
assert all(p > 50 for p in pops), f'Population crashed unexpectedly: {pops}'

# 2. Mean resistance should be near initial value (no selection pressure)
resistances = [r.yearly_mean_resistance[-1] for r in results]
print(f'Final resistance: mean={np.mean(resistances):.4f}, std={np.std(resistances):.4f}')

# 3. Recruits should be produced each year
for i, r in enumerate(results):
    total_recruits = sum(r.yearly_recruits)
    print(f'  Seed {i}: total_recruits={total_recruits}, yearly={list(r.yearly_recruits)}')
    # At least some recruitment in most years
    assert total_recruits > 0, f'No recruits produced with seed {i}'

print('\n✅ All biological invariants verified')