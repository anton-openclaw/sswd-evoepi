#!/usr/bin/env python3
"""Re-profile after optimization."""
import sys, time, json; sys.path.insert(0, '.')
from sswd_evoepi.model import run_coupled_simulation
from sswd_evoepi.config import default_config
from sswd_evoepi.perf import PerfMonitor

configs = [
    ('50 agents, 2yr', 50, 2),
    ('100 agents, 2yr', 100, 2),
    ('50 agents, 5yr', 50, 5),
]

results = []
for label, n, yrs in configs:
    cfg = default_config()
    perf = PerfMonitor(enabled=True)
    t0 = time.perf_counter()
    r = run_coupled_simulation(
        n_individuals=n, carrying_capacity=n, habitat_area=10000,
        T_celsius=14, salinity=30, phi_k=0.02, n_years=yrs,
        disease_year=999, seed=42, config=cfg, perf=perf)
    wall = time.perf_counter() - t0
    print(f'\n{label}: {wall:.3f}s')
    print(perf.report())
    results.append({'label': label, 'wall_s': round(wall, 4), 'components': perf.summary()})

# BASELINE COMPARISON (from pre-optimization):
print('\n' + '='*60)
print('BASELINE (pre-optimization): 50 agents, 2yr = 22.00s')
print(f'OPTIMIZED: 50 agents, 2yr = {results[0]["wall_s"]:.3f}s')
speedup = 22.0 / results[0]['wall_s'] if results[0]['wall_s'] > 0 else float('inf')
print(f'SPEEDUP: {speedup:.1f}x')
print('='*60)

with open('results/performance/post_optimization_profile.json', 'w') as f:
    json.dump(results, f, indent=2)