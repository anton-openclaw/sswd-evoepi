#!/usr/bin/env python3
"""Scaling test: measure performance at different population sizes."""
import sys, time, json; sys.path.insert(0, '.')
from sswd_evoepi.model import run_coupled_simulation
from sswd_evoepi.config import default_config
from sswd_evoepi.perf import PerfMonitor

configs = [
    ('50 agents, 5yr', 50, 5),
    ('100 agents, 5yr', 100, 5),
    ('200 agents, 5yr', 200, 5),
    ('500 agents, 5yr', 500, 5),
    ('50 agents, 20yr', 50, 20),
    ('200 agents, 20yr', 200, 20),
]

results = []
for label, n, yrs in configs:
    cfg = default_config()
    perf = PerfMonitor(enabled=True)
    t0 = time.perf_counter()
    r = run_coupled_simulation(
        n_individuals=n, carrying_capacity=n, habitat_area=max(10000, n*200),
        T_celsius=14, salinity=30, phi_k=0.02, n_years=yrs,
        disease_year=999, seed=42, config=cfg, perf=perf)
    wall = time.perf_counter() - t0
    print(f'{label}: {wall:.3f}s  pop={r.final_pop}')
    summary = perf.summary()
    # Show top 3 components
    top = sorted([(k,v) for k,v in summary.items() if k != '_total_s'], key=lambda x: -x[1].get('total_s',0))[:3]
    for name, data in top:
        print(f'  {name}: {data["total_s"]:.3f}s ({data["pct"]:.1f}%)')
    results.append({'label': label, 'wall_s': round(wall, 4), 'n': n, 'years': yrs, 'components': summary, 'final_pop': r.final_pop})
    sys.stdout.flush()

print('\nScaling summary:')
for r in results:
    print(f'  {r["label"]:<25} {r["wall_s"]:>8.3f}s')

# Estimate 500-agent 20yr
if len(results) >= 4:
    # Simple extrapolation
    r200_5 = next(r for r in results if r['n']==200 and r['years']==5)
    est_500_20 = r200_5['wall_s'] * (500/200)**1.5 * (20/5)
    print(f'\nEstimate 500 agents 20yr: ~{est_500_20:.0f}s ({est_500_20/60:.1f}min)')

with open('results/performance/scaling_data.json', 'w') as f:
    json.dump(results, f, indent=2)