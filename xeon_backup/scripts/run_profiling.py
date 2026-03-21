#!/usr/bin/env python3
"""Component-level profiling using the built-in PerfMonitor.

Runs multiple configurations to isolate where time is actually spent.
"""
import sys, time, json
sys.path.insert(0, '.')

import numpy as np
from sswd_evoepi.model import run_coupled_simulation
from sswd_evoepi.config import default_config
from sswd_evoepi.perf import PerfMonitor


def run_config(label, n_ind, n_years, spawning, disease_year, seed=42):
    """Run a single configuration and return timing data."""
    cfg = default_config()
    if not spawning:
        cfg.spawning = None

    perf = PerfMonitor(enabled=True)

    t0 = time.perf_counter()
    result = run_coupled_simulation(
        n_individuals=n_ind,
        carrying_capacity=n_ind,
        habitat_area=10000.0,
        T_celsius=14.0,
        salinity=30.0,
        phi_k=0.02,
        n_years=n_years,
        disease_year=disease_year,
        initial_infected=5,
        seed=seed,
        config=cfg,
        perf=perf,
    )
    wall = time.perf_counter() - t0

    summary = perf.summary()
    print(f"\n{'='*60}")
    print(f" {label}")
    print(f" Wall time: {wall:.3f}s | Final pop: {result.final_pop}")
    print(f"{'='*60}")
    print(perf.report())

    return {'label': label, 'wall_s': round(wall, 4), 'components': summary,
            'final_pop': result.final_pop, 'n_ind': n_ind, 'n_years': n_years,
            'spawning': spawning, 'disease_year': disease_year}


def main():
    results = []

    # Config 1: Baseline — no spawning, no disease
    results.append(run_config(
        "Baseline (no spawning, no disease)", 50, 2, False, 999))

    # Config 2: Disease only
    results.append(run_config(
        "Disease only (year 1)", 50, 2, False, 1))

    # Config 3: Spawning only — THIS IS THE BOTTLENECK
    results.append(run_config(
        "Spawning only (no disease)", 50, 2, True, 999))

    # Config 4: Everything
    results.append(run_config(
        "Full (spawning + disease)", 50, 2, True, 1))

    # Config 5: Larger population, spawning
    results.append(run_config(
        "Spawning 200 agents 2yr", 200, 2, True, 999))

    # Config 6: Longer sim, spawning
    results.append(run_config(
        "Spawning 50 agents 5yr", 50, 5, True, 999))

    # Save results
    with open('results/performance/profiling_data.json', 'w') as f:
        json.dump(results, f, indent=2)

    print("\n\n" + "="*60)
    print(" SUMMARY")
    print("="*60)
    for r in results:
        print(f"  {r['label']:<40} {r['wall_s']:>8.3f}s")
    print()

    # Identify biggest component across spawning runs
    spawning_runs = [r for r in results if r['spawning']]
    if spawning_runs:
        print("Top components in spawning runs:")
        for r in spawning_runs:
            comps = r['components']
            top = sorted(
                [(k, v) for k, v in comps.items() if k != '_total_s'],
                key=lambda x: -x[1].get('total_s', 0)
            )[:3]
            print(f"  {r['label']}:")
            for name, data in top:
                print(f"    {name}: {data['total_s']:.3f}s ({data['pct']:.1f}%)")


if __name__ == '__main__':
    main()
