#!/usr/bin/env python3
"""Phase 4: Final validation benchmarks after optimization.

1. Benchmark 5 spatial runs (3-node, 20yr, PE on, disease yr 3)
2. SA timing projections
3. Trajectory validation (10yr, no disease, daily recording)
"""
import os
import sys
import time
import json
from pathlib import Path

project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root))
os.chdir(project_root)

import numpy as np

from sswd_evoepi.config import SimulationConfig, default_config
from sswd_evoepi.spatial import NodeDefinition, build_network
from sswd_evoepi.model import run_spatial_simulation, run_coupled_simulation


def make_3node_network(seed=42):
    """Create 3-node benchmark network (Sitka / Howe Sound / Monterey)."""
    node_defs = [
        NodeDefinition(
            node_id=0, name="Sitka",
            lat=57.05, lon=-135.33, subregion="AK-SE",
            habitat_area=500_000.0, carrying_capacity=5000, is_fjord=False,
            mean_sst=8.5, sst_amplitude=3.5, sst_trend=0.02,
            salinity=32.0, flushing_rate=0.5,
        ),
        NodeDefinition(
            node_id=1, name="Howe Sound",
            lat=49.38, lon=-123.25, subregion="BC-SW",
            habitat_area=300_000.0, carrying_capacity=5000, is_fjord=True,
            sill_depth=40.0, mean_sst=10.0, sst_amplitude=3.0, sst_trend=0.02,
            salinity=28.0, flushing_rate=0.02,
        ),
        NodeDefinition(
            node_id=2, name="Monterey",
            lat=36.62, lon=-121.90, subregion="CA-C",
            habitat_area=400_000.0, carrying_capacity=5000, is_fjord=False,
            mean_sst=13.0, sst_amplitude=2.5, sst_trend=0.03,
            salinity=33.5, flushing_rate=0.5,
        ),
    ]
    return build_network(node_defs, seed=seed)


def benchmark_spatial():
    """Benchmark 5 spatial runs with different seeds."""
    seeds = [42, 123, 999, 7, 314]
    times = []

    cfg = default_config()
    cfg.pathogen_evolution.enabled = True

    for seed in seeds:
        network = make_3node_network(seed=seed)

        t0 = time.perf_counter()
        result = run_spatial_simulation(
            network=network,
            n_years=20,
            disease_year=3,
            initial_infected_per_node=5,
            seed=seed,
            config=cfg,
        )
        elapsed = time.perf_counter() - t0
        times.append(elapsed)

        final_pop = result.final_total_pop
        print(f"  Seed {seed:>3d}: {elapsed:.1f}s  (final pop: {final_pop})")

    mean_t = np.mean(times)
    min_t = np.min(times)
    max_t = np.max(times)

    print(f"\n  Mean: {mean_t:.1f}s  Min: {min_t:.1f}s  Max: {max_t:.1f}s")
    print(f"  vs pre-optimization: 107s → {mean_t:.1f}s ({(mean_t/107 - 1)*100:+.1f}%)")
    print(f"  vs pre-daily-mortality: 51s → {mean_t:.1f}s ({(mean_t/51 - 1)*100:+.1f}%)")

    return {'seeds': seeds, 'times': [round(t, 1) for t in times],
            'mean': round(mean_t, 1), 'min': round(min_t, 1), 'max': round(max_t, 1)}


def trajectory_validation():
    """10yr coupled sim, no disease, daily recording — check smoothness."""
    cfg = default_config()
    cfg.disease.enabled = False

    result = run_coupled_simulation(
        n_individuals=5000,
        carrying_capacity=5000,
        n_years=10,
        disease_year=999,  # effectively no disease
        seed=42,
        record_daily=True,
        config=cfg,
    )

    daily_pop = result.daily_pop
    if daily_pop is None or len(daily_pop) == 0:
        print("  ERROR: No daily records!")
        return {'passed': False, 'error': 'no daily records'}

    pops = np.array(daily_pop)

    # Compute day-over-day changes
    diffs = np.diff(pops)
    pct_changes = np.abs(diffs) / np.maximum(pops[:-1], 1) * 100

    max_drop_pct = float(np.max(pct_changes))
    max_drop_day = int(np.argmax(pct_changes))

    # Check year boundaries specifically
    year_boundary_days = [365 * y for y in range(1, 10)]
    year_boundary_drops = []
    for d in year_boundary_days:
        if d < len(pct_changes):
            year_boundary_drops.append(float(pct_changes[d]))

    max_year_boundary = max(year_boundary_drops) if year_boundary_drops else 0.0

    passed = max_drop_pct < 3.0

    print(f"  Max single-day drop: {max_drop_pct:.2f}% (day {max_drop_day}) — {'PASS' if passed else 'FAIL'} (limit: 3%)")
    print(f"  Max year-boundary change: {max_year_boundary:.2f}%")
    print(f"  Pop range: {int(pops.min())} - {int(pops.max())}")
    print(f"  Final pop: {int(pops[-1])}")

    return {
        'passed': passed,
        'max_drop_pct': round(max_drop_pct, 2),
        'max_drop_day': max_drop_day,
        'max_year_boundary_pct': round(max_year_boundary, 2),
        'pop_range': [int(pops.min()), int(pops.max())],
        'final_pop': int(pops[-1]),
    }


def sa_projections(mean_runtime):
    """Compute SA timing projections for Round 3."""
    n_params = 39
    n_cores = 8

    # Morris: (p+1) * r  where r = number of trajectories
    morris_traj = 21
    morris_runs = (n_params + 1) * morris_traj  # 40 * 21 = 840

    # Sobol: (2p+2) * N
    sobol_256_runs = (2 * n_params + 2) * 256  # 80 * 256 = 20,480
    sobol_128_runs = (2 * n_params + 2) * 128  # 80 * 128 = 10,240

    morris_hours = (morris_runs * mean_runtime) / n_cores / 3600
    sobol_256_hours = (sobol_256_runs * mean_runtime) / n_cores / 3600
    sobol_128_hours = (sobol_128_runs * mean_runtime) / n_cores / 3600

    print(f"\n  {'Method':<20s} {'Runs':>8s} {'Time/run':>10s} {'Total (8 cores)':>16s}")
    print(f"  {'─'*20} {'─'*8} {'─'*10} {'─'*16}")
    print(f"  {'Morris (21 traj)':<20s} {morris_runs:>8d} {mean_runtime:>9.1f}s {morris_hours:>12.1f}h")
    print(f"  {'Sobol N=128':<20s} {sobol_128_runs:>8d} {mean_runtime:>9.1f}s {sobol_128_hours:>12.1f}h")
    print(f"  {'Sobol N=256':<20s} {sobol_256_runs:>8d} {mean_runtime:>9.1f}s {sobol_256_hours:>12.1f}h")

    return {
        'n_params': n_params,
        'n_cores': n_cores,
        'morris': {'runs': morris_runs, 'hours': round(morris_hours, 1)},
        'sobol_128': {'runs': sobol_128_runs, 'hours': round(sobol_128_hours, 1)},
        'sobol_256': {'runs': sobol_256_runs, 'hours': round(sobol_256_hours, 1)},
    }


if __name__ == '__main__':
    print("=" * 60)
    print("PHASE 4: Final Validation & SA Timing")
    print(f"  {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)

    print("\n1. Spatial Benchmarks (3-node, 20yr, PE on, disease yr 3)")
    print("-" * 60)
    bench = benchmark_spatial()

    print("\n2. SA Timing Projections (39 params, 8 cores)")
    print("-" * 60)
    sa = sa_projections(bench['mean'])

    print("\n3. Trajectory Validation (10yr, no disease, K=5000)")
    print("-" * 60)
    traj = trajectory_validation()

    # Save results
    results = {
        'benchmark': bench,
        'sa_projections': sa,
        'trajectory_validation': traj,
        'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    }

    os.makedirs('results/continuous_mortality', exist_ok=True)
    with open('results/continuous_mortality/phase4_results.json', 'w') as f:
        json.dump(results, f, indent=2)

    print(f"\n{'=' * 60}")
    print("Results saved to results/continuous_mortality/phase4_results.json")
    print(f"{'=' * 60}")
