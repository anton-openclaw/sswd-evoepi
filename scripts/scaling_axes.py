#!/usr/bin/env python3
"""Measure scaling behavior along each independent axis of the model.

Axes:
  1. N — agents per node
  2. T — simulation years  
  3. L — genetic loci (estimated, hardcoded at 52)
  4. Spawning on/off overhead
  5. Disease on/off overhead

Each axis is swept while holding others constant.
"""
import sys, time, json
sys.path.insert(0, '.')

import numpy as np
from sswd_evoepi.model import run_coupled_simulation
from sswd_evoepi.config import default_config
from sswd_evoepi.perf import PerfMonitor


def timed_run(n, years, spawning=True, disease_year=999, seed=42, timeout=120):
    """Run simulation and return wall time + component breakdown."""
    cfg = default_config()
    if not spawning:
        cfg.spawning = None
    
    perf = PerfMonitor(enabled=True)
    t0 = time.perf_counter()
    try:
        r = run_coupled_simulation(
            n_individuals=n, carrying_capacity=n,
            habitat_area=max(10000, n * 200),
            T_celsius=14.0, salinity=30.0, phi_k=0.02,
            n_years=years, disease_year=disease_year,
            initial_infected=min(5, n // 10),
            seed=seed, config=cfg, perf=perf,
        )
        wall = time.perf_counter() - t0
        return {
            'wall_s': round(wall, 4),
            'final_pop': r.final_pop,
            'components': perf.summary(),
            'ok': True,
        }
    except Exception as e:
        wall = time.perf_counter() - t0
        return {'wall_s': round(wall, 4), 'error': str(e), 'ok': False}


def axis_N():
    """Sweep population size."""
    print("\n=== AXIS 1: Population Size (N) ===")
    print(f"{'N':>6} {'Time (s)':>10} {'Pop':>6} {'Top Component':>25}")
    print("-" * 55)
    
    results = []
    for n in [25, 50, 100, 200, 500, 1000, 2000]:
        r = timed_run(n, years=5, spawning=True)
        if r['ok']:
            # Find top component
            comps = r['components']
            top = max(
                [(k, v) for k, v in comps.items() if k != '_total_s'],
                key=lambda x: x[1].get('total_s', 0),
                default=('none', {'total_s': 0, 'pct': 0})
            )
            print(f"{n:>6} {r['wall_s']:>10.3f} {r['final_pop']:>6} {top[0]:>20} ({top[1]['pct']:.0f}%)")
        else:
            print(f"{n:>6} {r['wall_s']:>10.3f}   ERR  {r.get('error','')[:30]}")
        results.append({'n': n, **r})
        sys.stdout.flush()
    
    return results


def axis_T():
    """Sweep simulation duration."""
    print("\n=== AXIS 2: Simulation Duration (T years) ===")
    print(f"{'T':>6} {'Time (s)':>10} {'Pop':>6}")
    print("-" * 30)
    
    results = []
    for t in [1, 2, 5, 10, 20, 50]:
        r = timed_run(200, years=t, spawning=True)
        if r['ok']:
            print(f"{t:>6} {r['wall_s']:>10.3f} {r['final_pop']:>6}")
        else:
            print(f"{t:>6} {r['wall_s']:>10.3f}   ERR")
        results.append({'years': t, **r})
        sys.stdout.flush()
    
    return results


def axis_spawning_overhead():
    """Compare spawning on vs off at different scales."""
    print("\n=== AXIS 4: Spawning Overhead ===")
    print(f"{'N':>6} {'No Spawn':>10} {'Spawn':>10} {'Overhead':>10}")
    print("-" * 42)
    
    results = []
    for n in [50, 100, 200, 500, 1000]:
        r_off = timed_run(n, years=5, spawning=False)
        r_on = timed_run(n, years=5, spawning=True)
        overhead = r_on['wall_s'] / r_off['wall_s'] if r_off['wall_s'] > 0 else float('inf')
        print(f"{n:>6} {r_off['wall_s']:>10.3f} {r_on['wall_s']:>10.3f} {overhead:>9.1f}×")
        results.append({'n': n, 'no_spawn': r_off['wall_s'], 'spawn': r_on['wall_s'], 'overhead': round(overhead, 2)})
        sys.stdout.flush()
    
    return results


def axis_disease_overhead():
    """Compare disease on vs off."""
    print("\n=== AXIS 5: Disease Overhead ===")
    print(f"{'N':>6} {'No Disease':>10} {'Disease':>10} {'Overhead':>10}")
    print("-" * 42)
    
    results = []
    for n in [50, 200, 500]:
        r_off = timed_run(n, years=10, spawning=True, disease_year=999)
        r_on = timed_run(n, years=10, spawning=True, disease_year=3)
        overhead = r_on['wall_s'] / r_off['wall_s'] if r_off['wall_s'] > 0 else float('inf')
        print(f"{n:>6} {r_off['wall_s']:>10.3f} {r_on['wall_s']:>10.3f} {overhead:>9.1f}×")
        results.append({'n': n, 'no_disease': r_off['wall_s'], 'disease': r_on['wall_s'], 'overhead': round(overhead, 2)})
        sys.stdout.flush()
    
    return results


def main():
    all_results = {}
    
    all_results['axis_N'] = axis_N()
    all_results['axis_T'] = axis_T()
    all_results['axis_spawning'] = axis_spawning_overhead()
    all_results['axis_disease'] = axis_disease_overhead()
    
    # Save raw data
    with open('results/performance/scaling_axes_data.json', 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    
    # Summary
    print("\n" + "=" * 60)
    print(" SCALING SUMMARY")
    print("=" * 60)
    
    # Fit power law to N scaling
    n_data = [(r['n'], r['wall_s']) for r in all_results['axis_N'] if r['ok']]
    if len(n_data) >= 3:
        ns, ts = zip(*n_data)
        log_ns = np.log(ns)
        log_ts = np.log(ts)
        slope, intercept = np.polyfit(log_ns, log_ts, 1)
        print(f"  N scaling: O(N^{slope:.2f})")
    
    # Fit power law to T scaling
    t_data = [(r['years'], r['wall_s']) for r in all_results['axis_T'] if r['ok']]
    if len(t_data) >= 3:
        ts_val, times = zip(*t_data)
        log_ts = np.log(ts_val)
        log_times = np.log(times)
        slope_t, _ = np.polyfit(log_ts, log_times, 1)
        print(f"  T scaling: O(T^{slope_t:.2f})")
    
    # Projections
    print("\n  Projections (spawning enabled, no disease):")
    if n_data:
        # Use largest measured point to project
        n_ref, t_ref = n_data[-1]
        for target_n, target_t in [(500, 20), (1000, 20), (500, 50), (1000, 50), (2000, 20)]:
            est = t_ref * (target_n / n_ref) ** slope * (target_t / 5) ** slope_t
            unit = 's' if est < 60 else 'min' if est < 3600 else 'hr'
            val = est if est < 60 else est/60 if est < 3600 else est/3600
            print(f"    {target_n:>5} agents × {target_t:>3} yr: ~{val:.1f} {unit}")
    
    print(f"\n  Data saved: results/performance/scaling_axes_data.json")


if __name__ == '__main__':
    main()
