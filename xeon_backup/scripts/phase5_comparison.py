#!/usr/bin/env python3
"""Phase 5: Comparison run — 5-node 20yr spatial sim with disease.

Same config as final_5node_20yr to compare timing + trajectory smoothness.
"""
import time
import json
import numpy as np
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sswd_evoepi.config import default_config
from sswd_evoepi.spatial import build_network


def run_comparison(seed=42):
    """Run the 5-node 20yr spatial sim and return results + timing."""
    config = default_config()

    # 5-node network: Sitka, Howe Sound, SJI, Newport, Monterey
    from sswd_evoepi.spatial import NodeDefinition
    nodes = [
        NodeDefinition(node_id=0, name="Sitka", lat=57.05, lon=-135.33,
                       subregion="Alaska", carrying_capacity=5000,
                       mean_sst=8.5, salinity=31.0, habitat_area=100.0,
                       flushing_rate=0.15, is_fjord=False),
        NodeDefinition(node_id=1, name="Howe Sound", lat=49.38, lon=-123.25,
                       subregion="BC", carrying_capacity=5000,
                       mean_sst=10.5, salinity=28.0, habitat_area=80.0,
                       flushing_rate=0.08, is_fjord=True),
        NodeDefinition(node_id=2, name="SJI", lat=48.53, lon=-123.02,
                       subregion="WA", carrying_capacity=5000,
                       mean_sst=10.0, salinity=30.0, habitat_area=120.0,
                       flushing_rate=0.12, is_fjord=False),
        NodeDefinition(node_id=3, name="Newport", lat=44.63, lon=-124.05,
                       subregion="OR", carrying_capacity=5000,
                       mean_sst=11.5, salinity=33.0, habitat_area=90.0,
                       flushing_rate=0.20, is_fjord=False),
        NodeDefinition(node_id=4, name="Monterey", lat=36.60, lon=-121.90,
                       subregion="CA", carrying_capacity=5000,
                       mean_sst=13.0, salinity=33.5, habitat_area=110.0,
                       flushing_rate=0.18, is_fjord=False),
    ]
    network = build_network(
        nodes,
        alpha_self_fjord=config.spatial.alpha_self_fjord,
        alpha_self_open=config.spatial.alpha_self_open,
    )

    from sswd_evoepi.model import run_spatial_simulation
    t0 = time.perf_counter()
    result = run_spatial_simulation(
        network=network,
        n_years=20,
        disease_year=3,
        initial_infected_per_node=5,
        seed=seed,
        config=config,
        progress_callback=lambda y, n: print(f"  Year {y+1}/{n}"),
    )
    elapsed = time.perf_counter() - t0
    return result, elapsed


def main():
    print("=" * 60)
    print("Phase 5: 5-node 20yr comparison run (seed=42)")
    print("=" * 60)

    result, elapsed = run_comparison(seed=42)

    print(f"\n{'='*60}")
    print(f"Runtime: {elapsed:.1f}s (baseline was ~51s)")
    print(f"Initial total pop: {result.initial_total_pop}")
    print(f"Final total pop:   {result.final_total_pop}")
    print()

    # Per-node final populations
    for i, name in enumerate(result.node_names):
        final_pop = result.yearly_pop[i, -1]
        initial_pop = result.yearly_pop[i, 0] if result.yearly_pop.shape[1] > 0 else 0
        pct = (final_pop / initial_pop * 100) if initial_pop > 0 else 0
        print(f"  {name:15s}: {final_pop:6d} ({pct:5.1f}% of initial)")

    # Save results
    os.makedirs("results/continuous_mortality", exist_ok=True)
    summary = {
        "seed": 42,
        "runtime_s": round(elapsed, 2),
        "baseline_runtime_s": 51.0,
        "initial_total_pop": int(result.initial_total_pop),
        "final_total_pop": int(result.final_total_pop),
        "node_results": {},
    }
    for i, name in enumerate(result.node_names):
        summary["node_results"][name] = {
            "yearly_pop": result.yearly_pop[i].tolist(),
            "yearly_natural_deaths": result.yearly_natural_deaths[i].tolist(),
            "yearly_disease_deaths": result.yearly_disease_deaths[i].tolist(),
            "final_pop": int(result.yearly_pop[i, -1]),
        }

    with open("results/continuous_mortality/comparison_5node_20yr.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nResults saved to results/continuous_mortality/comparison_5node_20yr.json")

    # Multi-seed timing
    print(f"\n{'='*60}")
    print("Multi-seed timing (5 seeds)...")
    print("=" * 60)
    timings = [elapsed]  # seed=42 already done
    for s in [123, 456, 789, 1000]:
        print(f"  Running seed={s}...")
        _, t = run_comparison(seed=s)
        timings.append(t)
        print(f"    → {t:.1f}s")

    mean_t = np.mean(timings)
    min_t = np.min(timings)
    max_t = np.max(timings)
    print(f"\n  Mean: {mean_t:.1f}s | Range: [{min_t:.1f}s, {max_t:.1f}s]")

    summary["multi_seed_timings"] = {
        "seeds": [42, 123, 456, 789, 1000],
        "times_s": [round(t, 2) for t in timings],
        "mean_s": round(mean_t, 2),
        "min_s": round(min_t, 2),
        "max_s": round(max_t, 2),
    }
    with open("results/continuous_mortality/comparison_5node_20yr.json", "w") as f:
        json.dump(summary, f, indent=2)

    return summary


if __name__ == "__main__":
    main()
