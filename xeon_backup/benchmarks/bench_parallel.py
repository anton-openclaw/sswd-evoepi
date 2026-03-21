#!/usr/bin/env python3
"""Benchmark parallel node processing.

Runs a 5-node, 3-year simulation at workers=1,2,4,8
and reports wall-clock times.
"""

import time
import sys
import os

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from sswd_evoepi.spatial import build_network, NodeDefinition
from sswd_evoepi.config import default_config
from sswd_evoepi.model import run_spatial_simulation


def make_network(n_nodes=5):
    """Create a test network."""
    defs = []
    for i in range(n_nodes):
        defs.append(NodeDefinition(
            node_id=i,
            name=f"bench_{i}",
            lat=48.0 + i * 0.5,
            lon=-123.0 + i * 0.5,
            subregion="bench",
            habitat_area=1e4,
            carrying_capacity=100,
            mean_sst=12.0,
            sst_amplitude=4.0,
            sst_trend=0.0,
            salinity=30.0,
            flushing_rate=0.05,
            is_fjord=False,
        ))
    cfg = default_config()
    return build_network(defs, D_L=cfg.spatial.D_L, D_P=cfg.spatial.D_P,
                         r_total=cfg.spatial.r_total)


def benchmark(n_nodes=5, n_years=3, workers_list=None, seed=42):
    """Run benchmark across different worker counts."""
    if workers_list is None:
        workers_list = [1, 2, 4, 8]

    results = {}
    for w in workers_list:
        network = make_network(n_nodes)
        config = default_config()
        config.simulation.parallel_workers = w
        config.movement.enabled = True
        config.movement.spatial_transmission = True

        t0 = time.perf_counter()
        result = run_spatial_simulation(
            network=network,
            n_years=n_years,
            disease_year=1,
            initial_infected_per_node=5,
            seed=seed,
            config=config,
        )
        elapsed = time.perf_counter() - t0

        results[w] = {
            'elapsed': elapsed,
            'final_pop': result.final_total_pop,
            'total_disease_deaths': int(result.yearly_disease_deaths.sum()),
        }
        print(f"  workers={w:2d}  time={elapsed:6.2f}s  "
              f"pop={result.final_total_pop}  "
              f"disease_deaths={int(result.yearly_disease_deaths.sum())}")

    return results


if __name__ == "__main__":
    print(f"Benchmark: 5 nodes, 3 years, movement+disease enabled")
    print(f"{'='*60}")
    results = benchmark(n_nodes=5, n_years=3)

    print(f"\n{'='*60}")
    print("Summary:")
    serial_time = results[1]['elapsed']
    for w, r in results.items():
        speedup = serial_time / r['elapsed'] if r['elapsed'] > 0 else 0
        print(f"  workers={w:2d}: {r['elapsed']:6.2f}s  "
              f"speedup={speedup:.2f}x")
