#!/usr/bin/env python3
"""Run 5-node network simulation with CRW movement + spatial transmission.

Compare runtime and results against the mean-field baseline.
"""

import time
import sys
import numpy as np

sys.path.insert(0, '.')

from sswd_evoepi.spatial import make_5node_network
from sswd_evoepi.model import run_spatial_simulation
from sswd_evoepi.config import default_config
from sswd_evoepi.snapshots import SnapshotRecorder

N_YEARS = 15
SEED = 42
DISEASE_YEAR = 5

def run_sim(label, movement_enabled, spatial_tx, snapshot=False):
    """Run one simulation and report results."""
    config = default_config()
    config.movement.enabled = movement_enabled
    config.movement.spatial_transmission = spatial_tx

    network = make_5node_network()
    
    snap = None
    if snapshot:
        snap = SnapshotRecorder(
            enabled=True,
            interval_days=7,
            nodes=[n.definition.node_id for n in network.nodes],
        )

    print(f"\n{'='*60}")
    print(f"  {label}")
    print(f"  Movement: {movement_enabled}, Spatial TX: {spatial_tx}")
    print(f"  {N_YEARS} years, seed={SEED}, disease at year {DISEASE_YEAR}")
    print(f"{'='*60}")

    t0 = time.time()
    result = run_spatial_simulation(
        network=network,
        n_years=N_YEARS,
        disease_year=DISEASE_YEAR,
        initial_infected_per_node=5,
        seed=SEED,
        config=config,
        snapshot_recorder=snap,
        progress_callback=lambda y, n: print(f"  Year {y+1}/{n}", end='\r'),
    )
    elapsed = time.time() - t0
    print(f"\n  Runtime: {elapsed:.1f}s")

    # Report per-node results
    print(f"\n  {'Node':<15} {'Start':>6} {'End':>6} {'Crash%':>7} {'Mean r':>7}")
    print(f"  {'-'*42}")
    for i, name in enumerate(result.node_names):
        pop_start = result.yearly_pop[i, DISEASE_YEAR - 1]
        pop_end = result.yearly_pop[i, -1]
        crash_pct = (1 - pop_end / max(pop_start, 1)) * 100
        mean_r = result.yearly_mean_resistance[i, -1]
        print(f"  {name:<15} {pop_start:>6} {pop_end:>6} {crash_pct:>6.1f}% {mean_r:>7.4f}")

    total_start = sum(result.yearly_pop[i, DISEASE_YEAR - 1] for i in range(len(result.node_names)))
    total_end = sum(result.yearly_pop[i, -1] for i in range(len(result.node_names)))
    print(f"  {'TOTAL':<15} {total_start:>6} {total_end:>6} {(1-total_end/max(total_start,1))*100:>6.1f}%")

    if snap and len(snap.snapshots) > 0:
        print(f"\n  Snapshots: {len(snap.snapshots)} captures")

    return result, elapsed

# Run both configurations
print("SSWD-EvoEpi 5-Node Simulation: Movement Comparison")
print(f"Python {sys.version.split()[0]}")

r_mf, t_mf = run_sim("MEAN-FIELD (no movement)", False, False)
r_mv, t_mv = run_sim("CRW + SPATIAL TRANSMISSION", True, True, snapshot=True)

# Compare
print(f"\n{'='*60}")
print(f"  COMPARISON")
print(f"{'='*60}")
print(f"  Mean-field runtime:  {t_mf:.1f}s")
print(f"  Movement runtime:    {t_mv:.1f}s")
print(f"  Slowdown factor:     {t_mv/t_mf:.1f}x")

print(f"\n  Population comparison (final year):")
for i, name in enumerate(r_mf.node_names):
    mf_pop = r_mf.yearly_pop[i, -1]
    mv_pop = r_mv.yearly_pop[i, -1]
    mf_r = r_mf.yearly_mean_resistance[i, -1]
    mv_r = r_mv.yearly_mean_resistance[i, -1]
    print(f"  {name:<15} MF:{mf_pop:>5}  Mov:{mv_pop:>5}  Δr:{mv_r-mf_r:>+.4f}")

print("\nDone.")
