#!/usr/bin/env python3
"""Validate recruitment pipeline fix.

Compares population recovery dynamics with old vs new s0 values.
Runs a single 5-node, 10-year simulation with disease crash at year 1.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sswd_evoepi.config import default_config, SimulationConfig
from sswd_evoepi.spatial import NodeDefinition, build_network
from sswd_evoepi.model import run_spatial_simulation

def run_scenario(s0: float, label: str):
    """Run a 5-node sim with given s0 and report yearly populations."""
    config = default_config()
    config.population.settler_survival = s0
    config.simulation.seed = 42
    
    # 5 nodes spanning temperature gradient
    nodes = []
    for i, (lat, lon, name) in enumerate([
        (57.0, -135.0, "Alaska"),
        (52.0, -128.0, "BC-N"),
        (48.5, -123.0, "WA"),
        (44.0, -124.0, "OR"),
        (38.0, -123.0, "CA-N"),
    ]):
        nd = NodeDefinition(
            node_id=i, name=name, lat=lat, lon=lon, subregion=name,
            habitat_area=25000.0, carrying_capacity=5000,
            is_fjord=False, sill_depth=np.inf, flushing_rate=0.5,
            mean_sst=max(4.0, 25.0 - 0.35*lat),
            sst_amplitude=4.0, sst_trend=0.02, salinity=32.0,
            depth_range=(5.0, 60.0),
        )
        nodes.append(nd)
    
    network = build_network(nodes, D_L=400.0, D_P=15.0, seed=42)
    
    result = run_spatial_simulation(
        network=network, n_years=10, disease_year=1, seed=42, config=config,
    )
    
    print(f"\n=== {label} (s0={s0}) ===")
    print(f"{'Year':>5} {'Total Pop':>10} {'% of K':>8} {'Δ%':>8}")
    total_K = 5 * 5000
    for y in range(result.yearly_pop.shape[1]):
        pop = result.yearly_pop[:, y].sum()
        pct = pop / total_K * 100
        delta = (pop / max(result.yearly_pop[:, max(0,y-1)].sum(), 1) - 1) * 100 if y > 0 else 0
        print(f"{y:5d} {int(pop):10d} {pct:7.1f}% {delta:+7.1f}%")
    
    return result

if __name__ == "__main__":
    print("=== Recruitment Pipeline Validation ===")
    print("5-node network, K=5000/node, disease at year 1")
    
    # Run with new s0=0.001
    r_new = run_scenario(0.001, "NEW (s0=0.001)")
    
    # Run with old s0=0.03 for comparison
    r_old = run_scenario(0.03, "OLD (s0=0.03)")
    
    print("\n=== Summary ===")
    total_K = 25000
    for label, r in [("OLD s0=0.03", r_old), ("NEW s0=0.001", r_new)]:
        final_pct = r.yearly_pop[:, -1].sum() / total_K * 100
        min_pct = min(r.yearly_pop[:, y].sum() for y in range(r.yearly_pop.shape[1])) / total_K * 100
        print(f"{label}: min={min_pct:.1f}%, final={final_pct:.1f}%")
