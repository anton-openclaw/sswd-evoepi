#!/usr/bin/env python3
"""Generate example visualizations: co-evolution + spatial (Part 4).

Runs simulations and produces all 16 PNGs.
"""

import sys
import os
import time

# Ensure the project root is on the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
from pathlib import Path

from sswd_evoepi.config import default_config, SimulationConfig
from sswd_evoepi.model import run_coupled_simulation, run_spatial_simulation
from sswd_evoepi.spatial import (
    build_network, NodeDefinition, get_5node_definitions,
)

# Import viz modules
from sswd_evoepi.viz.coevolution import (
    plot_virulence_trajectory,
    plot_coevolution_phase_portrait,
    plot_virulence_distribution_over_time,
    plot_tradeoff_curve,
    plot_R0_by_virulence,
    plot_virulence_vs_host_density,
    plot_strain_competition,
    plot_coevolution_multi_seed,
)
from sswd_evoepi.viz.spatial import (
    plot_network_map,
    plot_connectivity_heatmap,
    plot_north_south_gradient,
    plot_fjord_vs_open,
    plot_metapopulation_timeseries,
    plot_larval_flow_diagram,
    plot_spatial_epidemic_timeline,
    plot_node_fate_matrix,
)


def make_pe_config() -> SimulationConfig:
    """Config with pathogen evolution enabled."""
    cfg = default_config()
    cfg.pathogen_evolution.enabled = True
    cfg.pathogen_evolution.v_init = 0.5
    cfg.pathogen_evolution.sigma_v_mutation = 0.03
    return cfg


def make_3node_defs():
    """3-node network for spatial examples: Sitka, Howe Sound, Monterey."""
    all5 = get_5node_definitions()
    # Indices 0 (Sitka), 1 (Howe Sound), 4 (Monterey) — renumber
    nodes = []
    for new_id, old_idx in enumerate([0, 1, 4]):
        nd = all5[old_idx]
        nodes.append(NodeDefinition(
            node_id=new_id,
            name=nd.name,
            lat=nd.lat, lon=nd.lon,
            subregion=nd.subregion,
            habitat_area=nd.habitat_area,
            carrying_capacity=nd.carrying_capacity,
            is_fjord=nd.is_fjord,
            sill_depth=nd.sill_depth,
            flushing_rate=nd.flushing_rate,
            mean_sst=nd.mean_sst,
            sst_amplitude=nd.sst_amplitude,
            sst_trend=nd.sst_trend,
            salinity=nd.salinity,
            depth_range=nd.depth_range,
        ))
    return nodes


def run_single_coevo(seed=42, n_ind=500, n_years=20):
    """Single-node co-evolution run."""
    cfg = make_pe_config()
    print(f"  Running single-node co-evo: {n_ind} agents, {n_years}yr, seed={seed}...")
    t0 = time.time()
    result = run_coupled_simulation(
        n_individuals=n_ind,
        carrying_capacity=n_ind,
        n_years=n_years,
        disease_year=3,
        seed=seed,
        config=cfg,
    )
    print(f"  Done in {time.time() - t0:.1f}s. Final pop: {result.final_pop}")
    return result


def run_multi_seed_coevo(seeds=range(42, 47), n_ind=200, n_years=15):
    """Multi-seed co-evolution runs."""
    cfg = make_pe_config()
    results = {}
    for s in seeds:
        print(f"  Multi-seed co-evo: seed={s}...")
        t0 = time.time()
        r = run_coupled_simulation(
            n_individuals=n_ind,
            carrying_capacity=n_ind,
            n_years=n_years,
            disease_year=3,
            seed=s,
            config=cfg,
        )
        results[s] = r
        print(f"    Done in {time.time() - t0:.1f}s. Final pop: {r.final_pop}")
    return results


def run_spatial_3node(n_years=20):
    """3-node spatial simulation with PE enabled."""
    cfg = make_pe_config()
    node_defs = make_3node_defs()
    print(f"  Building 3-node network...")
    network = build_network(node_defs, seed=42)
    print(f"  Running 3-node spatial sim: {n_years}yr...")
    t0 = time.time()
    result = run_spatial_simulation(
        network=network,
        n_years=n_years,
        disease_year=3,
        initial_infected_per_node=5,
        seed=42,
        config=cfg,
    )
    print(f"  Done in {time.time() - t0:.1f}s. Final total pop: {result.final_total_pop}")
    return result, network


def main():
    # Output directories
    coevo_dir = Path('results/viz_examples/coevolution')
    spatial_dir = Path('results/viz_examples/spatial')
    coevo_dir.mkdir(parents=True, exist_ok=True)
    spatial_dir.mkdir(parents=True, exist_ok=True)

    cfg = make_pe_config()

    # ─── 1. TRADEOFF CURVE (no sim needed) ────────────────────────────
    print("\n=== Tradeoff curve (no sim) ===")
    plot_tradeoff_curve(
        cfg.pathogen_evolution, cfg.disease,
        T_celsius=15.0, v_current=0.5,
        save_path=str(coevo_dir / 'tradeoff_curve.png'),
    )
    print("  ✓ tradeoff_curve.png")

    # ─── 2. R₀ BY VIRULENCE (no sim needed) ──────────────────────────
    print("\n=== R₀ by virulence (no sim) ===")
    plot_R0_by_virulence(
        cfg.disease, cfg.pathogen_evolution,
        T_celsius=15.0,
        save_path=str(coevo_dir / 'R0_by_virulence.png'),
    )
    print("  ✓ R0_by_virulence.png")

    # ─── 3. SINGLE-NODE CO-EVO ────────────────────────────────────────
    print("\n=== Single-node co-evolution simulation ===")
    result = run_single_coevo()

    print("\n  Generating co-evolution plots...")
    plot_virulence_trajectory(
        result, disease_year=3,
        v_init=cfg.pathogen_evolution.v_init,
        save_path=str(coevo_dir / 'virulence_trajectory.png'),
    )
    print("  ✓ virulence_trajectory.png")

    plot_coevolution_phase_portrait(
        result, disease_year=3,
        save_path=str(coevo_dir / 'coevolution_phase_portrait.png'),
    )
    print("  ✓ coevolution_phase_portrait.png")

    plot_virulence_distribution_over_time(
        result, disease_year=3,
        save_path=str(coevo_dir / 'virulence_distribution_over_time.png'),
    )
    print("  ✓ virulence_distribution_over_time.png")

    plot_virulence_vs_host_density(
        result, disease_year=3,
        save_path=str(coevo_dir / 'virulence_vs_host_density.png'),
    )
    print("  ✓ virulence_vs_host_density.png")

    plot_strain_competition(
        result, disease_year=3,
        save_path=str(coevo_dir / 'strain_competition.png'),
    )
    print("  ✓ strain_competition.png")

    # ─── 4. MULTI-SEED CO-EVO ─────────────────────────────────────────
    print("\n=== Multi-seed co-evolution ===")
    multi_results = run_multi_seed_coevo()

    plot_coevolution_multi_seed(
        multi_results, disease_year=3,
        save_path=str(coevo_dir / 'coevolution_multi_seed.png'),
    )
    print("  ✓ coevolution_multi_seed.png")

    # ─── 5. SPATIAL SIMULATION ────────────────────────────────────────
    print("\n=== 3-node spatial simulation ===")
    spatial_result, network = run_spatial_3node()

    print("\n  Generating spatial plots...")

    # Network map — colour by final population fraction
    final_pop = spatial_result.yearly_pop[:, -1].astype(float)
    K = np.array([n.definition.carrying_capacity for n in network.nodes]).astype(float)
    pop_fraction = final_pop / K
    plot_network_map(
        network, metric_by_node=pop_fraction,
        title='Final Population Fraction by Node',
        metric_label='Pop / K',
        save_path=str(spatial_dir / 'network_map.png'),
    )
    print("  ✓ network_map.png")

    plot_connectivity_heatmap(
        network, matrix='C',
        save_path=str(spatial_dir / 'connectivity_heatmap.png'),
    )
    print("  ✓ connectivity_heatmap.png")

    lats = np.array([n.definition.lat for n in network.nodes])
    plot_north_south_gradient(
        spatial_result, metric='mortality', lats=lats,
        save_path=str(spatial_dir / 'north_south_gradient_mortality.png'),
    )
    print("  ✓ north_south_gradient_mortality.png")

    plot_fjord_vs_open(
        spatial_result, fjord_idx=1, disease_year=3,
        save_path=str(spatial_dir / 'fjord_vs_open.png'),
    )
    print("  ✓ fjord_vs_open.png")

    plot_metapopulation_timeseries(
        spatial_result, disease_year=3,
        save_path=str(spatial_dir / 'metapopulation_timeseries.png'),
    )
    print("  ✓ metapopulation_timeseries.png")

    plot_larval_flow_diagram(
        spatial_result, network,
        save_path=str(spatial_dir / 'larval_flow_diagram.png'),
    )
    print("  ✓ larval_flow_diagram.png")

    plot_spatial_epidemic_timeline(
        spatial_result, disease_year=3,
        save_path=str(spatial_dir / 'spatial_epidemic_timeline.png'),
    )
    print("  ✓ spatial_epidemic_timeline.png")

    plot_node_fate_matrix(
        spatial_result, disease_year=3,
        save_path=str(spatial_dir / 'node_fate_matrix.png'),
    )
    print("  ✓ node_fate_matrix.png")

    # ─── Summary ──────────────────────────────────────────────────────
    coevo_pngs = list(coevo_dir.glob('*.png'))
    spatial_pngs = list(spatial_dir.glob('*.png'))
    print(f"\n{'='*60}")
    print(f"  Co-evolution plots: {len(coevo_pngs)}")
    for p in sorted(coevo_pngs):
        print(f"    {p.name} ({p.stat().st_size / 1024:.0f} KB)")
    print(f"  Spatial plots: {len(spatial_pngs)}")
    for p in sorted(spatial_pngs):
        print(f"    {p.name} ({p.stat().st_size / 1024:.0f} KB)")
    print(f"  Total: {len(coevo_pngs) + len(spatial_pngs)} PNGs")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
