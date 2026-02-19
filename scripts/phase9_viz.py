#!/usr/bin/env python3
"""Phase 9: Generate settlement & spawning visualizations.

Runs a fresh 5-node 5-year simulation with daily spawning tracking enabled,
then generates all 9 new settlement/spawning visualizations.

Author: Anton 🔬
"""

import sys
import os
import time
import numpy as np
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from sswd_evoepi.model import (
    run_coupled_simulation,
    run_spatial_simulation,
    default_config,
    CoupledSimResult,
    SpatialSimResult,
)
from sswd_evoepi.spatial import make_5node_network

# Output directory
OUT_DIR = Path('results/continuous_settlement/viz')
OUT_DIR.mkdir(parents=True, exist_ok=True)


def run_spatial_sim():
    """Run 5-node 5-year spatial sim with spawning tracking."""
    print("Running 5-node 5-year spatial simulation...")
    t0 = time.time()
    config = default_config()
    network = make_5node_network(seed=42)
    result = run_spatial_simulation(
        network, n_years=5, disease_year=3, seed=42, config=config,
    )
    elapsed = time.time() - t0
    print(f"  Spatial sim done in {elapsed:.1f}s")
    print(f"  daily_spawning_counts shape: {result.daily_spawning_counts.shape}")
    print(f"  Total spawners tracked: {np.sum(result.daily_spawning_counts)}")
    for i in range(result.n_nodes):
        max_val = np.max(result.daily_spawning_counts[i])
        total = np.sum(result.daily_spawning_counts[i])
        print(f"    {result.node_names[i]}: peak {max_val}/day, total {total}")
    return result


def run_coupled_sim():
    """Run coupled sim for before/after comparison."""
    print("Running coupled K=5000 10-year simulation (daily recording)...")
    t0 = time.time()
    config = default_config()
    result = run_coupled_simulation(
        n_individuals=5000, carrying_capacity=5000,
        n_years=10, disease_year=3, seed=42,
        config=config, record_daily=True,
    )
    elapsed = time.time() - t0
    print(f"  Coupled sim done in {elapsed:.1f}s")
    print(f"  daily_spawning_counts shape: {result.daily_spawning_counts.shape}")
    print(f"  Total spawners: {np.sum(result.daily_spawning_counts)}")
    print(f"  Max spawners/day: {np.max(result.daily_spawning_counts)}")
    return result


def generate_all_viz(spatial_result, coupled_result):
    """Generate all 9 settlement/spawning visualizations."""
    from sswd_evoepi.viz.settlement import (
        plot_settlement_timing_heatmap,
        plot_settlement_spread,
        plot_pld_temperature_curve,
        plot_daily_recruitment,
        plot_before_after_epidemic,
        plot_spawning_intensity,
        plot_spawning_heatmap,
        plot_spawning_density_dependence,
        plot_spawning_cascade,
    )

    figures = []

    # 1. Settlement timing heatmap
    print("  [1/9] Settlement timing heatmap...")
    fig = plot_settlement_timing_heatmap(
        spatial_result, save_path=str(OUT_DIR / 'settlement_timing_heatmap.png'))
    figures.append('settlement_timing_heatmap.png')

    # 2. Settlement spread comparison
    print("  [2/9] Settlement spread comparison...")
    fig = plot_settlement_spread(
        spatial_result, save_path=str(OUT_DIR / 'settlement_spread.png'))
    figures.append('settlement_spread.png')

    # 3. PLD vs temperature curve
    print("  [3/9] PLD temperature curve...")
    fig = plot_pld_temperature_curve(
        save_path=str(OUT_DIR / 'pld_temperature_curve.png'))
    figures.append('pld_temperature_curve.png')

    # 4. Daily recruitment timeseries
    print("  [4/9] Daily recruitment timeseries...")
    fig = plot_daily_recruitment(
        spatial_result, save_path=str(OUT_DIR / 'daily_recruitment.png'))
    figures.append('daily_recruitment.png')

    # 5. Before/after epidemic curve (KEY FIGURE)
    print("  [5/9] Before/after epidemic curve...")
    fig = plot_before_after_epidemic(
        coupled_result, save_path=str(OUT_DIR / 'before_after_epidemic.png'))
    figures.append('before_after_epidemic.png')

    # 6. Daily spawning intensity
    print("  [6/9] Spawning intensity...")
    fig = plot_spawning_intensity(
        spatial_result, save_path=str(OUT_DIR / 'spawning_intensity.png'))
    figures.append('spawning_intensity.png')

    # 7. Spawning event heatmap
    print("  [7/9] Spawning heatmap...")
    fig = plot_spawning_heatmap(
        spatial_result, save_path=str(OUT_DIR / 'spawning_heatmap.png'))
    figures.append('spawning_heatmap.png')

    # 8. Spawning density dependence (KEY FIGURE)
    print("  [8/9] Spawning density dependence...")
    fig = plot_spawning_density_dependence(
        spatial_result, save_path=str(OUT_DIR / 'spawning_density_dependence.png'))
    figures.append('spawning_density_dependence.png')

    # 9. Spawning cascade timeline
    print("  [9/9] Spawning cascade timeline...")
    # Find node with highest K for best cascade signal
    best_node = int(np.argmax(spatial_result.node_K))
    fig = plot_spawning_cascade(
        spatial_result, node_idx=best_node,
        save_path=str(OUT_DIR / 'spawning_cascade.png'))
    figures.append('spawning_cascade.png')

    return figures


def main():
    print("=" * 60)
    print("Phase 9: Settlement & Spawning Visualizations")
    print("=" * 60)

    # Run simulations
    spatial_result = run_spatial_sim()
    coupled_result = run_coupled_sim()

    # Generate all visualizations
    print("\nGenerating visualizations...")
    figures = generate_all_viz(spatial_result, coupled_result)

    # Verify all outputs
    print(f"\n{'=' * 60}")
    print(f"Generated {len(figures)} figures in {OUT_DIR}/")
    all_ok = True
    for f in figures:
        path = OUT_DIR / f
        if path.exists():
            size_kb = path.stat().st_size / 1024
            print(f"  ✅ {f} ({size_kb:.1f} KB)")
        else:
            print(f"  ❌ {f} MISSING")
            all_ok = False

    if all_ok:
        print(f"\n✅ All {len(figures)} figures generated successfully!")
    else:
        print(f"\n❌ Some figures missing!")
        sys.exit(1)


if __name__ == '__main__':
    main()
