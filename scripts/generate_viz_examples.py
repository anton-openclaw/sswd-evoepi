#!/usr/bin/env python3
"""Generate example population & demographics visualizations.

Runs single-node and spatial simulations, then calls every function
in sswd_evoepi.viz.population and saves PNGs.
"""

import matplotlib
matplotlib.use('Agg')

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import time
from pathlib import Path

import numpy as np

OUT_DIR = Path(__file__).resolve().parent.parent / 'results' / 'viz_examples' / 'population'
OUT_DIR.mkdir(parents=True, exist_ok=True)


def main():
    t0 = time.time()

    # ── 1. Single-node simulation ────────────────────────────────────
    print("Running single-node simulation (500 agents, 20 yr, disease yr 3) ...")
    from sswd_evoepi.config import default_config
    from sswd_evoepi.model import run_coupled_simulation, initialize_population, make_effect_sizes

    cfg = default_config()
    if hasattr(cfg, 'pathogen_evolution'):
        cfg.pathogen_evolution.enabled = True

    result = run_coupled_simulation(
        n_individuals=500,
        carrying_capacity=500,
        n_years=20,
        disease_year=3,
        seed=42,
        config=cfg,
        record_daily=True,
    )
    print(f"  Single-node done: initial={result.initial_pop}, final={result.final_pop}, "
          f"min={result.min_pop} (year {result.min_pop_year})")

    # We also need agent snapshots for the pyramid and sex-ratio plots.
    # Re-run a short sim and grab the agents array at the end.
    print("Running short sim for agent snapshot ...")
    from sswd_evoepi.model import (
        DAYS_PER_YEAR, annual_growth_and_aging, annual_natural_mortality,
    )
    # Just use a fresh init for the pyramid
    effect_sizes = make_effect_sizes(cfg.genetics.effect_size_seed)
    agents_snap, geno_snap = initialize_population(
        n_individuals=500,
        max_agents=1500,
        habitat_area=10000.0,
        effect_sizes=effect_sizes,
        pop_cfg=cfg.population,
        rng=np.random.default_rng(42),
        genetics_cfg=cfg.genetics,
    )

    # ── 2. Spatial simulation (3-node) ───────────────────────────────
    print("Running 3-node spatial simulation (K=5000 each, 20 yr) ...")
    from scripts.sensitivity.spatial_runner import build_network
    from sswd_evoepi.model import run_spatial_simulation

    network = build_network(rng_seed=42)
    spatial_result = run_spatial_simulation(
        network=network,
        n_years=20,
        disease_year=3,
        initial_infected_per_node=5,
        seed=42,
        config=cfg,
    )
    print(f"  Spatial done: initial={spatial_result.initial_total_pop}, "
          f"final={spatial_result.final_total_pop}")

    # ── 3. Generate all plots ────────────────────────────────────────
    from sswd_evoepi.viz.population import (
        plot_population_trajectory,
        plot_stage_composition,
        plot_cause_of_death_breakdown,
        plot_age_size_pyramid,
        plot_population_heatmap,
        plot_recruitment_timeseries,
        plot_survival_curves,
        plot_sex_ratio_over_time,
        plot_density_dependence,
        plot_node_comparison_bars,
    )
    import matplotlib.pyplot as plt

    plots = []

    # 1. Population trajectory
    p = str(OUT_DIR / '01_population_trajectory.png')
    print(f"  Generating {Path(p).name} ...")
    plot_population_trajectory(result, carrying_capacity=500, disease_year=3, save_path=p)
    plots.append(p)

    # 2. Stage composition
    p = str(OUT_DIR / '02_stage_composition.png')
    print(f"  Generating {Path(p).name} ...")
    plot_stage_composition(result, save_path=p)
    plots.append(p)

    # 3. Cause of death
    p = str(OUT_DIR / '03_cause_of_death.png')
    print(f"  Generating {Path(p).name} ...")
    plot_cause_of_death_breakdown(result, save_path=p)
    plots.append(p)

    # 4. Age-size pyramid
    p = str(OUT_DIR / '04_age_size_pyramid.png')
    print(f"  Generating {Path(p).name} ...")
    plot_age_size_pyramid(agents_snap, title='Initial Population Pyramid', save_path=p)
    plots.append(p)

    # 5. Population heatmap (spatial)
    p = str(OUT_DIR / '05_population_heatmap.png')
    print(f"  Generating {Path(p).name} ...")
    plot_population_heatmap(spatial_result, save_path=p)
    plots.append(p)

    # 6. Recruitment timeseries
    p = str(OUT_DIR / '06_recruitment_timeseries.png')
    print(f"  Generating {Path(p).name} ...")
    plot_recruitment_timeseries(result, save_path=p)
    plots.append(p)

    # 7. Survival curves
    p = str(OUT_DIR / '07_survival_curves.png')
    print(f"  Generating {Path(p).name} ...")
    plot_survival_curves(result, save_path=p)
    plots.append(p)

    # 8. Sex ratio
    p = str(OUT_DIR / '08_sex_ratio.png')
    print(f"  Generating {Path(p).name} ...")
    plot_sex_ratio_over_time(result, save_path=p)
    plots.append(p)

    # 9. Density dependence
    p = str(OUT_DIR / '09_density_dependence.png')
    print(f"  Generating {Path(p).name} ...")
    plot_density_dependence(result, carrying_capacity=500, save_path=p)
    plots.append(p)

    # 10. Node comparison bars (spatial)
    p = str(OUT_DIR / '10_node_comparison.png')
    print(f"  Generating {Path(p).name} ...")
    plot_node_comparison_bars(spatial_result, save_path=p)
    plots.append(p)

    # ── 4. Verify all files exist ────────────────────────────────────
    print("\n--- Verification ---")
    all_ok = True
    for p in plots:
        exists = os.path.isfile(p)
        size = os.path.getsize(p) if exists else 0
        status = f"OK ({size:,} bytes)" if exists else "MISSING"
        print(f"  {Path(p).name}: {status}")
        if not exists:
            all_ok = False

    elapsed = time.time() - t0
    print(f"\nTotal time: {elapsed:.1f}s")
    print(f"All {len(plots)} plots generated: {'✅ SUCCESS' if all_ok else '❌ FAILURES'}")

    return 0 if all_ok else 1


if __name__ == '__main__':
    sys.exit(main())
