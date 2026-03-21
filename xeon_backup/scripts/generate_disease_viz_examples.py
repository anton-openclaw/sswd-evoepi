#!/usr/bin/env python3
"""Generate example disease & epidemiology visualizations.

Runs single-node and spatial simulations, then calls every function
in sswd_evoepi.viz.disease and saves PNGs.
"""

import matplotlib
matplotlib.use('Agg')

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import time
from pathlib import Path

import numpy as np

OUT_DIR = Path(__file__).resolve().parent.parent / 'results' / 'viz_examples' / 'disease'
OUT_DIR.mkdir(parents=True, exist_ok=True)


def main():
    t0 = time.time()

    # ── 1. Single-node simulation (500 agents, 20yr, disease yr 3) ───
    print("Running single-node simulation (500 agents, 20 yr, disease yr 3, PE on) ...")
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
          f"disease_deaths={result.total_disease_deaths}")

    # Agent snapshot for FoI distribution and recovery plots
    print("Running short sim for agent snapshot with partial disease ...")
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

    # Run a short epidemic to get some recovered agents
    from sswd_evoepi.disease import (
        run_single_node_epidemic, NodeDiseaseState, daily_disease_update,
        sample_stage_duration, arrhenius, K_SHAPE_E
    )
    from sswd_evoepi.types import DiseaseState

    # Seed a few infections
    rng = np.random.default_rng(42)
    alive_idx = np.where(agents_snap['alive'])[0]
    n_infect = min(10, len(alive_idx))
    infect_idx = rng.choice(alive_idx, size=n_infect, replace=False)
    mu_EI1 = arrhenius(cfg.disease.mu_EI1_ref, cfg.disease.Ea_EI1, 15.0)
    for idx in infect_idx:
        agents_snap['disease_state'][idx] = DiseaseState.E
        agents_snap['disease_timer'][idx] = sample_stage_duration(mu_EI1, K_SHAPE_E, rng)

    # Run 90 days of disease to get some through the pipeline
    node_state = NodeDiseaseState(node_id=0, vibrio_concentration=5000.0)
    for day in range(90):
        node_state = daily_disease_update(
            agents=agents_snap, node_state=node_state,
            T_celsius=15.0, salinity=30.0, phi_k=0.02,
            dispersal_input=0.0, day=day, cfg=cfg.disease, rng=rng,
        )
    
    n_alive = int(np.sum(agents_snap['alive']))
    n_recovered = int(np.sum(agents_snap['alive'] & (agents_snap['disease_state'] == DiseaseState.R)))
    print(f"  Agent snapshot: {n_alive} alive, {n_recovered} recovered, "
          f"vibrio={node_state.vibrio_concentration:.0f}")

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

    # ── 3. Generate all disease plots ────────────────────────────────
    from sswd_evoepi.viz.disease import (
        plot_epidemic_curve,
        plot_vibrio_concentration,
        plot_force_of_infection_distribution,
        plot_R0_over_time,
        plot_disease_mortality_by_node,
        plot_epidemic_wave_timing,
        plot_compartment_flow_sankey,
        plot_shedding_timeseries,
        plot_disease_state_heatmap,
        plot_immunosuppression_overlap,
        plot_recovery_vs_resistance,
        plot_cfr_over_time,
    )

    plots = []

    # 1. Epidemic curve (THE SIGNATURE PLOT)
    p = str(OUT_DIR / '01_epidemic_curve.png')
    print(f"  Generating {Path(p).name} ...")
    plot_epidemic_curve(result, disease_year=3, save_path=p)
    plots.append(p)

    # 2. Vibrio concentration
    p = str(OUT_DIR / '02_vibrio_concentration.png')
    print(f"  Generating {Path(p).name} ...")
    plot_vibrio_concentration(result, epidemic_threshold=1000.0, disease_year=3, save_path=p)
    plots.append(p)

    # 3. Force of infection distribution
    p = str(OUT_DIR / '03_force_of_infection.png')
    print(f"  Generating {Path(p).name} ...")
    plot_force_of_infection_distribution(
        agents_snap, vibrio_concentration=5000.0, salinity=30.0, save_path=p
    )
    plots.append(p)

    # 4. R₀ over time
    p = str(OUT_DIR / '04_R0_over_time.png')
    print(f"  Generating {Path(p).name} ...")
    plot_R0_over_time(result, T_celsius=15.0, salinity=30.0, phi_k=0.02,
                      disease_year=3, save_path=p)
    plots.append(p)

    # 5. Disease mortality by node (spatial)
    p = str(OUT_DIR / '05_disease_mortality_by_node.png')
    print(f"  Generating {Path(p).name} ...")
    plot_disease_mortality_by_node(spatial_result, save_path=p)
    plots.append(p)

    # 6. Epidemic wave timing (spatial)
    p = str(OUT_DIR / '06_epidemic_wave_timing.png')
    print(f"  Generating {Path(p).name} ...")
    plot_epidemic_wave_timing(spatial_result, save_path=p)
    plots.append(p)

    # 7. Compartment flow Sankey
    p = str(OUT_DIR / '07_compartment_flow_sankey.png')
    print(f"  Generating {Path(p).name} ...")
    plot_compartment_flow_sankey(result, save_path=p)
    plots.append(p)

    # 8. Shedding timeseries
    p = str(OUT_DIR / '08_shedding_timeseries.png')
    print(f"  Generating {Path(p).name} ...")
    plot_shedding_timeseries(result, T_celsius=15.0, disease_year=3, save_path=p)
    plots.append(p)

    # 9. Disease state heatmap (spatial)
    p = str(OUT_DIR / '09_disease_state_heatmap.png')
    print(f"  Generating {Path(p).name} ...")
    plot_disease_state_heatmap(spatial_result, save_path=p)
    plots.append(p)

    # 10. Immunosuppression overlap
    p = str(OUT_DIR / '10_immunosuppression_overlap.png')
    print(f"  Generating {Path(p).name} ...")
    plot_immunosuppression_overlap(T_celsius_range=(7.0, 18.0), latitude=48.0,
                                    save_path=p)
    plots.append(p)

    # 11. Recovery vs resistance
    p = str(OUT_DIR / '11_recovery_vs_resistance.png')
    print(f"  Generating {Path(p).name} ...")
    plot_recovery_vs_resistance(agents=agents_snap, rho_rec=0.05, save_path=p)
    plots.append(p)

    # 12. CFR over time
    p = str(OUT_DIR / '12_cfr_over_time.png')
    print(f"  Generating {Path(p).name} ...")
    plot_cfr_over_time(result, disease_year=3, save_path=p)
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
