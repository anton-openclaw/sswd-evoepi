#!/usr/bin/env python3
"""Phase 14 patch: Generate viz that need CoupledSimResult or network."""

import sys, os, time, traceback
import numpy as np
from pathlib import Path

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))
os.chdir(ROOT)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from sswd_evoepi.model import (
    run_coupled_simulation,
    run_spatial_simulation,
    default_config,
)
from sswd_evoepi.spatial import make_5node_network

FINAL_VIZ = ROOT / 'results' / 'continuous_settlement' / 'final_viz'
FINAL_VIZ.mkdir(parents=True, exist_ok=True)

succeeded = []
failed = []

def save(name, func, *args, **kwargs):
    path = str(FINAL_VIZ / name)
    try:
        func(*args, save_path=path, **kwargs)
        size = os.path.getsize(path)
        succeeded.append((name, size))
        print(f"  ✅ {name} ({size/1024:.0f} KB)")
        plt.close('all')
    except Exception as e:
        failed.append((name, str(e)))
        print(f"  ❌ {name}: {e}")
        plt.close('all')


def main():
    t0 = time.time()
    config = default_config()

    # Coupled sim for single-site viz
    print("Running coupled sim K=5000, 10yr, daily recording...")
    coupled = run_coupled_simulation(
        n_individuals=5000, carrying_capacity=5000,
        n_years=10, disease_year=3, seed=42,
        config=config, record_daily=True,
    )
    print(f"  Done. daily_pop shape: {coupled.daily_pop.shape if coupled.daily_pop is not None else 'None'}")

    # Network for spatial viz
    network = make_5node_network(seed=42)

    # Spatial sim
    print("Running 5-node 20yr spatial sim...")
    spatial = run_spatial_simulation(
        network, n_years=20, disease_year=3, seed=42, config=config,
    )
    print(f"  Done.")

    # ─── Population viz (CoupledSimResult) ───
    from sswd_evoepi.viz.population import (
        plot_population_trajectory,
        plot_sex_ratio_over_time,
        plot_survival_curves,
        plot_recruitment_timeseries,
        plot_density_dependence,
        plot_cause_of_death_breakdown,
    )
    save('10_population_trajectory.png', plot_population_trajectory, coupled)
    save('13_sex_ratio.png', plot_sex_ratio_over_time, coupled)
    save('14_survival_curves.png', plot_survival_curves, coupled)
    save('15_recruitment_timeseries.png', plot_recruitment_timeseries, coupled)
    save('16_density_dependence.png', plot_density_dependence, coupled)
    save('19_cause_of_death.png', plot_cause_of_death_breakdown, coupled)

    # ─── Disease viz (CoupledSimResult) ───
    from sswd_evoepi.viz.disease import (
        plot_epidemic_curve,
        plot_cfr_over_time,
        plot_R0_over_time,
        plot_shedding_timeseries,
        plot_vibrio_concentration,
        plot_compartment_flow_sankey,
        plot_immunosuppression_overlap,
    )
    save('20_epidemic_curve.png', plot_epidemic_curve, coupled)
    save('22_cfr_over_time.png', plot_cfr_over_time, coupled)
    save('23_R0_over_time.png', plot_R0_over_time, coupled)
    save('24_shedding_timeseries.png', plot_shedding_timeseries, coupled)
    save('25_vibrio_concentration.png', plot_vibrio_concentration, coupled)
    save('31_compartment_flow_sankey.png', plot_compartment_flow_sankey, coupled)

    # ─── Genetics viz (CoupledSimResult) ───
    from sswd_evoepi.viz.genetics import (
        plot_resistance_trajectory,
        plot_additive_variance_over_time,
        plot_ef1a_dynamics,
    )
    save('32_resistance_trajectory.png', plot_resistance_trajectory, coupled)
    save('37_additive_variance.png', plot_additive_variance_over_time, coupled)
    save('38_ef1a_dynamics.png', plot_ef1a_dynamics, coupled)

    # ─── Coevolution viz (CoupledSimResult) ───
    from sswd_evoepi.viz.coevolution import (
        plot_virulence_trajectory,
        plot_coevolution_phase_portrait,
        plot_virulence_distribution_over_time,
        plot_virulence_vs_host_density,
    )
    save('44_virulence_trajectory.png', plot_virulence_trajectory, coupled)
    save('46_coevo_phase_portrait.png', plot_coevolution_phase_portrait, coupled)
    save('47_virulence_distribution.png', plot_virulence_distribution_over_time, coupled)
    save('48_virulence_vs_density.png', plot_virulence_vs_host_density, coupled)

    # ─── Spatial viz (need network object) ───
    from sswd_evoepi.viz.spatial import (
        plot_connectivity_heatmap,
        plot_larval_flow_diagram,
        plot_network_map,
    )
    save('51_connectivity_heatmap.png', plot_connectivity_heatmap, network)
    save('56_larval_flow_diagram.png', plot_larval_flow_diagram, spatial, network)
    save('57_network_map.png', plot_network_map, network)

    # ─── Dashboards (CoupledSimResult) ───
    from sswd_evoepi.viz.dashboards import (
        plot_simulation_dashboard,
        plot_evolutionary_rescue_assessment,
        plot_model_validation_panel,
    )
    save('58_simulation_dashboard.png', plot_simulation_dashboard, coupled)
    save('60_evolutionary_rescue.png', plot_evolutionary_rescue_assessment, coupled)
    save('61_validation_panel.png', plot_model_validation_panel, coupled)

    # ─── Tradeoff curve (needs config) ───
    from sswd_evoepi.viz.coevolution import plot_tradeoff_curve
    save('45_tradeoff_curve.png', plot_tradeoff_curve, coupled, config.disease)

    total = time.time() - t0
    print(f"\n{'='*60}")
    print(f"  Patch: {len(succeeded)} succeeded, {len(failed)} failed — {total:.1f}s")
    print(f"{'='*60}")
    if failed:
        for n, e in failed:
            print(f"    ❌ {n}: {e}")


if __name__ == '__main__':
    main()
