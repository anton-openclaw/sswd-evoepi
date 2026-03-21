#!/usr/bin/env python3
"""Phase 14: Generate ALL visualizations + final report data.

Runs a fresh 5-node 20yr simulation with all new features, then generates
every visualization in the library.

Author: Anton 🔬
"""

import sys
import os
import time
import json
import traceback
import numpy as np
from pathlib import Path

# Add project root to path
ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))
os.chdir(ROOT)

from sswd_evoepi.model import (
    SimulationConfig,
    run_spatial_simulation,
    run_coupled_simulation,
    default_config,
)
from sswd_evoepi.spatial import load_node_definitions_yaml, build_network, make_5node_network

# Output
FINAL_VIZ = ROOT / 'results' / 'continuous_settlement' / 'final_viz'
FINAL_VIZ.mkdir(parents=True, exist_ok=True)

import matplotlib
matplotlib.use('Agg')

succeeded = []
failed = []

def save(name, func, *args, **kwargs):
    """Try to generate a figure, track success/failure."""
    path = str(FINAL_VIZ / name)
    try:
        func(*args, save_path=path, **kwargs)
        size = os.path.getsize(path)
        succeeded.append((name, size))
        print(f"  ✅ {name} ({size/1024:.0f} KB)")
    except Exception as e:
        failed.append((name, str(e)))
        print(f"  ❌ {name}: {e}")
        traceback.print_exc()


def main():
    t0 = time.time()
    print("=" * 70)
    print("  Phase 14: Comprehensive Visualization Generation")
    print("=" * 70)

    # ── 1. Run simulations ─────────────────────────────────────────
    print("\n[1/4] Running 5-node 20yr spatial simulation (seed=42)...")
    config = default_config()
    network = make_5node_network(seed=42)
    t_sim = time.time()
    spatial_result = run_spatial_simulation(
        network, n_years=20, disease_year=3, seed=42, config=config,
    )
    print(f"  Spatial sim: {time.time()-t_sim:.1f}s")

    print("\n[1b] Running coupled sim K=5000, 10yr (daily recording)...")
    t_sim = time.time()
    coupled_result = run_coupled_simulation(
        n_individuals=5000, carrying_capacity=5000,
        n_years=10, disease_year=3, seed=42,
        config=config, record_daily=True,
    )
    print(f"  Coupled sim: {time.time()-t_sim:.1f}s")

    # Save summary data
    summary = {
        'n_nodes': spatial_result.n_nodes,
        'node_names': list(spatial_result.node_names),
        'node_K': [int(k) for k in spatial_result.node_K],
        'initial_total_pop': int(spatial_result.initial_total_pop),
        'final_total_pop': int(spatial_result.final_total_pop),
    }
    for i in range(spatial_result.n_nodes):
        name = spatial_result.node_names[i]
        yp = spatial_result.yearly_pop[i]
        init = yp[0] if yp[0] > 0 else spatial_result.node_K[i]
        crash = (1 - yp[-1] / init) * 100 if init > 0 else 0
        summary[f'{name}_crash_pct'] = round(float(crash), 1)
        summary[f'{name}_final_pop'] = int(yp[-1])
        if hasattr(spatial_result, 'yearly_mean_resistance'):
            r = spatial_result.yearly_mean_resistance[i]
            summary[f'{name}_resist_shift'] = round(float(r[-1] - r[0]), 4)
        if hasattr(spatial_result, 'yearly_mean_virulence') and spatial_result.yearly_mean_virulence is not None:
            v = spatial_result.yearly_mean_virulence[i]
            summary[f'{name}_mean_virulence'] = round(float(np.nanmean(v[3:])), 4)
    with open(FINAL_VIZ / 'simulation_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"\n  Summary: {json.dumps(summary, indent=2)}")

    # ── 2. Settlement visualizations ───────────────────────────────
    print("\n[2/4] Settlement visualizations...")
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

    save('01_settlement_timing_heatmap.png', plot_settlement_timing_heatmap, spatial_result)
    save('02_settlement_spread.png', plot_settlement_spread, spatial_result)
    save('03_pld_temperature_curve.png', plot_pld_temperature_curve)
    save('04_daily_recruitment.png', plot_daily_recruitment, spatial_result)
    save('05_before_after_epidemic.png', plot_before_after_epidemic, coupled_result)
    save('06_spawning_intensity.png', plot_spawning_intensity, spatial_result)
    save('07_spawning_heatmap.png', plot_spawning_heatmap, spatial_result)
    save('08_spawning_density_dependence.png', plot_spawning_density_dependence, spatial_result)
    best_node = int(np.argmax(spatial_result.node_K))
    save('09_spawning_cascade.png', plot_spawning_cascade, spatial_result, node_idx=best_node)

    # ── 3. Population, disease, genetics, coevolution, spatial ─────
    print("\n[3/4] Full visualization library...")

    from sswd_evoepi.viz.population import (
        plot_population_trajectory,
        plot_population_heatmap,
        plot_stage_composition,
        plot_sex_ratio_over_time,
        plot_survival_curves,
        plot_recruitment_timeseries,
        plot_density_dependence,
        plot_node_comparison_bars,
        plot_age_size_pyramid,
        plot_cause_of_death_breakdown,
    )
    save('10_population_trajectory.png', plot_population_trajectory, spatial_result)
    save('11_population_heatmap.png', plot_population_heatmap, spatial_result)
    save('12_stage_composition.png', plot_stage_composition, spatial_result)
    save('13_sex_ratio.png', plot_sex_ratio_over_time, spatial_result)
    save('14_survival_curves.png', plot_survival_curves, spatial_result)
    save('15_recruitment_timeseries.png', plot_recruitment_timeseries, spatial_result)
    save('16_density_dependence.png', plot_density_dependence, spatial_result)
    save('17_node_comparison_bars.png', plot_node_comparison_bars, spatial_result)
    save('18_age_size_pyramid.png', plot_age_size_pyramid, spatial_result)
    save('19_cause_of_death.png', plot_cause_of_death_breakdown, spatial_result)

    from sswd_evoepi.viz.disease import (
        plot_epidemic_curve,
        plot_disease_state_heatmap,
        plot_cfr_over_time,
        plot_R0_over_time,
        plot_shedding_timeseries,
        plot_vibrio_concentration,
        plot_disease_mortality_by_node,
        plot_epidemic_wave_timing,
        plot_recovery_vs_resistance,
        plot_force_of_infection_distribution,
        plot_immunosuppression_overlap,
        plot_compartment_flow_sankey,
    )
    save('20_epidemic_curve.png', plot_epidemic_curve, spatial_result)
    save('21_disease_state_heatmap.png', plot_disease_state_heatmap, spatial_result)
    save('22_cfr_over_time.png', plot_cfr_over_time, spatial_result)
    save('23_R0_over_time.png', plot_R0_over_time, spatial_result)
    save('24_shedding_timeseries.png', plot_shedding_timeseries, spatial_result)
    save('25_vibrio_concentration.png', plot_vibrio_concentration, spatial_result)
    save('26_disease_mortality_by_node.png', plot_disease_mortality_by_node, spatial_result)
    save('27_epidemic_wave_timing.png', plot_epidemic_wave_timing, spatial_result)
    save('28_recovery_vs_resistance.png', plot_recovery_vs_resistance, spatial_result)
    save('29_force_of_infection.png', plot_force_of_infection_distribution, spatial_result)
    save('30_immunosuppression_overlap.png', plot_immunosuppression_overlap, spatial_result)
    save('31_compartment_flow_sankey.png', plot_compartment_flow_sankey, spatial_result)

    from sswd_evoepi.viz.genetics import (
        plot_resistance_trajectory,
        plot_resistance_distribution,
        plot_allele_frequency_spaghetti,
        plot_selection_differential,
        plot_heritability_over_time,
        plot_additive_variance_over_time,
        plot_ef1a_dynamics,
        plot_genotype_phenotype_map,
        plot_resistance_by_node_violin,
        plot_locus_effect_size_distribution,
        plot_genetic_drift_null,
        plot_beta_init_visualization,
    )
    save('32_resistance_trajectory.png', plot_resistance_trajectory, spatial_result)
    save('33_resistance_distribution.png', plot_resistance_distribution, spatial_result)
    save('34_allele_freq_spaghetti.png', plot_allele_frequency_spaghetti, spatial_result)
    save('35_selection_differential.png', plot_selection_differential, spatial_result)
    save('36_heritability.png', plot_heritability_over_time, spatial_result)
    save('37_additive_variance.png', plot_additive_variance_over_time, spatial_result)
    save('38_ef1a_dynamics.png', plot_ef1a_dynamics, spatial_result)
    save('39_genotype_phenotype_map.png', plot_genotype_phenotype_map, spatial_result)
    save('40_resistance_violin.png', plot_resistance_by_node_violin, spatial_result)
    save('41_locus_effect_sizes.png', plot_locus_effect_size_distribution, spatial_result)
    save('42_genetic_drift_null.png', plot_genetic_drift_null, spatial_result)
    save('43_beta_init.png', plot_beta_init_visualization)

    from sswd_evoepi.viz.coevolution import (
        plot_virulence_trajectory,
        plot_tradeoff_curve,
        plot_coevolution_phase_portrait,
        plot_virulence_distribution_over_time,
        plot_virulence_vs_host_density,
        plot_R0_by_virulence,
        plot_coevolution_multi_seed,
        plot_strain_competition,
    )
    save('44_virulence_trajectory.png', plot_virulence_trajectory, spatial_result)
    save('45_tradeoff_curve.png', plot_tradeoff_curve, spatial_result)
    save('46_coevo_phase_portrait.png', plot_coevolution_phase_portrait, spatial_result)
    save('47_virulence_distribution.png', plot_virulence_distribution_over_time, spatial_result)
    save('48_virulence_vs_density.png', plot_virulence_vs_host_density, spatial_result)
    save('49_R0_by_virulence.png', plot_R0_by_virulence, spatial_result)
    # multi-seed and strain competition need multiple results
    # Skip for now — they require running multiple seeds

    from sswd_evoepi.viz.spatial import (
        plot_metapopulation_timeseries,
        plot_connectivity_heatmap,
        plot_north_south_gradient,
        plot_fjord_vs_open,
        plot_spatial_epidemic_timeline,
        plot_node_fate_matrix,
        plot_larval_flow_diagram,
        plot_network_map,
    )
    save('50_metapopulation_timeseries.png', plot_metapopulation_timeseries, spatial_result)
    save('51_connectivity_heatmap.png', plot_connectivity_heatmap, spatial_result)
    save('52_north_south_gradient.png', plot_north_south_gradient, spatial_result)
    save('53_fjord_vs_open.png', plot_fjord_vs_open, spatial_result)
    save('54_spatial_epidemic_timeline.png', plot_spatial_epidemic_timeline, spatial_result)
    save('55_node_fate_matrix.png', plot_node_fate_matrix, spatial_result)
    save('56_larval_flow_diagram.png', plot_larval_flow_diagram, spatial_result)
    save('57_network_map.png', plot_network_map, spatial_result)

    from sswd_evoepi.viz.dashboards import (
        plot_simulation_dashboard,
        plot_spatial_dashboard,
        plot_sensitivity_tornado,
        plot_sensitivity_heatmap,
        plot_evolutionary_rescue_assessment,
        plot_model_validation_panel,
    )
    save('58_simulation_dashboard.png', plot_simulation_dashboard, spatial_result)
    save('59_spatial_dashboard.png', plot_spatial_dashboard, spatial_result)
    save('60_evolutionary_rescue.png', plot_evolutionary_rescue_assessment, spatial_result)
    save('61_validation_panel.png', plot_model_validation_panel, spatial_result)

    # ── 4. Spawning viz (Phase 13 type) ────────────────────────────
    print("\n[4/4] Spawning dynamics viz (static + animations)...")
    from sswd_evoepi.viz.spawning_viz import (
        run_simulation_with_snapshots,
        plot_spawning_event_profile,
        plot_spawning_participation,
        plot_readiness_cascade,
        plot_spawning_vs_density,
        plot_spawning_before_after,
        create_spawning_animation,
        create_comparison_animation,
    )

    # High density, new params
    print("  Running snapshot sim (high density, new params)...")
    snaps_high, _ = run_simulation_with_snapshots(
        n_individuals=500, carrying_capacity=500, habitat_area=10000.0,
        T_celsius=10.5, n_years=3, seed=42,
        snapshot_start_doy=None, snapshot_end_doy=None, snapshot_year=1,
        disease_year=None,
    )
    print(f"    {len(snaps_high)} snapshots")

    # Low density
    print("  Running snapshot sim (low density)...")
    snaps_low, _ = run_simulation_with_snapshots(
        n_individuals=50, carrying_capacity=500, habitat_area=10000.0,
        T_celsius=10.5, n_years=3, seed=42,
        snapshot_start_doy=None, snapshot_end_doy=None, snapshot_year=1,
        disease_year=None,
    )

    # Old params
    print("  Running snapshot sim (old params)...")
    old_overrides = {
        'female_max_bouts': 1, 'male_refractory_days': 21,
        'induction_male_to_female': 0.30, 'induction_female_to_male': 0.50,
        'readiness_induction_prob': 0.0, 'gravity_enabled': False,
    }
    snaps_old, _ = run_simulation_with_snapshots(
        n_individuals=500, carrying_capacity=500, habitat_area=10000.0,
        T_celsius=10.5, n_years=3, seed=42,
        snapshot_start_doy=None, snapshot_end_doy=None, snapshot_year=1,
        disease_year=None, spawning_overrides=old_overrides,
    )

    save('62_spawning_event_profile.png', plot_spawning_event_profile, snaps_high)
    save('63_readiness_cascade.png', plot_readiness_cascade, snaps_high)
    save('64_spawning_density_scatter.png', plot_spawning_vs_density,
         [('High density (500/K=500)', snaps_high, 500),
          ('Low density (50/K=500)', snaps_low, 500)])
    save('65_spawning_before_after.png', plot_spawning_before_after, snaps_old, snaps_high)

    # Window around peak
    hab_side = np.sqrt(10000.0)
    if snaps_high:
        peak_idx = max(range(len(snaps_high)),
                       key=lambda i: snaps_high[i].n_spawning_today)
        peak_doy = snaps_high[peak_idx].doy

        def doy_dist(a, b):
            d = abs(a - b)
            return min(d, 365 - d)

        snaps_high_w = [s for s in snaps_high if doy_dist(s.doy, peak_doy) <= 30]
        snaps_low_w = [s for s in snaps_low if doy_dist(s.doy, peak_doy) <= 30]
        snaps_old_w = [s for s in snaps_old if doy_dist(s.doy, peak_doy) <= 30]
    else:
        snaps_high_w = snaps_high
        snaps_low_w = snaps_low
        snaps_old_w = snaps_old

    # GIFs
    print("  Generating animated GIFs...")
    save('66_spawning_highdensity.gif', create_spawning_animation,
         snaps_high_w, hab_side=hab_side, fps=4, title_prefix='HIGH DENSITY | ', dpi=80)
    save('67_spawning_lowdensity.gif', create_spawning_animation,
         snaps_low_w, hab_side=hab_side, fps=4, title_prefix='LOW DENSITY | ', dpi=80)
    save('68_spawning_param_comparison.gif', create_comparison_animation,
         snaps_old_w, snaps_high_w, hab_side=hab_side, fps=4,
         label_a='Old Params', label_b='Phase 12 Params',
         max_female_bouts_a=1, max_male_bouts_a=3,
         max_female_bouts_b=2, max_male_bouts_b=3, dpi=80)
    save('69_spawning_full_season.gif', create_spawning_animation,
         snaps_high, hab_side=hab_side, fps=8, title_prefix='FULL SEASON | ', dpi=80)

    # ── Summary ────────────────────────────────────────────────────
    total_time = time.time() - t0
    print(f"\n{'='*70}")
    print(f"  COMPLETE: {len(succeeded)} succeeded, {len(failed)} failed — {total_time:.1f}s")
    print(f"{'='*70}")

    if failed:
        print(f"\n  FAILED:")
        for name, err in failed:
            print(f"    ❌ {name}: {err}")

    # Save manifest
    manifest = {
        'total_time': round(total_time, 1),
        'succeeded': len(succeeded),
        'failed': len(failed),
        'figures': [{'name': n, 'size_bytes': s} for n, s in succeeded],
        'errors': [{'name': n, 'error': e} for n, e in failed],
    }
    with open(FINAL_VIZ / 'manifest.json', 'w') as f:
        json.dump(manifest, f, indent=2)

    print(f"\n  Output: {FINAL_VIZ}")
    return len(failed) == 0


if __name__ == '__main__':
    ok = main()
    sys.exit(0 if ok else 1)
