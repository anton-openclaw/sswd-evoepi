#!/usr/bin/env python3
"""Generate dashboard & composite view example plots.

Runs actual simulations and loads SA data to produce all 10+ dashboard PNGs.
Output: results/viz_examples/dashboards/
"""

import json
import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
from sswd_evoepi.config import default_config, SimulationConfig
from sswd_evoepi.model import run_coupled_simulation, run_spatial_simulation
from sswd_evoepi.spatial import make_5node_network, NodeDefinition, build_network

OUT = Path('results/viz_examples/dashboards')
OUT.mkdir(parents=True, exist_ok=True)

# ═══════════════════════════════════════════════════════════════════════
# HELPER: make a 3-node network
# ═══════════════════════════════════════════════════════════════════════

def make_3node_network(seed=42):
    """Sitka (open, cold), Howe Sound (fjord), Monterey (open, warm)."""
    nodes = [
        NodeDefinition(
            node_id=0, name='Sitka, AK', lat=57.06, lon=-135.34,
            subregion='AK-SE', habitat_area=5000, carrying_capacity=500,
            is_fjord=False, mean_sst=8.5, sst_amplitude=3.5, salinity=32.0,
            flushing_rate=0.5,
        ),
        NodeDefinition(
            node_id=1, name='Howe Sound, BC', lat=49.52, lon=-123.25,
            subregion='BC-S', habitat_area=2000, carrying_capacity=400,
            is_fjord=True, sill_depth=60.0, mean_sst=10.0, sst_amplitude=3.0,
            salinity=28.0, flushing_rate=0.02,
        ),
        NodeDefinition(
            node_id=2, name='Monterey, CA', lat=36.62, lon=-121.90,
            subregion='CA-C', habitat_area=3500, carrying_capacity=500,
            is_fjord=False, mean_sst=13.5, sst_amplitude=2.5, salinity=33.5,
            flushing_rate=0.5,
        ),
    ]
    return build_network(nodes, seed=seed)


# ═══════════════════════════════════════════════════════════════════════
# 1. Single-node simulation (PE on)
# ═══════════════════════════════════════════════════════════════════════
print('[1/7] Single-node simulation (500 agents, 20yr, PE on, seed=42)...')
cfg_pe = default_config()
cfg_pe.pathogen_evolution.enabled = True
cfg_pe.pathogen_evolution.v_init = 0.5

result_pe = run_coupled_simulation(
    n_individuals=500, carrying_capacity=500,
    n_years=20, disease_year=3, seed=42,
    record_daily=True, config=cfg_pe,
)
print(f'  Pop: {result_pe.initial_pop} → {result_pe.final_pop}, '
      f'disease deaths: {result_pe.total_disease_deaths}')

# ═══════════════════════════════════════════════════════════════════════
# 2. Single-node no PE
# ═══════════════════════════════════════════════════════════════════════
print('[2/7] Single-node simulation (no PE)...')
cfg_nope = default_config()
cfg_nope.pathogen_evolution.enabled = False

result_nope = run_coupled_simulation(
    n_individuals=500, carrying_capacity=500,
    n_years=20, disease_year=3, seed=42,
    record_daily=True, config=cfg_nope,
)
print(f'  Pop: {result_nope.initial_pop} → {result_nope.final_pop}')

# ═══════════════════════════════════════════════════════════════════════
# 3. PE on v=0.8 (higher virulence)
# ═══════════════════════════════════════════════════════════════════════
print('[3/7] Single-node simulation (PE on, v=0.8)...')
cfg_pe8 = default_config()
cfg_pe8.pathogen_evolution.enabled = True
cfg_pe8.pathogen_evolution.v_init = 0.8

result_pe8 = run_coupled_simulation(
    n_individuals=500, carrying_capacity=500,
    n_years=20, disease_year=3, seed=42,
    record_daily=True, config=cfg_pe8,
)
print(f'  Pop: {result_pe8.initial_pop} → {result_pe8.final_pop}')

# ═══════════════════════════════════════════════════════════════════════
# 4. Spatial simulation (3-node, PE on)
# ═══════════════════════════════════════════════════════════════════════
print('[4/7] Spatial simulation (3-node, 20yr, PE on)...')
net3 = make_3node_network(seed=42)
spatial_result = run_spatial_simulation(
    network=net3, n_years=20, disease_year=3,
    initial_infected_per_node=5, seed=42, config=cfg_pe,
)
print(f'  Total pop: {spatial_result.initial_total_pop} → {spatial_result.final_total_pop}')

# ═══════════════════════════════════════════════════════════════════════
# 5. Load SA data
# ═══════════════════════════════════════════════════════════════════════
print('[5/7] Loading sensitivity analysis data...')
sobol_path = Path('results/sensitivity/sobol_indices.json')
morris_path = Path('results/sensitivity/morris_screening.json')

sobol_indices = None
morris_results = None

if sobol_path.exists():
    with open(sobol_path) as f:
        sobol_indices = json.load(f)
    print(f'  Sobol: {len(sobol_indices["param_names"])} params, '
          f'{len(sobol_indices["metric_names"])} metrics')

if morris_path.exists():
    with open(morris_path) as f:
        morris_results = json.load(f)
    print(f'  Morris: {morris_results["n_params"]} params, '
          f'{morris_results["n_runs"]} runs')

# ═══════════════════════════════════════════════════════════════════════
# 6. Generate all dashboard plots
# ═══════════════════════════════════════════════════════════════════════
print('[6/7] Generating dashboard plots...')

from sswd_evoepi.viz.dashboards import (
    plot_simulation_dashboard,
    plot_spatial_dashboard,
    plot_scenario_comparison,
    plot_sensitivity_tornado,
    plot_sensitivity_heatmap,
    plot_parameter_interaction_web,
    plot_evolutionary_rescue_assessment,
    plot_conservation_scenario_matrix,
    plot_model_validation_panel,
    plot_parameter_space_exploration,
)

plots_generated = []

# 1. Simulation dashboard (PE on)
print('  → simulation_dashboard_pe.png')
plot_simulation_dashboard(
    result_pe, title='Single-Node Dashboard (PE Enabled, v₀=0.5)',
    disease_year=3, carrying_capacity=500,
    save_path=str(OUT / 'simulation_dashboard_pe.png'),
)
plots_generated.append(('simulation_dashboard_pe.png',
    'Master 2×3 dashboard: population, epidemic, resistance, deaths, virulence, co-evolution (PE on, v₀=0.5)'))

# 2. Simulation dashboard (no PE)
print('  → simulation_dashboard_nope.png')
plot_simulation_dashboard(
    result_nope, title='Single-Node Dashboard (No Pathogen Evolution)',
    disease_year=3, carrying_capacity=500,
    save_path=str(OUT / 'simulation_dashboard_nope.png'),
)
plots_generated.append(('simulation_dashboard_nope.png',
    'Master dashboard without PE — Vibrio concentration panel instead of virulence'))

# 3. Spatial dashboard
print('  → spatial_dashboard.png')
plot_spatial_dashboard(
    spatial_result, title='3-Node Spatial Dashboard (Sitka / Howe Sound / Monterey)',
    network=net3,
    save_path=str(OUT / 'spatial_dashboard.png'),
)
plots_generated.append(('spatial_dashboard.png',
    'Multi-node 2×3 overview: per-node trajectories, mortality bar, network map, resistance, N-S gradient, fjord vs open'))

# 4. Scenario comparison (3 configs)
print('  → scenario_comparison.png')
scenarios = {
    'No PE': result_nope,
    'PE v₀=0.5': result_pe,
    'PE v₀=0.8': result_pe8,
}
plot_scenario_comparison(
    scenarios, title='Scenario Comparison: Effect of Pathogen Evolution',
    save_path=str(OUT / 'scenario_comparison.png'),
)
plots_generated.append(('scenario_comparison.png',
    'Side-by-side 4×3 grid: population, resistance, virulence, deaths for 3 PE scenarios'))

# 5-7. Sensitivity analysis plots
if sobol_indices is not None:
    # Tornado for pop_crash_pct
    print('  → sensitivity_tornado_popcrash.png')
    plot_sensitivity_tornado(
        sobol_indices, metric='pop_crash_pct', top_n=15,
        save_path=str(OUT / 'sensitivity_tornado_popcrash.png'),
    )
    plots_generated.append(('sensitivity_tornado_popcrash.png',
        'Sobol S1/ST tornado plot for population crash metric, top 15 params, color by category'))

    # Tornado for resistance_shift_mean
    print('  → sensitivity_tornado_resistance.png')
    plot_sensitivity_tornado(
        sobol_indices, metric='resistance_shift_mean', top_n=15,
        save_path=str(OUT / 'sensitivity_tornado_resistance.png'),
    )
    plots_generated.append(('sensitivity_tornado_resistance.png',
        'Sobol tornado for mean resistance shift metric'))

    # Heatmap
    print('  → sensitivity_heatmap.png')
    plot_sensitivity_heatmap(
        sobol_indices, top_n=20,
        save_path=str(OUT / 'sensitivity_heatmap.png'),
    )
    plots_generated.append(('sensitivity_heatmap.png',
        'Sobol ST heatmap: 20 params × 14 metrics, category sidebar, top 5 cells highlighted'))

    # Interaction web for pop_crash_pct
    print('  → interaction_web_popcrash.png')
    plot_parameter_interaction_web(
        sobol_indices, metric='pop_crash_pct', top_n=12,
        save_path=str(OUT / 'interaction_web_popcrash.png'),
    )
    plots_generated.append(('interaction_web_popcrash.png',
        'Parameter interaction matrix: approximate pairwise interaction strength for pop crash'))

# 8. Evolutionary rescue assessment
print('  → evolutionary_rescue.png')
plot_evolutionary_rescue_assessment(
    result_pe, carrying_capacity=500, disease_year=3,
    save_path=str(OUT / 'evolutionary_rescue.png'),
)
plots_generated.append(('evolutionary_rescue.png',
    '4-panel rescue assessment: population vs threshold, Δr̄ rate, V_A depletion, generations to recovery'))

# 9. Conservation scenario matrix (using the 3 runs as stand-ins)
print('  → conservation_matrix.png')
conservation_scenarios = {
    'No Intervention': result_nope,
    'Captive Breeding (placeholder)': result_pe,
    'Assisted Gene Flow (placeholder)': result_pe8,
}
plot_conservation_scenario_matrix(
    conservation_scenarios,
    save_path=str(OUT / 'conservation_matrix.png'),
)
plots_generated.append(('conservation_matrix.png',
    '3×3 conservation comparison grid: no intervention vs captive breeding vs gene flow (placeholder data)'))

# 10. Model validation panel
print('  → model_validation.png')
plot_model_validation_panel(
    result_pe,
    observed_data={
        'hamilton_mortality_low': 0.86,
        'hamilton_mortality_high': 0.96,
        'sighting_year': 15,  # ~2025 if disease_year=3 ≈ 2013
    },
    disease_year=3,
    save_path=str(OUT / 'model_validation.png'),
)
plots_generated.append(('model_validation.png',
    'Model vs Hamilton 2021 mortality range, recovery timeline, decline pattern'))

# 11. Morris vs Sobol scatter
if sobol_indices is not None and morris_results is not None:
    print('  → morris_vs_sobol.png')
    plot_parameter_space_exploration(
        morris_results, sobol_indices, metric='pop_crash_pct',
        save_path=str(OUT / 'morris_vs_sobol.png'),
    )
    plots_generated.append(('morris_vs_sobol.png',
        'Morris μ* vs Sobol ST scatter: each point is a parameter, colored by category'))

    # Also do for resistance_shift_mean
    print('  → morris_vs_sobol_resistance.png')
    plot_parameter_space_exploration(
        morris_results, sobol_indices, metric='resistance_shift_mean',
        save_path=str(OUT / 'morris_vs_sobol_resistance.png'),
    )
    plots_generated.append(('morris_vs_sobol_resistance.png',
        'Morris vs Sobol scatter for resistance shift metric'))

# ═══════════════════════════════════════════════════════════════════════
# 7. Write INDEX.md
# ═══════════════════════════════════════════════════════════════════════
print(f'[7/7] Writing INDEX.md ({len(plots_generated)} plots)...')

# Collect ALL viz examples across all modules
index_path = Path('results/viz_examples/INDEX.md')
lines = ['# SSWD-EvoEpi Visualization Examples Index\n']
lines.append(f'Generated: 2026-02-18\n')
lines.append(f'Total modules: 6 (population, disease, genetics, coevolution, spatial, dashboards)\n\n')

# Scan all subdirs
for subdir in sorted(Path('results/viz_examples').iterdir()):
    if subdir.is_dir():
        pngs = sorted(subdir.glob('*.png'))
        if pngs:
            lines.append(f'## {subdir.name.title()} ({len(pngs)} plots)\n\n')
            for png in pngs:
                # Find description from generated list if dashboards
                desc = ''
                if subdir.name == 'dashboards':
                    for fname, d in plots_generated:
                        if fname == png.name:
                            desc = f' — {d}'
                            break
                lines.append(f'- `{png.name}`{desc}\n')
            lines.append('\n')

with open(index_path, 'w') as f:
    f.writelines(lines)

print(f'\n✅ Done! {len(plots_generated)} dashboard plots saved to {OUT}/')
print(f'   INDEX.md written to {index_path}')
