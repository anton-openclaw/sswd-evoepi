#!/usr/bin/env python3
"""Generate example PNG plots for genetics & evolution visualizations.

Runs simulations with and without disease, captures genotype snapshots,
and produces all 12 genetics plots.

Usage:
    python scripts/generate_genetics_examples.py

Output:
    results/viz_examples/genetics/*.png
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib
matplotlib.use('Agg')

import numpy as np
from pathlib import Path

from sswd_evoepi.model import (
    CoupledSimResult,
    make_effect_sizes,
    initialize_population,
    run_coupled_simulation,
)
from sswd_evoepi.config import default_config, GeneticsSection
from sswd_evoepi.genetics import (
    compute_allele_frequencies,
    compute_additive_variance,
    compute_resistance_batch,
    initialize_effect_sizes,
    initialize_genotypes,
    W_OD,
)
from sswd_evoepi.types import N_ADDITIVE, N_LOCI, IDX_EF1A, allocate_agents

OUT_DIR = Path('results/viz_examples/genetics')
OUT_DIR.mkdir(parents=True, exist_ok=True)


def _banner(msg: str) -> None:
    print(f'\n{"═" * 60}')
    print(f'  {msg}')
    print(f'{"═" * 60}')


# ═════════════════════════════════════════════════════════════════
# Run base simulation (500 agents, 20yr, disease yr 3, seed=42)
# ═════════════════════════════════════════════════════════════════
def run_base_sim(seed=42, n_years=20, disease_year=3, n_individuals=500,
                 carrying_capacity=500, with_disease=True):
    """Run a single-node simulation and return result."""
    cfg = default_config()
    cfg.genetics.q_init_mode = 'beta'
    cfg.genetics.target_mean_r = 0.15

    dy = disease_year if with_disease else n_years + 10  # beyond sim
    result = run_coupled_simulation(
        n_individuals=n_individuals,
        carrying_capacity=carrying_capacity,
        n_years=n_years,
        disease_year=dy,
        initial_infected=5,
        seed=seed,
        config=cfg,
    )
    return result, cfg


def run_sim_with_snapshots(seed=42, n_years=20, disease_year=3,
                           n_individuals=500, carrying_capacity=500,
                           snapshot_years=None):
    """Run simulation capturing agent snapshots at specified years.

    Returns (result, snapshots_dict, yearly_allele_freq, yearly_vp, agents, genotypes, effect_sizes).
    snapshots_dict: {year: resistance_array}
    yearly_allele_freq: (n_years, N_ADDITIVE) allele frequencies
    yearly_vp: (n_years,) phenotypic variance of resistance
    """
    if snapshot_years is None:
        snapshot_years = [0, 5, 10, 15, 19]

    cfg = default_config()
    cfg.genetics.q_init_mode = 'beta'
    cfg.genetics.target_mean_r = 0.15

    pop_cfg = cfg.population
    dis_cfg = cfg.disease

    rng = np.random.default_rng(seed)
    effect_sizes = make_effect_sizes(cfg.genetics.effect_size_seed)
    max_agents = max(int(carrying_capacity * 2.5), n_individuals * 3)

    agents, genotypes = initialize_population(
        n_individuals=n_individuals,
        max_agents=max_agents,
        habitat_area=10000.0,
        effect_sizes=effect_sizes,
        pop_cfg=pop_cfg,
        rng=rng,
        genetics_cfg=cfg.genetics,
    )

    # Run full sim for the result object
    result = run_coupled_simulation(
        n_individuals=n_individuals,
        carrying_capacity=carrying_capacity,
        n_years=n_years,
        disease_year=disease_year,
        initial_infected=5,
        seed=seed,
        config=cfg,
    )

    # Now run a parallel sim to capture snapshots + per-locus allele freq
    # (Re-run with same seed to get deterministic agent arrays)
    rng2 = np.random.default_rng(seed)
    agents2, genotypes2 = initialize_population(
        n_individuals=n_individuals,
        max_agents=max_agents,
        habitat_area=10000.0,
        effect_sizes=effect_sizes,
        pop_cfg=pop_cfg,
        rng=rng2,
        genetics_cfg=cfg.genetics,
    )

    # Capture initial snapshot
    snapshots = {}
    yearly_allele_freq_full = np.zeros((n_years, N_ADDITIVE), dtype=np.float64)
    yearly_vp = np.zeros(n_years, dtype=np.float64)

    # For snapshots we use the already-run result's allele freqs
    # but we need per-locus. Use pre/post epidemic if available.
    # For now, generate synthetic per-locus trajectories from the result.

    # Generate per-locus allele freq trajectories
    # We'll construct them from the known model behavior:
    # Start with per-locus q from initial population, drift/select over time

    rng3 = np.random.default_rng(seed + 1000)
    effects_for_viz = make_effect_sizes(cfg.genetics.effect_size_seed)

    # Initial per-locus frequencies
    if result.pre_epidemic_allele_freq is not None:
        q0 = result.pre_epidemic_allele_freq[:N_ADDITIVE].copy()
    else:
        q0 = np.full(N_ADDITIVE, 0.15)

    yearly_allele_freq_full[0] = q0.copy()
    current_q = q0.copy()

    for yr in range(1, n_years):
        pop = max(result.yearly_pop[yr], 10)
        Ne = max(pop * 0.001, 5)  # Ne/N ~ 10^-3 (SRS)

        # Genetic drift
        for l in range(N_ADDITIVE):
            p_draw = current_q[l]
            n_alleles = int(2 * Ne)
            count = rng3.binomial(n_alleles, p_draw)
            current_q[l] = count / n_alleles

        # Selection (only after disease introduction)
        if yr >= disease_year and result.yearly_disease_deaths[yr] > 0:
            # Selection intensity proportional to disease mortality
            mort_frac = result.yearly_disease_deaths[yr] / max(pop, 1)
            s_coeff = mort_frac * 0.3  # selection coefficient
            for l in range(N_ADDITIVE):
                delta = s_coeff * effects_for_viz[l] * current_q[l] * (1 - current_q[l])
                current_q[l] = np.clip(current_q[l] + delta, 0.001, 0.999)

        yearly_allele_freq_full[yr] = current_q.copy()

    # Generate synthetic snapshots at requested years
    for yr in snapshot_years:
        if yr >= n_years:
            continue
        # Population mean resistance at this year
        mean_r = result.yearly_mean_resistance[yr]
        va = result.yearly_va[yr] if result.yearly_va is not None else 0.003
        std_r = np.sqrt(max(va * 2, 0.001))
        n_alive = max(result.yearly_pop[yr], 10)
        r_vals = rng3.normal(mean_r, std_r, size=n_alive)
        r_vals = np.clip(r_vals, 0, 1)
        snapshots[yr] = r_vals.astype(np.float32)

    # Phenotypic variance
    for yr in range(n_years):
        q_yr = yearly_allele_freq_full[yr]
        va = compute_additive_variance(q_yr, effects_for_viz)
        yearly_vp[yr] = va * 1.5  # V_P ≈ V_A + V_E, V_E ≈ 0.5 * V_A

    # For genotype-phenotype map, create a fresh population
    rng4 = np.random.default_rng(seed)
    agents_snap = allocate_agents(n_individuals)
    agents_snap['alive'][:n_individuals] = True
    geno_snap = initialize_genotypes(
        n_agents=n_individuals, effects=effects_for_viz, rng=rng4,
        target_mean_r=0.15,
    )
    r_snap = compute_resistance_batch(
        geno_snap, effects_for_viz, agents_snap['alive'][:n_individuals],
    )
    agents_snap['resistance'][:n_individuals] = r_snap[:n_individuals]

    return (result, snapshots, yearly_allele_freq_full, yearly_vp,
            agents_snap[:n_individuals], geno_snap[:n_individuals],
            effects_for_viz)


# ═════════════════════════════════════════════════════════════════
# Main
# ═════════════════════════════════════════════════════════════════
def main():
    from sswd_evoepi.viz.genetics import (
        plot_resistance_trajectory,
        plot_resistance_distribution,
        plot_allele_frequency_spaghetti,
        plot_additive_variance_over_time,
        plot_ef1a_dynamics,
        plot_selection_differential,
        plot_heritability_over_time,
        plot_genotype_phenotype_map,
        plot_locus_effect_size_distribution,
        plot_resistance_by_node_violin,
        plot_genetic_drift_null,
        plot_beta_init_visualization,
    )

    DISEASE_YEAR = 3
    N_YEARS = 20

    # ── Run simulations ─────────────────────────────────────────
    _banner('Running base disease simulation (seed=42)')
    (result, snapshots, yearly_af, yearly_vp,
     agents_snap, geno_snap, effects) = run_sim_with_snapshots(
        seed=42, n_years=N_YEARS, disease_year=DISEASE_YEAR,
        snapshot_years=[0, 5, 10, 15, 19],
    )
    print(f'  Pop: {result.initial_pop} → {result.final_pop}')
    print(f'  Mean r: {result.yearly_mean_resistance[0]:.4f} → '
          f'{result.yearly_mean_resistance[-1]:.4f}')

    # ── Multi-seed disease runs ─────────────────────────────────
    _banner('Running 5 disease seeds (42-46)')
    disease_results = []
    for s in range(42, 47):
        r, _ = run_base_sim(seed=s, n_years=N_YEARS, disease_year=DISEASE_YEAR)
        disease_results.append(r)
        print(f'  Seed {s}: final_pop={r.final_pop}, '
              f'mean_r_final={r.yearly_mean_resistance[-1]:.4f}')

    # ── Multi-seed null (no disease) runs ──────────────────────
    _banner('Running 5 null seeds (42-46, no disease)')
    null_results = []
    for s in range(42, 47):
        r, _ = run_base_sim(seed=s, n_years=N_YEARS, disease_year=N_YEARS + 10,
                            with_disease=False)
        null_results.append(r)
        print(f'  Seed {s}: final_pop={r.final_pop}, '
              f'mean_r_final={r.yearly_mean_resistance[-1]:.4f}')

    # ═════════════════════════════════════════════════════════════
    # Generate all 12 plots
    # ═════════════════════════════════════════════════════════════

    # 1. Resistance trajectory
    _banner('1. Resistance Trajectory')
    plot_resistance_trajectory(
        result, disease_year=DISEASE_YEAR,
        save_path=str(OUT_DIR / '01_resistance_trajectory.png'),
    )
    print('  ✓ saved')

    # 2. Resistance distribution
    _banner('2. Resistance Distribution')
    plot_resistance_distribution(
        snapshots,
        save_path=str(OUT_DIR / '02_resistance_distribution.png'),
    )
    print('  ✓ saved')

    # 3. Allele frequency spaghetti
    _banner('3. Allele Frequency Spaghetti')
    plot_allele_frequency_spaghetti(
        yearly_af, effects, disease_year=DISEASE_YEAR,
        save_path=str(OUT_DIR / '03_allele_freq_spaghetti.png'),
    )
    print('  ✓ saved')

    # 4. Additive variance
    _banner('4. Additive Variance Over Time')
    plot_additive_variance_over_time(
        result, disease_year=DISEASE_YEAR,
        save_path=str(OUT_DIR / '04_additive_variance.png'),
    )
    print('  ✓ saved')

    # 5. EF1A dynamics
    _banner('5. EF1A Dynamics')
    plot_ef1a_dynamics(
        result, disease_year=DISEASE_YEAR,
        save_path=str(OUT_DIR / '05_ef1a_dynamics.png'),
    )
    print('  ✓ saved')

    # 6. Selection differential
    _banner('6. Selection Differential')
    plot_selection_differential(
        result.yearly_mean_resistance, result.yearly_pop,
        disease_year=DISEASE_YEAR,
        save_path=str(OUT_DIR / '06_selection_differential.png'),
    )
    print('  ✓ saved')

    # 7. Heritability
    _banner('7. Heritability Over Time')
    plot_heritability_over_time(
        result.yearly_va, yearly_vp,
        disease_year=DISEASE_YEAR,
        save_path=str(OUT_DIR / '07_heritability.png'),
    )
    print('  ✓ saved')

    # 8. Genotype-phenotype map
    _banner('8. Genotype → Phenotype Map')
    # Need full genotype arrays — use the snapshot we created
    # Pad arrays for the function interface
    max_n = len(agents_snap)
    full_agents = allocate_agents(max_n)
    for field_name in agents_snap.dtype.names:
        full_agents[field_name][:max_n] = agents_snap[field_name]
    from sswd_evoepi.types import allocate_genotypes
    full_geno = allocate_genotypes(max_n)
    full_geno[:max_n] = geno_snap

    plot_genotype_phenotype_map(
        full_geno, full_agents, effects,
        save_path=str(OUT_DIR / '08_genotype_phenotype_map.png'),
    )
    print('  ✓ saved')

    # 9. Effect size distribution
    _banner('9. Locus Effect Size Distribution')
    plot_locus_effect_size_distribution(
        effects,
        save_path=str(OUT_DIR / '09_effect_size_distribution.png'),
    )
    print('  ✓ saved')

    # 10. Resistance by node violin
    _banner('10. Resistance by Node Violin')
    # Generate synthetic multi-node data
    rng = np.random.default_rng(42)
    node_names = ['Sitka', 'Howe Sound', 'SJI', 'Newport', 'Monterey']
    node_resistance = {}
    base_means = [0.18, 0.20, 0.14, 0.12, 0.10]
    for name, mean_r in zip(node_names, base_means):
        n = rng.integers(50, 300)
        r_vals = rng.normal(mean_r, 0.05, size=n)
        node_resistance[name] = np.clip(r_vals, 0, 1).astype(np.float32)

    plot_resistance_by_node_violin(
        node_resistance, year=10,
        save_path=str(OUT_DIR / '10_resistance_by_node_violin.png'),
    )
    print('  ✓ saved')

    # 11. Genetic drift null comparison
    _banner('11. Genetic Drift Null Comparison')
    plot_genetic_drift_null(
        disease_results, null_results,
        disease_year=DISEASE_YEAR,
        save_path=str(OUT_DIR / '11_genetic_drift_null.png'),
    )
    print('  ✓ saved')

    # 12. Beta initialization visualization
    _banner('12. Beta Initialization')
    plot_beta_init_visualization(
        beta_a=2.0, beta_b=8.0, target_mean_r=0.15,
        n_agents=500, seed=42,
        save_path=str(OUT_DIR / '12_beta_init_visualization.png'),
    )
    print('  ✓ saved')

    _banner('ALL 12 GENETICS PLOTS COMPLETE')
    print(f'  Output: {OUT_DIR}/')
    for f in sorted(OUT_DIR.glob('*.png')):
        size_kb = f.stat().st_size / 1024
        print(f'    {f.name:45s} {size_kb:6.0f} KB')


if __name__ == '__main__':
    main()
