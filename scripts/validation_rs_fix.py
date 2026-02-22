#!/usr/bin/env python3
"""Validation with R→S fix: recovered stars return to susceptible.

Compares against baseline (permanent immunity) validation at K=5000.
Runs two simulations:
  1. R→S with sinusoidal SST (direct comparison to baseline)
  2. R→S with satellite SST (if data available)

New metrics enabled by R→S:
  - Total recovery events (can exceed unique infected count)
  - Reinfection tracking (inferred from recovery count vs population)
"""

import sys, os, time, json
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from sswd_evoepi.config import default_config
from sswd_evoepi.model import run_spatial_simulation
from sswd_evoepi.spatial import NodeDefinition, build_network
from sswd_evoepi.types import DiseaseState

# ═══════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════

K = 5000
YEARS = 20
DISEASE_YEAR = 3
SEED = 42

OUTDIR = "results/validation_rs_fix"
os.makedirs(OUTDIR, exist_ok=True)

# 5-node network (same as baseline)
NODE_DEFS = [
    NodeDefinition(node_id=0, name="Sitka", lat=57.06, lon=-135.34,
        subregion="AK-SE", habitat_area=333333.0, carrying_capacity=K,
        is_fjord=False, flushing_rate=0.8, mean_sst=8.0, sst_amplitude=3.5,
        sst_trend=0.015, salinity=32.0),
    NodeDefinition(node_id=1, name="Howe Sound", lat=49.52, lon=-123.25,
        subregion="SS", habitat_area=333333.0, carrying_capacity=K,
        is_fjord=True, sill_depth=30.0, flushing_rate=0.03, mean_sst=10.0,
        sst_amplitude=4.0, sst_trend=0.02, salinity=22.0),
    NodeDefinition(node_id=2, name="SJI", lat=48.53, lon=-123.02,
        subregion="SS", habitat_area=333333.0, carrying_capacity=K,
        is_fjord=False, flushing_rate=0.3, mean_sst=10.5, sst_amplitude=4.0,
        sst_trend=0.02, salinity=30.0),
    NodeDefinition(node_id=3, name="Newport", lat=44.63, lon=-124.05,
        subregion="OR", habitat_area=333333.0, carrying_capacity=K,
        is_fjord=False, flushing_rate=0.6, mean_sst=11.5, sst_amplitude=3.0,
        sst_trend=0.02, salinity=33.0),
    NodeDefinition(node_id=4, name="Monterey", lat=36.62, lon=-121.90,
        subregion="CA-CEN", habitat_area=333333.0, carrying_capacity=K,
        is_fjord=False, flushing_rate=0.4, mean_sst=13.0, sst_amplitude=2.5,
        sst_trend=0.025, salinity=33.5),
]

NAMES = [n.name for n in NODE_DEFS]


def build_net(cfg):
    """Build network from node definitions."""
    return build_network(
        NODE_DEFS,
        D_L=cfg.spatial.D_L,
        D_P=cfg.spatial.D_P,
        r_total=cfg.spatial.r_total,
        alpha_self_fjord=cfg.spatial.alpha_self_fjord,
        alpha_self_open=cfg.spatial.alpha_self_open,
    )


def run_and_extract(cfg, label, net):
    """Run simulation and extract per-node results."""
    print(f"\n{'='*70}")
    print(f"  {label}")
    print(f"  K={K:,}, {YEARS}yr, disease yr {DISEASE_YEAR}, seed={SEED}")
    print(f"  SST source: {cfg.simulation.sst_source}")
    print(f"  Config: rho_rec={cfg.disease.rho_rec}, "
          f"target_mean_r={cfg.genetics.target_mean_r}, "
          f"target_mean_t={cfg.genetics.target_mean_t}, "
          f"target_mean_c={cfg.genetics.target_mean_c}, "
          f"tau_max={cfg.disease.tau_max}")
    print(f"{'='*70}\n")

    t0 = time.time()
    res = run_spatial_simulation(net, n_years=YEARS, disease_year=DISEASE_YEAR,
                                  seed=SEED, config=cfg)
    elapsed = time.time() - t0
    print(f"\nRuntime: {elapsed:.1f}s ({elapsed/60:.1f} min)")

    # ── Extract per-node results ──────────────────────────────────────
    results = {}
    print(f"\n{'Node':<14} {'Init':>7} {'Min':>7} {'Final':>7} {'Crash%':>7} "
          f"{'Deaths':>7} {'Recov':>6} {'Δr':>8} {'Δt':>8} {'Δc':>8}")
    print("-" * 100)

    for i in range(5):
        pop_init = int(res.yearly_pop[i, 0])
        pop_final = int(res.yearly_pop[i, -1])
        pop_min = int(res.yearly_pop[i, DISEASE_YEAR:].min())
        min_yr = int(np.argmin(res.yearly_pop[i, DISEASE_YEAR:]) + DISEASE_YEAR)
        crash = 100 * (1 - pop_min / pop_init) if pop_init > 0 else 0

        deaths = int(res.yearly_disease_deaths[i].sum())
        recov = int(res.yearly_recoveries[i].sum()) if res.yearly_recoveries is not None else 0

        # Trait shifts (pre-disease baseline at year 2 vs final)
        dr = float(res.yearly_mean_resistance[i, -1] - res.yearly_mean_resistance[i, 2])
        dt = float(res.yearly_mean_tolerance[i, -1] - res.yearly_mean_tolerance[i, 2])
        dc = float(res.yearly_mean_recovery[i, -1] - res.yearly_mean_recovery[i, 2])

        results[NAMES[i]] = {
            'pop_init': pop_init, 'pop_final': pop_final, 'pop_min': pop_min,
            'min_yr': min_yr, 'crash_pct': round(crash, 1),
            'deaths': deaths, 'recoveries': recov,
            'delta_r': round(dr, 4), 'delta_t': round(dt, 4), 'delta_c': round(dc, 4),
            'r_final': round(float(res.yearly_mean_resistance[i, -1]), 4),
            't_final': round(float(res.yearly_mean_tolerance[i, -1]), 4),
            'c_final': round(float(res.yearly_mean_recovery[i, -1]), 4),
            'yearly_pop': res.yearly_pop[i].tolist(),
            'yearly_deaths': res.yearly_disease_deaths[i].tolist(),
            'yearly_recov': res.yearly_recoveries[i].tolist() if res.yearly_recoveries is not None else [],
            'yearly_mean_r': [round(float(x), 4) for x in res.yearly_mean_resistance[i]],
            'yearly_mean_t': [round(float(x), 4) for x in res.yearly_mean_tolerance[i]],
            'yearly_mean_c': [round(float(x), 4) for x in res.yearly_mean_recovery[i]],
        }

        print(f"{NAMES[i]:<14} {pop_init:>7,} {pop_min:>7,} {pop_final:>7,} "
              f"{crash:>6.1f}% {deaths:>7,} {recov:>6} "
              f"{dr:>+7.4f} {dt:>+7.4f} {dc:>+7.4f}")

    print()
    total_init = sum(r['pop_init'] for r in results.values())
    total_final = sum(r['pop_final'] for r in results.values())
    total_min = sum(r['pop_min'] for r in results.values())
    total_deaths = sum(r['deaths'] for r in results.values())
    total_recov = sum(r['recoveries'] for r in results.values())
    overall_crash = 100 * (1 - total_min / total_init) if total_init > 0 else 0
    recov_rate = total_recov / total_deaths * 100 if total_deaths > 0 else 0

    print(f"Total: {total_init:,} → {total_final:,}")
    print(f"Overall crash: {overall_crash:.1f}%")
    print(f"Disease deaths: {total_deaths:,}, Recoveries: {total_recov:,} "
          f"({recov_rate:.2f}% of infections)")

    return {
        'config': {
            'K': K, 'years': YEARS, 'disease_year': DISEASE_YEAR, 'seed': SEED,
            'sst_source': cfg.simulation.sst_source,
            'rho_rec': cfg.disease.rho_rec, 'tau_max': cfg.disease.tau_max,
            'target_mean_r': cfg.genetics.target_mean_r,
            'target_mean_t': cfg.genetics.target_mean_t,
            'target_mean_c': cfg.genetics.target_mean_c,
            'partition': f'{cfg.genetics.n_resistance}/{cfg.genetics.n_tolerance}/{cfg.genetics.n_recovery}',
        },
        'runtime_s': round(elapsed, 1),
        'per_node': results,
        'totals': {
            'initial': total_init, 'final': total_final,
            'overall_crash_pct': round(overall_crash, 1),
            'deaths': total_deaths, 'recoveries': total_recov,
            'recovery_rate_pct': round(recov_rate, 4),
        },
    }, res


def generate_figures(res_sin, res_sat, results_sin, results_sat, outdir):
    """Generate comparison figures."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
        plot_years = np.arange(YEARS)

        # ── Figure 1: Population trajectories (sinusoidal) ────────────
        fig, ax = plt.subplots(figsize=(12, 6))
        for i in range(5):
            ax.plot(plot_years, res_sin.yearly_pop[i], color=colors[i],
                    linewidth=2, label=NAMES[i])
        ax.axvline(DISEASE_YEAR, color='red', linestyle='--', alpha=0.5,
                   label='Epidemic onset')
        ax.set_xlabel('Year', fontsize=12)
        ax.set_ylabel('Population', fontsize=12)
        ax.set_title(f'R→S Population Trajectories (sinusoidal SST, K={K:,})', fontsize=14)
        ax.legend(fontsize=10)
        ax.set_ylim(bottom=0)
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(f'{outdir}/pop_trajectories_sinusoidal.png', dpi=150)
        plt.close(fig)
        print(f"Saved: {outdir}/pop_trajectories_sinusoidal.png")

        # ── Figure 2: Three-trait evolution (sinusoidal) ──────────────
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        trait_data = [
            (res_sin.yearly_mean_resistance, 'Resistance (r_i)', 'Mean resistance'),
            (res_sin.yearly_mean_tolerance, 'Tolerance (t_i)', 'Mean tolerance'),
            (res_sin.yearly_mean_recovery, 'Recovery (c_i)', 'Mean recovery'),
        ]
        for ax, (data, title, ylabel) in zip(axes, trait_data):
            for i in range(5):
                ax.plot(plot_years, data[i], color=colors[i], linewidth=2,
                        label=NAMES[i])
            ax.axvline(DISEASE_YEAR, color='red', linestyle='--', alpha=0.5)
            ax.set_xlabel('Year', fontsize=11)
            ax.set_ylabel(ylabel, fontsize=11)
            ax.set_title(title, fontsize=13)
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)
        fig.suptitle(f'R→S Trait Evolution (sinusoidal SST, K={K:,})', fontsize=14, y=1.02)
        fig.tight_layout()
        fig.savefig(f'{outdir}/trait_evolution_sinusoidal.png', dpi=150,
                    bbox_inches='tight')
        plt.close(fig)
        print(f"Saved: {outdir}/trait_evolution_sinusoidal.png")

        # ── Figure 3: Recovery events per year ────────────────────────
        fig, ax = plt.subplots(figsize=(12, 6))
        if res_sin.yearly_recoveries is not None:
            for i in range(5):
                ax.plot(plot_years, res_sin.yearly_recoveries[i], color=colors[i],
                        linewidth=2, label=NAMES[i])
        ax.axvline(DISEASE_YEAR, color='red', linestyle='--', alpha=0.5,
                   label='Epidemic onset')
        ax.set_xlabel('Year', fontsize=12)
        ax.set_ylabel('Recovery events', fontsize=12)
        ax.set_title(f'R→S Recovery Events per Year (K={K:,})', fontsize=14)
        ax.legend(fontsize=10)
        ax.set_ylim(bottom=0)
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(f'{outdir}/yearly_recoveries.png', dpi=150)
        plt.close(fig)
        print(f"Saved: {outdir}/yearly_recoveries.png")

        # ── Figure 4: Trait shift comparison (R→S sinusoidal) ─────────
        fig, ax = plt.subplots(figsize=(10, 6))
        x = np.arange(5)
        w = 0.25
        dr_vals = [results_sin['per_node'][n]['delta_r'] for n in NAMES]
        dt_vals = [results_sin['per_node'][n]['delta_t'] for n in NAMES]
        dc_vals = [results_sin['per_node'][n]['delta_c'] for n in NAMES]
        ax.bar(x - w, dr_vals, w, label='Δ Resistance', color='#1f77b4')
        ax.bar(x, dt_vals, w, label='Δ Tolerance', color='#ff7f0e')
        ax.bar(x + w, dc_vals, w, label='Δ Recovery', color='#2ca02c')
        ax.set_xticks(x)
        ax.set_xticklabels(NAMES, fontsize=10)
        ax.set_ylabel('Trait shift', fontsize=11)
        ax.set_title(f'R→S Evolutionary Response by Trait (sinusoidal SST, K={K:,})', fontsize=13)
        ax.legend(fontsize=10)
        ax.axhline(0, color='black', linewidth=0.5)
        ax.grid(True, alpha=0.3, axis='y')
        fig.tight_layout()
        fig.savefig(f'{outdir}/trait_shifts_sinusoidal.png', dpi=150)
        plt.close(fig)
        print(f"Saved: {outdir}/trait_shifts_sinusoidal.png")

        # ── Figure 5: Sinusoidal vs Satellite SST comparison ──────────
        if res_sat is not None:
            fig, axes = plt.subplots(1, 2, figsize=(16, 6))

            # Population comparison
            for i in range(5):
                axes[0].plot(plot_years, res_sin.yearly_pop[i], color=colors[i],
                            linewidth=2, linestyle='-', label=f'{NAMES[i]} (sin)')
                axes[0].plot(plot_years, res_sat.yearly_pop[i], color=colors[i],
                            linewidth=2, linestyle='--', label=f'{NAMES[i]} (sat)')
            axes[0].axvline(DISEASE_YEAR, color='red', linestyle='--', alpha=0.5)
            axes[0].set_xlabel('Year', fontsize=11)
            axes[0].set_ylabel('Population', fontsize=11)
            axes[0].set_title('Population: Sinusoidal vs Satellite SST', fontsize=13)
            axes[0].legend(fontsize=7, ncol=2)
            axes[0].set_ylim(bottom=0)
            axes[0].grid(True, alpha=0.3)

            # Recovery comparison
            if res_sin.yearly_recoveries is not None and res_sat.yearly_recoveries is not None:
                for i in range(5):
                    axes[1].plot(plot_years, res_sin.yearly_recoveries[i], color=colors[i],
                                linewidth=2, linestyle='-', label=f'{NAMES[i]} (sin)')
                    axes[1].plot(plot_years, res_sat.yearly_recoveries[i], color=colors[i],
                                linewidth=2, linestyle='--', label=f'{NAMES[i]} (sat)')
            axes[1].axvline(DISEASE_YEAR, color='red', linestyle='--', alpha=0.5)
            axes[1].set_xlabel('Year', fontsize=11)
            axes[1].set_ylabel('Recovery events', fontsize=11)
            axes[1].set_title('Recoveries: Sinusoidal vs Satellite SST', fontsize=13)
            axes[1].legend(fontsize=7, ncol=2)
            axes[1].set_ylim(bottom=0)
            axes[1].grid(True, alpha=0.3)

            fig.tight_layout()
            fig.savefig(f'{outdir}/sinusoidal_vs_satellite.png', dpi=150)
            plt.close(fig)
            print(f"Saved: {outdir}/sinusoidal_vs_satellite.png")

            # Satellite pop trajectories
            fig, ax = plt.subplots(figsize=(12, 6))
            for i in range(5):
                ax.plot(plot_years, res_sat.yearly_pop[i], color=colors[i],
                        linewidth=2, label=NAMES[i])
            ax.axvline(DISEASE_YEAR, color='red', linestyle='--', alpha=0.5,
                       label='Epidemic onset')
            ax.set_xlabel('Year', fontsize=12)
            ax.set_ylabel('Population', fontsize=12)
            ax.set_title(f'R→S Population Trajectories (satellite SST, K={K:,})', fontsize=14)
            ax.legend(fontsize=10)
            ax.set_ylim(bottom=0)
            ax.grid(True, alpha=0.3)
            fig.tight_layout()
            fig.savefig(f'{outdir}/pop_trajectories_satellite.png', dpi=150)
            plt.close(fig)
            print(f"Saved: {outdir}/pop_trajectories_satellite.png")

        print("\nAll figures saved.")

    except Exception as e:
        import traceback
        print(f"\nFigure generation error: {e}")
        traceback.print_exc()


# ═══════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    print("=" * 70)
    print("  R→S VALIDATION: Recovered stars return to susceptible")
    print("  Echinoderms lack adaptive immunity — reinfection possible")
    print("=" * 70)

    # ── Run 1: Sinusoidal SST (direct comparison to baseline) ─────────
    cfg_sin = default_config()
    cfg_sin.simulation.sst_source = 'sinusoidal'
    net_sin = build_net(cfg_sin)
    results_sin, res_sin = run_and_extract(cfg_sin, "RUN 1: R→S + Sinusoidal SST", net_sin)

    # ── Run 2: Satellite SST (if data available) ─────────────────────
    results_sat = None
    res_sat = None
    sst_data_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'sst')
    has_sst = all(
        os.path.exists(os.path.join(sst_data_dir, f"{n.replace(' ', '_')}_climatology.csv"))
        for n in NAMES
    )

    if has_sst:
        cfg_sat = default_config()
        cfg_sat.simulation.sst_source = 'satellite'
        cfg_sat.simulation.sst_data_dir = sst_data_dir
        net_sat = build_net(cfg_sat)
        results_sat, res_sat = run_and_extract(cfg_sat, "RUN 2: R→S + Satellite SST", net_sat)
    else:
        print("\n⚠ Satellite SST data not found for all nodes — skipping Run 2")

    # ── Generate figures ──────────────────────────────────────────────
    print("\n\n" + "=" * 70)
    print("  GENERATING FIGURES")
    print("=" * 70)
    generate_figures(res_sin, res_sat, results_sin, results_sat, OUTDIR)

    # ── Save JSON results ─────────────────────────────────────────────
    output = {
        'sinusoidal': results_sin,
    }
    if results_sat:
        output['satellite'] = results_sat

    with open(f'{OUTDIR}/results.json', 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved: {OUTDIR}/results.json")

    # ── Print summary ─────────────────────────────────────────────────
    print("\n\n" + "=" * 70)
    print("  SUMMARY")
    print("=" * 70)
    print(f"\nSinusoidal SST:")
    print(f"  Overall crash: {results_sin['totals']['overall_crash_pct']}%")
    print(f"  Total deaths:  {results_sin['totals']['deaths']:,}")
    print(f"  Total recov:   {results_sin['totals']['recoveries']:,} "
          f"({results_sin['totals']['recovery_rate_pct']:.2f}%)")
    if results_sat:
        print(f"\nSatellite SST:")
        print(f"  Overall crash: {results_sat['totals']['overall_crash_pct']}%")
        print(f"  Total deaths:  {results_sat['totals']['deaths']:,}")
        print(f"  Total recov:   {results_sat['totals']['recoveries']:,} "
              f"({results_sat['totals']['recovery_rate_pct']:.2f}%)")

    print("\nDone.")
