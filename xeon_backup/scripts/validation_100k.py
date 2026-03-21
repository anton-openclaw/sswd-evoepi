#!/usr/bin/env python3
"""5-node validation with K=100,000 per node.

Large-scale run to test scaling + recovery dynamics at realistic population sizes.
"""

import sys, os, time, json
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from sswd_evoepi.config import default_config
from sswd_evoepi.model import run_spatial_simulation
from sswd_evoepi.spatial import NodeDefinition, build_network

cfg = default_config()

K = 100000
YEARS = 20
DISEASE_YEAR = 3
SEED = 42

print(f"=== 100K Validation: K={K:,}, {YEARS}yr, disease yr {DISEASE_YEAR}, seed={SEED} ===")
print(f"Config: rho_rec={cfg.disease.rho_rec}, target_mean_r={cfg.genetics.target_mean_r}, "
      f"target_mean_t={cfg.genetics.target_mean_t}, target_mean_c={cfg.genetics.target_mean_c}, "
      f"tau_max={cfg.disease.tau_max}")
print(f"Partition: {cfg.genetics.n_resistance}R/{cfg.genetics.n_tolerance}T/{cfg.genetics.n_recovery}C")
print()

node_defs = [
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

names = [n.name for n in node_defs]
net = build_network(
    node_defs,
    D_L=cfg.spatial.D_L,
    D_P=cfg.spatial.D_P,
    r_total=cfg.spatial.r_total,
    alpha_self_fjord=cfg.spatial.alpha_self_fjord,
    alpha_self_open=cfg.spatial.alpha_self_open,
)

print(f"Network built: {len(node_defs)} nodes, K={K:,} each")
print(f"Total initial agents: ~{K*5:,}")
print()

t0 = time.time()
res = run_spatial_simulation(net, n_years=YEARS, disease_year=DISEASE_YEAR, seed=SEED, config=cfg)
elapsed = time.time() - t0

print(f"\nRuntime: {elapsed:.1f}s ({elapsed/60:.1f} min)")
print()

# ─── Results table ────────────────────────────────────────────────────

print(f"{'Node':<14} {'Init':>9} {'Min':>9} {'Final':>9} {'Crash%':>7} {'Deaths':>9} {'Recov':>7} {'Δr':>8} {'Δt':>8} {'Δc':>8}")
print("-" * 110)

results = {}
for i in range(5):
    pop_init = int(res.yearly_pop[i, 0])
    pop_final = int(res.yearly_pop[i, -1])
    pop_min = int(res.yearly_pop[i, DISEASE_YEAR:].min())
    min_yr = int(np.argmin(res.yearly_pop[i, DISEASE_YEAR:]) + DISEASE_YEAR)
    crash = 100 * (1 - pop_min / pop_init) if pop_init > 0 else 0

    deaths = int(res.yearly_disease_deaths[i].sum())
    recov = int(res.yearly_recoveries[i].sum()) if res.yearly_recoveries is not None else -1

    # Trait shifts (pre-disease baseline at year 2 vs final)
    dr = float(res.yearly_mean_resistance[i, -1] - res.yearly_mean_resistance[i, 2])
    dt = float(res.yearly_mean_tolerance[i, -1] - res.yearly_mean_tolerance[i, 2])
    dc = float(res.yearly_mean_recovery[i, -1] - res.yearly_mean_recovery[i, 2])

    results[names[i]] = {
        'pop_init': pop_init, 'pop_final': pop_final, 'pop_min': pop_min,
        'min_yr': min_yr, 'crash_pct': round(crash, 1),
        'deaths': deaths, 'recoveries': recov,
        'delta_r': round(dr, 4), 'delta_t': round(dt, 4), 'delta_c': round(dc, 4),
        'r_final': round(float(res.yearly_mean_resistance[i, -1]), 4),
        't_final': round(float(res.yearly_mean_tolerance[i, -1]), 4),
        'c_final': round(float(res.yearly_mean_recovery[i, -1]), 4),
    }

    print(f"{names[i]:<14} {pop_init:>9,} {pop_min:>9,} {pop_final:>9,} {crash:>6.1f}% {deaths:>9,} {recov:>7,} {dr:>+7.4f} {dt:>+7.4f} {dc:>+7.4f}")

print()
total_init = sum(r['pop_init'] for r in results.values())
total_final = sum(r['pop_final'] for r in results.values())
total_deaths = sum(r['deaths'] for r in results.values())
total_recov = sum(r['recoveries'] for r in results.values())
recov_rate = total_recov / total_deaths * 100 if total_deaths > 0 else 0
print(f"Total: {total_init:,} → {total_final:,} ({100*(1-total_final/total_init):.1f}% decline)")
print(f"Disease deaths: {total_deaths:,}, Recoveries: {total_recov:,} ({recov_rate:.2f}% of infections)")
print(f"Runtime: {elapsed:.1f}s ({elapsed/60:.1f} min)")

# ─── Memory usage ─────────────────────────────────────────────────────
import resource
max_rss_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
print(f"Peak memory: {max_rss_mb:.0f} MB")

# ─── Generate figures ─────────────────────────────────────────────────

outdir = "results/validation_100k"
os.makedirs(outdir, exist_ok=True)

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    plot_years = np.arange(YEARS)

    # Figure 1: Population trajectories
    fig, ax = plt.subplots(figsize=(12, 6))
    for i in range(5):
        ax.plot(plot_years, res.yearly_pop[i], color=colors[i], linewidth=2, label=names[i])
    ax.axvline(DISEASE_YEAR, color='red', linestyle='--', alpha=0.5, label='Epidemic onset')
    ax.set_xlabel('Year', fontsize=12)
    ax.set_ylabel('Population', fontsize=12)
    ax.set_title(f'Population Trajectories (K={K:,}, ρ_rec=0.05, mean_c=0.02)', fontsize=14)
    ax.legend(fontsize=10)
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(f'{outdir}/population_trajectories.png', dpi=150)
    plt.close(fig)
    print(f"\nSaved: {outdir}/population_trajectories.png")

    # Figure 2: Three-trait evolution (3 panels)
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    trait_data = [
        (res.yearly_mean_resistance, 'Resistance (r_i)', 'Mean resistance'),
        (res.yearly_mean_tolerance, 'Tolerance (t_i)', 'Mean tolerance'),
        (res.yearly_mean_recovery, 'Recovery (c_i)', 'Mean recovery'),
    ]
    for ax, (data, title, ylabel) in zip(axes, trait_data):
        for i in range(5):
            ax.plot(plot_years, data[i], color=colors[i], linewidth=2, label=names[i])
        ax.axvline(DISEASE_YEAR, color='red', linestyle='--', alpha=0.5)
        ax.set_xlabel('Year', fontsize=11)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_title(title, fontsize=13)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
    fig.suptitle(f'Three-Trait Evolution (K={K:,})', fontsize=14, y=1.02)
    fig.tight_layout()
    fig.savefig(f'{outdir}/trait_evolution.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {outdir}/trait_evolution.png")

    # Figure 3: Recovery events per year per node
    fig, ax = plt.subplots(figsize=(12, 6))
    if res.yearly_recoveries is not None:
        for i in range(5):
            ax.plot(plot_years, res.yearly_recoveries[i], color=colors[i], linewidth=2, label=names[i])
    ax.axvline(DISEASE_YEAR, color='red', linestyle='--', alpha=0.5, label='Epidemic onset')
    ax.set_xlabel('Year', fontsize=12)
    ax.set_ylabel('Recovery events', fontsize=12)
    ax.set_title(f'Disease Recoveries per Year (K={K:,})', fontsize=14)
    ax.legend(fontsize=10)
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(f'{outdir}/yearly_recoveries.png', dpi=150)
    plt.close(fig)
    print(f"Saved: {outdir}/yearly_recoveries.png")

    # Figure 4: Trait shift comparison
    fig, ax = plt.subplots(figsize=(10, 6))
    x = np.arange(5)
    w = 0.25
    dr_vals = [results[n]['delta_r'] for n in names]
    dt_vals = [results[n]['delta_t'] for n in names]
    dc_vals = [results[n]['delta_c'] for n in names]
    ax.bar(x - w, dr_vals, w, label='Δ Resistance', color='#1f77b4')
    ax.bar(x, dt_vals, w, label='Δ Tolerance', color='#ff7f0e')
    ax.bar(x + w, dc_vals, w, label='Δ Recovery', color='#2ca02c')
    ax.set_xticks(x)
    ax.set_xticklabels(names, fontsize=10)
    ax.set_ylabel('Trait shift', fontsize=11)
    ax.set_title(f'Evolutionary Response by Trait (K={K:,})', fontsize=13)
    ax.legend(fontsize=10)
    ax.axhline(0, color='black', linewidth=0.5)
    ax.grid(True, alpha=0.3, axis='y')
    fig.tight_layout()
    fig.savefig(f'{outdir}/trait_shifts.png', dpi=150)
    plt.close(fig)
    print(f"Saved: {outdir}/trait_shifts.png")

except Exception as e:
    import traceback
    print(f"Figure generation error: {e}")
    traceback.print_exc()

# Save JSON results
with open(f'{outdir}/results.json', 'w') as f:
    json.dump({
        'config': {
            'K': K, 'years': YEARS, 'disease_year': DISEASE_YEAR, 'seed': SEED,
            'rho_rec': cfg.disease.rho_rec, 'tau_max': cfg.disease.tau_max,
            'target_mean_r': cfg.genetics.target_mean_r,
            'target_mean_t': cfg.genetics.target_mean_t,
            'target_mean_c': cfg.genetics.target_mean_c,
            'partition': f'{cfg.genetics.n_resistance}/{cfg.genetics.n_tolerance}/{cfg.genetics.n_recovery}',
        },
        'runtime_s': round(elapsed, 1),
        'peak_memory_mb': round(max_rss_mb, 0),
        'per_node': results,
        'totals': {
            'initial': total_init, 'final': total_final,
            'deaths': total_deaths, 'recoveries': total_recov,
            'recovery_rate_pct': round(recov_rate, 4),
        }
    }, f, indent=2)
print(f"\nSaved: {outdir}/results.json")
print("\nDone.")
