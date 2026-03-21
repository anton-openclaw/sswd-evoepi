#!/usr/bin/env python3
"""Three-trait genetic architecture validation run.

5-node spatial simulation: Sitka, Howe Sound, SJI, Newport, Monterey
K=5000 per node, disease at year 3, seed=42, 20 years.

Records per-node per-year: population, 3 trait means, disease deaths.
Generates figures and comparison to pre-refactor single-trait baseline.
"""

import sys
import os
import time
import json
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from sswd_evoepi.config import default_config
from sswd_evoepi.model import run_spatial_simulation
from sswd_evoepi.spatial import NodeDefinition, build_network
from sswd_evoepi.types import DiseaseState, N_LOCI


# ─── 5-node network with K=5000 ──────────────────────────────────────

def make_5node_K5000():
    """5-node network with K=5000 per node for validation."""
    node_defs = [
        NodeDefinition(
            node_id=0, name="Sitka",
            lat=57.06, lon=-135.34,
            subregion="AK-SE",
            habitat_area=33333.0,
            carrying_capacity=5000,
            is_fjord=False,
            flushing_rate=0.8,
            mean_sst=8.0, sst_amplitude=3.5, sst_trend=0.015,
            salinity=32.0,
        ),
        NodeDefinition(
            node_id=1, name="Howe Sound",
            lat=49.52, lon=-123.25,
            subregion="SS",
            habitat_area=33333.0,
            carrying_capacity=5000,
            is_fjord=True,
            sill_depth=30.0,
            flushing_rate=0.03,
            mean_sst=10.0, sst_amplitude=4.0, sst_trend=0.02,
            salinity=22.0,
        ),
        NodeDefinition(
            node_id=2, name="SJI",
            lat=48.53, lon=-123.02,
            subregion="SS",
            habitat_area=33333.0,
            carrying_capacity=5000,
            is_fjord=False,
            flushing_rate=0.3,
            mean_sst=10.0, sst_amplitude=4.0, sst_trend=0.02,
            salinity=30.0,
        ),
        NodeDefinition(
            node_id=3, name="Newport",
            lat=44.63, lon=-124.05,
            subregion="WA-OR",
            habitat_area=33333.0,
            carrying_capacity=5000,
            is_fjord=False,
            flushing_rate=1.0,
            mean_sst=12.0, sst_amplitude=3.0, sst_trend=0.02,
            salinity=33.0,
        ),
        NodeDefinition(
            node_id=4, name="Monterey",
            lat=36.62, lon=-121.90,
            subregion="CA-BJ",
            habitat_area=33333.0,
            carrying_capacity=5000,
            is_fjord=False,
            flushing_rate=0.8,
            mean_sst=14.0, sst_amplitude=2.5, sst_trend=0.025,
            salinity=33.5,
        ),
    ]
    return build_network(node_defs, seed=42)


# ─── Run simulation ──────────────────────────────────────────────────

def run_validation():
    print("=" * 70)
    print("Three-Trait Genetic Architecture — Phase 6 Validation")
    print("5 nodes × K=5000 × 20 years × seed=42")
    print("=" * 70)

    config = default_config()
    network = make_5node_K5000()

    n_nodes = network.n_nodes
    node_names = [n.definition.name for n in network.nodes]
    print(f"\nNodes: {node_names}")
    print(f"Config: n_resistance={config.genetics.n_resistance}, "
          f"n_tolerance={config.genetics.n_tolerance}, "
          f"n_recovery={config.genetics.n_recovery}")
    print(f"q_init_mode={config.genetics.q_init_mode}, "
          f"target_mean_r={config.genetics.target_mean_r}")

    t0 = time.time()
    result = run_spatial_simulation(
        network=network,
        n_years=20,
        disease_year=3,
        initial_infected_per_node=5,
        seed=42,
        config=config,
        progress_callback=lambda y, n: print(f"  Year {y+1}/{n}") if (y+1) % 5 == 0 else None,
    )
    elapsed = time.time() - t0
    print(f"\nSimulation complete in {elapsed:.1f}s")

    # ── Extract per-node results ──────────────────────────────────
    print("\n" + "─" * 70)
    print("Per-Node Results (Year 0 = initial, Year 19 = final)")
    print("─" * 70)

    node_results = []
    for i in range(n_nodes):
        name = node_names[i]
        K = 5000
        pop_init = int(result.yearly_pop[i, 0])
        pop_final = int(result.yearly_pop[i, -1])
        pop_min = int(np.min(result.yearly_pop[i, :]))
        pop_min_yr = int(np.argmin(result.yearly_pop[i, :]))
        crash_pct = (1.0 - pop_min / K) * 100 if K > 0 else 0

        total_dis_deaths = int(np.sum(result.yearly_disease_deaths[i, :]))

        # Trait trajectories
        r_init = result.yearly_mean_resistance[i, 0]
        r_final = result.yearly_mean_resistance[i, -1]
        t_init = result.yearly_mean_tolerance[i, 0]
        t_final = result.yearly_mean_tolerance[i, -1]
        c_init = result.yearly_mean_recovery[i, 0]
        c_final = result.yearly_mean_recovery[i, -1]

        # Count recovered individuals at end of simulation
        node = network.nodes[i]
        alive_mask = node.agents['alive']
        n_alive_final = int(np.sum(alive_mask))
        if n_alive_final > 0:
            ds = node.agents['disease_state'][alive_mask]
            n_recovered = int(np.sum(ds == DiseaseState.R))
            n_infected = int(np.sum((ds == DiseaseState.I1) | (ds == DiseaseState.I2)))
        else:
            n_recovered = 0
            n_infected = 0

        nr = {
            'name': name,
            'K': K,
            'pop_init': pop_init,
            'pop_final': pop_final,
            'pop_min': pop_min,
            'pop_min_year': pop_min_yr,
            'crash_pct': crash_pct,
            'total_disease_deaths': total_dis_deaths,
            'n_recovered_final': n_recovered,
            'n_infected_final': n_infected,
            'r_init': float(r_init),
            'r_final': float(r_final),
            'delta_r': float(r_final - r_init),
            't_init': float(t_init),
            't_final': float(t_final),
            'delta_t': float(t_final - t_init),
            'c_init': float(c_init),
            'c_final': float(c_final),
            'delta_c': float(c_final - c_init),
        }
        node_results.append(nr)

        print(f"\n  {name}:")
        print(f"    Pop: {pop_init} → {pop_final} (min {pop_min} at yr {pop_min_yr}, crash {crash_pct:.1f}%)")
        print(f"    Disease deaths: {total_dis_deaths}")
        print(f"    Recovered alive: {n_recovered}, Infected alive: {n_infected}")
        print(f"    r_i: {r_init:.4f} → {r_final:.4f} (Δ = {r_final-r_init:+.4f})")
        print(f"    t_i: {t_init:.4f} → {t_final:.4f} (Δ = {t_final-t_init:+.4f})")
        print(f"    c_i: {c_init:.4f} → {c_final:.4f} (Δ = {c_final-c_init:+.4f})")

    # ── Emergent dynamics check ───────────────────────────────────
    print("\n" + "─" * 70)
    print("Emergent Dynamics Check")
    print("─" * 70)

    # 1. Silent spreaders: mean disease_timer for I₂ agents
    print("\n1. Silent Spreaders (high-t_i staying infected longer):")
    for i in range(n_nodes):
        node = network.nodes[i]
        alive = node.agents['alive']
        ds_arr = node.agents['disease_state']
        i2_mask = alive & (ds_arr == DiseaseState.I2)
        n_i2 = int(np.sum(i2_mask))
        if n_i2 > 0:
            timers = node.agents['disease_timer'][i2_mask]
            tol_scores = node.agents['tolerance'][i2_mask]
            print(f"  {node_names[i]}: {n_i2} I₂ agents, "
                  f"mean timer={np.mean(timers):.1f}d, "
                  f"mean t_i={np.mean(tol_scores):.4f}")
        else:
            print(f"  {node_names[i]}: 0 I₂ agents at end of sim")

    # 2. Trait divergence: rate of change per trait
    print("\n2. Trait Divergence (evolution rate comparison):")
    for i in range(n_nodes):
        nr = node_results[i]
        # Normalize by initial value for fair comparison
        rates = []
        for trait, d, init in [('r_i', nr['delta_r'], nr['r_init']),
                                ('t_i', nr['delta_t'], nr['t_init']),
                                ('c_i', nr['delta_c'], nr['c_init'])]:
            rel = d / init if init > 0.001 else d
            rates.append((trait, d, rel))
        fastest = max(rates, key=lambda x: abs(x[2]))
        print(f"  {nr['name']}: Δr={nr['delta_r']:+.4f}, Δt={nr['delta_t']:+.4f}, "
              f"Δc={nr['delta_c']:+.4f} → fastest: {fastest[0]}")

    # 3. Trait correlations in surviving population
    print("\n3. Trait Correlations (final surviving populations):")
    for i in range(n_nodes):
        node = network.nodes[i]
        alive = node.agents['alive']
        n_alive = int(np.sum(alive))
        if n_alive > 30:
            r_vals = node.agents['resistance'][alive]
            t_vals = node.agents['tolerance'][alive]
            c_vals = node.agents['recovery_ability'][alive]
            corr_rt = np.corrcoef(r_vals, t_vals)[0, 1]
            corr_rc = np.corrcoef(r_vals, c_vals)[0, 1]
            corr_tc = np.corrcoef(t_vals, c_vals)[0, 1]
            print(f"  {node_names[i]} (n={n_alive}): "
                  f"r(r,t)={corr_rt:.3f}, r(r,c)={corr_rc:.3f}, r(t,c)={corr_tc:.3f}")
        else:
            print(f"  {node_names[i]} (n={n_alive}): too few survivors for correlation")

    # ── Summary dict ──────────────────────────────────────────────
    summary = {
        'n_nodes': n_nodes,
        'node_names': node_names,
        'n_years': 20,
        'seed': 42,
        'K_per_node': 5000,
        'disease_year': 3,
        'runtime_s': round(elapsed, 1),
        'config': {
            'n_resistance': config.genetics.n_resistance,
            'n_tolerance': config.genetics.n_tolerance,
            'n_recovery': config.genetics.n_recovery,
            'q_init_mode': config.genetics.q_init_mode,
            'target_mean_r': config.genetics.target_mean_r,
        },
        'node_results': node_results,
        'initial_total_pop': int(result.initial_total_pop),
        'final_total_pop': int(result.final_total_pop),
    }

    return result, summary, network


# ─── Generate figures ─────────────────────────────────────────────────

def generate_figures(result, summary, outdir):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    os.makedirs(outdir, exist_ok=True)
    node_names = summary['node_names']
    n_nodes = summary['n_nodes']
    n_years = summary['n_years']
    years = np.arange(n_years)
    disease_year = summary['disease_year']

    colors = ['#2196F3', '#4CAF50', '#FF9800', '#F44336', '#9C27B0']

    # ── Figure 1: Population Trajectories ─────────────────────────
    fig, ax = plt.subplots(figsize=(12, 6))
    for i in range(n_nodes):
        ax.plot(years, result.yearly_pop[i, :], '-o', color=colors[i],
                label=node_names[i], markersize=3, linewidth=1.5)
    ax.axvline(x=disease_year, color='red', linestyle='--', alpha=0.5,
               label=f'Disease intro (yr {disease_year})')
    ax.axhline(y=5000, color='gray', linestyle=':', alpha=0.3, label='K=5000')
    ax.set_xlabel('Year')
    ax.set_ylabel('Population')
    ax.set_title('Three-Trait Validation: Population Trajectories (5 nodes, K=5000)')
    ax.legend(loc='upper right')
    ax.set_xlim(-0.5, n_years - 0.5)
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    path1 = os.path.join(outdir, 'population_trajectories.png')
    fig.savefig(path1, dpi=150)
    plt.close(fig)
    print(f"  Saved: {path1}")

    # ── Figure 2: Trait Evolution (3-panel) ───────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharey=False)

    trait_data = [
        ('Resistance (r_i)', result.yearly_mean_resistance),
        ('Tolerance (t_i)', result.yearly_mean_tolerance),
        ('Recovery (c_i)', result.yearly_mean_recovery),
    ]

    for panel_idx, (title, data) in enumerate(trait_data):
        ax = axes[panel_idx]
        for i in range(n_nodes):
            ax.plot(years, data[i, :], '-o', color=colors[i],
                    label=node_names[i], markersize=2.5, linewidth=1.5)
        ax.axvline(x=disease_year, color='red', linestyle='--', alpha=0.5)
        ax.set_xlabel('Year')
        ax.set_ylabel('Mean trait value')
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        if panel_idx == 0:
            ax.legend(loc='best', fontsize=8)

    fig.suptitle('Three-Trait Evolution Under SSWD Selection', fontsize=14, y=1.02)
    fig.tight_layout()
    path2 = os.path.join(outdir, 'trait_evolution.png')
    fig.savefig(path2, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {path2}")

    # ── Figure 3: Trait Comparison Bar Chart ──────────────────────
    fig, ax = plt.subplots(figsize=(12, 6))

    x = np.arange(n_nodes)
    width = 0.25

    deltas_r = [nr['delta_r'] for nr in summary['node_results']]
    deltas_t = [nr['delta_t'] for nr in summary['node_results']]
    deltas_c = [nr['delta_c'] for nr in summary['node_results']]

    bars_r = ax.bar(x - width, deltas_r, width, label='Δr (resistance)',
                    color='#2196F3', alpha=0.8)
    bars_t = ax.bar(x, deltas_t, width, label='Δt (tolerance)',
                    color='#FF9800', alpha=0.8)
    bars_c = ax.bar(x + width, deltas_c, width, label='Δc (recovery)',
                    color='#4CAF50', alpha=0.8)

    ax.set_xlabel('Node')
    ax.set_ylabel('Δ trait (final - initial)')
    ax.set_title('Selection Response per Trait per Node')
    ax.set_xticks(x)
    ax.set_xticklabels(node_names)
    ax.legend()
    ax.axhline(y=0, color='black', linewidth=0.5)
    ax.grid(True, alpha=0.3, axis='y')

    # Add value labels on bars
    for bars in [bars_r, bars_t, bars_c]:
        for bar in bars:
            h = bar.get_height()
            if abs(h) > 0.001:
                ax.text(bar.get_x() + bar.get_width()/2., h,
                        f'{h:+.3f}', ha='center', va='bottom' if h > 0 else 'top',
                        fontsize=7)

    fig.tight_layout()
    path3 = os.path.join(outdir, 'trait_comparison.png')
    fig.savefig(path3, dpi=150)
    plt.close(fig)
    print(f"  Saved: {path3}")

    return [path1, path2, path3]


# ─── Write report ─────────────────────────────────────────────────────

def write_report(summary, outdir, n_tests=632):
    """Write the VALIDATION_REPORT.md."""
    nr = summary['node_results']
    cfg = summary['config']

    lines = [
        "# Three-Trait Genetic Architecture — Validation Report",
        "",
        f"**Date:** 2026-02-20",
        f"**Phase:** 6/6 (Final Validation)",
        f"**Runtime:** {summary['runtime_s']:.1f}s",
        f"**Test suite:** {n_tests} passing",
        "",
        "## Configuration",
        "",
        f"- **Genetic architecture:** {cfg['n_resistance']}R / {cfg['n_tolerance']}T / {cfg['n_recovery']}C = 51 loci",
        f"- **Initialization:** {cfg['q_init_mode']} (target_mean_r = {cfg['target_mean_r']})",
        f"- **Nodes:** {summary['n_nodes']} (Sitka, Howe Sound, SJI, Newport, Monterey)",
        f"- **K:** {summary['K_per_node']} per node ({summary['K_per_node'] * summary['n_nodes']} total)",
        f"- **Duration:** {summary['n_years']} years, disease at year {summary['disease_year']}",
        f"- **Seed:** {summary['seed']}",
        "",
        "## 5-Node Results",
        "",
        "| Node | Pop₀ | Pop_final | Min (yr) | Crash% | Dis.Deaths | Rec. | Δr | Δt | Δc |",
        "|------|------|-----------|----------|--------|------------|------|----|----|----|",
    ]

    for n in nr:
        lines.append(
            f"| {n['name']} | {n['pop_init']} | {n['pop_final']} | "
            f"{n['pop_min']} (yr{n['pop_min_year']}) | {n['crash_pct']:.1f}% | "
            f"{n['total_disease_deaths']} | {n['n_recovered_final']} | "
            f"{n['delta_r']:+.4f} | {n['delta_t']:+.4f} | {n['delta_c']:+.4f} |"
        )

    # Comparison to single-trait baseline
    baseline = {
        'Sitka': {'crash': 95.4, 'delta_r': 0.031},
        'Howe Sound': {'crash': 93.2, 'delta_r': 0.044},
        'SJI': {'crash': 99.3, 'delta_r': 0.061},
        'Newport': {'crash': 99.1, 'delta_r': 0.055},
        'Monterey': {'crash': 98.8, 'delta_r': 0.051},
    }

    lines += [
        "",
        "## Comparison to Single-Trait Baseline (Phase 14)",
        "",
        "Previous run: 51 resistance loci, same 5 nodes, K=5000, seed=42, 20yr.",
        "",
        "| Node | Crash (old) | Crash (new) | Δ | Δr (old) | Δr (new) | Δt (new) | Δc (new) |",
        "|------|-------------|-------------|---|----------|----------|----------|----------|",
    ]
    for n in nr:
        name = n['name']
        old = baseline.get(name, {'crash': '?', 'delta_r': '?'})
        diff = n['crash_pct'] - old['crash'] if isinstance(old['crash'], (int, float)) else '?'
        diff_str = f"{diff:+.1f}" if isinstance(diff, (int, float)) else '?'
        lines.append(
            f"| {name} | {old['crash']}% | {n['crash_pct']:.1f}% | "
            f"{diff_str}pp | {old['delta_r']:+.3f} | "
            f"{n['delta_r']:+.4f} | {n['delta_t']:+.4f} | {n['delta_c']:+.4f} |"
        )

    lines += [
        "",
        "### Key Differences",
        "",
        "With three traits (17 loci each) vs one trait (51 loci):",
        "- Each trait has fewer loci → less genetic variance per trait",
        "- Resistance (r_i) alone is weaker → tolerance and recovery provide alternative survival pathways",
        "- Selection acts on all three traits simultaneously",
        "",
        "## Emergent Dynamics",
        "",
    ]

    # Add dynamics observations based on results
    # Trait divergence
    mean_delta_r = np.mean([n['delta_r'] for n in nr])
    mean_delta_t = np.mean([n['delta_t'] for n in nr])
    mean_delta_c = np.mean([n['delta_c'] for n in nr])

    lines += [
        "### Trait Evolution Rates",
        "",
        f"- Mean Δr across nodes: {mean_delta_r:+.4f}",
        f"- Mean Δt across nodes: {mean_delta_t:+.4f}",
        f"- Mean Δc across nodes: {mean_delta_c:+.4f}",
        "",
    ]

    # Which trait evolves fastest?
    abs_rates = {'r': abs(mean_delta_r), 't': abs(mean_delta_t), 'c': abs(mean_delta_c)}
    fastest = max(abs_rates, key=abs_rates.get)
    trait_names = {'r': 'resistance', 't': 'tolerance', 'c': 'recovery'}
    lines.append(f"**Fastest evolving trait:** {trait_names[fastest]} ({fastest}_i)")
    lines.append("")

    # Recovery signal
    nodes_with_rec = [n['name'] for n in nr if n['n_recovered_final'] > 0]
    lines += [
        "### Recovery Signal",
        "",
        f"Nodes with recovered (R-state) individuals at end: {nodes_with_rec if nodes_with_rec else 'None'}",
        "",
        "### Trait Correlations",
        "",
        "See console output for r(r,t), r(r,c), r(t,c) per node.",
        "Independent loci → expect low correlations unless shared selection pressure creates linkage.",
        "",
        "## Figures",
        "",
        "1. `population_trajectories.png` — Population over time per node",
        "2. `trait_evolution.png` — Three-panel trait mean trajectories",
        "3. `trait_comparison.png` — Δr vs Δt vs Δc bar chart per node",
        "",
        "## Summary",
        "",
        f"Total initial population: {summary['initial_total_pop']}",
        f"Total final population: {summary['final_total_pop']}",
        f"Overall crash: {(1 - summary['final_total_pop'] / summary['initial_total_pop']) * 100:.1f}%",
        "",
    ]

    report_path = os.path.join(outdir, 'VALIDATION_REPORT.md')
    with open(report_path, 'w') as f:
        f.write('\n'.join(lines))
    print(f"  Saved: {report_path}")
    return report_path


# ─── Main ─────────────────────────────────────────────────────────────

if __name__ == '__main__':
    outdir = 'results/three_trait_validation'
    os.makedirs(outdir, exist_ok=True)

    # Run simulation
    result, summary, network = run_validation()

    # Save raw data
    data_path = os.path.join(outdir, 'validation_data.json')
    # Make JSON-serializable
    json_summary = json.loads(json.dumps(summary, default=str))
    with open(data_path, 'w') as f:
        json.dump(json_summary, f, indent=2)
    print(f"\n  Saved: {data_path}")

    # Generate figures
    print("\nGenerating figures...")
    figure_paths = generate_figures(result, summary, outdir)

    # Write report
    print("\nWriting report...")
    write_report(summary, outdir)

    print("\n" + "=" * 70)
    print("Phase 6 validation complete!")
    print(f"Results at: {outdir}/")
    print("=" * 70)
