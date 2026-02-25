#!/usr/bin/env python3
"""Analysis and visualization for Monterey reintroduction experiment.

Generates 7 figures from experiment results JSON:

  1. Total metapopulation trajectory (all 19 scenarios, mean ± SE)
  2. Monterey release-site trajectory (population at CA-C-043)
  3. Final population by genetics × restoration level (heatmap)
  4. Genetic trait evolution at Monterey (r, t, c over time)
  5. Released individual survival curves
  6. Regional population recovery map (18 regions, final year)
  7. Dose-response: restoration level vs. final population

Usage:
    python3 experiments/analyze_reintroduction.py --results results/reintroduction/results_demo.json
    python3 experiments/analyze_reintroduction.py --results results/reintroduction/results_full.json

Authors: Anton 🔬 & Willem Weertman
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

# Project constants
START_YEAR = 2002
END_YEAR = 2050
EPIDEMIC_YEAR = 2013
RELEASE_YEAR = 2026

REGIONS = [
    'AK-WG', 'AK-AL', 'AK-EG', 'AK-PWS', 'AK-FN', 'AK-FS', 'AK-OC',
    'BC-N', 'BC-C',
    'JDF', 'SS-N', 'SS-S',
    'WA-O', 'OR',
    'CA-N', 'CA-C', 'CA-S',
    'BJ',
]

# Ordered genetic backgrounds (low → high resistance)
GENETICS_ORDER = [
    'pre_sswd', 'survivors_2019', 'bred_1gen',
    'bred_2gen', 'bred_5gen', 'optimal',
]
GENETICS_LABELS = {
    'pre_sswd': 'Pre-SSWD\n(r=0.15)',
    'survivors_2019': 'Survivors\n(r=0.18)',
    'bred_1gen': 'Bred 1gen\n(r=0.32)',
    'bred_2gen': 'Bred 2gen\n(r=0.44)',
    'bred_5gen': 'Bred 5gen\n(r=0.77)',
    'optimal': 'Optimal\n(r=1.00)',
}
GENETICS_SHORT = {
    'pre_sswd': 'Pre-SSWD',
    'survivors_2019': 'Survivors',
    'bred_1gen': '1 Gen',
    'bred_2gen': '2 Gen',
    'bred_5gen': '5 Gen',
    'optimal': 'Optimal',
}

RESTORATION_ORDER = ['partial', 'medium', 'full']
RESTORATION_LABELS = {'partial': '50', 'medium': '500', 'full': '5,000'}
RESTORATION_N = {'partial': 50, 'medium': 500, 'full': 5000}

# Color schemes
GENETICS_COLORS = {
    'pre_sswd': '#d62728',       # Red
    'survivors_2019': '#ff7f0e', # Orange
    'bred_1gen': '#2ca02c',      # Green
    'bred_2gen': '#1f77b4',      # Blue
    'bred_5gen': '#9467bd',      # Purple
    'optimal': '#17becf',        # Cyan
}

RESTORATION_COLORS = {
    'partial': '#fdae61',
    'medium': '#abdda4',
    'full': '#2b83ba',
}

# Matplotlib style
plt.rcParams.update({
    'figure.dpi': 150,
    'savefig.dpi': 200,
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
    'legend.fontsize': 8,
    'figure.facecolor': 'white',
    'axes.facecolor': '#fafafa',
    'axes.grid': True,
    'grid.alpha': 0.3,
})


# ═══════════════════════════════════════════════════════════════════════
# DATA LOADING
# ═══════════════════════════════════════════════════════════════════════


def load_results(path: str) -> List[Dict]:
    """Load experiment results JSON."""
    with open(path) as f:
        results = json.load(f)
    # Filter out errors
    valid = [r for r in results if 'error' not in r]
    errors = [r for r in results if 'error' in r]
    if errors:
        print(f"WARNING: {len(errors)} scenarios had errors:")
        for e in errors:
            print(f"  {e['scenario']} rep={e['replicate']}: {e['error']}")
    print(f"Loaded {len(valid)} valid results from {path}")
    return valid


def group_by_scenario(results: List[Dict]) -> Dict[str, List[Dict]]:
    """Group results by scenario name."""
    grouped = defaultdict(list)
    for r in results:
        grouped[r['scenario']].append(r)
    return dict(grouped)


def years_array() -> np.ndarray:
    """Return year labels for the simulation."""
    return np.arange(START_YEAR, END_YEAR)


def mean_se(values: List[np.ndarray]) -> Tuple[np.ndarray, np.ndarray]:
    """Compute mean and standard error across replicates."""
    arr = np.array(values, dtype=float)
    m = np.nanmean(arr, axis=0)
    se = np.nanstd(arr, axis=0) / np.sqrt(arr.shape[0])
    return m, se


# ═══════════════════════════════════════════════════════════════════════
# FIGURE 1: Total Metapopulation Trajectory
# ═══════════════════════════════════════════════════════════════════════


def fig1_total_trajectory(results: List[Dict], outdir: Path) -> None:
    """All 19 scenarios' total population over time."""
    grouped = group_by_scenario(results)
    years = years_array()

    fig, ax = plt.subplots(figsize=(12, 7))

    # Plot baseline
    if 'baseline' in grouped:
        pops = [r['yearly_total_pop'] for r in grouped['baseline'] if r.get('yearly_total_pop')]
        if pops:
            m, se = mean_se(pops)
            ax.plot(years[:len(m)], m, 'k-', linewidth=2.5, label='Baseline (no intervention)', zorder=10)
            ax.fill_between(years[:len(m)], m - se, m + se, color='gray', alpha=0.2)

    # Plot treatments grouped by restoration level
    for rest in RESTORATION_ORDER:
        for gen in GENETICS_ORDER:
            name = f'{rest}_{gen}'
            if name not in grouped:
                continue
            pops = [r['yearly_total_pop'] for r in grouped[name] if r.get('yearly_total_pop')]
            if not pops:
                continue
            m, se = mean_se(pops)
            color = GENETICS_COLORS[gen]
            ls = {'partial': ':', 'medium': '--', 'full': '-'}[rest]
            alpha = {'partial': 0.5, 'medium': 0.7, 'full': 1.0}[rest]
            label = f'{GENETICS_SHORT[gen]} ({RESTORATION_LABELS[rest]})'
            ax.plot(years[:len(m)], m, color=color, linestyle=ls, alpha=alpha,
                    linewidth=1.2, label=label)

    # Mark key events
    ax.axvline(EPIDEMIC_YEAR, color='red', linestyle='--', alpha=0.5, label='Epidemic onset (2013)')
    ax.axvline(RELEASE_YEAR, color='green', linestyle='--', alpha=0.5, label='Release (2026)')

    ax.set_xlabel('Year')
    ax.set_ylabel('Total Metapopulation (all 907 nodes)')
    ax.set_title('Fig 1. Metapopulation Trajectory Under Reintroduction Scenarios')
    ax.set_yscale('log')
    ax.set_xlim(START_YEAR, END_YEAR)
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f'{x:,.0f}'))
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(outdir / 'fig1_total_trajectory.png', bbox_inches='tight')
    plt.close(fig)
    print("  Fig 1: Total metapopulation trajectory ✓")


# ═══════════════════════════════════════════════════════════════════════
# FIGURE 2: Monterey Release-Site Trajectory
# ═══════════════════════════════════════════════════════════════════════


def fig2_monterey_trajectory(results: List[Dict], outdir: Path) -> None:
    """Population trajectory at CA-C-043 (Monterey)."""
    grouped = group_by_scenario(results)
    years = years_array()

    fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=True)

    for ax_i, rest in enumerate(RESTORATION_ORDER):
        ax = axes[ax_i]

        # Baseline
        if 'baseline' in grouped:
            pops = [r['monterey_yearly_pop'] for r in grouped['baseline'] if r.get('monterey_yearly_pop')]
            if pops:
                m, se = mean_se(pops)
                ax.plot(years[:len(m)], m, 'k-', linewidth=2, label='Baseline')

        # Treatments
        for gen in GENETICS_ORDER:
            name = f'{rest}_{gen}'
            if name not in grouped:
                continue
            pops = [r['monterey_yearly_pop'] for r in grouped[name] if r.get('monterey_yearly_pop')]
            if not pops:
                continue
            m, se = mean_se(pops)
            ax.plot(years[:len(m)], m, color=GENETICS_COLORS[gen], linewidth=1.5,
                    label=GENETICS_SHORT[gen])
            ax.fill_between(years[:len(m)], np.maximum(m - se, 0), m + se,
                            color=GENETICS_COLORS[gen], alpha=0.15)

        ax.axvline(RELEASE_YEAR, color='green', linestyle='--', alpha=0.5)
        ax.set_xlabel('Year')
        ax.set_title(f'{RESTORATION_LABELS[rest]} released')
        if ax_i == 0:
            ax.set_ylabel('Population at Monterey (CA-C-043)')
        ax.legend(fontsize=7)

    fig.suptitle('Fig 2. Monterey Release-Site Population Trajectory', fontsize=13)
    fig.tight_layout()
    fig.savefig(outdir / 'fig2_monterey_trajectory.png', bbox_inches='tight')
    plt.close(fig)
    print("  Fig 2: Monterey trajectory ✓")


# ═══════════════════════════════════════════════════════════════════════
# FIGURE 3: Final Population Heatmap (Genetics × Restoration)
# ═══════════════════════════════════════════════════════════════════════


def fig3_final_pop_heatmap(results: List[Dict], outdir: Path) -> None:
    """Heatmap of final (year 2050) total population."""
    grouped = group_by_scenario(results)

    # Build matrix
    matrix = np.full((len(GENETICS_ORDER), len(RESTORATION_ORDER)), np.nan)
    for gi, gen in enumerate(GENETICS_ORDER):
        for ri, rest in enumerate(RESTORATION_ORDER):
            name = f'{rest}_{gen}'
            if name not in grouped:
                continue
            finals = [r['final_total_pop'] for r in grouped[name]]
            matrix[gi, ri] = np.mean(finals)

    # Baseline for reference
    baseline_final = None
    if 'baseline' in grouped:
        baseline_final = np.mean([r['final_total_pop'] for r in grouped['baseline']])

    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(matrix, cmap='YlGnBu', aspect='auto')
    ax.set_xticks(range(len(RESTORATION_ORDER)))
    ax.set_xticklabels([RESTORATION_LABELS[r] for r in RESTORATION_ORDER])
    ax.set_yticks(range(len(GENETICS_ORDER)))
    ax.set_yticklabels([GENETICS_SHORT[g] for g in GENETICS_ORDER])
    ax.set_xlabel('Number Released')
    ax.set_ylabel('Genetic Background')

    # Annotate cells
    for gi in range(len(GENETICS_ORDER)):
        for ri in range(len(RESTORATION_ORDER)):
            val = matrix[gi, ri]
            if not np.isnan(val):
                color = 'white' if val > np.nanmax(matrix) * 0.6 else 'black'
                ax.text(ri, gi, f'{val:,.0f}', ha='center', va='center',
                        fontsize=9, color=color, fontweight='bold')

    cbar = fig.colorbar(im, ax=ax, label='Final Total Population (2050)')

    title = 'Fig 3. Final Metapopulation Size by Scenario'
    if baseline_final is not None:
        title += f'\nBaseline (no intervention): {baseline_final:,.0f}'
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(outdir / 'fig3_final_pop_heatmap.png', bbox_inches='tight')
    plt.close(fig)
    print("  Fig 3: Final population heatmap ✓")


# ═══════════════════════════════════════════════════════════════════════
# FIGURE 4: Genetic Trait Evolution at Monterey
# ═══════════════════════════════════════════════════════════════════════


def fig4_trait_evolution(results: List[Dict], outdir: Path) -> None:
    """Mean resistance, tolerance, recovery over time at Monterey."""
    grouped = group_by_scenario(results)
    years = years_array()

    # Focus on full restoration (most dramatic effect)
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    trait_keys = [
        ('monterey_yearly_mean_r', 'Resistance (r)', 0),
        ('monterey_yearly_mean_t', 'Tolerance (t)', 1),
        ('monterey_yearly_mean_c', 'Recovery (c)', 2),
    ]

    for key, label, ax_i in trait_keys:
        ax = axes[ax_i]

        # Baseline
        if 'baseline' in grouped:
            vals = [r[key] for r in grouped['baseline'] if r.get(key)]
            if vals:
                m, se = mean_se(vals)
                ax.plot(years[:len(m)], m, 'k-', linewidth=2, label='Baseline')

        for gen in GENETICS_ORDER:
            name = f'full_{gen}'
            if name not in grouped:
                continue
            vals = [r[key] for r in grouped[name] if r.get(key)]
            if not vals:
                continue
            m, se = mean_se(vals)
            ax.plot(years[:len(m)], m, color=GENETICS_COLORS[gen],
                    linewidth=1.5, label=GENETICS_SHORT[gen])
            ax.fill_between(years[:len(m)], m - se, m + se,
                            color=GENETICS_COLORS[gen], alpha=0.15)

        ax.axvline(EPIDEMIC_YEAR, color='red', linestyle='--', alpha=0.4)
        ax.axvline(RELEASE_YEAR, color='green', linestyle='--', alpha=0.4)
        ax.set_xlabel('Year')
        ax.set_ylabel(f'Mean {label}')
        ax.set_title(label)
        ax.legend(fontsize=7)

    fig.suptitle('Fig 4. Trait Evolution at Monterey (Full Restoration, 5000 Released)',
                 fontsize=13)
    fig.tight_layout()
    fig.savefig(outdir / 'fig4_trait_evolution.png', bbox_inches='tight')
    plt.close(fig)
    print("  Fig 4: Trait evolution ✓")


# ═══════════════════════════════════════════════════════════════════════
# FIGURE 5: Released Individual Survival
# ═══════════════════════════════════════════════════════════════════════


def fig5_released_survival(results: List[Dict], outdir: Path) -> None:
    """Survival of released individuals at Monterey over time."""
    grouped = group_by_scenario(results)
    years = years_array()

    fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=True)

    for ax_i, rest in enumerate(RESTORATION_ORDER):
        ax = axes[ax_i]
        n_released = RESTORATION_N[rest]

        for gen in GENETICS_ORDER:
            name = f'{rest}_{gen}'
            if name not in grouped:
                continue
            surv = [r['monterey_released_alive'] for r in grouped[name]
                    if r.get('monterey_released_alive')]
            if not surv:
                continue
            m, se = mean_se(surv)
            # Convert to fraction of released
            frac = m / n_released if n_released > 0 else m
            frac_se = se / n_released if n_released > 0 else se

            # Only plot from release year onward
            release_idx = RELEASE_YEAR - START_YEAR
            plot_years = years[release_idx:release_idx + len(frac) - release_idx]
            plot_frac = frac[release_idx:][:len(plot_years)]
            plot_se = frac_se[release_idx:][:len(plot_years)]

            ax.plot(plot_years, plot_frac, color=GENETICS_COLORS[gen],
                    linewidth=1.5, label=GENETICS_SHORT[gen])
            ax.fill_between(plot_years, np.maximum(plot_frac - plot_se, 0),
                            plot_frac + plot_se, color=GENETICS_COLORS[gen], alpha=0.15)

        ax.set_xlabel('Year')
        ax.set_title(f'{RESTORATION_LABELS[rest]} released')
        if ax_i == 0:
            ax.set_ylabel('Fraction of Released Still Alive')
        ax.set_ylim(-0.05, 1.05)
        ax.legend(fontsize=7)

    fig.suptitle('Fig 5. Survival of Released Individuals at Monterey', fontsize=13)
    fig.tight_layout()
    fig.savefig(outdir / 'fig5_released_survival.png', bbox_inches='tight')
    plt.close(fig)
    print("  Fig 5: Released survival ✓")


# ═══════════════════════════════════════════════════════════════════════
# FIGURE 6: Regional Recovery Map
# ═══════════════════════════════════════════════════════════════════════


def fig6_regional_recovery(results: List[Dict], outdir: Path) -> None:
    """Bar chart of final-year population by region for key scenarios."""
    grouped = group_by_scenario(results)

    # Compare: baseline, full_pre_sswd, full_bred_2gen, full_optimal
    compare = ['baseline', 'full_pre_sswd', 'full_bred_2gen', 'full_optimal']
    compare_labels = ['Baseline', 'Full Pre-SSWD', 'Full Bred-2gen', 'Full Optimal']
    compare_colors = ['gray', '#d62728', '#1f77b4', '#17becf']

    fig, ax = plt.subplots(figsize=(14, 6))

    x = np.arange(len(REGIONS))
    width = 0.2
    offsets = np.array([-1.5, -0.5, 0.5, 1.5]) * width

    for ci, (scenario, label, color) in enumerate(zip(compare, compare_labels, compare_colors)):
        if scenario not in grouped:
            continue

        # Average final-year regional pop across replicates
        region_finals = defaultdict(list)
        for r in grouped[scenario]:
            reg_pop = r.get('region_yearly_pop', {})
            for region in REGIONS:
                if region in reg_pop and len(reg_pop[region]) > 0:
                    region_finals[region].append(reg_pop[region][-1])

        values = []
        for region in REGIONS:
            if region in region_finals:
                values.append(np.mean(region_finals[region]))
            else:
                values.append(0)

        ax.bar(x + offsets[ci], values, width, label=label, color=color, alpha=0.8)

    ax.set_xticks(x)
    ax.set_xticklabels(REGIONS, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('Final Population (2050)')
    ax.set_title('Fig 6. Regional Population Recovery by Scenario')
    ax.legend()
    ax.set_yscale('symlog', linthresh=100)
    fig.tight_layout()
    fig.savefig(outdir / 'fig6_regional_recovery.png', bbox_inches='tight')
    plt.close(fig)
    print("  Fig 6: Regional recovery ✓")


# ═══════════════════════════════════════════════════════════════════════
# FIGURE 7: Dose-Response
# ═══════════════════════════════════════════════════════════════════════


def fig7_dose_response(results: List[Dict], outdir: Path) -> None:
    """Final population as function of restoration level, by genetics."""
    grouped = group_by_scenario(results)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Panel A: Total metapopulation
    ax = axes[0]
    for gen in GENETICS_ORDER:
        x_vals = []
        y_vals = []
        y_errs = []
        for rest in RESTORATION_ORDER:
            name = f'{rest}_{gen}'
            if name not in grouped:
                continue
            finals = [r['final_total_pop'] for r in grouped[name]]
            x_vals.append(RESTORATION_N[rest])
            y_vals.append(np.mean(finals))
            y_errs.append(np.std(finals) / np.sqrt(len(finals)))

        if x_vals:
            ax.errorbar(x_vals, y_vals, yerr=y_errs, marker='o',
                        color=GENETICS_COLORS[gen], label=GENETICS_SHORT[gen],
                        linewidth=1.5, capsize=3)

    # Baseline reference line
    if 'baseline' in grouped:
        bl = np.mean([r['final_total_pop'] for r in grouped['baseline']])
        ax.axhline(bl, color='gray', linestyle='--', alpha=0.7, label=f'Baseline ({bl:,.0f})')

    ax.set_xscale('log')
    ax.set_xlabel('Number Released')
    ax.set_ylabel('Final Total Population (2050)')
    ax.set_title('A) Total Metapopulation')
    ax.legend(fontsize=7)

    # Panel B: Monterey only
    ax = axes[1]
    for gen in GENETICS_ORDER:
        x_vals = []
        y_vals = []
        y_errs = []
        for rest in RESTORATION_ORDER:
            name = f'{rest}_{gen}'
            if name not in grouped:
                continue
            mont_finals = [
                r['monterey_yearly_pop'][-1] if r.get('monterey_yearly_pop') else 0
                for r in grouped[name]
            ]
            x_vals.append(RESTORATION_N[rest])
            y_vals.append(np.mean(mont_finals))
            y_errs.append(np.std(mont_finals) / np.sqrt(len(mont_finals)))

        if x_vals:
            ax.errorbar(x_vals, y_vals, yerr=y_errs, marker='o',
                        color=GENETICS_COLORS[gen], label=GENETICS_SHORT[gen],
                        linewidth=1.5, capsize=3)

    if 'baseline' in grouped:
        bl_mont = np.mean([
            r['monterey_yearly_pop'][-1] if r.get('monterey_yearly_pop') else 0
            for r in grouped['baseline']
        ])
        ax.axhline(bl_mont, color='gray', linestyle='--', alpha=0.7,
                    label=f'Baseline ({bl_mont:,.0f})')

    ax.set_xscale('log')
    ax.set_xlabel('Number Released')
    ax.set_ylabel('Final Monterey Population (2050)')
    ax.set_title('B) Monterey (CA-C-043)')
    ax.legend(fontsize=7)

    fig.suptitle('Fig 7. Dose-Response: Restoration Level vs. Final Population', fontsize=13)
    fig.tight_layout()
    fig.savefig(outdir / 'fig7_dose_response.png', bbox_inches='tight')
    plt.close(fig)
    print("  Fig 7: Dose-response ✓")


# ═══════════════════════════════════════════════════════════════════════
# SUMMARY STATISTICS
# ═══════════════════════════════════════════════════════════════════════


def compute_summary_stats(results: List[Dict]) -> Dict[str, Any]:
    """Compute summary statistics across all scenarios."""
    grouped = group_by_scenario(results)
    stats = {}

    # Baseline
    if 'baseline' in grouped:
        bl_finals = [r['final_total_pop'] for r in grouped['baseline']]
        bl_mont = [r['monterey_yearly_pop'][-1] for r in grouped['baseline']
                    if r.get('monterey_yearly_pop')]
        stats['baseline'] = {
            'final_pop_mean': np.mean(bl_finals),
            'final_pop_se': np.std(bl_finals) / np.sqrt(len(bl_finals)),
            'monterey_final_mean': np.mean(bl_mont) if bl_mont else 0,
            'n_replicates': len(bl_finals),
        }

    # Treatment scenarios
    for rest in RESTORATION_ORDER:
        for gen in GENETICS_ORDER:
            name = f'{rest}_{gen}'
            if name not in grouped:
                continue
            reps = grouped[name]
            finals = [r['final_total_pop'] for r in reps]
            mont_finals = [r['monterey_yearly_pop'][-1] for r in reps
                          if r.get('monterey_yearly_pop')]
            released_alive = [r.get('released_surviving', 0) for r in reps]
            elapsed = [r.get('elapsed_s', 0) for r in reps]

            stats[name] = {
                'restoration': rest,
                'genetics': gen,
                'n_released': RESTORATION_N[rest],
                'final_pop_mean': np.mean(finals),
                'final_pop_se': np.std(finals) / np.sqrt(len(finals)),
                'monterey_final_mean': np.mean(mont_finals) if mont_finals else 0,
                'released_surviving_mean': np.mean(released_alive),
                'n_replicates': len(reps),
                'mean_elapsed_s': np.mean(elapsed),
            }

    return stats


def write_summary_markdown(
    stats: Dict[str, Any],
    results: List[Dict],
    outdir: Path,
) -> None:
    """Write summary report as Markdown."""
    outfile = outdir / 'demo_summary.md'

    # Count total
    n_total = len(results)
    n_errors = len([r for r in results if 'error' in r])
    total_time = sum(r.get('elapsed_s', 0) for r in results if 'error' not in r)

    lines = [
        '# Reintroduction Experiment — Demo Summary',
        '',
        f'**Generated:** {__import__("datetime").datetime.now().strftime("%Y-%m-%d %H:%M")}',
        f'**Mode:** Demo (3 replicates per scenario)',
        f'**Network:** 907 nodes, 18 regions (Aleutians to Baja California)',
        f'**Carrying capacity:** K = 5,000 per node',
        f'**Timeline:** {START_YEAR}–{END_YEAR} ({END_YEAR - START_YEAR} years)',
        f'**Epidemic onset:** {EPIDEMIC_YEAR}',
        f'**Release:** {RELEASE_YEAR} at CA-C-043 (Monterey)',
        f'**Scenarios run:** {n_total} ({n_total - n_errors} successful, {n_errors} errors)',
        f'**Total compute time:** {total_time:.0f}s ({total_time/60:.1f} min)',
        '',
        '## Design',
        '',
        '| Factor | Levels |',
        '|--------|--------|',
        '| Restoration level | Partial (50), Medium (500), Full (5,000) |',
        '| Genetic background | Pre-SSWD, Survivors, Bred 1/2/5 gen, Optimal |',
        '| Baseline | No intervention |',
        '| **Total** | **19 scenarios × 3 replicates = 57 runs** |',
        '',
        '## Key Results',
        '',
    ]

    # Baseline
    if 'baseline' in stats:
        bl = stats['baseline']
        lines.extend([
            f'### Baseline (No Intervention)',
            f'- Final metapopulation: **{bl["final_pop_mean"]:,.0f}** ± {bl["final_pop_se"]:,.0f}',
            f'- Monterey final: **{bl["monterey_final_mean"]:,.0f}**',
            '',
        ])

    # Results table
    lines.extend([
        '### Treatment Scenarios (Final Metapopulation, Year 2050)',
        '',
        '| Genetics | 50 Released | 500 Released | 5,000 Released |',
        '|----------|------------|-------------|----------------|',
    ])

    for gen in GENETICS_ORDER:
        row = f'| {GENETICS_SHORT[gen]} |'
        for rest in RESTORATION_ORDER:
            name = f'{rest}_{gen}'
            if name in stats:
                s = stats[name]
                row += f' {s["final_pop_mean"]:,.0f} ± {s["final_pop_se"]:,.0f} |'
            else:
                row += ' — |'
        lines.append(row)

    lines.append('')

    # Monterey table
    lines.extend([
        '### Monterey (CA-C-043) Final Population',
        '',
        '| Genetics | 50 Released | 500 Released | 5,000 Released |',
        '|----------|------------|-------------|----------------|',
    ])

    for gen in GENETICS_ORDER:
        row = f'| {GENETICS_SHORT[gen]} |'
        for rest in RESTORATION_ORDER:
            name = f'{rest}_{gen}'
            if name in stats:
                row += f' {stats[name]["monterey_final_mean"]:,.0f} |'
            else:
                row += ' — |'
        lines.append(row)

    lines.append('')

    # Key findings
    lines.extend([
        '## Interpretation',
        '',
        '### Key Findings',
        '',
    ])

    # Find best and worst treatment
    treatment_stats = {k: v for k, v in stats.items() if k != 'baseline'}
    if treatment_stats:
        best_name = max(treatment_stats, key=lambda k: treatment_stats[k]['final_pop_mean'])
        worst_name = min(treatment_stats, key=lambda k: treatment_stats[k]['final_pop_mean'])
        best = treatment_stats[best_name]
        worst = treatment_stats[worst_name]

        baseline_pop = stats.get('baseline', {}).get('final_pop_mean', 0)
        best_gain = best['final_pop_mean'] - baseline_pop if baseline_pop else 0
        best_ratio = best['final_pop_mean'] / baseline_pop if baseline_pop > 0 else float('inf')

        lines.extend([
            f'1. **Best scenario:** {best_name} → {best["final_pop_mean"]:,.0f} '
            f'({best_ratio:.1f}× baseline)',
            f'2. **Worst treatment:** {worst_name} → {worst["final_pop_mean"]:,.0f}',
            f'3. **Net gain from best:** +{best_gain:,.0f} over baseline',
            '',
        ])

    # Genetics ranking
    lines.extend(['### Genetics Ranking (Full Restoration)', ''])
    for gen in GENETICS_ORDER:
        name = f'full_{gen}'
        if name in stats:
            s = stats[name]
            lines.append(f'- **{GENETICS_SHORT[gen]}**: {s["final_pop_mean"]:,.0f}')
    lines.append('')

    # Figures
    lines.extend([
        '## Figures',
        '',
        '| Figure | Description |',
        '|--------|-------------|',
        '| fig1_total_trajectory.png | Total metapopulation over time (all scenarios) |',
        '| fig2_monterey_trajectory.png | Population at CA-C-043 (Monterey) |',
        '| fig3_final_pop_heatmap.png | Final population heatmap (genetics × restoration) |',
        '| fig4_trait_evolution.png | Trait evolution at Monterey (r, t, c) |',
        '| fig5_released_survival.png | Survival curves of released individuals |',
        '| fig6_regional_recovery.png | Regional population recovery (18 regions) |',
        '| fig7_dose_response.png | Dose-response: restoration level vs outcome |',
        '',
    ])

    with open(outfile, 'w') as f:
        f.write('\n'.join(lines))
    print(f"  Summary: {outfile}")


# ═══════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════


def main():
    parser = argparse.ArgumentParser(description='Analyze reintroduction experiment results')
    parser.add_argument('--results', type=str, default='results/reintroduction/results_demo.json',
                        help='Path to results JSON')
    parser.add_argument('--outdir', type=str, default='results/reintroduction/figures',
                        help='Output directory for figures')
    args = parser.parse_args()

    # Resolve paths relative to project root
    project_root = Path(__file__).resolve().parent.parent
    results_path = project_root / args.results
    outdir = project_root / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Loading results from {results_path}")
    results = load_results(str(results_path))

    if not results:
        print("ERROR: No valid results to analyze!")
        sys.exit(1)

    print(f"\nGenerating figures in {outdir}/")

    fig1_total_trajectory(results, outdir)
    fig2_monterey_trajectory(results, outdir)
    fig3_final_pop_heatmap(results, outdir)
    fig4_trait_evolution(results, outdir)
    fig5_released_survival(results, outdir)
    fig6_regional_recovery(results, outdir)
    fig7_dose_response(results, outdir)

    print("\nComputing summary statistics...")
    stats = compute_summary_stats(results)
    write_summary_markdown(stats, results, outdir.parent)

    print("\n✅ Analysis complete!")


if __name__ == '__main__':
    main()
