#!/usr/bin/env python3
"""Generate all 14 figures for the Reintroduction Experiment 2050 report."""

import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import FancyBboxPatch
from pathlib import Path

# ─── Setup ────────────────────────────────────────────────────────────────
DATA_PATH = Path('/home/starbot/.openclaw/workspace/sswd-evoepi/reports/reintroduction_2050/data/data_summary.json')
OUT_DIR = Path('/home/starbot/.openclaw/workspace/sswd-evoepi/reports/reintroduction_2050/figures')
OUT_DIR.mkdir(parents=True, exist_ok=True)

with open(DATA_PATH) as f:
    D = json.load(f)

# ─── Constants ────────────────────────────────────────────────────────────
N_YEARS = len(D['baseline']['regions']['CA-C']['yearly_totals'])  # 38
YEARS = np.arange(2012, 2012 + N_YEARS)  # 2012-2049
RELEASE_IDX = 13  # 2025
EPIDEMIC_IDX = 1   # 2013

R_KEYS = ['R1', 'R2', 'R3', 'R4', 'R5']
D_KEYS = ['D1', 'D2', 'D3', 'D4']
R_LABELS = {
    'R1': 'R1: CA-2025 (r̄=0.31)',
    'R2': 'R2: Naive (r̄=0.15)',
    'R3': 'R3: Mid-epidemic (r̄=0.20)',
    'R4': 'R4: Selected (r̄=0.50)',
    'R5': 'R5: Immune (r̄=1.0)'
}
R_SHORT = {'R1': 'CA-2025\n(r̄=0.31)', 'R2': 'Naive\n(r̄=0.15)', 'R3': 'Mid-2019\n(r̄=0.20)',
            'R4': 'Selected\n(r̄=0.50)', 'R5': 'Immune\n(r̄=1.0)'}
R_RMEAN = {'R1': 0.31, 'R2': 0.15, 'R3': 0.20, 'R4': 0.50, 'R5': 1.0}
D_N = {'D1': 50, 'D2': 200, 'D3': 500, 'D4': 1000}
D_LABELS = {'D1': '50/site', 'D2': '200/site', 'D3': '500/site', 'D4': '1000/site (K)'}

# Colors
R_COLORS = {'R1': '#d62728', 'R2': '#ff7f0e', 'R3': '#DAA520', 'R4': '#4682B4', 'R5': '#2ca02c'}
D_COLORS = {'D1': '#9467bd', 'D2': '#1f77b4', 'D3': '#ff7f0e', 'D4': '#d62728'}

DESIGN = D['design']
REGIONS_ORDERED = DESIGN['regions_ordered']
KEY_REGIONS = ['CA-S', 'CA-C', 'CA-N', 'OR', 'JDF', 'AK-PWS', 'BC-N']

def scenario_key(loc, res, den):
    return f'REINTRO_{loc}_{res}_{den}'

def get_region_recovery(scenario_name, region):
    if scenario_name == 'baseline':
        return D['baseline']['regions'][region]['recovery_pct']
    return D['scenarios'][scenario_name]['regions'][region]['recovery_pct']

def get_yearly(scenario_name, region):
    if scenario_name == 'baseline':
        return np.array(D['baseline']['regions'][region]['yearly_totals'], dtype=float)
    return np.array(D['scenarios'][scenario_name]['regions'][region]['yearly_totals'], dtype=float)

def add_event_lines(ax):
    """Add epidemic and release vertical lines."""
    ax.axvline(YEARS[EPIDEMIC_IDX], color='red', ls='--', alpha=0.5, lw=1, label='Epidemic (2013)')
    ax.axvline(YEARS[RELEASE_IDX], color='green', ls='--', alpha=0.5, lw=1, label='Release (2025)')

def style_ax(ax, xlabel='', ylabel='', title=''):
    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)
    if title:
        ax.set_title(title, fontsize=11, fontweight='bold')
    ax.grid(True, alpha=0.3, color='gray')
    ax.tick_params(labelsize=9)

def savefig(fig, name):
    path = OUT_DIR / name
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"  ✓ {name}")


# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 1: Factorial Heatmap
# ═══════════════════════════════════════════════════════════════════════════
def fig_01():
    print("Fig 01: Factorial heatmap...")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    for ax, loc, title_suffix, metric_key in [
        (ax1, 'L1', 'L1: California Mean Recovery %', 'ca_total_pct'),
        (ax2, 'L2', 'L2: JDF Recovery %', 'jdf_pct')
    ]:
        hmap = D[f'{loc.lower()}_heatmap']
        data = np.zeros((5, 4))
        for i, r in enumerate(R_KEYS):
            for j, d in enumerate(D_KEYS):
                key = f'{r}_{d}'
                data[i, j] = hmap[key][metric_key]

        vmin, vmax = 0, max(data.max(), 1)
        im = ax.imshow(data, cmap='YlGn', aspect='auto', vmin=vmin, vmax=vmax)

        # Annotate cells
        for i in range(5):
            for j in range(4):
                val = data[i, j]
                color = 'white' if val > vmax * 0.7 else 'black'
                ax.text(j, i, f'{val:.1f}%', ha='center', va='center', fontsize=9, color=color, fontweight='bold')

        ax.set_xticks(range(4))
        ax.set_xticklabels([D_LABELS[d] for d in D_KEYS], fontsize=9)
        ax.set_yticks(range(5))
        ax.set_yticklabels([R_SHORT[r] for r in R_KEYS], fontsize=9)
        ax.set_title(title_suffix, fontsize=11, fontweight='bold')

    # Baseline annotation
    bl_ca = np.mean([D['baseline']['regions'][r]['recovery_pct'] for r in ['CA-S', 'CA-C', 'CA-N']])
    bl_jdf = D['baseline']['regions']['JDF']['recovery_pct']
    ax1.text(0.5, -0.12, f'Baseline CA mean: {bl_ca:.2f}%', transform=ax1.transAxes,
             ha='center', fontsize=9, style='italic', color='gray')
    ax2.text(0.5, -0.12, f'Baseline JDF: {bl_jdf:.2f}%', transform=ax2.transAxes,
             ha='center', fontsize=9, style='italic', color='gray')

    fig.colorbar(im, ax=[ax1, ax2], label='Recovery %', shrink=0.8)
    fig.suptitle('Recovery by Genetics × Density', fontsize=13, fontweight='bold', y=1.02)
    fig.tight_layout()
    savefig(fig, 'fig_01_factorial_heatmap.pdf')


# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 2: Baseline Comparison Bar Chart
# ═══════════════════════════════════════════════════════════════════════════
def fig_02():
    print("Fig 02: Baseline comparison...")
    fig, ax = plt.subplots(figsize=(12, 6))

    scenarios = {
        'Baseline': 'baseline',
        'L1: R5/D4 (Immune, CA)': scenario_key('L1', 'R5', 'D4'),
        'L2: R5/D4 (Immune, WA)': scenario_key('L2', 'R5', 'D4'),
    }
    colors = {'Baseline': '#888888', 'L1: R5/D4 (Immune, CA)': '#2ca02c', 'L2: R5/D4 (Immune, WA)': '#1f77b4'}

    x = np.arange(len(KEY_REGIONS))
    width = 0.25
    for i, (label, sname) in enumerate(scenarios.items()):
        vals = [get_region_recovery(sname, r) for r in KEY_REGIONS]
        bars = ax.bar(x + i * width - width, vals, width, label=label, color=colors[label], edgecolor='white', linewidth=0.5)
        # Add value labels on bars
        for bar, val in zip(bars, vals):
            if val > 0.5:
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                        f'{val:.1f}', ha='center', va='bottom', fontsize=7)

    ax.set_xticks(x)
    ax.set_xticklabels(KEY_REGIONS, fontsize=10)
    style_ax(ax, ylabel='Recovery %', title='Regional Recovery: Baseline vs Best Reintroduction Scenarios')
    ax.legend(fontsize=9, loc='upper left')
    fig.tight_layout()
    savefig(fig, 'fig_02_baseline_comparison.pdf')


# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 3: CA-C Trajectories (L1)
# ═══════════════════════════════════════════════════════════════════════════
def fig_03():
    print("Fig 03: CA-C trajectories (L1)...")
    fig, ax = plt.subplots(figsize=(10, 6))

    lines = [
        ('baseline', 'Baseline', 'black', '--', 2),
        (scenario_key('L1', 'R1', 'D2'), 'L1-R1/D2 (CA-2025, 200)', R_COLORS['R1'], '-', 1.5),
        (scenario_key('L1', 'R3', 'D2'), 'L1-R3/D2 (Mid-2019, 200)', R_COLORS['R3'], '-', 1.5),
        (scenario_key('L1', 'R4', 'D2'), 'L1-R4/D2 (Selected, 200)', R_COLORS['R4'], '-', 1.5),
        (scenario_key('L1', 'R5', 'D2'), 'L1-R5/D2 (Immune, 200)', R_COLORS['R5'], '-', 1.5),
        (scenario_key('L1', 'R5', 'D4'), 'L1-R5/D4 (Immune, 1000)', R_COLORS['R5'], '-', 2.5),
    ]

    for sname, label, color, ls, lw in lines:
        yt = get_yearly(sname, 'CA-C')
        ax.plot(YEARS, yt, color=color, ls=ls, lw=lw, label=label)

    add_event_lines(ax)
    style_ax(ax, xlabel='Year', ylabel='CA-C Total Population', title='CA-C Population Trajectories (L1 Releases)')
    ax.legend(fontsize=8, loc='upper right')
    ax.set_xlim(YEARS[0], YEARS[-1])
    fig.tight_layout()
    savefig(fig, 'fig_03_ca_trajectories.pdf')


# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 4: JDF Trajectories (L2)
# ═══════════════════════════════════════════════════════════════════════════
def fig_04():
    print("Fig 04: JDF trajectories (L2)...")
    fig, ax = plt.subplots(figsize=(10, 6))

    lines = [
        ('baseline', 'Baseline', 'black', '--', 2),
        (scenario_key('L2', 'R1', 'D2'), 'L2-R1/D2 (CA-2025, 200)', R_COLORS['R1'], '-', 1.5),
        (scenario_key('L2', 'R4', 'D2'), 'L2-R4/D2 (Selected, 200)', R_COLORS['R4'], '-', 1.5),
        (scenario_key('L2', 'R5', 'D2'), 'L2-R5/D2 (Immune, 200)', R_COLORS['R5'], '-', 1.5),
        (scenario_key('L2', 'R5', 'D4'), 'L2-R5/D4 (Immune, 1000)', R_COLORS['R5'], '-', 2.5),
    ]

    for sname, label, color, ls, lw in lines:
        yt = get_yearly(sname, 'JDF')
        ax.plot(YEARS, yt, color=color, ls=ls, lw=lw, label=label)

    add_event_lines(ax)
    style_ax(ax, xlabel='Year', ylabel='JDF Total Population', title='JDF Population Trajectories (L2 Releases)')
    ax.legend(fontsize=8, loc='upper right')
    ax.set_xlim(YEARS[0], YEARS[-1])
    fig.tight_layout()
    savefig(fig, 'fig_04_jdf_trajectories.pdf')


# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 5: Competitive Dilution (KEY FIGURE)
# ═══════════════════════════════════════════════════════════════════════════
def fig_05():
    print("Fig 05: Competitive dilution (KEY FIGURE)...")
    fig, ax = plt.subplots(figsize=(10, 6))

    densities = [D_N[d] for d in D_KEYS]
    baseline_cac = D['baseline']['regions']['CA-C']['recovery_pct']

    for r in R_KEYS:
        vals = []
        for d in D_KEYS:
            sname = scenario_key('L1', r, d)
            vals.append(get_region_recovery(sname, 'CA-C'))
        ax.plot(densities, vals, 'o-', color=R_COLORS[r], lw=2, markersize=8, label=R_LABELS[r])

    ax.axhline(baseline_cac, color='black', ls='--', lw=1.5, alpha=0.7, label=f'Baseline ({baseline_cac:.2f}%)')

    # Highlight the dilution zone
    ax.fill_between([400, 1050], -0.5, baseline_cac, alpha=0.08, color='red')
    ax.annotate('Competitive\nDilution Zone', xy=(700, baseline_cac * 0.3), fontsize=9,
                ha='center', color='red', style='italic', alpha=0.7)

    # Arrow showing decrease for R1
    r1_d2 = get_region_recovery(scenario_key('L1', 'R1', 'D2'), 'CA-C')
    r1_d4 = get_region_recovery(scenario_key('L1', 'R1', 'D4'), 'CA-C')
    ax.annotate('', xy=(1000, r1_d4), xytext=(200, r1_d2),
                arrowprops=dict(arrowstyle='->', color='red', lw=1.5, ls='--'))

    style_ax(ax, xlabel='Release Density (animals/site)', ylabel='CA-C Recovery %',
             title='Competitive Dilution: More Animals Can HARM Recovery')
    ax.legend(fontsize=8, loc='upper left', framealpha=0.9)
    ax.set_xticks(densities)
    ax.set_xticklabels([str(d) for d in densities])
    ax.set_ylim(bottom=-0.5)
    fig.tight_layout()
    savefig(fig, 'fig_05_competitive_dilution.pdf')


# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 6: Resistance Threshold
# ═══════════════════════════════════════════════════════════════════════════
def fig_06():
    print("Fig 06: Resistance threshold...")
    fig, ax = plt.subplots(figsize=(10, 6))

    r_means = [R_RMEAN[r] for r in R_KEYS]
    baseline_cac = D['baseline']['regions']['CA-C']['recovery_pct']

    for d_key in D_KEYS:
        vals = []
        for r in R_KEYS:
            sname = scenario_key('L1', r, d_key)
            vals.append(get_region_recovery(sname, 'CA-C'))
        ax.plot(r_means, vals, 'o-', color=D_COLORS[d_key], lw=2, markersize=8, label=D_LABELS[d_key])

    ax.axhline(baseline_cac, color='black', ls='--', lw=1.5, alpha=0.7, label=f'Baseline ({baseline_cac:.2f}%)')

    # Threshold annotation
    ax.axvspan(0.31, 0.50, alpha=0.1, color='orange')
    ax.annotate('Threshold\nZone', xy=(0.405, ax.get_ylim()[1] * 0.5), fontsize=9,
                ha='center', color='darkorange', style='italic')

    style_ax(ax, xlabel='Mean Resistance (r̄)', ylabel='CA-C Recovery %',
             title='Resistance Threshold for CA-C Recovery (L1)')
    ax.legend(fontsize=8, loc='upper left')
    ax.set_xticks(r_means)
    ax.set_xticklabels([f'{v:.2f}' for v in r_means])
    fig.tight_layout()
    savefig(fig, 'fig_06_resistance_threshold.pdf')


# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 7: Dose Response (two panels)
# ═══════════════════════════════════════════════════════════════════════════
def fig_07():
    print("Fig 07: Dose response...")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    densities = [D_N[d] for d in D_KEYS]

    # Panel a: L1 CA mean
    bl_ca = np.mean([D['baseline']['regions'][r]['recovery_pct'] for r in ['CA-S', 'CA-C', 'CA-N']])
    for r in R_KEYS:
        vals = []
        for d in D_KEYS:
            hm = D['l1_heatmap'][f'{r}_{d}']
            vals.append(hm['ca_total_pct'])
        ax1.plot(densities, vals, 'o-', color=R_COLORS[r], lw=2, markersize=7, label=R_LABELS[r])
    ax1.axhline(bl_ca, color='black', ls='--', lw=1.5, alpha=0.7, label=f'Baseline ({bl_ca:.2f}%)')
    ax1.set_xscale('log')
    ax1.set_xticks(densities)
    ax1.set_xticklabels([str(d) for d in densities])
    style_ax(ax1, xlabel='Release Density (log scale)', ylabel='CA Mean Recovery %',
             title='(a) L1: California Mean Recovery')
    ax1.legend(fontsize=7, loc='upper left')

    # Panel b: L2 JDF
    bl_jdf = D['baseline']['regions']['JDF']['recovery_pct']
    for r in R_KEYS:
        vals = []
        for d in D_KEYS:
            hm = D['l2_heatmap'][f'{r}_{d}']
            vals.append(hm['jdf_pct'])
        ax2.plot(densities, vals, 'o-', color=R_COLORS[r], lw=2, markersize=7, label=R_LABELS[r])
    ax2.axhline(bl_jdf, color='black', ls='--', lw=1.5, alpha=0.7, label=f'Baseline ({bl_jdf:.2f}%)')
    ax2.set_xscale('log')
    ax2.set_xticks(densities)
    ax2.set_xticklabels([str(d) for d in densities])
    style_ax(ax2, xlabel='Release Density (log scale)', ylabel='JDF Recovery %',
             title='(b) L2: JDF Recovery')
    ax2.legend(fontsize=7, loc='upper left')

    fig.tight_layout()
    savefig(fig, 'fig_07_dose_response.pdf')


# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 8: Spillover Heatmap
# ═══════════════════════════════════════════════════════════════════════════
def fig_08():
    print("Fig 08: Spillover heatmap...")
    fig, ax = plt.subplots(figsize=(12, 8))

    col_scenarios = [
        ('Baseline', 'baseline'),
        ('L1-R5/D4', scenario_key('L1', 'R5', 'D4')),
        ('L2-R5/D4', scenario_key('L2', 'R5', 'D4')),
        ('L1-R4/D2', scenario_key('L1', 'R4', 'D2')),
        ('L2-R4/D2', scenario_key('L2', 'R4', 'D2')),
    ]

    data = np.zeros((len(REGIONS_ORDERED), len(col_scenarios)))
    for j, (label, sname) in enumerate(col_scenarios):
        for i, reg in enumerate(REGIONS_ORDERED):
            data[i, j] = get_region_recovery(sname, reg)

    im = ax.imshow(data, cmap='YlGn', aspect='auto', vmin=0, vmax=data.max())

    for i in range(len(REGIONS_ORDERED)):
        for j in range(len(col_scenarios)):
            val = data[i, j]
            color = 'white' if val > data.max() * 0.65 else 'black'
            ax.text(j, i, f'{val:.1f}', ha='center', va='center', fontsize=8, color=color)

    ax.set_xticks(range(len(col_scenarios)))
    ax.set_xticklabels([c[0] for c in col_scenarios], fontsize=9, rotation=30, ha='right')
    ax.set_yticks(range(len(REGIONS_ORDERED)))
    ax.set_yticklabels(REGIONS_ORDERED, fontsize=9)
    style_ax(ax, title='Recovery % by Region and Scenario (Spillover Analysis)')
    fig.colorbar(im, ax=ax, label='Recovery %', shrink=0.7)
    fig.tight_layout()
    savefig(fig, 'fig_08_spillover_heatmap.pdf')


# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 9: AK-PWS Stability
# ═══════════════════════════════════════════════════════════════════════════
def fig_09():
    print("Fig 09: AK-PWS stability...")
    fig, ax = plt.subplots(figsize=(10, 6))

    bl_yt = get_yearly('baseline', 'AK-PWS')
    all_finals = []

    for sname in sorted(D['scenarios'].keys()):
        yt = get_yearly(sname, 'AK-PWS')
        ax.plot(YEARS, yt, color='gray', alpha=0.3, lw=0.8)
        all_finals.append(yt[-1])

    ax.plot(YEARS, bl_yt, color='black', lw=2.5, label='Baseline')

    # Stats band
    mean_final = np.mean(all_finals)
    std_final = np.std(all_finals)
    ax.fill_between(YEARS, bl_yt * 0.99, bl_yt * 1.01, alpha=0.15, color='blue', label='±1% band')

    add_event_lines(ax)
    style_ax(ax, xlabel='Year', ylabel='AK-PWS Total Population',
             title=f'AK-PWS Stability Across All 40 Scenarios (final mean={mean_final:.0f}±{std_final:.0f})')
    ax.legend(fontsize=9)
    ax.set_xlim(YEARS[0], YEARS[-1])

    # Recovery annotation
    ak_rec = D['baseline']['regions']['AK-PWS']['recovery_pct']
    ax.text(0.95, 0.05, f'Recovery: ~{ak_rec:.1f}% (invariant)', transform=ax.transAxes,
            ha='right', fontsize=10, style='italic', color='navy',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', alpha=0.8))

    fig.tight_layout()
    savefig(fig, 'fig_09_ak_stability.pdf')


# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 10: Regional Map
# ═══════════════════════════════════════════════════════════════════════════
def fig_10():
    print("Fig 10: Regional map...")
    # Try to load site coordinates from an NPZ file
    npz_path = Path('/home/starbot/.openclaw/workspace/sswd-evoepi/results/reintro_2050/REINTRO_L1_R5_D4/monthly_seed42.npz')
    bl_npz_path = Path('/home/starbot/.openclaw/workspace/sswd-evoepi/results/forecast_2050/F01_baseline/monthly_seed42.npz')

    try:
        if npz_path.exists():
            npz = np.load(npz_path, allow_pickle=True)
        elif bl_npz_path.exists():
            npz = np.load(bl_npz_path, allow_pickle=True)
        else:
            raise FileNotFoundError("No NPZ found")

        lats = npz['site_lats']
        lons = npz['site_lons']
        pops = npz['populations']  # (time, sites)

        # Final populations (last timestep)
        final_pop_bl = pops[-1, :]  # This is for L1_R5_D4 or baseline

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))

        # Load baseline populations
        if bl_npz_path.exists():
            bl_npz = np.load(bl_npz_path, allow_pickle=True)
            bl_final = bl_npz['populations'][-1, :]
        else:
            bl_final = final_pop_bl

        # Load best scenario
        best_npz = np.load(npz_path, allow_pickle=True)
        best_final = best_npz['populations'][-1, :]

        for ax, pop_data, title in [(ax1, bl_final, '(a) Baseline'), (ax2, best_final, '(b) L1-R5/D4 (Best)')]:
            # Color by recovery fraction (pop / K)
            K = 1000
            recovery_frac = np.clip(pop_data / K, 0, 1)
            sizes = np.clip(pop_data / 5, 1, 50)

            sc = ax.scatter(lons, lats, c=recovery_frac, s=sizes, cmap='YlGn',
                           vmin=0, vmax=1, alpha=0.7, edgecolors='gray', linewidth=0.3)
            ax.set_title(title, fontsize=11, fontweight='bold')
            ax.set_xlabel('Longitude', fontsize=10)
            ax.set_ylabel('Latitude', fontsize=10)
            ax.grid(True, alpha=0.2)

            # Add region labels for key regions
            # Approximate region centroids
            region_centroids = {
                'AK-PWS': (-147, 60.5), 'CA-C': (-122, 36.5), 'CA-N': (-124, 40),
                'CA-S': (-120, 34), 'JDF': (-123.5, 48.3), 'OR': (-124.5, 44),
                'BC-N': (-129, 53)
            }
            for reg, (rlon, rlat) in region_centroids.items():
                if rlon > lons.min() - 1 and rlon < lons.max() + 1:
                    ax.annotate(reg, (rlon, rlat), fontsize=7, ha='center', color='navy', fontweight='bold')

        fig.colorbar(sc, ax=[ax1, ax2], label='Recovery Fraction (pop/K)', shrink=0.7)
        fig.suptitle('Geographic Distribution of Recovery', fontsize=13, fontweight='bold', y=1.02)
        fig.tight_layout()
        savefig(fig, 'fig_10_regional_map.pdf')

    except Exception as e:
        print(f"  ⚠ Fig 10 error: {e}")
        # Fallback: simple schematic based on region recovery data
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        # Approximate lat/lon for each region (south to north)
        region_coords = {
            'BJ': (-117, 30), 'CA-S': (-120, 34), 'CA-C': (-122, 36.5),
            'CA-N': (-124, 40), 'OR': (-124.5, 44), 'WA-O': (-124.5, 47),
            'JDF': (-123.5, 48.3), 'SS-S': (-123, 49), 'SS-N': (-125, 50),
            'BC-C': (-128, 51.5), 'BC-N': (-130, 53), 'AK-FS': (-133, 56),
            'AK-FN': (-135, 58), 'AK-OC': (-150, 59), 'AK-PWS': (-147, 60.5),
            'AK-EG': (-140, 59), 'AK-WG': (-155, 57), 'AK-AL': (-165, 54)
        }

        for ax, sname, title in [(ax1, 'baseline', '(a) Baseline'), (ax2, scenario_key('L1', 'R5', 'D4'), '(b) L1-R5/D4')]:
            for reg in REGIONS_ORDERED:
                lon, lat = region_coords[reg]
                rec = get_region_recovery(sname, reg)
                color_val = rec / 30.0  # normalize
                ax.scatter(lon, lat, s=max(rec * 15, 20), c=[plt.cm.YlGn(min(color_val, 1))],
                          edgecolors='gray', linewidth=0.5)
                ax.annotate(f'{reg}\n{rec:.1f}%', (lon, lat), fontsize=6, ha='center', va='bottom')
            ax.set_title(title, fontsize=11, fontweight='bold')
            ax.set_xlabel('Longitude')
            ax.set_ylabel('Latitude')
            ax.grid(True, alpha=0.2)

        fig.suptitle('Geographic Distribution of Recovery', fontsize=13, fontweight='bold')
        fig.tight_layout()
        savefig(fig, 'fig_10_regional_map.pdf')


# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 11: Interaction Plot
# ═══════════════════════════════════════════════════════════════════════════
def fig_11():
    print("Fig 11: Interaction plot...")
    fig, ax = plt.subplots(figsize=(10, 6))

    densities = [D_N[d] for d in D_KEYS]
    baseline_cac = D['baseline']['regions']['CA-C']['recovery_pct']

    for r in R_KEYS:
        vals = []
        for d in D_KEYS:
            sname = scenario_key('L1', r, d)
            vals.append(get_region_recovery(sname, 'CA-C'))
        ax.plot(densities, vals, 'o-', color=R_COLORS[r], lw=2, markersize=8, label=R_LABELS[r])

    ax.axhline(baseline_cac, color='black', ls='--', lw=1.5, alpha=0.7, label=f'Baseline ({baseline_cac:.2f}%)')

    # Annotate crossover
    ax.annotate('Crossover:\nR4-R5 increase,\nR1-R3 decrease',
                xy=(750, 5), fontsize=9, ha='center', style='italic',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', alpha=0.8))

    style_ax(ax, xlabel='Release Density (animals/site)', ylabel='CA-C Recovery %',
             title='Genetics × Density Interaction (L1, CA-C)')
    ax.legend(fontsize=8, loc='upper left')
    ax.set_xticks(densities)
    fig.tight_layout()
    savefig(fig, 'fig_11_interaction_plot.pdf')


# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 12: CA Region Comparison (3-panel)
# ═══════════════════════════════════════════════════════════════════════════
def fig_12():
    print("Fig 12: CA region comparison...")
    fig, axes = plt.subplots(1, 3, figsize=(14, 5), sharey=True)

    ca_regions = ['CA-S', 'CA-C', 'CA-N']
    for ax, reg in zip(axes, ca_regions):
        baseline_val = get_region_recovery('baseline', reg)
        vals = [baseline_val]
        labels_list = ['Baseline']
        colors_list = ['#888888']

        for r in R_KEYS:
            sname = scenario_key('L1', r, 'D2')
            vals.append(get_region_recovery(sname, reg))
            labels_list.append(r)
            colors_list.append(R_COLORS[r])

        bars = ax.bar(range(len(vals)), vals, color=colors_list, edgecolor='white', linewidth=0.5)
        for bar, val in zip(bars, vals):
            if val > 0.01:
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                        f'{val:.1f}', ha='center', va='bottom', fontsize=8)

        ax.set_xticks(range(len(labels_list)))
        ax.set_xticklabels(labels_list, fontsize=8, rotation=45, ha='right')
        style_ax(ax, ylabel='Recovery %' if reg == 'CA-S' else '', title=reg)

    fig.suptitle('California Region Recovery at D2 (200/site, L1)', fontsize=12, fontweight='bold')
    fig.tight_layout()
    savefig(fig, 'fig_12_ca_region_comparison.pdf')


# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 13: Genetics Evolution
# ═══════════════════════════════════════════════════════════════════════════
def fig_13():
    print("Fig 13: Genetics evolution...")
    fig, ax = plt.subplots(figsize=(10, 6))

    # All scenarios share the same pre-release genetics (indices 0-12)
    # Use R2/D4 as "baseline" genetics trace since baseline data doesn't have genetics
    scenarios_for_genetics = [
        (scenario_key('L1', 'R2', 'D4'), 'R2/D4: Naive flood', R_COLORS['R2'], '-', 1.5),
        (scenario_key('L1', 'R1', 'D2'), 'R1/D2: CA-2025', R_COLORS['R1'], '-', 1.5),
        (scenario_key('L1', 'R4', 'D2'), 'R4/D2: Selected', R_COLORS['R4'], '-', 1.5),
        (scenario_key('L1', 'R5', 'D2'), 'R5/D2: Immune', R_COLORS['R5'], '-', 1.5),
        (scenario_key('L1', 'R5', 'D4'), 'R5/D4: Immune (K)', R_COLORS['R5'], '-', 2.5),
    ]

    for sname, label, color, ls, lw in scenarios_for_genetics:
        r_data = D['scenarios'][sname]['regions']['CA-C']['yearly_mean_resistance']
        ax.plot(YEARS[:len(r_data)], r_data, color=color, ls=ls, lw=lw, label=label)

    add_event_lines(ax)

    # Reference lines for input r̄ values
    for r_key, rmean in [('R2', 0.15), ('R4', 0.50), ('R5', 1.0)]:
        ax.axhline(rmean, color=R_COLORS[r_key], ls=':', alpha=0.3, lw=0.8)
        ax.text(YEARS[-1] + 0.5, rmean, f'r̄={rmean}', fontsize=7, color=R_COLORS[r_key], va='center')

    style_ax(ax, xlabel='Year', ylabel='Mean Resistance Trait (CA-C)',
             title='Genetic Evolution of Resistance at CA-C After Release')
    ax.legend(fontsize=8, loc='right')
    ax.set_xlim(YEARS[0], YEARS[-1] + 3)
    ax.set_ylim(0, 1.05)
    fig.tight_layout()
    savefig(fig, 'fig_13_genetics_evolution.pdf')


# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 14: Summary Table
# ═══════════════════════════════════════════════════════════════════════════
def fig_14():
    print("Fig 14: Summary table...")

    # Build table data
    headers = ['Scenario', 'Loc', 'Resistance', 'Density', 'CA-S%', 'CA-C%', 'CA-N%', 'JDF%', 'AK-PWS%', 'Overall%']
    rows = []

    # Baseline row
    bl = D['baseline']
    bl_vals = [
        'Baseline', '—', '—', '—',
        f"{bl['regions']['CA-S']['recovery_pct']:.1f}",
        f"{bl['regions']['CA-C']['recovery_pct']:.1f}",
        f"{bl['regions']['CA-N']['recovery_pct']:.1f}",
        f"{bl['regions']['JDF']['recovery_pct']:.1f}",
        f"{bl['regions']['AK-PWS']['recovery_pct']:.1f}",
        f"{bl['final_pop_pct']:.1f}"
    ]
    rows.append(bl_vals)

    # Numeric data for coloring
    numeric_data = []
    bl_numeric = [
        bl['regions']['CA-S']['recovery_pct'],
        bl['regions']['CA-C']['recovery_pct'],
        bl['regions']['CA-N']['recovery_pct'],
        bl['regions']['JDF']['recovery_pct'],
        bl['regions']['AK-PWS']['recovery_pct'],
        bl['final_pop_pct']
    ]
    numeric_data.append(bl_numeric)

    for loc in ['L1', 'L2']:
        for r in R_KEYS:
            for d_key in D_KEYS:
                sname = scenario_key(loc, r, d_key)
                s = D['scenarios'][sname]
                r_info = DESIGN['resistance'][r]
                d_info = DESIGN['density'][d_key]

                ca_s = s['regions']['CA-S']['recovery_pct']
                ca_c = s['regions']['CA-C']['recovery_pct']
                ca_n = s['regions']['CA-N']['recovery_pct']
                jdf = s['regions']['JDF']['recovery_pct']
                ak = s['regions']['AK-PWS']['recovery_pct']
                overall = s.get('final_pop_pct', 0)
                if 'final_pop_pct' not in s:
                    # Calculate overall as mean of all region recoveries
                    overall = np.mean([s['regions'][reg]['recovery_pct'] for reg in REGIONS_ORDERED])

                row = [
                    sname.replace('REINTRO_', ''),
                    loc,
                    f"{r} (r̄={r_info['r_mean']})",
                    f"{d_info['n']}",
                    f"{ca_s:.1f}", f"{ca_c:.1f}", f"{ca_n:.1f}",
                    f"{jdf:.1f}", f"{ak:.1f}", f"{overall:.1f}"
                ]
                rows.append(row)
                numeric_data.append([ca_s, ca_c, ca_n, jdf, ak, overall])

    # Create figure
    n_rows = len(rows)
    fig_height = max(12, n_rows * 0.3 + 2)
    fig, ax = plt.subplots(figsize=(16, fig_height))
    ax.axis('off')

    # Create table
    table = ax.table(cellText=rows, colLabels=headers, loc='center', cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(7)
    table.scale(1, 1.2)

    # Color-code numeric columns (CA-S, CA-C, CA-N, JDF, AK-PWS, Overall)
    num_col_start = 4  # columns 4-9 are numeric
    cmap = plt.cm.YlGn
    max_vals = [max(nd[i] for nd in numeric_data) for i in range(6)]

    for row_idx in range(n_rows):
        for col_idx in range(6):
            val = numeric_data[row_idx][col_idx]
            max_v = max_vals[col_idx] if max_vals[col_idx] > 0 else 1
            norm_val = min(val / max_v, 1.0)
            cell = table[row_idx + 1, num_col_start + col_idx]
            cell.set_facecolor(cmap(norm_val * 0.7))
            if norm_val > 0.6:
                cell.set_text_props(color='white')

    # Header styling
    for j in range(len(headers)):
        cell = table[0, j]
        cell.set_facecolor('#4472C4')
        cell.set_text_props(color='white', fontweight='bold')

    # Baseline row highlighting
    for j in range(len(headers)):
        cell = table[1, j]
        cell.set_edgecolor('black')
        cell.set_linewidth(1.5)

    ax.set_title('Summary of All 40 Reintroduction Scenarios + Baseline',
                 fontsize=13, fontweight='bold', pad=20)
    fig.tight_layout()
    savefig(fig, 'fig_14_summary_table.pdf')


# ═══════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.size': 10,
        'axes.titlesize': 11,
        'figure.facecolor': 'white',
    })

    print("=" * 60)
    print("Generating all 14 figures for Reintroduction Experiment 2050")
    print("=" * 60)

    errors = []
    for i, func in enumerate([fig_01, fig_02, fig_03, fig_04, fig_05, fig_06,
                                fig_07, fig_08, fig_09, fig_10, fig_11, fig_12,
                                fig_13, fig_14], 1):
        try:
            func()
        except Exception as e:
            import traceback
            err_msg = f"Figure {i:02d} FAILED: {e}"
            print(f"  ✗ {err_msg}")
            traceback.print_exc()
            errors.append(err_msg)

    print("\n" + "=" * 60)
    if errors:
        print(f"Completed with {len(errors)} error(s):")
        for e in errors:
            print(f"  - {e}")
    else:
        print("All 14 figures generated successfully!")
    print("=" * 60)
