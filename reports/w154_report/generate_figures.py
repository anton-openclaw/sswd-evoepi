#!/usr/bin/env python3
"""Generate all figures for the W154 calibration report."""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import numpy as np
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import matplotlib.ticker as mticker
from scipy.interpolate import UnivariateSpline
from sswd_evoepi.metrics import RECOVERY_TARGETS

plt.style.use('seaborn-v0_8-whitegrid')

# ── Load data ──────────────────────────────────────────────────────────────
seeds = [42, 123, 999]
seed_colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
seed_labels = ['Seed 42', 'Seed 123', 'Seed 999']

npz_data = {}
for s in seeds:
    npz_data[s] = np.load(f'data/monthly_seed{s}.npz', allow_pickle=True)

result_data = {}
for s in seeds:
    with open(f'data/result_seed{s}.json') as f:
        result_data[s] = json.load(f)

K = 5000

# Region ordering (N → S)
regions_order = ['AK-WG','AK-AL','AK-OC','AK-EG','AK-PWS','AK-FN','AK-FS',
                 'BC-N','BC-C','SS-N','SS-S','JDF','WA-O','OR','CA-N','CA-C','CA-S','BJ']

# CENTRALIZED: moved to sswd_evoepi.metrics
targets = RECOVERY_TARGETS
# targets = {
#     'AK-PWS': 0.50, 'AK-FN': 0.50, 'AK-FS': 0.20, 'BC-N': 0.20,
#     'SS-S': 0.05, 'JDF': 0.02, 'OR': 0.0025, 'CA-N': 0.001
# }

# Build site → region mapping
site_names = npz_data[42]['site_names']
site_lats = npz_data[42]['site_lats']
site_lons = npz_data[42]['site_lons']

def get_region(name):
    parts = str(name).split('-')
    if len(parts) >= 2:
        return parts[0] + '-' + parts[1]
    return str(name)

site_regions = np.array([get_region(n) for n in site_names])

# Region mean latitudes and site counts
region_lats = {}
region_site_counts = {}
for reg in regions_order:
    mask = site_regions == reg
    if mask.any():
        region_lats[reg] = np.mean(site_lats[mask])
        region_site_counts[reg] = int(mask.sum())
    else:
        region_lats[reg] = np.nan
        region_site_counts[reg] = 0

# Month → date mapping
n_months = 159  # 2012-01 to 2025-03
years_frac = np.array([2012 + (m // 12) + (m % 12) / 12.0 for m in range(n_months)])

# ── FIGURE 1: Regional Recovery Bar Chart ──────────────────────────────────
fig1, ax1 = plt.subplots(figsize=(14, 6))

bar_width = 0.22
x = np.arange(len(regions_order))

for i, s in enumerate(seeds):
    recoveries = []
    colors = []
    for reg in regions_order:
        rec = result_data[s]['region_recovery'].get(reg, 0.0)
        rec = max(rec, 1e-5)  # floor for log scale
        recoveries.append(rec)
        
        if reg in targets:
            ratio = rec / targets[reg]
            if 0.5 <= ratio <= 2.0:
                colors.append('#2ca02c')  # green - within 50%
            elif 0.25 <= ratio <= 4.0:
                colors.append('#ffbf00')  # yellow - within 2×
            else:
                colors.append('#d62728')  # red - >2× off
        else:
            colors.append('#7f7f7f')  # gray for no target
    
    bars = ax1.bar(x + (i - 1) * bar_width, recoveries, bar_width,
                   color=colors, alpha=0.8, edgecolor='white', linewidth=0.5)

# Overlay target diamonds
for j, reg in enumerate(regions_order):
    if reg in targets:
        ax1.scatter(j, targets[reg], marker='D', s=100, color='black',
                   zorder=10, edgecolors='white', linewidth=1.5)

ax1.set_yscale('log')
ax1.set_xticks(x)
ax1.set_xticklabels(regions_order, rotation=45, ha='right', fontsize=9)
ax1.set_ylabel('Recovery Fraction', fontsize=12)
ax1.set_title('W154 Regional Recovery: Modeled vs Observed', fontsize=14, fontweight='bold')
ax1.set_ylim(1e-5, 1.0)
ax1.axhline(y=1.0, color='gray', linestyle=':', alpha=0.3)

# Legend
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='D', color='w', markerfacecolor='black',
           markeredgecolor='white', markersize=10, label='Observed target'),
    plt.Rectangle((0, 0), 1, 1, fc='#2ca02c', alpha=0.8, label='Within 50% of target'),
    plt.Rectangle((0, 0), 1, 1, fc='#ffbf00', alpha=0.8, label='Within 2× of target'),
    plt.Rectangle((0, 0), 1, 1, fc='#d62728', alpha=0.8, label='>2× from target'),
    plt.Rectangle((0, 0), 1, 1, fc='#7f7f7f', alpha=0.8, label='No target'),
]
ax1.legend(handles=legend_elements, loc='upper right', fontsize=9, framealpha=0.9)

# Add N→S arrow annotation
ax1.annotate('N → S', xy=(0.02, 0.95), xycoords='axes fraction',
            fontsize=10, fontstyle='italic', color='gray')

plt.tight_layout()
fig1.savefig('figures/fig1_regional_recovery.png', dpi=150, bbox_inches='tight')
plt.close(fig1)
print("Figure 1 saved.")

# ── FIGURE 2: Population Time Series ──────────────────────────────────────
key_regions = ['AK-PWS', 'AK-FN', 'BC-N', 'SS-S', 'CA-N']
fig2, axes = plt.subplots(5, 1, figsize=(12, 15), sharex=True)

for panel_idx, reg in enumerate(key_regions):
    ax = axes[panel_idx]
    mask = site_regions == reg
    n_sites = mask.sum()
    
    all_pops = []
    all_inf = []
    for s in seeds:
        pops = npz_data[s]['populations'][:, mask].mean(axis=1) / K
        inf = npz_data[s]['infected'][:, mask].mean(axis=1) / K
        all_pops.append(pops)
        all_inf.append(inf)
        ax.plot(years_frac, pops, color=seed_colors[seeds.index(s)],
               alpha=0.4, linewidth=1, label=f'Seed {s}')
    
    # Mean across seeds
    mean_pop = np.mean(all_pops, axis=0)
    mean_inf = np.mean(all_inf, axis=0)
    ax.plot(years_frac, mean_pop, color='black', linewidth=2.5, label='Mean', zorder=5)
    
    # Shade infected fraction
    ax.fill_between(years_frac, 0, mean_inf, alpha=0.25, color='#d62728', label='Infected (mean)')
    
    # Target level
    if reg in targets:
        ax.axhline(y=targets[reg], color='black', linestyle='--', alpha=0.5, linewidth=1)
        ax.text(2024.5, targets[reg], f'Target: {targets[reg]*100:.1f}%',
               fontsize=8, va='bottom', ha='right', color='black', alpha=0.7)
    
    # Disease arrival year
    ax.axvline(x=2013, color='#d62728', linestyle=':', alpha=0.5, linewidth=1.5)
    if panel_idx == 0:
        ax.text(2013.1, ax.get_ylim()[1] * 0.9, 'Disease\narrival',
               fontsize=8, color='#d62728', alpha=0.7, va='top')
    
    ax.set_ylabel(f'{reg}\n(N/K)', fontsize=10)
    ax.set_ylim(bottom=0)
    if panel_idx == 0:
        ax.legend(loc='upper right', fontsize=8, ncol=3, framealpha=0.9)
    ax.set_title(f'{reg} ({n_sites} sites)', fontsize=11, fontweight='bold', loc='left')

axes[-1].set_xlabel('Year', fontsize=12)
axes[-1].set_xlim(2012, 2025)
fig2.suptitle('Population Dynamics in Key Regions', fontsize=14, fontweight='bold', y=1.01)
plt.tight_layout()
fig2.savefig('figures/fig2_timeseries.png', dpi=150, bbox_inches='tight')
plt.close(fig2)
print("Figure 2 saved.")

# ── FIGURE 3: Latitudinal Recovery Gradient ───────────────────────────────
fig3, ax3 = plt.subplots(figsize=(10, 6))

# Mean recovery per region across seeds
reg_mean_recovery = {}
for reg in regions_order:
    vals = [result_data[s]['region_recovery'].get(reg, 0.0) for s in seeds]
    reg_mean_recovery[reg] = np.mean(vals)

lats_plot = []
recs_plot = []
sizes_plot = []
labels_plot = []

for reg in regions_order:
    if region_lats[reg] is not np.nan and reg_mean_recovery[reg] > 0:
        lats_plot.append(region_lats[reg])
        recs_plot.append(max(reg_mean_recovery[reg], 1e-6))
        sizes_plot.append(region_site_counts[reg])
        labels_plot.append(reg)

lats_arr = np.array(lats_plot)
recs_arr = np.array(recs_plot)
sizes_arr = np.array(sizes_plot)

# Scale sizes for visibility
size_scale = sizes_arr / sizes_arr.max() * 300 + 30

ax3.scatter(lats_arr, recs_arr, s=size_scale, c='#1f77b4', alpha=0.6,
           edgecolors='white', linewidth=1, zorder=5)

# Labels
for i, lab in enumerate(labels_plot):
    offset_y = 1.3 if recs_arr[i] < 0.01 else 0.75
    ax3.annotate(lab, (lats_arr[i], recs_arr[i]), fontsize=7,
                ha='center', va='bottom', xytext=(0, 8),
                textcoords='offset points', color='#333333')

# Target values as red diamonds
for reg, targ in targets.items():
    if reg in region_lats and not np.isnan(region_lats[reg]):
        ax3.scatter(region_lats[reg], targ, marker='D', s=80, color='#d62728',
                   edgecolors='white', linewidth=1.5, zorder=10)

# Trend line (fit in log space)
sort_idx = np.argsort(lats_arr)
lats_sorted = lats_arr[sort_idx]
log_recs_sorted = np.log10(recs_arr[sort_idx])

# Use polynomial fit
z = np.polyfit(lats_sorted, log_recs_sorted, 2)
p = np.poly1d(z)
lat_smooth = np.linspace(lats_sorted.min() - 0.5, lats_sorted.max() + 0.5, 100)
ax3.plot(lat_smooth, 10**p(lat_smooth), color='#1f77b4', linestyle='--',
        linewidth=2, alpha=0.5, zorder=3)

ax3.set_yscale('log')
ax3.set_xlabel('Latitude (°N)', fontsize=12)
ax3.set_ylabel('Recovery Fraction', fontsize=12)
ax3.set_title('Recovery vs Latitude', fontsize=14, fontweight='bold')
ax3.set_ylim(1e-6, 1.0)

legend_elements3 = [
    plt.scatter([], [], s=100, c='#1f77b4', alpha=0.6, edgecolors='white', label='Modeled (mean of 3 seeds)'),
    Line2D([0], [0], marker='D', color='w', markerfacecolor='#d62728',
           markeredgecolor='white', markersize=10, label='Observed target'),
    Line2D([0], [0], color='#1f77b4', linestyle='--', linewidth=2, alpha=0.5, label='Trend (quadratic)'),
]
ax3.legend(handles=legend_elements3, loc='lower right', fontsize=9, framealpha=0.9)

plt.tight_layout()
fig3.savefig('figures/fig3_latitudinal_gradient.png', dpi=150, bbox_inches='tight')
plt.close(fig3)
print("Figure 3 saved.")

# ── FIGURE 4: Seed Variability ────────────────────────────────────────────
target_regions = ['AK-PWS', 'AK-FN', 'AK-FS', 'BC-N', 'SS-S', 'JDF', 'OR', 'CA-N']

fig4, ax4 = plt.subplots(figsize=(10, 5))

x4 = np.arange(len(target_regions))
for i, s in enumerate(seeds):
    vals = [result_data[s]['region_recovery'].get(reg, 0.0) for reg in target_regions]
    ax4.scatter(x4 + (i - 1) * 0.15, vals, s=80, color=seed_colors[i],
               edgecolors='white', linewidth=1, zorder=5, label=f'Seed {s}')

# Target markers
for j, reg in enumerate(target_regions):
    ax4.scatter(j, targets[reg], marker='D', s=100, color='black',
               edgecolors='white', linewidth=1.5, zorder=10)
    
    # Connect seeds with a line to show range
    seed_vals = [result_data[s]['region_recovery'].get(reg, 0.0) for s in seeds]
    ax4.plot([j, j], [min(seed_vals), max(seed_vals)], color='gray',
            linewidth=1, alpha=0.5, zorder=3)

ax4.set_yscale('log')
ax4.set_xticks(x4)
ax4.set_xticklabels(target_regions, rotation=45, ha='right', fontsize=10)
ax4.set_ylabel('Recovery Fraction', fontsize=12)
ax4.set_title('Stochastic Variability Across Seeds (Target Regions)', fontsize=14, fontweight='bold')
ax4.set_ylim(1e-4, 1.0)

legend_elements4 = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor=seed_colors[0],
           markeredgecolor='white', markersize=10, label='Seed 42'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor=seed_colors[1],
           markeredgecolor='white', markersize=10, label='Seed 123'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor=seed_colors[2],
           markeredgecolor='white', markersize=10, label='Seed 999'),
    Line2D([0], [0], marker='D', color='w', markerfacecolor='black',
           markeredgecolor='white', markersize=10, label='Observed target'),
]
ax4.legend(handles=legend_elements4, loc='upper right', fontsize=9, framealpha=0.9)

plt.tight_layout()
fig4.savefig('figures/fig4_seed_variability.png', dpi=150, bbox_inches='tight')
plt.close(fig4)
print("Figure 4 saved.")

print("\nAll figures generated successfully!")
