#!/usr/bin/env python3
"""Generate all 8 figures + data_summary.json for Parametric vs WOA23 report."""
import sys
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from collections import Counter
from pathlib import Path

# Setup paths
REPO = Path('/home/starbot/.openclaw/workspace/sswd-evoepi')
sys.path.insert(0, str(REPO))
from sswd_evoepi.viz.salinity import _add_coastline, REGION_COLORS, MONTH_LABELS

OUTDIR = REPO / 'reports' / 'salinity-comparison'
OUTDIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------------------------
npz = np.load('/home/starbot/.openclaw/workspace/salinity_validation/woa23_site_monthly.npz',
              allow_pickle=True)
woa = npz['woa_monthly']          # (896, 12)
woa_dir = npz['woa_direct']       # (896, 12) — NaN = land
p0 = npz['param_monthly_0']       # (896, 12)
p15 = npz['param_monthly_15']     # (896, 12)
names = npz['names']
lats = npz['lats']
lons = npz['lons']
regions = npz['regions']
fjord_depth = npz['fjord_depth_norm']
nn_dist = npz['nn_distance']

N = len(names)
direct_mask = ~np.isnan(woa_dir.astype(float)).any(axis=1)

dfo = pd.read_csv('/home/starbot/.openclaw/workspace/salinity_validation/dfo_monthly_climatology.csv')
encl = pd.read_csv(REPO / 'data' / 'nodes' / 'site_enclosedness.csv')

# Region groupings
REGION_GROUPS = {
    'Alaska': ['AK-AL', 'AK-WG', 'AK-OC', 'AK-EG', 'AK-PWS', 'AK-FN', 'AK-FS'],
    'British Columbia': ['BC-N', 'BC-C'],
    'Salish Sea / WA': ['SS-N', 'SS-S', 'JDF', 'WA-O'],
    'OR / California': ['OR', 'CA-N', 'CA-C', 'CA-S', 'BJ'],
}

def group_for(region):
    for g, regs in REGION_GROUPS.items():
        if region in regs:
            return g
    return 'Other'

groups = np.array([group_for(r) for r in regions])
GROUP_COLORS = {
    'Alaska': '#1f77b4',
    'British Columbia': '#ff7f0e',
    'Salish Sea / WA': '#2ca02c',
    'OR / California': '#d62728',
}

# Precompute stats
woa_ann = woa.mean(axis=1)
p0_ann = p0.mean(axis=1)
diff_monthly = p0 - woa
bias_monthly = float(diff_monthly.mean())
rmse_monthly = float(np.sqrt((diff_monthly**2).mean()))
max_err_monthly = float(np.abs(diff_monthly).max())
bias_annual = float((p0_ann - woa_ann).mean())
rmse_annual = float(np.sqrt(((p0_ann - woa_ann)**2).mean()))

woa_range = woa.max(axis=1) - woa.min(axis=1)
min_months = woa.argmin(axis=1)
min_month_counts = Counter(int(m) for m in min_months)

# Per-region stats
per_region = {}
for r in sorted(set(regions)):
    mask = regions == r
    d = diff_monthly[mask]
    per_region[str(r)] = {
        'n_sites': int(mask.sum()),
        'bias': round(float(d.mean()), 3),
        'rmse': round(float(np.sqrt((d**2).mean())), 3),
        'woa_annual_mean': round(float(woa_ann[mask].mean()), 2),
        'param_annual_mean': round(float(p0_ann[mask].mean()), 2),
        'woa_seasonal_range_mean': round(float(woa_range[mask].mean()), 2),
    }

# Per-group stats
per_group = {}
for g, regs in REGION_GROUPS.items():
    mask = np.isin(regions, regs)
    d = diff_monthly[mask]
    per_group[g] = {
        'n_sites': int(mask.sum()),
        'bias': round(float(d.mean()), 3),
        'rmse': round(float(np.sqrt((d**2).mean())), 3),
        'woa_annual_mean': round(float(woa_ann[mask].mean()), 2),
        'woa_seasonal_range_mean': round(float(woa_range[mask].mean()), 2),
    }

# ---------------------------------------------------------------------------
# DATA SUMMARY JSON
# ---------------------------------------------------------------------------
data_summary = {
    'n_sites': N,
    'n_direct': int(direct_mask.sum()),
    'n_nn_filled': int((~direct_mask).sum()),
    'max_nn_distance_deg': round(float(nn_dist.max()), 4),
    'mean_nn_distance_deg': round(float(nn_dist.mean()), 4),
    'monthly_bias_psu': round(bias_monthly, 4),
    'monthly_rmse_psu': round(rmse_monthly, 4),
    'monthly_max_error_psu': round(max_err_monthly, 4),
    'annual_bias_psu': round(bias_annual, 4),
    'annual_rmse_psu': round(rmse_annual, 4),
    'woa_seasonal_range_mean_psu': round(float(woa_range.mean()), 4),
    'woa_seasonal_range_max_psu': round(float(woa_range.max()), 4),
    'param_seasonal_range_psu': 0.0,
    'min_month_distribution': {MONTH_LABELS[k]: v for k, v in sorted(min_month_counts.items())},
    'top_min_months': [MONTH_LABELS[k] for k in sorted(min_month_counts, key=lambda x: -min_month_counts[x])[:4]],
    'per_region': per_region,
    'per_group': per_group,
    'woa_min_psu': round(float(woa.min()), 2),
    'woa_max_psu': round(float(woa.max()), 2),
    'param_min_psu': round(float(p0.min()), 2),
    'param_max_psu': round(float(p0.max()), 2),
    'lat_range': [round(float(lats.min()), 2), round(float(lats.max()), 2)],
    'lon_range': [round(float(lons.min()), 2), round(float(lons.max()), 2)],
    'dfo_stations': len(dfo),
}

with open(OUTDIR / 'data_summary.json', 'w') as f:
    json.dump(data_summary, f, indent=2)
print("✓ data_summary.json written")

# ===========================================================================
# FIGURE STYLE
# ===========================================================================
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
    'savefig.dpi': 200,
    'savefig.bbox': 'tight',
    'savefig.facecolor': 'white',
})

# ===========================================================================
# FIG 1: Scatter — Parametric vs WOA23 (annual mean), colored by group
# ===========================================================================
fig, ax = plt.subplots(figsize=(6, 5.5))
for g in ['Alaska', 'British Columbia', 'Salish Sea / WA', 'OR / California']:
    mask = groups == g
    ax.scatter(woa_ann[mask], p0_ann[mask], s=12, alpha=0.6,
               color=GROUP_COLORS[g], label=g, edgecolors='none')
lo, hi = 22, 35
ax.plot([lo, hi], [lo, hi], 'k--', lw=1, alpha=0.5, label='1:1 line')
ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
ax.set_xlabel('WOA23 Annual Mean Salinity (psu)')
ax.set_ylabel('Parametric Annual Mean Salinity (psu)')
ax.set_title('Parametric vs WOA23: Annual Mean Salinity')
ax.legend(fontsize=8, loc='lower right')
ax.set_aspect('equal')
# Add stats annotation
ax.text(0.05, 0.95, f'Bias = {bias_annual:.2f} psu\nRMSE = {rmse_annual:.2f} psu\nn = {N}',
        transform=ax.transAxes, va='top', fontsize=9,
        bbox=dict(boxstyle='round,pad=0.3', facecolor='wheat', alpha=0.8))
fig.savefig(OUTDIR / 'fig1_scatter_annual.pdf')
fig.savefig(OUTDIR / 'fig1_scatter_annual.png')
plt.close(fig)
print("✓ Fig 1 saved")

# ===========================================================================
# FIG 2: Bias map — geographic map of annual mean bias
# ===========================================================================
fig, ax = plt.subplots(figsize=(8, 8))
bias_site = p0_ann - woa_ann
vmax = 6
sc = ax.scatter(lons, lats, c=bias_site, cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                s=14, alpha=0.8, edgecolors='none', zorder=3)
ax.set_xlim(-180, -112)
ax.set_ylim(25, 62)
_add_coastline(ax, theme='light')
cb = fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.02)
cb.set_label('Parametric − WOA23 (psu)')
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_title('Annual Mean Salinity Bias (Parametric − WOA23)')
fig.savefig(OUTDIR / 'fig2_bias_map.pdf')
fig.savefig(OUTDIR / 'fig2_bias_map.png')
plt.close(fig)
print("✓ Fig 2 saved")

# ===========================================================================
# FIG 3: Seasonal timing histogram — month of min WOA23 salinity
# ===========================================================================
fig, ax = plt.subplots(figsize=(7, 4.5))
group_order = ['Alaska', 'British Columbia', 'Salish Sea / WA', 'OR / California']
bottom = np.zeros(12)
for g in group_order:
    mask = groups == g
    counts = np.array([int((min_months[mask] == m).sum()) for m in range(12)])
    ax.bar(np.arange(12), counts, bottom=bottom, color=GROUP_COLORS[g],
           label=g, edgecolor='white', linewidth=0.5)
    bottom += counts
ax.set_xticks(range(12))
ax.set_xticklabels(MONTH_LABELS)
ax.set_xlabel('Month of Minimum WOA23 Salinity')
ax.set_ylabel('Number of Sites')
ax.set_title('Timing of Salinity Minimum Across 896 Sites')
ax.legend(fontsize=8, loc='upper right')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
fig.savefig(OUTDIR / 'fig3_min_month_hist.pdf')
fig.savefig(OUTDIR / 'fig3_min_month_hist.png')
plt.close(fig)
print("✓ Fig 3 saved")

# ===========================================================================
# FIG 4: Regional seasonal profiles (4 panels)
# ===========================================================================
months = np.arange(12)
dfo_months_cols = ['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']

fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
for ax_i, (g, regs) in enumerate(REGION_GROUPS.items()):
    ax = axes.flat[ax_i]
    mask = np.isin(regions, regs)
    woa_g = woa[mask]
    p0_g = p0[mask]

    # WOA23 mean + range shading
    woa_mean = woa_g.mean(axis=0)
    woa_lo = np.percentile(woa_g, 10, axis=0)
    woa_hi = np.percentile(woa_g, 90, axis=0)
    ax.fill_between(months, woa_lo, woa_hi, alpha=0.2, color=GROUP_COLORS[g])
    ax.plot(months, woa_mean, '-', color=GROUP_COLORS[g], lw=2.5, label='WOA23 mean')

    # Parametric mean
    p0_mean = p0_g.mean(axis=0)
    ax.plot(months, p0_mean, '--', color='gray', lw=2, label='Parametric (fw=0)')

    # DFO stations in this region
    for _, row in dfo.iterrows():
        dlat, dlon = row['lat'], row['lon']
        # Check if this DFO station is in one of the regions in this group
        # Approximate: find nearest site
        dist = np.sqrt((lats - dlat)**2 + (lons - dlon)**2)
        nearest_idx = dist.argmin()
        if regions[nearest_idx] in regs:
            vals = row[dfo_months_cols].values.astype(float)
            ax.plot(months, vals, 'o-', ms=4, lw=1, alpha=0.7,
                    color='black', zorder=5)
            ax.annotate(row['station'].split(' ')[0], (11, vals[-1]),
                       fontsize=5, alpha=0.6)

    ax.set_title(g, fontsize=11)
    ax.set_xticks(range(12))
    ax.set_xticklabels(MONTH_LABELS, fontsize=7, rotation=45)
    ax.set_ylabel('Salinity (psu)')
    if ax_i == 0:
        ax.legend(fontsize=7, loc='lower left')

fig.suptitle('Regional Seasonal Salinity Profiles: WOA23 vs Parametric', fontsize=13, y=1.02)
fig.tight_layout()
fig.savefig(OUTDIR / 'fig4_regional_profiles.pdf')
fig.savefig(OUTDIR / 'fig4_regional_profiles.png')
plt.close(fig)
print("✓ Fig 4 saved")

# ===========================================================================
# FIG 5: NN fill distance map
# ===========================================================================
fig, ax = plt.subplots(figsize=(8, 8))
# Plot direct sites first
direct_idx = np.where(direct_mask)[0]
nn_idx = np.where(~direct_mask)[0]
ax.scatter(lons[direct_idx], lats[direct_idx], c='#2ecc71', s=12, alpha=0.7,
           edgecolors='none', zorder=3, label=f'Direct WOA23 (n={len(direct_idx)})')
sc = ax.scatter(lons[nn_idx], lats[nn_idx], c=nn_dist[nn_idx], cmap='YlOrRd',
                vmin=0, vmax=0.18, s=12, alpha=0.8, edgecolors='none', zorder=3)
ax.set_xlim(-180, -112)
ax.set_ylim(25, 62)
_add_coastline(ax, theme='light')
cb = fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.02)
cb.set_label('NN Distance (degrees)')
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_title(f'WOA23 Extraction: Direct vs Nearest-Neighbour Fill')
# Custom legend
handles = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#2ecc71',
           markersize=6, label=f'Direct (n={len(direct_idx)})'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#e74c3c',
           markersize=6, label=f'NN-filled (n={len(nn_idx)})'),
]
ax.legend(handles=handles, fontsize=9, loc='lower left')
fig.savefig(OUTDIR / 'fig5_nn_fill_map.pdf')
fig.savefig(OUTDIR / 'fig5_nn_fill_map.png')
plt.close(fig)
print("✓ Fig 5 saved")

# ===========================================================================
# FIG 6: Bias by fjord depth
# ===========================================================================
fig, ax = plt.subplots(figsize=(7, 5))
abs_bias_ann = np.abs(p0_ann - woa_ann)
for g in group_order:
    mask = groups == g
    ax.scatter(fjord_depth[mask], abs_bias_ann[mask], s=14, alpha=0.5,
               color=GROUP_COLORS[g], label=g, edgecolors='none')
ax.set_xlabel('Fjord Depth (normalised)')
ax.set_ylabel('|Parametric − WOA23| Annual Mean (psu)')
ax.set_title('Parametric Bias vs Fjord Depth')
ax.legend(fontsize=8)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
# Add trend line
from numpy.polynomial import polynomial as P
c = P.polyfit(fjord_depth, abs_bias_ann, 1)
x_fit = np.linspace(0, 1, 100)
y_fit = P.polyval(x_fit, c)
ax.plot(x_fit, y_fit, 'k--', lw=1.5, alpha=0.6)
corr = np.corrcoef(fjord_depth, abs_bias_ann)[0, 1]
ax.text(0.95, 0.95, f'r = {corr:.3f}', transform=ax.transAxes,
        ha='right', va='top', fontsize=9,
        bbox=dict(boxstyle='round,pad=0.3', facecolor='wheat', alpha=0.8))
fig.savefig(OUTDIR / 'fig6_bias_fjord.pdf')
fig.savefig(OUTDIR / 'fig6_bias_fjord.png')
plt.close(fig)
print("✓ Fig 6 saved")

# ===========================================================================
# FIG 7: DFO validation — comparison
# ===========================================================================
fig, axes = plt.subplots(3, 4, figsize=(14, 9), sharex=True)
for i, (_, row) in enumerate(dfo.iterrows()):
    ax = axes.flat[i]
    dfo_vals = row[dfo_months_cols].values.astype(float)

    # Find nearest model site
    dist = np.sqrt((lats - row['lat'])**2 + (lons - row['lon'])**2)
    nearest = dist.argmin()
    woa_vals = woa[nearest]
    p0_vals = p0[nearest]

    ax.plot(months, dfo_vals, 'ko-', ms=4, lw=1.5, label='DFO obs')
    ax.plot(months, woa_vals, 's-', ms=3, lw=1.2, color='#1f77b4', label='WOA23')
    ax.plot(months, p0_vals, 'd--', ms=3, lw=1.2, color='#d62728', label='Parametric')

    ax.set_title(row['station'], fontsize=8)
    ax.set_xticks(range(0, 12, 2))
    ax.set_xticklabels([MONTH_LABELS[m] for m in range(0, 12, 2)], fontsize=6)
    if i % 4 == 0:
        ax.set_ylabel('Salinity (psu)', fontsize=8)
    if i == 0:
        ax.legend(fontsize=6, loc='lower left')

fig.suptitle('DFO Lighthouse Stations: Observed vs WOA23 vs Parametric', fontsize=13, y=1.01)
fig.tight_layout()
fig.savefig(OUTDIR / 'fig7_dfo_validation.pdf')
fig.savefig(OUTDIR / 'fig7_dfo_validation.png')
plt.close(fig)
print("✓ Fig 7 saved")

# ===========================================================================
# FIG 8: Seasonal range map
# ===========================================================================
fig, ax = plt.subplots(figsize=(8, 8))
sc = ax.scatter(lons, lats, c=woa_range, cmap='viridis', vmin=0, vmax=8,
                s=14, alpha=0.8, edgecolors='none', zorder=3)
ax.set_xlim(-180, -112)
ax.set_ylim(25, 62)
_add_coastline(ax, theme='light')
cb = fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.02)
cb.set_label('WOA23 Seasonal Range (max − min, psu)')
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_title('WOA23 Seasonal Salinity Range Across 896 Sites')
ax.text(0.05, 0.05, f'Mean range: {woa_range.mean():.2f} psu\nMax range: {woa_range.max():.2f} psu',
        transform=ax.transAxes, fontsize=9,
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
fig.savefig(OUTDIR / 'fig8_seasonal_range_map.pdf')
fig.savefig(OUTDIR / 'fig8_seasonal_range_map.png')
plt.close(fig)
print("✓ Fig 8 saved")

# ===========================================================================
# Additional stats for data_summary.json
# ===========================================================================
# DFO comparison stats
dfo_rmse_woa = []
dfo_rmse_param = []
for _, row in dfo.iterrows():
    dist = np.sqrt((lats - row['lat'])**2 + (lons - row['lon'])**2)
    nearest = dist.argmin()
    dfo_vals = row[dfo_months_cols].values.astype(float)
    woa_vals = woa[nearest]
    p0_vals = p0[nearest]
    dfo_rmse_woa.append(np.sqrt(((woa_vals - dfo_vals)**2).mean()))
    dfo_rmse_param.append(np.sqrt(((p0_vals - dfo_vals)**2).mean()))

data_summary['dfo_rmse_woa_mean'] = round(float(np.mean(dfo_rmse_woa)), 3)
data_summary['dfo_rmse_param_mean'] = round(float(np.mean(dfo_rmse_param)), 3)
data_summary['dfo_rmse_woa_per_station'] = {row['station']: round(float(r), 3)
                                             for (_, row), r in zip(dfo.iterrows(), dfo_rmse_woa)}
data_summary['dfo_rmse_param_per_station'] = {row['station']: round(float(r), 3)
                                               for (_, row), r in zip(dfo.iterrows(), dfo_rmse_param)}

# Correlation between fjord depth and bias
data_summary['fjord_bias_correlation'] = round(float(np.corrcoef(fjord_depth, abs_bias_ann)[0, 1]), 4)

# Worst sites
worst_idx = np.argsort(abs_bias_ann)[-10:][::-1]
data_summary['worst_10_sites'] = [
    {'name': str(names[i]), 'region': str(regions[i]),
     'bias': round(float(p0_ann[i] - woa_ann[i]), 2),
     'woa_ann': round(float(woa_ann[i]), 2),
     'param_ann': round(float(p0_ann[i]), 2)}
    for i in worst_idx
]

# Best sites (smallest bias)
best_idx = np.argsort(abs_bias_ann)[:10]
data_summary['best_10_sites'] = [
    {'name': str(names[i]), 'region': str(regions[i]),
     'bias': round(float(p0_ann[i] - woa_ann[i]), 2),
     'woa_ann': round(float(woa_ann[i]), 2),
     'param_ann': round(float(p0_ann[i]), 2)}
    for i in best_idx
]

with open(OUTDIR / 'data_summary.json', 'w') as f:
    json.dump(data_summary, f, indent=2)
print("✓ data_summary.json updated with DFO + worst/best sites")
print("\nDone! All 8 figures + data_summary.json generated.")
