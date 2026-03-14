#!/usr/bin/env python3
"""Generate all 9 figures + data_summary.json for Parametric vs WOA23 report.

Figures match BRIEF.md specification exactly:
  1. Scatter: Parametric vs WOA23 annual mean
  2. Bias map (geographic)
  3. Seasonal timing map (geographic, colored by month of min salinity)
  4. Regional seasonal profiles (4-panel)
  5. Seasonal range comparison (scatter)
  6. NN fill quality map
  7. DFO validation (12-panel)
  8. Suppression comparison (side-by-side maps)
  9. Latitude gradient (3-panel)
"""
import sys, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.colors as mcolors
from collections import Counter
from pathlib import Path

# Setup
REPO = Path('/home/starbot/.openclaw/workspace/sswd-evoepi')
sys.path.insert(0, str(REPO))
from sswd_evoepi.viz.salinity import _add_coastline

OUTDIR = REPO / 'reports' / 'salinity-comparison'
FIGDIR = OUTDIR / 'figures'
FIGDIR.mkdir(parents=True, exist_ok=True)

MONTH_LABELS = ['Jan','Feb','Mar','Apr','May','Jun',
                'Jul','Aug','Sep','Oct','Nov','Dec']

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
dfo_months_cols = ['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']

# Region groupings
REGION_GROUPS = {
    'Alaska': ['AK-AL', 'AK-WG', 'AK-OC', 'AK-EG', 'AK-PWS', 'AK-FN', 'AK-FS'],
    'British Columbia': ['BC-N', 'BC-C'],
    'Salish Sea / WA': ['SS-N', 'SS-S', 'JDF', 'WA-O'],
    'OR / California': ['OR', 'CA-N', 'CA-C', 'CA-S', 'BJ'],
}

def group_for(r):
    for g, regs in REGION_GROUPS.items():
        if r in regs: return g
    return 'Other'

groups = np.array([group_for(r) for r in regions])
GROUP_COLORS = {
    'Alaska': '#1f77b4',
    'British Columbia': '#ff7f0e',
    'Salish Sea / WA': '#2ca02c',
    'OR / California': '#d62728',
}
GROUP_ORDER = ['Alaska', 'British Columbia', 'Salish Sea / WA', 'OR / California']

# Precompute stats
woa_ann = woa.mean(axis=1)
p0_ann = p0.mean(axis=1)
p15_ann = p15.mean(axis=1)
diff_monthly = p0 - woa
bias_monthly = float(diff_monthly.mean())
rmse_monthly = float(np.sqrt((diff_monthly**2).mean()))
max_err_monthly = float(np.abs(diff_monthly).max())
bias_annual = float((p0_ann - woa_ann).mean())
rmse_annual = float(np.sqrt(((p0_ann - woa_ann)**2).mean()))
abs_bias_ann = np.abs(p0_ann - woa_ann)

# R² calculation
ss_res = ((p0_ann - woa_ann)**2).sum()
ss_tot = ((woa_ann - woa_ann.mean())**2).sum()
r_squared = 1 - ss_res / ss_tot

woa_range = woa.max(axis=1) - woa.min(axis=1)
p0_range = p0.max(axis=1) - p0.min(axis=1)  # should be ~0
p15_range = p15.max(axis=1) - p15.min(axis=1)
min_months = woa.argmin(axis=1)
min_month_counts = Counter(int(m) for m in min_months)

# Salinity modifier function
def sal_modifier(S, s_min=10.0, s_full=28.0, eta=2.0):
    S = np.asarray(S, dtype=float)
    out = np.ones_like(S)
    mask_low = S <= s_min
    mask_mid = (S > s_min) & (S < s_full)
    out[mask_low] = 0.0
    frac = (S[mask_mid] - s_min) / (s_full - s_min)
    out[mask_mid] = frac ** eta
    return out

# Suppression: June values (month index 5)
woa_june = woa[:, 5]
p0_june = p0[:, 5]
sal_mod_woa_june = sal_modifier(woa_june)
sal_mod_p0_june = sal_modifier(p0_june)
suppression_woa_june = (1.0 - sal_mod_woa_june) * 100
suppression_p0_june = (1.0 - sal_mod_p0_june) * 100

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

# DFO comparison
dfo_rmse_woa = {}
dfo_rmse_param = {}
for _, row in dfo.iterrows():
    dist = np.sqrt((lats - row['lat'])**2 + (lons - row['lon'])**2)
    nearest = dist.argmin()
    dfo_vals = row[dfo_months_cols].values.astype(float)
    woa_vals = woa[nearest]
    p0_vals = p0[nearest]
    dfo_rmse_woa[row['station']] = round(float(np.sqrt(((woa_vals - dfo_vals)**2).mean())), 3)
    dfo_rmse_param[row['station']] = round(float(np.sqrt(((p0_vals - dfo_vals)**2).mean())), 3)

# Worst/best sites
worst_idx = np.argsort(abs_bias_ann)[-10:][::-1]
best_idx = np.argsort(abs_bias_ann)[:10]

# Fjord depth correlation
fjord_bias_corr = float(np.corrcoef(fjord_depth, abs_bias_ann)[0, 1])

# Suppression stats
n_suppressed_woa = int((suppression_woa_june > 1.0).sum())
n_suppressed_param = int((suppression_p0_june > 1.0).sum())
mean_supp_woa = float(suppression_woa_june.mean())
mean_supp_param = float(suppression_p0_june.mean())
max_supp_woa = float(suppression_woa_june.max())
max_supp_param = float(suppression_p0_june.max())

# Per-group suppression
supp_per_group = {}
for g, regs in REGION_GROUPS.items():
    mask = np.isin(regions, regs)
    supp_per_group[g] = {
        'mean_suppression_woa_pct': round(float(suppression_woa_june[mask].mean()), 2),
        'mean_suppression_param_pct': round(float(suppression_p0_june[mask].mean()), 2),
        'max_suppression_woa_pct': round(float(suppression_woa_june[mask].max()), 2),
        'max_suppression_param_pct': round(float(suppression_p0_june[mask].max()), 2),
        'woa_june_mean': round(float(woa_june[mask].mean()), 2),
        'param_june_mean': round(float(p0_june[mask].mean()), 2),
    }

# Latitude gradient stats
lat_corr_woa = float(np.corrcoef(lats, woa_june)[0, 1])
lat_corr_param = float(np.corrcoef(lats, p0_june)[0, 1])

# ---------------------------------------------------------------------------
# DATA SUMMARY JSON
# ---------------------------------------------------------------------------
data_summary = {
    'n_sites': N,
    'n_direct': int(direct_mask.sum()),
    'n_nn_filled': int((~direct_mask).sum()),
    'pct_nn_filled': round(100 * (~direct_mask).sum() / N, 1),
    'max_nn_distance_deg': round(float(nn_dist.max()), 4),
    'mean_nn_distance_deg': round(float(nn_dist.mean()), 4),
    'monthly_bias_psu': round(bias_monthly, 4),
    'monthly_rmse_psu': round(rmse_monthly, 4),
    'monthly_max_error_psu': round(max_err_monthly, 4),
    'annual_bias_psu': round(bias_annual, 4),
    'annual_rmse_psu': round(rmse_annual, 4),
    'annual_r_squared': round(float(r_squared), 4),
    'woa_seasonal_range_mean_psu': round(float(woa_range.mean()), 4),
    'woa_seasonal_range_max_psu': round(float(woa_range.max()), 3),
    'woa_seasonal_range_median_psu': round(float(np.median(woa_range)), 4),
    'param_seasonal_range_psu': round(float(p0_range.mean()), 4),
    'param_fw15_seasonal_range_mean_psu': round(float(p15_range.mean()), 4),
    'min_month_distribution': {MONTH_LABELS[k]: v for k, v in sorted(min_month_counts.items())},
    'top_min_months': [MONTH_LABELS[k] for k in sorted(min_month_counts,
                       key=lambda x: -min_month_counts[x])[:4]],
    'per_region': per_region,
    'per_group': per_group,
    'woa_min_psu': round(float(woa.min()), 2),
    'woa_max_psu': round(float(woa.max()), 2),
    'param_min_psu': round(float(p0.min()), 2),
    'param_max_psu': round(float(p0.max()), 2),
    'lat_range': [round(float(lats.min()), 2), round(float(lats.max()), 2)],
    'lon_range': [round(float(lons.min()), 2), round(float(lons.max()), 2)],
    'dfo_stations': len(dfo),
    'dfo_rmse_woa_mean': round(float(np.mean(list(dfo_rmse_woa.values()))), 3),
    'dfo_rmse_param_mean': round(float(np.mean(list(dfo_rmse_param.values()))), 3),
    'dfo_rmse_woa_per_station': dfo_rmse_woa,
    'dfo_rmse_param_per_station': dfo_rmse_param,
    'fjord_bias_correlation': round(fjord_bias_corr, 4),
    'suppression_june': {
        'n_suppressed_woa': n_suppressed_woa,
        'n_suppressed_param': n_suppressed_param,
        'mean_suppression_woa_pct': round(mean_supp_woa, 2),
        'mean_suppression_param_pct': round(mean_supp_param, 2),
        'max_suppression_woa_pct': round(max_supp_woa, 2),
        'max_suppression_param_pct': round(max_supp_param, 2),
        'per_group': supp_per_group,
    },
    'latitude_correlation': {
        'woa_june_vs_lat': round(lat_corr_woa, 4),
        'param_june_vs_lat': round(lat_corr_param, 4),
    },
    'worst_10_sites': [
        {'name': str(names[i]), 'region': str(regions[i]),
         'bias': round(float(p0_ann[i] - woa_ann[i]), 2),
         'woa_ann': round(float(woa_ann[i]), 2),
         'param_ann': round(float(p0_ann[i]), 2)}
        for i in worst_idx
    ],
    'best_10_sites': [
        {'name': str(names[i]), 'region': str(regions[i]),
         'bias': round(float(p0_ann[i] - woa_ann[i]), 2),
         'woa_ann': round(float(woa_ann[i]), 2),
         'param_ann': round(float(p0_ann[i]), 2)}
        for i in best_idx
    ],
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
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.facecolor': 'white',
})

# ===========================================================================
# FIG 1: Scatter — Parametric vs WOA23 annual mean
# ===========================================================================
print("Fig 1: scatter annual mean...")
fig, ax = plt.subplots(figsize=(6, 5.5))
for g in GROUP_ORDER:
    mask = groups == g
    ax.scatter(woa_ann[mask], p0_ann[mask], s=14, alpha=0.6,
               color=GROUP_COLORS[g], label=g, edgecolors='none')
lo, hi = 22, 35
ax.plot([lo, hi], [lo, hi], 'k--', lw=1, alpha=0.5, label='1:1 line')
ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
ax.set_xlabel('WOA23 Annual Mean Salinity (psu)')
ax.set_ylabel('Parametric Annual Mean Salinity (psu)')
ax.set_title('Parametric vs WOA23: Annual Mean Salinity')
ax.legend(fontsize=8, loc='lower right')
ax.set_aspect('equal')
ax.text(0.05, 0.95,
        f'Bias = {bias_annual:.2f} psu\nRMSE = {rmse_annual:.2f} psu\n$R^2$ = {r_squared:.3f}\nn = {N}',
        transform=ax.transAxes, va='top', fontsize=9,
        bbox=dict(boxstyle='round,pad=0.3', facecolor='wheat', alpha=0.8))
fig.savefig(FIGDIR / 'fig1_scatter_annual.pdf')
plt.close(fig)
print("  ✓ Fig 1")

# ===========================================================================
# FIG 2: Bias map
# ===========================================================================
print("Fig 2: bias map...")
fig, ax = plt.subplots(figsize=(8, 8))
bias_site = p0_ann - woa_ann
vmax = 6
sc = ax.scatter(lons, lats, c=bias_site, cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                s=14, alpha=0.8, edgecolors='none', zorder=3)
ax.set_xlim(-180, -112); ax.set_ylim(25, 62)
_add_coastline(ax, theme='light')
cb = fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.02)
cb.set_label('Parametric − WOA23 (psu)')
ax.set_xlabel('Longitude'); ax.set_ylabel('Latitude')
ax.set_title('Annual Mean Salinity Bias (Parametric − WOA23)')
fig.savefig(FIGDIR / 'fig2_bias_map.pdf')
plt.close(fig)
print("  ✓ Fig 2")

# ===========================================================================
# FIG 3: Seasonal timing MAP — geographic scatter colored by month of min
# ===========================================================================
print("Fig 3: seasonal timing map...")
# Create a cyclic colormap for months
month_cmap = plt.cm.hsv
fig, ax = plt.subplots(figsize=(8, 8))
sc = ax.scatter(lons, lats, c=min_months, cmap='hsv', vmin=-0.5, vmax=11.5,
                s=14, alpha=0.8, edgecolors='none', zorder=3)
ax.set_xlim(-180, -112); ax.set_ylim(25, 62)
_add_coastline(ax, theme='light')
cb = fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.02, ticks=range(12))
cb.ax.set_yticklabels(MONTH_LABELS, fontsize=8)
cb.set_label('Month of Minimum WOA23 Salinity')
ax.set_xlabel('Longitude'); ax.set_ylabel('Latitude')
ax.set_title('Timing of Salinity Minimum (WOA23)')
fig.savefig(FIGDIR / 'fig3_timing_map.pdf')
plt.close(fig)
print("  ✓ Fig 3")

# Also make the histogram (useful supplementary)
print("Fig 3b: timing histogram...")
fig, ax = plt.subplots(figsize=(7, 4.5))
bottom = np.zeros(12)
for g in GROUP_ORDER:
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
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
fig.savefig(FIGDIR / 'fig3b_timing_hist.pdf')
plt.close(fig)
print("  ✓ Fig 3b")

# ===========================================================================
# FIG 4: Regional seasonal profiles (4-panel)
# ===========================================================================
print("Fig 4: regional profiles...")
months_arr = np.arange(12)
fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
for ax_i, (g, regs) in enumerate(REGION_GROUPS.items()):
    ax = axes.flat[ax_i]
    mask = np.isin(regions, regs)
    woa_g = woa[mask]
    p0_g = p0[mask]
    p15_g = p15[mask]

    # WOA23 mean + range shading
    woa_mean = woa_g.mean(axis=0)
    woa_lo = np.percentile(woa_g, 10, axis=0)
    woa_hi = np.percentile(woa_g, 90, axis=0)
    ax.fill_between(months_arr, woa_lo, woa_hi, alpha=0.2, color=GROUP_COLORS[g])
    ax.plot(months_arr, woa_mean, '-', color=GROUP_COLORS[g], lw=2.5, label='WOA23 mean')

    # Parametric baseline (fw=0) — dashed
    p0_mean = p0_g.mean(axis=0)
    ax.plot(months_arr, p0_mean, '--', color='gray', lw=2, label='Parametric (fw=0)')

    # Parametric with fw_strength=15 — solid thin
    p15_mean = p15_g.mean(axis=0)
    ax.plot(months_arr, p15_mean, '-', color='gray', lw=1.5, alpha=0.7,
            label='Parametric (fw=15)')

    # DFO stations in this region
    for _, row in dfo.iterrows():
        dist = np.sqrt((lats - row['lat'])**2 + (lons - row['lon'])**2)
        nearest_idx = dist.argmin()
        if regions[nearest_idx] in regs:
            vals = row[dfo_months_cols].values.astype(float)
            ax.plot(months_arr, vals, 'o-', ms=3, lw=0.8, alpha=0.6,
                    color='black', zorder=5)

    ax.set_title(f'{g} (n={int(mask.sum())})', fontsize=11)
    ax.set_xticks(range(12))
    ax.set_xticklabels(MONTH_LABELS, fontsize=7, rotation=45)
    ax.set_ylabel('Salinity (psu)')
    if ax_i == 0:
        ax.legend(fontsize=7, loc='lower left')

fig.suptitle('Regional Seasonal Salinity Profiles: WOA23 vs Parametric', fontsize=13, y=1.02)
fig.tight_layout()
fig.savefig(FIGDIR / 'fig4_regional_profiles.pdf')
plt.close(fig)
print("  ✓ Fig 4")

# ===========================================================================
# FIG 5: Seasonal range comparison scatter
# ===========================================================================
print("Fig 5: seasonal range comparison...")
fig, axes = plt.subplots(1, 2, figsize=(11, 5))

# Panel a: WOA23 range vs parametric baseline range (fw=0)
ax = axes[0]
for g in GROUP_ORDER:
    mask = groups == g
    ax.scatter(woa_range[mask], p0_range[mask], s=14, alpha=0.6,
               color=GROUP_COLORS[g], label=g, edgecolors='none')
ax.plot([0, 9], [0, 9], 'k--', lw=1, alpha=0.3)
ax.set_xlabel('WOA23 Seasonal Range (psu)')
ax.set_ylabel('Parametric (fw=0) Seasonal Range (psu)')
ax.set_title('(a) Parametric baseline: zero seasonality')
ax.set_xlim(-0.2, 9); ax.set_ylim(-0.2, 9)
ax.legend(fontsize=7, loc='upper left')
ax.text(0.5, 0.5, 'All sites = 0.0 psu',
        transform=ax.transAxes, ha='center', va='center',
        fontsize=14, color='red', alpha=0.6, fontweight='bold')

# Panel b: WOA23 range vs parametric fw=15 range
ax = axes[1]
for g in GROUP_ORDER:
    mask = groups == g
    ax.scatter(woa_range[mask], p15_range[mask], s=14, alpha=0.6,
               color=GROUP_COLORS[g], label=g, edgecolors='none')
ax.plot([0, 9], [0, 9], 'k--', lw=1, alpha=0.3)
ax.set_xlabel('WOA23 Seasonal Range (psu)')
ax.set_ylabel('Parametric (fw=15) Seasonal Range (psu)')
ax.set_title('(b) With melt pulse (fw=15)')
ax.set_xlim(-0.2, 9); ax.set_ylim(-0.2, 9)
ax.legend(fontsize=7, loc='upper left')

fig.suptitle('Seasonal Range: WOA23 vs Parametric Model', fontsize=13)
fig.tight_layout()
fig.savefig(FIGDIR / 'fig5_seasonal_range.pdf')
plt.close(fig)
print("  ✓ Fig 5")

# ===========================================================================
# FIG 6: NN fill quality map
# ===========================================================================
print("Fig 6: NN fill map...")
fig, ax = plt.subplots(figsize=(8, 8))
direct_idx = np.where(direct_mask)[0]
nn_idx = np.where(~direct_mask)[0]
ax.scatter(lons[direct_idx], lats[direct_idx], c='#2ecc71', s=12, alpha=0.7,
           edgecolors='none', zorder=3, label=f'Direct WOA23 (n={len(direct_idx)})')
sc = ax.scatter(lons[nn_idx], lats[nn_idx], c=nn_dist[nn_idx], cmap='YlOrRd',
                vmin=0, vmax=0.18, s=12, alpha=0.8, edgecolors='none', zorder=3)
ax.set_xlim(-180, -112); ax.set_ylim(25, 62)
_add_coastline(ax, theme='light')
cb = fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.02)
cb.set_label('NN Distance (degrees)')
ax.set_xlabel('Longitude'); ax.set_ylabel('Latitude')
ax.set_title('WOA23 Extraction: Direct vs Nearest-Neighbour Fill')
handles = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#2ecc71',
           markersize=6, label=f'Direct (n={len(direct_idx)})'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#e74c3c',
           markersize=6, label=f'NN-filled (n={len(nn_idx)})'),
]
ax.legend(handles=handles, fontsize=9, loc='lower left')
fig.savefig(FIGDIR / 'fig6_nn_fill_map.pdf')
plt.close(fig)
print("  ✓ Fig 6")

# ===========================================================================
# FIG 7: DFO validation (12-panel)
# ===========================================================================
print("Fig 7: DFO validation...")
fig, axes = plt.subplots(3, 4, figsize=(14, 9), sharex=True)
for i, (_, row) in enumerate(dfo.iterrows()):
    ax = axes.flat[i]
    dfo_vals = row[dfo_months_cols].values.astype(float)
    dist = np.sqrt((lats - row['lat'])**2 + (lons - row['lon'])**2)
    nearest = dist.argmin()
    woa_vals = woa[nearest]
    p0_vals = p0[nearest]

    ax.plot(months_arr, dfo_vals, 'ko-', ms=4, lw=1.5, label='DFO obs')
    ax.plot(months_arr, woa_vals, 's-', ms=3, lw=1.2, color='#1f77b4', label='WOA23')
    ax.plot(months_arr, p0_vals, 'd--', ms=3, lw=1.2, color='#d62728', label='Parametric')

    rmse_w = dfo_rmse_woa[row['station']]
    rmse_p = dfo_rmse_param[row['station']]
    ax.set_title(f"{row['station']}\nWOA: {rmse_w:.1f} / Par: {rmse_p:.1f}", fontsize=7)
    ax.set_xticks(range(0, 12, 2))
    ax.set_xticklabels([MONTH_LABELS[m] for m in range(0, 12, 2)], fontsize=6)
    if i % 4 == 0:
        ax.set_ylabel('Salinity (psu)', fontsize=8)
    if i == 0:
        ax.legend(fontsize=6, loc='lower left')

fig.suptitle('DFO Lighthouse Stations: Observed vs WOA23 vs Parametric', fontsize=13, y=1.01)
fig.tight_layout()
fig.savefig(FIGDIR / 'fig7_dfo_validation.pdf')
plt.close(fig)
print("  ✓ Fig 7")

# ===========================================================================
# FIG 8: Suppression comparison — side-by-side maps
# ===========================================================================
print("Fig 8: suppression comparison...")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7))

vmax_s = 50  # percent suppression
for ax, supp, title_str in [
    (ax1, suppression_p0_june, 'Parametric Baseline'),
    (ax2, suppression_woa_june, 'WOA23 Observational')]:
    sc = ax.scatter(lons, lats, c=supp, cmap='YlOrRd', vmin=0, vmax=vmax_s,
                    s=14, alpha=0.8, edgecolors='none', zorder=3)
    ax.set_xlim(-180, -112); ax.set_ylim(25, 62)
    _add_coastline(ax, theme='light')
    ax.set_xlabel('Longitude'); ax.set_ylabel('Latitude')
    ax.set_title(f'June Salinity Suppression\n{title_str}')
    cb = fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.02)
    cb.set_label('Transmission Suppression (%)')
    # Annotate stats
    ax.text(0.05, 0.05,
            f'Mean: {supp.mean():.1f}%\nMax: {supp.max():.1f}%\nn>1%: {int((supp > 1).sum())}',
            transform=ax.transAxes, fontsize=8, va='bottom',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

fig.suptitle('June Disease Suppression: Parametric vs WOA23-Informed', fontsize=14, y=1.02)
fig.tight_layout()
fig.savefig(FIGDIR / 'fig8_suppression_comparison.pdf')
plt.close(fig)
print("  ✓ Fig 8")

# ===========================================================================
# FIG 9: Latitude gradient (3-panel)
# ===========================================================================
print("Fig 9: latitude gradient...")
fig, axes = plt.subplots(1, 3, figsize=(14, 5), sharey=True)

# Panel a: WOA23 June salinity vs latitude
ax = axes[0]
for g in GROUP_ORDER:
    mask = groups == g
    ax.scatter(lats[mask], woa_june[mask], s=10, alpha=0.5,
               color=GROUP_COLORS[g], label=g, edgecolors='none')
ax.set_xlabel('Latitude (°N)'); ax.set_ylabel('June Salinity (psu)')
ax.set_title('(a) WOA23 June Salinity')
ax.legend(fontsize=6, loc='lower left')
ax.axhline(28, color='red', ls=':', lw=1, alpha=0.5)
ax.text(27, 27.5, 'Suppression threshold (28 psu)', fontsize=7, color='red', alpha=0.7)

# Panel b: Parametric June salinity vs latitude
ax = axes[1]
for g in GROUP_ORDER:
    mask = groups == g
    ax.scatter(lats[mask], p0_june[mask], s=10, alpha=0.5,
               color=GROUP_COLORS[g], edgecolors='none')
ax.set_xlabel('Latitude (°N)')
ax.set_title('(b) Parametric June Salinity')
ax.axhline(28, color='red', ls=':', lw=1, alpha=0.5)
# Show the parametric formula line
lat_line = np.linspace(26, 62, 100)
param_line = 31.32 + 0.054 * (lat_line - 50)
ax.plot(lat_line, param_line, 'k-', lw=1.5, alpha=0.5, label='Formula')
ax.legend(fontsize=7, loc='lower left')

# Panel c: Difference (parametric - WOA23) vs latitude
ax = axes[2]
diff_june = p0_june - woa_june
for g in GROUP_ORDER:
    mask = groups == g
    ax.scatter(lats[mask], diff_june[mask], s=10, alpha=0.5,
               color=GROUP_COLORS[g], edgecolors='none')
ax.axhline(0, color='black', ls='-', lw=0.5, alpha=0.3)
ax.set_xlabel('Latitude (°N)')
ax.set_title('(c) Difference (Param − WOA23)')
ax.set_ylabel('Δ Salinity (psu)')

fig.suptitle('Latitude–Salinity Gradient in June', fontsize=14, y=1.02)
fig.tight_layout()
fig.savefig(FIGDIR / 'fig9_latitude_gradient.pdf')
plt.close(fig)
print("  ✓ Fig 9")

# ===========================================================================
# SUPPLEMENTARY: Seasonal range map + fjord depth scatter
# ===========================================================================
print("Supplementary: seasonal range map...")
fig, ax = plt.subplots(figsize=(8, 8))
sc = ax.scatter(lons, lats, c=woa_range, cmap='viridis', vmin=0, vmax=8,
                s=14, alpha=0.8, edgecolors='none', zorder=3)
ax.set_xlim(-180, -112); ax.set_ylim(25, 62)
_add_coastline(ax, theme='light')
cb = fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.02)
cb.set_label('WOA23 Seasonal Range (max − min, psu)')
ax.set_xlabel('Longitude'); ax.set_ylabel('Latitude')
ax.set_title('WOA23 Seasonal Salinity Range Across 896 Sites')
ax.text(0.05, 0.05, f'Mean range: {woa_range.mean():.2f} psu\nMax range: {woa_range.max():.2f} psu',
        transform=ax.transAxes, fontsize=9,
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
fig.savefig(FIGDIR / 'figS1_seasonal_range_map.pdf')
plt.close(fig)
print("  ✓ Fig S1")

print("\n✓ All figures generated in", FIGDIR)
print("✓ data_summary.json updated")
