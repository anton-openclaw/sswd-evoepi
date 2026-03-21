#!/usr/bin/env python3
"""Phase 1: Generate all figures and data_summary.json for the salinity-v2 report.

Produces 8+ publication-quality PDF figures and a JSON file with every
number that appears in the report text.
"""

import sys
import os
import json
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker

# Add project root
sys.path.insert(0, '/home/starbot/.openclaw/workspace/sswd-evoepi')

from sswd_evoepi.spatial import NodeDefinition
from sswd_evoepi.salinity import (
    compute_salinity_array, ocean_baseline, freshwater_melt_pulse,
    latitude_melt_factor, DAYS_PER_YEAR, _PEAK_DAY, _monthly_to_daily,
)
from sswd_evoepi.viz.salinity import (
    plot_suppression_map, plot_suppression_regional_zoom,
    plot_salinity_heatmap, plot_depression_heatmap,
    plot_regional_salinity_profiles, plot_fw_strength_sensitivity,
    plot_depth_exp_comparison, plot_suppression_monthly_panels,
    plot_fjord_depth_by_region, plot_salinity_regional_profiles,
    REGION_COLORS, MONTH_LABELS, MONTH_STARTS,
    _sal_mod, _sal_mod_array,
    FRESH_WATER, SALT_WATER, SUPPRESSION, NEUTRAL, FJORD_COLOR, OPEN_COLOR,
)
from sswd_evoepi.viz.style import (
    themed_figure, save_figure, LIGHT_TEXT, LIGHT_GRID, LIGHT_PANEL,
)

WORKDIR = '/home/starbot/.openclaw/workspace/sswd-evoepi/reports/salinity-v2'
FIGDIR = os.path.join(WORKDIR, 'figures')
os.makedirs(FIGDIR, exist_ok=True)

# ── Load nodes ──────────────────────────────────────────────────────
def load_nodes():
    with open('/home/starbot/.openclaw/workspace/sswd-evoepi/data/nodes/all_sites.json') as f:
        sites = json.load(f)
    enc = {}
    with open('/home/starbot/.openclaw/workspace/sswd-evoepi/data/nodes/site_enclosedness.csv') as f:
        for row in csv.DictReader(f):
            enc[row['name']] = float(row.get('fjord_depth_norm', 0.0))
    nodes = []
    for i, s in enumerate(sites):
        nd = NodeDefinition(
            node_id=i, name=s['name'],
            lat=float(s.get('lat') or s.get('latitude')),
            lon=float(s.get('lon') or s.get('longitude')),
            subregion=s.get('region', ''),
            habitat_area=1000.0, carrying_capacity=5000,
            fjord_depth_norm=enc.get(s['name'], 0.0),
        )
        nodes.append(nd)
    return nodes

# ── Load validation data ────────────────────────────────────────────
def load_npz():
    d = np.load('/home/starbot/.openclaw/workspace/salinity_validation/woa23_site_monthly.npz',
                allow_pickle=True)
    return {k: d[k] for k in d.keys()}

def load_dfo():
    rows = []
    with open('/home/starbot/.openclaw/workspace/salinity_validation/dfo_monthly_climatology.csv') as f:
        for row in csv.DictReader(f):
            row['lat'] = float(row['lat'])
            row['lon'] = float(row['lon'])
            for m in MONTH_LABELS:
                ml = m.lower()
                row[ml] = float(row[ml]) if ml in row else float(row[m])
            rows.append(row)
    return rows

# ── Compute stats for data_summary ──────────────────────────────────
def compute_stats(npz, nodes, dfo):
    summary = {}

    names = npz['names']
    woa = npz['woa_monthly']        # (896, 12) — WOA23 values at each site
    param0 = npz['param_monthly_0']  # (896, 12) — parametric baseline (fw=0)
    param15 = npz['param_monthly_15']
    nn_dist = npz['nn_distance']
    woa_direct = npz['woa_direct']   # (896, 12) — NaN where NN-filled
    regions = npz['regions']
    lats = npz['lats']
    fd_norm = npz['fjord_depth_norm']

    N = len(names)
    summary['n_sites'] = int(N)

    # Direct vs NN fill
    # woa_direct has actual values where direct, NaN where NN-filled
    has_direct = ~np.isnan(woa_direct[:, 0])
    n_direct = int(has_direct.sum())
    n_nn = int(N - n_direct)
    summary['n_direct'] = n_direct
    summary['n_nn_fill'] = n_nn
    summary['max_nn_distance_deg'] = float(np.nanmax(nn_dist))
    summary['max_nn_distance_km'] = round(float(np.nanmax(nn_dist)) * 111.0, 1)

    # Annual means
    woa_annual = woa.mean(axis=1)
    param_annual = param0.mean(axis=1)  # parametric baseline (constant, so all same per site)

    # Parametric vs WOA23 comparison
    diff = woa_annual - param_annual
    summary['mean_bias_psu'] = round(float(np.mean(diff)), 2)
    summary['rmse_psu'] = round(float(np.sqrt(np.mean(diff**2))), 2)
    summary['max_error_psu'] = round(float(np.max(np.abs(diff))), 2)

    # Seasonal range
    woa_range = woa.max(axis=1) - woa.min(axis=1)
    param_range = param0.max(axis=1) - param0.min(axis=1)
    summary['woa_seasonal_range_mean'] = round(float(woa_range.mean()), 2)
    summary['param_seasonal_range_mean'] = round(float(param_range.mean()), 2)

    # Month of minimum
    woa_min_month = np.argmin(woa, axis=1)  # 0-indexed month
    month_counts = {}
    for m in range(12):
        month_counts[MONTH_LABELS[m]] = int((woa_min_month == m).sum())
    summary['month_of_minimum'] = month_counts

    # Top months
    sorted_months = sorted(month_counts.items(), key=lambda x: -x[1])
    summary['top_min_months'] = {m: c for m, c in sorted_months[:5]}

    # Regional stats
    unique_regions = sorted(set(regions))
    regional_bias = {}
    regional_rmse = {}
    regional_n = {}
    for r in unique_regions:
        mask = regions == r
        rd = diff[mask]
        regional_bias[r] = round(float(rd.mean()), 2)
        regional_rmse[r] = round(float(np.sqrt(np.mean(rd**2))), 2)
        regional_n[r] = int(mask.sum())
    summary['regional_bias'] = regional_bias
    summary['regional_rmse'] = regional_rmse
    summary['regional_n'] = regional_n

    # DFO validation
    dfo_obs = []
    dfo_pred = []
    dfo_stations = []
    for row in dfo:
        lat = row['lat']
        lon = row['lon']
        station = row['station']
        dfo_stations.append(station)
        for mi, m in enumerate(MONTH_LABELS):
            obs_val = row[m.lower()]
            # Get WOA23 prediction at this lat/lon
            # Use the salinity module to compute
            syn_node = NodeDefinition(
                node_id=0, name=station, lat=lat, lon=lon,
                subregion='DFO', habitat_area=1000, carrying_capacity=100,
                fjord_depth_norm=0.0,
            )
            sal = compute_salinity_array([syn_node], fw_strength=0.0)
            start = MONTH_STARTS[mi]
            end = MONTH_STARTS[mi + 1] if mi < 11 else 365
            pred_val = float(sal[0, start:end].mean())
            dfo_obs.append(obs_val)
            dfo_pred.append(pred_val)

    dfo_obs = np.array(dfo_obs)
    dfo_pred = np.array(dfo_pred)
    dfo_resid = dfo_pred - dfo_obs
    summary['dfo_n_stations'] = len(dfo)
    summary['dfo_n_points'] = len(dfo_obs)
    summary['dfo_mean_bias'] = round(float(dfo_resid.mean()), 2)
    summary['dfo_rmse'] = round(float(np.sqrt(np.mean(dfo_resid**2))), 2)
    summary['dfo_max_error'] = round(float(np.max(np.abs(dfo_resid))), 2)
    summary['dfo_correlation'] = round(float(np.corrcoef(dfo_obs, dfo_pred)[0, 1]), 3)

    # Also compute parametric baseline DFO stats
    dfo_param_pred = []
    for row in dfo:
        lat = row['lat']
        baseline = ocean_baseline(lat)
        for mi in range(12):
            dfo_param_pred.append(baseline)
    dfo_param_pred = np.array(dfo_param_pred)
    dfo_param_resid = dfo_param_pred - dfo_obs
    summary['dfo_param_mean_bias'] = round(float(dfo_param_resid.mean()), 2)
    summary['dfo_param_rmse'] = round(float(np.sqrt(np.mean(dfo_param_resid**2))), 2)
    summary['dfo_param_correlation'] = round(float(np.corrcoef(dfo_obs, dfo_param_pred)[0, 1]), 3)

    # Store DFO arrays for figure generation
    summary['_dfo_obs'] = dfo_obs.tolist()
    summary['_dfo_pred'] = dfo_pred.tolist()
    summary['_dfo_param_pred'] = dfo_param_pred.tolist()
    summary['_dfo_stations'] = dfo_stations

    # Latitude statistics
    summary['lat_range'] = [round(float(lats.min()), 2), round(float(lats.max()), 2)]
    summary['n_regions'] = len(unique_regions)
    summary['regions_list'] = list(unique_regions)

    # WOA23 extremes
    summary['woa_min_salinity'] = round(float(woa.min()), 2)
    summary['woa_max_salinity'] = round(float(woa.max()), 2)
    summary['woa_annual_mean_min'] = round(float(woa_annual.min()), 2)
    summary['woa_annual_mean_max'] = round(float(woa_annual.max()), 2)

    return summary


# ── Figure 1: Regional zoom suppression map ─────────────────────────
def fig_regional_zoom(nodes):
    print("  Fig 1: Regional zoom suppression map...")
    fig = plot_suppression_regional_zoom(
        nodes, fw_strength=15.0, fw_depth_exp=1.0,
        month=5, theme='light',
        save_path=os.path.join(FIGDIR, 'fig_regional_zoom.pdf'),
    )
    plt.close(fig)

# ── Figure 2: Full coast suppression map ────────────────────────────
def fig_full_coast(nodes):
    print("  Fig 2: Full coast suppression map...")
    fig = plot_suppression_map(
        nodes, fw_strength=15.0, fw_depth_exp=1.0,
        month=5, theme='light',
        save_path=os.path.join(FIGDIR, 'fig_full_coast.pdf'),
    )
    plt.close(fig)

# ── Figure 3: WOA23 vs parametric scatter ───────────────────────────
def fig_scatter_comparison(npz):
    print("  Fig 3: WOA23 vs parametric scatter...")
    woa = npz['woa_monthly']
    param = npz['param_monthly_0']
    regions = npz['regions']

    woa_annual = woa.mean(axis=1)
    param_annual = param.mean(axis=1)

    unique_regions = sorted(set(regions))
    # Group into macro-regions for cleaner colours
    macro_map = {}
    for r in unique_regions:
        if r.startswith('AK'):
            macro_map[r] = 'Alaska'
        elif r.startswith('BC'):
            macro_map[r] = 'British Columbia'
        elif r.startswith('SS') or r == 'JDF':
            macro_map[r] = 'Salish Sea'
        elif r in ('WA-O', 'OR'):
            macro_map[r] = 'WA/OR'
        elif r.startswith('CA') or r == 'BJ':
            macro_map[r] = 'CA/Baja'
        else:
            macro_map[r] = 'Other'

    macro_colors = {
        'Alaska': '#1f77b4',
        'British Columbia': '#ff7f0e',
        'Salish Sea': '#2ca02c',
        'WA/OR': '#d62728',
        'CA/Baja': '#9467bd',
        'Other': '#7f7f7f',
    }

    fig, ax = themed_figure(theme='light', figsize=(9, 8))

    for macro_name, color in macro_colors.items():
        mask = np.array([macro_map.get(r, 'Other') == macro_name for r in regions])
        if mask.any():
            ax.scatter(param_annual[mask], woa_annual[mask],
                      c=color, s=20, alpha=0.6, label=macro_name, zorder=3)

    # 1:1 line
    lo = min(param_annual.min(), woa_annual.min()) - 1
    hi = max(param_annual.max(), woa_annual.max()) + 1
    ax.plot([lo, hi], [lo, hi], 'k--', linewidth=1, alpha=0.5, label='1:1 line')

    # Stats annotation
    diff = woa_annual - param_annual
    rmse = np.sqrt(np.mean(diff**2))
    bias = np.mean(diff)
    ax.text(0.03, 0.97,
            f'N = {len(woa_annual)}\nBias = {bias:.2f} psu\nRMSE = {rmse:.2f} psu',
            transform=ax.transAxes, fontsize=10, va='top',
            bbox=dict(boxstyle='round', facecolor='white', edgecolor='#cccccc', alpha=0.9))

    ax.set_xlabel('Parametric baseline salinity (psu)', fontsize=12)
    ax.set_ylabel('WOA23 annual mean salinity (psu)', fontsize=12)
    ax.set_title('WOA23 vs Parametric Baseline — Annual Mean', fontsize=14, fontweight='bold')
    ax.legend(fontsize=9, loc='lower right', framealpha=0.9)
    ax.set_aspect('equal')
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.grid(True, alpha=0.3)

    save_figure(fig, os.path.join(FIGDIR, 'fig_scatter_comparison.pdf'))
    plt.close(fig)

# ── Figure 4: Seasonal timing histogram ─────────────────────────────
def fig_timing_histogram(npz):
    print("  Fig 4: Seasonal timing histogram...")
    woa = npz['woa_monthly']
    regions = npz['regions']

    min_month = np.argmin(woa, axis=1)  # 0-indexed

    # Macro regions
    macro_labels = []
    for r in regions:
        if r.startswith('AK'):
            macro_labels.append('Alaska')
        elif r.startswith('BC'):
            macro_labels.append('BC')
        elif r.startswith('SS') or r == 'JDF':
            macro_labels.append('Salish Sea')
        elif r in ('WA-O', 'OR'):
            macro_labels.append('WA/OR')
        else:
            macro_labels.append('CA/Baja')
    macro_labels = np.array(macro_labels)

    macro_order = ['Alaska', 'BC', 'Salish Sea', 'WA/OR', 'CA/Baja']
    macro_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']

    fig, ax = themed_figure(theme='light', figsize=(10, 6))

    bottoms = np.zeros(12)
    for macro, color in zip(macro_order, macro_colors):
        mask = macro_labels == macro
        counts = np.array([(min_month[mask] == m).sum() for m in range(12)])
        ax.bar(range(12), counts, bottom=bottoms, color=color, label=macro,
               edgecolor='white', linewidth=0.5)
        bottoms += counts

    ax.set_xticks(range(12))
    ax.set_xticklabels(MONTH_LABELS, fontsize=10)
    ax.set_xlabel('Month of Minimum WOA23 Salinity', fontsize=12)
    ax.set_ylabel('Number of Sites', fontsize=12)
    ax.set_title('Seasonal Timing of Salinity Minimum (WOA23)', fontsize=14, fontweight='bold')
    ax.legend(fontsize=10, loc='upper right', framealpha=0.9)
    ax.grid(axis='y', alpha=0.3)

    save_figure(fig, os.path.join(FIGDIR, 'fig_timing_histogram.pdf'))
    plt.close(fig)

# ── Figure 5: Example site profiles ─────────────────────────────────
def fig_site_profiles(nodes, npz):
    print("  Fig 5: Example site profiles...")
    woa = npz['woa_monthly']
    regions = npz['regions']
    names = npz['names']
    fd = npz['fjord_depth_norm']
    lats = npz['lats']

    # Pick 4 contrasting sites
    # AK fjord: high lat, high fd
    # BC coast: mid lat, low fd
    # Salish Sea: SS region
    # CA open: low lat, low fd
    picks = {}
    for i in range(len(names)):
        r = str(regions[i])
        lat = lats[i]
        fdn = fd[i]
        if 'AK fjord' not in picks and r.startswith('AK') and fdn > 0.6 and lat > 56:
            picks['AK fjord'] = i
        elif 'BC coast' not in picks and r.startswith('BC') and fdn < 0.15 and 49 < lat < 53:
            picks['BC coast'] = i
        elif 'Salish Sea' not in picks and (r.startswith('SS') or r == 'JDF') and fdn > 0.2:
            picks['Salish Sea'] = i
        elif 'CA open' not in picks and r.startswith('CA') and fdn < 0.15 and lat < 38:
            picks['CA open'] = i
        if len(picks) >= 4:
            break

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    days = np.arange(365)

    fig, axes = themed_figure(theme='light', nrows=2, ncols=2, figsize=(14, 10))

    for ax, (label, idx), color in zip(axes.flat, picks.items(), colors):
        monthly_woa = woa[idx]
        lat = lats[idx]
        fdn = fd[idx]
        site_name = str(names[idx])

        # WOA23 daily interpolated
        daily_woa = _monthly_to_daily(monthly_woa)

        # Parametric baseline (constant)
        param_base = ocean_baseline(lat)
        param_line = np.full(365, param_base)

        # Combined model at fw=15
        node = NodeDefinition(
            node_id=0, name=site_name, lat=lat, lon=float(npz['lons'][idx]),
            subregion=str(regions[idx]), habitat_area=1000, carrying_capacity=5000,
            fjord_depth_norm=fdn,
        )
        sal_combined = compute_salinity_array([node], fw_strength=15.0)

        # Plot
        ax.plot(days, daily_woa, color=color, linewidth=2, label='WOA23 baseline')
        ax.axhline(param_base, color='#888888', linewidth=1.5, linestyle='--', label=f'Parametric ({param_base:.1f})')
        ax.plot(days, sal_combined[0], color=color, linewidth=2, linestyle=':', label='Combined (fw=15)')

        # Shading between WOA and combined
        ax.fill_between(days, daily_woa, sal_combined[0], alpha=0.15, color=color)

        ax.set_xticks(MONTH_STARTS)
        ax.set_xticklabels(MONTH_LABELS, fontsize=8)
        ax.set_ylabel('Salinity (psu)', fontsize=10)
        ax.set_title(f'{label}: {site_name}\n({lat:.1f}°N, fd={fdn:.2f})',
                     fontsize=11, fontweight='bold')
        ax.legend(fontsize=8, loc='lower left', framealpha=0.9)
        ax.grid(True, alpha=0.2)

        # Reference lines
        ax.axhline(28, color='#cccccc', linewidth=0.8, linestyle=':')
        ax.axhline(10, color='#cccccc', linewidth=0.8, linestyle=':')

    fig.suptitle('Example Site Salinity Profiles: WOA23 vs Parametric vs Combined',
                 fontsize=14, fontweight='bold', y=1.02)
    fig.tight_layout()

    save_figure(fig, os.path.join(FIGDIR, 'fig_site_profiles.pdf'))
    plt.close(fig)

    return picks  # return for data_summary

# ── Figure 6: Regional salinity profiles ────────────────────────────
def fig_regional_profiles(nodes):
    print("  Fig 6: Regional salinity profiles...")
    fig = plot_salinity_regional_profiles(
        nodes, fw_strength=15.0, fw_depth_exp=1.0,
        theme='light',
        save_path=os.path.join(FIGDIR, 'fig_regional_profiles.pdf'),
    )
    plt.close(fig)

# ── Figure 7: DFO validation scatter ────────────────────────────────
def fig_dfo_validation(dfo, summary):
    print("  Fig 7: DFO validation scatter...")
    dfo_obs = np.array(summary['_dfo_obs'])
    dfo_pred = np.array(summary['_dfo_pred'])
    dfo_param = np.array(summary['_dfo_param_pred'])
    stations = summary['_dfo_stations']

    n_stations = len(dfo)

    fig, (ax1, ax2) = themed_figure(theme='light', nrows=1, ncols=2, figsize=(14, 7))

    # Panel 1: WOA23 vs observed
    station_colors = plt.cm.tab10(np.linspace(0, 1, n_stations))
    for si, station_row in enumerate(dfo):
        station_name = station_row['station']
        start = si * 12
        end = (si + 1) * 12
        ax1.scatter(dfo_obs[start:end], dfo_pred[start:end],
                   c=[station_colors[si]], s=40, alpha=0.7,
                   label=station_name, zorder=3)

    lo = min(dfo_obs.min(), dfo_pred.min()) - 1
    hi = max(dfo_obs.max(), dfo_pred.max()) + 1
    ax1.plot([lo, hi], [lo, hi], 'k--', linewidth=1, alpha=0.5)
    ax1.set_xlabel('DFO Observed Salinity (psu)', fontsize=12)
    ax1.set_ylabel('WOA23 Predicted Salinity (psu)', fontsize=12)
    ax1.set_title('WOA23 vs DFO Lighthouse Stations', fontsize=13, fontweight='bold')
    ax1.set_aspect('equal')
    ax1.set_xlim(lo, hi)
    ax1.set_ylim(lo, hi)

    resid = dfo_pred - dfo_obs
    ax1.text(0.03, 0.97,
             f'N = {len(dfo_obs)} ({n_stations} stations × 12 months)\n'
             f'Bias = {resid.mean():.2f} psu\n'
             f'RMSE = {np.sqrt(np.mean(resid**2)):.2f} psu\n'
             f'r = {np.corrcoef(dfo_obs, dfo_pred)[0,1]:.3f}',
             transform=ax1.transAxes, fontsize=9, va='top',
             bbox=dict(boxstyle='round', facecolor='white', edgecolor='#ccc', alpha=0.9))
    ax1.legend(fontsize=7, loc='lower right', framealpha=0.9, ncol=2)
    ax1.grid(True, alpha=0.3)

    # Panel 2: Parametric vs observed
    for si, station_row in enumerate(dfo):
        station_name = station_row['station']
        start = si * 12
        end = (si + 1) * 12
        ax2.scatter(dfo_obs[start:end], dfo_param[start:end],
                   c=[station_colors[si]], s=40, alpha=0.7,
                   label=station_name, zorder=3)

    lo2 = min(dfo_obs.min(), dfo_param.min()) - 1
    hi2 = max(dfo_obs.max(), dfo_param.max()) + 1
    ax2.plot([lo2, hi2], [lo2, hi2], 'k--', linewidth=1, alpha=0.5)
    ax2.set_xlabel('DFO Observed Salinity (psu)', fontsize=12)
    ax2.set_ylabel('Parametric Predicted Salinity (psu)', fontsize=12)
    ax2.set_title('Parametric Baseline vs DFO', fontsize=13, fontweight='bold')
    ax2.set_aspect('equal')
    ax2.set_xlim(lo2, hi2)
    ax2.set_ylim(lo2, hi2)

    resid_p = dfo_param - dfo_obs
    ax2.text(0.03, 0.97,
             f'Bias = {resid_p.mean():.2f} psu\n'
             f'RMSE = {np.sqrt(np.mean(resid_p**2)):.2f} psu\n'
             f'r = {np.corrcoef(dfo_obs, dfo_param)[0,1]:.3f}',
             transform=ax2.transAxes, fontsize=9, va='top',
             bbox=dict(boxstyle='round', facecolor='white', edgecolor='#ccc', alpha=0.9))
    ax2.legend(fontsize=7, loc='lower right', framealpha=0.9, ncol=2)
    ax2.grid(True, alpha=0.3)

    fig.suptitle('Model Validation Against DFO Lighthouse Observations',
                 fontsize=14, fontweight='bold', y=1.02)
    fig.tight_layout()
    save_figure(fig, os.path.join(FIGDIR, 'fig_dfo_validation.pdf'))
    plt.close(fig)

# ── Figure 8: fw_strength sensitivity ───────────────────────────────
def fig_fw_sensitivity(nodes):
    print("  Fig 8: fw_strength sensitivity...")
    fig = plot_fw_strength_sensitivity(
        nodes, fw_depth_exp=1.0, theme='light',
        save_path=os.path.join(FIGDIR, 'fig_fw_sensitivity.pdf'),
    )
    plt.close(fig)

# ═══════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    print("Loading data...")
    nodes = load_nodes()
    npz = load_npz()
    dfo = load_dfo()

    print(f"Loaded {len(nodes)} nodes, {len(dfo)} DFO stations")

    print("Computing statistics...")
    summary = compute_stats(npz, nodes, dfo)

    print("Generating figures...")
    fig_regional_zoom(nodes)
    fig_full_coast(nodes)
    fig_scatter_comparison(npz)
    fig_timing_histogram(npz)
    site_picks = fig_site_profiles(nodes, npz)
    fig_regional_profiles(nodes)
    fig_dfo_validation(dfo, summary)
    fig_fw_sensitivity(nodes)

    # Record site picks
    summary['example_sites'] = {label: str(npz['names'][idx]) for label, idx in site_picks.items()}

    # Remove internal arrays before saving JSON
    for key in list(summary.keys()):
        if key.startswith('_'):
            del summary[key]

    # Save data_summary.json
    json_path = os.path.join(WORKDIR, 'data_summary.json')
    with open(json_path, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"\nSaved data_summary.json to {json_path}")

    print("\nAll figures generated in:", FIGDIR)
    print("Done!")
