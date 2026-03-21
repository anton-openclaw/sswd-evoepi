#!/usr/bin/env python3
"""Generate all salinity mechanism figures and data_summary.json.

Uses the full 896-node network from load_sites() + build_node_defs().
"""

import sys, json, os, csv
import numpy as np

# Setup paths
PROJECT_ROOT = '/home/starbot/.openclaw/workspace/sswd-evoepi'
sys.path.insert(0, PROJECT_ROOT)
os.chdir(PROJECT_ROOT)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ── STEP 1: Load 896-node network ──────────────────────────────────
print("STEP 1: Loading 896-node network...")
from experiments.calibration_runner import load_sites, build_node_defs
sites = load_sites()
node_defs = build_node_defs(sites, K=5000)
print(f"  Loaded {len(node_defs)} nodes")

from sswd_evoepi.salinity import (
    compute_salinity_array, ocean_baseline, latitude_melt_factor,
    freshwater_melt_pulse, _PEAK_DAY, DAYS_PER_YEAR,
)

# ── Compute salinity for multiple fw_strength values ───────────────
print("  Computing salinity arrays for fw_strength sweep...")
fw_values = [0, 5, 8, 10, 12, 15, 18, 20, 25]
sal_arrays = {}
for fw in fw_values:
    sal_arrays[fw] = compute_salinity_array(node_defs, fw_strength=float(fw))
    print(f"    fw={fw}: shape={sal_arrays[fw].shape}, "
          f"min={sal_arrays[fw].min():.2f}, max={sal_arrays[fw].max():.2f}")

# Salinity modifier helper
def sal_mod_scalar(s, s_min=10.0, s_full=28.0, eta=2.0):
    if s >= s_full: return 1.0
    if s <= s_min: return 0.0
    return ((s - s_min) / (s_full - s_min)) ** eta

def sal_mod_array(arr, s_min=10.0, s_full=28.0, eta=2.0):
    x = np.clip((arr - s_min) / (s_full - s_min), 0.0, 1.0)
    return x ** eta


# ── STEP 2: Load DFO data ──────────────────────────────────────────
print("\nSTEP 2: Loading DFO lighthouse climatology...")
dfo_path = 'salinity_validation/dfo_monthly_climatology.csv'
dfo_stations = []
month_names = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
with open(dfo_path) as f:
    reader = csv.DictReader(f)
    for row in reader:
        months_data = [float(row[m]) for m in month_names]
        dfo_stations.append({
            'station': row['station'],
            'latitude': float(row['latitude']),
            'longitude': float(row['longitude']),
            'monthly_salinity': months_data,
        })
print(f"  Loaded {len(dfo_stations)} DFO stations")


# ── STEP 3: Compute summary statistics ─────────────────────────────
print("\nSTEP 3: Computing summary statistics...")

# Per-region stats
region_stats = {}
for i, nd in enumerate(node_defs):
    r = nd.subregion
    if r not in region_stats:
        region_stats[r] = {
            'lats': [], 'fd_norms': [], 'indices': [],
        }
    region_stats[r]['lats'].append(nd.lat)
    region_stats[r]['fd_norms'].append(nd.fjord_depth_norm)
    region_stats[r]['indices'].append(i)

# fw=15 reference
sal_15 = sal_arrays[15]
sal_0 = sal_arrays[0]

# June suppression (month 5, days 151-181)
june_start, june_end = 151, 181
june_sal_15 = sal_15[:, june_start:june_end].mean(axis=1)
june_sal_0 = sal_0[:, june_start:june_end].mean(axis=1)
june_depression = june_sal_0 - june_sal_15
june_sal_mod = sal_mod_array(june_sal_15)
june_suppression = (1.0 - june_sal_mod) * 100

# Peak-day metrics
peak_sal_15 = sal_15[:, _PEAK_DAY]
peak_sal_0 = sal_0[:, _PEAK_DAY]
peak_depression = peak_sal_0 - peak_sal_15
peak_sal_mod = sal_mod_array(peak_sal_15)
peak_suppression = (1.0 - peak_sal_mod) * 100

region_summary = {}
for r, data in sorted(region_stats.items()):
    idxs = data['indices']
    fd_arr = np.array(data['fd_norms'])
    lat_arr = np.array(data['lats'])
    region_summary[r] = {
        'n_nodes': len(idxs),
        'mean_latitude': round(float(lat_arr.mean()), 2),
        'mean_fjord_depth_norm': round(float(fd_arr.mean()), 4),
        'n_nodes_fd_gt_0.5': int((fd_arr > 0.5).sum()),
        'fw15_june_mean_depression_psu': round(float(june_depression[idxs].mean()), 3),
        'fw15_june_mean_suppression_pct': round(float(june_suppression[idxs].mean()), 2),
        'fw15_peak_mean_depression_psu': round(float(peak_depression[idxs].mean()), 3),
        'fw15_peak_mean_sal_mod': round(float(peak_sal_mod[idxs].mean()), 4),
        'fw15_peak_mean_suppression_pct': round(float(peak_suppression[idxs].mean()), 2),
    }

# Model predictions at DFO latitudes (open-coast fd=0)
from sswd_evoepi.spatial import NodeDefinition
dfo_model_comparison = []
for st in dfo_stations:
    lat = st['latitude']
    syn = NodeDefinition(
        node_id=0, name=st['station'], lat=lat, lon=st['longitude'],
        subregion='DFO', habitat_area=1000, carrying_capacity=100,
        fjord_depth_norm=0.0,
    )
    model_sal = compute_salinity_array([syn], fw_strength=15.0)
    # Monthly means
    month_starts = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]
    model_monthly = []
    for m in range(12):
        s = month_starts[m]
        e = month_starts[m+1] if m < 11 else 365
        model_monthly.append(round(float(model_sal[0, s:e].mean()), 2))
    
    residuals = [round(model_monthly[m] - st['monthly_salinity'][m], 2) for m in range(12)]
    dfo_model_comparison.append({
        'station': st['station'],
        'latitude': lat,
        'observed_monthly': st['monthly_salinity'],
        'model_monthly': model_monthly,
        'residuals': residuals,
        'mean_abs_residual': round(float(np.mean(np.abs(residuals))), 3),
    })

# Asymmetry metrics
ak_mask = np.array([nd.lat > 55 for nd in node_defs])
ca_mask = np.array([nd.lat < 40 for nd in node_defs])
ak_mean_sup = float(peak_suppression[ak_mask].mean()) if ak_mask.any() else 0
ca_mean_sup = float(peak_suppression[ca_mask].mean()) if ca_mask.any() else 0
asymmetry_ratio = ak_mean_sup / ca_mean_sup if ca_mean_sup > 0.01 else float('inf')

asymmetry_metrics = {
    'ak_n_nodes': int(ak_mask.sum()),
    'ca_n_nodes': int(ca_mask.sum()),
    'ak_mean_peak_suppression_pct': round(ak_mean_sup, 2),
    'ca_mean_peak_suppression_pct': round(ca_mean_sup, 2),
    'asymmetry_ratio': round(asymmetry_ratio, 1) if asymmetry_ratio != float('inf') else 'inf',
    'ak_mean_peak_depression_psu': round(float(peak_depression[ak_mask].mean()), 3) if ak_mask.any() else 0,
    'ca_mean_peak_depression_psu': round(float(peak_depression[ca_mask].mean()), 3) if ca_mask.any() else 0,
    'ak_mean_fjord_depth_norm': round(float(np.array([nd.fjord_depth_norm for nd in node_defs])[ak_mask].mean()), 4) if ak_mask.any() else 0,
    'ca_mean_fjord_depth_norm': round(float(np.array([nd.fjord_depth_norm for nd in node_defs])[ca_mask].mean()), 4) if ca_mask.any() else 0,
}

print(f"  Computed stats for {len(region_summary)} regions")
print(f"  AK mean suppression: {ak_mean_sup:.1f}%, CA: {ca_mean_sup:.1f}%, ratio: {asymmetry_ratio:.1f}x")


# ── STEP 4: Generate all figures ────────────────────────────────────
print("\nSTEP 4: Generating figures...")

from sswd_evoepi.viz.salinity import (
    plot_mechanism_components,
    plot_sal_mod_transfer,
    plot_fjord_depth_by_region,
    plot_regional_salinity_profiles,
    plot_latitude_asymmetry,
    plot_fw_strength_sensitivity,
    plot_regional_suppression_bars,
    plot_suppression_monthly_panels,
    plot_dfo_validation,
    plot_depression_heatmap,
    plot_depth_exp_comparison,
)

figdir = '/home/starbot/.openclaw/workspace/sswd-evoepi/reports/salinity/figures/'

# 1. Mechanism components
print("  1/11 Mechanism components...")
fig = plot_mechanism_components(save_path=f'{figdir}fig_mechanism.pdf')
plt.close('all')

# 2. sal_mod transfer function with annotated sites
print("  2/11 sal_mod transfer function...")
# Find example sites: AK-FN fjord and CA-N coast
ak_fn_fjord_idx = None
ca_n_coast_idx = None
for i, nd in enumerate(node_defs):
    if nd.subregion == 'AK-FN' and nd.fjord_depth_norm > 0.6 and ak_fn_fjord_idx is None:
        ak_fn_fjord_idx = i
    if nd.subregion == 'CA-N' and nd.fjord_depth_norm < 0.1 and ca_n_coast_idx is None:
        ca_n_coast_idx = i
    if ak_fn_fjord_idx is not None and ca_n_coast_idx is not None:
        break

annotate = {}
if ak_fn_fjord_idx is not None:
    s = float(sal_15[ak_fn_fjord_idx, _PEAK_DAY])
    annotate[f'AK-FN fjord ({node_defs[ak_fn_fjord_idx].lat:.0f}°N)'] = s
if ca_n_coast_idx is not None:
    s = float(sal_15[ca_n_coast_idx, _PEAK_DAY])
    annotate[f'CA-N coast ({node_defs[ca_n_coast_idx].lat:.0f}°N)'] = s

# Add a mid-range BC example
bc_mid_idx = None
for i, nd in enumerate(node_defs):
    if nd.subregion == 'BC-N' and 0.3 < nd.fjord_depth_norm < 0.6:
        bc_mid_idx = i
        break
if bc_mid_idx is not None:
    s = float(sal_15[bc_mid_idx, _PEAK_DAY])
    annotate[f'BC-N ({node_defs[bc_mid_idx].lat:.0f}°N)'] = s

fig = plot_sal_mod_transfer(annotate_sites=annotate if annotate else None,
                            save_path=f'{figdir}fig_sal_mod.pdf')
plt.close('all')

# 3. Fjord depth by region
print("  3/11 Fjord depth by region...")
fig = plot_fjord_depth_by_region(node_defs, save_path=f'{figdir}fig_fjord_depth.pdf')
plt.close('all')

# 4. Regional salinity profiles (fw=15)
print("  4/11 Regional salinity profiles...")
fig = plot_regional_salinity_profiles(node_defs, fw_strength=15.0,
                                      save_path=f'{figdir}fig_regional_profiles.pdf')
plt.close('all')

# 5. Latitude asymmetry (fw=15)
print("  5/11 Latitude asymmetry...")
fig = plot_latitude_asymmetry(node_defs, fw_strength=15.0,
                              save_path=f'{figdir}fig_asymmetry.pdf')
plt.close('all')

# 6. fw_strength sensitivity
print("  6/11 fw_strength sensitivity...")
fig = plot_fw_strength_sensitivity(node_defs, save_path=f'{figdir}fig_sensitivity.pdf')
plt.close('all')

# 7. Regional suppression bars (fw=15, June)
print("  7/11 Regional suppression bars...")
fig = plot_regional_suppression_bars(node_defs, fw_strength=15.0,
                                     save_path=f'{figdir}fig_suppression_bars.pdf')
plt.close('all')

# 8. Monthly suppression panels
print("  8/11 Monthly suppression panels...")
fig = plot_suppression_monthly_panels(node_defs, fw_strength=15.0,
                                      save_path=f'{figdir}fig_monthly_panels.pdf')
plt.close('all')

# 9. DFO validation
print("  9/11 DFO validation...")
fig = plot_dfo_validation(dfo_path, save_path=f'{figdir}fig_dfo_validation.pdf')
plt.close('all')

# 10. Depression heatmap — top 50 most depressed
print("  10/11 Depression heatmap (top 50)...")
max_dep = (sal_0 - sal_15).max(axis=1)
top50_idx = np.argsort(max_dep)[-50:]
top50_nodes = [node_defs[i] for i in top50_idx]
fig = plot_depression_heatmap(top50_nodes, fw_strength=15.0,
                              save_path=f'{figdir}fig_depression_heatmap.pdf')
plt.close('all')

# 11. Depth exponent comparison
print("  11/11 Depth exponent comparison...")
fig = plot_depth_exp_comparison(node_defs, fw_strength=15.0,
                                save_path=f'{figdir}fig_depth_exp.pdf')
plt.close('all')


# ── STEP 5: Write data_summary.json ────────────────────────────────
print("\nSTEP 5: Writing data_summary.json...")

# Top-50 most depressed sites detail
top50_details = []
for i in np.argsort(max_dep)[-50:][::-1]:  # most depressed first
    nd = node_defs[i]
    top50_details.append({
        'name': nd.name,
        'subregion': nd.subregion,
        'latitude': round(nd.lat, 3),
        'longitude': round(nd.lon, 3),
        'fjord_depth_norm': round(nd.fjord_depth_norm, 4),
        'peak_depression_psu': round(float(max_dep[i]), 3),
        'june_sal_mod': round(float(sal_mod_scalar(float(sal_15[i, _PEAK_DAY]))), 4),
        'peak_salinity_psu': round(float(sal_15[i, _PEAK_DAY]), 2),
    })

# fw_strength sweep summary
fw_sweep = {}
for fw in fw_values:
    sal = sal_arrays[fw]
    pk_sal = sal[:, _PEAK_DAY]
    pk_mod = sal_mod_array(pk_sal)
    pk_sup = (1.0 - pk_mod) * 100
    fw_sweep[str(fw)] = {
        'global_mean_peak_salinity': round(float(pk_sal.mean()), 2),
        'global_min_peak_salinity': round(float(pk_sal.min()), 2),
        'global_mean_peak_suppression_pct': round(float(pk_sup.mean()), 2),
        'global_max_peak_suppression_pct': round(float(pk_sup.max()), 2),
        'n_nodes_sup_gt_10pct': int((pk_sup > 10).sum()),
        'n_nodes_sup_gt_50pct': int((pk_sup > 50).sum()),
    }

data_summary = {
    'n_total_nodes': len(node_defs),
    'regions': region_summary,
    'asymmetry': asymmetry_metrics,
    'fw_strength_sweep': fw_sweep,
    'dfo_validation': dfo_model_comparison,
    'top50_most_depressed': top50_details,
    'model_parameters': {
        'ocean_baseline_formula': 'S_ocean(lat) = 31.32 + 0.054 * (lat - 50.0)',
        'melt_pulse': 'cos(2π(day-166)/365), clipped to [0,1]',
        'latitude_melt_factor': 'clip((lat-35)/25, 0, 1)',
        'peak_day': 166,
        's_min': 10.0,
        's_full': 28.0,
        'eta': 2.0,
        's_floor': 5.0,
    },
}

summary_path = '/home/starbot/.openclaw/workspace/sswd-evoepi/reports/salinity/data_summary.json'
with open(summary_path, 'w') as f:
    json.dump(data_summary, f, indent=2, default=str)
print(f"  Wrote {summary_path}")


# ── STEP 6: Verify ─────────────────────────────────────────────────
print("\nSTEP 6: Verification...")
import subprocess
result = subprocess.run(['ls', '-la', f'{figdir}'], capture_output=True, text=True)
print(result.stdout)

# Check all expected files
expected = [
    'fig_mechanism.pdf', 'fig_sal_mod.pdf', 'fig_fjord_depth.pdf',
    'fig_regional_profiles.pdf', 'fig_asymmetry.pdf', 'fig_sensitivity.pdf',
    'fig_suppression_bars.pdf', 'fig_monthly_panels.pdf', 'fig_dfo_validation.pdf',
    'fig_depression_heatmap.pdf', 'fig_depth_exp.pdf',
]
missing = []
for fname in expected:
    fpath = os.path.join(figdir, fname)
    if not os.path.exists(fpath):
        missing.append(fname)
    else:
        size = os.path.getsize(fpath)
        if size == 0:
            missing.append(f"{fname} (EMPTY)")

if missing:
    print(f"  ❌ Missing/empty: {missing}")
else:
    print(f"  ✅ All {len(expected)} figures generated successfully")

print("\n✅ DONE - All figures and data_summary.json generated.")
