#!/usr/bin/env python3
"""Generate all figures for W115-W124 T_vbnc_min calibration report."""

import json
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as mpatches

BASE = '/home/starbot/.openclaw/workspace/sswd-evoepi/results/calibration'
FIGDIR = '/home/starbot/.openclaw/workspace/sswd-evoepi/reports/tvbnc_min_w115_w124/figures'
os.makedirs(FIGDIR, exist_ok=True)

RUNS = ['W115','W116','W117','W118','W119','W120','W121','W122','W123','W124']
PREV_RUNS = ['W105','W106','W107','W108','W109','W110']

TARGETS = {
    'AK-PWS': 0.50, 'AK-FN': 0.50, 'AK-FS': 0.20, 'BC-N': 0.20,
    'SS-S': 0.05, 'JDF': 0.02, 'OR': 0.0025, 'CA-N': 0.001
}
CALIB_REGIONS = list(TARGETS.keys())
KEY_REGIONS = ['AK-PWS', 'BC-N', 'SS-S', 'CA-N']
ALL_REGIONS_ORDERED = [
    'AK-WG','AK-AL','AK-PWS','AK-OC','AK-EG','AK-FN','AK-FS',
    'BC-N','BC-C','SS-N','SS-S','JDF','WA-O','OR','CA-N','CA-C','CA-S','BJ'
]

# Colors
COLORS = {
    'AK-PWS': '#1f77b4', 'AK-FN': '#17becf', 'AK-FS': '#7fcdbb',
    'BC-N': '#2ca02c', 'SS-S': '#ff7f0e', 'JDF': '#d62728',
    'OR': '#9467bd', 'CA-N': '#e377c2', 'CA-S': '#8c564b', 'CA-C': '#bcbd22'
}

def load_run(run_name):
    path = f'{BASE}/{run_name}/combined_results.json'
    with open(path) as f:
        return json.load(f)

def load_all():
    data = {}
    for w in RUNS:
        data[w] = load_run(w)
    for w in PREV_RUNS:
        p = f'{BASE}/{w}/combined_results.json'
        if os.path.exists(p):
            data[w] = load_run(w)
    return data

data = load_all()

# ============================================================
# Figure 1: RMSE ranking bar chart
# ============================================================
def fig1_rmse_ranking():
    fig, ax = plt.subplots(figsize=(10, 5))
    
    # Current round
    rmses = [(w, data[w]['mean_rmse_log']) for w in RUNS]
    rmses.sort(key=lambda x: x[1])
    
    names = [r[0] for r in rmses]
    vals = [r[1] for r in rmses]
    
    # Color by T_vbnc_min
    colors = []
    for w, _ in rmses:
        p = data[w]['param_overrides']
        t_min = p.get('disease.T_vbnc_min', 9)
        v_evo = p.get('disease.virulence_evolution', True)
        if not v_evo:
            colors.append('#e74c3c')  # red = control
        elif t_min == 10:
            colors.append('#3498db')  # blue = T_min=10
        else:
            colors.append('#2ecc71')  # green = T_min=9
    
    bars = ax.bar(names, vals, color=colors, edgecolor='black', linewidth=0.5)
    
    # Add previous best
    prev_best = min((data[w]['mean_rmse_log'], w) for w in PREV_RUNS if w in data)
    ax.axhline(prev_best[0], color='gray', linestyle='--', linewidth=1.5, label=f'Best previous ({prev_best[1]}, RMSE={prev_best[0]:.3f})')
    
    # Value labels
    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                f'{val:.3f}', ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    # Legend
    patches = [
        mpatches.Patch(color='#2ecc71', label='T_vbnc_min=9°C, v_evo ON'),
        mpatches.Patch(color='#3498db', label='T_vbnc_min=10°C, v_evo ON'),
        mpatches.Patch(color='#e74c3c', label='T_vbnc_min=9°C, v_evo OFF (control)'),
    ]
    ax.legend(handles=patches + [plt.Line2D([0],[0], color='gray', linestyle='--', label=f'Best previous ({prev_best[1]})')],
              loc='upper right', fontsize=9)
    
    ax.set_ylabel('RMSE (log-scale)', fontsize=12)
    ax.set_title('W115–W124: RMSE Ranking (Lower = Better)', fontsize=14, fontweight='bold')
    ax.set_ylim(0, max(vals) * 1.2)
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    fig.savefig(f'{FIGDIR}/fig1_rmse_ranking.pdf', dpi=150)
    plt.close()
    print("Fig1 done")

# ============================================================
# Figure 2: Per-region recovery bar chart (best runs vs targets)
# ============================================================
def fig2_recovery_bars():
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Pick best 3 + controls
    best_runs = ['W117', 'W121', 'W118', 'W116']  # sorted by RMSE
    x = np.arange(len(CALIB_REGIONS))
    width = 0.15
    
    # Targets
    target_vals = [TARGETS[r]*100 for r in CALIB_REGIONS]
    ax.bar(x - 2*width, target_vals, width, label='Target', color='#27ae60', alpha=0.8, edgecolor='black', linewidth=0.5)
    
    run_colors = ['#3498db', '#e74c3c', '#9b59b6', '#f39c12']
    for i, w in enumerate(best_runs):
        rec = data[w]['results'][0]['region_recovery']
        vals = [rec.get(r, 0)*100 for r in CALIB_REGIONS]
        p = data[w]['param_overrides']
        t_min = p.get('disease.T_vbnc_min', 9)
        v_evo = p.get('disease.virulence_evolution', True)
        label = f'{w} (T_min={t_min}°C, v_evo={"ON" if v_evo else "OFF"})'
        ax.bar(x + (i-1)*width, vals, width, label=label, color=run_colors[i], alpha=0.8, edgecolor='black', linewidth=0.5)
    
    ax.set_xticks(x)
    ax.set_xticklabels(CALIB_REGIONS, fontsize=10, rotation=45, ha='right')
    ax.set_ylabel('Recovery Fraction (%)', fontsize=12)
    ax.set_title('Recovery by Region: Best Runs vs Targets', fontsize=14, fontweight='bold')
    ax.legend(fontsize=9, loc='upper right')
    ax.set_yscale('log')
    ax.set_ylim(0.05, 100)
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    fig.savefig(f'{FIGDIR}/fig2_recovery_bars.pdf', dpi=150)
    plt.close()
    print("Fig2 done")

# ============================================================
# Figure 3: Recovery time series for 4 key regions
# ============================================================
def fig3_timeseries():
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()
    
    best_runs = ['W117', 'W121', 'W123', 'W115']
    run_styles = {
        'W117': ('#3498db', '-', 'W117 (T_min=10, v_evo ON)'),
        'W121': ('#e74c3c', '--', 'W121 (T_min=9, v_evo OFF)'),
        'W123': ('#9b59b6', '-.', 'W123 (T_min=10, v_evo ON, seed 123)'),
        'W115': ('#2ecc71', ':', 'W115 (T_min=9, v_evo ON)'),
    }
    
    years = list(range(1, 14))
    
    for i, reg in enumerate(KEY_REGIONS):
        ax = axes[i]
        target = TARGETS.get(reg, None)
        
        for w in best_runs:
            rd = data[w]['results'][0]['region_details']
            if reg in rd:
                yt = rd[reg].get('yearly_totals', [])
                peak = rd[reg].get('peak_pop', 1)
                pcts = [100 * v / peak if peak > 0 else 0 for v in yt]
                color, ls, label = run_styles[w]
                ax.plot(years[:len(pcts)], pcts, color=color, linestyle=ls, linewidth=2, label=label, marker='o', markersize=3)
        
        if target:
            ax.axhline(target * 100, color='#27ae60', linestyle='--', linewidth=1.5, alpha=0.7, label=f'Target ({target*100:.1f}%)')
        
        # La Niña shading (2020-2021 = years 8-9)
        ax.axvspan(10, 11, alpha=0.1, color='blue', label='La Niña cooling' if i == 0 else None)
        
        ax.set_title(reg, fontsize=14, fontweight='bold')
        ax.set_xlabel('Year', fontsize=10)
        ax.set_ylabel('% of Peak Population', fontsize=10)
        ax.set_xlim(1, 13)
        ax.set_ylim(0, 105)
        ax.grid(alpha=0.3)
        if i == 0:
            ax.legend(fontsize=7, loc='lower left')
    
    fig.suptitle('Recovery Trajectories: Key Regions (Year 1–13)', fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()
    fig.savefig(f'{FIGDIR}/fig3_timeseries.pdf', dpi=150, bbox_inches='tight')
    plt.close()
    print("Fig3 done")

# ============================================================
# Figure 4: T_vbnc_min effect comparison (9 vs 10 vs 4-6)
# ============================================================
def fig4_tvbnc_effect():
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    
    # Group: T_min=9 (W115,W116), T_min=10 (W117,W118), T_min=4-6 previous (W105,W106)
    groups = {
        'T_vbnc_min = 4°C\n(W105–W106, previous)': ['W105', 'W106'],
        'T_vbnc_min = 9°C\n(W115–W116)': ['W115', 'W116'],
        'T_vbnc_min = 10°C\n(W117–W118)': ['W117', 'W118'],
    }
    group_colors = ['#95a5a6', '#2ecc71', '#3498db']
    
    # Panel A: Mean recovery per region
    ax = axes[0]
    x = np.arange(len(CALIB_REGIONS))
    width = 0.25
    for gi, (gname, gruns) in enumerate(groups.items()):
        mean_rec = []
        for reg in CALIB_REGIONS:
            vals = [data[w]['results'][0]['region_recovery'].get(reg, 0) * 100 for w in gruns if w in data]
            mean_rec.append(np.mean(vals) if vals else 0)
        ax.bar(x + (gi-1)*width, mean_rec, width, label=gname, color=group_colors[gi], edgecolor='black', linewidth=0.5)
    
    # Targets
    for xi, reg in enumerate(CALIB_REGIONS):
        ax.plot(xi, TARGETS[reg]*100, '*', color='gold', markersize=12, zorder=5)
    
    ax.set_xticks(x)
    ax.set_xticklabels(CALIB_REGIONS, rotation=45, ha='right', fontsize=9)
    ax.set_ylabel('Recovery (%)', fontsize=11)
    ax.set_title('(a) Recovery by T_vbnc_min', fontsize=12, fontweight='bold')
    ax.legend(fontsize=7, loc='upper right')
    ax.set_yscale('log')
    ax.set_ylim(0.05, 100)
    ax.grid(axis='y', alpha=0.3)
    
    # Panel B: AK-PWS time series comparison
    ax = axes[1]
    years = list(range(1, 14))
    for gi, (gname, gruns) in enumerate(groups.items()):
        for w in gruns:
            if w not in data:
                continue
            rd = data[w]['results'][0]['region_details']
            if 'AK-PWS' in rd:
                yt = rd['AK-PWS']['yearly_totals']
                peak = rd['AK-PWS']['peak_pop']
                pcts = [100*v/peak for v in yt]
                ax.plot(years[:len(pcts)], pcts, color=group_colors[gi], linewidth=1.5, alpha=0.7)
    
    ax.axhline(50, color='#27ae60', linestyle='--', alpha=0.5, label='Target 50%')
    ax.set_title('(b) AK-PWS Trajectory', fontsize=12, fontweight='bold')
    ax.set_xlabel('Year')
    ax.set_ylabel('% of Peak')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    
    # Panel C: RMSE comparison
    ax = axes[2]
    for gi, (gname, gruns) in enumerate(groups.items()):
        rmses = [data[w]['mean_rmse_log'] for w in gruns if w in data]
        ax.bar(gi, np.mean(rmses), color=group_colors[gi], edgecolor='black', linewidth=0.5, yerr=np.std(rmses) if len(rmses)>1 else 0, capsize=5)
        ax.text(gi, np.mean(rmses) + 0.02, f'{np.mean(rmses):.3f}', ha='center', fontsize=11, fontweight='bold')
    
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(['4°C\n(prev)', '9°C', '10°C'], fontsize=10)
    ax.set_ylabel('RMSE', fontsize=11)
    ax.set_title('(c) RMSE by T_vbnc_min', fontsize=12, fontweight='bold')
    ax.set_ylim(0, 1.0)
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    fig.savefig(f'{FIGDIR}/fig4_tvbnc_effect.pdf', dpi=150)
    plt.close()
    print("Fig4 done")

# ============================================================
# Figure 5: Spatial pathogen evolution (T_vbnc and v_local by region)
# ============================================================
def fig5_pathogen_evolution():
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Top row: T_vbnc by region for W115 (T_min=9) and W117 (T_min=10)
    for col, (w, tmin) in enumerate([('W115', 9), ('W117', 10)]):
        ax = axes[0][col]
        rd = data[w]['results'][0]['region_details']
        regions = [r for r in ALL_REGIONS_ORDERED if r in rd]
        tvbnc_vals = [rd[r].get('final_mean_T_vbnc', 0) for r in regions]
        
        colors = ['#3498db' if r.startswith('AK') else '#2ecc71' if r.startswith('BC') else '#f39c12' if r.startswith('SS') or r == 'JDF' else '#e74c3c' for r in regions]
        
        bars = ax.bar(range(len(regions)), tvbnc_vals, color=colors, edgecolor='black', linewidth=0.5)
        ax.axhline(tmin, color='red', linestyle='--', linewidth=2, label=f'T_vbnc_min = {tmin}°C')
        ax.set_xticks(range(len(regions)))
        ax.set_xticklabels(regions, rotation=60, ha='right', fontsize=7)
        ax.set_ylabel('Final T_vbnc (°C)', fontsize=10)
        ax.set_title(f'({"a" if col==0 else "b"}) T_vbnc Adaptation — {w} (T_min={tmin}°C)', fontsize=11, fontweight='bold')
        ax.legend(fontsize=9)
        ax.set_ylim(8, 13)
        ax.grid(axis='y', alpha=0.3)
    
    # Bottom row: v_local by region for W115 and W117
    for col, (w, tmin) in enumerate([('W115', 9), ('W117', 10)]):
        ax = axes[1][col]
        rd = data[w]['results'][0]['region_details']
        regions = [r for r in ALL_REGIONS_ORDERED if r in rd]
        vlocal_vals = [rd[r].get('final_mean_v_local', 0) for r in regions]
        
        colors = ['#3498db' if r.startswith('AK') else '#2ecc71' if r.startswith('BC') else '#f39c12' if r.startswith('SS') or r == 'JDF' else '#e74c3c' for r in regions]
        
        bars = ax.bar(range(len(regions)), vlocal_vals, color=colors, edgecolor='black', linewidth=0.5)
        ax.set_xticks(range(len(regions)))
        ax.set_xticklabels(regions, rotation=60, ha='right', fontsize=7)
        ax.set_ylabel('Final v_local', fontsize=10)
        ax.set_title(f'({"c" if col==0 else "d"}) Virulence Adaptation — {w} (T_min={tmin}°C)', fontsize=11, fontweight='bold')
        ax.set_ylim(0, 0.35)
        ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    fig.savefig(f'{FIGDIR}/fig5_pathogen_evolution.pdf', dpi=150)
    plt.close()
    print("Fig5 done")

# ============================================================
# Figure 6: Control comparison (v_evo OFF vs ON)
# ============================================================
def fig6_control_comparison():
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    
    # W115/W116 (v_evo ON, T_min=9) vs W121/W122 (v_evo OFF, T_min=9)
    on_runs = ['W115', 'W116']
    off_runs = ['W121', 'W122']
    
    # Panel A: Recovery comparison
    ax = axes[0]
    x = np.arange(len(CALIB_REGIONS))
    width = 0.3
    
    on_rec = []
    off_rec = []
    for reg in CALIB_REGIONS:
        on_vals = [data[w]['results'][0]['region_recovery'].get(reg, 0)*100 for w in on_runs]
        off_vals = [data[w]['results'][0]['region_recovery'].get(reg, 0)*100 for w in off_runs]
        on_rec.append(np.mean(on_vals))
        off_rec.append(np.mean(off_vals))
    
    ax.bar(x - width/2, on_rec, width, label='v_evo ON (W115–116)', color='#3498db', edgecolor='black', linewidth=0.5)
    ax.bar(x + width/2, off_rec, width, label='v_evo OFF (W121–122)', color='#e74c3c', edgecolor='black', linewidth=0.5)
    
    for xi, reg in enumerate(CALIB_REGIONS):
        ax.plot(xi, TARGETS[reg]*100, '*', color='gold', markersize=12, zorder=5)
    
    ax.set_xticks(x)
    ax.set_xticklabels(CALIB_REGIONS, rotation=45, ha='right', fontsize=9)
    ax.set_ylabel('Recovery (%)', fontsize=11)
    ax.set_title('(a) Recovery: v_evo ON vs OFF', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9)
    ax.set_yscale('log')
    ax.set_ylim(0.05, 100)
    ax.grid(axis='y', alpha=0.3)
    
    # Panel B: v_local comparison
    ax = axes[1]
    regions = [r for r in ALL_REGIONS_ORDERED if r in data['W115']['results'][0]['region_details']]
    
    on_vlocal = []
    off_vlocal = []
    for reg in regions:
        on_vals = [data[w]['results'][0]['region_details'][reg].get('final_mean_v_local', 0) for w in on_runs]
        off_vals = [data[w]['results'][0]['region_details'][reg].get('final_mean_v_local', 0) for w in off_runs]
        on_vlocal.append(np.mean(on_vals))
        off_vlocal.append(np.mean(off_vals))
    
    x2 = np.arange(len(regions))
    ax.bar(x2 - width/2, on_vlocal, width, label='v_evo ON', color='#3498db', edgecolor='black', linewidth=0.5)
    ax.bar(x2 + width/2, off_vlocal, width, label='v_evo OFF (fixed 0.5)', color='#e74c3c', edgecolor='black', linewidth=0.5)
    ax.set_xticks(x2)
    ax.set_xticklabels(regions, rotation=60, ha='right', fontsize=7)
    ax.set_ylabel('v_local', fontsize=11)
    ax.set_title('(b) Final Virulence', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(axis='y', alpha=0.3)
    
    # Panel C: RMSE
    ax = axes[2]
    on_rmses = [data[w]['mean_rmse_log'] for w in on_runs]
    off_rmses = [data[w]['mean_rmse_log'] for w in off_runs]
    
    ax.bar(0, np.mean(on_rmses), color='#3498db', edgecolor='black', linewidth=0.5, yerr=np.std(on_rmses), capsize=5)
    ax.bar(1, np.mean(off_rmses), color='#e74c3c', edgecolor='black', linewidth=0.5, yerr=np.std(off_rmses), capsize=5)
    ax.text(0, np.mean(on_rmses)+0.02, f'{np.mean(on_rmses):.3f}', ha='center', fontsize=12, fontweight='bold')
    ax.text(1, np.mean(off_rmses)+0.02, f'{np.mean(off_rmses):.3f}', ha='center', fontsize=12, fontweight='bold')
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['v_evo ON\n(W115–116)', 'v_evo OFF\n(W121–122)'], fontsize=11)
    ax.set_ylabel('RMSE', fontsize=11)
    ax.set_title('(c) RMSE Comparison', fontsize=12, fontweight='bold')
    ax.set_ylim(0, 1.0)
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    fig.savefig(f'{FIGDIR}/fig6_control_comparison.pdf', dpi=150)
    plt.close()
    print("Fig6 done")

# ============================================================
# Figure 7: All 18 regions, best run vs targets
# ============================================================
def fig7_all_regions():
    fig, ax = plt.subplots(figsize=(14, 6))
    
    w = 'W117'  # best RMSE
    rd = data[w]['results'][0]['region_details']
    regions = [r for r in ALL_REGIONS_ORDERED if r in rd]
    
    rec_vals = [rd[r].get('recovery_frac', 0)*100 for r in regions]
    
    colors = ['#3498db' if r.startswith('AK') else '#2ecc71' if r.startswith('BC') else '#f39c12' if 'SS' in r or r == 'JDF' or r == 'WA-O' else '#e74c3c' for r in regions]
    
    bars = ax.bar(range(len(regions)), rec_vals, color=colors, edgecolor='black', linewidth=0.5, label=f'{w} (RMSE={data[w]["mean_rmse_log"]:.3f})')
    
    # Targets as stars
    for xi, reg in enumerate(regions):
        if reg in TARGETS:
            ax.plot(xi, TARGETS[reg]*100, '*', color='gold', markersize=15, zorder=5)
    
    ax.set_xticks(range(len(regions)))
    ax.set_xticklabels(regions, rotation=45, ha='right', fontsize=9)
    ax.set_ylabel('Recovery Fraction (%)', fontsize=12)
    ax.set_title(f'All 18 Regions — {w} (Best RMSE={data[w]["mean_rmse_log"]:.3f})', fontsize=14, fontweight='bold')
    ax.set_yscale('log')
    ax.set_ylim(0.001, 100)
    ax.grid(axis='y', alpha=0.3)
    
    # Region group labels
    ax.text(3, 80, 'Alaska', fontsize=11, ha='center', color='#3498db', fontweight='bold')
    ax.text(8, 80, 'BC/Salish', fontsize=11, ha='center', color='#2ecc71', fontweight='bold')
    ax.text(13, 80, 'South', fontsize=11, ha='center', color='#e74c3c', fontweight='bold')
    
    patches = [mpatches.Patch(color='gold', label='Calibration target')]
    ax.legend(handles=patches, fontsize=10, loc='upper right')
    
    plt.tight_layout()
    fig.savefig(f'{FIGDIR}/fig7_all_regions.pdf', dpi=150)
    plt.close()
    print("Fig7 done")

# ============================================================
# Figure 8: T_vbnc adaptation distance from floor
# ============================================================
def fig8_tvbnc_adaptation_distance():
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Compare T_vbnc adaptation: T_min=4 (prev) vs T_min=9 vs T_min=10
    comparisons = [
        ('W106', 4, 'T_min=4°C (W106)', '#95a5a6'),
        ('W115', 9, 'T_min=9°C (W115)', '#2ecc71'),
        ('W117', 10, 'T_min=10°C (W117)', '#3498db'),
    ]
    
    ax = axes[0]
    x = np.arange(len(ALL_REGIONS_ORDERED))
    width = 0.25
    
    for gi, (w, tmin, label, color) in enumerate(comparisons):
        rd = data[w]['results'][0]['region_details']
        tvbnc_vals = [rd[r].get('final_mean_T_vbnc', tmin) - tmin for r in ALL_REGIONS_ORDERED if r in rd]
        regions = [r for r in ALL_REGIONS_ORDERED if r in rd]
        ax.bar(np.arange(len(regions)) + (gi-1)*width, tvbnc_vals, width, label=label, color=color, edgecolor='black', linewidth=0.5)
    
    ax.set_xticks(np.arange(len(regions)))
    ax.set_xticklabels(regions, rotation=60, ha='right', fontsize=7)
    ax.set_ylabel('ΔT_vbnc above floor (°C)', fontsize=11)
    ax.set_title('(a) T_vbnc Adaptation Above Floor', fontsize=12, fontweight='bold')
    ax.legend(fontsize=8)
    ax.grid(axis='y', alpha=0.3)
    
    # Panel B: Did AK hit the floor?
    ax = axes[1]
    ak_regions = [r for r in ALL_REGIONS_ORDERED if r.startswith('AK')]
    for gi, (w, tmin, label, color) in enumerate(comparisons):
        rd = data[w]['results'][0]['region_details']
        tvbnc_vals = [rd[r].get('final_mean_T_vbnc', tmin) for r in ak_regions if r in rd]
        ax.bar(np.arange(len(ak_regions)) + (gi-1)*width, tvbnc_vals, width, label=label, color=color, edgecolor='black', linewidth=0.5)
    
    ax.axhline(9, color='#2ecc71', linestyle='--', linewidth=2, alpha=0.7, label='Floor = 9°C')
    ax.axhline(10, color='#3498db', linestyle='--', linewidth=2, alpha=0.7, label='Floor = 10°C')
    
    ax.set_xticks(np.arange(len(ak_regions)))
    ax.set_xticklabels(ak_regions, fontsize=9)
    ax.set_ylabel('Final T_vbnc (°C)', fontsize=11)
    ax.set_title('(b) Alaska T_vbnc: Did It Hit the Floor?', fontsize=12, fontweight='bold')
    ax.legend(fontsize=8)
    ax.set_ylim(3, 12)
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    fig.savefig(f'{FIGDIR}/fig8_tvbnc_floor.pdf', dpi=150)
    plt.close()
    print("Fig8 done")

# ============================================================
# Run all
# ============================================================
fig1_rmse_ranking()
fig2_recovery_bars()
fig3_timeseries()
fig4_tvbnc_effect()
fig5_pathogen_evolution()
fig6_control_comparison()
fig7_all_regions()
fig8_tvbnc_adaptation_distance()
print("\nAll figures generated!")
