#!/usr/bin/env python3
"""Generate all 6 figures for the W05-W16 wavefront calibration report."""

import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import cm
import os

RESULTS_DIR = '/home/starbot/.openclaw/workspace/sswd-evoepi/results/calibration'
FIG_DIR = '/home/starbot/.openclaw/workspace/sswd-evoepi/reports/wavefront_calibration_r2/figures'

# Parameter configs
CONFIGS = {
    'W05': {'k_vbnc': 1.5, 'threshold': 100, 'P_env_max': 2000},
    'W06': {'k_vbnc': 1.5, 'threshold': 100, 'P_env_max': 5000},
    'W07': {'k_vbnc': 1.5, 'threshold': 200, 'P_env_max': 2000},
    'W08': {'k_vbnc': 1.5, 'threshold': 200, 'P_env_max': 5000},
    'W09': {'k_vbnc': 1.5, 'threshold': 500, 'P_env_max': 2000},
    'W10': {'k_vbnc': 1.5, 'threshold': 500, 'P_env_max': 5000},
    'W11': {'k_vbnc': 2.0, 'threshold': 100, 'P_env_max': 2000},
    'W12': {'k_vbnc': 2.0, 'threshold': 100, 'P_env_max': 5000},
    'W13': {'k_vbnc': 2.0, 'threshold': 200, 'P_env_max': 2000},
    'W14': {'k_vbnc': 2.0, 'threshold': 200, 'P_env_max': 5000},
    'W15': {'k_vbnc': 2.0, 'threshold': 500, 'P_env_max': 2000},
    'W16': {'k_vbnc': 2.0, 'threshold': 500, 'P_env_max': 5000},
}

RECOVERY_REGIONS = ['AK-PWS', 'AK-FN', 'AK-FS', 'BC-N', 'SS-S', 'JDF', 'OR', 'CA-N']
SST = {'AK-PWS': 8.4, 'AK-FN': 8.7, 'AK-FS': 9.0, 'BC-N': 10.1, 'SS-S': 10.1, 'JDF': 10.0, 'OR': 11.4, 'CA-N': 11.6}

TIMING_REGIONS = ['CA-S', 'CA-C', 'CA-N', 'OR', 'WA-O', 'JDF', 'SS-S', 'SS-N', 'BC-C', 'BC-N',
                  'AK-FS', 'AK-FN', 'AK-PWS', 'AK-EG', 'AK-OC', 'AK-WG', 'AK-AL']

def load_data():
    """Load all W05-W16 and W01-W04 results."""
    data = {}
    # W05-W16
    for w in range(5, 17):
        name = f'W{w:02d}'
        with open(f'{RESULTS_DIR}/W05-W16/{name}.json') as f:
            data[name] = json.load(f)
    # W01-W04
    for w in range(1, 5):
        name = f'W{w:02d}'
        with open(f'{RESULTS_DIR}/{name}/result_seed42.json') as f:
            data[name] = json.load(f)
    return data

def fig1_recovery_bar(data):
    """Recovery fractions: model vs targets, best round (W15) highlighted."""
    fig, ax = plt.subplots(figsize=(10, 5.5))
    
    targets = [data['W15']['scoring']['per_region'][r]['target_pct'] for r in RECOVERY_REGIONS]
    
    # Show W15 (best) plus a few others for comparison
    rounds_to_show = ['W05', 'W09', 'W11', 'W13', 'W15']
    
    x = np.arange(len(RECOVERY_REGIONS))
    width = 0.12
    n = len(rounds_to_show) + 1  # +1 for targets
    
    # Plot targets
    bars_target = ax.bar(x - width * n/2, targets, width, label='Target', color='black', alpha=0.3, edgecolor='black')
    
    colors = ['#7fcdbb', '#41b6c4', '#1d91c0', '#225ea8', '#c51b7d']
    for i, rnd in enumerate(rounds_to_show):
        vals = [data[rnd]['scoring']['per_region'][r]['actual_pct'] for r in RECOVERY_REGIONS]
        edgecolor = 'red' if rnd == 'W15' else 'none'
        linewidth = 2.0 if rnd == 'W15' else 0.5
        bars = ax.bar(x - width * n/2 + width * (i + 1), vals, width, 
                      label=rnd + (' (best)' if rnd == 'W15' else ''),
                      color=colors[i], edgecolor=edgecolor, linewidth=linewidth)
    
    ax.set_xticks(x)
    ax.set_xticklabels(RECOVERY_REGIONS, fontsize=9)
    ax.set_ylabel('Recovery (%)', fontsize=11)
    ax.set_title('Recovery Fractions: Model vs Targets (P_env_max=2000 rounds)', fontsize=12)
    ax.legend(fontsize=8, loc='upper right')
    ax.set_ylim(0, 100)
    
    # Add a horizontal note
    ax.axhline(y=50, color='gray', linestyle=':', alpha=0.3)
    
    plt.tight_layout()
    fig.savefig(f'{FIG_DIR}/fig1_recovery_bars.pdf', bbox_inches='tight')
    plt.close()
    print("Fig 1 saved")

def fig2_rmse_bars(data):
    """RMSE bar chart for all 12 rounds, colored by P_env_max."""
    fig, ax = plt.subplots(figsize=(10, 5))
    
    rounds = [f'W{w:02d}' for w in range(5, 17)]
    rmses = [data[r]['scoring']['rmse_log'] for r in rounds]
    colors = ['#2166ac' if CONFIGS[r]['P_env_max'] == 2000 else '#b2182b' for r in rounds]
    
    bars = ax.bar(range(len(rounds)), rmses, color=colors, edgecolor='black', linewidth=0.5)
    
    # Highlight best
    best_idx = np.argmin(rmses)
    bars[best_idx].set_edgecolor('gold')
    bars[best_idx].set_linewidth(3)
    
    ax.set_xticks(range(len(rounds)))
    ax.set_xticklabels(rounds, rotation=45, fontsize=9)
    ax.set_ylabel('RMSE (log scale)', fontsize=11)
    ax.set_title('Recovery RMSE by Calibration Round', fontsize=12)
    
    # Custom legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='#2166ac', label='P_env_max=2000'),
                       Patch(facecolor='#b2182b', label='P_env_max=5000'),
                       Patch(facecolor='white', edgecolor='gold', linewidth=2, label='Best (W15)')]
    ax.legend(handles=legend_elements, fontsize=9)
    
    # Add value labels
    for i, (bar, val) in enumerate(zip(bars, rmses)):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                f'{val:.2f}', ha='center', va='bottom', fontsize=7)
    
    ax.set_ylim(0, max(rmses) * 1.15)
    plt.tight_layout()
    fig.savefig(f'{FIG_DIR}/fig2_rmse_bars.pdf', bbox_inches='tight')
    plt.close()
    print("Fig 2 saved")

def fig3_heatmap(data):
    """Parameter sensitivity heatmap: k_vbnc × threshold × P_env_max showing RMSE."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), sharey=True)
    
    thresholds = [100, 200, 500]
    k_vals = [1.5, 2.0]
    p_vals = [2000, 5000]
    
    for pidx, pval in enumerate(p_vals):
        ax = axes[pidx]
        grid = np.zeros((len(k_vals), len(thresholds)))
        
        for ki, kv in enumerate(k_vals):
            for ti, tv in enumerate(thresholds):
                # Find matching round
                for rname, cfg in CONFIGS.items():
                    if cfg['k_vbnc'] == kv and cfg['threshold'] == tv and cfg['P_env_max'] == pval:
                        grid[ki, ti] = data[rname]['scoring']['rmse_log']
                        break
        
        im = ax.imshow(grid, cmap='RdYlGn_r', aspect='auto', vmin=0.35, vmax=1.7)
        
        # Labels
        ax.set_xticks(range(len(thresholds)))
        ax.set_xticklabels(thresholds)
        ax.set_yticks(range(len(k_vals)))
        ax.set_yticklabels(k_vals)
        ax.set_xlabel('Activation Threshold', fontsize=10)
        if pidx == 0:
            ax.set_ylabel('k_vbnc', fontsize=10)
        ax.set_title(f'P_env_max = {pval}', fontsize=11)
        
        # Annotate
        for ki in range(len(k_vals)):
            for ti in range(len(thresholds)):
                val = grid[ki, ti]
                color = 'white' if val > 1.0 else 'black'
                ax.text(ti, ki, f'{val:.3f}', ha='center', va='center', fontsize=11, fontweight='bold', color=color)
    
    fig.colorbar(im, ax=axes, label='RMSE (log)', shrink=0.8)
    fig.suptitle('Parameter Sensitivity: RMSE across k_vbnc × Threshold × P_env_max', fontsize=12, y=1.02)
    plt.tight_layout()
    fig.savefig(f'{FIG_DIR}/fig3_heatmap.pdf', bbox_inches='tight')
    plt.close()
    print("Fig 3 saved")

def fig4_recovery_gradient(data):
    """Recovery fraction vs mean SST for each round, with target points."""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Targets
    target_sst = [SST[r] for r in RECOVERY_REGIONS]
    target_rec = [data['W15']['scoring']['per_region'][r]['target_pct'] for r in RECOVERY_REGIONS]
    
    # Plot each P_env_max=2000 round as a line
    rounds_2k = ['W05', 'W07', 'W09', 'W11', 'W13', 'W15']
    rounds_5k = ['W06', 'W08', 'W10', 'W12', 'W14', 'W16']
    
    cmap_2k = plt.cm.Blues(np.linspace(0.3, 0.9, len(rounds_2k)))
    cmap_5k = plt.cm.Reds(np.linspace(0.3, 0.9, len(rounds_5k)))
    
    for i, rnd in enumerate(rounds_2k):
        ssts = [SST[r] for r in RECOVERY_REGIONS]
        recs = [data[rnd]['scoring']['per_region'][r]['actual_pct'] for r in RECOVERY_REGIONS]
        lw = 2.5 if rnd == 'W15' else 1.0
        alpha = 1.0 if rnd == 'W15' else 0.5
        ax.plot(ssts, recs, 'o-', color=cmap_2k[i], label=f'{rnd} (2K)', linewidth=lw, alpha=alpha, markersize=4)
    
    for i, rnd in enumerate(rounds_5k):
        ssts = [SST[r] for r in RECOVERY_REGIONS]
        recs = [data[rnd]['scoring']['per_region'][r]['actual_pct'] for r in RECOVERY_REGIONS]
        ax.plot(ssts, recs, 's--', color=cmap_5k[i], label=f'{rnd} (5K)', linewidth=0.8, alpha=0.5, markersize=3)
    
    # Targets
    ax.scatter(target_sst, target_rec, marker='*', s=200, color='black', zorder=10, label='Targets')
    
    # Region labels for targets
    for r in RECOVERY_REGIONS:
        tgt = data['W15']['scoring']['per_region'][r]['target_pct']
        ax.annotate(r, (SST[r], tgt), textcoords="offset points", xytext=(5, 5), fontsize=7, color='black')
    
    ax.set_xlabel('Mean SST (°C)', fontsize=11)
    ax.set_ylabel('Recovery (%)', fontsize=11)
    ax.set_title('Recovery Gradient: Recovery vs Sea Surface Temperature', fontsize=12)
    ax.legend(fontsize=7, loc='upper right', ncol=2)
    ax.set_ylim(-2, 100)
    ax.invert_xaxis()  # Cold (high recovery) on left
    
    plt.tight_layout()
    fig.savefig(f'{FIG_DIR}/fig4_recovery_gradient.pdf', bbox_inches='tight')
    plt.close()
    print("Fig 4 saved")

def fig5_timing(data):
    """Arrival timing comparison: W15 vs targets."""
    fig, ax = plt.subplots(figsize=(12, 5.5))
    
    regions_show = TIMING_REGIONS
    
    targets = []
    w15_actual = []
    for r in regions_show:
        tr = data['W15']['arrival_timing']['per_region'][r]
        targets.append(tr['target_months'])
        act = tr['actual_months']
        if act is None:
            w15_actual.append(np.nan)
        else:
            w15_actual.append(act)
    
    x = np.arange(len(regions_show))
    width = 0.35
    
    bars1 = ax.bar(x - width/2, targets, width, label='Target', color='#2166ac', alpha=0.7, edgecolor='black', linewidth=0.5)
    
    # For W15, mark unreached regions differently
    colors_w15 = []
    for i, r in enumerate(regions_show):
        if np.isnan(w15_actual[i]):
            colors_w15.append('#d9d9d9')
        else:
            colors_w15.append('#b2182b')
    
    # Replace NaN with penalty value for display
    w15_display = [v if not np.isnan(v) else 0 for v in w15_actual]
    bars2 = ax.bar(x + width/2, w15_display, width, label='W15 (best)', color=colors_w15, alpha=0.7, edgecolor='black', linewidth=0.5)
    
    # Mark unreached with X
    for i, r in enumerate(regions_show):
        if np.isnan(w15_actual[i]):
            ax.text(x[i] + width/2, 2, '✗', ha='center', va='bottom', fontsize=14, color='red', fontweight='bold')
            ax.text(x[i] + width/2, 6, 'N/R', ha='center', va='bottom', fontsize=7, color='red')
    
    ax.set_xticks(x)
    ax.set_xticklabels(regions_show, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('Arrival Time (months)', fontsize=11)
    ax.set_title('Wavefront Arrival Timing: W15 vs Targets', fontsize=12)
    ax.legend(fontsize=9)
    ax.set_ylim(0, 50)
    
    # Add vertical lines to separate groups
    for pos in [0.5, 2.5, 4.5, 9.5, 11.5, 14.5]:
        ax.axvline(x=pos, color='gray', linestyle=':', alpha=0.3)
    
    plt.tight_layout()
    fig.savefig(f'{FIG_DIR}/fig5_timing.pdf', bbox_inches='tight')
    plt.close()
    print("Fig 5 saved")

def fig6_before_after(data):
    """Before/after comparison: W01-W04 average vs W15."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # --- Panel A: Recovery ---
    ax = axes[0]
    
    # W01-W04 average
    w_old_avg = {}
    for r in RECOVERY_REGIONS:
        vals = [data[f'W{w:02d}']['scoring']['per_region'][r]['actual_pct'] for w in range(1, 5)]
        w_old_avg[r] = np.mean(vals)
    
    targets = [data['W15']['scoring']['per_region'][r]['target_pct'] for r in RECOVERY_REGIONS]
    old_vals = [w_old_avg[r] for r in RECOVERY_REGIONS]
    new_vals = [data['W15']['scoring']['per_region'][r]['actual_pct'] for r in RECOVERY_REGIONS]
    
    x = np.arange(len(RECOVERY_REGIONS))
    width = 0.25
    
    ax.bar(x - width, targets, width, label='Target', color='black', alpha=0.3, edgecolor='black')
    ax.bar(x, old_vals, width, label='W01-W04 avg', color='#fc8d59', edgecolor='black', linewidth=0.5)
    ax.bar(x + width, new_vals, width, label='W15 (best)', color='#91bfdb', edgecolor='black', linewidth=0.5)
    
    ax.set_xticks(x)
    ax.set_xticklabels(RECOVERY_REGIONS, fontsize=8, rotation=45)
    ax.set_ylabel('Recovery (%)', fontsize=10)
    ax.set_title('A) Recovery: Before vs After', fontsize=11)
    ax.legend(fontsize=8)
    ax.set_ylim(0, 100)
    
    # --- Panel B: Timing ---
    ax = axes[1]
    
    # Subset of timing regions for clarity
    timing_subset = ['CA-S', 'CA-N', 'OR', 'WA-O', 'JDF', 'SS-S', 'BC-N', 'AK-FS', 'AK-PWS', 'AK-WG']
    
    targets_t = []
    old_t = []
    new_t = []
    labels_t = []
    
    for r in timing_subset:
        tr_new = data['W15']['arrival_timing']['per_region'][r]
        # Average old timing
        old_vals_t = []
        for w in range(1, 5):
            tr_old = data[f'W{w:02d}']['arrival_timing']['per_region'][r]
            act = tr_old['actual_months']
            if act is not None:
                old_vals_t.append(act)
        
        targets_t.append(tr_new['target_months'])
        old_t.append(np.mean(old_vals_t) if old_vals_t else np.nan)
        act_new = tr_new['actual_months']
        new_t.append(act_new if act_new is not None else np.nan)
        labels_t.append(r)
    
    x2 = np.arange(len(labels_t))
    
    ax.bar(x2 - width, targets_t, width, label='Target', color='black', alpha=0.3, edgecolor='black')
    
    # Old timing - replace NaN for display
    old_display = [v if not np.isnan(v) else 0 for v in old_t]
    old_colors = ['#fc8d59' if not np.isnan(v) else '#d9d9d9' for v in old_t]
    ax.bar(x2, old_display, width, label='W01-W04 avg', color=old_colors, edgecolor='black', linewidth=0.5)
    
    new_display = [v if not np.isnan(v) else 0 for v in new_t]
    new_colors = ['#91bfdb' if not np.isnan(v) else '#d9d9d9' for v in new_t]
    ax.bar(x2 + width, new_display, width, label='W15 (best)', color=new_colors, edgecolor='black', linewidth=0.5)
    
    # Mark unreached
    for i, r in enumerate(labels_t):
        if np.isnan(new_t[i]):
            ax.text(x2[i] + width, 2, '✗', ha='center', va='bottom', fontsize=12, color='red')
    
    ax.set_xticks(x2)
    ax.set_xticklabels(labels_t, fontsize=8, rotation=45)
    ax.set_ylabel('Arrival Time (months)', fontsize=10)
    ax.set_title('B) Timing: Before vs After', fontsize=11)
    ax.legend(fontsize=8)
    ax.set_ylim(0, 50)
    
    plt.tight_layout()
    fig.savefig(f'{FIG_DIR}/fig6_before_after.pdf', bbox_inches='tight')
    plt.close()
    print("Fig 6 saved")


if __name__ == '__main__':
    data = load_data()
    fig1_recovery_bar(data)
    fig2_rmse_bars(data)
    fig3_heatmap(data)
    fig4_recovery_gradient(data)
    fig5_timing(data)
    fig6_before_after(data)
    print("\nAll figures generated!")
