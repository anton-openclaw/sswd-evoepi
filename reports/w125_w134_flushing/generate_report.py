#!/usr/bin/env python3
"""Generate W125-W134 flushing analysis report with figures."""

import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

BASE = Path('/home/starbot/.openclaw/workspace/sswd-evoepi')
RESULTS = BASE / 'results' / 'calibration'
FIGURES = BASE / 'reports' / 'w125_w134_flushing' / 'figures'
FIGURES.mkdir(parents=True, exist_ok=True)

# CENTRALIZED: moved to sswd_evoepi.metrics
from sswd_evoepi.metrics import RECOVERY_TARGETS
TARGETS = RECOVERY_TARGETS
# TARGETS = {
#     'AK-PWS': 0.50, 'AK-FN': 0.50, 'AK-FS': 0.20, 'BC-N': 0.20,
#     'SS-S': 0.05, 'JDF': 0.02, 'OR': 0.0025, 'CA-N': 0.001
# }
TARGET_ORDER = list(TARGETS.keys())

CONFIGS = {
    'W125': {'label': 's₀=0.0005', 'param': 'settler_survival', 'val': 0.0005},
    'W126': {'label': 's₀=0.001 (base)', 'param': 'settler_survival', 'val': 0.001},
    'W127': {'label': 's₀=0.002', 'param': 'settler_survival', 'val': 0.002},
    'W128': {'label': 'α=0.15', 'param': 'alpha_env', 'val': 0.15},
    'W129': {'label': 'α=0.25', 'param': 'alpha_env', 'val': 0.25},
    'W130': {'label': 'α=0.30', 'param': 'alpha_env', 'val': 0.30},
    'W131': {'label': 'floor=300', 'param': 'P_env_floor', 'val': 300},
    'W132': {'label': 'floor=800', 'param': 'P_env_floor', 'val': 800},
}

def load_results():
    """Load combined_results.json for each run."""
    data = {}
    for run_id in CONFIGS:
        path = RESULTS / run_id / 'combined_results.json'
        if path.exists():
            data[run_id] = json.loads(path.read_text())
    # Also load W117 baseline
    w117_path = RESULTS / 'W117' / 'combined_results.json'
    if w117_path.exists():
        data['W117'] = json.loads(w117_path.read_text())
    return data

def load_monthly(run_id, seed):
    """Load monthly NPZ data."""
    path = RESULTS / run_id / f'monthly_seed{seed}.npz'
    if path.exists():
        return np.load(str(path), allow_pickle=True)
    return None

def get_regional_recovery(data_npz):
    """Compute per-region recovery fraction timeseries from monthly NPZ."""
    pops = data_npz['populations']  # (T, 896)
    site_names = data_npz['site_names']
    K = int(data_npz['K'])
    
    # Map sites to scored regions
    region_indices = {}
    for i, name in enumerate(site_names):
        name_str = str(name)
        # Extract region: everything before last dash-number
        parts = name_str.rsplit('-', 1)
        if len(parts) == 2:
            reg = parts[0]
            if reg in TARGETS:
                region_indices.setdefault(reg, []).append(i)
    
    # Monthly recovery = sum(pop_sites_in_region) / (n_sites * K)
    T = pops.shape[0]
    recovery = {}
    for reg, indices in region_indices.items():
        idx = np.array(indices)
        reg_pops = pops[:, idx].sum(axis=1)
        recovery[reg] = reg_pops / (len(idx) * K)
    
    return recovery

def get_seed_averaged_regional(run_data):
    """Get seed-averaged per-region actual values from scoring."""
    avgs = {reg: [] for reg in TARGET_ORDER}
    for r in run_data.get('results', []):
        scoring = r.get('scoring', {}).get('per_region', {})
        for reg in TARGET_ORDER:
            if reg in scoring:
                avgs[reg].append(scoring[reg].get('actual', 0))
    return {reg: np.mean(vals) if vals else 0 for reg, vals in avgs.items()}


# ============================================================
# FIGURE 1: RMSE Comparison Bar Chart
# ============================================================
def fig1_rmse_comparison(data):
    fig, ax = plt.subplots(figsize=(10, 5))
    
    runs = list(CONFIGS.keys())
    rmses = [data[r]['mean_rmse_log'] for r in runs]
    labels = [CONFIGS[r]['label'] for r in runs]
    colors = plt.cm.Set2(np.linspace(0, 1, len(runs)))
    
    bars = ax.bar(range(len(runs)), rmses, color=colors, edgecolor='black', linewidth=0.5)
    
    # Baseline lines
    if 'W117' in data:
        w117_rmse = data['W117']['mean_rmse_log']
        ax.axhline(w117_rmse, color='red', linestyle='--', linewidth=2, 
                    label=f'W117 baseline (φ=0.5): {w117_rmse:.3f}')
    
    # Best of batch
    best_idx = np.argmin(rmses)
    bars[best_idx].set_edgecolor('gold')
    bars[best_idx].set_linewidth(3)
    
    ax.set_xticks(range(len(runs)))
    ax.set_xticklabels(labels, rotation=30, ha='right')
    ax.set_ylabel('RMSE (log-scale)')
    ax.set_title('W125-W132: Per-Site Flushing Calibration\n(all worse than W117 uniform φ=0.5)')
    ax.legend(fontsize=11)
    
    # Annotate bars
    for i, (r, v) in enumerate(zip(runs, rmses)):
        ax.text(i, v + 0.01, f'{v:.3f}', ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    ax.set_ylim(0, max(rmses) * 1.15)
    plt.tight_layout()
    fig.savefig(FIGURES / 'fig1_rmse_comparison.png', dpi=150)
    plt.close()
    print("Figure 1: RMSE comparison ✓")


# ============================================================
# FIGURE 2: Regional Recovery Comparison — W129 vs W117
# ============================================================
def fig2_regional_comparison(data):
    fig, ax = plt.subplots(figsize=(12, 6))
    
    x = np.arange(len(TARGET_ORDER))
    width = 0.3
    
    # W129 (best with flushing)
    w129_vals = get_seed_averaged_regional(data['W129'])
    w129_bars = [w129_vals[r] * 100 for r in TARGET_ORDER]
    
    # W117 (best without flushing)
    w117_vals = get_seed_averaged_regional(data['W117'])
    w117_bars = [w117_vals[r] * 100 for r in TARGET_ORDER]
    
    # Targets
    target_bars = [TARGETS[r] * 100 for r in TARGET_ORDER]
    
    ax.bar(x - width, target_bars, width, label='Target', color='green', alpha=0.6, edgecolor='black')
    ax.bar(x, w117_bars, width, label='W117 (uniform φ=0.5)', color='steelblue', alpha=0.8, edgecolor='black')
    ax.bar(x + width, w129_bars, width, label='W129 (per-site φ)', color='coral', alpha=0.8, edgecolor='black')
    
    ax.set_xticks(x)
    ax.set_xticklabels(TARGET_ORDER, rotation=30, ha='right')
    ax.set_ylabel('Recovery fraction (%)')
    ax.set_title('Regional Recovery: Per-Site Flushing (W129) vs Uniform (W117)')
    ax.legend(fontsize=11)
    ax.set_yscale('symlog', linthresh=0.5)
    
    plt.tight_layout()
    fig.savefig(FIGURES / 'fig2_regional_comparison.png', dpi=150)
    plt.close()
    print("Figure 2: Regional comparison ✓")


# ============================================================
# FIGURE 3: The Mechanism — Flushing Rate vs Recovery Change
# ============================================================
def fig3_flushing_mechanism(data):
    """Show that per-site flushing rate correlates with recovery CHANGE from W117."""
    # Load enclosedness data for phi values
    enc_path = BASE / 'data' / 'nodes' / 'site_enclosedness.json'
    if not enc_path.exists():
        print("Figure 3: SKIP (no enclosedness data)")
        return
    
    enc_list = json.loads(enc_path.read_text())
    
    # Compute mean phi per scored region
    region_phi = {}
    for info in enc_list:
        site_id = info.get('name', '')
        parts = site_id.rsplit('-', 1)
        reg = parts[0] if len(parts) == 2 else site_id
        if reg in TARGETS:
            region_phi.setdefault(reg, []).append(info.get('flushing_rate', 0.5))
    region_phi_mean = {r: np.mean(v) for r, v in region_phi.items()}
    
    w129_vals = get_seed_averaged_regional(data['W129'])
    w117_vals = get_seed_averaged_regional(data['W117'])
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Left: phi vs W129 recovery
    for reg in TARGET_ORDER:
        phi = region_phi_mean.get(reg, 0.5)
        w129_r = w129_vals[reg] * 100
        w117_r = w117_vals[reg] * 100
        target_r = TARGETS[reg] * 100
        
        ax1.scatter(phi, w129_r, s=100, zorder=5, color='coral', edgecolor='black')
        ax1.annotate(reg, (phi, w129_r), textcoords="offset points", 
                     xytext=(5, 5), fontsize=9)
    
    ax1.set_xlabel('Mean flushing rate φ (per-site)')
    ax1.set_ylabel('Recovery fraction (%) — W129')
    ax1.set_title('Per-Site Flushing Rate vs Recovery')
    ax1.axvline(0.5, color='gray', linestyle=':', label='Old uniform φ=0.5')
    ax1.legend()
    
    # Right: phi vs CHANGE in recovery (W129 - W117)
    for reg in TARGET_ORDER:
        phi = region_phi_mean.get(reg, 0.5)
        delta = (w129_vals[reg] - w117_vals[reg]) * 100
        color = 'red' if delta > 0 else 'blue'
        ax2.scatter(phi, delta, s=100, zorder=5, color=color, edgecolor='black')
        ax2.annotate(reg, (phi, delta), textcoords="offset points", 
                     xytext=(5, 5), fontsize=9)
    
    ax2.axhline(0, color='black', linewidth=0.5)
    ax2.axvline(0.5, color='gray', linestyle=':', label='Old uniform φ=0.5')
    ax2.set_xlabel('Mean flushing rate φ (per-site)')
    ax2.set_ylabel('Δ Recovery (W129 - W117) in pp')
    ax2.set_title('How Flushing Changed Recovery\n(red=increased, blue=decreased)')
    ax2.legend()
    
    plt.tight_layout()
    fig.savefig(FIGURES / 'fig3_flushing_mechanism.png', dpi=150)
    plt.close()
    print("Figure 3: Flushing mechanism ✓")


# ============================================================
# FIGURE 4: Recovery Trajectories (W129 best run)
# ============================================================
def fig4_trajectories(data):
    seeds = [42, 123, 999]
    fig, axes = plt.subplots(3, 3, figsize=(15, 12))
    axes = axes.flatten()
    
    regions_to_plot = TARGET_ORDER + ['ALL']
    
    for idx, reg in enumerate(regions_to_plot):
        if idx >= 9:
            break
        ax = axes[idx]
        
        # W129 seeds
        for seed in seeds:
            npz = load_monthly('W129', seed)
            if npz is None:
                continue
            recovery = get_regional_recovery(npz)
            if reg == 'ALL':
                # Overall
                pops = npz['populations']
                K = int(npz['K'])
                total = pops.sum(axis=1) / (896 * K)
                months = np.arange(len(total))
                years = months / 12.0
                ax.plot(years, total * 100, color='coral', alpha=0.3, linewidth=1)
            elif reg in recovery:
                months = np.arange(len(recovery[reg]))
                years = months / 12.0
                ax.plot(years, recovery[reg] * 100, color='coral', alpha=0.3, linewidth=1)
        
        # W117 seeds for comparison
        for seed in seeds:
            npz = load_monthly('W117', seed)
            if npz is None:
                continue
            recovery = get_regional_recovery(npz)
            if reg == 'ALL':
                pops = npz['populations']
                K = int(npz['K'])
                total = pops.sum(axis=1) / (896 * K)
                months = np.arange(len(total))
                years = months / 12.0
                ax.plot(years, total * 100, color='steelblue', alpha=0.3, linewidth=1)
            elif reg in recovery:
                months = np.arange(len(recovery[reg]))
                years = months / 12.0
                ax.plot(years, recovery[reg] * 100, color='steelblue', alpha=0.3, linewidth=1)
        
        # Target line
        if reg in TARGETS:
            ax.axhline(TARGETS[reg] * 100, color='green', linestyle='--', linewidth=1.5, alpha=0.7)
        
        ax.set_title(reg, fontweight='bold')
        ax.set_xlabel('Year')
        ax.set_ylabel('Pop/K (%)')
        ax.set_ylim(bottom=0)
    
    # Legend in last subplot
    if len(regions_to_plot) < 9:
        axes[-1].axis('off')
    
    # Custom legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color='coral', linewidth=2, label='W129 (per-site φ)'),
        Line2D([0], [0], color='steelblue', linewidth=2, label='W117 (uniform φ=0.5)'),
        Line2D([0], [0], color='green', linestyle='--', linewidth=1.5, label='Target'),
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=3, fontsize=12,
               bbox_to_anchor=(0.5, -0.02))
    
    fig.suptitle('Recovery Trajectories: W129 (per-site flushing) vs W117 (uniform)', 
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig(FIGURES / 'fig4_trajectories.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Figure 4: Recovery trajectories ✓")


# ============================================================
# FIGURE 5: Regional Heatmap (all runs)
# ============================================================
def fig5_heatmap(data):
    runs = list(CONFIGS.keys())
    matrix = np.zeros((len(runs), len(TARGET_ORDER)))
    
    for i, run_id in enumerate(runs):
        vals = get_seed_averaged_regional(data[run_id])
        for j, reg in enumerate(TARGET_ORDER):
            matrix[i, j] = vals[reg] * 100
    
    fig, ax = plt.subplots(figsize=(12, 6))
    im = ax.imshow(matrix, cmap='RdYlGn', aspect='auto', vmin=0, vmax=30)
    
    # Annotate cells
    for i in range(len(runs)):
        for j in range(len(TARGET_ORDER)):
            val = matrix[i, j]
            target = TARGETS[TARGET_ORDER[j]] * 100
            # Color text based on whether close to target
            color = 'white' if val > 15 else 'black'
            ax.text(j, i, f'{val:.1f}%', ha='center', va='center', 
                    fontsize=9, color=color, fontweight='bold')
    
    # Target row
    ax.set_xticks(range(len(TARGET_ORDER)))
    ax.set_xticklabels([f'{r}\n(t={TARGETS[r]*100:.1f}%)' for r in TARGET_ORDER], fontsize=9)
    ax.set_yticks(range(len(runs)))
    ax.set_yticklabels([f'{r} ({CONFIGS[r]["label"]})' for r in runs])
    
    plt.colorbar(im, ax=ax, label='Recovery %')
    ax.set_title('Regional Recovery Heatmap (Year 13, seed-averaged)')
    plt.tight_layout()
    fig.savefig(FIGURES / 'fig5_heatmap.png', dpi=150)
    plt.close()
    print("Figure 5: Heatmap ✓")


# ============================================================
# FIGURE 6: Diagnosis — Where RMSE Gets Worse
# ============================================================
def fig6_diagnosis(data):
    """Per-region log-squared-error comparison W129 vs W117."""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    x = np.arange(len(TARGET_ORDER))
    width = 0.35
    
    # Compute per-region mean log_sq_error
    def get_lse(run_data):
        lse = {reg: [] for reg in TARGET_ORDER}
        for r in run_data.get('results', []):
            scoring = r.get('scoring', {}).get('per_region', {})
            for reg in TARGET_ORDER:
                if reg in scoring:
                    lse[reg].append(scoring[reg].get('log_sq_error', 0))
        return {reg: np.mean(vals) if vals else 0 for reg, vals in lse.items()}
    
    w117_lse = get_lse(data['W117'])
    w129_lse = get_lse(data['W129'])
    
    w117_bars = [w117_lse[r] for r in TARGET_ORDER]
    w129_bars = [w129_lse[r] for r in TARGET_ORDER]
    
    ax.bar(x - width/2, w117_bars, width, label='W117 (uniform φ=0.5)', color='steelblue', alpha=0.8, edgecolor='black')
    ax.bar(x + width/2, w129_bars, width, label='W129 (per-site φ)', color='coral', alpha=0.8, edgecolor='black')
    
    # Mark where W129 is worse
    for i in range(len(TARGET_ORDER)):
        if w129_bars[i] > w117_bars[i]:
            ax.annotate('↑WORSE', (i + width/2, w129_bars[i]), 
                        textcoords="offset points", xytext=(0, 5),
                        ha='center', fontsize=8, color='red', fontweight='bold')
        else:
            ax.annotate('↓better', (i + width/2, w129_bars[i]), 
                        textcoords="offset points", xytext=(0, 5),
                        ha='center', fontsize=8, color='green')
    
    ax.set_xticks(x)
    ax.set_xticklabels(TARGET_ORDER, rotation=30, ha='right')
    ax.set_ylabel('Log-squared error (lower = better)')
    ax.set_title('Per-Region Error: Where Per-Site Flushing Hurts')
    ax.legend(fontsize=11)
    plt.tight_layout()
    fig.savefig(FIGURES / 'fig6_diagnosis.png', dpi=150)
    plt.close()
    print("Figure 6: Diagnosis ✓")


# ============================================================
# FIGURE 7: Alpha_env sensitivity with flushing
# ============================================================
def fig7_alpha_sensitivity(data):
    fig, ax = plt.subplots(figsize=(8, 5))
    
    # W126=baseline(alpha=0.2), W128=0.15, W129=0.25, W130=0.3
    alphas = [0.15, 0.20, 0.25, 0.30]
    run_ids = ['W128', 'W126', 'W129', 'W130']
    rmses = [data[r]['mean_rmse_log'] for r in run_ids]
    
    ax.plot(alphas, rmses, 'o-', color='coral', linewidth=2, markersize=10, markeredgecolor='black')
    
    for a, rmse, rid in zip(alphas, rmses, run_ids):
        ax.annotate(f'{rid}\n{rmse:.3f}', (a, rmse), textcoords="offset points",
                    xytext=(0, 12), ha='center', fontsize=10, fontweight='bold')
    
    ax.axhline(0.603, color='steelblue', linestyle='--', linewidth=2,
               label='W117 baseline (uniform φ, RMSE=0.603)')
    
    ax.set_xlabel('α_env (host amplification rate)')
    ax.set_ylabel('RMSE')
    ax.set_title('Alpha_env Sensitivity Under Per-Site Flushing')
    ax.legend()
    ax.set_xlim(0.12, 0.33)
    plt.tight_layout()
    fig.savefig(FIGURES / 'fig7_alpha_sensitivity.png', dpi=150)
    plt.close()
    print("Figure 7: Alpha sensitivity ✓")


def main():
    data = load_results()
    print(f"Loaded {len(data)} runs: {list(data.keys())}")
    
    if 'W117' not in data:
        print("WARNING: W117 baseline not available for comparison!")
    
    fig1_rmse_comparison(data)
    fig2_regional_comparison(data)
    fig3_flushing_mechanism(data)
    fig4_trajectories(data)
    fig5_heatmap(data)
    fig6_diagnosis(data)
    fig7_alpha_sensitivity(data)
    
    print("\nAll figures generated in:", FIGURES)


if __name__ == '__main__':
    main()
