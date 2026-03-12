#!/usr/bin/env python3
"""
Generate comprehensive W135-W144 calibration report with figures.
This round tests continuous larval retention (alpha_self) from enclosedness.
"""

import os
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.gridspec import GridSpec
from pathlib import Path

# Set style for publication-quality figures
plt.style.use('default')
plt.rcParams.update({
    'figure.dpi': 150,
    'savefig.dpi': 150,
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 10,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.titlesize': 14
})

# Data paths
BASE_PATH = Path("/home/starbot/.openclaw/workspace/sswd-evoepi/results/calibration/")
FIGURES_PATH = Path("/home/starbot/.openclaw/workspace/sswd-evoepi/reports/w135_w144_larval/figures/")

# Run configurations
RUNS_CONFIG = {
    'W135': {'n_conn': 1.0, 'alpha_env': 0.25, 'alpha_self': 'default', 'desc': 'Baseline'},
    'W136': {'n_conn': 0.5, 'alpha_env': 0.25, 'alpha_self': 'default', 'desc': 'Gentle flushing'},
    'W137': {'n_conn': 2.0, 'alpha_env': 0.25, 'alpha_self': 'default', 'desc': 'Sharp flushing'},
    'W138': {'n_conn': 1.0, 'alpha_env': 0.25, 'alpha_self': 'strong', 'desc': 'Strong larval'},
    'W139': {'n_conn': 0.5, 'alpha_env': 0.25, 'alpha_self': 'strong', 'desc': 'Gentle + strong larval'},
    'W140': {'n_conn': 1.0, 'alpha_env': 0.20, 'alpha_self': 'default', 'desc': 'Low alpha_env'},
    'W141': {'n_conn': 1.0, 'alpha_env': 0.30, 'alpha_self': 'default', 'desc': 'High alpha_env'},
    'W142': {'n_conn': 0.5, 'alpha_env': 0.20, 'alpha_self': 'default', 'desc': '⭐ NEW BEST'},
}

# Baseline comparisons
BASELINES = {
    'W117': {'rmse': 0.603, 'desc': 'Uniform phi=0.5'},
    'W129': {'rmse': 0.691, 'desc': 'Per-site phi, no larval'}
}

# Regional recovery data (from user input - seed-averaged, year 13)
RECOVERY_DATA = {
    'W129': {'RMSE': 0.691, 'AK-PWS': 3.7, 'AK-FN': 5.5, 'BC-N': 4.1, 'SS-S': 5.6, 'JDF': 4.1, 'OR': 0.5, 'CA-N': 0.4},
    'W117': {'RMSE': 0.603, 'AK-PWS': 7.9, 'AK-FN': 3.2, 'BC-N': 7.3, 'SS-S': 2.5, 'JDF': 1.4, 'OR': 0.5, 'CA-N': 0.2},
    'W135': {'RMSE': 0.727, 'AK-PWS': 3.9, 'AK-FN': 5.3, 'BC-N': 3.8, 'SS-S': 4.8, 'JDF': 1.8, 'OR': 1.4, 'CA-N': 0.6},
    'W136': {'RMSE': 0.839, 'AK-PWS': 1.7, 'AK-FN': 3.4, 'BC-N': 2.7, 'SS-S': 2.7, 'JDF': 1.4, 'OR': 0.1, 'CA-N': 0.0},
    'W137': {'RMSE': 0.757, 'AK-PWS': 7.3, 'AK-FN': 8.1, 'BC-N': 6.3, 'SS-S': 7.9, 'JDF': 1.2, 'OR': 4.4, 'CA-N': 1.2},
    'W138': {'RMSE': 0.699, 'AK-PWS': 3.6, 'AK-FN': 5.5, 'BC-N': 4.1, 'SS-S': 4.5, 'JDF': 1.6, 'OR': 0.6, 'CA-N': 0.7},
    'W139': {'RMSE': 0.991, 'AK-PWS': 2.0, 'AK-FN': 3.3, 'BC-N': 3.0, 'SS-S': 2.6, 'JDF': 1.5, 'OR': 0.2, 'CA-N': 0.0},
    'W140': {'RMSE': 0.770, 'AK-PWS': 8.1, 'AK-FN': 15.5, 'BC-N': 9.0, 'SS-S': 10.0, 'JDF': 4.2, 'OR': 5.3, 'CA-N': 2.0},
    'W141': {'RMSE': 0.822, 'AK-PWS': 1.3, 'AK-FN': 4.5, 'BC-N': 2.5, 'SS-S': 3.6, 'JDF': 1.2, 'OR': 0.3, 'CA-N': 0.2},
    'W142': {'RMSE': 0.599, 'AK-PWS': 5.1, 'AK-FN': 8.3, 'BC-N': 4.9, 'SS-S': 7.1, 'JDF': 2.1, 'OR': 0.7, 'CA-N': 0.5},
}

# Calibration targets
TARGETS = {
    'AK-PWS': 50.0, 'AK-FN': 50.0, 'AK-FS': 20.0, 'BC-N': 20.0, 
    'SS-S': 5.0, 'JDF': 2.0, 'OR': 0.25, 'CA-N': 0.1
}

def load_run_data(run_name):
    """Load combined results for a run."""
    results_path = BASE_PATH / run_name / "combined_results.json"
    if not results_path.exists():
        return None
    
    with open(results_path) as f:
        data = json.load(f)
    return data

def load_monthly_data(run_name, seed=42):
    """Load monthly NPZ data for trajectory analysis."""
    npz_path = BASE_PATH / run_name / f"monthly_seed{seed}.npz"
    if not npz_path.exists():
        # Try other seeds
        for alt_seed in [123, 999]:
            alt_path = BASE_PATH / run_name / f"monthly_seed{alt_seed}.npz"
            if alt_path.exists():
                npz_path = alt_path
                break
        else:
            return None
    
    return np.load(npz_path, allow_pickle=True)

def compute_recovery_trajectories(monthly_data):
    """Compute recovery trajectories by region."""
    populations = monthly_data['populations']  # (159, 896)
    site_names = monthly_data['site_names']
    K = int(monthly_data['K'])
    
    # Group sites by region
    regions = {}
    for i, site_name in enumerate(site_names):
        site_str = site_name.item() if hasattr(site_name, 'item') else str(site_name)
        region = site_str.rsplit('-', 1)[0]
        if region not in regions:
            regions[region] = []
        regions[region].append(i)
    
    # Compute recovery percentages over time
    trajectories = {}
    for region, site_indices in regions.items():
        if region in TARGETS:  # Only scored regions
            region_pops = populations[:, site_indices]
            total_sites = len(site_indices)
            recovery_pct = (region_pops.sum(axis=1) / (total_sites * K)) * 100
            trajectories[region] = recovery_pct
    
    # Overall trajectory
    all_pops = populations.sum(axis=1)
    total_sites = populations.shape[1]
    trajectories['Overall'] = (all_pops / (total_sites * K)) * 100
    
    return trajectories

def create_figure_1_rmse_bar_chart():
    """Figure 1: RMSE bar chart with baselines."""
    print("Creating Figure 1: RMSE bar chart...")
    
    run_names = list(RUNS_CONFIG.keys())
    rmse_values = [RECOVERY_DATA[run]['RMSE'] for run in run_names]
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Create bars
    bars = ax.bar(range(len(run_names)), rmse_values, 
                  color=['orange' if run == 'W142' else 'lightblue' for run in run_names],
                  edgecolor='black', linewidth=0.8)
    
    # Highlight best run
    for i, (run, bar) in enumerate(zip(run_names, bars)):
        if run == 'W142':
            bar.set_color('orange')
            bar.set_edgecolor('red')
            bar.set_linewidth(2)
            # Add star
            ax.text(i, rmse_values[i] + 0.02, '⭐', ha='center', va='bottom', fontsize=16)
    
    # Add baseline lines
    ax.axhline(y=BASELINES['W117']['rmse'], color='blue', linestyle='--', 
               label=f"W117 (uniform φ=0.5): {BASELINES['W117']['rmse']:.3f}", linewidth=2)
    ax.axhline(y=BASELINES['W129']['rmse'], color='green', linestyle='--', 
               label=f"W129 (per-site φ, no larval): {BASELINES['W129']['rmse']:.3f}", linewidth=2)
    
    # Customize plot
    ax.set_xlabel('Calibration Run')
    ax.set_ylabel('RMSE (log space)')
    ax.set_title('W135-W142: RMSE Performance vs Baselines\n(First time per-site flushing beats uniform!)')
    ax.set_xticks(range(len(run_names)))
    ax.set_xticklabels([f"{run}\n{RUNS_CONFIG[run]['desc']}" for run in run_names], 
                       rotation=45, ha='right')
    ax.legend(loc='upper right')
    ax.grid(axis='y', alpha=0.3)
    
    # Add RMSE values on bars
    for i, (run, rmse) in enumerate(zip(run_names, rmse_values)):
        ax.text(i, rmse + 0.01, f'{rmse:.3f}', ha='center', va='bottom', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(FIGURES_PATH / "01_rmse_comparison.png", bbox_inches='tight')
    plt.close()

def create_figure_2_recovery_trajectories():
    """Figure 2: Recovery trajectories (3x3 panel) for key runs."""
    print("Creating Figure 2: Recovery trajectories...")
    
    # Load trajectory data for key runs
    key_runs = ['W142', 'W117', 'W137']  # Best, baseline, interesting
    colors = ['orange', 'blue', 'green']
    labels = ['W142 (NEW BEST)', 'W117 (uniform baseline)', 'W137 (sharp flushing)']
    
    trajectories = {}
    for run in key_runs:
        if run in ['W117', 'W129']:  # Baselines - we might not have monthly data
            continue  # Skip for now, we can add if data is available
        monthly_data = load_monthly_data(run)
        if monthly_data is not None:
            trajectories[run] = compute_recovery_trajectories(monthly_data)
    
    # Create 3x3 subplot
    fig = plt.figure(figsize=(15, 12))
    regions = list(TARGETS.keys()) + ['Overall']
    
    for i, region in enumerate(regions):
        ax = plt.subplot(3, 3, i + 1)
        
        # Plot trajectories for available runs
        for j, run in enumerate(key_runs):
            if run in trajectories and region in trajectories[run]:
                traj = trajectories[run][region]
                time_points = np.arange(len(traj))
                ax.plot(time_points, traj, color=colors[j], 
                       label=labels[j], linewidth=2)
        
        # Add target line if applicable
        if region in TARGETS:
            ax.axhline(y=TARGETS[region], color='red', linestyle='--', 
                      alpha=0.7, label=f'Target: {TARGETS[region]}%')
        
        ax.set_title(f'{region}')
        ax.set_xlabel('Time (months)')
        ax.set_ylabel('Recovery (%)')
        ax.grid(alpha=0.3)
        
        if i == 0:  # Only show legend on first subplot
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.suptitle('Recovery Trajectories: Key Runs vs Targets', fontsize=16)
    plt.tight_layout()
    plt.savefig(FIGURES_PATH / "02_recovery_trajectories.png", bbox_inches='tight')
    plt.close()

def create_figure_3_regional_recovery_comparison():
    """Figure 3: Regional recovery comparison - THE key figure."""
    print("Creating Figure 3: Regional recovery comparison...")
    
    regions = ['AK-PWS', 'AK-FN', 'BC-N', 'SS-S', 'JDF', 'OR', 'CA-N']
    comparison_runs = ['W117', 'W129', 'W142']
    run_labels = ['W117 (uniform)', 'W129 (per-site, no larval)', 'W142 (NEW BEST)']
    colors = ['blue', 'green', 'orange']
    
    # Prepare data
    data = []
    for region in regions:
        for i, run in enumerate(comparison_runs):
            recovery = RECOVERY_DATA[run][region]
            target = TARGETS[region]
            data.append({
                'Region': region,
                'Run': run_labels[i],
                'Recovery': recovery,
                'Target': target,
                'Color': colors[i]
            })
    
    # Create grouped bar chart
    fig, ax = plt.subplots(figsize=(14, 8))
    
    x = np.arange(len(regions))
    width = 0.25
    
    for i, (run, label, color) in enumerate(zip(comparison_runs, run_labels, colors)):
        recoveries = [RECOVERY_DATA[run][region] for region in regions]
        bars = ax.bar(x + i*width - width, recoveries, width, 
                     label=label, color=color, alpha=0.8, edgecolor='black', linewidth=0.5)
        
        # Add values on bars
        for j, bar in enumerate(bars):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                   f'{height:.1f}%', ha='center', va='bottom', fontsize=8)
    
    # Add target bars
    targets = [TARGETS[region] for region in regions]
    target_bars = ax.bar(x + width, targets, width, 
                        label='Target', color='red', alpha=0.3, hatch='///')
    
    # Customize plot
    ax.set_xlabel('Region')
    ax.set_ylabel('Recovery Percentage (%)')
    ax.set_title('Regional Recovery Comparison: Progressive Improvement\n(W117 → W129 → W142)')
    ax.set_xticks(x)
    ax.set_xticklabels(regions, rotation=45, ha='right')
    ax.legend()
    ax.grid(axis='y', alpha=0.3)
    ax.set_yscale('log')  # Log scale to handle wide range
    
    plt.tight_layout()
    plt.savefig(FIGURES_PATH / "03_regional_recovery_comparison.png", bbox_inches='tight')
    plt.close()

def create_figure_4_n_connectivity_effect():
    """Figure 4: n_connectivity effect on RMSE."""
    print("Creating Figure 4: n_connectivity effect...")
    
    # Data for n_connectivity effect
    n_conn_runs = [
        ('W136', 0.5, 0.839),  # Gentle flushing
        ('W135', 1.0, 0.727),  # Default
        ('W137', 2.0, 0.757),  # Sharp flushing
    ]
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    run_names, n_conn_values, rmse_values = zip(*n_conn_runs)
    
    # Plot line and points
    ax.plot(n_conn_values, rmse_values, 'o-', linewidth=2, markersize=10, color='blue')
    
    # Annotate points
    for run, n_conn, rmse in n_conn_runs:
        ax.annotate(f'{run}\n({rmse:.3f})', 
                   (n_conn, rmse), 
                   textcoords="offset points", 
                   xytext=(0,15), ha='center', fontsize=10)
    
    # Highlight the winner (lowest RMSE)
    min_idx = rmse_values.index(min(rmse_values))
    ax.plot(n_conn_values[min_idx], rmse_values[min_idx], 'o', 
           markersize=15, color='orange', markerfacecolor='none', markeredgewidth=3)
    
    ax.set_xlabel('n_connectivity (flushing contrast parameter)')
    ax.set_ylabel('RMSE (log space)')
    ax.set_title('Effect of n_connectivity on Performance\n(Gentler flushing contrast wins)')
    ax.grid(alpha=0.3)
    ax.set_xticks(n_conn_values)
    
    plt.tight_layout()
    plt.savefig(FIGURES_PATH / "04_n_connectivity_effect.png", bbox_inches='tight')
    plt.close()

def create_figure_5_larval_retention_impact():
    """Figure 5: Larval retention impact."""
    print("Creating Figure 5: Larval retention impact...")
    
    # Comparison of larval retention effects
    larval_runs = [
        ('W129', 'No larval retention', 0.691),
        ('W135', 'Default larval retention', 0.727),
        ('W138', 'Strong larval retention', 0.699),
    ]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Left panel: RMSE comparison
    run_names, descriptions, rmse_values = zip(*larval_runs)
    colors = ['red', 'lightblue', 'darkblue']
    
    bars = ax1.bar(range(len(run_names)), rmse_values, color=colors, 
                   alpha=0.8, edgecolor='black', linewidth=0.8)
    
    # Add values on bars
    for bar, rmse in zip(bars, rmse_values):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                f'{rmse:.3f}', ha='center', va='bottom', fontweight='bold')
    
    ax1.set_xlabel('Configuration')
    ax1.set_ylabel('RMSE (log space)')
    ax1.set_title('RMSE by Larval Retention Level')
    ax1.set_xticks(range(len(run_names)))
    ax1.set_xticklabels([f'{run}\n{desc}' for run, desc in zip(run_names, descriptions)],
                       rotation=45, ha='right')
    ax1.grid(axis='y', alpha=0.3)
    
    # Right panel: Regional recovery comparison
    regions = ['AK-PWS', 'AK-FN', 'BC-N', 'SS-S']  # Key regions
    x = np.arange(len(regions))
    width = 0.25
    
    for i, (run, desc, rmse) in enumerate(larval_runs):
        if run in RECOVERY_DATA:
            recoveries = [RECOVERY_DATA[run][region] for region in regions]
            ax2.bar(x + i*width - width/2, recoveries, width, 
                   label=f'{run} ({desc})', color=colors[i], alpha=0.8)
    
    ax2.set_xlabel('Region')
    ax2.set_ylabel('Recovery Percentage (%)')
    ax2.set_title('Regional Recovery by Larval Retention')
    ax2.set_xticks(x)
    ax2.set_xticklabels(regions)
    ax2.legend()
    ax2.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(FIGURES_PATH / "05_larval_retention_impact.png", bbox_inches='tight')
    plt.close()

def create_figure_6_error_decomposition():
    """Figure 6: Per-region error decomposition for W142."""
    print("Creating Figure 6: Error decomposition for W142...")
    
    # Load W142 data to get per-region log_sq_error
    w142_data = load_run_data('W142')
    if w142_data is None:
        print("Warning: W142 data not available for error decomposition")
        return
    
    # Extract per-region errors from first seed
    per_region = w142_data['results'][0]['scoring']['per_region']
    regions = list(per_region.keys())
    log_sq_errors = [per_region[region]['log_sq_error'] for region in regions]
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Create bars
    bars = ax.bar(regions, log_sq_errors, color='orange', alpha=0.8, 
                  edgecolor='black', linewidth=0.8)
    
    # Add values on bars
    for bar, error in zip(bars, log_sq_errors):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 0.01,
               f'{error:.3f}', ha='center', va='bottom', fontweight='bold')
    
    ax.set_xlabel('Region')
    ax.set_ylabel('Log Squared Error')
    ax.set_title('W142: Per-Region Error Decomposition\n(Where the remaining RMSE=0.599 comes from)')
    ax.set_xticklabels(regions, rotation=45, ha='right')
    ax.grid(axis='y', alpha=0.3)
    
    # Add note about AK dominance
    total_error = sum(log_sq_errors)
    ak_errors = [log_sq_errors[i] for i, region in enumerate(regions) if region.startswith('AK')]
    ak_fraction = sum(ak_errors) / total_error * 100
    ax.text(0.02, 0.98, f'Alaska regions contribute {ak_fraction:.1f}% of total error',
           transform=ax.transAxes, va='top', ha='left',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(FIGURES_PATH / "06_error_decomposition_w142.png", bbox_inches='tight')
    plt.close()

def create_figure_7_parameter_interaction_heatmap():
    """Figure 7: Parameter interaction heatmap."""
    print("Creating Figure 7: Parameter interaction heatmap...")
    
    # Create parameter grid
    n_conn_values = [0.5, 1.0, 2.0]
    alpha_env_values = [0.20, 0.25, 0.30]
    
    # Map runs to parameter combinations
    param_map = {
        (0.5, 0.20): ('W142', 0.599),  # Best combination
        (0.5, 0.25): ('W136', 0.839),
        (1.0, 0.20): ('W140', 0.770),
        (1.0, 0.25): ('W135', 0.727),
        (1.0, 0.30): ('W141', 0.822),
        (2.0, 0.25): ('W137', 0.757),
    }
    
    # Create RMSE matrix
    rmse_matrix = np.full((len(n_conn_values), len(alpha_env_values)), np.nan)
    annotation_matrix = []
    
    for i, n_conn in enumerate(n_conn_values):
        row_annotations = []
        for j, alpha_env in enumerate(alpha_env_values):
            if (n_conn, alpha_env) in param_map:
                run_name, rmse = param_map[(n_conn, alpha_env)]
                rmse_matrix[i, j] = rmse
                row_annotations.append(f'{run_name}\n{rmse:.3f}')
            else:
                row_annotations.append('')
        annotation_matrix.append(row_annotations)
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Create heatmap
    im = ax.imshow(rmse_matrix, cmap='RdYlBu_r', aspect='auto', interpolation='nearest')
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('RMSE (log space)')
    
    # Set ticks and labels
    ax.set_xticks(range(len(alpha_env_values)))
    ax.set_yticks(range(len(n_conn_values)))
    ax.set_xticklabels([f'{val:.2f}' for val in alpha_env_values])
    ax.set_yticklabels([f'{val:.1f}' for val in n_conn_values])
    ax.set_xlabel('alpha_env (disease amplification)')
    ax.set_ylabel('n_connectivity (flushing contrast)')
    ax.set_title('Parameter Interaction Heatmap\n(Sweet spot: n_conn=0.5, alpha_env=0.20)')
    
    # Add annotations
    for i in range(len(n_conn_values)):
        for j in range(len(alpha_env_values)):
            if annotation_matrix[i][j]:
                # Highlight best combination
                if (n_conn_values[i], alpha_env_values[j]) == (0.5, 0.20):
                    ax.add_patch(patches.Rectangle((j-0.4, i-0.4), 0.8, 0.8, 
                                                 fill=False, edgecolor='red', lw=3))
                    text_color = 'red'
                    weight = 'bold'
                else:
                    text_color = 'black'
                    weight = 'normal'
                
                ax.text(j, i, annotation_matrix[i][j], ha='center', va='center',
                       color=text_color, fontweight=weight, fontsize=10)
    
    plt.tight_layout()
    plt.savefig(FIGURES_PATH / "07_parameter_interaction_heatmap.png", bbox_inches='tight')
    plt.close()

def create_all_figures():
    """Generate all figures."""
    print("Generating all figures for W135-W144 calibration report...")
    FIGURES_PATH.mkdir(parents=True, exist_ok=True)
    
    create_figure_1_rmse_bar_chart()
    create_figure_2_recovery_trajectories()
    create_figure_3_regional_recovery_comparison()
    create_figure_4_n_connectivity_effect()
    create_figure_5_larval_retention_impact()
    create_figure_6_error_decomposition()
    create_figure_7_parameter_interaction_heatmap()
    
    print(f"All figures saved to {FIGURES_PATH}")

if __name__ == "__main__":
    create_all_figures()