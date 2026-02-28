#!/usr/bin/env python3
"""
SSWD-EvoEpi Calibration Analysis Script
Generates figures and analysis for the first calibration batch.
"""

import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Configuration
DATA_DIR = Path("/home/starbot/.openclaw/workspace/sswd-evoepi/data/calibration")
OUTPUT_DIR = Path("/home/starbot/.openclaw/workspace/sswd-evoepi/reports/calibration_r1/figures")
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

# Target regions and their expert estimates
TARGETS = {
    'AK-PWS': 0.50,  # 50%
    'AK-FN': 0.50,   # 50%
    'AK-FS': 0.20,   # 20%
    'BC-N': 0.20,    # 20%
    'SS-S': 0.05,    # 5%
    'JDF': 0.02,     # 2%
    'OR': 0.0025,    # 0.25%
    'CA-N': 0.001    # 0.1%
}

REGION_ORDER = list(TARGETS.keys())

# Round descriptions
ROUND_DESCRIPTIONS = {
    0: "Baseline (defaults)",
    1: "K_half=30,000",
    2: "K_half=200,000", 
    3: "P_env_max=100",
    4: "P_env_max=2,000",
    5: "Combo (K_half=200K + P_env_max=100 + a_exposure=0.3)"
}

def load_results():
    """Load all available results from JSON files."""
    results = {}
    
    for round_num in range(6):  # Rounds 0-5
        round_dir = DATA_DIR / f"round_{round_num:02d}"
        if not round_dir.exists():
            continue
            
        round_results = {}
        seeds = [42, 123, 999]
        
        for seed in seeds:
            result_file = round_dir / f"result_seed{seed}.json"
            if result_file.exists():
                try:
                    with open(result_file, 'r') as f:
                        data = json.load(f)
                    
                    # Check if run completed successfully
                    if 'early_stop' in data:
                        print(f"Round {round_num:02d}, Seed {seed}: Early stop - {data['early_stop']}")
                        continue
                        
                    # Check if we have regional data
                    if not data.get('scoring', {}).get('per_region', {}):
                        print(f"Round {round_num:02d}, Seed {seed}: No regional scoring data")
                        continue
                        
                    round_results[seed] = data
                    
                except Exception as e:
                    print(f"Error loading round {round_num:02d}, seed {seed}: {e}")
                    continue
        
        if round_results:
            results[round_num] = round_results
            print(f"Loaded Round {round_num:02d}: {len(round_results)} seeds")
    
    return results

def create_summary_df(results):
    """Create a summary DataFrame from results."""
    rows = []
    
    for round_num, round_data in results.items():
        for seed, data in round_data.items():
            row = {
                'round': round_num,
                'seed': seed,
                'description': ROUND_DESCRIPTIONS.get(round_num, f"Round {round_num}"),
                'rmse_log': data['scoring']['rmse_log'],
                'within_2x': data['scoring']['within_2x'],
                'within_5x': data['scoring']['within_5x'],
                'wall_time_hours': data['wall_time_seconds'] / 3600
            }
            
            # Add regional data
            for region in REGION_ORDER:
                if region in data['scoring']['per_region']:
                    row[f'{region}_target'] = data['scoring']['per_region'][region]['target']
                    row[f'{region}_actual'] = data['scoring']['per_region'][region]['actual']
                    row[f'{region}_actual_pct'] = data['scoring']['per_region'][region]['actual_pct']
                    row[f'{region}_log_error'] = data['scoring']['per_region'][region]['log_error']
                    row[f'{region}_within_2x'] = data['scoring']['per_region'][region]['within_2x']
                    row[f'{region}_within_5x'] = data['scoring']['per_region'][region]['within_5x']
                else:
                    row[f'{region}_target'] = TARGETS[region]
                    row[f'{region}_actual'] = np.nan
                    row[f'{region}_actual_pct'] = np.nan
                    row[f'{region}_log_error'] = np.nan
                    row[f'{region}_within_2x'] = False
                    row[f'{region}_within_5x'] = False
            
            rows.append(row)
    
    return pd.DataFrame(rows)

def plot_grouped_bar_chart(df, output_path):
    """Create grouped bar chart of target vs actual recovery."""
    # Calculate means across seeds
    mean_data = []
    
    for round_num in sorted(df['round'].unique()):
        round_data = df[df['round'] == round_num]
        row = {'round': round_num, 'description': ROUND_DESCRIPTIONS.get(round_num, f"Round {round_num}")}
        
        for region in REGION_ORDER:
            target_col = f'{region}_target'
            actual_col = f'{region}_actual_pct'
            
            if target_col in round_data.columns and not round_data[target_col].isna().all():
                row[f'{region}_target'] = round_data[target_col].iloc[0]  # Target is constant
                row[f'{region}_actual'] = round_data[actual_col].mean() if not round_data[actual_col].isna().all() else np.nan
            else:
                row[f'{region}_target'] = TARGETS[region] * 100
                row[f'{region}_actual'] = np.nan
        
        mean_data.append(row)
    
    mean_df = pd.DataFrame(mean_data)
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(14, 8))
    
    x = np.arange(len(REGION_ORDER))
    width = 0.12
    
    colors = plt.cm.Set3(np.linspace(0, 1, len(mean_df)))
    
    # Plot targets once
    target_values = [TARGETS[region] * 100 for region in REGION_ORDER]
    ax.bar(x - width * len(mean_df) / 2, target_values, width, 
           label='Target', color='black', alpha=0.7, hatch='///')
    
    # Plot actual values for each round
    for i, (_, row) in enumerate(mean_df.iterrows()):
        actual_values = [row[f'{region}_actual'] for region in REGION_ORDER]
        offset = width * (i - len(mean_df) / 2 + 0.5)
        
        ax.bar(x + offset, actual_values, width, 
               label=f"R{row['round']:02d}: {row['description'][:20]}", 
               color=colors[i], alpha=0.8)
    
    ax.set_xlabel('Region')
    ax.set_ylabel('Recovery Percentage (%)')
    ax.set_title('Target vs Actual Recovery by Region and Round\n(Seed-averaged)')
    ax.set_yscale('log')
    ax.set_ylim(0.01, 100)
    ax.set_xticks(x)
    ax.set_xticklabels(REGION_ORDER, rotation=45)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def plot_heatmap(df, output_path):
    """Create heatmap of log10(actual/target) ratios."""
    # Calculate mean ratios across seeds
    heatmap_data = []
    
    for round_num in sorted(df['round'].unique()):
        round_data = df[df['round'] == round_num]
        row_data = []
        
        for region in REGION_ORDER:
            target_col = f'{region}_target'
            actual_col = f'{region}_actual'
            
            if target_col in round_data.columns and not round_data[actual_col].isna().all():
                target = round_data[target_col].iloc[0]
                actual = round_data[actual_col].mean()
                
                if not np.isnan(actual) and target > 0 and actual > 0:
                    ratio = np.log10(actual / target)
                else:
                    ratio = np.nan
            else:
                ratio = np.nan
            
            row_data.append(ratio)
        
        heatmap_data.append(row_data)
    
    # Create the heatmap
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Create custom colormap (blue-white-red)
    colors = ['#2166ac', '#4393c3', '#92c5de', '#d1e5f0', 
              'white', 
              '#fdbf6f', '#ff7f00', '#e31a1c', '#b10026']
    n_bins = 256
    cmap = LinearSegmentedColormap.from_list('custom', colors, N=n_bins)
    
    # Plot heatmap
    im = ax.imshow(heatmap_data, cmap=cmap, aspect='auto', vmin=-2, vmax=2)
    
    # Add text annotations
    for i in range(len(heatmap_data)):
        for j in range(len(REGION_ORDER)):
            if not np.isnan(heatmap_data[i][j]):
                text = f'{heatmap_data[i][j]:.2f}'
                ax.text(j, i, text, ha="center", va="center", fontweight='bold')
    
    ax.set_xticks(range(len(REGION_ORDER)))
    ax.set_xticklabels(REGION_ORDER)
    ax.set_yticks(range(len(df['round'].unique())))
    ax.set_yticklabels([f"R{r:02d}" for r in sorted(df['round'].unique())])
    ax.set_xlabel('Region')
    ax.set_ylabel('Round')
    ax.set_title('Model Performance Heatmap\nlog₁₀(Actual/Target)\nGreen=Match, Red=Overshoot, Blue=Undershoot')
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('log₁₀(Actual/Target)')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def plot_rmse_comparison(df, output_path):
    """Create RMSE comparison bar chart."""
    # Calculate RMSE stats per round
    rmse_stats = df.groupby(['round', 'description'])['rmse_log'].agg(['mean', 'std']).reset_index()
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    x = np.arange(len(rmse_stats))
    bars = ax.bar(x, rmse_stats['mean'], yerr=rmse_stats['std'], 
                  capsize=5, color=plt.cm.Set2(np.arange(len(rmse_stats))))
    
    ax.set_xlabel('Round')
    ax.set_ylabel('RMSE(log₁₀)')
    ax.set_title('Root Mean Square Error Comparison\n(Mean ± Std across seeds)')
    ax.set_xticks(x)
    ax.set_xticklabels([f"R{row['round']:02d}\n{row['description'][:15]}" 
                       for _, row in rmse_stats.iterrows()], rotation=45)
    ax.grid(True, alpha=0.3)
    
    # Add value labels on bars
    for i, bar in enumerate(bars):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + rmse_stats.iloc[i]['std'],
                f'{height:.3f}', ha='center', va='bottom', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def plot_population_trajectories(results, output_path):
    """Plot population trajectories for the best round (Round 02)."""
    # Use Round 02 (K_half=200,000) as it appears to be the best performer
    best_round = 2
    
    if best_round not in results:
        print(f"Best round {best_round} not available for trajectory plot")
        return
    
    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    axes = axes.flatten()
    
    for i, region in enumerate(REGION_ORDER):
        ax = axes[i]
        
        # Collect trajectories for all seeds
        trajectories = []
        
        for seed, data in results[best_round].items():
            if region in data.get('region_details', {}):
                yearly_totals = data['region_details'][region].get('yearly_totals', [])
                if yearly_totals:
                    trajectories.append(yearly_totals)
        
        if trajectories:
            # Convert to numpy array for easier manipulation
            max_years = max(len(traj) for traj in trajectories)
            padded_trajectories = []
            
            for traj in trajectories:
                padded = traj + [np.nan] * (max_years - len(traj))
                padded_trajectories.append(padded)
            
            trajectories_array = np.array(padded_trajectories)
            years = np.arange(max_years)
            
            # Plot individual seed trajectories as thin lines
            colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
            seeds = list(results[best_round].keys())
            
            for j, seed in enumerate(seeds):
                if j < len(trajectories):
                    ax.plot(years, trajectories[j], color=colors[j], alpha=0.6, linewidth=1,
                           label=f'Seed {seed}')
            
            # Plot mean trajectory as thick line
            mean_traj = np.nanmean(trajectories_array, axis=0)
            ax.plot(years, mean_traj, 'black', linewidth=3, label='Mean')
        
        ax.set_title(f'{region}')
        ax.set_xlabel('Year')
        ax.set_ylabel('Population')
        ax.grid(True, alpha=0.3)
        
        if i == 0:  # Only add legend to first subplot
            ax.legend()
    
    fig.suptitle(f'Population Trajectories - Round {best_round:02d} (K_half=200,000)', fontsize=16)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def plot_scatter_plot(df, output_path):
    """Create scatter plot of actual vs target recovery."""
    fig, ax = plt.subplots(figsize=(10, 8))
    
    colors = plt.cm.Set3(np.linspace(0, 1, len(df['round'].unique())))
    
    for i, round_num in enumerate(sorted(df['round'].unique())):
        round_data = df[df['round'] == round_num]
        
        targets = []
        actuals = []
        
        for region in REGION_ORDER:
            target_col = f'{region}_target'
            actual_col = f'{region}_actual_pct'
            
            if target_col in round_data.columns:
                for _, row in round_data.iterrows():
                    if not np.isnan(row[actual_col]):
                        targets.append(row[target_col] * 100)  # Convert to percentage
                        actuals.append(row[actual_col])
        
        if targets and actuals:
            ax.scatter(targets, actuals, color=colors[i], alpha=0.7, s=50,
                      label=f"R{round_num:02d}: {ROUND_DESCRIPTIONS.get(round_num, '')[:15]}")
    
    # Add reference lines
    lims = [0.01, 100]
    ax.plot(lims, lims, 'k-', alpha=0.75, linewidth=2, label='1:1 line')
    ax.plot(lims, [2*x for x in lims], 'r--', alpha=0.5, label='2× bounds')
    ax.plot(lims, [0.5*x for x in lims], 'r--', alpha=0.5)
    ax.plot(lims, [5*x for x in lims], 'b--', alpha=0.5, label='5× bounds')
    ax.plot(lims, [0.2*x for x in lims], 'b--', alpha=0.5)
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel('Target Recovery (%)')
    ax.set_ylabel('Actual Recovery (%)')
    ax.set_title('Model Performance: Actual vs Target Recovery\n(All regions, all seeds)')
    ax.grid(True, alpha=0.3)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def create_summary_table(df):
    """Create summary statistics table."""
    summary_rows = []
    
    for round_num in sorted(df['round'].unique()):
        round_data = df[df['round'] == round_num]
        
        row = {
            'Round': f"R{round_num:02d}",
            'Description': ROUND_DESCRIPTIONS.get(round_num, f"Round {round_num}"),
            'Seeds': len(round_data),
            'RMSE(log)': f"{round_data['rmse_log'].mean():.3f} ± {round_data['rmse_log'].std():.3f}",
            'Within 2×': f"{round_data['within_2x'].mean():.2f}",
            'Within 5×': f"{round_data['within_5x'].mean():.2f}",
            'Runtime (h)': f"{round_data['wall_time_hours'].mean():.1f} ± {round_data['wall_time_hours'].std():.1f}"
        }
        
        summary_rows.append(row)
    
    return pd.DataFrame(summary_rows)

def main():
    """Main analysis function."""
    print("Loading calibration results...")
    results = load_results()
    
    if not results:
        print("No results found!")
        return
    
    print(f"Found data for {len(results)} rounds")
    
    # Create summary DataFrame
    df = create_summary_df(results)
    print(f"Created summary with {len(df)} rows")
    
    # Generate figures
    print("Creating Figure 1: Grouped bar chart...")
    plot_grouped_bar_chart(df, OUTPUT_DIR / "figure1_grouped_bars.pdf")
    
    print("Creating Figure 2: Performance heatmap...")
    plot_heatmap(df, OUTPUT_DIR / "figure2_heatmap.pdf")
    
    print("Creating Figure 3: RMSE comparison...")
    plot_rmse_comparison(df, OUTPUT_DIR / "figure3_rmse.pdf")
    
    print("Creating Figure 4: Population trajectories...")
    plot_population_trajectories(results, OUTPUT_DIR / "figure4_trajectories.pdf")
    
    print("Creating Figure 5: Scatter plot...")
    plot_scatter_plot(df, OUTPUT_DIR / "figure5_scatter.pdf")
    
    # Create summary table
    print("Creating summary table...")
    summary_table = create_summary_table(df)
    summary_table.to_csv(OUTPUT_DIR / "summary_table.csv", index=False)
    print(summary_table.to_string(index=False))
    
    print(f"\nAll figures saved to: {OUTPUT_DIR}")

if __name__ == "__main__":
    main()