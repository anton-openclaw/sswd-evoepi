#!/usr/bin/env python3
"""
Generate population trajectory figures for all calibration rounds.
Creates both multi-panel and individual round figures.
"""

import json
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages

# Configuration
TARGET_REGIONS = ['AK-PWS', 'AK-FN', 'AK-FS', 'BC-N', 'SS-S', 'JDF', 'OR', 'CA-N']
ROUNDS = ['round_00', 'round_01', 'round_02', 'round_03', 'round_04']
ROUND_NAMES = {
    'round_00': 'R00: Defaults',
    'round_01': 'R01: K_half=30K',
    'round_02': 'R02: K_half=200K', 
    'round_03': 'R03: P_env_max=100',
    'round_04': 'R04: P_env_max=2000'
}

SEEDS = [42, 123, 999]

# Use tab10 colormap for consistent colors
colors = plt.cm.tab10(np.linspace(0, 1, len(TARGET_REGIONS)))
REGION_COLORS = dict(zip(TARGET_REGIONS, colors))

def load_data(results_dir):
    """Load all result data from calibration rounds."""
    data = {}
    
    for round_name in ROUNDS:
        data[round_name] = {}
        round_dir = os.path.join(results_dir, 'calibration', round_name)
        
        for seed in SEEDS:
            result_file = os.path.join(round_dir, f'result_seed{seed}.json')
            if not os.path.exists(result_file):
                print(f"Warning: Missing {result_file}")
                continue
                
            with open(result_file, 'r') as f:
                result = json.load(f)
                
            # Handle early-stopped runs with empty per_region
            if not result.get('scoring', {}).get('per_region', {}):
                print(f"Warning: Skipping {round_name} seed {seed} - empty per_region")
                continue
                
            data[round_name][seed] = result
    
    return data

def calculate_recovery_fractions(data):
    """Calculate recovery fractions for all regions and seeds."""
    recovery_data = {}
    target_data = {}
    
    for round_name in ROUNDS:
        recovery_data[round_name] = {}
        target_data[round_name] = {}
        
        if round_name not in data:
            continue
            
        for region in TARGET_REGIONS:
            recovery_data[round_name][region] = []
            target_data[round_name][region] = None
            
            for seed in SEEDS:
                if seed not in data[round_name]:
                    continue
                    
                result = data[round_name][seed]
                
                # Get yearly totals for this region
                if region not in result.get('region_details', {}):
                    continue
                    
                yearly_totals = result['region_details'][region]['yearly_totals']
                if len(yearly_totals) < 13:
                    print(f"Warning: {round_name} seed {seed} region {region} has incomplete data")
                    continue
                
                # Calculate recovery fractions (population / initial population)
                initial_pop = yearly_totals[0]
                if initial_pop == 0:
                    continue
                    
                fractions = [pop / initial_pop for pop in yearly_totals]
                recovery_data[round_name][region].append(fractions)
                
                # Get target fraction (only need to do this once per region per round)
                if target_data[round_name][region] is None:
                    if region in result.get('scoring', {}).get('per_region', {}):
                        target_pct = result['scoring']['per_region'][region]['target_pct']
                        target_data[round_name][region] = target_pct / 100.0
                    
    return recovery_data, target_data

def create_multi_panel_figure(recovery_data, target_data, output_file):
    """Create multi-panel figure with all 5 rounds."""
    fig, axes = plt.subplots(1, 5, figsize=(20, 4))
    years = np.arange(13)
    
    for i, round_name in enumerate(ROUNDS):
        ax = axes[i]
        
        if round_name not in recovery_data:
            ax.text(0.5, 0.5, f'No data\nfor {ROUND_NAMES[round_name]}', 
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_title(ROUND_NAMES[round_name])
            continue
        
        # Plot individual seed trajectories (thin lines)
        for region in TARGET_REGIONS:
            if region not in recovery_data[round_name]:
                continue
                
            trajectories = recovery_data[round_name][region]
            if not trajectories:
                continue
                
            color = REGION_COLORS[region]
            
            for trajectory in trajectories:
                ax.plot(years, trajectory, color=color, alpha=0.3, linewidth=0.5)
            
            # Plot mean trajectory (thick line)
            if trajectories:
                mean_trajectory = np.mean(trajectories, axis=0)
                ax.plot(years, mean_trajectory, color=color, linewidth=2, 
                       label=region if i == 0 else "")  # Only label in first panel
            
            # Add target line
            if region in target_data[round_name] and target_data[round_name][region] is not None:
                target_val = target_data[round_name][region]
                ax.axhline(y=target_val, color=color, linestyle='--', alpha=0.7, linewidth=1)
        
        ax.set_yscale('log')
        ax.set_ylim(0.001, 2.0)
        ax.set_xlabel('Year')
        if i == 0:
            ax.set_ylabel('Recovery Fraction')
        ax.set_title(ROUND_NAMES[round_name], fontsize=10)
        ax.grid(True, alpha=0.3)
        
        # Only add legend to first panel
        if i == 0:
            ax.legend(bbox_to_anchor=(0, -0.3), loc='upper left', ncol=2, fontsize=8)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved multi-panel figure: {output_file}")

def create_individual_figures(recovery_data, target_data, output_dir):
    """Create individual figures for each round."""
    years = np.arange(13)
    
    for i, round_name in enumerate(ROUNDS):
        if round_name not in recovery_data:
            continue
            
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Plot trajectories for each region
        legend_elements = []
        for region in TARGET_REGIONS:
            if region not in recovery_data[round_name]:
                continue
                
            trajectories = recovery_data[round_name][region]
            if not trajectories:
                continue
                
            color = REGION_COLORS[region]
            
            # Plot individual seed trajectories (thin lines)
            for trajectory in trajectories:
                ax.plot(years, trajectory, color=color, alpha=0.3, linewidth=0.8)
            
            # Plot mean trajectory (thick line)
            mean_trajectory = np.mean(trajectories, axis=0)
            ax.plot(years, mean_trajectory, color=color, linewidth=3, label=region)
            
            # Add target line
            if region in target_data[round_name] and target_data[round_name][region] is not None:
                target_val = target_data[round_name][region]
                ax.axhline(y=target_val, color=color, linestyle='--', alpha=0.8, linewidth=2)
                
                # Add to legend 
                legend_elements.append(
                    mpatches.Patch(color=color, label=f'{region} (target: {target_val:.1%})')
                )
        
        ax.set_yscale('log')
        ax.set_ylim(0.001, 2.0)
        ax.set_xlabel('Year', fontsize=14)
        ax.set_ylabel('Recovery Fraction', fontsize=14)
        ax.set_title(f'{ROUND_NAMES[round_name]} - Population Recovery Trajectories', fontsize=16)
        ax.grid(True, alpha=0.3)
        
        # Add legend with target values
        if legend_elements:
            ax.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left')
        
        # Save individual figure
        output_file = os.path.join(output_dir, f'figure4{"abcde"[i]}_round{round_name[-2:]}.pdf')
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved individual figure: {output_file}")

def main():
    # Set working directory
    os.chdir('/home/starbot/.openclaw/workspace/sswd-evoepi')
    
    # Load all data
    print("Loading data from calibration rounds...")
    data = load_data('results')
    
    # Calculate recovery fractions
    print("Calculating recovery fractions...")
    recovery_data, target_data = calculate_recovery_fractions(data)
    
    # Create output directory
    output_dir = 'reports/calibration_r1/figures'
    os.makedirs(output_dir, exist_ok=True)
    
    # Create multi-panel figure
    print("Creating multi-panel figure...")
    multi_panel_file = os.path.join(output_dir, 'figure4_trajectories.pdf')
    create_multi_panel_figure(recovery_data, target_data, multi_panel_file)
    
    # Create individual figures
    print("Creating individual round figures...")
    create_individual_figures(recovery_data, target_data, output_dir)
    
    print("All trajectory figures completed!")

if __name__ == "__main__":
    main()