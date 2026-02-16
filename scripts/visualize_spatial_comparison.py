#!/usr/bin/env python3
"""Generate spatial validation visualizations comparing distance methods.

Creates publication-quality plots comparing Haversine vs overwater distance 
methods for the SSWD-EvoEpi model.

Generates:
  a) distance_comparison_heatmap.png - Side-by-side distance matrices
  b) connectivity_matrices.png - C and D matrix comparison
  c) epidemic_spread_comparison.png - Population trajectory comparison
  d) epidemic_wavefront_map.png - Geographic spread visualization
  e) dispersal_kernel_comparison.png - Distance kernel analysis

Dark theme with scientific color palette.
"""

import json
import sys
from pathlib import Path
from typing import Dict, Any, Tuple, List, Optional

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

# ═══════════════════════════════════════════════════════════════════════
# CONFIGURATION & STYLING
# ═══════════════════════════════════════════════════════════════════════

COMPARISON_DIR = project_root / "results" / "distance_method_comparison"
OUTPUT_DIR = project_root / "results" / "spatial_validation"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Dark theme colors
BG_COLOR = '#1a1a2e'
TEXT_COLOR = '#eeeeff'
ACCENT_COLOR = '#00ff66'  # Circuit board green
NEON_BLUE = '#16537e'
NEON_ORANGE = '#ff6b35'
NEON_PURPLE = '#7209b7'

# Set dark theme globally
plt.style.use('dark_background')
plt.rcParams.update({
    'figure.facecolor': BG_COLOR,
    'axes.facecolor': BG_COLOR,
    'text.color': TEXT_COLOR,
    'axes.labelcolor': TEXT_COLOR,
    'xtick.color': TEXT_COLOR,
    'ytick.color': TEXT_COLOR,
    'axes.edgecolor': TEXT_COLOR,
    'grid.color': '#333344',
    'font.size': 10,
    'axes.titlesize': 12,
    'figure.titlesize': 14,
})

# Create custom colormap
def create_neon_cmap(base_color: str, name: str) -> LinearSegmentedColormap:
    """Create a neon-style colormap."""
    colors = ['#000011', base_color, '#ffffff']
    return LinearSegmentedColormap.from_list(name, colors, N=256)

DISTANCE_CMAP = create_neon_cmap('#16537e', 'neon_blue')
CONNECTIVITY_CMAP = create_neon_cmap('#00ff66', 'neon_green')
EPIDEMIC_CMAP = create_neon_cmap('#ff6b35', 'neon_orange')

# Node definitions for 5-node network
NODE_NAMES = ['Sitka', 'Howe Sound', 'San Juan Islands', 'Newport', 'Monterey']
NODE_COORDS = {
    'Sitka': (-135.33, 57.05),
    'Howe Sound': (-123.16, 49.57), 
    'San Juan Islands': (-123.02, 48.52),
    'Newport': (-124.05, 44.63),
    'Monterey': (-121.89, 36.60)
}

# ═══════════════════════════════════════════════════════════════════════
# DATA LOADING
# ═══════════════════════════════════════════════════════════════════════

def load_scenario_data(scenario_name: str) -> Tuple[Dict, np.ndarray, Optional[np.ndarray], Optional[np.ndarray]]:
    """Load data for a single scenario."""
    scenario_dir = COMPARISON_DIR / scenario_name
    
    # Load metadata
    with open(scenario_dir / "metadata.json", 'r') as f:
        metadata = json.load(f)
    
    # Load arrays
    data = np.load(scenario_dir / "results.npz")
    yearly_pop = data['yearly_pop']
    C_matrix = data.get('C_matrix')
    D_matrix = data.get('D_matrix')
    distances = data.get('distances')
    
    return metadata, yearly_pop, C_matrix, D_matrix, distances

def load_all_scenarios() -> Dict[str, Tuple]:
    """Load data for all available scenarios."""
    scenarios = {}
    
    for scenario_dir in COMPARISON_DIR.iterdir():
        if scenario_dir.is_dir() and (scenario_dir / "results.npz").exists():
            scenario_name = scenario_dir.name
            scenarios[scenario_name] = load_scenario_data(scenario_name)
    
    return scenarios

# ═══════════════════════════════════════════════════════════════════════
# PLOT A: DISTANCE COMPARISON HEATMAP
# ═══════════════════════════════════════════════════════════════════════

def plot_distance_comparison(scenarios: Dict[str, Tuple], save_path: Path):
    """Create side-by-side distance matrix heatmaps."""
    fig = plt.figure(figsize=(14, 6))
    gs = GridSpec(1, 2, figure=fig, wspace=0.3)
    
    # Extract distance data (from connectivity - inverse relationship)
    haversine_distances = None
    overwater_distances = None
    
    for name, (metadata, yearly_pop, C_matrix, D_matrix, distances) in scenarios.items():
        if 'Haversine' in name and haversine_distances is None:
            # Convert connectivity to approximate distances (inverse log relationship)
            if C_matrix is not None:
                # Higher connectivity = shorter distance
                haversine_distances = -np.log(C_matrix + 1e-10) * 100  # Scale for visibility
        elif 'Overwater' in name and overwater_distances is None:
            if distances is not None:
                overwater_distances = distances
            elif C_matrix is not None:
                overwater_distances = -np.log(C_matrix + 1e-10) * 100
    
    # Plot Haversine
    ax1 = fig.add_subplot(gs[0, 0])
    if haversine_distances is not None:
        im1 = ax1.imshow(haversine_distances, cmap=DISTANCE_CMAP, aspect='equal')
        ax1.set_title('Haversine × 1.5 Distances', color=TEXT_COLOR, fontweight='bold')
        
        # Annotate with values
        for i in range(len(NODE_NAMES)):
            for j in range(len(NODE_NAMES)):
                if i != j:
                    dist_val = haversine_distances[i, j]
                    ax1.text(j, i, f'{dist_val:.0f}', ha='center', va='center',
                            color='white' if dist_val > haversine_distances.max()/2 else 'black',
                            fontsize=8, fontweight='bold')
        
        plt.colorbar(im1, ax=ax1, label='Distance (scaled)', shrink=0.8)
    else:
        ax1.text(0.5, 0.5, 'Haversine data\nnot available', ha='center', va='center',
                transform=ax1.transAxes, color=TEXT_COLOR, fontsize=12)
    
    ax1.set_xticks(range(len(NODE_NAMES)))
    ax1.set_yticks(range(len(NODE_NAMES)))
    ax1.set_xticklabels([name.replace(' ', '\n') for name in NODE_NAMES], rotation=45, ha='right')
    ax1.set_yticklabels(NODE_NAMES)
    
    # Plot Overwater
    ax2 = fig.add_subplot(gs[0, 1])
    if overwater_distances is not None:
        im2 = ax2.imshow(overwater_distances, cmap=DISTANCE_CMAP, aspect='equal')
        ax2.set_title('Overwater Distances', color=TEXT_COLOR, fontweight='bold')
        
        # Annotate with values
        for i in range(len(NODE_NAMES)):
            for j in range(len(NODE_NAMES)):
                if i != j:
                    dist_val = overwater_distances[i, j]
                    ax2.text(j, i, f'{dist_val:.0f}', ha='center', va='center',
                            color='white' if dist_val > overwater_distances.max()/2 else 'black',
                            fontsize=8, fontweight='bold')
        
        plt.colorbar(im2, ax=ax2, label='Distance (km)', shrink=0.8)
    else:
        ax2.text(0.5, 0.5, 'Overwater data\nnot available', ha='center', va='center',
                transform=ax2.transAxes, color=TEXT_COLOR, fontsize=12)
    
    ax2.set_xticks(range(len(NODE_NAMES)))
    ax2.set_yticks(range(len(NODE_NAMES)))
    ax2.set_xticklabels([name.replace(' ', '\n') for name in NODE_NAMES], rotation=45, ha='right')
    ax2.set_yticklabels(NODE_NAMES)
    
    plt.suptitle('Distance Matrix Comparison', color=TEXT_COLOR, fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight', facecolor=BG_COLOR)
    plt.close()
    
    print(f"Saved: {save_path}")

# ═══════════════════════════════════════════════════════════════════════
# PLOT B: CONNECTIVITY MATRICES
# ═══════════════════════════════════════════════════════════════════════

def plot_connectivity_matrices(scenarios: Dict[str, Tuple], save_path: Path):
    """Create side-by-side C and D matrix comparison."""
    fig = plt.figure(figsize=(16, 8))
    gs = GridSpec(2, 2, figure=fig, hspace=0.3, wspace=0.3)
    
    # Extract matrices
    haversine_C, haversine_D = None, None
    overwater_C, overwater_D = None, None
    
    for name, (metadata, yearly_pop, C_matrix, D_matrix, distances) in scenarios.items():
        if 'Haversine_1.5x' in name:
            haversine_C, haversine_D = C_matrix, D_matrix
        elif 'Overwater' in name:
            overwater_C, overwater_D = C_matrix, D_matrix
    
    matrices = [
        (haversine_C, "Haversine × 1.5\nLarval Connectivity (C)", (0, 0)),
        (overwater_C, "Overwater\nLarval Connectivity (C)", (0, 1)),
        (haversine_D, "Haversine × 1.5\nPathogen Dispersal (D)", (1, 0)),
        (overwater_D, "Overwater\nPathogen Dispersal (D)", (1, 1)),
    ]
    
    for matrix, title, pos in matrices:
        ax = fig.add_subplot(gs[pos[0], pos[1]])
        
        if matrix is not None:
            # Use log scale for better visualization
            matrix_log = np.log10(matrix + 1e-12)
            im = ax.imshow(matrix_log, cmap=CONNECTIVITY_CMAP, aspect='equal')
            
            ax.set_title(title, color=TEXT_COLOR, fontweight='bold')
            
            # Annotate with original values (scientific notation)
            for i in range(len(NODE_NAMES)):
                for j in range(len(NODE_NAMES)):
                    if i != j and matrix[i, j] > 1e-8:  # Only show significant values
                        val = matrix[i, j]
                        ax.text(j, i, f'{val:.2e}', ha='center', va='center',
                               color='white', fontsize=7, fontweight='bold')
            
            plt.colorbar(im, ax=ax, label='log₁₀(connectivity)', shrink=0.8)
        else:
            ax.text(0.5, 0.5, 'Matrix data\nnot available', ha='center', va='center',
                   transform=ax.transAxes, color=TEXT_COLOR, fontsize=12)
            ax.set_title(title, color=TEXT_COLOR, fontweight='bold')
        
        ax.set_xticks(range(len(NODE_NAMES)))
        ax.set_yticks(range(len(NODE_NAMES)))
        ax.set_xticklabels([name.replace(' ', '\n') for name in NODE_NAMES], rotation=45, ha='right')
        ax.set_yticklabels(NODE_NAMES)
    
    plt.suptitle('Connectivity Matrix Comparison', color=TEXT_COLOR, fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight', facecolor=BG_COLOR)
    plt.close()
    
    print(f"Saved: {save_path}")

# ═══════════════════════════════════════════════════════════════════════
# PLOT C: EPIDEMIC SPREAD COMPARISON
# ═══════════════════════════════════════════════════════════════════════

def plot_epidemic_spread_comparison(scenarios: Dict[str, Tuple], save_path: Path):
    """Create 5-panel plot showing population trajectories."""
    fig, axes = plt.subplots(1, 5, figsize=(20, 4), sharey=True)
    
    # Colors for different scenarios
    colors = {'Haversine_1.5x': NEON_BLUE, 'Overwater': ACCENT_COLOR, 'Haversine_1.0x': NEON_PURPLE}
    styles = {'Haversine_1.5x': '--', 'Overwater': '-', 'Haversine_1.0x': '-.'}
    
    disease_year = None
    years = None
    
    for node_idx in range(5):
        ax = axes[node_idx]
        node_name = NODE_NAMES[node_idx]
        
        for scenario_name, (metadata, yearly_pop, C_matrix, D_matrix, distances) in scenarios.items():
            if disease_year is None:
                disease_year = metadata['disease_year']
                years = np.arange(metadata['n_years'])
            
            # Plot population trajectory
            pop_trajectory = yearly_pop[node_idx, :]
            color = colors.get(scenario_name, TEXT_COLOR)
            style = styles.get(scenario_name, '-')
            
            ax.plot(years, pop_trajectory, style, color=color, 
                   label=scenario_name.replace('_', ' '), linewidth=2, alpha=0.8)
        
        # Shade epidemic period
        ax.axvspan(disease_year, years[-1], alpha=0.2, color=NEON_ORANGE, 
                  label='Epidemic Period' if node_idx == 0 else "")
        
        ax.set_title(f'{node_name}', color=TEXT_COLOR, fontweight='bold')
        ax.set_xlabel('Year')
        ax.grid(True, alpha=0.3)
        
        if node_idx == 0:
            ax.set_ylabel('Population')
            ax.legend(loc='upper right', fontsize=8)
    
    plt.suptitle('Population Trajectories: Distance Method Comparison', 
                color=TEXT_COLOR, fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight', facecolor=BG_COLOR)
    plt.close()
    
    print(f"Saved: {save_path}")

# ═══════════════════════════════════════════════════════════════════════
# PLOT D: EPIDEMIC WAVEFRONT MAP
# ═══════════════════════════════════════════════════════════════════════

def analyze_epidemic_arrival(yearly_pop: np.ndarray, disease_year: int) -> Dict[int, int]:
    """Determine epidemic arrival time at each node."""
    arrival_times = {}
    
    for node_id in range(yearly_pop.shape[0]):
        # Look for significant population decline after disease year
        pre_pop = yearly_pop[node_id, disease_year]
        post_pops = yearly_pop[node_id, disease_year:]
        
        # Find first year with >10% decline
        decline_years = np.where((pre_pop - post_pops) / pre_pop > 0.1)[0]
        
        if len(decline_years) > 0:
            arrival_times[node_id] = disease_year + decline_years[0]
        else:
            arrival_times[node_id] = -1  # Never arrived
    
    return arrival_times

def plot_epidemic_wavefront_map(scenarios: Dict[str, Tuple], save_path: Path):
    """Create geographic map showing epidemic spread patterns."""
    fig = plt.figure(figsize=(16, 8))
    gs = GridSpec(1, 2, figure=fig, wspace=0.3)
    
    scenario_pairs = [
        ('Haversine_1.5x', 'Haversine × 1.5'),
        ('Overwater', 'Overwater Distances')
    ]
    
    for idx, (scenario_name, display_name) in enumerate(scenario_pairs):
        if scenario_name not in scenarios:
            continue
            
        ax = fig.add_subplot(gs[0, idx])
        
        metadata, yearly_pop, C_matrix, D_matrix, distances = scenarios[scenario_name]
        disease_year = metadata['disease_year']
        
        # Get epidemic arrival times
        arrival_times = analyze_epidemic_arrival(yearly_pop, disease_year)
        
        # Plot coastline (simplified)
        coast_lons = [-140, -135, -130, -125, -122, -120]
        coast_lats = [60, 58, 54, 48, 40, 35]
        ax.plot(coast_lons, coast_lats, '-', color=ACCENT_COLOR, linewidth=3, alpha=0.6)
        
        # Plot nodes with color coding by arrival time
        max_arrival = max([t for t in arrival_times.values() if t > 0])
        min_arrival = min([t for t in arrival_times.values() if t > 0])
        
        for node_idx, node_name in enumerate(NODE_NAMES):
            lon, lat = NODE_COORDS[node_name]
            arrival_time = arrival_times[node_idx]
            
            if arrival_time > 0:
                # Color by arrival time
                normalized_time = (arrival_time - min_arrival) / (max_arrival - min_arrival)
                color = plt.cm.plasma(normalized_time)
                
                ax.scatter(lon, lat, s=200, c=[color], edgecolor='white', 
                          linewidth=2, zorder=10)
                ax.annotate(f'{node_name}\nYear {arrival_time}', 
                           (lon, lat), xytext=(5, 5), textcoords='offset points',
                           color=TEXT_COLOR, fontsize=9, fontweight='bold',
                           bbox=dict(boxstyle='round,pad=0.3', facecolor=BG_COLOR, 
                                   alpha=0.8, edgecolor='white'))
            else:
                # No epidemic arrival
                ax.scatter(lon, lat, s=200, c='gray', marker='x', linewidth=3, zorder=10)
                ax.annotate(f'{node_name}\nNo arrival', 
                           (lon, lat), xytext=(5, 5), textcoords='offset points',
                           color=TEXT_COLOR, fontsize=9, fontweight='bold',
                           bbox=dict(boxstyle='round,pad=0.3', facecolor=BG_COLOR, 
                                   alpha=0.8, edgecolor='gray'))
        
        # Draw arrows showing spread direction (from earliest to latest)
        sorted_arrivals = sorted([(t, i) for i, t in arrival_times.items() if t > 0])
        for i in range(len(sorted_arrivals) - 1):
            from_node = sorted_arrivals[i][1]
            to_node = sorted_arrivals[i + 1][1]
            
            from_lon, from_lat = NODE_COORDS[NODE_NAMES[from_node]]
            to_lon, to_lat = NODE_COORDS[NODE_NAMES[to_node]]
            
            ax.annotate('', xy=(to_lon, to_lat), xytext=(from_lon, from_lat),
                       arrowprops=dict(arrowstyle='->', color=NEON_ORANGE, 
                                     lw=2, alpha=0.7))
        
        ax.set_title(f'Epidemic Spread\n{display_name}', color=TEXT_COLOR, fontweight='bold')
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
        ax.grid(True, alpha=0.3)
        ax.set_xlim(-142, -118)
        ax.set_ylim(34, 62)
    
    plt.suptitle('Epidemic Wavefront Comparison', color=TEXT_COLOR, fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight', facecolor=BG_COLOR)
    plt.close()
    
    print(f"Saved: {save_path}")

# ═══════════════════════════════════════════════════════════════════════
# PLOT E: DISPERSAL KERNEL COMPARISON
# ═══════════════════════════════════════════════════════════════════════

def plot_dispersal_kernel_comparison(scenarios: Dict[str, Tuple], save_path: Path):
    """Plot exponential dispersal kernel with node pairs marked."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Dispersal parameters (from model config)
    lambda_larvae = 200  # km (larval dispersal scale)
    lambda_vibrio = 50   # km (pathogen dispersal scale)
    
    # Distance range for kernel
    distances = np.linspace(0, 3000, 1000)
    
    # Exponential kernels
    larval_kernel = np.exp(-distances / lambda_larvae)
    vibrio_kernel = np.exp(-distances / lambda_vibrio)
    
    # Plot larval kernel
    ax1.plot(distances, larval_kernel, '-', color=ACCENT_COLOR, linewidth=3, 
            label=f'Larval (λ={lambda_larvae}km)')
    ax1.set_title('Larval Dispersal Kernel', color=TEXT_COLOR, fontweight='bold')
    ax1.set_xlabel('Distance (km)')
    ax1.set_ylabel('Connectivity')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # Plot vibrio kernel  
    ax2.plot(distances, vibrio_kernel, '-', color=NEON_ORANGE, linewidth=3,
            label=f'Pathogen (λ={lambda_vibrio}km)')
    ax2.set_title('Pathogen Dispersal Kernel', color=TEXT_COLOR, fontweight='bold')
    ax2.set_xlabel('Distance (km)')
    ax2.set_ylabel('Connectivity')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    # Mark node pairs on both kernels
    colors = {'Haversine_1.5x': NEON_BLUE, 'Overwater': ACCENT_COLOR}
    markers = {'Haversine_1.5x': 'o', 'Overwater': 's'}
    
    for scenario_name, (metadata, yearly_pop, C_matrix, D_matrix, distances_matrix) in scenarios.items():
        if scenario_name not in ['Haversine_1.5x', 'Overwater']:
            continue
            
        color = colors[scenario_name]
        marker = markers[scenario_name]
        
        if distances_matrix is not None:
            # Extract unique distances (exclude diagonal)
            unique_distances = []
            for i in range(5):
                for j in range(5):
                    if i != j:
                        unique_distances.append(distances_matrix[i, j])
            
            unique_distances = sorted(set(unique_distances))
            
            # Mark on larval kernel
            for dist in unique_distances:
                kernel_val = np.exp(-dist / lambda_larvae)
                ax1.scatter([dist], [kernel_val], c=color, marker=marker, s=60, 
                           edgecolor='white', linewidth=1, alpha=0.8, 
                           label=scenario_name.replace('_', ' ') if dist == unique_distances[0] else "")
            
            # Mark on vibrio kernel
            for dist in unique_distances:
                kernel_val = np.exp(-dist / lambda_vibrio)
                ax2.scatter([dist], [kernel_val], c=color, marker=marker, s=60, 
                           edgecolor='white', linewidth=1, alpha=0.8,
                           label=scenario_name.replace('_', ' ') if dist == unique_distances[0] else "")
    
    ax1.legend()
    ax2.legend()
    
    plt.suptitle('Dispersal Kernel Analysis', color=TEXT_COLOR, fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight', facecolor=BG_COLOR)
    plt.close()
    
    print(f"Saved: {save_path}")

# ═══════════════════════════════════════════════════════════════════════
# MAIN VISUALIZATION GENERATOR
# ═══════════════════════════════════════════════════════════════════════

def main():
    print("=" * 72)
    print("SPATIAL VALIDATION VISUALIZATION GENERATOR")
    print("=" * 72)
    print()
    
    # Check if comparison data exists
    if not COMPARISON_DIR.exists():
        print(f"ERROR: Comparison directory not found: {COMPARISON_DIR}")
        print("Run scripts/compare_distance_methods.py first.")
        return 1
    
    # Load all scenario data
    print("Loading scenario data...")
    scenarios = load_all_scenarios()
    
    if not scenarios:
        print(f"ERROR: No scenario data found in {COMPARISON_DIR}")
        return 1
    
    print(f"Found {len(scenarios)} scenarios:")
    for name in scenarios.keys():
        print(f"  - {name}")
    print()
    
    # Generate visualizations
    print("Generating visualizations...")
    
    # A) Distance comparison heatmap
    plot_distance_comparison(scenarios, OUTPUT_DIR / "distance_comparison_heatmap.png")
    
    # B) Connectivity matrices
    plot_connectivity_matrices(scenarios, OUTPUT_DIR / "connectivity_matrices.png")
    
    # C) Epidemic spread comparison  
    plot_epidemic_spread_comparison(scenarios, OUTPUT_DIR / "epidemic_spread_comparison.png")
    
    # D) Epidemic wavefront map
    plot_epidemic_wavefront_map(scenarios, OUTPUT_DIR / "epidemic_wavefront_map.png")
    
    # E) Dispersal kernel comparison
    plot_dispersal_kernel_comparison(scenarios, OUTPUT_DIR / "dispersal_kernel_comparison.png")
    
    print()
    print("All visualizations saved to:", OUTPUT_DIR)
    return 0

if __name__ == "__main__":
    sys.exit(main())