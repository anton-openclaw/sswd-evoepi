#!/usr/bin/env python3
"""
Create timing comparison visualizations.
"""

import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

def load_timing_data():
    """Load timing comparison data."""
    data_file = Path('results/performance/focused_timing_comparison.json')
    if data_file.exists():
        with open(data_file) as f:
            return json.load(f)
    return None

def create_comparison_plot(data):
    """Create before/after comparison plot."""
    results = data['comparison_results']
    valid_results = [r for r in results if not np.isnan(r.get('overhead_factor', np.nan))]
    
    if len(valid_results) == 0:
        print("No valid timing data for plotting")
        return None
    
    # Set up the plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle('SSWD-EvoEpi Spawning System: Before/After Optimization', fontsize=14, fontweight='bold')
    
    # Extract data for plotting
    labels = [r['label'].title() for r in valid_results]
    n_agents = [r['n_agents'] for r in valid_results]
    spawn_times = [r['spawning_time_s'] for r in valid_results]
    nospawn_times = [r['no_spawning_time_s'] for r in valid_results]
    overhead_factors = [r['overhead_factor'] for r in valid_results]
    
    x = np.arange(len(labels))
    width = 0.35
    
    # Plot 1: Absolute timing comparison
    bars1 = ax1.bar(x - width/2, spawn_times, width, label='With Spawning', color='#2E86AB', alpha=0.8)
    bars2 = ax1.bar(x + width/2, nospawn_times, width, label='Without Spawning', color='#A23B72', alpha=0.8)
    
    ax1.set_xlabel('Population Size')
    ax1.set_ylabel('Runtime (seconds)')
    ax1.set_title('Absolute Runtime Comparison')
    ax1.set_xticks(x)
    ax1.set_xticklabels([f"{l}\n({n} agents)" for l, n in zip(labels, n_agents)])
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Add value labels on bars
    for bar in bars1:
        height = bar.get_height()
        ax1.annotate(f'{height:.1f}s',
                    xy=(bar.get_x() + bar.get_width()/2, height),
                    xytext=(0, 3), textcoords="offset points",
                    ha='center', va='bottom', fontsize=9)
    
    for bar in bars2:
        height = bar.get_height()
        ax1.annotate(f'{height:.1f}s',
                    xy=(bar.get_x() + bar.get_width()/2, height),
                    xytext=(0, 3), textcoords="offset points",
                    ha='center', va='bottom', fontsize=9)
    
    # Plot 2: Spawning overhead factor
    bars3 = ax2.bar(x, overhead_factors, color='#F18F01', alpha=0.8)
    ax2.axhline(y=1.0, color='black', linestyle='--', alpha=0.5, label='No overhead')
    
    ax2.set_xlabel('Population Size')
    ax2.set_ylabel('Overhead Factor (×)')
    ax2.set_title('Spawning System Overhead')
    ax2.set_xticks(x)
    ax2.set_xticklabels([f"{l}\n({n} agents)" for l, n in zip(labels, n_agents)])
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    # Add overhead factor labels
    for i, (bar, factor) in enumerate(zip(bars3, overhead_factors)):
        height = bar.get_height()
        ax2.annotate(f'{factor:.2f}×',
                    xy=(bar.get_x() + bar.get_width()/2, height),
                    xytext=(0, 3), textcoords="offset points",
                    ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    plt.tight_layout()
    return fig

def create_scaling_analysis(data):
    """Create scaling analysis plot."""
    results = data['comparison_results']
    valid_results = [r for r in results if not np.isnan(r.get('overhead_factor', np.nan))]
    
    if len(valid_results) < 2:
        return None
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle('Spawning System Performance Scaling Analysis', fontsize=14, fontweight='bold')
    
    n_agents = np.array([r['n_agents'] for r in valid_results])
    spawn_times = np.array([r['spawning_time_s'] for r in valid_results])
    overhead_s = np.array([r['overhead_s'] for r in valid_results])
    
    # Plot 1: Runtime scaling
    ax1.plot(n_agents, spawn_times, 'o-', color='#2E86AB', linewidth=2, markersize=8, label='With Spawning')
    ax1.set_xlabel('Population Size (agents)')
    ax1.set_ylabel('Runtime (seconds)')
    ax1.set_title('Runtime Scaling')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # Add data point labels
    for n, t in zip(n_agents, spawn_times):
        ax1.annotate(f'{t:.1f}s', (n, t), xytext=(5, 5), textcoords='offset points')
    
    # Plot 2: Overhead per agent
    overhead_per_agent = overhead_s / n_agents * 1000  # ms per agent
    ax2.plot(n_agents, overhead_per_agent, 'o-', color='#F18F01', linewidth=2, markersize=8)
    ax2.set_xlabel('Population Size (agents)')
    ax2.set_ylabel('Spawning Overhead (ms per agent)')
    ax2.set_title('Per-Agent Overhead Scaling')
    ax2.grid(True, alpha=0.3)
    
    # Add trend analysis
    if len(n_agents) > 2:
        z = np.polyfit(n_agents, overhead_per_agent, 1)
        p = np.poly1d(z)
        ax2.plot(n_agents, p(n_agents), "--", alpha=0.8, color='red', 
                label=f'Trend: {z[0]:.3f} ms/agent² + {z[1]:.1f}')
        ax2.legend()
    
    # Add data point labels
    for n, o in zip(n_agents, overhead_per_agent):
        ax2.annotate(f'{o:.1f}', (n, o), xytext=(5, 5), textcoords='offset points')
    
    plt.tight_layout()
    return fig

def main():
    print("Creating timing comparison visualizations...")
    
    # Load data
    data = load_timing_data()
    if data is None:
        print("No timing data found. Skipping visualization.")
        return
    
    print(f"Loaded data with {len(data['comparison_results'])} comparisons")
    
    # Create plots
    fig1 = create_comparison_plot(data)
    if fig1:
        fig1.savefig('results/performance/timing_comparison.png', dpi=300, bbox_inches='tight', 
                     facecolor='white', edgecolor='none')
        print("Saved: results/performance/timing_comparison.png")
    
    fig2 = create_scaling_analysis(data)
    if fig2:
        fig2.savefig('results/performance/spawning_scaling_analysis.png', dpi=300, bbox_inches='tight',
                     facecolor='white', edgecolor='none')
        print("Saved: results/performance/spawning_scaling_analysis.png")
    
    # Print summary
    if data['summary']['successful_comparisons'] > 0:
        print(f"\nSummary:")
        print(f"  Successful comparisons: {data['summary']['successful_comparisons']}")
        print(f"  Mean overhead: {data['summary']['mean_overhead_s']:.3f}s")
        print(f"  Mean overhead factor: {data['summary']['mean_overhead_factor']:.2f}×")
    
    print("Visualization complete!")

if __name__ == '__main__':
    main()