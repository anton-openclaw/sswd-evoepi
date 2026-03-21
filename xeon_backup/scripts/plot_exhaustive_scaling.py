#!/usr/bin/env python3
"""
Comprehensive scaling analysis plots from all collected axis data.
Creates two figures:
1. scaling_axes_comprehensive.png - 8-panel analysis
2. scaling_projections.png - projections for key scenarios
"""

import json
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.optimize import curve_fit
import os

# Dark theme colors
BG_COLOR = '#1a1a2e'
PANEL_COLOR = '#16213e'
PRIMARY = '#e94560'
SECONDARY = '#00d2ff'
TERTIARY = '#533483'
ACCENT = '#f5a623'
WHITE = '#ffffff'
GRAY = '#8892b0'

def power_law(x, a, b):
    """Power law: y = a * x^b"""
    return a * (x ** b)

def load_all_data():
    """Load all axis data files"""
    data = {}
    for f in glob.glob('results/performance/axis_*.json'):
        name = f.split('axis_')[1].replace('.json', '')
        with open(f) as fh:
            data[name] = json.load(fh)
    return data

def fit_power_law(x_data, y_data, name=""):
    """Fit power law and return parameters + equation string"""
    try:
        # Filter out zeros/negatives for log space
        mask = (np.array(x_data) > 0) & (np.array(y_data) > 0)
        if np.sum(mask) < 3:
            return None, f"No fit ({name})"
        
        x_clean = np.array(x_data)[mask]
        y_clean = np.array(y_data)[mask]
        
        popt, _ = curve_fit(power_law, x_clean, y_clean, p0=[1, 1])
        a, b = popt
        
        # R² calculation
        y_pred = power_law(x_clean, a, b)
        ss_res = np.sum((y_clean - y_pred) ** 2)
        ss_tot = np.sum((y_clean - np.mean(y_clean)) ** 2)
        r_squared = 1 - (ss_res / ss_tot)
        
        eq_str = f"T = {a:.3f}·N^{b:.2f}\nR² = {r_squared:.3f}"
        return (a, b, r_squared), eq_str
    except:
        return None, f"Fit failed ({name})"

def create_comprehensive_figure(data):
    """Create 8-panel comprehensive scaling analysis"""
    fig, axes = plt.subplots(4, 2, figsize=(16, 20), facecolor=BG_COLOR)
    fig.suptitle('SSWD-EvoEpi: Comprehensive Scaling Analysis', 
                 fontsize=20, color=WHITE, y=0.98)
    
    # Panel A: N scaling (single_N)
    ax = axes[0, 0]
    ax.set_facecolor(PANEL_COLOR)
    
    if 'single_N' in data:
        n_vals = [d['n'] for d in data['single_N']]
        t_vals = [d['wall_s'] for d in data['single_N']]
        
        ax.loglog(n_vals, t_vals, 'o-', color=PRIMARY, markersize=8, linewidth=2, label='Measured')
        
        # Fit power law
        fit_params, eq_str = fit_power_law(n_vals, t_vals, "N scaling")
        if fit_params:
            x_fit = np.logspace(np.log10(min(n_vals)), np.log10(max(n_vals)), 100)
            y_fit = power_law(x_fit, *fit_params[:2])
            ax.loglog(x_fit, y_fit, '--', color=SECONDARY, linewidth=2, alpha=0.8)
            ax.text(0.05, 0.95, eq_str, transform=ax.transAxes, color=WHITE, 
                   fontsize=10, verticalalignment='top', fontfamily='monospace')
    else:
        ax.text(0.5, 0.5, 'No Data\n(single_N)', ha='center', va='center', 
               color=GRAY, fontsize=16, transform=ax.transAxes)
    
    ax.set_xlabel('Population Size (N)', color=WHITE)
    ax.set_ylabel('Wall Time (s)', color=WHITE)
    ax.set_title('A. Population Size Scaling', color=WHITE, fontweight='bold')
    ax.tick_params(colors=WHITE)
    ax.grid(True, alpha=0.3)
    
    # Panel B: T scaling (single_T)
    ax = axes[0, 1]
    ax.set_facecolor(PANEL_COLOR)
    
    if 'single_T' in data:
        t_vals = [d['years'] for d in data['single_T']]
        wall_vals = [d['wall_s'] for d in data['single_T']]
        
        ax.loglog(t_vals, wall_vals, 's-', color=PRIMARY, markersize=8, linewidth=2, label='Measured')
        
        # Fit power law
        fit_params, eq_str = fit_power_law(t_vals, wall_vals, "T scaling")
        if fit_params:
            x_fit = np.logspace(np.log10(min(t_vals)), np.log10(max(t_vals)), 100)
            y_fit = power_law(x_fit, *fit_params[:2])
            ax.loglog(x_fit, y_fit, '--', color=SECONDARY, linewidth=2, alpha=0.8)
            ax.text(0.05, 0.95, eq_str, transform=ax.transAxes, color=WHITE, 
                   fontsize=10, verticalalignment='top', fontfamily='monospace')
    else:
        ax.text(0.5, 0.5, 'No Data\n(single_T)', ha='center', va='center', 
               color=GRAY, fontsize=16, transform=ax.transAxes)
    
    ax.set_xlabel('Simulation Years (T)', color=WHITE)
    ax.set_ylabel('Wall Time (s)', color=WHITE)
    ax.set_title('B. Temporal Scaling', color=WHITE, fontweight='bold')
    ax.tick_params(colors=WHITE)
    ax.grid(True, alpha=0.3)
    
    # Panel C: K scaling (spatial_K)
    ax = axes[1, 0]
    ax.set_facecolor(PANEL_COLOR)
    
    if 'spatial_K' in data:
        k_vals = [d['k'] for d in data['spatial_K']]
        wall_vals = [d['wall_s'] for d in data['spatial_K']]
        
        ax.loglog(k_vals, wall_vals, '^-', color=PRIMARY, markersize=8, linewidth=2, label='Measured')
        
        # Fit power law
        fit_params, eq_str = fit_power_law(k_vals, wall_vals, "K scaling")
        if fit_params:
            x_fit = np.logspace(np.log10(min(k_vals)), np.log10(max(k_vals)), 100)
            y_fit = power_law(x_fit, *fit_params[:2])
            ax.loglog(x_fit, y_fit, '--', color=SECONDARY, linewidth=2, alpha=0.8)
            ax.text(0.05, 0.95, eq_str, transform=ax.transAxes, color=WHITE, 
                   fontsize=10, verticalalignment='top', fontfamily='monospace')
    else:
        ax.text(0.5, 0.5, 'No Data\n(spatial_K)', ha='center', va='center', 
               color=GRAY, fontsize=16, transform=ax.transAxes)
    
    ax.set_xlabel('Spatial Nodes (K)', color=WHITE)
    ax.set_ylabel('Wall Time (s)', color=WHITE)
    ax.set_title('C. Spatial Scaling', color=WHITE, fontweight='bold')
    ax.tick_params(colors=WHITE)
    ax.grid(True, alpha=0.3)
    
    # Panel D: N×T heatmap (NxT_matrix)
    ax = axes[1, 1]
    ax.set_facecolor(PANEL_COLOR)
    
    if 'NxT_matrix' in data:
        # Extract unique N and T values
        n_vals = sorted(list(set([d['n'] for d in data['NxT_matrix']])))
        t_vals = sorted(list(set([d['years'] for d in data['NxT_matrix']])))
        
        # Create matrix
        matrix = np.full((len(t_vals), len(n_vals)), np.nan)
        for d in data['NxT_matrix']:
            i = t_vals.index(d['years'])
            j = n_vals.index(d['n'])
            matrix[i, j] = d['wall_s']
        
        # Create heatmap
        im = ax.imshow(matrix, cmap='plasma', aspect='auto', origin='lower')
        
        # Set ticks
        ax.set_xticks(range(len(n_vals)))
        ax.set_xticklabels([str(n) for n in n_vals])
        ax.set_yticks(range(len(t_vals)))
        ax.set_yticklabels([str(t) for t in t_vals])
        
        # Colorbar
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label('Wall Time (s)', color=WHITE)
        cbar.ax.tick_params(colors=WHITE)
        
        ax.set_xlabel('Population Size (N)', color=WHITE)
        ax.set_ylabel('Simulation Years (T)', color=WHITE)
    else:
        ax.text(0.5, 0.5, 'No Data\n(NxT_matrix)', ha='center', va='center', 
               color=GRAY, fontsize=16, transform=ax.transAxes)
    
    ax.set_title('D. N×T Parameter Space', color=WHITE, fontweight='bold')
    ax.tick_params(colors=WHITE)
    
    # Panel E: Spawning overhead (spawning_overhead)
    ax = axes[2, 0]
    ax.set_facecolor(PANEL_COLOR)
    
    if 'spawning_overhead' in data:
        n_vals = [d['n'] for d in data['spawning_overhead']]
        no_spawn = [d['no_spawn_s'] for d in data['spawning_overhead']]
        spawn = [d['spawn_s'] for d in data['spawning_overhead']]
        overhead = [d['overhead'] for d in data['spawning_overhead']]
        
        x = np.arange(len(n_vals))
        width = 0.35
        
        ax.bar(x - width/2, no_spawn, width, label='No Spawning', color=SECONDARY, alpha=0.8)
        ax.bar(x + width/2, spawn, width, label='With Spawning', color=PRIMARY, alpha=0.8)
        
        # Overhead percentages on top
        for i, (oh, sp) in enumerate(zip(overhead, spawn)):
            ax.text(i + width/2, sp + 0.1, f'+{oh:.1f}%', ha='center', color=WHITE, fontsize=8)
        
        ax.set_xticks(x)
        ax.set_xticklabels([str(n) for n in n_vals])
        ax.set_xlabel('Population Size (N)', color=WHITE)
        ax.set_ylabel('Wall Time (s)', color=WHITE)
        ax.legend(facecolor=PANEL_COLOR, edgecolor=WHITE, labelcolor=WHITE)
    else:
        ax.text(0.5, 0.5, 'No Data\n(spawning_overhead)', ha='center', va='center', 
               color=GRAY, fontsize=16, transform=ax.transAxes)
    
    ax.set_title('E. Spawning Overhead', color=WHITE, fontweight='bold')
    ax.tick_params(colors=WHITE)
    ax.grid(True, alpha=0.3, axis='y')
    
    # Panel F: K×N stress test (KxN_stress)
    ax = axes[2, 1]
    ax.set_facecolor(PANEL_COLOR)
    
    if 'KxN_stress' in data:
        # Group by K values
        k_groups = {}
        for d in data['KxN_stress']:
            k = d['k']
            if k not in k_groups:
                k_groups[k] = {'total': [], 'wall_s': []}
            k_groups[k]['total'].append(d['total'])
            k_groups[k]['wall_s'].append(d['wall_s'])
        
        colors = [PRIMARY, SECONDARY, TERTIARY, ACCENT]
        for i, (k, group) in enumerate(sorted(k_groups.items())):
            color = colors[i % len(colors)]
            ax.loglog(group['total'], group['wall_s'], 'o-', 
                     color=color, label=f'K={k}', markersize=6)
        
        ax.set_xlabel('Total Agents (K×N)', color=WHITE)
        ax.set_ylabel('Wall Time (s)', color=WHITE)
        ax.legend(facecolor=PANEL_COLOR, edgecolor=WHITE, labelcolor=WHITE)
    else:
        ax.text(0.5, 0.5, 'No Data\n(KxN_stress)', ha='center', va='center', 
               color=GRAY, fontsize=16, transform=ax.transAxes)
    
    ax.set_title('F. Spatial×Population Stress', color=WHITE, fontweight='bold')
    ax.tick_params(colors=WHITE)
    ax.grid(True, alpha=0.3)
    
    # Panel G: Parallel efficiency (parallel)
    ax = axes[3, 0]
    ax.set_facecolor(PANEL_COLOR)
    
    if 'parallel' in data:
        procs = [d['n_procs'] for d in data['parallel']]
        throughput = [d['throughput'] for d in data['parallel']]
        efficiency = [d['efficiency_pct'] for d in data['parallel']]
        
        ax2 = ax.twinx()
        
        # Throughput bars
        bars = ax.bar(procs, throughput, color=PRIMARY, alpha=0.7, label='Throughput (runs/s)')
        ax.set_xlabel('Number of Processes', color=WHITE)
        ax.set_ylabel('Throughput (runs/s)', color=PRIMARY)
        
        # Efficiency line
        line = ax2.plot(procs, efficiency, 'o-', color=SECONDARY, linewidth=2, 
                       markersize=8, label='Efficiency %')
        ax2.set_ylabel('Parallel Efficiency (%)', color=SECONDARY)
        ax2.axhline(y=100, color=WHITE, linestyle='--', alpha=0.5)
        
        # Combine legends
        lines1, labels1 = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax.legend(lines1 + lines2, labels1 + labels2, 
                 facecolor=PANEL_COLOR, edgecolor=WHITE, labelcolor=WHITE)
        
        ax.tick_params(colors=WHITE)
        ax2.tick_params(colors=WHITE)
    else:
        ax.text(0.5, 0.5, 'No Data\n(parallel)', ha='center', va='center', 
               color=GRAY, fontsize=16, transform=ax.transAxes)
    
    ax.set_title('G. Parallel Scaling', color=WHITE, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    
    # Panel H: Memory usage vs total agents
    ax = axes[3, 1]
    ax.set_facecolor(PANEL_COLOR)
    
    # Combine memory data from all sources
    mem_data = []
    for source_name, source_data in data.items():
        for d in source_data:
            if 'mem_total_mb' in d and 'n' in d:
                total_agents = d.get('total', d['n'])  # Use total if available, else n
                mem_data.append((total_agents, d['mem_total_mb']))
            elif 'mem_delta_mb' in d and 'n' in d:
                total_agents = d.get('total', d['n'])
                # Estimate total from delta (rough)
                mem_data.append((total_agents, d['mem_delta_mb'] + 500))  # Base + delta
    
    if mem_data:
        agents, memory = zip(*sorted(mem_data))
        ax.loglog(agents, memory, 'o', color=PRIMARY, markersize=6, alpha=0.7)
        
        # 31GB line
        ax.axhline(y=31*1024, color=ACCENT, linestyle='--', linewidth=2, 
                  label='31GB RAM Limit', alpha=0.8)
        
        ax.set_xlabel('Total Agents', color=WHITE)
        ax.set_ylabel('Memory Usage (MB)', color=WHITE)
        ax.legend(facecolor=PANEL_COLOR, edgecolor=WHITE, labelcolor=WHITE)
    else:
        ax.text(0.5, 0.5, 'No Memory Data\nAvailable', ha='center', va='center', 
               color=GRAY, fontsize=16, transform=ax.transAxes)
    
    ax.set_title('H. Memory Scaling', color=WHITE, fontweight='bold')
    ax.tick_params(colors=WHITE)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.96])
    return fig

def create_projections_figure(data):
    """Create projections figure for key scenarios"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 8), facecolor=BG_COLOR)
    fig.suptitle('SSWD-EvoEpi: Runtime Projections for Key Scenarios', 
                 fontsize=16, color=WHITE, y=0.95)
    
    # Extract scaling parameters from single_N and single_T
    n_params = None
    t_params = None
    
    if 'single_N' in data:
        n_vals = [d['n'] for d in data['single_N']]
        wall_vals = [d['wall_s'] for d in data['single_N']]
        n_params, _ = fit_power_law(n_vals, wall_vals)
    
    if 'single_T' in data:
        t_vals = [d['years'] for d in data['single_T']]
        wall_vals = [d['wall_s'] for d in data['single_T']]
        t_params, _ = fit_power_law(t_vals, wall_vals)
    
    # Panel A: Runtime projections
    ax1.set_facecolor(PANEL_COLOR)
    
    scenarios = [
        ('Sensitivity Analysis\n1000×500×20yr', 1000*500*20, 1000, 500, 20),
        ('Full Coastline\n150×500×20yr', 150*500*20, 150, 500, 20),
        ('Monte Carlo\n10000×200×10yr', 10000*200*10, 10000, 200, 10),
        ('Century Run\n500×100yr', 500*100, 1, 500, 100)
    ]
    
    seq_times = []
    par_times = []
    
    for name, complexity, runs, n, t in scenarios:
        # Estimate single run time
        if n_params and t_params:
            # Combine N and T scaling: T_total ≈ T_n * T_t / T_baseline
            base_n, base_t = 500, 20  # Reference point
            t_n = power_law(n, *n_params[:2]) / power_law(base_n, *n_params[:2])
            t_t = power_law(t, *t_params[:2]) / power_law(base_t, *t_params[:2])
            single_run_s = t_n * t_t * 6.2  # 6.2s from 5-node reference
        else:
            # Fallback rough estimate
            single_run_s = (n/500) * (t/20) * 6.2
        
        total_seq_s = runs * single_run_s
        total_par_s = total_seq_s / 16  # 16 cores
        
        seq_times.append(total_seq_s)
        par_times.append(total_par_s)
    
    # Convert to hours/minutes
    seq_hours = [t/3600 for t in seq_times]
    par_hours = [t/3600 for t in par_times]
    
    x = np.arange(len(scenarios))
    width = 0.35
    
    bars1 = ax1.bar(x - width/2, seq_hours, width, label='Sequential', 
                    color=PRIMARY, alpha=0.8)
    bars2 = ax1.bar(x + width/2, par_hours, width, label='16-Core Parallel', 
                    color=SECONDARY, alpha=0.8)
    
    # Add value labels
    for i, (seq_h, par_h) in enumerate(zip(seq_hours, par_hours)):
        if seq_h < 1:
            seq_label = f'{seq_h*60:.0f}m'
        elif seq_h < 24:
            seq_label = f'{seq_h:.1f}h'
        else:
            seq_label = f'{seq_h/24:.1f}d'
            
        if par_h < 1:
            par_label = f'{par_h*60:.0f}m'
        elif par_h < 24:
            par_label = f'{par_h:.1f}h'
        else:
            par_label = f'{par_h/24:.1f}d'
            
        ax1.text(i - width/2, seq_h + 0.1, seq_label, ha='center', color=WHITE, fontsize=8)
        ax1.text(i + width/2, par_h + 0.1, par_label, ha='center', color=WHITE, fontsize=8)
    
    ax1.set_xticks(x)
    ax1.set_xticklabels([s[0] for s in scenarios], rotation=15, ha='right')
    ax1.set_ylabel('Runtime (hours)', color=WHITE)
    ax1.set_yscale('log')
    ax1.legend(facecolor=PANEL_COLOR, edgecolor=WHITE, labelcolor=WHITE)
    ax1.set_title('A. Projected Runtimes', color=WHITE, fontweight='bold')
    ax1.tick_params(colors=WHITE)
    ax1.grid(True, alpha=0.3, axis='y')
    
    # Panel B: Memory budget
    ax2.set_facecolor(PANEL_COLOR)
    
    # Memory scaling from collected data
    mem_agents = []
    mem_usage = []
    for source_name, source_data in data.items():
        for d in source_data:
            if 'mem_total_mb' in d and 'n' in d:
                total_agents = d.get('total', d['n'])
                mem_agents.append(total_agents)
                mem_usage.append(d['mem_total_mb'])
    
    if mem_agents:
        # Fit memory scaling
        mem_params, _ = fit_power_law(mem_agents, mem_usage, "Memory")
        
        # Project memory usage for scenarios
        scenario_names = ['Sensitivity\n1K runs', 'Full Coast\n150 nodes', 
                         'Monte Carlo\n10K runs', 'Century\n100yr']
        scenario_agents = [500, 500, 200, 500]  # Per-run agent counts
        
        proj_mem = []
        for n in scenario_agents:
            if mem_params:
                mem_mb = power_law(n, *mem_params[:2])
                proj_mem.append(mem_mb)
            else:
                proj_mem.append(n * 2)  # Rough estimate: 2MB per agent
        
        x = np.arange(len(scenario_names))
        bars = ax2.bar(x, proj_mem, color=TERTIARY, alpha=0.8)
        
        # 31GB limit line
        ax2.axhline(y=31*1024, color=ACCENT, linestyle='--', linewidth=2, 
                   label='31GB RAM Limit', alpha=0.8)
        
        # Memory labels
        for i, mem in enumerate(proj_mem):
            if mem < 1024:
                label = f'{mem:.0f}MB'
            else:
                label = f'{mem/1024:.1f}GB'
            ax2.text(i, mem + 100, label, ha='center', color=WHITE, fontsize=9)
        
        ax2.set_xticks(x)
        ax2.set_xticklabels(scenario_names)
        ax2.set_ylabel('Memory Usage (MB)', color=WHITE)
        ax2.legend(facecolor=PANEL_COLOR, edgecolor=WHITE, labelcolor=WHITE)
    else:
        ax2.text(0.5, 0.5, 'No Memory\nData Available', ha='center', va='center', 
               color=GRAY, fontsize=16, transform=ax2.transAxes)
    
    ax2.set_title('B. Memory Budget Analysis', color=WHITE, fontweight='bold')
    ax2.tick_params(colors=WHITE)
    ax2.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    return fig

def main():
    """Generate comprehensive scaling analysis"""
    print("Loading performance data...")
    data = load_all_data()
    
    print(f"Found {len(data)} axis datasets:")
    for name, dataset in data.items():
        print(f"  {name}: {len(dataset)} points")
    
    # Ensure output directory
    os.makedirs('results/performance', exist_ok=True)
    
    # Generate comprehensive figure
    print("\nGenerating comprehensive scaling figure...")
    fig1 = create_comprehensive_figure(data)
    fig1.savefig('results/performance/scaling_axes_comprehensive.png', 
                 dpi=300, bbox_inches='tight', facecolor=BG_COLOR)
    plt.close(fig1)
    print("Saved: scaling_axes_comprehensive.png")
    
    # Generate projections figure
    print("Generating projections figure...")
    fig2 = create_projections_figure(data)
    fig2.savefig('results/performance/scaling_projections.png', 
                 dpi=300, bbox_inches='tight', facecolor=BG_COLOR)
    plt.close(fig2)
    print("Saved: scaling_projections.png")
    
    print("\nScaling analysis figures complete!")

if __name__ == '__main__':
    main()