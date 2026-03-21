#!/usr/bin/env python3
"""
SST Visualization Suite for SSWD-EvoEpi
Generates publication-quality temperature figures from OISST observations
and CMIP6 SSP2-4.5 projections for all 11 stepping-stone nodes.

Output: results/sst_figures/
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.lines import Line2D
from pathlib import Path

# Paths
BASE = Path(__file__).parent.parent
SST_DIR = BASE / "data" / "sst"
PROJ_DIR = SST_DIR / "projections"
OUT_DIR = BASE / "results" / "sst_figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Node info: name, latitude (for ordering/coloring)
NODES = [
    ("Sitka", 57.05),
    ("Ketchikan", 55.34),
    ("Haida_Gwaii", 53.25),
    ("Bella_Bella", 52.16),
    ("Howe_Sound", 49.38),
    ("SJI", 48.53),
    ("Westport", 46.89),
    ("Newport", 44.63),
    ("Crescent_City", 41.74),
    ("Fort_Bragg", 39.43),
    ("Monterey", 36.60),
]

# Style
plt.rcParams.update({
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
    'font.size': 11,
    'axes.titlesize': 13,
    'axes.labelsize': 12,
    'legend.fontsize': 9,
    'figure.dpi': 150,
    'savefig.dpi': 200,
})


def load_observed(node_name):
    """Load OISST monthly observations."""
    path = SST_DIR / f"{node_name}_monthly.csv"
    if not path.exists():
        return None
    df = pd.read_csv(path)
    df['date'] = pd.to_datetime(df[['year', 'month']].assign(day=15))
    return df


def load_projection(node_name, scenario='ssp245'):
    """Load CMIP6 projection (bias-corrected ensemble mean)."""
    path = PROJ_DIR / f"{node_name}_{scenario}_monthly.csv"
    if not path.exists():
        return None
    df = pd.read_csv(path)
    df['date'] = pd.to_datetime(df[['year', 'month']].assign(day=15))
    return df


def load_climatology(node_name):
    """Load daily climatology."""
    path = SST_DIR / f"{node_name}_climatology.csv"
    if not path.exists():
        return None
    return pd.read_csv(path)


def get_lat_color(lat, lats):
    """Map latitude to color on a cool-warm gradient."""
    norm = (lat - min(lats)) / (max(lats) - min(lats))
    return cm.coolwarm(norm)


def fig1_all_nodes_timeseries():
    """Figure 1: Full timeseries (2002-2100) for all 11 nodes, observed + projected."""
    fig, ax = plt.subplots(figsize=(14, 7))
    lats = [n[1] for n in NODES]
    
    for name, lat in NODES:
        color = get_lat_color(lat, lats)
        obs = load_observed(name)
        proj = load_projection(name)
        
        if obs is not None:
            # Annual mean for cleaner plot
            annual = obs.groupby('year')['sst'].mean().reset_index()
            ax.plot(annual['year'], annual['sst'], color=color, linewidth=1.5, alpha=0.9)
        
        if proj is not None:
            annual = proj.groupby('year')['sst'].mean().reset_index()
            ax.plot(annual['year'], annual['sst'], color=color, linewidth=1.5, 
                    alpha=0.6, linestyle='--')
    
    # Transition line
    ax.axvline(x=2025.5, color='gray', linestyle=':', alpha=0.5, linewidth=1)
    ax.text(2025.5, ax.get_ylim()[1] * 0.98, ' Observed → Projected', 
            fontsize=9, color='gray', va='top')
    
    # SSWD onset
    ax.axvline(x=2013, color='red', linestyle=':', alpha=0.4, linewidth=1)
    ax.text(2013.2, ax.get_ylim()[0] + 0.3, 'SSWD\nonset', fontsize=8, color='red', alpha=0.6)
    
    # Legend: latitude gradient
    legend_nodes = [NODES[0], NODES[5], NODES[-1]]
    handles = [Line2D([0], [0], color=get_lat_color(lat, lats), linewidth=2, label=f"{name} ({lat}°N)")
               for name, lat in legend_nodes]
    handles.append(Line2D([0], [0], color='black', linewidth=1.5, label='Observed'))
    handles.append(Line2D([0], [0], color='black', linewidth=1.5, linestyle='--', alpha=0.6, label='SSP2-4.5'))
    ax.legend(handles=handles, loc='upper left', framealpha=0.9)
    
    ax.set_xlabel('Year')
    ax.set_ylabel('Annual Mean SST (°C)')
    ax.set_title('Sea Surface Temperature: Observations (2002-2025) & CMIP6 SSP2-4.5 Projections (2026-2100)')
    ax.set_xlim(2002, 2100)
    ax.grid(True, alpha=0.2)
    
    fig.savefig(OUT_DIR / "fig1_all_nodes_timeseries.png")
    plt.close(fig)
    print("  ✓ fig1_all_nodes_timeseries.png")


def fig2_latitude_heatmap():
    """Figure 2: Latitude × Year heatmap of annual mean SST."""
    years_obs = range(2002, 2026)
    years_proj = range(2026, 2101)
    all_years = list(years_obs) + list(years_proj)
    
    sst_matrix = np.full((len(NODES), len(all_years)), np.nan)
    
    for i, (name, lat) in enumerate(NODES):
        obs = load_observed(name)
        proj = load_projection(name)
        
        if obs is not None:
            annual = obs.groupby('year')['sst'].mean()
            for yr in years_obs:
                if yr in annual.index:
                    sst_matrix[i, yr - 2002] = annual[yr]
        
        if proj is not None:
            annual = proj.groupby('year')['sst'].mean()
            for yr in years_proj:
                if yr in annual.index:
                    sst_matrix[i, yr - 2002] = annual[yr]
    
    fig, ax = plt.subplots(figsize=(16, 6))
    im = ax.imshow(sst_matrix, aspect='auto', cmap='RdYlBu_r',
                   extent=[2002, 2100, len(NODES) - 0.5, -0.5],
                   interpolation='nearest')
    
    ax.set_yticks(range(len(NODES)))
    ax.set_yticklabels([f"{name} ({lat}°N)" for name, lat in NODES], fontsize=10)
    
    # Transition line
    ax.axvline(x=2025.5, color='white', linestyle='-', linewidth=2, alpha=0.8)
    
    cbar = fig.colorbar(im, ax=ax, shrink=0.8, label='Annual Mean SST (°C)')
    
    ax.set_xlabel('Year')
    ax.set_title('SST Across Latitude: OISST Observations + CMIP6 SSP2-4.5 Projections')
    
    fig.savefig(OUT_DIR / "fig2_latitude_heatmap.png")
    plt.close(fig)
    print("  ✓ fig2_latitude_heatmap.png")


def fig3_seasonal_cycle():
    """Figure 3: Monthly climatology by node — seasonal amplitude and timing."""
    fig, ax = plt.subplots(figsize=(12, 7))
    lats = [n[1] for n in NODES]
    months = np.arange(1, 13)
    month_labels = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']
    
    for name, lat in NODES:
        obs = load_observed(name)
        if obs is None:
            continue
        color = get_lat_color(lat, lats)
        monthly_mean = obs.groupby('month')['sst'].mean()
        monthly_std = obs.groupby('month')['sst'].std()
        
        ax.plot(months, monthly_mean.values, color=color, linewidth=2, label=f"{name}")
        ax.fill_between(months, 
                        (monthly_mean - monthly_std).values,
                        (monthly_mean + monthly_std).values,
                        color=color, alpha=0.1)
    
    ax.set_xticks(months)
    ax.set_xticklabels(month_labels)
    ax.set_xlabel('Month')
    ax.set_ylabel('SST (°C)')
    ax.set_title('Seasonal SST Cycle by Node (2002-2025 OISST, mean ± 1 SD)')
    ax.legend(loc='upper left', ncol=2, fontsize=9, framealpha=0.9)
    ax.grid(True, alpha=0.2)
    
    fig.savefig(OUT_DIR / "fig3_seasonal_cycle.png")
    plt.close(fig)
    print("  ✓ fig3_seasonal_cycle.png")


def fig4_warming_trend():
    """Figure 4: Observed warming trends (2002-2025) and projected (2026-2100) per node."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
    lats = [n[1] for n in NODES]
    
    obs_trends = []
    proj_trends = []
    
    for name, lat in NODES:
        obs = load_observed(name)
        proj = load_projection(name)
        
        if obs is not None:
            annual = obs.groupby('year')['sst'].mean().reset_index()
            if len(annual) > 2:
                z = np.polyfit(annual['year'], annual['sst'], 1)
                obs_trends.append((name, lat, z[0] * 10))  # °C per decade
        
        if proj is not None:
            annual = proj.groupby('year')['sst'].mean().reset_index()
            if len(annual) > 2:
                z = np.polyfit(annual['year'], annual['sst'], 1)
                proj_trends.append((name, lat, z[0] * 10))
    
    # Observed trends
    ax = axes[0]
    for name, lat, trend in obs_trends:
        color = get_lat_color(lat, lats)
        ax.barh(name, trend, color=color, alpha=0.8, edgecolor='white', linewidth=0.5)
    ax.axvline(x=0, color='gray', linewidth=0.5)
    ax.set_xlabel('Warming Rate (°C / decade)')
    ax.set_title('Observed (2002-2025)')
    ax.grid(True, alpha=0.2, axis='x')
    
    # Projected trends
    ax = axes[1]
    for name, lat, trend in proj_trends:
        color = get_lat_color(lat, lats)
        ax.barh(name, trend, color=color, alpha=0.8, edgecolor='white', linewidth=0.5)
    ax.axvline(x=0, color='gray', linewidth=0.5)
    ax.set_xlabel('Warming Rate (°C / decade)')
    ax.set_title('Projected SSP2-4.5 (2026-2100)')
    ax.grid(True, alpha=0.2, axis='x')
    
    fig.suptitle('Warming Trends by Node', fontsize=14, y=1.02)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig4_warming_trends.png")
    plt.close(fig)
    print("  ✓ fig4_warming_trends.png")


def fig5_transition_detail():
    """Figure 5: Detailed view of observed-to-projected transition (2020-2035) for 3 key nodes."""
    key_nodes = ["Sitka", "SJI", "Monterey"]
    fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
    
    for ax, name in zip(axes, key_nodes):
        lat = dict(NODES)[name]
        obs = load_observed(name)
        proj = load_projection(name)
        
        if obs is not None:
            mask = obs['year'] >= 2020
            monthly = obs[mask]
            ax.plot(monthly['date'], monthly['sst'], color='#2166ac', linewidth=1.2, 
                    label='OISST Observed', alpha=0.9)
        
        if proj is not None:
            mask = proj['year'] <= 2035
            monthly = proj[mask]
            ax.plot(monthly['date'], monthly['sst'], color='#b2182b', linewidth=1.2,
                    label='CMIP6 SSP2-4.5 (bias-corrected)', alpha=0.8)
        
        # Transition
        ax.axvline(x=pd.Timestamp('2025-12-15'), color='gray', linestyle=':', alpha=0.6)
        
        ax.set_ylabel('SST (°C)')
        ax.set_title(f'{name} ({lat}°N)', fontsize=12)
        ax.legend(loc='upper right', fontsize=9)
        ax.grid(True, alpha=0.2)
    
    axes[-1].set_xlabel('Date')
    fig.suptitle('Observed → Projected Transition Detail (2020-2035)', fontsize=14)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig5_transition_detail.png")
    plt.close(fig)
    print("  ✓ fig5_transition_detail.png")


def fig6_summer_max_trend():
    """Figure 6: Summer maximum SST over time — critical for disease dynamics."""
    fig, ax = plt.subplots(figsize=(14, 7))
    lats = [n[1] for n in NODES]
    
    for name, lat in NODES:
        color = get_lat_color(lat, lats)
        
        obs = load_observed(name)
        proj = load_projection(name)
        
        # Summer = months 6-9
        summer_max_obs = []
        summer_max_proj = []
        
        if obs is not None:
            summer = obs[obs['month'].isin([6, 7, 8, 9])]
            annual_max = summer.groupby('year')['sst'].max().reset_index()
            summer_max_obs = annual_max
            ax.plot(annual_max['year'], annual_max['sst'], color=color, 
                    linewidth=1.5, alpha=0.9)
        
        if proj is not None:
            summer = proj[proj['month'].isin([6, 7, 8, 9])]
            annual_max = summer.groupby('year')['sst'].max().reset_index()
            summer_max_proj = annual_max
            ax.plot(annual_max['year'], annual_max['sst'], color=color,
                    linewidth=1.5, alpha=0.6, linestyle='--')
    
    ax.axvline(x=2025.5, color='gray', linestyle=':', alpha=0.5)
    ax.axhline(y=20, color='red', linestyle=':', alpha=0.3, linewidth=2)
    ax.text(2003, 20.2, 'T_ref = 20°C (Prentice disease calibration)', 
            fontsize=9, color='red', alpha=0.5)
    
    # Legend
    legend_nodes = [NODES[0], NODES[5], NODES[-1]]
    handles = [Line2D([0], [0], color=get_lat_color(lat, lats), linewidth=2, 
               label=f"{name} ({lat}°N)") for name, lat in legend_nodes]
    handles.append(Line2D([0], [0], color='black', linewidth=1.5, label='Observed'))
    handles.append(Line2D([0], [0], color='black', linewidth=1.5, linestyle='--', alpha=0.6, label='SSP2-4.5'))
    ax.legend(handles=handles, loc='upper left', framealpha=0.9)
    
    ax.set_xlabel('Year')
    ax.set_ylabel('Peak Summer SST (°C)')
    ax.set_title('Summer Maximum SST (Jun-Sep): Disease Pressure Indicator')
    ax.set_xlim(2002, 2100)
    ax.grid(True, alpha=0.2)
    
    fig.savefig(OUT_DIR / "fig6_summer_max_trend.png")
    plt.close(fig)
    print("  ✓ fig6_summer_max_trend.png")


def fig7_disease_window():
    """Figure 7: Days per year above disease-critical temperature thresholds."""
    thresholds = [12, 15, 18]  # °C — relevant for Vibrio dynamics
    
    fig, axes = plt.subplots(1, 3, figsize=(16, 6), sharey=True)
    lats = [n[1] for n in NODES]
    
    for ax, thresh in zip(axes, thresholds):
        for name, lat in NODES:
            color = get_lat_color(lat, lats)
            obs = load_observed(name)
            proj = load_projection(name)
            
            days_above = []
            years = []
            
            # For monthly data, estimate days above threshold
            # Approximate: if monthly mean > threshold, count ~30 days; 
            # if within 2°C, interpolate
            for df, yr_range in [(obs, range(2002, 2026)), (proj, range(2026, 2101))]:
                if df is None:
                    continue
                for yr in yr_range:
                    yr_data = df[df['year'] == yr]
                    if len(yr_data) == 0:
                        continue
                    # Count months where mean SST exceeds threshold
                    months_above = (yr_data['sst'] >= thresh).sum()
                    # Rough: 30 days per month
                    days_above.append(months_above * 30.4)
                    years.append(yr)
            
            if years:
                style = '-' if years[0] < 2026 else '--'
                # Split into obs and proj for different line styles
                obs_mask = [y < 2026 for y in years]
                proj_mask = [y >= 2026 for y in years]
                
                obs_yrs = [y for y, m in zip(years, obs_mask) if m]
                obs_days = [d for d, m in zip(days_above, obs_mask) if m]
                proj_yrs = [y for y, m in zip(years, proj_mask) if m]
                proj_days = [d for d, m in zip(days_above, proj_mask) if m]
                
                if obs_yrs:
                    ax.plot(obs_yrs, obs_days, color=color, linewidth=1.2, alpha=0.8)
                if proj_yrs:
                    ax.plot(proj_yrs, proj_days, color=color, linewidth=1.2, 
                            alpha=0.5, linestyle='--')
        
        ax.axvline(x=2025.5, color='gray', linestyle=':', alpha=0.4)
        ax.set_xlabel('Year')
        ax.set_title(f'Days ≥ {thresh}°C')
        ax.grid(True, alpha=0.2)
        ax.set_xlim(2002, 2100)
    
    axes[0].set_ylabel('Estimated Days per Year')
    
    fig.suptitle('Disease Window: Days Above Temperature Thresholds', fontsize=14)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig7_disease_window.png")
    plt.close(fig)
    print("  ✓ fig7_disease_window.png")


if __name__ == "__main__":
    print("Generating SST visualizations...")
    print()
    
    fig1_all_nodes_timeseries()
    fig2_latitude_heatmap()
    fig3_seasonal_cycle()
    fig4_warming_trend()
    fig5_transition_detail()
    fig6_summer_max_trend()
    fig7_disease_window()
    
    print()
    print(f"All figures saved to {OUT_DIR}/")
