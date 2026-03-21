#!/usr/bin/env python3
"""Plot observed SST (2013-2025) alongside CMIP6 projections (2026-2050).

Validates the transition at the 2025/2026 boundary is smooth.
Outputs: results/sst_validation/observed_vs_projected.png
"""

import csv
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
from datetime import datetime


def load_monthly(path):
    """Load a monthly CSV → dict of (year, month) → sst."""
    data = {}
    with open(path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            yr = int(row['year'])
            mo = int(row['month'])
            data[(yr, mo)] = float(row['sst'])
    return data


def monthly_to_arrays(data, year_start, year_end):
    """Convert monthly dict to parallel date/sst arrays."""
    dates = []
    ssts = []
    for yr in range(year_start, year_end + 1):
        for mo in range(1, 13):
            if (yr, mo) in data:
                dates.append(datetime(yr, mo, 15))
                ssts.append(data[(yr, mo)])
    return dates, ssts


def main():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    obs_dir = os.path.join(base_dir, 'data', 'sst')
    proj_dir = os.path.join(base_dir, 'data', 'sst', 'projections')
    out_dir = os.path.join(base_dir, 'results', 'sst_validation')
    os.makedirs(out_dir, exist_ok=True)

    nodes = ['Sitka', 'SJI', 'Monterey']
    scenarios = ['ssp245', 'ssp585']

    fig, axes = plt.subplots(3, 1, figsize=(14, 12), sharex=True)
    fig.suptitle('Observed SST (OISST) vs CMIP6 Projections\n11-Node Stepping-Stone Network',
                 fontsize=14, fontweight='bold')

    colors_obs = '#2166AC'
    colors_proj = {'ssp245': '#E67E22', 'ssp585': '#C0392B'}
    labels_proj = {'ssp245': 'SSP2-4.5 (CMIP6 ensemble)', 'ssp585': 'SSP5-8.5 (placeholder)'}

    for idx, node in enumerate(nodes):
        ax = axes[idx]

        # Load observations
        obs_path = os.path.join(obs_dir, f'{node}_monthly.csv')
        obs_data = load_monthly(obs_path)
        obs_dates, obs_ssts = monthly_to_arrays(obs_data, 2013, 2025)

        # Plot observations
        ax.plot(obs_dates, obs_ssts, color=colors_obs, linewidth=1.2,
                label='OISST Observations', zorder=3)

        # Load and plot projections
        for scenario in scenarios:
            # Try exact, then placeholder
            proj_path = os.path.join(proj_dir, f'{node}_{scenario}_monthly.csv')
            if not os.path.exists(proj_path):
                proj_path = os.path.join(proj_dir, f'{node}_{scenario}_placeholder_monthly.csv')
            if not os.path.exists(proj_path):
                continue

            proj_data = load_monthly(proj_path)
            proj_dates, proj_ssts = monthly_to_arrays(proj_data, 2026, 2050)

            ax.plot(proj_dates, proj_ssts, color=colors_proj[scenario],
                    linewidth=1.0, alpha=0.85, label=labels_proj[scenario], zorder=2)

        # Mark the transition boundary
        ax.axvline(datetime(2026, 1, 1), color='gray', linestyle='--',
                   linewidth=0.8, alpha=0.7, zorder=1)
        ax.text(datetime(2026, 3, 1), ax.get_ylim()[0] if idx > 0 else 4,
                '← obs | proj →', fontsize=8, color='gray', va='bottom')

        # Compute and annotate the transition jump for ssp245
        if (2025, 12) in obs_data:
            last_obs = obs_data[(2025, 12)]
            proj245_path = os.path.join(proj_dir, f'{node}_ssp245_monthly.csv')
            if os.path.exists(proj245_path):
                proj245 = load_monthly(proj245_path)
                if (2026, 1) in proj245:
                    first_proj = proj245[(2026, 1)]
                    jump = first_proj - last_obs
                    ax.annotate(f'Δ = {jump:+.2f}°C',
                               xy=(datetime(2026, 1, 15), first_proj),
                               xytext=(datetime(2027, 6, 1), first_proj + 1.5),
                               fontsize=8, color=colors_proj['ssp245'],
                               arrowprops=dict(arrowstyle='->', color=colors_proj['ssp245'],
                                               lw=0.8))

        ax.set_ylabel('SST (°C)', fontsize=11)
        ax.set_title(f'{node}', fontsize=12, fontweight='bold', loc='left')
        ax.grid(True, alpha=0.3)
        if idx == 0:
            ax.legend(loc='upper left', fontsize=9, framealpha=0.9)

    axes[-1].xaxis.set_major_locator(mdates.YearLocator(5))
    axes[-1].xaxis.set_minor_locator(mdates.YearLocator(1))
    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    axes[-1].set_xlabel('Year', fontsize=11)

    plt.tight_layout()
    out_path = os.path.join(out_dir, 'observed_vs_projected.png')
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f'Saved: {out_path}')

    # Print transition statistics
    print('\nTransition Statistics (Dec 2025 → Jan 2026):')
    print('-' * 50)
    for node in nodes:
        obs_data = load_monthly(os.path.join(obs_dir, f'{node}_monthly.csv'))
        if (2025, 12) not in obs_data:
            print(f'  {node}: no Dec 2025 observation')
            continue
        last_obs = obs_data[(2025, 12)]
        for scenario in scenarios:
            proj_path = os.path.join(proj_dir, f'{node}_{scenario}_monthly.csv')
            if not os.path.exists(proj_path):
                proj_path = os.path.join(proj_dir, f'{node}_{scenario}_placeholder_monthly.csv')
            if not os.path.exists(proj_path):
                continue
            proj_data = load_monthly(proj_path)
            if (2026, 1) in proj_data:
                first_proj = proj_data[(2026, 1)]
                jump = first_proj - last_obs
                print(f'  {node} ({scenario}): {last_obs:.2f} → {first_proj:.2f} (Δ = {jump:+.2f}°C)')


if __name__ == '__main__':
    main()
