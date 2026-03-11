#!/usr/bin/env python3
"""Fig 1: Coast-wide population trajectories (2012-2050), all 4 scenarios."""
import numpy as np
import matplotlib.pyplot as plt
from fig_style import (apply_style, load_summary, years_from_yearly,
                       SCENARIO_COLORS, SCENARIO_LABELS, REGIONS_NS, savefig, FIG_DIR)

def main():
    apply_style()
    data = load_summary()

    fig, ax = plt.subplots(figsize=(7, 4))

    for sc in ['F01', 'F02', 'F03', 'F04']:
        sc_data = data[sc]['regions']
        # Sum across all regions for each year
        n_years = len(sc_data[REGIONS_NS[0]]['yearly_pop_mean'])
        years = years_from_yearly(n_years)

        # We need per-seed totals to compute std. Use the summary yearly_pop_mean/std per region.
        # Since we only have mean per region, sum means for coastwide mean.
        # For std, we need to go back to raw data. Use monthly npz instead.
        total_mean = np.zeros(n_years)
        for reg in REGIONS_NS:
            total_mean += np.array(sc_data[reg]['yearly_pop_mean'])

        # For std bands, load per-seed
        from fig_style import load_monthly_regional, SEEDS
        seed_coastwide = []
        for seed in SEEDS:
            regional = load_monthly_regional(sc, seed)
            monthly_years = regional['_years']
            # Aggregate to yearly by taking population at each year mark
            coastwide_monthly = np.zeros(len(monthly_years))
            for reg in REGIONS_NS:
                coastwide_monthly += regional[reg]
            # Resample to yearly (take value at each year boundary)
            yearly_vals = []
            for y in years:
                mask = (monthly_years >= y) & (monthly_years < y + 1)
                if mask.any():
                    yearly_vals.append(coastwide_monthly[mask].mean())
                else:
                    yearly_vals.append(0)
            seed_coastwide.append(yearly_vals)

        seed_arr = np.array(seed_coastwide) / 1e6  # millions
        mean = seed_arr.mean(axis=0)
        std = seed_arr.std(axis=0)

        color = SCENARIO_COLORS[sc]
        label = SCENARIO_LABELS[sc]
        ax.plot(years, mean, color=color, label=label, linewidth=1.5)
        ax.fill_between(years, mean - std, mean + std, color=color, alpha=0.15)

    # Mark calibration → forecast transition
    ax.axvline(2025, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
    ax.text(2025.3, ax.get_ylim()[1] * 0.95, 'Forecast →', fontsize=7, color='gray',
            va='top')

    ax.set_xlabel('Year')
    ax.set_ylabel('Total population (millions)')
    ax.set_title('Coast-wide population trajectory')
    ax.legend(loc='upper right', frameon=True, framealpha=0.9)
    ax.set_xlim(2012, 2050)
    ax.set_ylim(bottom=0)

    savefig(fig, 'fig_coastwide_population.png')

if __name__ == '__main__':
    main()
