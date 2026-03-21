#!/usr/bin/env python3
"""Fig 10: Population heatmap — time × regions, F01 only."""
import numpy as np
import matplotlib.pyplot as plt
from fig_style import (apply_style, load_monthly_regional, SEEDS,
                       REGIONS_NS, savefig)

def main():
    apply_style()

    # Load all seeds and average
    all_data = []
    for seed in SEEDS:
        regional = load_monthly_regional('F01', seed)
        years = regional['_years']
        mat = np.array([regional[reg].astype(float) for reg in REGIONS_NS])
        all_data.append(mat)
    all_data = np.array(all_data)  # (3, 18, T)
    mean_data = all_data.mean(axis=0)  # (18, T)

    # Normalize each region by its peak (row-wise)
    peak = mean_data.max(axis=1, keepdims=True)
    peak[peak == 0] = 1
    frac = mean_data / peak

    fig, ax = plt.subplots(figsize=(14, 6))

    # Use imshow with extent for correct axis mapping
    extent = [years[0], years[-1], len(REGIONS_NS) - 0.5, -0.5]
    im = ax.imshow(frac, aspect='auto', cmap='RdYlGn', vmin=0, vmax=1,
                   extent=extent, interpolation='nearest')

    ax.set_yticks(np.arange(len(REGIONS_NS)))
    ax.set_yticklabels(REGIONS_NS)
    ax.set_xlabel('Year')
    ax.set_title('Population fraction of peak (F01)')

    # Calibration/forecast line
    ax.axvline(2025, color='black', linestyle='--', linewidth=0.8, alpha=0.7)

    cbar = fig.colorbar(im, ax=ax, shrink=0.8, label='Fraction of peak')

    fig.tight_layout()
    savefig(fig, 'fig_heatmap_population.png')

if __name__ == '__main__':
    main()
