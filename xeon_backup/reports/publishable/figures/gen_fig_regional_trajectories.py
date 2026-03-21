#!/usr/bin/env python3
"""Fig 3: Population trajectories for 6 key regions, all scenarios overlaid."""
import numpy as np
import matplotlib.pyplot as plt
from fig_style import (apply_style, load_monthly_regional, SEEDS,
                       SCENARIO_COLORS, SCENARIO_LABELS, KEY_REGIONS, savefig)

def main():
    apply_style()

    fig, axes = plt.subplots(2, 3, figsize=(14, 7), sharex=True)
    axes = axes.flatten()

    for ax_i, reg in enumerate(KEY_REGIONS):
        ax = axes[ax_i]
        for sc in ['F01', 'F02', 'F03', 'F04']:
            seed_series = []
            for seed in SEEDS:
                regional = load_monthly_regional(sc, seed)
                seed_series.append(regional[reg].astype(float))
            years = load_monthly_regional(sc, SEEDS[0])['_years']
            arr = np.array(seed_series)
            mean = arr.mean(axis=0)
            std = arr.std(axis=0)
            color = SCENARIO_COLORS[sc]
            ax.plot(years, mean, color=color, linewidth=1.2,
                    label=SCENARIO_LABELS[sc] if ax_i == 0 else None)
            ax.fill_between(years, mean - std, mean + std, color=color, alpha=0.12)

        ax.axvline(2025, color='gray', linestyle='--', linewidth=0.6, alpha=0.6)
        ax.set_title(reg, fontsize=10, fontweight='bold')
        ax.set_ylim(bottom=0)
        if ax_i >= 3:
            ax.set_xlabel('Year')
        if ax_i % 3 == 0:
            ax.set_ylabel('Population')

    # Single legend
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=4, frameon=True,
              bbox_to_anchor=(0.5, -0.02))
    fig.suptitle('Regional population trajectories', fontsize=12, y=1.01)
    fig.tight_layout()
    savefig(fig, 'fig_regional_trajectories.png')

if __name__ == '__main__':
    main()
