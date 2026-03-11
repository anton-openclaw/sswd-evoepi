#!/usr/bin/env python3
"""Fig 6: Additive genetic variance (Va) for resistance over time, 6 key regions, F01."""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from fig_style import (apply_style, load_summary, years_from_yearly,
                       KEY_REGIONS, savefig)

def main():
    apply_style()
    data = load_summary()

    fig, ax = plt.subplots(figsize=(7, 4.5))

    cmap = cm.get_cmap('viridis', len(KEY_REGIONS))
    for i, reg in enumerate(KEY_REGIONS):
        rd = data['F01']['regions'][reg]
        va = rd['yearly_va_resistance_mean']
        years = years_from_yearly(len(va))
        ax.plot(years, va, color=cmap(i), linewidth=1.5, label=reg,
                marker='s', markersize=2)

    ax.axvline(2025, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
    ax.set_xlabel('Year')
    ax.set_ylabel('Additive genetic variance ($V_a$) for resistance')
    ax.set_title('Genetic diversity erosion (F01)')
    ax.legend(loc='best', frameon=True, framealpha=0.9)
    ax.set_xlim(2012, 2050)
    ax.set_ylim(bottom=0)

    fig.tight_layout()
    savefig(fig, 'fig_genetic_diversity.png')

if __name__ == '__main__':
    main()
