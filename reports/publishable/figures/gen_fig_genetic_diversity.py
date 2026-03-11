#!/usr/bin/env python3
"""Fig 6: Additive genetic variance (Va) for resistance over time, all 18 regions, F01."""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from fig_style import (apply_style, load_summary, years_from_yearly,
                       REGIONS_NS, savefig)

def main():
    apply_style()
    data = load_summary()

    fig, ax = plt.subplots(figsize=(10, 5.5))

    cmap = cm.get_cmap('turbo', len(REGIONS_NS))
    for i, reg in enumerate(REGIONS_NS):
        rd = data['F01']['regions'][reg]
        va = rd.get('yearly_va_resistance_mean', [])
        if not va:
            continue
        years = years_from_yearly(len(va))
        ax.plot(years, va, color=cmap(i), linewidth=1.2, label=reg, alpha=0.85)

    ax.axvline(2025, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
    ax.set_xlabel('Year')
    ax.set_ylabel('Additive genetic variance ($V_a$) for resistance')
    ax.set_title('Genetic diversity erosion by region (F01)')
    ax.legend(loc='center left', bbox_to_anchor=(1.01, 0.5), frameon=True,
              framealpha=0.9, fontsize=7, ncol=1)
    ax.set_xlim(2012, 2050)
    ax.set_ylim(bottom=0)

    fig.tight_layout()
    savefig(fig, 'fig_genetic_diversity.png')

if __name__ == '__main__':
    main()
