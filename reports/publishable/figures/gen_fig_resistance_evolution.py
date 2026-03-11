#!/usr/bin/env python3
"""Fig 5: Mean host resistance (r) evolution 2012-2050 for 6 key regions, F01."""
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
        years = years_from_yearly(len(rd['yearly_resistance_mean']))
        ax.plot(years, rd['yearly_resistance_mean'], color=cmap(i),
                linewidth=1.5, label=reg, marker='o', markersize=2)

    ax.axvline(2025, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
    ax.text(2013.5, ax.get_ylim()[1] * 0.97, 'Calibration', fontsize=7, color='gray')
    ax.text(2026, ax.get_ylim()[1] * 0.97, 'Forecast', fontsize=7, color='gray')

    ax.set_xlabel('Year')
    ax.set_ylabel('Mean host resistance ($\\bar{r}$)')
    ax.set_title('Host resistance evolution (F01)')
    ax.legend(loc='best', frameon=True, framealpha=0.9)
    ax.set_xlim(2012, 2050)

    fig.tight_layout()
    savefig(fig, 'fig_resistance_evolution_forecast.png')

if __name__ == '__main__':
    main()
