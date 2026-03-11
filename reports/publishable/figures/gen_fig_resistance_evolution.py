#!/usr/bin/env python3
"""Fig 5: Mean host resistance (r) evolution 2012-2050 for all 18 regions, F01."""
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
        r_data = rd.get('yearly_resistance_mean', [])
        if not r_data:
            continue
        years = years_from_yearly(len(r_data))
        ax.plot(years, r_data, color=cmap(i),
                linewidth=1.2, label=reg, alpha=0.85)

    ax.axvline(2025, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
    ax.text(2013.5, ax.get_ylim()[1] * 0.97, 'Calibration', fontsize=7, color='gray')
    ax.text(2026, ax.get_ylim()[1] * 0.97, 'Forecast', fontsize=7, color='gray')

    ax.set_xlabel('Year')
    ax.set_ylabel('Mean host resistance ($\\bar{r}$)')
    ax.set_title('Host resistance evolution by region (F01)')
    ax.legend(loc='center left', bbox_to_anchor=(1.01, 0.5), frameon=True,
              framealpha=0.9, fontsize=7, ncol=1)
    ax.set_xlim(2012, 2050)

    fig.tight_layout()
    savefig(fig, 'fig_resistance_evolution_forecast.png')

if __name__ == '__main__':
    main()
