#!/usr/bin/env python3
"""Fig 9: Recovery % vs latitude scatter, all 4 scenarios with F01 trend line."""
import numpy as np
import matplotlib.pyplot as plt
from fig_style import (apply_style, load_summary, SCENARIO_COLORS,
                       SCENARIO_LABELS, REGIONS_NS, savefig)

def main():
    apply_style()
    data = load_summary()

    markers = {'F01': 'o', 'F02': 's', 'F03': '^', 'F04': 'D'}

    fig, ax = plt.subplots(figsize=(7, 5))

    for sc in ['F01', 'F02', 'F03', 'F04']:
        lats, recs = [], []
        for reg in REGIONS_NS:
            rd = data[sc]['regions'][reg]
            lats.append(rd['lat'])
            recs.append(rd['recovery_pct'])
        ax.scatter(lats, recs, c=SCENARIO_COLORS[sc], marker=markers[sc],
                   s=40, label=SCENARIO_LABELS[sc], alpha=0.8, edgecolors='white',
                   linewidths=0.3, zorder=3)

    # Trend line for F01
    lats_f01 = [data['F01']['regions'][r]['lat'] for r in REGIONS_NS]
    recs_f01 = [data['F01']['regions'][r]['recovery_pct'] for r in REGIONS_NS]
    z = np.polyfit(lats_f01, recs_f01, 1)
    p = np.poly1d(z)
    lat_range = np.linspace(min(lats_f01) - 1, max(lats_f01) + 1, 100)
    ax.plot(lat_range, p(lat_range), '--', color=SCENARIO_COLORS['F01'],
            linewidth=1, alpha=0.6, label=f'F01 trend (slope={z[0]:.2f})')

    ax.set_xlabel('Latitude (°N)')
    ax.set_ylabel('Recovery in 2050 (% of peak)')
    ax.set_title('Recovery vs. latitude')
    ax.legend(loc='best', frameon=True, framealpha=0.9, fontsize=7)
    ax.set_ylim(bottom=0)

    fig.tight_layout()
    savefig(fig, 'fig_recovery_vs_latitude.png')

if __name__ == '__main__':
    main()
