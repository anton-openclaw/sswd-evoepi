#!/usr/bin/env python3
"""Fig 8: Evolution comparison — F01 vs F02 vs F03 recovery %, grouped bars."""
import numpy as np
import matplotlib.pyplot as plt
from fig_style import (apply_style, load_summary, SCENARIO_COLORS,
                       SCENARIO_LABELS, REGIONS_NS, savefig)

def main():
    apply_style()
    data = load_summary()

    scenarios = ['F01', 'F02', 'F03']
    n_regions = len(REGIONS_NS)
    n_sc = len(scenarios)
    bar_width = 0.25
    x = np.arange(n_regions)

    fig, ax = plt.subplots(figsize=(14, 5))

    for i, sc in enumerate(scenarios):
        vals = [data[sc]['regions'][reg]['recovery_pct'] for reg in REGIONS_NS]
        offset = (i - n_sc / 2 + 0.5) * bar_width
        ax.bar(x + offset, vals, bar_width, color=SCENARIO_COLORS[sc],
               label=SCENARIO_LABELS[sc], edgecolor='white', linewidth=0.3)

    ax.set_xticks(x)
    ax.set_xticklabels(REGIONS_NS, rotation=45, ha='right')
    ax.set_ylabel('Recovery in 2050 (% of peak)')
    ax.set_title('Effect of pathogen evolution on recovery')
    ax.legend(loc='upper left', frameon=True, framealpha=0.9)

    fig.tight_layout()
    savefig(fig, 'fig_evo_comparison.png')

if __name__ == '__main__':
    main()
