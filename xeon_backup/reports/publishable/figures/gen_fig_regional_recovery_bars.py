#!/usr/bin/env python3
"""Fig 2: Regional recovery (%) grouped bar chart — 18 regions × 4 scenarios."""
import numpy as np
import matplotlib.pyplot as plt
from fig_style import (apply_style, load_summary, SCENARIO_COLORS,
                       SCENARIO_LABELS, REGIONS_NS, savefig)

def main():
    apply_style()
    data = load_summary()

    scenarios = ['F01', 'F02', 'F03', 'F04']
    n_regions = len(REGIONS_NS)
    n_sc = len(scenarios)
    bar_width = 0.18
    x = np.arange(n_regions)

    fig, ax = plt.subplots(figsize=(14, 5))

    for i, sc in enumerate(scenarios):
        vals = []
        for reg in REGIONS_NS:
            vals.append(data[sc]['regions'][reg]['recovery_pct'])
        offset = (i - n_sc / 2 + 0.5) * bar_width
        bars = ax.bar(x + offset, vals, bar_width, color=SCENARIO_COLORS[sc],
                      label=SCENARIO_LABELS[sc], edgecolor='white', linewidth=0.3)

    ax.set_xticks(x)
    ax.set_xticklabels(REGIONS_NS, rotation=45, ha='right')
    ax.set_ylabel('Recovery in 2050 (% of peak)')
    ax.set_title('Regional recovery by scenario')
    ax.legend(loc='upper left', frameon=True, framealpha=0.9)
    ax.set_xlim(-0.5, n_regions - 0.5)

    # Add a subtle horizontal line at 100% for reference
    ax.axhline(100, color='gray', linestyle=':', linewidth=0.5, alpha=0.5)

    fig.tight_layout()
    savefig(fig, 'fig_regional_recovery_bars.png')

if __name__ == '__main__':
    main()
