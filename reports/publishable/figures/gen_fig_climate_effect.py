#!/usr/bin/env python3
"""Fig 4: Diverging bars — climate effect (F04−F01 recovery %) per region."""
import numpy as np
import matplotlib.pyplot as plt
from fig_style import apply_style, load_summary, REGIONS_NS, savefig

def main():
    apply_style()
    data = load_summary()
    ce = data['comparisons']['climate_effect']

    diffs = [ce[r]['diff_pct'] for r in REGIONS_NS]
    colors = ['#2ca02c' if d >= 0 else '#d62728' for d in diffs]

    fig, ax = plt.subplots(figsize=(7, 6))
    y = np.arange(len(REGIONS_NS))
    ax.barh(y, diffs, color=colors, edgecolor='white', linewidth=0.3, height=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels(REGIONS_NS)
    ax.invert_yaxis()  # N at top
    ax.axvline(0, color='black', linewidth=0.5)
    ax.set_xlabel('Δ Recovery % (SSP5-8.5 − SSP2-4.5)')
    ax.set_title('Climate effect on recovery')

    # Annotations
    ax.text(0.98, 0.02, 'Green = warming helps\nRed = warming hurts',
            transform=ax.transAxes, fontsize=7, ha='right', va='bottom',
            color='gray', style='italic')

    fig.tight_layout()
    savefig(fig, 'fig_climate_effect.png')

if __name__ == '__main__':
    main()
