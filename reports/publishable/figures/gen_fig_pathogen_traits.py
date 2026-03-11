#!/usr/bin/env python3
"""Fig 7: Final (2050) pathogen traits — T_vbnc and v_local for all 18 regions, F01."""
import numpy as np
import matplotlib.pyplot as plt
from fig_style import apply_style, load_summary, REGIONS_NS, savefig

def main():
    apply_style()
    data = load_summary()

    t_vbnc = []
    v_local = []
    for reg in REGIONS_NS:
        rd = data['F01']['regions'][reg]
        t_vbnc.append(rd['final_T_vbnc_mean'])
        v_local.append(rd['final_v_local_mean'])

    x = np.arange(len(REGIONS_NS))

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 6), sharex=True)

    # Panel A: T_vbnc
    ax1.bar(x, t_vbnc, color='tab:blue', edgecolor='white', linewidth=0.3)
    ax1.set_ylabel('$T_{VBNC}$ (°C)')
    ax1.set_title('(a) Pathogen thermal threshold')

    # Panel B: v_local
    ax2.bar(x, v_local, color='tab:red', edgecolor='white', linewidth=0.3)
    ax2.set_ylabel('$v_{local}$')
    ax2.set_title('(b) Pathogen virulence')
    ax2.set_xticks(x)
    ax2.set_xticklabels(REGIONS_NS, rotation=45, ha='right')

    fig.suptitle('Final pathogen traits by region (F01, 2050)', fontsize=11, y=1.01)
    fig.tight_layout()
    savefig(fig, 'fig_pathogen_traits_final.png')

if __name__ == '__main__':
    main()
