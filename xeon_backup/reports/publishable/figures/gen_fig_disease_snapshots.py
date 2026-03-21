#!/usr/bin/env python3
"""Fig 22: Four snapshots of range-wide disease simulation from F01 seed42.

Proper map extent covering all 896 sites, extinct sites transparent,
colorbar on the right, clean publication styling.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from fig_style import apply_style, savefig, NPZ_DIR


def main():
    apply_style()
    d = np.load(NPZ_DIR / 'F01' / 'monthly_seed42.npz', allow_pickle=True)
    pops = d['populations']      # (T, 896)
    infected = d['infected']     # (T, 896)
    lats = d['site_lats']
    lons = d['site_lons']
    sim_days = d['sim_days']
    K = int(d['K'])

    years = 2012 + sim_days / 365.25

    # Compute map extent from actual site coordinates with padding
    lon_pad, lat_pad = 2.0, 1.5
    lon_min, lon_max = lons.min() - lon_pad, lons.max() + lon_pad
    lat_min, lat_max = lats.min() - lat_pad, lats.max() + lat_pad

    # Four snapshot timepoints
    snapshots = [
        ('Pre-disease (2012)', 0),
        ('Peak crash (mid-2014)', np.argmin(np.abs(years - 2014.5))),
        ('Mid-recovery (2030)', np.argmin(np.abs(years - 2030))),
        ('Final state (2050)', len(years) - 1),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    for ax, (title, tidx) in zip(axes.flat, snapshots):
        pop = pops[tidx]
        inf = infected[tidx]
        alive_mask = pop > 0
        extinct_mask = ~alive_mask

        # Infection prevalence for alive sites
        prevalence = np.zeros_like(pop, dtype=float)
        prevalence[alive_mask] = inf[alive_mask] / pop[alive_mask]

        # Population fraction of K for sizing
        frac_alive = pop / K

        # Background
        ax.set_facecolor('#e8ecf1')

        # Plot EXTINCT sites first (small, transparent gray x)
        if extinct_mask.any():
            ax.scatter(lons[extinct_mask], lats[extinct_mask],
                      s=4, c='#999999', marker='.', alpha=0.15, zorder=1)

        # Plot ALIVE sites: size by population, color by prevalence
        if alive_mask.any():
            sizes = 8 + 60 * (frac_alive[alive_mask] ** 0.5)
            sc = ax.scatter(lons[alive_mask], lats[alive_mask],
                           s=sizes, c=prevalence[alive_mask],
                           cmap='RdYlGn_r', vmin=0, vmax=1,
                           edgecolors='#333333', linewidths=0.2, alpha=0.9,
                           zorder=2)

        ax.set_xlim(lon_min, lon_max)
        ax.set_ylim(lat_min, lat_max)
        ax.set_title(title, fontsize=11, fontweight='bold')
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
        ax.set_aspect('equal', adjustable='box')

        # Stats annotation
        total_pop = pop.sum()
        total_inf = inf.sum()
        n_alive = alive_mask.sum()
        prev = total_inf / total_pop if total_pop > 0 else 0
        ax.text(0.02, 0.02,
                f'N = {total_pop:,.0f} ({n_alive}/896 sites alive)\n'
                f'Prevalence: {prev:.1%}',
                transform=ax.transAxes, fontsize=7,
                verticalalignment='bottom',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                         edgecolor='#cccccc', alpha=0.9))

    # Shared colorbar on the right side
    fig.subplots_adjust(right=0.88, top=0.92, hspace=0.30, wspace=0.25)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.70])
    cbar = fig.colorbar(sc, cax=cbar_ax)
    cbar.set_label('Infection prevalence', fontsize=10)

    # Legend for extinct sites
    legend_elements = [
        Line2D([0], [0], marker='.', color='w', markerfacecolor='#999999',
               markersize=5, alpha=0.4, label='Extinct site'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='green',
               markeredgecolor='#333333', markersize=8, label='Alive (healthy)'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='red',
               markeredgecolor='#333333', markersize=8, label='Alive (infected)'),
    ]
    axes[0, 0].legend(handles=legend_elements, loc='upper left', fontsize=7,
                      framealpha=0.9, edgecolor='#cccccc')

    fig.suptitle('Disease simulation snapshots — F01 baseline (seed 42)',
                 fontsize=13, fontweight='bold')
    savefig(fig, 'fig_disease_snapshots.png')


if __name__ == '__main__':
    main()
