#!/usr/bin/env python3
"""Fig 22: Four snapshots of range-wide disease simulation from F01 seed42."""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
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

    # Four snapshot timepoints
    snapshots = [
        ('Pre-disease (2012)', 0),                           # t=0
        ('Peak crash (2014)', np.argmin(np.abs(years - 2014.5))),
        ('Partial recovery (2030)', np.argmin(np.abs(years - 2030))),
        ('Final state (2050)', len(years) - 1),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    for ax, (title, tidx) in zip(axes.flat, snapshots):
        pop = pops[tidx]
        inf = infected[tidx]
        frac_alive = pop / K  # fraction of K
        frac_infected = np.where(pop > 0, inf / pop, 0)

        # Size by population fraction
        sizes = 5 + 80 * (frac_alive ** 0.5)  # sqrt scale for visibility
        # Color by infection prevalence
        colors = frac_infected

        # Background: coastline hint
        ax.set_facecolor('#f0f4f8')

        # Plot all sites
        sc = ax.scatter(lons, lats, s=sizes, c=colors,
                       cmap='RdYlGn_r', vmin=0, vmax=1,
                       edgecolors='gray', linewidths=0.3, alpha=0.85,
                       zorder=2)

        # Mark extinct sites
        extinct = pop == 0
        if extinct.any():
            ax.scatter(lons[extinct], lats[extinct], s=8, c='black',
                      marker='x', linewidths=0.5, alpha=0.4, zorder=1)

        ax.set_title(title, fontsize=11, fontweight='bold')
        ax.set_xlim(-180, -115)
        ax.set_ylim(28, 62)
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')

        # Stats annotation
        total_pop = pop.sum()
        total_inf = inf.sum()
        n_alive = (pop > 0).sum()
        prev = total_inf / total_pop if total_pop > 0 else 0
        ax.text(0.02, 0.02,
                f'Pop: {total_pop:,.0f} ({n_alive}/896 sites)\n'
                f'Prevalence: {prev:.1%}',
                transform=ax.transAxes, fontsize=7,
                verticalalignment='bottom',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

    # Shared colorbar
    cbar = fig.colorbar(sc, ax=axes, shrink=0.6, pad=0.04,
                        label='Infection prevalence (fraction of alive)')

    fig.suptitle('Disease simulation snapshots — F01 baseline (seed 42)',
                 fontsize=13)
    fig.subplots_adjust(top=0.93, hspace=0.25, wspace=0.25)
    savefig(fig, 'fig_disease_snapshots.png')

if __name__ == '__main__':
    main()
