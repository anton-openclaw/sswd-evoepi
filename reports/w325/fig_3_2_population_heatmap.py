#!/usr/bin/env python3
"""Fig 3.2: Population heatmap — fraction of K over time × region."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from fig_common import *

plt = setup_style()

data = load_npz()
populations = data["populations"]
sim_days = data["sim_days"]
site_names = data["site_names"]
regions_per_site = get_regions_per_site(site_names)
site_counts = get_region_site_counts(site_names)

years = START_YEAR + sim_days / 365.25

pop_reg = aggregate_to_regions(populations, regions_per_site)

# Fraction of K
pop_frac = np.zeros_like(pop_reg)
for j, reg in enumerate(COASTLINE_ORDER):
    n_sites = site_counts.get(reg, 1)
    pop_frac[:, j] = pop_reg[:, j] / (n_sites * K_PER_SITE)

fig, ax = plt.subplots(figsize=(14, 6))

extent = [years[0], years[-1], -0.5, len(COASTLINE_ORDER) - 0.5]
im = ax.imshow(pop_frac.T, aspect="auto", origin="lower", extent=extent,
               cmap="YlGn", vmin=0, vmax=1.0, interpolation="nearest")

ax.set_yticks(range(len(COASTLINE_ORDER)))
ax.set_yticklabels(COASTLINE_ORDER, fontsize=9)
ax.set_xlabel("Year")
ax.set_ylabel("Region (S → N)")
ax.set_title("W325: Population as Fraction of K by Region Over Time", fontsize=13, fontweight="bold")

ax.axvline(2013.5, color="black", ls="--", lw=1.2, alpha=0.7)
ax.text(2013.6, len(COASTLINE_ORDER) - 1.5, "SSWD onset", fontsize=8, va="top", ha="left",
        color="black", alpha=0.8)

cb = fig.colorbar(im, ax=ax, label="Population / K", shrink=0.8, pad=0.02)

fig.savefig(f"{FIGDIR}/fig_3_2_population_heatmap.pdf")
plt.close()
print("✓ fig_3_2_population_heatmap.pdf")
