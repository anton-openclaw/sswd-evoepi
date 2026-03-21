#!/usr/bin/env python3
"""Fig 3.1: Infection prevalence heatmap — Time × Region (S→N)."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from fig_common import *

plt = setup_style()

data = load_npz()
populations = data["populations"]
infected = data["infected"]
sim_days = data["sim_days"]
site_names = data["site_names"]
regions_per_site = get_regions_per_site(site_names)

years = START_YEAR + sim_days / 365.25

# Aggregate to regions
pop_reg = aggregate_to_regions(populations, regions_per_site)
inf_reg = aggregate_to_regions(infected, regions_per_site)

# Prevalence = infected / population (avoid div by zero)
with np.errstate(divide='ignore', invalid='ignore'):
    prevalence = np.where(pop_reg > 0, inf_reg / pop_reg, 0.0)

fig, ax = plt.subplots(figsize=(14, 6))

# Use imshow with aspect='auto' — prevalence.T is (18 regions, 159 timesteps)
extent = [years[0], years[-1], -0.5, len(COASTLINE_ORDER) - 0.5]
im = ax.imshow(prevalence.T, aspect="auto", origin="lower", extent=extent,
               cmap="Reds", vmin=0, vmax=0.5, interpolation="nearest")

ax.set_yticks(range(len(COASTLINE_ORDER)))
ax.set_yticklabels(COASTLINE_ORDER, fontsize=9)
ax.set_xlabel("Year")
ax.set_ylabel("Region (S → N)")
ax.set_title("W325: Infection Prevalence by Region Over Time", fontsize=13, fontweight="bold")

# Annotate SSWD onset
ax.axvline(2013.5, color="black", ls="--", lw=1.2, alpha=0.7)
ax.text(2013.6, len(COASTLINE_ORDER) - 1.5, "SSWD onset", fontsize=8, va="top", ha="left",
        color="black", alpha=0.8)

cb = fig.colorbar(im, ax=ax, label="Infection Prevalence", shrink=0.8, pad=0.02)

fig.savefig(f"{FIGDIR}/fig_3_1_infection_heatmap.pdf")
plt.close()
print("✓ fig_3_1_infection_heatmap.pdf")
