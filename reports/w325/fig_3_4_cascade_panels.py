#!/usr/bin/env python3
"""Fig 3.4: 4-panel cascade — population & infected by latitude band."""
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

pop_reg = aggregate_to_regions(populations, regions_per_site)
inf_reg = aggregate_to_regions(infected, regions_per_site)

# Define latitude bands
bands = {
    "Alaska": ["AK-FS", "AK-FN", "AK-OC", "AK-PWS", "AK-EG", "AK-WG", "AK-AL"],
    "BC / Salish Sea": ["BC-C", "BC-N", "SS-S", "SS-N"],
    "WA / OR": ["WA-O", "JDF", "OR"],
    "CA / BJ": ["BJ", "CA-S", "CA-C", "CA-N"],
}

fig, axes = plt.subplots(4, 1, figsize=(12, 12), sharex=True)

for ax, (band_name, regs) in zip(axes, bands.items()):
    # Sum across regions in band
    idxs = [COASTLINE_ORDER.index(r) for r in regs if r in COASTLINE_ORDER]
    pop_band = pop_reg[:, idxs].sum(axis=1)
    inf_band = inf_reg[:, idxs].sum(axis=1)
    
    color_pop = "#2166ac"
    color_inf = "#b2182b"
    
    ax.plot(years, pop_band, color=color_pop, lw=2, label="Population")
    ax.fill_between(years, 0, pop_band, color=color_pop, alpha=0.15)
    ax.set_ylabel("Population", color=color_pop)
    ax.tick_params(axis="y", labelcolor=color_pop)
    
    ax2 = ax.twinx()
    ax2.plot(years, inf_band, color=color_inf, lw=1.5, ls="--", label="Infected")
    ax2.fill_between(years, 0, inf_band, color=color_inf, alpha=0.1)
    ax2.set_ylabel("Infected", color=color_inf)
    ax2.tick_params(axis="y", labelcolor=color_inf)
    
    ax.set_title(band_name, fontsize=11, fontweight="bold")
    ax.axvline(2013.5, color="gray", ls=":", lw=0.8, alpha=0.6)
    ax.grid(True, alpha=0.2)

axes[-1].set_xlabel("Year")
fig.suptitle("W325: Population & Infection Cascades by Latitude Band", fontsize=13, fontweight="bold", y=0.98)
fig.tight_layout(rect=[0, 0, 1, 0.96])

fig.savefig(f"{FIGDIR}/fig_3_4_cascade_panels.pdf")
plt.close()
print("✓ fig_3_4_cascade_panels.pdf")
