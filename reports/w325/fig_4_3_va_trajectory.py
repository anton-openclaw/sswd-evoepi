#!/usr/bin/env python3
"""Fig 4.3: Additive variance (Va) for resistance over time."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from fig_common import *

plt = setup_style()
import matplotlib.cm as cm

result = load_result()
rd = result["region_details"]
data = load_npz()
region_lats = get_region_lats(data["site_names"], data["site_lats"])

fig, ax = plt.subplots(figsize=(11, 6))

all_lats = [region_lats[r] for r in COASTLINE_ORDER if r in region_lats]
norm = plt.Normalize(min(all_lats), max(all_lats))
cmap = cm.get_cmap("coolwarm")
years = YEARS

for reg in COASTLINE_ORDER:
    if reg not in rd:
        continue
    va = rd[reg]["yearly_va_resistance"]
    if len(va) != len(years):
        continue
    lat = region_lats.get(reg, 45)
    color = cmap(norm(lat))
    # Thicker for regions with big Va crashes
    va_drop = (va[0] - min(va)) / va[0] if va[0] > 0 else 0
    lw = 2.0 if va_drop > 0.3 else 0.8
    alpha = 1.0 if va_drop > 0.3 else 0.4
    ax.plot(years, va, color=color, lw=lw, alpha=alpha)
    if va_drop > 0.3:
        ax.annotate(reg, (years[-1], va[-1]), fontsize=7, ha="left",
                    xytext=(3, 0), textcoords="offset points")

ax.axvline(2013, color="gray", ls=":", lw=0.8, alpha=0.6)
ax.set_xlabel("Year")
ax.set_ylabel("Additive Variance (Va) — Resistance")
ax.set_title("W325: Resistance Va Trajectories — Selective Sweep Signatures", fontsize=13, fontweight="bold")
ax.grid(True, alpha=0.3)

sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cb = fig.colorbar(sm, ax=ax, label="Latitude (°N)", shrink=0.7, pad=0.02)

fig.savefig(f"{FIGDIR}/fig_4_3_va_trajectory.pdf")
plt.close()
print("✓ fig_4_3_va_trajectory.pdf")
