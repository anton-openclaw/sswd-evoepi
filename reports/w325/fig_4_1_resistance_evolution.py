#!/usr/bin/env python3
"""Fig 4.1: Resistance evolution over time per region."""
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

focal = ["AK-PWS", "JDF", "CA-N"]
years = YEARS

for reg in COASTLINE_ORDER:
    if reg not in rd:
        continue
    d = rd[reg]
    res = d["yearly_mean_resistance"]
    if len(res) != len(years):
        continue
    lat = region_lats.get(reg, 45)
    color = cmap(norm(lat))
    is_focal = reg in focal
    lw = 2.5 if is_focal else 0.8
    alpha = 1.0 if is_focal else 0.4
    ax.plot(years, res, color=color, lw=lw, alpha=alpha, label=reg if is_focal else None)
    
    # ±1 SD band for focal regions
    if is_focal:
        va = d["yearly_va_resistance"]
        sd = [np.sqrt(v) if v > 0 else 0 for v in va]
        upper = [r + s for r, s in zip(res, sd)]
        lower = [r - s for r, s in zip(res, sd)]
        ax.fill_between(years, lower, upper, color=color, alpha=0.15)

ax.axvline(2013, color="gray", ls=":", lw=0.8, alpha=0.6)
ax.set_xlabel("Year")
ax.set_ylabel("Mean Resistance")
ax.set_title("W325: Resistance Evolution by Region", fontsize=13, fontweight="bold")
ax.legend(loc="upper left", fontsize=9)
ax.grid(True, alpha=0.3)

# Colorbar for latitude
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cb = fig.colorbar(sm, ax=ax, label="Latitude (°N)", shrink=0.7, pad=0.02)

fig.savefig(f"{FIGDIR}/fig_4_1_resistance_evolution.pdf")
plt.close()
print("✓ fig_4_1_resistance_evolution.pdf")
