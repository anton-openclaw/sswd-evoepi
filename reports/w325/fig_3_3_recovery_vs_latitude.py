#!/usr/bin/env python3
"""Fig 3.3: Recovery fraction vs latitude scatter."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from fig_common import *

plt = setup_style()
import matplotlib.cm as cm

result = load_result()
data = load_npz()
region_lats = get_region_lats(data["site_names"], data["site_lats"])
recovery = result["region_recovery"]

fig, ax = plt.subplots(figsize=(10, 6))

# Get latitude range for colormap normalization
all_lats = [region_lats[r] for r in COASTLINE_ORDER if r in region_lats]
norm = plt.Normalize(min(all_lats), max(all_lats))
cmap = cm.get_cmap("coolwarm")

# Plot all regions as scatter
for reg in COASTLINE_ORDER:
    if reg not in region_lats or reg not in recovery:
        continue
    lat = region_lats[reg]
    rec = recovery[reg]
    color = cmap(norm(lat))
    is_scored = reg in SCORED_REGIONS
    marker = "o" if not is_scored else "s"
    size = 100 if is_scored else 50
    edgecolor = "black" if is_scored else "gray"
    ax.scatter(lat, rec, c=[color], s=size, marker=marker, edgecolors=edgecolor, linewidths=1.2, zorder=3)
    if is_scored or rec > 0.15:
        offset = (0.3, 0.005) if rec > 0.01 else (0.3, 0.0005)
        ax.annotate(reg, (lat, rec), fontsize=7, ha="left", va="bottom",
                    xytext=(5, 3), textcoords="offset points")

# Plot targets as red diamonds
for reg, target in TARGETS.items():
    if reg in region_lats:
        lat = region_lats[reg]
        ax.scatter(lat, target, marker="D", c="red", s=60, edgecolors="darkred",
                   linewidths=0.8, zorder=4, alpha=0.7)

ax.set_yscale("log")
ax.set_xlabel("Latitude (°N)")
ax.set_ylabel("Recovery Fraction (log scale)")
ax.set_title("W325: Recovery Fraction vs Latitude", fontsize=13, fontweight="bold")
ax.set_ylim(1e-6, 1.0)

# Legend
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='s', color='w', markerfacecolor='gray', markersize=10, label='Scored region (actual)'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=8, label='Unscored region'),
    Line2D([0], [0], marker='D', color='w', markerfacecolor='red', markersize=8, label='Target'),
]
ax.legend(handles=legend_elements, loc="upper left", fontsize=9)
ax.grid(True, alpha=0.3)

fig.savefig(f"{FIGDIR}/fig_3_3_recovery_vs_latitude.pdf")
plt.close()
print("✓ fig_3_3_recovery_vs_latitude.pdf")
