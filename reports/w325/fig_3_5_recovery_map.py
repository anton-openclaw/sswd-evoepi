#!/usr/bin/env python3
"""Fig 3.5: Recovery % bar chart by region (S→N coastline order)."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from fig_common import *

plt = setup_style()
import matplotlib.cm as cm

result = load_result()
recovery = result["region_recovery"]
data = load_npz()
region_lats = get_region_lats(data["site_names"], data["site_lats"])

fig, ax = plt.subplots(figsize=(12, 5))

x = np.arange(len(COASTLINE_ORDER))
rec_vals = [recovery.get(r, 0) * 100 for r in COASTLINE_ORDER]

# Color by coolwarm based on position (S→N)
cmap = cm.get_cmap("coolwarm")
colors = [cmap(i / (len(COASTLINE_ORDER) - 1)) for i in range(len(COASTLINE_ORDER))]

bars = ax.bar(x, rec_vals, color=colors, edgecolor="gray", linewidth=0.5)

# Highlight scored regions with thicker edges
for i, reg in enumerate(COASTLINE_ORDER):
    if reg in SCORED_REGIONS:
        bars[i].set_edgecolor("black")
        bars[i].set_linewidth(1.5)

# Target lines for scored regions
for reg, target in TARGETS.items():
    if reg in COASTLINE_ORDER:
        idx = COASTLINE_ORDER.index(reg)
        ax.plot([idx - 0.4, idx + 0.4], [target * 100, target * 100],
                color="red", lw=2, ls="--", zorder=5)

ax.set_xticks(x)
ax.set_xticklabels(COASTLINE_ORDER, rotation=45, ha="right", fontsize=9)
ax.set_ylabel("Recovery (%)")
ax.set_title("W325: Recovery Fraction by Region (S→N)", fontsize=13, fontweight="bold")
ax.set_ylim(0, max(rec_vals) * 1.15)
ax.grid(axis="y", alpha=0.3)

# Legend
from matplotlib.lines import Line2D
legend = [
    Line2D([0], [0], color="red", ls="--", lw=2, label="Calibration target"),
]
ax.legend(handles=legend, loc="upper left", fontsize=9)

fig.savefig(f"{FIGDIR}/fig_3_5_recovery_map.pdf")
plt.close()
print("✓ fig_3_5_recovery_map.pdf")
