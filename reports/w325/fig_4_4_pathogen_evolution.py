#!/usr/bin/env python3
"""Fig 4.4: Pathogen cold-adaptation — final_mean_T_vbnc by region."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from fig_common import *

plt = setup_style()
import matplotlib.cm as cm

result = load_result()
rd = result["region_details"]

fig, ax = plt.subplots(figsize=(12, 5))

x = np.arange(len(COASTLINE_ORDER))
t_vbnc = []
for reg in COASTLINE_ORDER:
    if reg in rd:
        t_vbnc.append(rd[reg]["final_mean_T_vbnc"])
    else:
        t_vbnc.append(np.nan)

# Color by magma
cmap = cm.get_cmap("magma")
norm = plt.Normalize(min(t for t in t_vbnc if not np.isnan(t)),
                     max(t for t in t_vbnc if not np.isnan(t)))
colors = [cmap(norm(v)) if not np.isnan(v) else "gray" for v in t_vbnc]

bars = ax.bar(x, t_vbnc, color=colors, edgecolor="gray", linewidth=0.5)

# Mark initial T_vbnc_min = 10°C
ax.axhline(10.0, color="red", ls="--", lw=1.2, alpha=0.7, label="T_vbnc_min = 10°C")

ax.set_xticks(x)
ax.set_xticklabels(COASTLINE_ORDER, rotation=45, ha="right", fontsize=9)
ax.set_ylabel("Final Mean T_vbnc (°C)")
ax.set_title("W325: Pathogen Cold-Adaptation by Region", fontsize=13, fontweight="bold")
ax.legend(loc="upper right", fontsize=9)
ax.grid(axis="y", alpha=0.3)

# Add colorbar
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
fig.colorbar(sm, ax=ax, label="T_vbnc (°C)", shrink=0.7, pad=0.02)

fig.savefig(f"{FIGDIR}/fig_4_4_pathogen_evolution.pdf")
plt.close()
print("✓ fig_4_4_pathogen_evolution.pdf")
