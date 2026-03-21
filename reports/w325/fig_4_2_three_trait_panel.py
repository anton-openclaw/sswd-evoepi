#!/usr/bin/env python3
"""Fig 4.2: Three-trait panel — resistance, tolerance, recovery over time."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from fig_common import *

plt = setup_style()

result = load_result()
rd = result["region_details"]

# Representative regions spanning the coast
reps = ["AK-PWS", "AK-FS", "JDF", "CA-N"]
colors = {"AK-PWS": "#1f77b4", "AK-FS": "#ff7f0e", "JDF": "#2ca02c", "CA-N": "#d62728"}
years = YEARS

traits = [
    ("yearly_mean_resistance", "Mean Resistance", "Resistance"),
    ("yearly_mean_tolerance", "Mean Tolerance", "Tolerance"),
    ("yearly_mean_recovery", "Mean Recovery Rate", "Recovery"),
]

fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)

for ax, (key, ylabel, title) in zip(axes, traits):
    for reg in reps:
        if reg in rd:
            vals = rd[reg][key]
            if len(vals) == len(years):
                ax.plot(years, vals, color=colors[reg], lw=2, label=reg)
    ax.axvline(2013, color="gray", ls=":", lw=0.8, alpha=0.6)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=8, ncol=2)

axes[-1].set_xlabel("Year")
fig.suptitle("W325: Three-Trait Evolutionary Response", fontsize=13, fontweight="bold", y=0.98)
fig.tight_layout(rect=[0, 0, 1, 0.96])

fig.savefig(f"{FIGDIR}/fig_4_2_three_trait_panel.pdf")
plt.close()
print("✓ fig_4_2_three_trait_panel.pdf")
