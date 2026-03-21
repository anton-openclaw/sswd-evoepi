#!/usr/bin/env python3
"""Fig 5.1: Grouped bar chart — target vs actual recovery for scored regions."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from fig_common import *

plt = setup_style()

result = load_result()
scoring = result["scoring"]["per_region"]

fig, ax = plt.subplots(figsize=(10, 6))

regions = SCORED_COASTLINE
x = np.arange(len(regions))
width = 0.35

targets = [TARGETS[r] * 100 for r in regions]
actuals = [scoring[r]["actual_pct"] for r in regions]

bars1 = ax.bar(x - width/2, targets, width, label="Target", color="#999999", edgecolor="gray")
bars2 = ax.bar(x + width/2, actuals, width, label="Actual (W325)", color="#2166ac", edgecolor="gray")

ax.set_yscale("log")
ax.set_xticks(x)
ax.set_xticklabels(regions, rotation=45, ha="right", fontsize=10)
ax.set_ylabel("Recovery (%)")
ax.set_title("W325: Target vs Actual Recovery (Scored Regions)", fontsize=13, fontweight="bold")
ax.legend(fontsize=10)
ax.grid(axis="y", alpha=0.3)

# Annotate within_2x
for i, reg in enumerate(regions):
    within = scoring[reg]["within_2x"]
    symbol = "✓" if within else "✗"
    color = "green" if within else "red"
    ypos = max(targets[i], actuals[i]) * 1.5
    ax.text(x[i], ypos, symbol, ha="center", fontsize=12, color=color, fontweight="bold")

ax.set_ylim(0.005, 200)
fig.tight_layout()

fig.savefig(f"{FIGDIR}/fig_5_1_recovery_vs_targets.pdf")
plt.close()
print("✓ fig_5_1_recovery_vs_targets.pdf")
