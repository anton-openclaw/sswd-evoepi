#!/usr/bin/env python3
"""Fig 5.2: Horizontal bar chart — per-region contribution to RMSLE²."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from fig_common import *

plt = setup_style()

result = load_result()
scoring = result["scoring"]["per_region"]

fig, ax = plt.subplots(figsize=(9, 5))

regions = SCORED_COASTLINE
log_sq_errors = [scoring[r]["log_sq_error"] for r in regions]

# Color: higher error = more red
cmap = plt.cm.get_cmap("RdYlGn_r")
max_err = max(log_sq_errors)
colors = [cmap(e / max_err) for e in log_sq_errors]

y = np.arange(len(regions))
bars = ax.barh(y, log_sq_errors, color=colors, edgecolor="gray", linewidth=0.5)

ax.set_yticks(y)
ax.set_yticklabels(regions, fontsize=10)
ax.set_xlabel("(log error)²")
ax.set_title(f"W325: Per-Region Contribution to RMSLE² (total RMSLE = {result['scoring']['rmsle']:.3f})",
             fontsize=12, fontweight="bold")
ax.grid(axis="x", alpha=0.3)

# Annotate values
for i, (reg, err) in enumerate(zip(regions, log_sq_errors)):
    ax.text(err + max_err * 0.02, i, f"{err:.4f}", va="center", fontsize=9)

fig.tight_layout()

fig.savefig(f"{FIGDIR}/fig_5_2_scoring_breakdown.pdf")
plt.close()
print("✓ fig_5_2_scoring_breakdown.pdf")
