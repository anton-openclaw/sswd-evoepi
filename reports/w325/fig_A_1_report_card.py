#!/usr/bin/env python3
"""Fig A.1: Report card — 2×2 composite summary figure."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from fig_common import *

plt = setup_style()
import matplotlib.cm as cm

result = load_result()
rd = result["region_details"]
scoring = result["scoring"]["per_region"]
recovery = result["region_recovery"]

data = load_npz()
populations = data["populations"]
infected = data["infected"]
sim_days = data["sim_days"]
site_names = data["site_names"]
regions_per_site = get_regions_per_site(site_names)
site_counts = get_region_site_counts(site_names)

years_ts = START_YEAR + sim_days / 365.25

pop_reg = aggregate_to_regions(populations, regions_per_site)
inf_reg = aggregate_to_regions(infected, regions_per_site)

fig, axes = plt.subplots(2, 2, figsize=(16, 10))

# (a) Infection heatmap
ax = axes[0, 0]
with np.errstate(divide='ignore', invalid='ignore'):
    prevalence = np.where(pop_reg > 0, inf_reg / pop_reg, 0.0)
extent = [years_ts[0], years_ts[-1], -0.5, len(COASTLINE_ORDER) - 0.5]
im = ax.imshow(prevalence.T, aspect="auto", origin="lower", extent=extent,
               cmap="Reds", vmin=0, vmax=0.5, interpolation="nearest")
ax.set_yticks(range(len(COASTLINE_ORDER)))
ax.set_yticklabels(COASTLINE_ORDER, fontsize=7)
ax.set_xlabel("Year", fontsize=9)
ax.set_title("(a) Infection Prevalence", fontsize=11, fontweight="bold")
ax.axvline(2013.5, color="black", ls="--", lw=0.8, alpha=0.6)
fig.colorbar(im, ax=ax, shrink=0.7, pad=0.02)

# (b) Population heatmap
ax = axes[0, 1]
pop_frac = np.zeros_like(pop_reg)
for j, reg in enumerate(COASTLINE_ORDER):
    n_sites = site_counts.get(reg, 1)
    pop_frac[:, j] = pop_reg[:, j] / (n_sites * K_PER_SITE)
im = ax.imshow(pop_frac.T, aspect="auto", origin="lower", extent=extent,
               cmap="YlGn", vmin=0, vmax=1.0, interpolation="nearest")
ax.set_yticks(range(len(COASTLINE_ORDER)))
ax.set_yticklabels(COASTLINE_ORDER, fontsize=7)
ax.set_xlabel("Year", fontsize=9)
ax.set_title("(b) Population / K", fontsize=11, fontweight="bold")
ax.axvline(2013.5, color="black", ls="--", lw=0.8, alpha=0.6)
fig.colorbar(im, ax=ax, shrink=0.7, pad=0.02)

# (c) Recovery bar chart
ax = axes[1, 0]
x = np.arange(len(COASTLINE_ORDER))
rec_vals = [recovery.get(r, 0) * 100 for r in COASTLINE_ORDER]
cmap_cw = cm.get_cmap("coolwarm")
colors = [cmap_cw(i / (len(COASTLINE_ORDER) - 1)) for i in range(len(COASTLINE_ORDER))]
bars = ax.bar(x, rec_vals, color=colors, edgecolor="gray", linewidth=0.3)
for i, reg in enumerate(COASTLINE_ORDER):
    if reg in SCORED_REGIONS:
        bars[i].set_edgecolor("black")
        bars[i].set_linewidth(1.2)
        if reg in TARGETS:
            ax.plot([i - 0.35, i + 0.35], [TARGETS[reg]*100]*2, color="red", lw=1.5, ls="--")
ax.set_xticks(x)
ax.set_xticklabels(COASTLINE_ORDER, rotation=45, ha="right", fontsize=7)
ax.set_ylabel("Recovery (%)", fontsize=9)
ax.set_title("(c) Recovery by Region", fontsize=11, fontweight="bold")
ax.grid(axis="y", alpha=0.3)

# (d) RMSLE breakdown
ax = axes[1, 1]
regions_scored = SCORED_COASTLINE
log_sq = [scoring[r]["log_sq_error"] for r in regions_scored]
cmap_err = cm.get_cmap("RdYlGn_r")
max_e = max(log_sq) if log_sq else 1
colors_err = [cmap_err(e / max_e) for e in log_sq]
y = np.arange(len(regions_scored))
ax.barh(y, log_sq, color=colors_err, edgecolor="gray", linewidth=0.5)
ax.set_yticks(y)
ax.set_yticklabels(regions_scored, fontsize=9)
ax.set_xlabel("(log error)²", fontsize=9)
ax.set_title(f"(d) RMSLE Breakdown (total = {result['scoring']['rmsle']:.3f})", fontsize=11, fontweight="bold")
ax.grid(axis="x", alpha=0.3)
for i, e in enumerate(log_sq):
    ax.text(e + max_e * 0.02, i, f"{e:.3f}", va="center", fontsize=8)

fig.suptitle("W325 Report Card — RMSLE = 0.348", fontsize=14, fontweight="bold", y=0.99)
fig.tight_layout(rect=[0, 0, 1, 0.97])

fig.savefig(f"{FIGDIR}/fig_A_1_report_card.pdf")
plt.close()
print("✓ fig_A_1_report_card.pdf")
