#!/usr/bin/env python3
"""Generate map-based figures for SSWD-EvoEpi progress report."""

import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────
BASE = Path("/home/starbot/.openclaw/workspace/sswd-evoepi")
SITES_JSON = BASE / "data/nodes/all_sites.json"
NPZ_FILE = BASE / "results/calibration/W142/monthly_seed42.npz"
OUT_DIR = BASE / "reports/progress_report/figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Load data ──────────────────────────────────────────────────────────
with open(SITES_JSON) as f:
    sites = json.load(f)

data = np.load(NPZ_FILE, allow_pickle=True)
populations = data["populations"]   # (159, 896)
infected = data["infected"]         # (159, 896)
site_lats = data["site_lats"]       # (896,)
site_lons = data["site_lons"]       # (896,)
K = int(data["K"])

# Region list (canonical order for legend)
REGIONS = [
    "AK-AL", "AK-WG", "AK-PWS", "AK-EG", "AK-OC", "AK-FN", "AK-FS",
    "BC-N", "BC-C", "SS-N", "SS-S", "JDF", "WA-O", "OR",
    "CA-N", "CA-C", "CA-S", "BJ",
]
region_to_idx = {r: i for i, r in enumerate(REGIONS)}
site_regions = [site["region"] for site in sites]
site_region_idx = np.array([region_to_idx[r] for r in site_regions])

# ── Projection ─────────────────────────────────────────────────────────
proj = ccrs.LambertConformal(central_longitude=-135, central_latitude=45)
data_crs = ccrs.PlateCarree()

# ── Style constants ────────────────────────────────────────────────────
LAND_COLOR = "#F5E6CC"      # tan
OCEAN_COLOR = "#F0F4F8"     # light blue-white
COAST_COLOR = "#555555"

def setup_ax(ax):
    """Common axis styling."""
    ax.set_facecolor(OCEAN_COLOR)
    ax.add_feature(cfeature.LAND.with_scale("50m"), facecolor=LAND_COLOR, edgecolor="none")
    ax.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=0.5, edgecolor=COAST_COLOR)
    ax.set_extent([-180, -110, 20, 62], crs=data_crs)


# ══════════════════════════════════════════════════════════════════════
# Figure 1: Network Map
# ══════════════════════════════════════════════════════════════════════
print("Generating fig_network_map.png ...")

cmap20 = plt.cm.get_cmap("tab20", 20)
region_colors = [cmap20(i) for i in range(len(REGIONS))]

fig, ax = plt.subplots(figsize=(9, 12), subplot_kw={"projection": proj})
setup_ax(ax)

for i, region in enumerate(REGIONS):
    mask = site_region_idx == i
    if not mask.any():
        continue
    ax.scatter(
        site_lons[mask], site_lats[mask],
        s=12, color=region_colors[i], edgecolors="k", linewidths=0.2,
        label=region, transform=data_crs, zorder=5,
    )

ax.legend(
    loc="lower left", fontsize=7, ncol=2,
    framealpha=0.9, edgecolor="#ccc",
    title="Region", title_fontsize=8,
    markerscale=1.5,
)
ax.set_title("SSWD-EvoEpi: 896-Site Network", fontsize=14, fontweight="bold", pad=12)

fig.savefig(OUT_DIR / "fig_network_map.png", dpi=150, bbox_inches="tight",
            facecolor="white", edgecolor="none")
plt.close(fig)
print("  ✓ fig_network_map.png saved")


# ══════════════════════════════════════════════════════════════════════
# Figure 2: Disease Snapshots
# ══════════════════════════════════════════════════════════════════════
print("Generating fig_disease_snapshots.png ...")

# Colormap: green → yellow → red
disease_cmap = mcolors.LinearSegmentedColormap.from_list(
    "disease", ["#22AA22", "#FFDD00", "#DD2222"], N=256
)

# Timepoints
frames = [
    (11,  "Pre-disease (Dec 2012)"),
    (35,  "Peak crash (Dec 2014)"),
    (59,  "La Niña recovery (Dec 2016)"),
    (158, "Final state (Mar 2025)"),
]

fig, axes = plt.subplots(2, 2, figsize=(14, 12),
                         subplot_kw={"projection": proj})
axes = axes.flatten()

MAX_MARKER = 100  # pts²

for ax, (frame_idx, title) in zip(axes, frames):
    setup_ax(ax)
    ax.set_title(title, fontsize=12, fontweight="bold", pad=8)

    pop = populations[frame_idx].astype(float)
    inf = infected[frame_idx].astype(float)

    # Infected fraction
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = np.where(pop > 0, inf / pop, 0.0)
    frac = np.clip(frac, 0, 1)

    # Marker size proportional to pop/K
    sizes = (pop / K) * MAX_MARKER
    sizes = np.clip(sizes, 0, MAX_MARKER)

    # Alpha: extinct sites fully transparent
    alphas = np.where(pop > 0, 0.85, 0.0)

    # Sort so larger/redder dots render on top
    order = np.argsort(frac)
    sc = ax.scatter(
        site_lons[order], site_lats[order],
        s=sizes[order],
        c=frac[order],
        cmap=disease_cmap, vmin=0, vmax=1,
        alpha=alphas[order],
        edgecolors="k", linewidths=0.15,
        transform=data_crs, zorder=5,
    )

# Shared colorbar
cbar_ax = fig.add_axes([0.25, 0.04, 0.50, 0.018])
cbar = fig.colorbar(
    plt.cm.ScalarMappable(norm=mcolors.Normalize(0, 1), cmap=disease_cmap),
    cax=cbar_ax, orientation="horizontal",
)
cbar.set_label("Infected fraction", fontsize=11)
cbar.ax.tick_params(labelsize=9)

fig.suptitle("Disease Spread Over Time", fontsize=15, fontweight="bold", y=0.97)
fig.subplots_adjust(hspace=0.08, wspace=0.05, top=0.93, bottom=0.08)

fig.savefig(OUT_DIR / "fig_disease_snapshots.png", dpi=150, bbox_inches="tight",
            facecolor="white", edgecolor="none")
plt.close(fig)
print("  ✓ fig_disease_snapshots.png saved")

# ── Report ─────────────────────────────────────────────────────────────
print("\nFiles created:")
for f in sorted(OUT_DIR.glob("fig_*.png")):
    size_kb = f.stat().st_size / 1024
    print(f"  {f.name}: {size_kb:.0f} KB")
