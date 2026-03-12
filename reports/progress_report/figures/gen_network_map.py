#!/usr/bin/env python3
"""Generate SSWD-EvoEpi network map with Salish Sea inset panel."""

import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from pathlib import Path

# --- Paths ---
BASE = Path("/home/starbot/.openclaw/workspace/sswd-evoepi")
SITES_JSON = BASE / "data/nodes/all_sites.json"
DIST_NPZ = BASE / "results/overwater/distance_matrix.npz"
OUT_PNG = BASE / "reports/progress_report/figures/fig_network_map.png"

# --- Load data ---
with open(SITES_JSON) as f:
    sites = json.load(f)

npz = np.load(DIST_NPZ, allow_pickle=True)
distances = npz["distances"]  # (896, 896)
npz_names = npz["names"]      # (896,)
npz_regions = npz["regions"]  # (896,)

# Build arrays from JSON (preserving order) - handle both lat/lon and latitude/longitude keys
names = np.array([s["name"] for s in sites])
lats = np.array([s.get("lat", s.get("latitude")) for s in sites], dtype=float)
lons = np.array([s.get("lon", s.get("longitude")) for s in sites], dtype=float)
regions = np.array([s["region"] for s in sites])

# --- Region order and colormap ---
REGION_ORDER = [
    "AK-AL", "AK-WG", "AK-PWS", "AK-EG", "AK-OC", "AK-FN", "AK-FS",
    "BC-N", "BC-C", "SS-N", "SS-S", "JDF", "WA-O", "OR",
    "CA-N", "CA-C", "CA-S", "BJ"
]

cmap = plt.cm.tab20
region_colors = {r: cmap(i / len(REGION_ORDER)) for i, r in enumerate(REGION_ORDER)}
site_colors = np.array([region_colors[r] for r in regions])

# --- Projection ---
proj = ccrs.LambertConformal(central_longitude=-135, central_latitude=45)

# --- Style ---
LAND_COLOR = "#F5E6CC"
OCEAN_COLOR = "#F0F4F8"
COAST_COLOR = "#888888"

# --- Create figure with manual axes positioning ---
fig = plt.figure(figsize=(14, 10), dpi=150)

# Main panel: left 68% width
ax_main = fig.add_axes([0.02, 0.05, 0.62, 0.88], projection=proj)
# Inset panel: right 30% width
ax_inset = fig.add_axes([0.66, 0.05, 0.33, 0.88], projection=proj)

# ==================== MAIN PANEL ====================
ax_main.set_facecolor(OCEAN_COLOR)
ax_main.add_feature(cfeature.LAND.with_scale("50m"), facecolor=LAND_COLOR, edgecolor="none")
ax_main.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=0.5, edgecolor=COAST_COLOR)
ax_main.add_feature(cfeature.BORDERS.with_scale("50m"), linewidth=0.3, edgecolor="#BBBBBB", linestyle="--")

# Set extent for full network (Alaska to Baja)
ax_main.set_extent([-180, -115, 24, 62], crs=ccrs.PlateCarree())

# Plot sites by region (so legend entries are grouped)
for region in REGION_ORDER:
    mask = regions == region
    ax_main.scatter(
        lons[mask], lats[mask],
        s=12, c=[region_colors[region]], label=region,
        transform=ccrs.PlateCarree(), zorder=5,
        edgecolors="black", linewidths=0.2, alpha=0.85
    )

# Draw red dashed rectangle for Salish Sea zoom area
zoom_lon = [-125.5, -122.0]
zoom_lat = [47.0, 50.5]

# Create rectangle in PlateCarree and transform
rect_lons = [zoom_lon[0], zoom_lon[1], zoom_lon[1], zoom_lon[0], zoom_lon[0]]
rect_lats = [zoom_lat[0], zoom_lat[0], zoom_lat[1], zoom_lat[1], zoom_lat[0]]
ax_main.plot(
    rect_lons, rect_lats,
    transform=ccrs.PlateCarree(),
    color="red", linewidth=1.5, linestyle="--", zorder=10
)

# Legend
legend = ax_main.legend(
    loc="lower left", ncol=2, fontsize=6.5,
    markerscale=1.5, framealpha=0.9, edgecolor="#CCCCCC",
    title="Regions (N→S)", title_fontsize=7.5,
    borderpad=0.6, columnspacing=0.8, handletextpad=0.3
)

ax_main.set_title("SSWD-EvoEpi: 896-Site Network", fontsize=13, fontweight="bold", pad=10)

# ==================== INSET PANEL ====================
ax_inset.set_facecolor(OCEAN_COLOR)
ax_inset.add_feature(cfeature.LAND.with_scale("50m"), facecolor=LAND_COLOR, edgecolor="none")
ax_inset.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=0.5, edgecolor=COAST_COLOR)

# Set extent for Salish Sea
ax_inset.set_extent([zoom_lon[0], zoom_lon[1], zoom_lat[0], zoom_lat[1]], crs=ccrs.PlateCarree())

# Find sites within the zoom bounds
in_zoom = (
    (lons >= zoom_lon[0]) & (lons <= zoom_lon[1]) &
    (lats >= zoom_lat[0]) & (lats <= zoom_lat[1])
)
zoom_indices = np.where(in_zoom)[0]
print(f"Sites in Salish Sea zoom: {len(zoom_indices)}")
print(f"  Regions present: {sorted(set(regions[in_zoom]))}")

# Draw connectivity lines for pairs within zoom where distance < 150 km
# Need to map site indices between JSON order and NPZ order
# Build name->npz_index mapping
npz_name_to_idx = {n: i for i, n in enumerate(npz_names)}

# Get NPZ indices for zoom sites
zoom_npz_indices = []
for si in zoom_indices:
    name = names[si]
    if name in npz_name_to_idx:
        zoom_npz_indices.append((si, npz_name_to_idx[name]))
    else:
        print(f"  Warning: {name} not found in distance matrix")

print(f"  Matched to distance matrix: {len(zoom_npz_indices)}")

# Draw lines
line_count = 0
for i in range(len(zoom_npz_indices)):
    si_i, ni_i = zoom_npz_indices[i]
    for j in range(i + 1, len(zoom_npz_indices)):
        si_j, ni_j = zoom_npz_indices[j]
        dist = distances[ni_i, ni_j]
        if 0 < dist < 150:
            ax_inset.plot(
                [lons[si_i], lons[si_j]], [lats[si_i], lats[si_j]],
                transform=ccrs.PlateCarree(),
                color="gray", alpha=0.15, linewidth=0.5, zorder=3
            )
            line_count += 1

print(f"  Connectivity lines drawn: {line_count}")

# Plot sites in zoom area (larger dots)
for region in REGION_ORDER:
    mask = in_zoom & (regions == region)
    if mask.any():
        ax_inset.scatter(
            lons[mask], lats[mask],
            s=30, c=[region_colors[region]], label=region,
            transform=ccrs.PlateCarree(), zorder=5,
            edgecolors="black", linewidths=0.3, alpha=0.9
        )

ax_inset.set_title("Salish Sea Detail", fontsize=11, fontweight="bold", pad=10)

# Add note about connectivity threshold
ax_inset.text(
    0.5, 0.02, "Lines show connections < 150 km (overwater)",
    transform=ax_inset.transAxes, fontsize=7,
    ha="center", va="bottom",
    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#CCCCCC", alpha=0.85)
)

# --- Save ---
OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT_PNG, dpi=150, bbox_inches="tight", facecolor="white")
plt.close(fig)

file_size = OUT_PNG.stat().st_size
print(f"\nSaved: {OUT_PNG}")
print(f"File size: {file_size:,} bytes ({file_size / 1024 / 1024:.1f} MB)")
