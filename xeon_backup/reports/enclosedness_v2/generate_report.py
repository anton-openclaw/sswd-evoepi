#!/usr/bin/env python3
"""
Generate figures for the SSWD-EvoEpi site enclosedness analysis with fjord depth metric.
"""

import json
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import warnings
warnings.filterwarnings('ignore')

# Try to import cartopy and geopandas for mapping
try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    CARTOPY_AVAILABLE = True
except ImportError:
    CARTOPY_AVAILABLE = False
    print("Cartopy not available, using matplotlib for maps")

try:
    import geopandas as gpd
    GEOPANDAS_AVAILABLE = True
except ImportError:
    GEOPANDAS_AVAILABLE = False
    print("GeoPandas not available, mapping will be simplified")

# Set up paths
DATA_PATH = Path("../../data/nodes/site_enclosedness.json")
FIGURES_PATH = Path("figures")
FIGURES_PATH.mkdir(exist_ok=True)

# Load data
print("Loading site enclosedness data...")
with open(DATA_PATH, 'r') as f:
    sites_data = json.load(f)

# Convert to DataFrame
df = pd.DataFrame(sites_data)

print(f"Loaded {len(df)} sites across {df['region'].nunique()} regions")

# Set up plotting parameters
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.style.use('default')

# CENTRALIZED: moved to sswd_evoepi.results (canonical 18-region ordering)
# NOTE: This file uses a different, more granular region list than the canonical one.
# Keeping the local version since it includes sub-regions not in the canonical ordering.
from sswd_evoepi.results import REGION_ORDER as _CANONICAL_REGION_ORDER
REGION_ORDER = [
    "AK-FN", "AK-PWS", "AK-AL", "AK-BS", "AK-SE", 
    "BC-N", "BC-C", "BC-HG", "BC-S", 
    "WA-N", "WA-PS", "SS-S", "JDF", "WA-C", "WA-S", 
    "OR", "CA-N", "CA-S"
]

# Ensure all regions in data are in our ordering
missing_regions = set(df['region'].unique()) - set(REGION_ORDER)
if missing_regions:
    print(f"Warning: Found regions not in ordering: {missing_regions}")
    REGION_ORDER.extend(sorted(missing_regions))

df['region'] = pd.Categorical(df['region'], categories=REGION_ORDER, ordered=True)

def create_basemap(ax, extent=None):
    """Create a coastline basemap for the Pacific Northwest region."""
    if extent is None:
        # Default extent covering Alaska to Baja California
        extent = [-180, -110, 30, 70]
    
    if CARTOPY_AVAILABLE:
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_extent(extent, crs=ccrs.PlateCarree())
        ax.add_feature(cfeature.LAND, color='lightgray', alpha=0.7)
        ax.add_feature(cfeature.OCEAN, color='lightblue', alpha=0.3)
        ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
        ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
    else:
        # Fallback to simple matplotlib
        if GEOPANDAS_AVAILABLE:
            try:
                # Try to load the shapefile if available
                shapefile_path = Path("../../data/shorelines/ne_10m_land/ne_10m_land.shp")
                if shapefile_path.exists():
                    land = gpd.read_file(shapefile_path)
                    land.plot(ax=ax, color='lightgray', alpha=0.7)
            except Exception as e:
                print(f"Warning: Could not load shapefile: {e}")
        
        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2], extent[3])
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
        ax.grid(True, alpha=0.3)
    
    return ax

# Figure 1: Map of Fjord Depth
print("Creating Figure 1: Fjord Depth Map...")
fig, ax = plt.subplots(1, 1, figsize=(12, 8))
if CARTOPY_AVAILABLE:
    ax = plt.axes(projection=ccrs.PlateCarree())
ax = create_basemap(ax)

scatter = ax.scatter(df['lon'], df['lat'], c=df['fjord_depth_km'], 
                    cmap='viridis', s=8, alpha=0.8, 
                    transform=ccrs.PlateCarree() if CARTOPY_AVAILABLE else None)
cbar = plt.colorbar(scatter, ax=ax, shrink=0.6, aspect=20)
cbar.set_label('Fjord Depth (km)', rotation=270, labelpad=15)

plt.title('Site Fjord Depth: Over-water Distance to Open Coast\n'
          '(896 sites, SSWD-EvoEpi)', fontsize=14, fontweight='bold', pad=20)
plt.tight_layout()
plt.savefig(FIGURES_PATH / "01_fjord_depth_map.png", dpi=300, bbox_inches='tight')
plt.close()

# Figure 2: Map of Updated Flushing Rate
print("Creating Figure 2: Flushing Rate Map...")
fig, ax = plt.subplots(1, 1, figsize=(12, 8))
if CARTOPY_AVAILABLE:
    ax = plt.axes(projection=ccrs.PlateCarree())
ax = create_basemap(ax)

# Use RdYlBu_r colormap (blue=low/retained, red=high/flushed)
scatter = ax.scatter(df['lon'], df['lat'], c=df['flushing_rate'], 
                    cmap='RdYlBu_r', s=8, alpha=0.8, vmin=0, vmax=1,
                    transform=ccrs.PlateCarree() if CARTOPY_AVAILABLE else None)
cbar = plt.colorbar(scatter, ax=ax, shrink=0.6, aspect=20)
cbar.set_label('Flushing Rate (φ)', rotation=270, labelpad=15)

plt.title('Site Flushing Rate: Updated with Fjord Depth\n'
          '(Blue=retained, Red=flushed)', fontsize=14, fontweight='bold', pad=20)
plt.tight_layout()
plt.savefig(FIGURES_PATH / "02_flushing_rate_map.png", dpi=300, bbox_inches='tight')
plt.close()

# Figure 3: Map of Enclosedness (combined)
print("Creating Figure 3: Enclosedness Map...")
fig, ax = plt.subplots(1, 1, figsize=(12, 8))
if CARTOPY_AVAILABLE:
    ax = plt.axes(projection=ccrs.PlateCarree())
ax = create_basemap(ax)

scatter = ax.scatter(df['lon'], df['lat'], c=df['enclosedness_combined'], 
                    cmap='plasma', s=8, alpha=0.8, vmin=0, vmax=1,
                    transform=ccrs.PlateCarree() if CARTOPY_AVAILABLE else None)
cbar = plt.colorbar(scatter, ax=ax, shrink=0.6, aspect=20)
cbar.set_label('Enclosedness (combined)', rotation=270, labelpad=15)

plt.title('Site Enclosedness: Tortuosity + Ray-casting\n'
          '(Higher values = more enclosed)', fontsize=14, fontweight='bold', pad=20)
plt.tight_layout()
plt.savefig(FIGURES_PATH / "03_enclosedness_map.png", dpi=300, bbox_inches='tight')
plt.close()

# Figure 4: Scatter plot - Fjord Depth vs Enclosedness
print("Creating Figure 4: Fjord Depth vs Enclosedness...")
fig, ax = plt.subplots(1, 1, figsize=(10, 8))

# Create color map for regions
region_colors = plt.cm.tab20(np.linspace(0, 1, len(REGION_ORDER)))
region_color_map = dict(zip(REGION_ORDER, region_colors))

for region in REGION_ORDER:
    if region in df['region'].values:
        region_data = df[df['region'] == region]
        ax.scatter(region_data['enclosedness_combined'], region_data['fjord_depth_km'],
                  c=[region_color_map[region]], label=region, alpha=0.7, s=20)

ax.set_xlabel('Enclosedness (combined)', fontsize=12)
ax.set_ylabel('Fjord Depth (km)', fontsize=12)
ax.set_title('Fjord Depth vs Enclosedness by Region\n'
             'Deep fjords = high enclosedness + high depth', 
             fontsize=14, fontweight='bold')
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8, ncol=1)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(FIGURES_PATH / "04_depth_vs_enclosedness.png", dpi=300, bbox_inches='tight')
plt.close()

# Figure 5: Box plot - Flushing Rate by Region
print("Creating Figure 5: Flushing Rate by Region...")
fig, ax = plt.subplots(1, 1, figsize=(12, 8))

# Create box plot data
box_data = [df[df['region'] == region]['flushing_rate'].values 
           for region in REGION_ORDER if region in df['region'].values]
box_labels = [region for region in REGION_ORDER if region in df['region'].values]

bp = ax.boxplot(box_data, labels=box_labels, patch_artist=True)

# Color boxes with a gradient
colors = plt.cm.viridis(np.linspace(0, 1, len(bp['boxes'])))
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)

ax.set_ylabel('Flushing Rate (φ)', fontsize=12)
ax.set_xlabel('Region', fontsize=12)
ax.set_title('Flushing Rate Distribution by Region\n'
             '(ordered roughly North → South)', fontsize=14, fontweight='bold')
ax.tick_params(axis='x', rotation=45)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(FIGURES_PATH / "05_flushing_by_region.png", dpi=300, bbox_inches='tight')
plt.close()

# Figure 6: Box plot - Fjord Depth by Region
print("Creating Figure 6: Fjord Depth by Region...")
fig, ax = plt.subplots(1, 1, figsize=(12, 8))

# Create box plot data
box_data = [df[df['region'] == region]['fjord_depth_km'].values 
           for region in REGION_ORDER if region in df['region'].values]

bp = ax.boxplot(box_data, labels=box_labels, patch_artist=True)

# Color boxes with a gradient
colors = plt.cm.plasma(np.linspace(0, 1, len(bp['boxes'])))
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)

ax.set_ylabel('Fjord Depth (km)', fontsize=12)
ax.set_xlabel('Region', fontsize=12)
ax.set_title('Fjord Depth Distribution by Region\n'
             '(over-water distance to open coast)', fontsize=14, fontweight='bold')
ax.tick_params(axis='x', rotation=45)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(FIGURES_PATH / "06_depth_by_region.png", dpi=300, bbox_inches='tight')
plt.close()

# Figure 7: Before vs After fjord depth comparison
print("Creating Figure 7: Before vs After fjord depth...")
fig, ax = plt.subplots(1, 1, figsize=(10, 8))

# Calculate v1 flushing rate (without fjord depth)
# φ = 0.8*(1-encl) + 0.03*encl = 0.8 - 0.77*encl
df['flushing_rate_v1'] = 0.8 * (1 - df['enclosedness_combined']) + 0.03 * df['enclosedness_combined']

# Current flushing rate is v2 (with fjord depth)
scatter = ax.scatter(df['flushing_rate_v1'], df['flushing_rate'], 
                    c=df['fjord_depth_km'], cmap='viridis', alpha=0.6, s=15)

# Add 1:1 line
lims = [0, 1]
ax.plot(lims, lims, 'k--', alpha=0.75, zorder=0, label='1:1 line')

ax.set_xlabel('Flushing Rate v1 (without fjord depth)', fontsize=12)
ax.set_ylabel('Flushing Rate v2 (with fjord depth)', fontsize=12)
ax.set_title('Impact of Adding Fjord Depth Metric\n'
             'Points above line: depth increased flushing (fjord mouths)\n'
             'Points below line: depth decreased flushing', 
             fontsize=14, fontweight='bold')

cbar = plt.colorbar(scatter, ax=ax)
cbar.set_label('Fjord Depth (km)', rotation=270, labelpad=15)

ax.legend()
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)

plt.tight_layout()
plt.savefig(FIGURES_PATH / "07_before_after_depth.png", dpi=300, bbox_inches='tight')
plt.close()

# Figure 8: Mean Flushing Rate by Region (horizontal bar chart)
print("Creating Figure 8: Mean Flushing Rate by Region...")
fig, ax = plt.subplots(1, 1, figsize=(10, 8))

# Calculate regional statistics
regional_stats = df.groupby('region')['flushing_rate'].agg(['mean', 'std']).reset_index()
regional_stats = regional_stats.set_index('region').reindex(REGION_ORDER).dropna()

# Create horizontal bar plot
y_pos = np.arange(len(regional_stats))
bars = ax.barh(y_pos, regional_stats['mean'], 
               xerr=regional_stats['std'], capsize=3, alpha=0.7)

# Color bars with a gradient
colors = plt.cm.RdYlBu_r(np.linspace(0, 1, len(bars)))
for bar, color in zip(bars, colors):
    bar.set_color(color)

ax.set_yticks(y_pos)
ax.set_yticklabels(regional_stats.index)
ax.set_xlabel('Mean Flushing Rate (φ) ± SD', fontsize=12)
ax.set_ylabel('Region', fontsize=12)
ax.set_title('Regional Mean Flushing Rates\n'
             '(with fjord depth adjustment)', fontsize=14, fontweight='bold')
ax.grid(True, alpha=0.3, axis='x')

# Add values on bars
for i, (mean_val, std_val) in enumerate(zip(regional_stats['mean'], regional_stats['std'])):
    ax.text(mean_val + std_val + 0.01, i, f'{mean_val:.3f}', 
            va='center', ha='left', fontsize=9)

plt.tight_layout()
plt.savefig(FIGURES_PATH / "08_mean_flushing_by_region.png", dpi=300, bbox_inches='tight')
plt.close()

print("\n" + "="*50)
print("FIGURE GENERATION COMPLETE")
print("="*50)
print(f"Generated 8 figures in: {FIGURES_PATH.absolute()}")

# Print summary statistics for the report
print("\nSUMMARY STATISTICS:")
print(f"Total sites: {len(df)}")
print(f"Regions: {df['region'].nunique()}")
print(f"Flushing rate range: {df['flushing_rate'].min():.3f} - {df['flushing_rate'].max():.3f}")
print(f"Fjord depth range: {df['fjord_depth_km'].min():.1f} - {df['fjord_depth_km'].max():.1f} km")

print("\nKEY REGIONAL RESULTS:")
for region in ["AK-FN", "SS-S", "JDF", "AK-PWS", "OR", "CA-N", "CA-S"]:
    if region in df['region'].values:
        region_data = df[df['region'] == region]
        mean_depth = region_data['fjord_depth_km'].mean()
        mean_phi = region_data['flushing_rate'].mean()
        print(f"{region}: mean depth {mean_depth:.0f}km, φ={mean_phi:.3f}")

print("\nFigures ready for LaTeX compilation!")