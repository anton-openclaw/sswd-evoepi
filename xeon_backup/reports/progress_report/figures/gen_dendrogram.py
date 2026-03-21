#!/usr/bin/env python3
"""Generate dendrogram + resistance divergence figure for SSWD-EvoEpi progress report."""

import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import pdist

# --- Load data ---
with open('/home/starbot/.openclaw/workspace/sswd-evoepi/results/calibration/W142/combined_results.json') as f:
    data = json.load(f)

rd = data['results'][0]['region_details']

# Regions in latitude order (north→south)
REGIONS = ['AK-AL', 'AK-WG', 'AK-PWS', 'AK-EG', 'AK-OC', 'AK-FN', 'AK-FS',
           'BC-N', 'BC-C', 'SS-N', 'SS-S', 'JDF', 'WA-O', 'OR', 'CA-N', 'CA-C', 'CA-S', 'BJ']

# Approximate latitudes for colormap
LATITUDES = {
    'AK-AL': 58.0, 'AK-WG': 57.5, 'AK-PWS': 60.5, 'AK-EG': 59.5,
    'AK-OC': 57.0, 'AK-FN': 56.5, 'AK-FS': 55.5,
    'BC-N': 53.0, 'BC-C': 50.5,
    'SS-N': 49.5, 'SS-S': 49.0, 'JDF': 48.3, 'WA-O': 47.5,
    'OR': 44.5, 'CA-N': 41.0, 'CA-C': 37.5, 'CA-S': 34.0, 'BJ': 30.0
}

years = list(range(2012, 2025))  # 13 years

# --- Build trait matrix for dendrogram ---
trait_names = ['resistance', 'tolerance', 'recovery', 'T_vbnc', 'v_local']
trait_matrix = []
for reg in REGIONS:
    r = rd[reg]
    vec = [
        r['yearly_mean_resistance'][-1],   # final resistance
        r['yearly_mean_tolerance'][-1],     # final tolerance
        r['yearly_mean_recovery'][-1],      # final recovery
        r['final_mean_T_vbnc'],             # pathogen thermal threshold
        r['final_mean_v_local'],            # pathogen virulence
    ]
    trait_matrix.append(vec)

trait_matrix = np.array(trait_matrix)

# Standardize (z-score)
means = trait_matrix.mean(axis=0)
stds = trait_matrix.std(axis=0)
stds[stds == 0] = 1  # avoid division by zero
trait_z = (trait_matrix - means) / stds

# --- Clustering ---
Z = linkage(trait_z, method='ward')

# Choose threshold for ~3-4 clusters
max_d = 0.7 * Z[-1, 2]  # 70% of max distance as starting point
clusters = fcluster(Z, t=max_d, criterion='distance')
n_clusters = len(set(clusters))
# Adjust if needed
if n_clusters < 3:
    max_d = 0.5 * Z[-1, 2]
elif n_clusters > 5:
    max_d = 0.8 * Z[-1, 2]

# --- Set up figure ---
plt.style.use('seaborn-v0_8-whitegrid')
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7), dpi=150,
                                 gridspec_kw={'width_ratios': [1, 1.1]})

# --- Left panel: Dendrogram ---
# Color palette for clusters
cluster_colors = ['#2166ac', '#d6604d', '#4daf4a', '#984ea3', '#ff7f00']

from scipy.cluster.hierarchy import set_link_color_palette
set_link_color_palette(cluster_colors)

dn = dendrogram(
    Z,
    labels=REGIONS,
    leaf_rotation=45,
    leaf_font_size=9,
    ax=ax1,
    color_threshold=max_d,
    above_threshold_color='#666666',
)

ax1.set_title('Regional Clustering by Evolved Traits', fontsize=13, fontweight='bold', pad=12)
ax1.set_ylabel('Ward distance (standardized traits)', fontsize=10)
ax1.set_xlabel('')
ax1.tick_params(axis='x', labelsize=9)
ax1.axhline(y=max_d, color='gray', linestyle='--', linewidth=0.8, alpha=0.6)
ax1.text(ax1.get_xlim()[1] * 0.98, max_d * 1.03, f'cluster threshold',
         ha='right', va='bottom', fontsize=8, color='gray', style='italic')

# Remove grid on x
ax1.grid(axis='x', visible=False)

# --- Right panel: Resistance divergence over time ---
lat_vals = np.array([LATITUDES[r] for r in REGIONS])
norm = Normalize(vmin=lat_vals.min(), vmax=lat_vals.max())
cmap = cm.coolwarm_r  # blue=high lat (north), red=low lat (south)

for reg in REGIONS:
    r = rd[reg]
    resistance = r['yearly_mean_resistance']
    lat = LATITUDES[reg]
    color = cmap(norm(lat))
    ax2.plot(years, resistance, color=color, linewidth=1.3, alpha=0.85)

# Add region labels at the end of each line
for reg in REGIONS:
    r = rd[reg]
    resistance = r['yearly_mean_resistance']
    lat = LATITUDES[reg]
    color = cmap(norm(lat))
    ax2.annotate(reg, xy=(2024, resistance[-1]), xytext=(2024.15, resistance[-1]),
                 fontsize=5.5, color=color, va='center', fontweight='bold',
                 annotation_clip=False)

ax2.set_title('Resistance Divergence from Common Starting Point',
              fontsize=13, fontweight='bold', pad=12)
ax2.set_xlabel('Year', fontsize=10)
ax2.set_ylabel('Mean host resistance ($r$)', fontsize=10)
ax2.set_xlim(2012, 2025.5)

# Colorbar
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax2, shrink=0.8, pad=0.02, aspect=25)
cbar.set_label('Latitude (°N)', fontsize=10)
cbar.ax.tick_params(labelsize=8)

# Add subtle annotation
ax2.annotate('common\nancestral\nvalue', xy=(2012, resistance[0]),
             xytext=(2012.5, 0.13),
             fontsize=7, color='gray', style='italic',
             arrowprops=dict(arrowstyle='->', color='gray', lw=0.8),
             ha='left')

# --- Final touches ---
fig.tight_layout(w_pad=3)

outpath = '/home/starbot/.openclaw/workspace/sswd-evoepi/reports/progress_report/figures/fig_evolution_dendrogram.png'
fig.savefig(outpath, dpi=150, bbox_inches='tight', facecolor='white')
plt.close()

import os
size_kb = os.path.getsize(outpath) / 1024
print(f'Saved: {outpath}')
print(f'Size: {size_kb:.0f} KB')
