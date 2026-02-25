"""Fig 3: Multi-metric heatmap — top 20 params × 23 metrics, ST values."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from load_data import *
import numpy as np

plt = setup_style()
d = load_sobol()
params = get_param_names(d)
short = get_short_names(params)
metrics = get_metric_names(d)
sr = d['sobol_results']

# Build ST matrix: params × metrics
n_p = len(params)
n_m = len(metrics)
ST_mat = np.zeros((n_p, n_m))
for j, m in enumerate(metrics):
    ST_mat[:, j] = np.array(sr[m]['ST'])

# Top 20 params by max ST across any metric
max_ST = ST_mat.max(axis=1)
top20_idx = np.argsort(max_ST)[-20:][::-1]

ST_sub = ST_mat[top20_idx, :]
short_sub = [short[i] for i in top20_idx]

# Simple metric short names
metric_short = [m.replace('_', '\n') if len(m) > 15 else m for m in metrics]

# Cluster columns (metrics) by similarity
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

if ST_sub.shape[1] > 2:
    col_linkage = linkage(pdist(ST_sub.T, 'euclidean'), method='ward')
    col_order = leaves_list(col_linkage)
else:
    col_order = np.arange(n_m)

if ST_sub.shape[0] > 2:
    row_linkage = linkage(pdist(ST_sub, 'euclidean'), method='ward')
    row_order = leaves_list(row_linkage)
else:
    row_order = np.arange(len(top20_idx))

ST_clustered = ST_sub[row_order][:, col_order]
short_clustered = [short_sub[i] for i in row_order]
metric_clustered = [metrics[i] for i in col_order]

fig, ax = plt.subplots(figsize=(14, 8))
# Clip negative to 0 for display
ST_disp = np.clip(ST_clustered, 0, None)
im = ax.imshow(ST_disp, aspect='auto', cmap='YlOrRd', interpolation='nearest')

ax.set_xticks(np.arange(n_m))
ax.set_xticklabels(metric_clustered, rotation=55, ha='right', fontsize=7.5)
ax.set_yticks(np.arange(len(short_clustered)))
ax.set_yticklabels(short_clustered, fontsize=8)

# Add text annotations for high values
for i in range(ST_disp.shape[0]):
    for j in range(ST_disp.shape[1]):
        v = ST_disp[i, j]
        if v > 0.05:
            ax.text(j, i, f'{v:.2f}', ha='center', va='center', fontsize=5.5,
                    color='white' if v > 0.15 else 'black', fontweight='bold')

cbar = plt.colorbar(im, ax=ax, shrink=0.7, label='$S_T$')
ax.set_title('Multi-Metric Sensitivity: Top 20 Parameters × 23 Metrics ($S_T$)',
             fontsize=13, fontweight='bold')

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'fig03_heatmap.png'))
plt.close()
print("Fig 3 done")
