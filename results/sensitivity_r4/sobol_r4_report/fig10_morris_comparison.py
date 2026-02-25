"""Fig 10: Morris vs Sobol rank comparison for pop_crash_pct."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from load_data import *
import numpy as np

plt = setup_style()
d = load_sobol()
params = get_param_names(d)
short = get_short_names(params)
groups = get_param_groups(params)
sr = d['sobol_results']['pop_crash_pct']

ST = np.array(sr['ST'])

# Sobol ST ranking (1 = highest)
sobol_order = np.argsort(ST)[::-1]
sobol_rank = np.zeros(len(params), dtype=int)
for rank, idx in enumerate(sobol_order):
    sobol_rank[idx] = rank + 1

# Morris ranks
morris_ranks = []
sobol_ranks = []
labels = []
grps = []
for i, p in enumerate(params):
    if p in MORRIS_RANKING:
        morris_ranks.append(MORRIS_RANKING[p])
        sobol_ranks.append(sobol_rank[i])
        labels.append(short[i])
        grps.append(groups[i])

morris_ranks = np.array(morris_ranks)
sobol_ranks = np.array(sobol_ranks)

fig, ax = plt.subplots(figsize=(8, 8))

# Diagonal = perfect agreement
ax.plot([0, 48], [0, 48], 'k--', alpha=0.3, linewidth=1, label='Perfect agreement')

# Color by group
for i in range(len(labels)):
    c = GROUP_COLORS.get(grps[i], '#999')
    ax.scatter(morris_ranks[i], sobol_ranks[i], c=c, s=50, edgecolors='white', linewidth=0.5, zorder=3)

# Highlight big rank changes (>10 positions)
for i in range(len(labels)):
    delta = abs(morris_ranks[i] - sobol_ranks[i])
    if delta >= 10 or sobol_ranks[i] <= 5 or morris_ranks[i] <= 5:
        ax.annotate(labels[i], (morris_ranks[i], sobol_ranks[i]),
                    textcoords='offset points', xytext=(5, 5),
                    fontsize=7, fontweight='bold' if delta >= 10 else 'normal',
                    color='red' if delta >= 10 else '#333',
                    arrowprops=dict(arrowstyle='-', color='#999', lw=0.5))

ax.set_xlabel('Morris $\\mu^*$ Rank')
ax.set_ylabel('Sobol $S_T$ Rank')
ax.set_title('Morris vs Sobol Ranking Comparison (pop_crash_pct)', fontsize=12, fontweight='bold')

# Shade disagreement zones
ax.fill_between([0, 48], [10, 58], [48, 48], alpha=0.03, color='blue', label='Sobol ranks higher')
ax.fill_between([0, 48], [0, 0], [-10, 38], alpha=0.03, color='red', label='Morris ranks higher')

ax.set_xlim(0, 48)
ax.set_ylim(0, 48)
ax.invert_yaxis()
ax.invert_xaxis()

# Make it: rank 1 at top-left
ax.invert_yaxis()
ax.invert_xaxis()

from matplotlib.patches import Patch
handles = [Patch(facecolor=GROUP_COLORS[g], label=GROUP_LABELS[g]) for g in GROUP_ORDER]
handles.insert(0, plt.Line2D([0], [0], linestyle='--', color='k', alpha=0.3, label='Perfect agreement'))
ax.legend(handles=handles, loc='lower right', fontsize=7, framealpha=0.9)

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'fig10_morris_comparison.png'))
plt.close()
print("Fig 10 done")
