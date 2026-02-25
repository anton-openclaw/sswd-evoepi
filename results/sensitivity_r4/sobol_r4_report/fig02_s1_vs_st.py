"""Fig 2: S1 vs ST scatter for pop_crash_pct."""
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

S1 = np.array(sr['S1'])
ST = np.array(sr['ST'])
colors = [GROUP_COLORS.get(g, '#999') for g in groups]

fig, ax = plt.subplots(figsize=(8, 7))

# Diagonal line
mx = max(ST.max(), 0.5)
ax.plot([0, mx], [0, mx], 'k--', alpha=0.3, linewidth=1, label='$S_1 = S_T$ (no interactions)')

for i in range(len(params)):
    ax.scatter(S1[i], ST[i], c=colors[i], s=50, edgecolors='white', linewidth=0.5, zorder=3)

# Label top 10 by ST
top10 = np.argsort(ST)[-10:]
from matplotlib.offsetbox import AnnotationBbox
for i in top10:
    ax.annotate(short[i], (S1[i], ST[i]),
                textcoords='offset points', xytext=(6, 4),
                fontsize=7.5, fontweight='bold', color='#333',
                arrowprops=dict(arrowstyle='-', color='#999', lw=0.5))

ax.set_xlabel('First-order index ($S_1$)')
ax.set_ylabel('Total-order index ($S_T$)')
ax.set_title('$S_1$ vs $S_T$: Interaction Detection (pop_crash_pct)', fontsize=12, fontweight='bold')

# Shaded region above diagonal = interactions
ax.fill_between([0, mx], [0, mx], [mx, mx], alpha=0.05, color='red', label='Interaction region')

from matplotlib.patches import Patch
handles = [Patch(facecolor=GROUP_COLORS[g], label=GROUP_LABELS[g]) for g in GROUP_ORDER]
handles.insert(0, plt.Line2D([0], [0], linestyle='--', color='k', alpha=0.3, label='No interactions'))
ax.legend(handles=handles, loc='upper left', fontsize=8, framealpha=0.9)

ax.set_xlim(left=min(S1.min()-0.02, -0.02))
ax.set_ylim(bottom=0)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'fig02_s1_vs_st.png'))
plt.close()
print("Fig 2 done")
