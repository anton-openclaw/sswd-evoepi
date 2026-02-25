"""Fig 1: Global Parameter Ranking — ST for pop_crash_pct, all 47 params."""
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
ST_conf = np.array(sr['ST_conf'])

# Sort by ST
idx = np.argsort(ST)
ST_s = ST[idx]
ST_conf_s = ST_conf[idx]
short_s = [short[i] for i in idx]
groups_s = [groups[i] for i in idx]
colors = [GROUP_COLORS.get(g, '#999') for g in groups_s]

fig, ax = plt.subplots(figsize=(8, 12))
y = np.arange(len(ST_s))
ax.barh(y, ST_s, xerr=ST_conf_s, color=colors, edgecolor='white', linewidth=0.3, height=0.75,
        error_kw=dict(elinewidth=0.8, capsize=2, capthick=0.8, color='#333'))
ax.set_yticks(y)
ax.set_yticklabels(short_s, fontsize=7.5)
ax.set_xlabel('Total-order Sobol index ($S_T$)')
ax.set_title('Global Parameter Ranking: $S_T$ for Population Crash (%)', fontsize=13, fontweight='bold')

# Legend
from matplotlib.patches import Patch
handles = [Patch(facecolor=GROUP_COLORS[g], label=GROUP_LABELS[g]) for g in GROUP_ORDER]
ax.legend(handles=handles, loc='lower right', fontsize=8, framealpha=0.9)

ax.set_xlim(left=0)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'fig01_global_ranking.png'))
plt.close()
print("Fig 1 done")
