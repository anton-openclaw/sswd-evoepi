"""Fig 7: Disease dynamics panel — 2×2."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from load_data import *
import numpy as np

plt = setup_style()
d = load_sobol()
params = get_param_names(d)
short = get_short_names(params)
groups = get_param_groups(params)
sr = d['sobol_results']

dis_metrics = ['pop_crash_pct', 'peak_mortality', 'time_to_nadir', 'disease_death_fraction']
dis_titles = ['Population Crash (%)', 'Peak Mortality', 'Time to Nadir', 'Disease Death Fraction']

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

for ax, metric, title in zip(axes.flat, dis_metrics, dis_titles):
    ST = np.array(sr[metric]['ST'])
    ST_conf = np.array(sr[metric]['ST_conf'])
    top10 = np.argsort(ST)[-10:][::-1]

    y = np.arange(10)
    colors = [GROUP_COLORS.get(groups[i], '#999') for i in top10]
    ax.barh(y, ST[top10], xerr=ST_conf[top10], color=colors, edgecolor='white', linewidth=0.3,
            height=0.7, error_kw=dict(elinewidth=0.7, capsize=2, capthick=0.7, color='#333'))
    ax.set_yticks(y)
    ax.set_yticklabels([short[i] for i in top10], fontsize=8)
    ax.set_xlabel('$S_T$', fontsize=9)
    ax.set_title(title, fontsize=10, fontweight='bold')
    ax.invert_yaxis()
    ax.set_xlim(left=0)

from matplotlib.patches import Patch
handles = [Patch(facecolor=GROUP_COLORS[g], label=GROUP_LABELS[g]) for g in GROUP_ORDER]
fig.legend(handles=handles, loc='lower center', ncol=6, fontsize=8, framealpha=0.9,
           bbox_to_anchor=(0.5, -0.02))

plt.suptitle('Disease Dynamics: Top 10 Parameters by $S_T$', fontsize=13, fontweight='bold')
plt.tight_layout(rect=[0, 0.03, 1, 0.96])
plt.savefig(os.path.join(FIG_DIR, 'fig07_disease.png'))
plt.close()
print("Fig 7 done")
