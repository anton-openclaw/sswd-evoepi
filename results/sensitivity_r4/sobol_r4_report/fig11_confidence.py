"""Fig 11: Confidence intervals — dot plot for top 15 params (pop_crash_pct)."""
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
S1_conf = np.array(sr['S1_conf'])
ST = np.array(sr['ST'])
ST_conf = np.array(sr['ST_conf'])

top15 = np.argsort(ST)[-15:][::-1]

fig, ax = plt.subplots(figsize=(10, 7))
y = np.arange(len(top15))

# ST with confidence
ax.errorbar(ST[top15], y - 0.15, xerr=ST_conf[top15], fmt='o', color='#e74c3c',
            markersize=6, capsize=3, capthick=1, linewidth=1.2, label='$S_T \\pm$ conf', zorder=3)

# S1 with confidence
ax.errorbar(S1[top15], y + 0.15, xerr=S1_conf[top15], fmt='s', color='#3498db',
            markersize=5, capsize=3, capthick=1, linewidth=1.2, label='$S_1 \\pm$ conf', zorder=3)

ax.set_yticks(y)
ax.set_yticklabels([short[i] for i in top15], fontsize=9)
ax.axvline(x=0, color='#999', linewidth=0.5, linestyle='-')
ax.set_xlabel('Sobol Index Value')
ax.set_title('Confidence Intervals: Top 15 Parameters (pop_crash_pct)', fontsize=12, fontweight='bold')
ax.legend(loc='lower right', fontsize=9)
ax.invert_yaxis()

# Add a subtle grid
ax.grid(True, axis='x', alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'fig11_confidence.png'))
plt.close()
print("Fig 11 done")
