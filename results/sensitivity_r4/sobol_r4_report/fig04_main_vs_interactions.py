"""Fig 4: Main effects vs interactions — stacked bar for top 15 params."""
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

# Top 15 by ST
top15 = np.argsort(ST)[-15:][::-1]
S1_top = np.clip(S1[top15], 0, None)  # clip negative S1 to 0
ST_top = ST[top15]
interaction = np.clip(ST_top - S1_top, 0, None)
short_top = [short[i] for i in top15]
groups_top = [groups[i] for i in top15]

fig, ax = plt.subplots(figsize=(10, 6))
y = np.arange(len(top15))

ax.barh(y, S1_top, color='#3498db', edgecolor='white', linewidth=0.3, label='$S_1$ (main effect)', height=0.7)
ax.barh(y, interaction, left=S1_top, color='#e67e22', edgecolor='white', linewidth=0.3, label='$S_T - S_1$ (interactions)', height=0.7)

ax.set_yticks(y)
ax.set_yticklabels(short_top, fontsize=9)
ax.set_xlabel('Sobol Index')
ax.set_title('Main Effects vs Interactions: Top 15 Parameters (pop_crash_pct)', fontsize=12, fontweight='bold')
ax.legend(loc='lower right', fontsize=9)
ax.invert_yaxis()
ax.set_xlim(left=0)

# Add percentage labels
for i in range(len(top15)):
    total = ST_top[i]
    if total > 0.01:
        pct = S1_top[i] / total * 100
        ax.text(total + 0.005, i, f'{pct:.0f}% main', fontsize=7, va='center', color='#555')

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'fig04_main_vs_interactions.png'))
plt.close()
print("Fig 4 done")
