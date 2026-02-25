"""Fig 5: Parameter group contributions — grouped bar for pop_crash_pct."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from load_data import *
import numpy as np

plt = setup_style()
d = load_sobol()
params = get_param_names(d)
groups = get_param_groups(params)
sr = d['sobol_results']['pop_crash_pct']
ST = np.array(sr['ST'])

# Sum ST by group
group_totals = {}
group_counts = {}
for g, st in zip(groups, ST):
    group_totals[g] = group_totals.get(g, 0) + max(st, 0)
    group_counts[g] = group_counts.get(g, 0) + 1

# Sort by total
sorted_groups = sorted(group_totals.keys(), key=lambda g: group_totals[g], reverse=True)

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Left: total ST by group
ax = axes[0]
vals = [group_totals[g] for g in sorted_groups]
colors = [GROUP_COLORS.get(g, '#999') for g in sorted_groups]
labels = [GROUP_LABELS.get(g, g) for g in sorted_groups]
bars = ax.bar(range(len(sorted_groups)), vals, color=colors, edgecolor='white', width=0.7)
ax.set_xticks(range(len(sorted_groups)))
ax.set_xticklabels(labels, rotation=30, ha='right', fontsize=9)
ax.set_ylabel('Sum of $S_T$')
ax.set_title('Total $S_T$ by Parameter Group', fontsize=11, fontweight='bold')
for bar, v, g in zip(bars, vals, sorted_groups):
    ax.text(bar.get_x() + bar.get_width()/2, v + 0.01, f'{v:.3f}\n({group_counts[g]} params)',
            ha='center', fontsize=8, color='#333')

# Right: mean ST by group
ax = axes[1]
means = [group_totals[g] / group_counts[g] for g in sorted_groups]
bars = ax.bar(range(len(sorted_groups)), means, color=colors, edgecolor='white', width=0.7)
ax.set_xticks(range(len(sorted_groups)))
ax.set_xticklabels(labels, rotation=30, ha='right', fontsize=9)
ax.set_ylabel('Mean $S_T$ per parameter')
ax.set_title('Mean $S_T$ by Parameter Group', fontsize=11, fontweight='bold')
for bar, v in zip(bars, means):
    ax.text(bar.get_x() + bar.get_width()/2, v + 0.002, f'{v:.4f}',
            ha='center', fontsize=8, color='#333')

plt.suptitle('Parameter Group Contributions to Population Crash', fontsize=13, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'fig05_group_contributions.png'))
plt.close()
print("Fig 5 done")
