"""Fig 12: Calibration priority matrix — ST vs range width."""
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
range_widths = get_range_widths(params)

fig, ax = plt.subplots(figsize=(10, 8))

# Quadrant thresholds
st_thresh = 0.02
rw_thresh = np.median(range_widths)

# Background quadrant shading
ax.axvline(x=st_thresh, color='#999', linewidth=0.8, linestyle='--', alpha=0.5)
ax.axhline(y=rw_thresh, color='#999', linewidth=0.8, linestyle='--', alpha=0.5)

# Quadrant labels
mx_st = max(ST) * 1.1
mx_rw = max(range_widths) * 1.05
ax.text(mx_st * 0.55, mx_rw * 0.85, 'MUST CALIBRATE', fontsize=11, fontweight='bold',
        color='red', alpha=0.4, ha='center')
ax.text(mx_st * 0.55, rw_thresh * 0.3, 'WELL CONSTRAINED\n(verify range)', fontsize=9,
        color='#2ecc71', alpha=0.5, ha='center')
ax.text(st_thresh * 0.4, mx_rw * 0.85, 'WIDE BUT\nUNIMPORTANT', fontsize=9,
        color='#f39c12', alpha=0.5, ha='center')
ax.text(st_thresh * 0.4, rw_thresh * 0.3, 'FIX AT\nNOMINAL', fontsize=9,
        color='#3498db', alpha=0.5, ha='center')

colors = [GROUP_COLORS.get(g, '#999') for g in groups]
for i in range(len(params)):
    ax.scatter(ST[i], range_widths[i], c=colors[i], s=60, edgecolors='white', linewidth=0.5, zorder=3)

# Label important ones
for i in range(len(params)):
    if ST[i] > 0.03 or (ST[i] > 0.015 and range_widths[i] > rw_thresh):
        ax.annotate(short[i], (ST[i], range_widths[i]),
                    textcoords='offset points', xytext=(5, 5),
                    fontsize=7, fontweight='bold', color='#333',
                    arrowprops=dict(arrowstyle='-', color='#bbb', lw=0.5))

ax.set_xlabel('$S_T$ (Sensitivity Importance)', fontsize=11)
ax.set_ylabel('Normalized Range Width (range / midpoint)', fontsize=11)
ax.set_title('Calibration Priority Matrix', fontsize=13, fontweight='bold')

from matplotlib.patches import Patch
handles = [Patch(facecolor=GROUP_COLORS[g], label=GROUP_LABELS[g]) for g in GROUP_ORDER]
ax.legend(handles=handles, loc='upper left', fontsize=8, framealpha=0.9)

ax.set_xlim(left=-0.01)
ax.set_ylim(bottom=0)

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'fig12_calibration_priority.png'))
plt.close()
print("Fig 12 done")
