#!/usr/bin/env python3
"""Generate progress report figures for SSWD-EvoEpi."""

import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# Load data
with open('/home/starbot/.openclaw/workspace/sswd-evoepi/results/calibration/W142/combined_results.json') as f:
    data = json.load(f)

rd = data['results'][0]['region_details']

# Config
YEARS = list(range(2012, 2025))  # 13 years
REGIONS = ['AK-PWS', 'AK-FN', 'BC-N', 'SS-S', 'CA-N', 'CA-C']
COLORS = {
    'AK-PWS': '#1b9e77',
    'AK-FN':  '#d95f02',
    'BC-N':   '#7570b3',
    'SS-S':   '#e7298a',
    'CA-N':   '#66a61e',
    'CA-C':   '#e6ab02',
}
OUTDIR = '/home/starbot/.openclaw/workspace/sswd-evoepi/reports/progress_report/figures/'

# ─── Figure 1: All three traits ───────────────────────────────────────────────
plt.style.use('seaborn-v0_8-whitegrid')
fig, axes = plt.subplots(3, 1, figsize=(10, 10), dpi=150, sharex=True)

traits = [
    ('yearly_mean_resistance', 'Mean resistance (r)', 0.15),
    ('yearly_mean_tolerance',  'Mean tolerance (t)',  0.10),
    ('yearly_mean_recovery',   'Mean recovery (c)',   0.02),
]

for ax, (key, ylabel, init_val) in zip(axes, traits):
    for reg in REGIONS:
        vals = rd[reg][key]
        ax.plot(YEARS, vals, color=COLORS[reg], label=reg, linewidth=1.8, marker='o', markersize=3)
    ax.axhline(init_val, color='gray', linestyle='--', alpha=0.5, linewidth=1)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.tick_params(labelsize=10)

axes[0].set_title('Host Trait Evolution Across Regions', fontsize=14, fontweight='bold', pad=10)
axes[-1].set_xlabel('Year', fontsize=12)
axes[-1].set_xticks(YEARS)
axes[-1].set_xticklabels(YEARS, rotation=45, ha='right')

# Shared legend at bottom
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', ncol=6, fontsize=10,
           bbox_to_anchor=(0.5, -0.02), frameon=True, fancybox=True)

fig.tight_layout(rect=[0, 0.03, 1, 1])
path1 = OUTDIR + 'fig_all_traits_evolution.png'
fig.savefig(path1, bbox_inches='tight')
plt.close(fig)
print(f'Saved: {path1}')

# ─── Figure 2: Survival probability ──────────────────────────────────────────
fig2, ax2 = plt.subplots(figsize=(10, 6), dpi=150)

for reg in REGIONS:
    vals = rd[reg]['yearly_mean_resistance']
    ax2.plot(YEARS, vals, color=COLORS[reg], label=reg, linewidth=2, marker='o', markersize=4)

# Reference line
ax2.axhline(0.15, color='gray', linestyle='--', linewidth=1.2, alpha=0.7)
ax2.annotate('Initial population mean', xy=(2012.3, 0.153), fontsize=9, color='gray')

# Annotation
ax2.annotate(
    'Even after 13 years of selection,\nr = 0.25 means 75% chance\nof infection per exposure',
    xy=(2020, 0.25), xytext=(2015.5, 0.40),
    fontsize=9, fontstyle='italic',
    arrowprops=dict(arrowstyle='->', color='#555555', lw=1.2),
    bbox=dict(boxstyle='round,pad=0.4', facecolor='#ffffcc', edgecolor='#cccccc', alpha=0.9),
    color='#333333',
)

ax2.set_ylim(0, 0.5)
ax2.set_xlim(2011.5, 2024.5)
ax2.set_xticks(YEARS)
ax2.set_xticklabels(YEARS, rotation=45, ha='right')
ax2.set_xlabel('Year', fontsize=12)
ax2.set_ylabel('Probability of avoiding infection per exposure', fontsize=12)
ax2.set_title('Host Immune Defense Over Time', fontsize=14, fontweight='bold')
ax2.legend(loc='upper left', fontsize=10, frameon=True, fancybox=True)
ax2.tick_params(labelsize=10)

fig2.tight_layout()
path2 = OUTDIR + 'fig_survival_probability.png'
fig2.savefig(path2, bbox_inches='tight')
plt.close(fig2)
print(f'Saved: {path2}')

# Report sizes
import os
for p in [path1, path2]:
    sz = os.path.getsize(p)
    print(f'  {p}: {sz/1024:.1f} KB')
