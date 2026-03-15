#!/usr/bin/env python3
"""Generate all figures for the W75-W84 calibration report."""

import json
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import os
from sswd_evoepi.metrics import RECOVERY_TARGETS

# ── Configuration ──────────────────────────────────────────────────────────
REPORT_DIR = '/home/starbot/.openclaw/workspace/sswd-evoepi/reports/khalf_sweep_w75_w84'
FIG_DIR = os.path.join(REPORT_DIR, 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

YEARS = list(range(2012, 2025))  # 13 years
TARGET_REGIONS = ['AK-PWS', 'AK-FN', 'AK-FS', 'BC-N', 'SS-S', 'JDF', 'OR', 'CA-N']
# CENTRALIZED: moved to sswd_evoepi.metrics
TARGETS = RECOVERY_TARGETS
# TARGETS = {'AK-PWS': 0.50, 'AK-FN': 0.50, 'AK-FS': 0.20, 'BC-N': 0.20,
#            'SS-S': 0.05, 'JDF': 0.02, 'OR': 0.0025, 'CA-N': 0.001}
TARGET_PCT = {k: v*100 for k, v in TARGETS.items()}

REGION_COLORS = {
    'AK-PWS': '#1f77b4', 'AK-FN': '#4a90d9', 'AK-FS': '#7eb3e8',
    'BC-N': '#2ca02c', 'SS-S': '#ff7f0e', 'JDF': '#d62728',
    'OR': '#9467bd', 'CA-N': '#e377c2'
}

# All regions in the model
ALL_REGIONS = ['AK-AL', 'AK-WG', 'AK-OC', 'AK-EG', 'AK-PWS', 'AK-FN', 'AK-FS',
               'BC-N', 'BC-C', 'SS-N', 'SS-S', 'JDF', 'WA-O', 'OR', 'CA-N', 'CA-C', 'CA-S']

ALL_REGION_COLORS = {
    'AK-AL': '#08306b', 'AK-WG': '#08519c', 'AK-OC': '#2171b5',
    'AK-EG': '#4292c6', 'AK-PWS': '#1f77b4', 'AK-FN': '#4a90d9',
    'AK-FS': '#7eb3e8', 'BC-N': '#2ca02c', 'BC-C': '#74c476',
    'SS-N': '#fdae6b', 'SS-S': '#ff7f0e', 'JDF': '#d62728',
    'WA-O': '#e6550d', 'OR': '#9467bd', 'CA-N': '#e377c2',
    'CA-C': '#bcbd22', 'CA-S': '#17becf'
}

RUN_IDS = [71, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84]
RUN_PARAMS = {
    71:  {'K_half': 400, 'Pmax': 2000, 'floor': 100, 'label': 'W71 (K½=400K)'},
    75:  {'K_half': 600, 'Pmax': 2000, 'floor': 100, 'label': 'W75 (K½=600K, Pm=2000)'},
    76:  {'K_half': 800, 'Pmax': 2000, 'floor': 100, 'label': 'W76 (K½=800K, Pm=2000)'},
    77:  {'K_half': 1000, 'Pmax': 2000, 'floor': 100, 'label': 'W77 (K½=1M, Pm=2000)'},
    78:  {'K_half': 600, 'Pmax': 1500, 'floor': 100, 'label': 'W78 (K½=600K, Pm=1500)'},
    79:  {'K_half': 800, 'Pmax': 1500, 'floor': 100, 'label': 'W79 (K½=800K, Pm=1500)'},
    80:  {'K_half': 1000, 'Pmax': 1500, 'floor': 100, 'label': 'W80 (K½=1M, Pm=1500)'},
    81:  {'K_half': 600, 'Pmax': 1000, 'floor': 100, 'label': 'W81 (K½=600K, Pm=1000)'},
    82:  {'K_half': 800, 'Pmax': 1000, 'floor': 100, 'label': 'W82 (K½=800K, Pm=1000)'},
    83:  {'K_half': 800, 'Pmax': 2000, 'floor': 75,  'label': 'W83 (K½=800K, Pm=2000, f=75)'},
    84:  {'K_half': 800, 'Pmax': 1500, 'floor': 75,  'label': 'W84 (K½=800K, Pm=1500, f=75)'},
}

# ── Load data ──────────────────────────────────────────────────────────────
def load_run(wid):
    with open(f'/tmp/W{wid}.json') as f:
        data = json.load(f)
    result = data['results'][0]
    scoring = result['scoring']
    details = result['region_details']
    
    # Extract per-region recovery fractions (population / peak)
    recovery_ts = {}
    for region, info in details.items():
        pops = info['yearly_totals']
        peak = max(pops[:3])  # max of first 3 years
        if peak > 0:
            recovery_ts[region] = [p / peak for p in pops]
        else:
            recovery_ts[region] = [0.0] * len(pops)
    
    # Per-target-region actual_pct
    actual_pct = {}
    for reg in TARGET_REGIONS:
        if reg in scoring['per_region']:
            actual_pct[reg] = scoring['per_region'][reg]['actual_pct']
    
    return {
        'scoring': scoring,
        'details': details,
        'recovery_ts': recovery_ts,
        'actual_pct': actual_pct,
        'rmse': scoring['rmse_log'],
        'within_2x': scoring['within_2x'],
        'within_5x': scoring['within_5x'],
    }

runs = {}
for wid in RUN_IDS:
    runs[wid] = load_run(wid)
    print(f"W{wid}: RMSE={runs[wid]['rmse']:.3f}, " + 
          ", ".join(f"{r}={runs[wid]['actual_pct'].get(r, 0):.1f}%" for r in TARGET_REGIONS))

# ── Print summary table ──────────────────────────────────────────────────
print("\n=== SUMMARY TABLE ===")
print(f"{'Run':<6} {'K½':>6} {'Pmax':>5} {'floor':>5} | " + " ".join(f"{r:>7}" for r in TARGET_REGIONS) + f" | {'RMSE':>6}")
print("-" * 120)
for wid in RUN_IDS:
    p = RUN_PARAMS[wid]
    vals = " ".join(f"{runs[wid]['actual_pct'].get(r, 0):>6.1f}%" for r in TARGET_REGIONS)
    print(f"W{wid:<4} {p['K_half']:>5}K {p['Pmax']:>5} {p['floor']:>5} | {vals} | {runs[wid]['rmse']:>6.3f}")
print("Targets: " + " ".join(f"{TARGET_PCT[r]:>6.1f}%" for r in TARGET_REGIONS))


# ═══════════════════════════════════════════════════════════════════════════
# FIG 1: 2×2 key runs recovery time series
# ═══════════════════════════════════════════════════════════════════════════
def plot_recovery_panel(ax, wid, title):
    """Plot recovery time series for target regions on one axis."""
    data = runs[wid]
    for reg in TARGET_REGIONS:
        if reg in data['recovery_ts']:
            ts = data['recovery_ts'][reg]
            ax.plot(YEARS, ts, color=REGION_COLORS[reg], linewidth=1.5, label=reg)
    # Add target lines
    for reg in TARGET_REGIONS:
        ax.axhline(y=TARGETS[reg], color=REGION_COLORS[reg], linestyle=':', alpha=0.3, linewidth=0.8)
    ax.set_title(title, fontsize=10, fontweight='bold')
    ax.set_xlim(2012, 2024)
    ax.set_ylim(0, 1.1)
    ax.set_ylabel('Recovery fraction')
    ax.set_xlabel('Year')
    ax.grid(True, alpha=0.3)

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
plot_recovery_panel(axes[0, 0], 71, 'W71 — K½=400K (Reference)')
plot_recovery_panel(axes[0, 1], 75, 'W75 — K½=600K')
plot_recovery_panel(axes[1, 0], 76, 'W76 — K½=800K')
plot_recovery_panel(axes[1, 1], 77, 'W77 — K½=1M')

# Add shared legend
handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', ncol=8, fontsize=8, bbox_to_anchor=(0.5, -0.02))
fig.suptitle('Recovery Time Series: Key K½ Comparison', fontsize=13, fontweight='bold')
plt.tight_layout(rect=[0, 0.04, 1, 0.96])
plt.savefig(os.path.join(FIG_DIR, 'fig1_key_runs_recovery.png'), dpi=200, bbox_inches='tight')
plt.close()
print("✓ Fig 1 saved")


# ═══════════════════════════════════════════════════════════════════════════
# FIG 2: All 10 runs recovery time series grid (+ W71 reference)
# ═══════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(3, 4, figsize=(20, 13))
all_wids = [71, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84]

for idx, wid in enumerate(all_wids):
    row, col = divmod(idx, 4)
    ax = axes[row, col]
    p = RUN_PARAMS[wid]
    title = f"W{wid}: K½={p['K_half']}K"
    if wid == 71:
        title += " [REF]"
    else:
        title += f", Pm={p['Pmax']}, f={p['floor']}"
    plot_recovery_panel(ax, wid, title)
    # Add RMSE annotation
    ax.text(0.98, 0.98, f"RMSE={runs[wid]['rmse']:.2f}", transform=ax.transAxes,
            fontsize=8, ha='right', va='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

# Hide unused subplot
axes[2, 3].set_visible(False)

handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', ncol=8, fontsize=9, bbox_to_anchor=(0.5, -0.01))
fig.suptitle('All Runs: Recovery Time Series (W71 + W75–W84)', fontsize=14, fontweight='bold')
plt.tight_layout(rect=[0, 0.03, 1, 0.96])
plt.savefig(os.path.join(FIG_DIR, 'fig2_all_runs_grid.png'), dpi=200, bbox_inches='tight')
plt.close()
print("✓ Fig 2 saved")


# ═══════════════════════════════════════════════════════════════════════════
# FIG 3: K_half effect — same floor=100, Pmax=2000, varying K_half
# ═══════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
khalf_runs = [(75, 'K½=600K'), (76, 'K½=800K'), (77, 'K½=1M')]

for ax, (wid, label) in zip(axes, khalf_runs):
    data = runs[wid]
    for reg in TARGET_REGIONS:
        if reg in data['recovery_ts']:
            ts = data['recovery_ts'][reg]
            ax.plot(YEARS, ts, color=REGION_COLORS[reg], linewidth=2, label=reg)
            # Mark final value
            ax.annotate(f'{ts[-1]*100:.0f}%', xy=(2024, ts[-1]), fontsize=7,
                       color=REGION_COLORS[reg], ha='left', va='center',
                       xytext=(5, 0), textcoords='offset points')
    for reg in TARGET_REGIONS:
        ax.axhline(y=TARGETS[reg], color=REGION_COLORS[reg], linestyle=':', alpha=0.3, linewidth=0.8)
    ax.set_title(f'{label} (RMSE={runs[wid]["rmse"]:.2f})', fontsize=11, fontweight='bold')
    ax.set_xlim(2012, 2024.8)
    ax.set_ylim(0, 1.1)
    ax.set_ylabel('Recovery fraction')
    ax.set_xlabel('Year')
    ax.grid(True, alpha=0.3)

handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', ncol=8, fontsize=9, bbox_to_anchor=(0.5, -0.04))
fig.suptitle('K½ Effect on Recovery (Pmax=2000, floor=100)', fontsize=13, fontweight='bold')
plt.tight_layout(rect=[0, 0.06, 1, 0.95])
plt.savefig(os.path.join(FIG_DIR, 'fig3_khalf_effect.png'), dpi=200, bbox_inches='tight')
plt.close()
print("✓ Fig 3 saved")


# ═══════════════════════════════════════════════════════════════════════════
# FIG 4: Recovery fraction bar chart — W71 vs W76 vs targets
# ═══════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(12, 6))
x = np.arange(len(TARGET_REGIONS))
width = 0.25

target_vals = [TARGET_PCT[r] for r in TARGET_REGIONS]
w71_vals = [runs[71]['actual_pct'].get(r, 0) for r in TARGET_REGIONS]
w76_vals = [runs[76]['actual_pct'].get(r, 0) for r in TARGET_REGIONS]

bars1 = ax.bar(x - width, target_vals, width, label='Target', color='#333333', alpha=0.6)
bars2 = ax.bar(x, w71_vals, width, label='W71 (K½=400K)', color='#1f77b4', alpha=0.8)
bars3 = ax.bar(x + width, w76_vals, width, label='W76 (K½=800K)', color='#ff7f0e', alpha=0.8)

# Add value labels
for bars in [bars1, bars2, bars3]:
    for bar in bars:
        height = bar.get_height()
        if height > 0.5:
            ax.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                    f'{height:.1f}%', ha='center', va='bottom', fontsize=7, rotation=45)

ax.set_ylabel('Recovery (%)', fontsize=11)
ax.set_title('Recovery Comparison: W71 vs W76 vs Targets', fontsize=13, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(TARGET_REGIONS, fontsize=10)
ax.legend(fontsize=10)
ax.set_yscale('log')
ax.set_ylim(0.05, 100)
ax.grid(True, alpha=0.3, axis='y')

# Add annotations for gradient ratio
gradient_w71 = runs[71]['actual_pct']['AK-PWS'] / runs[71]['actual_pct']['CA-N']
gradient_w76 = runs[76]['actual_pct']['AK-PWS'] / runs[76]['actual_pct']['CA-N']
ax.text(0.02, 0.98, f'Gradient (AK-PWS/CA-N):\nW71: {gradient_w71:.0f}×\nW76: {gradient_w76:.1f}×',
        transform=ax.transAxes, fontsize=9, va='top',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'fig4_recovery_comparison.png'), dpi=200, bbox_inches='tight')
plt.close()
print("✓ Fig 4 saved")


# ═══════════════════════════════════════════════════════════════════════════
# FIG 5: P_env_max independence proof — overlay W75/W78/W81
# ═══════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

# Group: K_half=600K with Pmax=2000/1500/1000
group1 = [(75, 'Pmax=2000', '-'), (78, 'Pmax=1500', '--'), (81, 'Pmax=1000', ':')]
group2 = [(76, 'Pmax=2000', '-'), (79, 'Pmax=1500', '--'), (82, 'Pmax=1000', ':')]
group3 = [(77, 'Pmax=2000', '-'), (80, 'Pmax=1500', '--')]

groups = [
    (group1, 'K½=600K: Pmax Independence'),
    (group2, 'K½=800K: Pmax Independence'),
    (group3, 'K½=1M: Pmax Independence'),
]

for ax, (group, title) in zip(axes, groups):
    for reg in ['AK-PWS', 'SS-S', 'CA-N']:
        for wid, plabel, ls in group:
            ts = runs[wid]['recovery_ts'].get(reg, [0]*13)
            ax.plot(YEARS, ts, color=REGION_COLORS[reg], linewidth=2, linestyle=ls,
                    label=f'{reg} W{wid} ({plabel})')
    ax.set_title(title, fontsize=10, fontweight='bold')
    ax.set_xlim(2012, 2024)
    ax.set_ylim(0, 1.1)
    ax.set_ylabel('Recovery fraction')
    ax.set_xlabel('Year')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=6, loc='upper left')

fig.suptitle('P_env_max Has ZERO Effect: Lines Overlap Completely', fontsize=13, fontweight='bold', color='red')
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig(os.path.join(FIG_DIR, 'fig5_pmax_independence.png'), dpi=200, bbox_inches='tight')
plt.close()
print("✓ Fig 5 saved")


# ═══════════════════════════════════════════════════════════════════════════
# FIG 6: RMSE comparison across all runs
# ═══════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(12, 5))
wids = [71, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84]
rmse_vals = [runs[w]['rmse'] for w in wids]
labels = [f"W{w}\n{RUN_PARAMS[w]['label'].split('(')[1].rstrip(')')}" for w in wids]

colors = []
for w in wids:
    if w == 71:
        colors.append('#2ca02c')  # reference green
    elif runs[w]['rmse'] == min(rmse_vals[1:]):  # best of sweep
        colors.append('#1f77b4')
    elif runs[w]['rmse'] == max(rmse_vals):
        colors.append('#d62728')
    else:
        colors.append('#7f7f7f')

bars = ax.bar(range(len(wids)), rmse_vals, color=colors, alpha=0.8, edgecolor='black', linewidth=0.5)

for i, (val, wid) in enumerate(zip(rmse_vals, wids)):
    ax.text(i, val + 0.02, f'{val:.3f}', ha='center', va='bottom', fontsize=8, fontweight='bold')

ax.set_xticks(range(len(wids)))
ax.set_xticklabels(labels, fontsize=7)
ax.set_ylabel('RMSE (log scale)', fontsize=11)
ax.set_title('RMSE Comparison Across All Runs', fontsize=13, fontweight='bold')
ax.grid(True, alpha=0.3, axis='y')

# Add legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor='#2ca02c', label='W71 Reference'),
                   Patch(facecolor='#1f77b4', label='Best of sweep'),
                   Patch(facecolor='#d62728', label='Worst'),
                   Patch(facecolor='#7f7f7f', label='Other')]
ax.legend(handles=legend_elements, fontsize=9, loc='upper right')

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'fig6_rmse_comparison.png'), dpi=200, bbox_inches='tight')
plt.close()
print("✓ Fig 6 saved")


# ═══════════════════════════════════════════════════════════════════════════
# FIG 7: Disease deaths + recruits time series (W76)
# ═══════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(2, 1, figsize=(14, 9), sharex=True)

for reg in TARGET_REGIONS:
    if reg in runs[76]['details']:
        info = runs[76]['details'][reg]
        deaths = info.get('yearly_disease_deaths', [0]*13)
        recruits = info.get('yearly_recruits', [0]*13)
        axes[0].plot(YEARS, [d/1000 for d in deaths], color=REGION_COLORS[reg], 
                     linewidth=1.5, label=reg, marker='o', markersize=3)
        axes[1].plot(YEARS, [r/1000 for r in recruits], color=REGION_COLORS[reg],
                     linewidth=1.5, label=reg, marker='s', markersize=3)

axes[0].set_title('W76 — Disease Deaths per Year (thousands)', fontsize=11, fontweight='bold')
axes[0].set_ylabel('Deaths (×1000)')
axes[0].legend(fontsize=8, ncol=4, loc='upper right')
axes[0].grid(True, alpha=0.3)

axes[1].set_title('W76 — Recruits per Year (thousands)', fontsize=11, fontweight='bold')
axes[1].set_ylabel('Recruits (×1000)')
axes[1].set_xlabel('Year')
axes[1].legend(fontsize=8, ncol=4, loc='upper right')
axes[1].grid(True, alpha=0.3)

fig.suptitle('W76 (K½=800K): Disease Dynamics', fontsize=13, fontweight='bold')
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig(os.path.join(FIG_DIR, 'fig7_disease_dynamics.png'), dpi=200, bbox_inches='tight')
plt.close()
print("✓ Fig 7 saved")


# ═══════════════════════════════════════════════════════════════════════════
# FIG 8: Full 18-region recovery comparison (W71 vs W76)
# ═══════════════════════════════════════════════════════════════════════════
# Get all regions present in both runs
common_regions = sorted(set(runs[71]['recovery_ts'].keys()) & set(runs[76]['recovery_ts'].keys()))
# Order by geography (north to south roughly)
ordered_regions = [r for r in ALL_REGIONS if r in common_regions]

fig, axes = plt.subplots(2, 1, figsize=(16, 10))

# W71
ax = axes[0]
for reg in ordered_regions:
    ts = runs[71]['recovery_ts'][reg]
    ax.plot(YEARS, ts, color=ALL_REGION_COLORS.get(reg, '#999999'), linewidth=1.5, label=reg)
ax.set_title('W71 (K½=400K) — All Regions', fontsize=11, fontweight='bold')
ax.set_ylabel('Recovery fraction')
ax.set_xlim(2012, 2024)
ax.set_ylim(0, 1.1)
ax.legend(fontsize=7, ncol=6, loc='upper left')
ax.grid(True, alpha=0.3)

# W76
ax = axes[1]
for reg in ordered_regions:
    ts = runs[76]['recovery_ts'][reg]
    ax.plot(YEARS, ts, color=ALL_REGION_COLORS.get(reg, '#999999'), linewidth=1.5, label=reg)
ax.set_title('W76 (K½=800K) — All Regions', fontsize=11, fontweight='bold')
ax.set_ylabel('Recovery fraction')
ax.set_xlabel('Year')
ax.set_xlim(2012, 2024)
ax.set_ylim(0, 1.1)
ax.legend(fontsize=7, ncol=6, loc='upper left')
ax.grid(True, alpha=0.3)

fig.suptitle('Full Region Comparison: W71 vs W76', fontsize=13, fontweight='bold')
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig(os.path.join(FIG_DIR, 'fig8_full_region_comparison.png'), dpi=200, bbox_inches='tight')
plt.close()
print("✓ Fig 8 saved")

# ═══════════════════════════════════════════════════════════════════════════
# Save summary data for LaTeX table
# ═══════════════════════════════════════════════════════════════════════════
summary = {}
for wid in RUN_IDS:
    summary[f'W{wid}'] = {
        'K_half': RUN_PARAMS[wid]['K_half'],
        'Pmax': RUN_PARAMS[wid]['Pmax'],
        'floor': RUN_PARAMS[wid]['floor'],
        'rmse': runs[wid]['rmse'],
        'within_2x': runs[wid]['within_2x'],
        'within_5x': runs[wid]['within_5x'],
        'actual_pct': runs[wid]['actual_pct'],
    }

with open(os.path.join(REPORT_DIR, 'summary_data.json'), 'w') as f:
    json.dump(summary, f, indent=2)

print("\n✓ All figures saved to", FIG_DIR)
print("✓ Summary data saved to summary_data.json")
