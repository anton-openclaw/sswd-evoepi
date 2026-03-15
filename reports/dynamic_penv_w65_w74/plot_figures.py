#!/usr/bin/env python3
"""Generate figures for Dynamic P_env calibration report (W65-W74)."""

import json
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.ticker as mticker
from sswd_evoepi.metrics import RECOVERY_TARGETS

# Load data
with open("/tmp/w65_w74_all_data.json") as f:
    all_data = json.load(f)

OUT = "/home/starbot/.openclaw/workspace/sswd-evoepi/reports/dynamic_penv_w65_w74/figures"
import os
os.makedirs(OUT, exist_ok=True)

# CENTRALIZED: moved to sswd_evoepi.metrics
targets = RECOVERY_TARGETS
# targets = {"AK-PWS": 0.50, "AK-FN": 0.50, "AK-FS": 0.20, "BC-N": 0.20,
#            "SS-S": 0.05, "JDF": 0.02, "OR": 0.0025, "CA-N": 0.001}

scored_regions = ["AK-PWS", "AK-FN", "AK-FS", "BC-N", "SS-S", "JDF", "OR", "CA-N"]
years = np.arange(2012, 2025)  # 13 years

# Color scheme for regions (north=blue, south=red)
region_colors = {
    "AK-PWS": "#1f77b4", "AK-FN": "#4a90d9", "AK-FS": "#7eb3e8",
    "BC-N": "#2ca02c", "SS-S": "#ff7f0e", "JDF": "#d62728",
    "OR": "#9467bd", "CA-N": "#e377c2"
}

# ══════════════════════════════════════════════════════════════════════
# FIGURE 1: Recovery time series for key runs (2×2 panel)
# ══════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharex=True, sharey=True)
fig.suptitle("Recovery Trajectories: Dynamic P_env vs Static P_env", fontsize=14, fontweight='bold')

key_runs = [
    ("W62", "W62 (Static P_env) — Best static run"),
    ("W71", "W71 (floor=100, δ=0.02) — Best RMSE"),
    ("W68", "W68 (floor=50, δ=0.02)"),
    ("W74", "W74 (floor=50, δ=0.05, α=0.2)"),
]

for ax, (wname, title) in zip(axes.flat, key_runs):
    d = all_data[wname]
    for reg in scored_regions:
        if reg not in d["regions"]:
            continue
        yt = d["regions"][reg]["yearly_totals"]
        if not yt or len(yt) < 2:
            continue
        peak = max(yt[:3]) if len(yt) >= 3 else max(yt)
        if peak == 0:
            continue
        frac = [p / peak for p in yt[:13]]
        ax.plot(years[:len(frac)], frac, '-o', markersize=3, linewidth=1.5,
                color=region_colors.get(reg, 'gray'), label=reg)
    
    # Target markers on right edge
    for reg in scored_regions:
        ax.axhline(y=targets[reg], color=region_colors.get(reg, 'gray'),
                   linestyle=':', alpha=0.3, linewidth=0.8)
    
    rmse = d.get("rmse", "?")
    ax.set_title(f"{title}\nRMSE={rmse:.3f}" if isinstance(rmse, float) else title, fontsize=10)
    ax.set_ylim(-0.05, 1.05)
    ax.set_ylabel("Population / Pre-crash Peak")
    ax.set_xlabel("Year")
    ax.grid(True, alpha=0.3)
    ax.axvspan(2020.5, 2021.5, color='lightblue', alpha=0.3, label='La Niña' if ax == axes[0,0] else None)

# Single legend
handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', ncol=5, fontsize=8,
           bbox_to_anchor=(0.5, -0.02))
plt.tight_layout(rect=[0, 0.04, 1, 0.96])
plt.savefig(f"{OUT}/fig1_recovery_timeseries.png", dpi=150, bbox_inches='tight')
plt.close()
print("Fig 1 done")

# ══════════════════════════════════════════════════════════════════════
# FIGURE 2: All W65-W74 recovery time series (3×3 grid + W74)
# ══════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(3, 4, figsize=(18, 12), sharex=True, sharey=True)
fig.suptitle("Dynamic P_env: All 10 Runs — Recovery Trajectories by Region", fontsize=14, fontweight='bold')

all_runs = [
    ("W65", "f=25 δ=0.02"), ("W66", "f=25 δ=0.05"), ("W67", "f=25 δ=0.1"),
    ("W68", "f=50 δ=0.02"), ("W69", "f=50 δ=0.05"), ("W70", "f=50 δ=0.1"),
    ("W71", "f=100 δ=0.02"), ("W72", "f=100 δ=0.05"), ("W73", "f=100 δ=0.1"),
    ("W74", "f=50 δ=0.05 α=0.2"),
]

for idx, (wname, label) in enumerate(all_runs):
    row, col = divmod(idx, 4) if idx < 9 else (2, 3)
    if idx < 3:
        row, col = 0, idx
    elif idx < 6:
        row, col = 1, idx - 3
    elif idx < 9:
        row, col = 2, idx - 6
    else:
        row, col = 2, 3
    
    ax = axes[row, col]
    d = all_data[wname]
    
    for reg in scored_regions:
        if reg not in d["regions"]:
            continue
        yt = d["regions"][reg]["yearly_totals"]
        if not yt or len(yt) < 2:
            continue
        peak = max(yt[:3]) if len(yt) >= 3 else max(yt)
        if peak == 0:
            continue
        frac = [p / peak for p in yt[:13]]
        ax.plot(years[:len(frac)], frac, '-o', markersize=2, linewidth=1.2,
                color=region_colors.get(reg, 'gray'), label=reg if idx == 0 else None)
    
    rmse = d.get("rmse", "?")
    rmse_str = f" RMSE={rmse:.3f}" if isinstance(rmse, float) else ""
    ax.set_title(f"{wname} ({label}){rmse_str}", fontsize=9)
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, alpha=0.3)

# Hide unused subplot (row 0 col 3, row 1 col 3)
axes[0, 3].set_visible(False)
axes[1, 3].set_visible(False)

# Add legend
handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper right', fontsize=9, bbox_to_anchor=(0.98, 0.95))

for ax in axes[-1, :]:
    ax.set_xlabel("Year")
for ax in axes[:, 0]:
    ax.set_ylabel("Pop / Peak")

plt.tight_layout()
plt.savefig(f"{OUT}/fig2_all_runs_timeseries.png", dpi=150, bbox_inches='tight')
plt.close()
print("Fig 2 done")

# ══════════════════════════════════════════════════════════════════════
# FIGURE 3: δ_env effect — same floor, varying decay rate
# ══════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=True)
fig.suptitle("Effect of δ_env (Decay Rate) on Recovery — floor=100", fontsize=13, fontweight='bold')

delta_runs = [("W71", "δ=0.02 (35d ½-life)"), ("W72", "δ=0.05 (14d ½-life)"), ("W73", "δ=0.1 (7d ½-life)")]
focus_regions = ["AK-PWS", "BC-N", "SS-S", "CA-N"]

for ax, (wname, label) in zip(axes, delta_runs):
    d = all_data[wname]
    for reg in focus_regions:
        if reg not in d["regions"]:
            continue
        yt = d["regions"][reg]["yearly_totals"]
        if not yt:
            continue
        peak = max(yt[:3]) if len(yt) >= 3 else max(yt)
        if peak == 0:
            continue
        frac = [p / peak for p in yt[:13]]
        ax.plot(years[:len(frac)], frac, '-o', markersize=4, linewidth=2,
                color=region_colors.get(reg, 'gray'), label=reg)
        # Target line
        if reg in targets:
            ax.axhline(y=targets[reg], color=region_colors.get(reg, 'gray'),
                       linestyle='--', alpha=0.4, linewidth=1)
    
    rmse = d.get("rmse", "?")
    ax.set_title(f"{wname}: {label}\nRMSE={rmse:.3f}" if isinstance(rmse, float) else label, fontsize=10)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel("Year")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc='upper right')

axes[0].set_ylabel("Population / Pre-crash Peak")
plt.tight_layout()
plt.savefig(f"{OUT}/fig3_delta_effect.png", dpi=150, bbox_inches='tight')
plt.close()
print("Fig 3 done")

# ══════════════════════════════════════════════════════════════════════
# FIGURE 4: Recovery fraction bar chart — W62 vs W71 vs targets
# ══════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(12, 6))

x = np.arange(len(scored_regions))
width = 0.25

# Target values
target_vals = [targets[r] * 100 for r in scored_regions]

# W62 (static)
w62_vals = [all_data["W62"]["recovery"].get(r, 0) * 100 for r in scored_regions]

# W71 (dynamic best)
w71_vals = [all_data["W71"]["recovery"].get(r, 0) * 100 for r in scored_regions]

bars1 = ax.bar(x - width, target_vals, width, label='Target', color='#2ecc71', alpha=0.8, edgecolor='black', linewidth=0.5)
bars2 = ax.bar(x, w62_vals, width, label='W62 (Static P_env)', color='#3498db', alpha=0.8, edgecolor='black', linewidth=0.5)
bars3 = ax.bar(x + width, w71_vals, width, label='W71 (Dynamic P_env)', color='#e74c3c', alpha=0.8, edgecolor='black', linewidth=0.5)

ax.set_xlabel('Region (North → South)', fontsize=12)
ax.set_ylabel('Recovery Fraction (%)', fontsize=12)
ax.set_title('Static vs Dynamic P_env: Recovery by Region\n(W62: RMSE=0.757 vs W71: RMSE=0.692)', fontsize=13, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(scored_regions, fontsize=10)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3, axis='y')
ax.set_yscale('symlog', linthresh=1)
ax.set_ylim(0, 100)

plt.tight_layout()
plt.savefig(f"{OUT}/fig4_static_vs_dynamic.png", dpi=150, bbox_inches='tight')
plt.close()
print("Fig 4 done")

# ══════════════════════════════════════════════════════════════════════
# FIGURE 5: Resistance evolution comparison (W62 vs W71)
# ══════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
fig.suptitle("Resistance Evolution: Static vs Dynamic P_env", fontsize=13, fontweight='bold')

for ax, (wname, title) in zip(axes, [("W62", "W62 (Static)"), ("W71", "W71 (Dynamic)")]):
    d = all_data[wname]
    for reg in ["AK-PWS", "BC-N", "SS-S", "JDF", "OR", "CA-N"]:
        if reg not in d["regions"]:
            continue
        yr = d["regions"][reg].get("yearly_mean_resistance", [])
        if not yr:
            continue
        ax.plot(years[:len(yr)], yr[:13], '-o', markersize=3, linewidth=1.5,
                color=region_colors.get(reg, 'gray'), label=reg)
    
    ax.axhline(y=0.15, color='gray', linestyle='--', alpha=0.5, linewidth=1, label='Initial (0.15)')
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("Year")
    ax.set_ylabel("Mean Resistance (r̄)")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc='upper left')

plt.tight_layout()
plt.savefig(f"{OUT}/fig5_resistance_evolution.png", dpi=150, bbox_inches='tight')
plt.close()
print("Fig 5 done")

# ══════════════════════════════════════════════════════════════════════
# FIGURE 6: RMSE heatmap — floor × delta
# ══════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(8, 6))

floors = [25, 50, 100]
deltas = [0.02, 0.05, 0.1]
rmse_grid = np.zeros((3, 3))
labels_grid = [["" for _ in range(3)] for _ in range(3)]

run_map = {
    (25, 0.02): "W65", (25, 0.05): "W66", (25, 0.1): "W67",
    (50, 0.02): "W68", (50, 0.05): "W69", (50, 0.1): "W70",
    (100, 0.02): "W71", (100, 0.05): "W72", (100, 0.1): "W73",
}

for i, f in enumerate(floors):
    for j, d in enumerate(deltas):
        wname = run_map[(f, d)]
        rmse = all_data[wname].get("rmse", 2.0)
        rmse_grid[i, j] = rmse if isinstance(rmse, (int, float)) else 2.0
        labels_grid[i][j] = f"{wname}\n{rmse_grid[i,j]:.3f}"

im = ax.imshow(rmse_grid, cmap='RdYlGn_r', aspect='auto', vmin=0.6, vmax=1.6)
cbar = plt.colorbar(im, ax=ax, label='RMSE (log-space)')

for i in range(3):
    for j in range(3):
        color = 'white' if rmse_grid[i, j] > 1.2 else 'black'
        ax.text(j, i, labels_grid[i][j], ha='center', va='center', fontsize=10, color=color, fontweight='bold')

ax.set_xticks(range(3))
ax.set_xticklabels([f"δ={d}\n({int(0.693/d)}d ½-life)" for d in deltas], fontsize=10)
ax.set_yticks(range(3))
ax.set_yticklabels([f"floor={f}" for f in floors], fontsize=10)
ax.set_xlabel("Decay Rate (δ_env)", fontsize=12)
ax.set_ylabel("Community Floor (P_env_floor)", fontsize=12)
ax.set_title("RMSE Heatmap: floor × δ_env\n(lower = better; W71 is best at 0.692)", fontsize=13, fontweight='bold')

plt.tight_layout()
plt.savefig(f"{OUT}/fig6_rmse_heatmap.png", dpi=150, bbox_inches='tight')
plt.close()
print("Fig 6 done")

# ══════════════════════════════════════════════════════════════════════
# FIGURE 7: Disease deaths + recruits time series (W71)
# ══════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True)
fig.suptitle("W71: Disease Deaths & Recruitment by Region", fontsize=13, fontweight='bold')

d = all_data["W71"]
focus = ["AK-PWS", "BC-N", "SS-S", "JDF", "OR", "CA-N"]

for reg in focus:
    if reg not in d["regions"]:
        continue
    dd = d["regions"][reg].get("yearly_disease_deaths", [])
    rec = d["regions"][reg].get("yearly_recruits", [])
    if dd:
        axes[0].plot(years[:len(dd)], [x/1000 for x in dd[:13]], '-o', markersize=3,
                     linewidth=1.5, color=region_colors.get(reg, 'gray'), label=reg)
    if rec:
        axes[1].plot(years[:len(rec)], [x/1000 for x in rec[:13]], '-o', markersize=3,
                     linewidth=1.5, color=region_colors.get(reg, 'gray'), label=reg)

axes[0].set_ylabel("Disease Deaths (thousands)")
axes[0].set_title("Disease Deaths per Year")
axes[0].legend(fontsize=8, loc='upper right')
axes[0].grid(True, alpha=0.3)

axes[1].set_ylabel("Recruits (thousands)")
axes[1].set_xlabel("Year")
axes[1].set_title("Recruitment per Year")
axes[1].legend(fontsize=8, loc='upper right')
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f"{OUT}/fig7_deaths_recruits.png", dpi=150, bbox_inches='tight')
plt.close()
print("Fig 7 done")

# ══════════════════════════════════════════════════════════════════════
# FIGURE 8: Full 18-region recovery comparison (W62 vs W71)
# ══════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(16, 6))

all_regions = ["AK-WG", "AK-AL", "AK-EG", "AK-PWS", "AK-OC", "AK-FN", "AK-FS",
               "BC-N", "BC-C", "JDF", "SS-N", "SS-S", "WA-O", "OR", "CA-N", "CA-C", "CA-S", "BJ"]

x = np.arange(len(all_regions))
width = 0.35

w62_vals = [all_data["W62"]["recovery"].get(r, 0) * 100 for r in all_regions]
w71_vals = [all_data["W71"]["recovery"].get(r, 0) * 100 for r in all_regions]

ax.bar(x - width/2, w62_vals, width, label='W62 (Static)', color='#3498db', alpha=0.8)
ax.bar(x + width/2, w71_vals, width, label='W71 (Dynamic)', color='#e74c3c', alpha=0.8)

# Add target markers for scored regions
for i, reg in enumerate(all_regions):
    if reg in targets:
        ax.plot(i, targets[reg] * 100, 'k*', markersize=12, zorder=5)

ax.set_xticks(x)
ax.set_xticklabels(all_regions, rotation=45, ha='right', fontsize=9)
ax.set_ylabel("Recovery Fraction (%)")
ax.set_title("All 18 Regions: Static vs Dynamic P_env\n(★ = calibration target)", fontsize=13, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig(f"{OUT}/fig8_all_regions.png", dpi=150, bbox_inches='tight')
plt.close()
print("Fig 8 done")

print(f"\nAll figures saved to {OUT}/")
