#!/usr/bin/env python3
"""Extract data from W01-W04 calibration results and generate all figures."""

import json
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# ============================================================
# Configuration
# ============================================================
BASE = "/home/starbot/.openclaw/workspace/sswd-evoepi/results/calibration"
FIG_DIR = "/home/starbot/.openclaw/workspace/sswd-evoepi/reports/wavefront_calibration/figures"
os.makedirs(FIG_DIR, exist_ok=True)

ROUNDS = ["W01", "W02", "W03", "W04"]
SEEDS = [42, 123, 999]

# Round parameters
ROUND_PARAMS = {
    "W01": {"D_P": 50, "D_P_max_range": 175},
    "W02": {"D_P": 100, "D_P_max_range": 350},
    "W03": {"D_P": 150, "D_P_max_range": 525},
    "W04": {"D_P": 200, "D_P_max_range": 700},
}

# Recovery targets (south to north)
RECOVERY_REGIONS = ["CA-N", "OR", "JDF", "SS-S", "BC-N", "AK-FS", "AK-FN", "AK-PWS"]
RECOVERY_TARGETS = {
    "CA-N": 0.001, "OR": 0.0025, "JDF": 0.02, "SS-S": 0.05,
    "BC-N": 0.2, "AK-FS": 0.2, "AK-FN": 0.5, "AK-PWS": 0.5
}

# Arrival timing targets (south to north)
TIMING_REGIONS = [
    "CA-S", "CA-C", "CA-N", "OR", "WA-O", "JDF", "SS-S", "SS-N",
    "BC-C", "BC-N", "AK-FS", "AK-FN", "AK-PWS", "AK-EG", "AK-OC",
    "AK-WG", "AK-AL"
]
TIMING_TARGETS = {
    "CA-S": 0, "CA-C": 6, "CA-N": 6, "OR": 15, "WA-O": 15,
    "JDF": 26, "SS-S": 26, "SS-N": 26, "BC-C": 26, "BC-N": 26,
    "AK-FS": 26, "AK-FN": 26, "AK-PWS": 33, "AK-EG": 33, "AK-OC": 33,
    "AK-WG": 42, "AK-AL": 42
}

# ============================================================
# Load data
# ============================================================
def load_seed_result(rnd, seed):
    path = os.path.join(BASE, rnd, f"result_seed{seed}.json")
    with open(path) as f:
        return json.load(f)

# Store all data
all_data = {}
for rnd in ROUNDS:
    all_data[rnd] = {}
    for seed in SEEDS:
        all_data[rnd][seed] = load_seed_result(rnd, seed)

# ============================================================
# Extract metrics
# ============================================================

# Recovery fractions per round per seed
recovery_data = {}  # {round: {region: [seed42, seed123, seed999]}}
for rnd in ROUNDS:
    recovery_data[rnd] = {}
    for region in RECOVERY_REGIONS:
        vals = []
        for seed in SEEDS:
            d = all_data[rnd][seed]
            vals.append(d["scoring"]["per_region"][region]["actual"])
        recovery_data[rnd][region] = vals

# Arrival timing per round per seed
timing_data = {}  # {round: {region: [seed42, seed123, seed999]}}
for rnd in ROUNDS:
    timing_data[rnd] = {}
    for region in TIMING_REGIONS:
        vals = []
        for seed in SEEDS:
            d = all_data[rnd][seed]
            vals.append(d["arrival_timing"]["per_region"][region]["actual_months"])
        timing_data[rnd][region] = vals

# RMSE per round per seed
rmse_data = {}  # {round: [seed42, seed123, seed999]}
for rnd in ROUNDS:
    vals = []
    for seed in SEEDS:
        d = all_data[rnd][seed]
        vals.append(d["scoring"]["rmse_log"])
    rmse_data[rnd] = vals

# Print summary for verification
print("=" * 70)
print("SUMMARY OF EXTRACTED DATA")
print("=" * 70)

for rnd in ROUNDS:
    print(f"\n--- {rnd} (D_P={ROUND_PARAMS[rnd]['D_P']}km) ---")
    print(f"  RMSE (log): {[f'{v:.4f}' for v in rmse_data[rnd]]}  mean={np.mean(rmse_data[rnd]):.4f}")
    print(f"  Recovery fractions (mean across seeds):")
    for region in RECOVERY_REGIONS:
        vals = recovery_data[rnd][region]
        target = RECOVERY_TARGETS[region]
        print(f"    {region}: {np.mean(vals)*100:.1f}% (target {target*100:.1f}%)")
    print(f"  Arrival timing (mean across seeds):")
    for region in TIMING_REGIONS:
        vals = timing_data[rnd][region]
        target = TIMING_TARGETS[region]
        print(f"    {region}: {np.mean(vals):.1f} mo (target {target} mo)")

# ============================================================
# Color scheme
# ============================================================
ROUND_COLORS = {
    "W01": "#1f77b4",  # blue
    "W02": "#ff7f0e",  # orange
    "W03": "#2ca02c",  # green
    "W04": "#d62728",  # red
}
TARGET_COLOR = "#333333"

plt.rcParams.update({
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 150,
})

# ============================================================
# Figure 1: Recovery fractions — model vs target (grouped bar)
# ============================================================
fig, ax = plt.subplots(figsize=(10, 5))

n_regions = len(RECOVERY_REGIONS)
n_groups = len(ROUNDS) + 1  # 4 rounds + target
bar_width = 0.15
x = np.arange(n_regions)

# Plot target bars
target_vals = [RECOVERY_TARGETS[r] * 100 for r in RECOVERY_REGIONS]
ax.bar(x - 2*bar_width, target_vals, bar_width, color=TARGET_COLOR, alpha=0.7, label='Target', edgecolor='black', linewidth=0.5)

for i, rnd in enumerate(ROUNDS):
    means = [np.mean(recovery_data[rnd][r]) * 100 for r in RECOVERY_REGIONS]
    stds = [np.std(recovery_data[rnd][r]) * 100 for r in RECOVERY_REGIONS]
    ax.bar(x + (i-1)*bar_width, means, bar_width, 
           color=ROUND_COLORS[rnd], alpha=0.85, 
           label=f'{rnd} (D_P={ROUND_PARAMS[rnd]["D_P"]}km)',
           edgecolor='black', linewidth=0.3,
           yerr=stds, capsize=2, error_kw={'linewidth': 0.8})

ax.set_xlabel('Region (South → North)')
ax.set_ylabel('Recovery Fraction (%)')
ax.set_title('Recovery Fractions: Model vs Target by Region')
ax.set_xticks(x)
ax.set_xticklabels(RECOVERY_REGIONS, rotation=45, ha='right')
ax.legend(loc='upper left', framealpha=0.9)
ax.set_ylim(0, 55)
ax.grid(axis='y', alpha=0.3)
plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "fig1_recovery_bars.pdf"), bbox_inches='tight')
plt.close(fig)
print("\nSaved fig1_recovery_bars.pdf")

# ============================================================
# Figure 2: Wavefront arrival timing — model vs target (grouped bar)
# ============================================================
fig, ax = plt.subplots(figsize=(14, 5.5))

n_regions = len(TIMING_REGIONS)
bar_width = 0.15
x = np.arange(n_regions)

# Plot target bars
target_vals = [TIMING_TARGETS[r] for r in TIMING_REGIONS]
ax.bar(x - 2*bar_width, target_vals, bar_width, color=TARGET_COLOR, alpha=0.7, label='Target', edgecolor='black', linewidth=0.5)

for i, rnd in enumerate(ROUNDS):
    means = [np.mean(timing_data[rnd][r]) for r in TIMING_REGIONS]
    stds = [np.std(timing_data[rnd][r]) for r in TIMING_REGIONS]
    ax.bar(x + (i-1)*bar_width, means, bar_width,
           color=ROUND_COLORS[rnd], alpha=0.85,
           label=f'{rnd} (D_P={ROUND_PARAMS[rnd]["D_P"]}km)',
           edgecolor='black', linewidth=0.3,
           yerr=stds, capsize=2, error_kw={'linewidth': 0.8})

ax.set_xlabel('Region (South → North)')
ax.set_ylabel('Arrival Time (months after origin)')
ax.set_title('Wavefront Arrival Timing: Model vs Target by Region')
ax.set_xticks(x)
ax.set_xticklabels(TIMING_REGIONS, rotation=45, ha='right')
ax.legend(loc='upper left', framealpha=0.9)
ax.grid(axis='y', alpha=0.3)
plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "fig2_timing_bars.pdf"), bbox_inches='tight')
plt.close(fig)
print("Saved fig2_timing_bars.pdf")

# ============================================================
# Figure 3: Arrival timing progression (line plot — the key figure)
# ============================================================
fig, ax = plt.subplots(figsize=(12, 6))

x = np.arange(len(TIMING_REGIONS))

# Target line
target_vals = [TIMING_TARGETS[r] for r in TIMING_REGIONS]
ax.plot(x, target_vals, 'k-o', linewidth=2.5, markersize=7, label='Target', zorder=5)

for rnd in ROUNDS:
    means = [np.mean(timing_data[rnd][r]) for r in TIMING_REGIONS]
    ax.plot(x, means, '-s', color=ROUND_COLORS[rnd], linewidth=1.8, markersize=5,
            label=f'{rnd} (D_P={ROUND_PARAMS[rnd]["D_P"]}km)', alpha=0.85)

ax.fill_between(x, 0, [TIMING_TARGETS[r] for r in TIMING_REGIONS], alpha=0.08, color='gray')

ax.set_xlabel('Region (South → North)')
ax.set_ylabel('Arrival Time (months after origin)')
ax.set_title('Wavefront Propagation: Model vs Target\n(All rounds arrive too fast — gap widens with distance)')
ax.set_xticks(x)
ax.set_xticklabels(TIMING_REGIONS, rotation=45, ha='right')
ax.legend(loc='upper left', framealpha=0.9)
ax.grid(alpha=0.3)
ax.set_ylim(0, 48)

# Add annotation
ax.annotate('Wavefront arrives\n~20 months too early\nin Alaska', 
            xy=(14, 17), xytext=(10, 38),
            arrowprops=dict(arrowstyle='->', color='red', lw=1.5),
            fontsize=10, color='red', ha='center',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', edgecolor='red', alpha=0.8))

plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "fig3_timing_progression.pdf"), bbox_inches='tight')
plt.close(fig)
print("Saved fig3_timing_progression.pdf")

# ============================================================
# Figure 4: Recovery fraction scatter (log-log, model vs target)
# ============================================================
fig, ax = plt.subplots(figsize=(7, 7))

for rnd in ROUNDS:
    targets = []
    actuals = []
    for region in RECOVERY_REGIONS:
        t = RECOVERY_TARGETS[region] * 100
        a = np.mean(recovery_data[rnd][region]) * 100
        targets.append(t)
        actuals.append(a)
    ax.scatter(targets, actuals, s=60, color=ROUND_COLORS[rnd], alpha=0.8, edgecolors='black', linewidths=0.5,
               label=f'{rnd} (D_P={ROUND_PARAMS[rnd]["D_P"]}km)', zorder=3)

# Perfect fit line
lims = [0.05, 70]
ax.plot(lims, lims, 'k--', linewidth=1.5, alpha=0.5, label='Perfect fit')

# 2x and 5x bands
ax.fill_between(lims, [l/2 for l in lims], [l*2 for l in lims], alpha=0.08, color='green', label='2× band')
ax.fill_between(lims, [l/5 for l in lims], [l*5 for l in lims], alpha=0.05, color='blue', label='5× band')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Target Recovery (%)')
ax.set_ylabel('Model Recovery (%)')
ax.set_title('Recovery Fraction: Model vs Target (log-log)')
ax.set_xlim(0.05, 70)
ax.set_ylim(0.05, 70)
ax.legend(loc='upper left', framealpha=0.9)
ax.grid(True, alpha=0.3, which='both')
ax.set_aspect('equal')
plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "fig4_recovery_scatter.pdf"), bbox_inches='tight')
plt.close(fig)
print("Saved fig4_recovery_scatter.pdf")

# ============================================================
# Figure 5: Per-seed variability
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

# 5a: Recovery per seed
ax = axes[0]
seed_colors = {42: '#e377c2', 123: '#7f7f7f', 999: '#bcbd22'}
seed_markers = {42: 'o', 123: 's', 999: '^'}

for rnd_idx, rnd in enumerate(ROUNDS):
    for seed in SEEDS:
        vals = []
        for region in RECOVERY_REGIONS:
            vals.append(recovery_data[rnd][region][SEEDS.index(seed)] * 100)
        offset = (rnd_idx - 1.5) * 0.08
        ax.scatter(np.arange(len(RECOVERY_REGIONS)) + offset, vals, 
                   color=ROUND_COLORS[rnd], marker=seed_markers[seed], s=25, alpha=0.7,
                   edgecolors='black', linewidths=0.3)

# Custom legend
from matplotlib.lines import Line2D
legend_elements = []
for rnd in ROUNDS:
    legend_elements.append(Line2D([0], [0], marker='o', color='w', markerfacecolor=ROUND_COLORS[rnd],
                                  markersize=8, label=rnd))
for seed in SEEDS:
    legend_elements.append(Line2D([0], [0], marker=seed_markers[seed], color='w', markerfacecolor='gray',
                                  markersize=8, label=f'Seed {seed}'))

ax.set_xticks(range(len(RECOVERY_REGIONS)))
ax.set_xticklabels(RECOVERY_REGIONS, rotation=45, ha='right')
ax.set_ylabel('Recovery Fraction (%)')
ax.set_title('Per-Seed Recovery Fractions')
ax.legend(handles=legend_elements, loc='upper left', fontsize=8, ncol=2)
ax.grid(axis='y', alpha=0.3)

# 5b: Timing per seed (select key regions)
ax = axes[1]
key_timing_regions = ["CA-C", "OR", "JDF", "SS-S", "BC-N", "AK-FS", "AK-PWS", "AK-WG", "AK-AL"]

for rnd_idx, rnd in enumerate(ROUNDS):
    for seed in SEEDS:
        vals = []
        for region in key_timing_regions:
            vals.append(timing_data[rnd][region][SEEDS.index(seed)])
        offset = (rnd_idx - 1.5) * 0.08
        ax.scatter(np.arange(len(key_timing_regions)) + offset, vals,
                   color=ROUND_COLORS[rnd], marker=seed_markers[seed], s=25, alpha=0.7,
                   edgecolors='black', linewidths=0.3)

# Target line
target_vals = [TIMING_TARGETS[r] for r in key_timing_regions]
ax.scatter(range(len(key_timing_regions)), target_vals, marker='x', color='black', s=80, linewidths=2, label='Target', zorder=5)

ax.set_xticks(range(len(key_timing_regions)))
ax.set_xticklabels(key_timing_regions, rotation=45, ha='right')
ax.set_ylabel('Arrival Time (months)')
ax.set_title('Per-Seed Arrival Timing')
ax.legend(handles=legend_elements + [Line2D([0], [0], marker='x', color='black', markersize=10, linestyle='None', label='Target')],
          loc='upper left', fontsize=8, ncol=2)
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "fig5_per_seed.pdf"), bbox_inches='tight')
plt.close(fig)
print("Saved fig5_per_seed.pdf")

# ============================================================
# Figure 6: RMSE comparison
# ============================================================
fig, ax = plt.subplots(figsize=(7, 4.5))

x = np.arange(len(ROUNDS))
means = [np.mean(rmse_data[rnd]) for rnd in ROUNDS]
stds = [np.std(rmse_data[rnd]) for rnd in ROUNDS]

bars = ax.bar(x, means, 0.5, color=[ROUND_COLORS[rnd] for rnd in ROUNDS], 
              edgecolor='black', linewidth=0.5, alpha=0.85,
              yerr=stds, capsize=5, error_kw={'linewidth': 1.2})

for i, (m, s) in enumerate(zip(means, stds)):
    ax.text(i, m + s + 0.01, f'{m:.3f}', ha='center', va='bottom', fontsize=10, fontweight='bold')

ax.set_xlabel('Calibration Round')
ax.set_ylabel('RMSE (log-scale recovery)')
ax.set_title('Recovery RMSE Across Wavefront Rounds')
ax.set_xticks(x)
ax.set_xticklabels([f'{rnd}\n(D_P={ROUND_PARAMS[rnd]["D_P"]}km)' for rnd in ROUNDS])
ax.set_ylim(0, 1.4)
ax.grid(axis='y', alpha=0.3)
plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "fig6_rmse.pdf"), bbox_inches='tight')
plt.close(fig)
print("Saved fig6_rmse.pdf")

# ============================================================
# Output data summary for LaTeX tables
# ============================================================
print("\n" + "=" * 70)
print("DATA FOR LATEX TABLES")
print("=" * 70)

print("\n--- Recovery Fractions (mean ± std across seeds, %) ---")
header = f"{'Region':<8} {'Target':>8}"
for rnd in ROUNDS:
    header += f"  {rnd:>12}"
print(header)
for region in RECOVERY_REGIONS:
    row = f"{region:<8} {RECOVERY_TARGETS[region]*100:>8.2f}"
    for rnd in ROUNDS:
        m = np.mean(recovery_data[rnd][region]) * 100
        s = np.std(recovery_data[rnd][region]) * 100
        row += f"  {m:>5.1f}±{s:>4.1f}"
    print(row)

print("\n--- Arrival Timing (mean ± std across seeds, months) ---")
header = f"{'Region':<8} {'Target':>8}"
for rnd in ROUNDS:
    header += f"  {rnd:>12}"
print(header)
for region in TIMING_REGIONS:
    row = f"{region:<8} {TIMING_TARGETS[region]:>8.0f}"
    for rnd in ROUNDS:
        m = np.mean(timing_data[rnd][region])
        s = np.std(timing_data[rnd][region])
        row += f"  {m:>5.1f}±{s:>4.1f}"
    print(row)

print("\n--- RMSE Summary ---")
for rnd in ROUNDS:
    vals = rmse_data[rnd]
    print(f"{rnd}: seeds={[f'{v:.4f}' for v in vals]}, mean={np.mean(vals):.4f}, std={np.std(vals):.4f}")

print("\n--- MAE Timing (months) ---")
for rnd in ROUNDS:
    maes = []
    for seed in SEEDS:
        d = all_data[rnd][seed]
        maes.append(d["arrival_timing"]["mae_months"])
    print(f"{rnd}: seeds={maes}, mean={np.mean(maes):.2f}")

print("\nAll figures saved to:", FIG_DIR)
