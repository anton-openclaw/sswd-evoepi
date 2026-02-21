#!/usr/bin/env python3
"""
Morris R4 Analysis: Figures + Report
47 parameters, 23 metrics, 11-node network with pathogen evolution + three-trait genetics
"""

import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from pathlib import Path
import textwrap

# Paths
BASE = Path(__file__).resolve().parents[2]
R4_DIR = BASE / 'results' / 'sensitivity_r4'
R3_DIR = BASE / 'results' / 'sensitivity_r3'
FIG_DIR = R4_DIR / 'figures'
FIG_DIR.mkdir(parents=True, exist_ok=True)

# ── Load R4 ──
with open(R4_DIR / 'morris_results.json') as f:
    r4 = json.load(f)

param_names = r4['param_names']
metric_names = r4['metric_names']
n_params = len(param_names)
n_metrics = len(metric_names)

# Build mu_star and sigma matrices: shape (n_params, n_metrics)
mu_star = np.zeros((n_params, n_metrics))
sigma = np.zeros((n_params, n_metrics))
mu_star_conf = np.zeros((n_params, n_metrics))
y_std = np.zeros(n_metrics)

for j, m in enumerate(metric_names):
    r = r4['results'][m]
    mu_star[:, j] = r['mu_star']
    sigma[:, j] = r['sigma']
    if 'mu_star_conf' in r:
        mu_star_conf[:, j] = r['mu_star_conf']
    if 'y_std' in r:
        y_std[j] = r['y_std']

# Normalize mu_star per metric (divide by metric std or max) for cross-metric comparison
# Use y_std if available, else max of mu_star per metric
norm_factors = np.where(y_std > 0, y_std, mu_star.max(axis=0))
norm_factors = np.where(norm_factors > 0, norm_factors, 1.0)
mu_star_norm = mu_star / norm_factors[np.newaxis, :]

# Mean normalized mu_star across all metrics
mean_mu_star_norm = mu_star_norm.mean(axis=1)

# Ranking
rank_order = np.argsort(-mean_mu_star_norm)
ranked_params = [param_names[i] for i in rank_order]
ranked_values = mean_mu_star_norm[rank_order]

# Interaction ratio: sigma / mu_star (mean across metrics, only where mu_star > threshold)
with np.errstate(divide='ignore', invalid='ignore'):
    ratio = sigma / mu_star
    ratio = np.where(np.isfinite(ratio) & (mu_star > 0.1), ratio, np.nan)
mean_ratio = np.nanmean(ratio, axis=1)

# ── Load R3 for comparison ──
with open(R3_DIR / 'morris_results.json') as f:
    r3 = json.load(f)

r3_metrics = list(r3.keys())
r3_param_names = r3[r3_metrics[0]]['names']
n_r3_params = len(r3_param_names)

# Build R3 mu_star matrix
r3_mu_star = np.zeros((n_r3_params, len(r3_metrics)))
r3_y_std = np.zeros(len(r3_metrics))
for j, m in enumerate(r3_metrics):
    r3_mu_star[:, j] = r3[m]['mu_star']
    # R3 doesn't have y_std, normalize by max
r3_norm_factors = r3_mu_star.max(axis=0)
r3_norm_factors = np.where(r3_norm_factors > 0, r3_norm_factors, 1.0)
r3_mu_star_norm = r3_mu_star / r3_norm_factors[np.newaxis, :]
r3_mean_norm = r3_mu_star_norm.mean(axis=1)
r3_rank_order = np.argsort(-r3_mean_norm)
r3_ranked_params = [r3_param_names[i] for i in r3_rank_order]

# Short param names for display
def short_name(p):
    return p.split('.')[-1]

# ── Figure 1: Top-20 bar chart ──
fig, ax = plt.subplots(figsize=(12, 8))
top_n = min(20, n_params)
colors = plt.cm.viridis(np.linspace(0.2, 0.9, top_n))

# Group by module
module_colors = {
    'disease': '#e74c3c',
    'population': '#3498db',
    'genetics': '#2ecc71',
    'spawning': '#f39c12',
    'spatial': '#9b59b6',
    'pathogen_evolution': '#e67e22',
}

bar_colors = []
for i in range(top_n):
    p = ranked_params[i]
    module = p.split('.')[0]
    bar_colors.append(module_colors.get(module, '#95a5a6'))

bars = ax.barh(range(top_n-1, -1, -1), ranked_values[:top_n], color=bar_colors, edgecolor='white', linewidth=0.5)
ax.set_yticks(range(top_n-1, -1, -1))
ax.set_yticklabels([short_name(p) for p in ranked_params[:top_n]], fontsize=10)
ax.set_xlabel('Mean Normalized μ* (across 23 metrics)', fontsize=12)
ax.set_title('Morris R4: Top 20 Parameters by Global Sensitivity\n(47 params, 23 metrics, 11-node network)', fontsize=13, fontweight='bold')

# Legend for modules
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=c, label=m.replace('_', ' ').title()) for m, c in module_colors.items()]
ax.legend(handles=legend_elements, loc='lower right', fontsize=9, title='Module')

# Add value labels
for i in range(top_n):
    ax.text(ranked_values[i] + 0.005, top_n-1-i, f'{ranked_values[i]:.3f}', va='center', fontsize=8)

ax.set_xlim(0, ranked_values[0] * 1.15)
plt.tight_layout()
plt.savefig(FIG_DIR / 'morris_r4_top20.png', dpi=150)
plt.close()
print("✓ Figure 1: Top-20 bar chart saved")

# ── Figure 2: σ/μ* interaction ratio plot ──
fig, ax = plt.subplots(figsize=(12, 8))
# Plot mu_star_mean vs sigma_mean for each parameter
mu_star_mean = mu_star.mean(axis=1)
sigma_mean = sigma.mean(axis=1)

scatter_colors = [module_colors.get(p.split('.')[0], '#95a5a6') for p in param_names]

ax.scatter(mu_star_mean, sigma_mean, c=scatter_colors, s=80, alpha=0.7, edgecolors='black', linewidth=0.5, zorder=5)

# Label top parameters
top_indices = rank_order[:15]
for i in top_indices:
    ax.annotate(short_name(param_names[i]),
                (mu_star_mean[i], sigma_mean[i]),
                xytext=(5, 5), textcoords='offset points', fontsize=7,
                fontweight='bold' if i in rank_order[:5] else 'normal')

# Reference lines: σ/μ* = 0.5 and 1.0
max_mu = mu_star_mean.max() * 1.1
ax.plot([0, max_mu], [0, 0.5 * max_mu], '--', color='gray', alpha=0.5, label='σ/μ* = 0.5')
ax.plot([0, max_mu], [0, 1.0 * max_mu], '-.', color='gray', alpha=0.5, label='σ/μ* = 1.0')

ax.set_xlabel('Mean μ* (absolute, across metrics)', fontsize=12)
ax.set_ylabel('Mean σ (across metrics)', fontsize=12)
ax.set_title('Morris R4: μ* vs σ — Interaction Detection\n(σ/μ* > 1.0 = strong interactions/nonlinearity)', fontsize=13, fontweight='bold')
ax.legend(handles=legend_elements + [
    plt.Line2D([0], [0], linestyle='--', color='gray', label='σ/μ* = 0.5 (moderate)'),
    plt.Line2D([0], [0], linestyle='-.', color='gray', label='σ/μ* = 1.0 (strong)'),
], loc='upper left', fontsize=8)
plt.tight_layout()
plt.savefig(FIG_DIR / 'morris_r4_interaction.png', dpi=150)
plt.close()
print("✓ Figure 2: Interaction ratio plot saved")

# ── Figure 3: Per-metric heatmap ──
fig, ax = plt.subplots(figsize=(18, 14))

# Use normalized mu_star for heatmap
# Sort params by mean ranking
sorted_mu_star_norm = mu_star_norm[rank_order, :]

# Truncate labels
short_params = [short_name(p) for p in ranked_params]
short_metrics = [m.replace('_', ' ') for m in metric_names]

im = ax.imshow(sorted_mu_star_norm, aspect='auto', cmap='YlOrRd', interpolation='nearest')
ax.set_xticks(range(n_metrics))
ax.set_xticklabels(short_metrics, rotation=45, ha='right', fontsize=8)
ax.set_yticks(range(n_params))
ax.set_yticklabels(short_params, fontsize=8)
ax.set_xlabel('Metric', fontsize=11)
ax.set_ylabel('Parameter (ranked by global μ*)', fontsize=11)
ax.set_title('Morris R4: Normalized μ* Heatmap\n(Parameter × Metric sensitivity)', fontsize=13, fontweight='bold')

cbar = plt.colorbar(im, ax=ax, shrink=0.8)
cbar.set_label('Normalized μ*', fontsize=10)

# Add horizontal lines to separate top-10
ax.axhline(9.5, color='white', linewidth=2)
ax.text(-1.5, 4.5, 'Top 10', fontsize=9, fontweight='bold', color='red', va='center', ha='right')

plt.tight_layout()
plt.savefig(FIG_DIR / 'morris_r4_heatmap.png', dpi=150)
plt.close()
print("✓ Figure 3: Heatmap saved")

# ── Figure 4: R3 vs R4 comparison ──
# Find common parameters
r3_rank_dict = {p: i+1 for i, p in enumerate(r3_ranked_params)}
r4_rank_dict = {p: i+1 for i, p in enumerate(ranked_params)}

# Common params
common_params = [p for p in param_names if p in r3_rank_dict]
new_params = [p for p in param_names if p not in r3_rank_dict]

fig, ax = plt.subplots(figsize=(10, 12))

# Plot rank changes for common params
r3_ranks = [r3_rank_dict[p] for p in common_params]
r4_ranks = [r4_rank_dict[p] for p in common_params]

# Sort by R4 rank
sort_idx = np.argsort(r4_ranks)
common_sorted = [common_params[i] for i in sort_idx]
r3_sorted = [r3_ranks[i] for i in sort_idx]
r4_sorted = [r4_ranks[i] for i in sort_idx]

y_pos = range(len(common_sorted))
for i, (p, r3r, r4r) in enumerate(zip(common_sorted, r3_sorted, r4_sorted)):
    delta = r3r - r4r  # positive = moved up (better rank in R4)
    color = '#2ecc71' if delta > 3 else ('#e74c3c' if delta < -3 else '#95a5a6')
    arrow = '↑' if delta > 0 else ('↓' if delta < 0 else '→')
    ax.plot([r3r, r4r], [i, i], '-', color=color, linewidth=2, alpha=0.7)
    ax.plot(r3r, i, 'o', color='#3498db', markersize=6, zorder=5)
    ax.plot(r4r, i, 's', color='#e74c3c', markersize=6, zorder=5)
    label = f'{short_name(p)} ({arrow}{abs(delta)})'
    ax.text(max(r3r, r4r) + 0.5, i, label, fontsize=7, va='center')

ax.set_yticks([])
ax.set_xlabel('Rank (lower = more important)', fontsize=12)
ax.set_title('Morris R3→R4: Parameter Rank Changes\n(43 common params; R4 adds 4 new: n_tolerance, target_mean_t, target_mean_c, tau_max)',
             fontsize=11, fontweight='bold')
ax.invert_xaxis()

# Legend
from matplotlib.lines import Line2D
legend_elems = [
    Line2D([0], [0], marker='o', color='#3498db', linestyle='', markersize=8, label='R3 rank'),
    Line2D([0], [0], marker='s', color='#e74c3c', linestyle='', markersize=8, label='R4 rank'),
    Line2D([0], [0], color='#2ecc71', linewidth=2, label='Moved up >3'),
    Line2D([0], [0], color='#e74c3c', linewidth=2, label='Moved down >3'),
    Line2D([0], [0], color='#95a5a6', linewidth=2, label='Stable (±3)'),
]
ax.legend(handles=legend_elems, loc='lower left', fontsize=9)

plt.tight_layout()
plt.savefig(FIG_DIR / 'morris_r4_vs_r3.png', dpi=150)
plt.close()
print("✓ Figure 4: R3 vs R4 comparison saved")

# ── Figure 5: Module-level sensitivity breakdown ──
fig, ax = plt.subplots(figsize=(10, 6))
modules = {}
for i, p in enumerate(param_names):
    mod = p.split('.')[0]
    if mod not in modules:
        modules[mod] = []
    modules[mod].append(mean_mu_star_norm[i])

mod_names = sorted(modules.keys(), key=lambda m: -np.mean(modules[m]))
mod_means = [np.mean(modules[m]) for m in mod_names]
mod_maxes = [np.max(modules[m]) for m in mod_names]
mod_counts = [len(modules[m]) for m in mod_names]

x = range(len(mod_names))
bars = ax.bar(x, mod_means, color=[module_colors.get(m, '#95a5a6') for m in mod_names],
              edgecolor='white', linewidth=0.5)
ax.bar(x, mod_maxes, color='none', edgecolor='black', linewidth=1.5, linestyle='--')

ax.set_xticks(x)
ax.set_xticklabels([m.replace('_', '\n').title() for m in mod_names], fontsize=10)
ax.set_ylabel('Normalized μ*', fontsize=12)
ax.set_title('Morris R4: Module-Level Sensitivity\n(solid = mean, dashed = max parameter)', fontsize=13, fontweight='bold')

for i, (mean, mx, cnt) in enumerate(zip(mod_means, mod_maxes, mod_counts)):
    ax.text(i, mx + 0.01, f'n={cnt}', ha='center', fontsize=9)

plt.tight_layout()
plt.savefig(FIG_DIR / 'morris_r4_modules.png', dpi=150)
plt.close()
print("✓ Figure 5: Module breakdown saved")

# ── Figure 6: New R4 params focus (three-trait genetics) ──
# R4 new params: n_tolerance, target_mean_t, target_mean_c, tau_max
# Also n_resistance replaced n_additive
r4_new_params = ['genetics.n_tolerance', 'genetics.target_mean_t', 'genetics.target_mean_c', 'genetics.tau_max']
# Also genetics.n_resistance (renamed from n_additive), genetics.target_mean_r (existed before)

trait_params = [p for p in param_names if 'genetics' in p]
trait_indices = [param_names.index(p) for p in trait_params]

fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# Left: bar chart of genetic params
trait_vals = mean_mu_star_norm[trait_indices]
sort_idx = np.argsort(-trait_vals)
axes[0].barh(range(len(trait_params)-1, -1, -1),
             trait_vals[sort_idx],
             color='#2ecc71', edgecolor='white')
axes[0].set_yticks(range(len(trait_params)-1, -1, -1))
axes[0].set_yticklabels([short_name(trait_params[i]) for i in sort_idx], fontsize=10)
axes[0].set_xlabel('Mean Normalized μ*', fontsize=11)
axes[0].set_title('Genetic Architecture Parameters', fontsize=12, fontweight='bold')

# Highlight new R4 params
for j, i in enumerate(sort_idx):
    p = trait_params[i]
    if p in r4_new_params:
        axes[0].get_children()[len(trait_params)-1 - j].set_facecolor('#e74c3c')
        axes[0].get_children()[len(trait_params)-1 - j].set_edgecolor('black')

from matplotlib.patches import Patch
axes[0].legend(handles=[
    Patch(facecolor='#2ecc71', label='Existing (R3)'),
    Patch(facecolor='#e74c3c', label='New in R4'),
], loc='lower right', fontsize=9)

# Right: heatmap of genetic params across key metrics
key_metrics = ['resistance_shift_mean', 'tolerance_shift_mean', 'recovery_shift_mean',
               'evolutionary_rescue_index', 'pop_crash_pct', 'recovery_rate']
key_metric_idx = [metric_names.index(m) for m in key_metrics if m in metric_names]
key_metric_names = [m for m in key_metrics if m in metric_names]

trait_heatmap = mu_star_norm[np.ix_(trait_indices, key_metric_idx)]
im = axes[1].imshow(trait_heatmap, aspect='auto', cmap='YlOrRd')
axes[1].set_xticks(range(len(key_metric_names)))
axes[1].set_xticklabels([m.replace('_', '\n') for m in key_metric_names], fontsize=8, rotation=45, ha='right')
axes[1].set_yticks(range(len(trait_params)))
axes[1].set_yticklabels([short_name(p) for p in trait_params], fontsize=9)
axes[1].set_title('Genetic Params × Key Metrics', fontsize=12, fontweight='bold')
plt.colorbar(im, ax=axes[1], shrink=0.8, label='Normalized μ*')

plt.suptitle('Morris R4: Three-Trait Genetic Architecture Analysis', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(FIG_DIR / 'morris_r4_genetics.png', dpi=150, bbox_inches='tight')
plt.close()
print("✓ Figure 6: Genetics focus saved")

# ── Figure 7: Pathogen evolution params ──
pe_params = [p for p in param_names if 'pathogen_evolution' in p]
pe_indices = [param_names.index(p) for p in pe_params]

fig, ax = plt.subplots(figsize=(10, 5))
pe_vals = mean_mu_star_norm[pe_indices]
pe_ratios = mean_ratio[pe_indices]

x = range(len(pe_params))
bars = ax.bar(x, pe_vals, color='#e67e22', edgecolor='white', alpha=0.8, label='Mean norm μ*')
ax2 = ax.twinx()
ax2.plot(x, pe_ratios, 'ko-', markersize=8, label='σ/μ* ratio')
ax2.axhline(1.0, color='gray', linestyle='--', alpha=0.5)

ax.set_xticks(x)
ax.set_xticklabels([short_name(p) for p in pe_params], fontsize=10, rotation=30)
ax.set_ylabel('Mean Normalized μ*', fontsize=11)
ax2.set_ylabel('σ/μ* Interaction Ratio', fontsize=11)
ax.set_title('Morris R4: Pathogen Evolution Parameters\n(new in R4 — virulence trade-offs)', fontsize=13, fontweight='bold')

lines1, labels1 = ax.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax.legend(lines1 + lines2, labels1 + labels2, loc='upper right', fontsize=9)

plt.tight_layout()
plt.savefig(FIG_DIR / 'morris_r4_pathogen_evolution.png', dpi=150)
plt.close()
print("✓ Figure 7: Pathogen evolution params saved")


# ══════════════════════════════════════════════════════════
# Generate Report Data
# ══════════════════════════════════════════════════════════

print("\n" + "="*60)
print("SUMMARY FOR REPORT")
print("="*60)

print(f"\nR4: {n_params} params, {n_metrics} metrics, {r4['total_runs']} runs")
print(f"R3: {n_r3_params} params, {len(r3_metrics)} metrics")
print(f"New params in R4: {new_params}")

# R3→R4 mapping for genetics.n_additive → genetics.n_resistance
print("\nParameter name changes: genetics.n_additive → genetics.n_resistance")

print("\n── Top 47 Parameters Ranked ──")
for i, (p, v) in enumerate(zip(ranked_params, ranked_values)):
    r3r = r3_rank_dict.get(p, '—')
    if p == 'genetics.n_resistance' and 'genetics.n_additive' in r3_rank_dict:
        r3r = r3_rank_dict['genetics.n_additive']
    ratio_val = mean_ratio[param_names.index(p)]
    print(f"  {i+1:>2}. {short_name(p):<30s} μ*={v:.4f}  σ/μ*={ratio_val:.2f}  R3={r3r}")

print("\n── Biggest Rank Changes (R3→R4) ──")
changes = []
for p in common_params:
    r3r = r3_rank_dict[p]
    r4r = r4_rank_dict[p]
    changes.append((p, r3r, r4r, r3r - r4r))
# Also handle n_additive → n_resistance
if 'genetics.n_additive' in r3_rank_dict and 'genetics.n_resistance' in r4_rank_dict:
    r3r = r3_rank_dict['genetics.n_additive']
    r4r = r4_rank_dict['genetics.n_resistance']
    changes.append(('genetics.n_resistance (was n_additive)', r3r, r4r, r3r - r4r))

changes.sort(key=lambda x: -abs(x[3]))
print("  Movers (positive = improved rank):")
for p, r3r, r4r, delta in changes[:15]:
    arrow = '↑' if delta > 0 else '↓'
    print(f"    {short_name(p):<30s} R3=#{r3r} → R4=#{r4r}  ({arrow}{abs(delta)})")

print("\n── High Interaction Params (σ/μ* > 1.0) ──")
for i in rank_order:
    r = mean_ratio[i]
    if r > 1.0:
        print(f"    {short_name(param_names[i]):<30s} σ/μ*={r:.2f}")

print("\n── Module Summary ──")
for m in mod_names:
    vals = modules[m]
    print(f"  {m:<25s} n={len(vals):>2}  mean={np.mean(vals):.4f}  max={np.max(vals):.4f}")

# Save full ranking as JSON for the report
ranking_data = []
for i, idx in enumerate(rank_order):
    p = param_names[idx]
    r3r = r3_rank_dict.get(p, None)
    if p == 'genetics.n_resistance' and r3r is None and 'genetics.n_additive' in r3_rank_dict:
        r3r = r3_rank_dict['genetics.n_additive']
    ranking_data.append({
        'rank': i + 1,
        'parameter': p,
        'short_name': short_name(p),
        'module': p.split('.')[0],
        'mean_norm_mu_star': float(ranked_values[i]),
        'sigma_mu_star_ratio': float(mean_ratio[idx]),
        'r3_rank': r3r,
        'rank_change': (r3r - (i+1)) if r3r else None,
    })

with open(R4_DIR / 'morris_r4_ranking.json', 'w') as f:
    json.dump(ranking_data, f, indent=2)

print(f"\nAll figures saved to {FIG_DIR}")
print(f"Ranking data saved to {R4_DIR / 'morris_r4_ranking.json'}")
