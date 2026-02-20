#!/usr/bin/env python3
"""Regenerate all 7 Morris R3 figures from corrected data (v2 - degenerate metrics fixed)."""

import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import cm
import os

RESULTS_DIR = os.path.join(os.path.dirname(__file__), '../../results/sensitivity_r3')
FIG_DIR = os.path.join(RESULTS_DIR, 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

with open(os.path.join(RESULTS_DIR, 'morris_results.json')) as f:
    data = json.load(f)

metrics = list(data.keys())
names = data[metrics[0]]['names']
n_params = len(names)
n_metrics = len(metrics)

# Short param names for display
short_names = [n.split('.')[-1] for n in names]

# Module colors
MODULE_COLORS = {
    'disease': '#e74c3c',
    'population': '#3498db',
    'genetics': '#2ecc71',
    'spawning': '#f39c12',
    'pathogen_evolution': '#9b59b6',
    'spatial': '#1abc9c',
}

def get_module(name):
    parts = name.split('.')
    if len(parts) >= 2:
        mod = parts[0]
        if mod == 'pathogen_evolution':
            return 'pathogen_evolution'
        return mod
    return 'other'

param_colors = [MODULE_COLORS.get(get_module(n), '#95a5a6') for n in names]

# Compute normalized mu* matrix (n_params x n_metrics)
norm_matrix = np.zeros((n_params, n_metrics))
raw_matrix = np.zeros((n_params, n_metrics))
sigma_matrix = np.zeros((n_params, n_metrics))
for j, metric in enumerate(metrics):
    mu_stars = np.array(data[metric]['mu_star'])
    sigmas = np.array(data[metric]['sigma'])
    raw_matrix[:, j] = mu_stars
    sigma_matrix[:, j] = sigmas
    max_ms = max(mu_stars) if max(mu_stars) > 0 else 1
    norm_matrix[:, j] = mu_stars / max_ms

# Global ranking by mean normalized mu*
mean_norm = norm_matrix.mean(axis=1)
global_order = np.argsort(-mean_norm)

# ============================================================
# Figure 1: Global Ranking Bar Chart
# ============================================================
print("Figure 1: Global ranking...")
fig, ax = plt.subplots(figsize=(12, 10))
y_pos = np.arange(n_params)
sorted_vals = mean_norm[global_order]
sorted_names = [short_names[i] for i in global_order]
sorted_colors = [param_colors[i] for i in global_order]

bars = ax.barh(y_pos, sorted_vals, color=sorted_colors, edgecolor='white', linewidth=0.5)
ax.set_yticks(y_pos)
ax.set_yticklabels(sorted_names, fontsize=8)
ax.invert_yaxis()
ax.set_xlabel('Mean Normalized μ* (across 20 metrics)', fontsize=11)
ax.set_title('SA Round 3 — Morris Global Parameter Ranking (v2, 880 runs)', fontsize=13, fontweight='bold')

# Add module legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=c, label=m.replace('_', ' ').title()) 
                   for m, c in MODULE_COLORS.items()]
ax.legend(handles=legend_elements, loc='lower right', fontsize=8, framealpha=0.9)

# Mark threshold
ax.axvline(x=0.033, color='red', linestyle='--', alpha=0.5, label='5% elimination threshold')
ax.text(0.04, n_params - 2, 'All 43 above threshold', fontsize=8, color='red', alpha=0.7)

plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'fig1_global_ranking.png'), dpi=150, bbox_inches='tight')
plt.close()

# ============================================================
# Figure 2: Heatmap (normalized mu*)
# ============================================================
print("Figure 2: Heatmap...")
# Group metrics by category
metric_groups = {
    'Population': ['pop_crash_pct', 'final_pop_frac', 'recovery', 'extinction', 'peak_mortality', 'time_to_nadir', 'total_disease_deaths'],
    'Evolutionary': ['resistance_shift_mean', 'resistance_shift_max', 'va_retention_mean', 'ef1a_shift_mean', 'evolutionary_rescue_index'],
    'Spatial': ['n_extinct_nodes', 'north_south_mortality_gradient', 'fjord_protection_effect'],
    'Pathogen Evo': ['mean_final_virulence', 'virulence_shift'],
    'Demographic': ['disease_death_fraction', 'spawning_participation', 'mean_recruitment_rate'],
}
# Flatten in group order
metric_order = []
group_boundaries = []
for gname, gmetrics in metric_groups.items():
    group_boundaries.append((len(metric_order), gname))
    for m in gmetrics:
        if m in metrics:
            metric_order.append(m)

metric_indices = [metrics.index(m) for m in metric_order]
heatmap_data = norm_matrix[global_order][:, metric_indices]

fig, ax = plt.subplots(figsize=(14, 12))
im = ax.imshow(heatmap_data, aspect='auto', cmap='YlOrRd', vmin=0, vmax=1)

ax.set_xticks(range(len(metric_order)))
ax.set_xticklabels([m.replace('_', '\n') for m in metric_order], fontsize=6, rotation=45, ha='right')
ax.set_yticks(range(n_params))
ax.set_yticklabels([short_names[i] for i in global_order], fontsize=7)

# Group separators
for start, gname in group_boundaries:
    ax.axvline(x=start - 0.5, color='black', linewidth=1.5)
    ax.text(start + 0.3, -1.5, gname, fontsize=7, fontweight='bold', rotation=0)

plt.colorbar(im, ax=ax, label='Normalized μ*', shrink=0.6)
ax.set_title('SA R3 Morris — Parameter × Metric Heatmap (v2)', fontsize=13, fontweight='bold')
plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'fig2_heatmap.png'), dpi=150, bbox_inches='tight')
plt.close()

# ============================================================
# Figure 3: μ* vs σ scatter for key metrics  
# ============================================================
print("Figure 3: Scatter...")
key_metrics = ['pop_crash_pct', 'resistance_shift_mean', 'mean_final_virulence', 'fjord_protection_effect',
               'disease_death_fraction', 'spawning_participation']
fig, axes = plt.subplots(2, 3, figsize=(16, 10))
axes = axes.flatten()

for idx, metric in enumerate(key_metrics):
    ax = axes[idx]
    j = metrics.index(metric)
    mu_stars = raw_matrix[:, j]
    sigmas = sigma_matrix[:, j]
    
    for i in range(n_params):
        ax.scatter(mu_stars[i], sigmas[i], color=param_colors[i], s=30, alpha=0.7, edgecolor='white', linewidth=0.3)
    
    # Label top 5
    top5 = np.argsort(-mu_stars)[:5]
    for i in top5:
        ax.annotate(short_names[i], (mu_stars[i], sigmas[i]), fontsize=6,
                     xytext=(3, 3), textcoords='offset points')
    
    # σ = μ* line
    max_val = max(max(mu_stars), max(sigmas)) * 1.1
    if max_val > 0:
        ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.3, linewidth=0.8)
        ax.text(max_val * 0.6, max_val * 0.55, 'σ = μ*', fontsize=7, alpha=0.4, rotation=40)
    
    ax.set_xlabel('μ*', fontsize=9)
    ax.set_ylabel('σ', fontsize=9)
    ax.set_title(metric.replace('_', ' ').title(), fontsize=10, fontweight='bold')

fig.suptitle('SA R3 Morris — μ* vs σ (v2): Above diagonal = interaction-dominated', fontsize=12, fontweight='bold')
plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'fig3_scatter.png'), dpi=150, bbox_inches='tight')
plt.close()

# ============================================================
# Figure 4: Category-specific top 10
# ============================================================
print("Figure 4: Category top 10...")
categories = {
    'Population Outcomes': ['pop_crash_pct', 'final_pop_frac', 'recovery', 'extinction', 'peak_mortality', 'time_to_nadir', 'total_disease_deaths'],
    'Evolutionary Outcomes': ['resistance_shift_mean', 'resistance_shift_max', 'va_retention_mean', 'ef1a_shift_mean', 'evolutionary_rescue_index'],
    'Spatial Patterns': ['n_extinct_nodes', 'north_south_mortality_gradient', 'fjord_protection_effect'],
    'Pathogen Evolution': ['mean_final_virulence', 'virulence_shift'],
    'Demographic/Spawning': ['disease_death_fraction', 'spawning_participation', 'mean_recruitment_rate'],
}

fig, axes = plt.subplots(2, 3, figsize=(16, 10))
axes = axes.flatten()

for idx, (cat_name, cat_metrics) in enumerate(categories.items()):
    ax = axes[idx]
    cat_indices = [metrics.index(m) for m in cat_metrics if m in metrics]
    cat_norm = norm_matrix[:, cat_indices].mean(axis=1)
    cat_order = np.argsort(-cat_norm)[:10]
    
    y_pos = np.arange(10)
    ax.barh(y_pos, cat_norm[cat_order], color=[param_colors[i] for i in cat_order])
    ax.set_yticks(y_pos)
    ax.set_yticklabels([short_names[i] for i in cat_order], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel('Mean Norm μ*', fontsize=9)
    ax.set_title(cat_name, fontsize=10, fontweight='bold')

# Remove unused subplot
axes[5].set_visible(False)
fig.suptitle('SA R3 Morris — Category-Specific Top 10 (v2)', fontsize=12, fontweight='bold')
plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'fig4_category_top10.png'), dpi=150, bbox_inches='tight')
plt.close()

# ============================================================
# Figure 5: Rank shift comparison (R1 Sobol vs R3 Morris v2)
# ============================================================
print("Figure 5: Rank shift...")
# R1 Sobol final rankings (from memory)
r1_sobol = {
    'mu_I2D_ref': 1, 'susceptibility_multiplier': 2, 'a_exposure': 3, 'sigma_2_eff': 4,
    'n_additive': 5, 'mu_EI1_ref': 6, 'sigma_D': 7, 'K_half': 8, 'T_ref': 9,
    'L_min_repro': 10, 'sigma_1_eff': 11, 'mu_I1I2_ref': 12, 'k_growth': 13,
    'rho_rec': 14, 'alpha_srs': 15, 'senescence_age': 16, 'gamma_fert': 17,
    'F0': 18, 'settler_survival': 19, 'P_env_max': 20, 's_min': 21,
    'D_L': 22, 'T_vbnc': 23
}

# R3 v2 rankings
r3_ranks = {}
for rank, idx in enumerate(global_order):
    sname = short_names[idx]
    r3_ranks[sname] = rank + 1

# Find shared params
shared = [p for p in r1_sobol if p in r3_ranks]

fig, ax = plt.subplots(figsize=(10, 8))
for param in shared:
    r1_rank = r1_sobol[param]
    r3_rank = r3_ranks[param]
    shift = r1_rank - r3_rank  # positive = improved in R3
    color = '#2ecc71' if shift > 3 else '#e74c3c' if shift < -3 else '#95a5a6'
    ax.plot([0, 1], [r1_rank, r3_rank], color=color, linewidth=1.5, alpha=0.7)
    ax.scatter(0, r1_rank, color=color, s=40, zorder=5)
    ax.scatter(1, r3_rank, color=color, s=40, zorder=5)
    ax.text(-0.05, r1_rank, param, fontsize=7, ha='right', va='center')
    ax.text(1.05, r3_rank, param, fontsize=7, ha='left', va='center')

ax.set_xticks([0, 1])
ax.set_xticklabels(['Sobol R1 (23 params)', 'Morris R3v2 (43 params)'], fontsize=11)
ax.set_ylabel('Rank', fontsize=11)
ax.invert_yaxis()
ax.set_title('Parameter Rank Shift: Sobol R1 → Morris R3v2', fontsize=13, fontweight='bold')

# Legend
from matplotlib.lines import Line2D
legend_lines = [
    Line2D([0], [0], color='#2ecc71', linewidth=2, label='Rose >3 ranks'),
    Line2D([0], [0], color='#e74c3c', linewidth=2, label='Fell >3 ranks'),
    Line2D([0], [0], color='#95a5a6', linewidth=2, label='Stable (±3)'),
]
ax.legend(handles=legend_lines, loc='lower right', fontsize=9)
plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'fig5_rank_shift.png'), dpi=150, bbox_inches='tight')
plt.close()

# ============================================================
# Figure 6: Absolute effects for 6 key metrics
# ============================================================
print("Figure 6: Absolute effects...")
abs_metrics = ['pop_crash_pct', 'total_disease_deaths', 'resistance_shift_mean',
               'mean_final_virulence', 'fjord_protection_effect', 'disease_death_fraction']
abs_units = ['% crash', 'deaths', 'Δr̄', 'virulence', 'Δ survival', 'fraction']

fig, axes = plt.subplots(2, 3, figsize=(16, 10))
axes = axes.flatten()

for idx, (metric, unit) in enumerate(zip(abs_metrics, abs_units)):
    ax = axes[idx]
    j = metrics.index(metric)
    mu_stars = raw_matrix[:, j]
    top10 = np.argsort(-mu_stars)[:10]
    
    y_pos = np.arange(10)
    ax.barh(y_pos, mu_stars[top10], color=[param_colors[i] for i in top10])
    ax.set_yticks(y_pos)
    ax.set_yticklabels([short_names[i] for i in top10], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel(f'μ* ({unit})', fontsize=9)
    ax.set_title(metric.replace('_', ' ').title(), fontsize=10, fontweight='bold')

fig.suptitle('SA R3 Morris — Absolute Effects Top 10 (v2)', fontsize=12, fontweight='bold')
plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'fig6_absolute_effects.png'), dpi=150, bbox_inches='tight')
plt.close()

# ============================================================
# Figure 7: Interaction ratio (σ/μ*) heatmap for top 15 params × key metrics
# ============================================================
print("Figure 7: Interaction ratio...")
top15 = global_order[:15]
key_met = ['pop_crash_pct', 'resistance_shift_mean', 'mean_final_virulence', 
           'fjord_protection_effect', 'total_disease_deaths', 'disease_death_fraction',
           'spawning_participation', 'evolutionary_rescue_index']
key_met_idx = [metrics.index(m) for m in key_met if m in metrics]

ratio_data = np.zeros((len(top15), len(key_met_idx)))
for i, pi in enumerate(top15):
    for j_idx, mj in enumerate(key_met_idx):
        ms = raw_matrix[pi, mj]
        sig = sigma_matrix[pi, mj]
        if ms > 1e-10:
            ratio_data[i, j_idx] = sig / ms
        else:
            ratio_data[i, j_idx] = 0

fig, ax = plt.subplots(figsize=(12, 8))
im = ax.imshow(ratio_data, aspect='auto', cmap='RdYlBu_r', vmin=0, vmax=3)
ax.set_xticks(range(len(key_met_idx)))
ax.set_xticklabels([key_met[i].replace('_', '\n') for i in range(len(key_met_idx))], fontsize=7, rotation=45, ha='right')
ax.set_yticks(range(len(top15)))
ax.set_yticklabels([short_names[i] for i in top15], fontsize=8)

# Add text values
for i in range(len(top15)):
    for j in range(len(key_met_idx)):
        val = ratio_data[i, j]
        color = 'white' if val > 2 else 'black'
        ax.text(j, i, f'{val:.1f}', ha='center', va='center', fontsize=6, color=color)

plt.colorbar(im, ax=ax, label='σ/μ* (>1 = interaction-dominated)', shrink=0.7)
ax.set_title('SA R3 Morris — Interaction Ratio σ/μ* for Top 15 Parameters (v2)', fontsize=12, fontweight='bold')
plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'fig7_interaction_ratio.png'), dpi=150, bbox_inches='tight')
plt.close()

print("\nAll 7 figures regenerated successfully!")
print(f"Output: {FIG_DIR}/fig[1-7]_*.png")
