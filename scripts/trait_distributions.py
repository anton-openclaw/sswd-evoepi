#!/usr/bin/env python3
"""Visualize trait distributions and breeding/screening implications.

Questions answered:
1. What do r, t, c distributions look like in the initial population?
2. What fraction of the population exceeds each trait level?
3. Can an individual be maxed on all three traits?
4. How many individuals must you screen to find "good genetics"?

Usage:
    python3 scripts/trait_distributions.py
    
Output: results/trait_distributions/*.png
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from sswd_evoepi.genetics import (
    initialize_trait_effect_sizes,
    initialize_genotypes_three_trait,
    compute_trait_batch,
    trait_slices,
    N_RESISTANCE_DEFAULT,
    N_TOLERANCE_DEFAULT,
    N_RECOVERY_DEFAULT,
)

# ── Parameters ───────────────────────────────────────────────────────
N_POP = 100_000       # Large population to get smooth distributions
SEED = 42
TARGET_R = 0.15
TARGET_T = 0.10
TARGET_C = 0.02        # Current model default

# Disease parameters (for survival calculations)
TAU_MAX = 0.85
RHO_REC = 0.05
MU_I2D_REF = 0.563       # I2→D rate at 20°C (d⁻¹)
MU_I1I2_REF = 0.434      # I1→I2 rate at 20°C
C_EARLY_THRESH = 0.5     # Min c_i for I1 early recovery

# Reference population sizes for ratio context
POP_SIZES = [1_000, 5_000, 10_000, 50_000, 100_000]

# ── Initialize ───────────────────────────────────────────────────────
rng = np.random.default_rng(SEED)

# Effect sizes
effects_r = initialize_trait_effect_sizes(rng, N_RESISTANCE_DEFAULT, 1.0)
effects_t = initialize_trait_effect_sizes(rng, N_TOLERANCE_DEFAULT, 1.0)
effects_c = initialize_trait_effect_sizes(rng, N_RECOVERY_DEFAULT, 1.0)

# Genotypes
geno = initialize_genotypes_three_trait(
    N_POP, effects_r, effects_t, effects_c, rng,
    target_mean_r=TARGET_R,
    target_mean_t=TARGET_T,
    target_mean_c=TARGET_C,
)

# Compute trait scores
alive = np.ones(N_POP, dtype=bool)
res_s, tol_s, rec_s = trait_slices(
    N_RESISTANCE_DEFAULT, N_TOLERANCE_DEFAULT, N_RECOVERY_DEFAULT
)

r_scores = compute_trait_batch(geno, effects_r, alive, res_s)
t_scores = compute_trait_batch(geno, effects_t, alive, tol_s)
c_scores = compute_trait_batch(geno, effects_c, alive, rec_s)


def survival_probability(r, t, c):
    """P(survive if infected) given trait values."""
    eff_rate = MU_I2D_REF * (1 - t * TAU_MAX)
    eff_rate = max(eff_rate, MU_I2D_REF * 0.05)
    mean_i2 = 1 / eff_rate
    p_rec_daily = RHO_REC * c
    p_rec_i2 = 1 - (1 - p_rec_daily) ** mean_i2
    if c > C_EARLY_THRESH:
        p_early = RHO_REC * 2.0 * (c - C_EARLY_THRESH)
        mean_i1 = 1.0 / MU_I1I2_REF
        p_rec_i1 = 1 - (1 - p_early) ** mean_i1
    else:
        p_rec_i1 = 0
    return p_rec_i1 + (1 - p_rec_i1) * p_rec_i2


surv_probs = np.array([survival_probability(r, t, c)
                        for r, t, c in zip(r_scores, t_scores, c_scores)])
breeding_value = r_scores + (1 - r_scores) * surv_probs

# ── Output directory ─────────────────────────────────────────────────
outdir = Path(__file__).resolve().parent.parent / "results" / "trait_distributions"
outdir.mkdir(parents=True, exist_ok=True)

traits = [
    (r_scores, "Resistance (r)", "#2196F3", TARGET_R),
    (t_scores, "Tolerance (t)", "#4CAF50", TARGET_T),
    (c_scores, "Recovery (c)", "#FF9800", TARGET_C),
]

# ══════════════════════════════════════════════════════════════════════
# FIGURE 1: Trait distributions as fraction of population
# ══════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 3, figsize=(16, 5.5))

for ax, (scores, label, color, target) in zip(axes, traits):
    bins = np.linspace(0, max(0.6, np.max(scores) * 1.15), 60)
    counts, edges = np.histogram(scores, bins=bins)
    fractions = counts / N_POP
    centers = (edges[:-1] + edges[1:]) / 2
    width = edges[1] - edges[0]

    ax.bar(centers, fractions * 100, width=width * 0.9, color=color, alpha=0.7,
           edgecolor='white', linewidth=0.3)
    ax.axvline(np.mean(scores), color='black', linestyle='-', linewidth=2,
               label=f'Mean = {np.mean(scores):.3f}')
    ax.axvline(np.percentile(scores, 99), color='red', linestyle='--', linewidth=1.5,
               label=f'99th %ile = {np.percentile(scores, 99):.3f}')
    ax.axvline(np.max(scores), color='darkred', linestyle=':', linewidth=1.5,
               label=f'Max = {np.max(scores):.3f}')
    ax.set_xlabel(label, fontsize=12)
    ax.set_ylabel("% of population" if ax == axes[0] else "")
    ax.set_xlim(0, max(0.6, np.max(scores) * 1.2))
    ax.legend(fontsize=9)
    ax.set_title(f"N = {N_POP:,}", fontsize=10)

fig.suptitle("Trait Distributions — Fraction of Population", fontsize=14, fontweight='bold')
plt.tight_layout()
fig.savefig(outdir / "fig1_trait_distributions.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved fig1_trait_distributions.png")

# ══════════════════════════════════════════════════════════════════════
# FIGURE 2: Cumulative "1 in N" exceedance (ratio view)
# ══════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 3, figsize=(16, 5.5))

for ax, (scores, label, color, target) in zip(axes, traits):
    sorted_scores = np.sort(scores)[::-1]
    # Fraction of population with trait >= threshold
    thresholds = np.linspace(0, np.max(scores) * 1.05, 500)
    fractions = np.array([np.mean(scores >= t) for t in thresholds])

    ax.semilogy(thresholds, fractions, color=color, linewidth=2.5)

    # Reference ratio lines
    for ratio, lbl, alpha in [(1/10, "1 in 10", 0.6), (1/100, "1 in 100", 0.45),
                               (1/1000, "1 in 1,000", 0.3), (1/10000, "1 in 10,000", 0.2)]:
        ax.axhline(ratio, color='gray', linestyle=':', alpha=alpha)
        ax.text(ax.get_xlim()[0] + 0.005 if ax.get_xlim()[0] > 0 else 0.005,
                ratio * 1.15, lbl, fontsize=7, color='gray', alpha=0.8)

    ax.set_xlabel(label, fontsize=12)
    ax.set_ylabel("Fraction of population ≥ threshold" if ax == axes[0] else "")
    ax.set_ylim(1e-5, 1.1)
    ax.grid(True, alpha=0.2)

    # Right y-axis: "1 in N" labels
    ax2 = ax.twinx()
    ax2.set_yscale('log')
    ax2.set_ylim(ax.get_ylim())
    ratio_ticks = [1, 0.1, 0.01, 0.001, 0.0001, 0.00001]
    ratio_labels = ["1:1", "1:10", "1:100", "1:1K", "1:10K", "1:100K"]
    ax2.set_yticks(ratio_ticks)
    ax2.set_yticklabels(ratio_labels, fontsize=8)
    if ax == axes[2]:
        ax2.set_ylabel("Ratio (1 in N)", fontsize=10)

fig.suptitle("How Rare Are High-Trait Individuals? (Fraction of Population)",
             fontsize=14, fontweight='bold')
plt.tight_layout()
fig.savefig(outdir / "fig2_exceedance_curves.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved fig2_exceedance_curves.png")

# ══════════════════════════════════════════════════════════════════════
# FIGURE 3: Screening effort — framed as ratios + absolute numbers
# ══════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(15, 6.5))

# Left: How many to screen (absolute)
ax = axes[0]
thresholds_r = np.arange(0.15, 0.56, 0.005)
p_exceed = np.array([np.mean(r_scores >= t) for t in thresholds_r])

n_needed_95 = np.where(p_exceed > 0, np.ceil(np.log(0.05) / np.log(1 - p_exceed)), np.inf)
n_needed_50 = np.where(p_exceed > 0, np.ceil(np.log(0.50) / np.log(1 - p_exceed)), np.inf)

ax.semilogy(thresholds_r, n_needed_95, color='#D32F2F', linewidth=2.5, label='95% confidence')
ax.semilogy(thresholds_r, n_needed_50, color='#1976D2', linewidth=2.5, label='50% confidence')

key_thresholds = [
    (0.20, "r=0.20\n(modest)"),
    (0.30, "r=0.30\n(good)"),
    (0.40, "r=0.40\n(strong)"),
    (0.50, "r=0.50\n(elite)"),
]
for thresh, label in key_thresholds:
    idx = np.argmin(np.abs(thresholds_r - thresh))
    if n_needed_95[idx] < 1e7:
        ax.annotate(f"{label}\nn≈{int(n_needed_95[idx]):,}",
                    xy=(thresh, n_needed_95[idx]),
                    xytext=(thresh + 0.015, n_needed_95[idx] * 3),
                    fontsize=9, ha='left',
                    arrowprops=dict(arrowstyle='->', color='gray', lw=1))

for n_ref, lbl in [(50, "50 (field trip)"), (200, "200 (survey)"),
                    (1000, "1,000 (major effort)")]:
    ax.axhline(n_ref, color='orange', linestyle=':', alpha=0.4)
    ax.text(0.155, n_ref * 1.15, lbl, fontsize=8, color='orange', alpha=0.7)

ax.set_xlabel("Minimum Resistance Threshold", fontsize=12)
ax.set_ylabel("Number of Individuals to Screen", fontsize=12)
ax.set_title("Screening Effort (absolute)", fontsize=11)
ax.legend(fontsize=11)
ax.set_ylim(1, 1e6)
ax.set_xlim(0.15, 0.55)
ax.grid(True, alpha=0.3)

# Right: Population ratios — what fraction has r >= threshold?
ax = axes[1]
thresholds_fine = np.linspace(0.10, 0.55, 200)
frac_exceed = np.array([np.mean(r_scores >= t) for t in thresholds_fine])

ax.plot(thresholds_fine, frac_exceed * 100, color='#2196F3', linewidth=2.5)
ax.fill_between(thresholds_fine, 0, frac_exceed * 100, color='#2196F3', alpha=0.1)

# Annotate key points
for thresh, pct_label in [(0.15, None), (0.20, "1 in 4"), (0.25, "1 in 9"),
                           (0.30, "1 in 25"), (0.35, "1 in 84"),
                           (0.40, "1 in 345"), (0.50, "1 in 10,000")]:
    frac = np.mean(r_scores >= thresh)
    if frac > 0 and pct_label:
        ax.plot(thresh, frac * 100, 'o', color='#D32F2F', markersize=6, zorder=5)
        ax.annotate(f"r≥{thresh:.2f}: {frac*100:.1f}%\n({pct_label})",
                    xy=(thresh, frac * 100),
                    xytext=(thresh + 0.02, frac * 100 + 3),
                    fontsize=8, ha='left',
                    arrowprops=dict(arrowstyle='->', color='gray', lw=0.8))

ax.set_xlabel("Resistance Threshold", fontsize=12)
ax.set_ylabel("% of population at or above threshold", fontsize=12)
ax.set_title("Population Fraction Exceeding Resistance Threshold", fontsize=11)
ax.set_xlim(0.10, 0.55)
ax.set_ylim(0, 55)
ax.grid(True, alpha=0.3)

fig.suptitle("Resistance Screening: Effort & Population Ratios", fontsize=14, fontweight='bold')
plt.tight_layout()
fig.savefig(outdir / "fig3_screening_effort.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved fig3_screening_effort.png")

# ══════════════════════════════════════════════════════════════════════
# FIGURE 4: Joint trait distribution — can you be good at everything?
# ══════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

pairs = [
    (r_scores, t_scores, "Resistance", "Tolerance"),
    (r_scores, c_scores, "Resistance", "Recovery"),
    (t_scores, c_scores, "Tolerance", "Recovery"),
]

for ax, (x, y, xlabel, ylabel) in zip(axes, pairs):
    idx = rng.choice(N_POP, size=min(5000, N_POP), replace=False)
    ax.scatter(x[idx], y[idx], alpha=0.15, s=3, c='steelblue')
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    corr = np.corrcoef(x, y)[0, 1]
    ax.set_title(f"ρ = {corr:.3f}", fontsize=10)
    combined = x + y
    best = np.argmax(combined)
    ax.scatter([x[best]], [y[best]], color='red', s=80, zorder=5,
               edgecolors='black', linewidth=1.5, label='Best combined')
    ax.legend(fontsize=9)

fig.suptitle("Joint Trait Distributions — Are Traits Correlated?", fontsize=14, fontweight='bold')
plt.tight_layout()
fig.savefig(outdir / "fig4_joint_distributions.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved fig4_joint_distributions.png")

# ══════════════════════════════════════════════════════════════════════
# FIGURE 5: Ratio table visualization — expected count per pop size
# ══════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(14, 8))

# For each threshold, show expected count in populations of different sizes
r_thresholds = [0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]
pop_fracs = [np.mean(r_scores >= t) for t in r_thresholds]

# Create grouped bar chart
x = np.arange(len(r_thresholds))
width = 0.15
colors = plt.cm.plasma(np.linspace(0.15, 0.85, len(POP_SIZES)))

for i, (N, col) in enumerate(zip(POP_SIZES, colors)):
    counts = [frac * N for frac in pop_fracs]
    bars = ax.bar(x + i * width, counts, width, label=f'N={N:,}', color=col, alpha=0.8)
    # Label bars with counts
    for bar, count in zip(bars, counts):
        if count >= 1:
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                    f'{int(count):,}', ha='center', va='bottom', fontsize=6, rotation=45)
        elif count > 0:
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                    f'{count:.1f}', ha='center', va='bottom', fontsize=6, rotation=45)

ax.set_xlabel("Resistance Threshold (r ≥)", fontsize=13)
ax.set_ylabel("Expected Number of Individuals", fontsize=13)
ax.set_title("Expected Count of Resistant Individuals by Population Size",
             fontsize=14, fontweight='bold')
ax.set_xticks(x + width * (len(POP_SIZES) - 1) / 2)
ax.set_xticklabels([f'r ≥ {t:.2f}\n({p*100:.1f}%)' for t, p in zip(r_thresholds, pop_fracs)],
                    fontsize=9)
ax.legend(fontsize=10, title="Population Size")
ax.set_yscale('log')
ax.set_ylim(0.1, N_POP * 1.5)
ax.grid(True, alpha=0.2, axis='y')

# Add ratio annotation at top
for i, (t, f) in enumerate(zip(r_thresholds, pop_fracs)):
    if f > 0:
        ratio = int(1 / f) if f > 0 else "∞"
        ax.text(x[i] + width * (len(POP_SIZES) - 1) / 2, ax.get_ylim()[1] * 0.7,
                f'1 in {ratio:,}', ha='center', fontsize=8, color='gray', style='italic')

plt.tight_layout()
fig.savefig(outdir / "fig5_population_counts.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved fig5_population_counts.png")

# ══════════════════════════════════════════════════════════════════════
# FIGURE 6: Survival distributions + breeding value
# ══════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

ax = axes[0]
bins = np.linspace(0, np.max(surv_probs) * 1.2, 60)
counts, edges = np.histogram(surv_probs * 100, bins=80)
fracs = counts / N_POP * 100
centers = (edges[:-1] + edges[1:]) / 2
ax.bar(centers, fracs, width=centers[1] - centers[0], color='#7B1FA2', alpha=0.7,
       edgecolor='white', linewidth=0.3)
ax.axvline(np.mean(surv_probs) * 100, color='black', linewidth=2,
           label=f'Mean = {np.mean(surv_probs)*100:.2f}%')
ax.axvline(np.percentile(surv_probs, 99) * 100, color='red', linestyle='--', linewidth=1.5,
           label=f'99th %ile = {np.percentile(surv_probs, 99)*100:.2f}%')
ax.axvline(np.max(surv_probs) * 100, color='darkred', linestyle=':', linewidth=1.5,
           label=f'Max = {np.max(surv_probs)*100:.2f}%')
ax.set_xlabel("P(survive | infected) (%)", fontsize=12)
ax.set_ylabel("% of population", fontsize=12)
ax.set_title("Disease Survival Probability\n(tolerance + recovery only)", fontsize=11)
ax.legend(fontsize=10)

ax = axes[1]
counts, edges = np.histogram(breeding_value * 100, bins=80)
fracs = counts / N_POP * 100
centers = (edges[:-1] + edges[1:]) / 2
ax.bar(centers, fracs, width=centers[1] - centers[0], color='#00695C', alpha=0.7,
       edgecolor='white', linewidth=0.3)
ax.axvline(np.mean(breeding_value) * 100, color='black', linewidth=2,
           label=f'Mean = {np.mean(breeding_value)*100:.1f}%')
ax.axvline(np.percentile(breeding_value, 99) * 100, color='red', linestyle='--', linewidth=1.5,
           label=f'99th %ile = {np.percentile(breeding_value, 99)*100:.1f}%')
ax.axvline(np.max(breeding_value) * 100, color='darkred', linestyle=':', linewidth=1.5,
           label=f'Max = {np.max(breeding_value)*100:.1f}%')
ax.set_xlabel("P(survive one exposure event) (%)", fontsize=12)
ax.set_ylabel("% of population", fontsize=12)
ax.set_title("Combined Breeding Value\nr + (1−r) × P(survive|infected)", fontsize=11)
ax.legend(fontsize=10)

fig.suptitle("Disease Outcome Distributions (as % of population)", fontsize=14, fontweight='bold')
plt.tight_layout()
fig.savefig(outdir / "fig6_survival_distributions.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved fig6_survival_distributions.png")

# ══════════════════════════════════════════════════════════════════════
# FIGURE 7: Landrace screening — expected best individual vs sample size
# ══════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
sample_sizes = np.array([10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000])
n_sims = 2000

ax = axes[0]
for scores, label, color in [(r_scores, "Resistance", "#2196F3"),
                               (t_scores, "Tolerance", "#4CAF50"),
                               (c_scores, "Recovery", "#FF9800")]:
    medians, q25, q75 = [], [], []
    for n in sample_sizes:
        max_vals = np.array([np.max(rng.choice(scores, size=n, replace=False))
                             for _ in range(n_sims)])
        medians.append(np.median(max_vals))
        q25.append(np.percentile(max_vals, 25))
        q75.append(np.percentile(max_vals, 75))
    medians, q25, q75 = map(np.array, (medians, q25, q75))
    ax.semilogx(sample_sizes, medians, 'o-', color=color, linewidth=2, label=label, markersize=5)
    ax.fill_between(sample_sizes, q25, q75, color=color, alpha=0.15)

ax.set_xlabel("Number of Individuals Screened", fontsize=12)
ax.set_ylabel("Expected Best Trait Value Found", fontsize=12)
ax.set_title("Expected Best Individual vs Screening Effort", fontsize=11)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)
ax.set_xticks(sample_sizes)
ax.set_xticklabels([str(n) for n in sample_sizes], fontsize=8, rotation=45)

ax = axes[1]
medians_bv, q25_bv, q75_bv = [], [], []
for n in sample_sizes:
    max_bv = np.array([np.max(rng.choice(breeding_value, size=n, replace=False))
                        for _ in range(n_sims)])
    medians_bv.append(np.median(max_bv))
    q25_bv.append(np.percentile(max_bv, 25))
    q75_bv.append(np.percentile(max_bv, 75))
medians_bv, q25_bv, q75_bv = map(np.array, (medians_bv, q25_bv, q75_bv))

ax.semilogx(sample_sizes, medians_bv * 100, 'o-', color='#00695C', linewidth=2.5, markersize=6)
ax.fill_between(sample_sizes, q25_bv * 100, q75_bv * 100, color='#00695C', alpha=0.15)
ax.set_xlabel("Number of Individuals Screened", fontsize=12)
ax.set_ylabel("Expected Best Breeding Value (%)", fontsize=12)
ax.set_title("Expected Best Survival per Exposure vs Screening", fontsize=11)
ax.grid(True, alpha=0.3)
ax.set_xticks(sample_sizes)
ax.set_xticklabels([str(n) for n in sample_sizes], fontsize=8, rotation=45)

fig.suptitle("Landrace Screening: Diminishing Returns", fontsize=14, fontweight='bold')
plt.tight_layout()
fig.savefig(outdir / "fig7_landrace_screening.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved fig7_landrace_screening.png")

# ══════════════════════════════════════════════════════════════════════
# SUMMARY TABLE (console + markdown)
# ══════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("SUMMARY STATISTICS")
print("=" * 70)

for scores, label, target in [(r_scores, "Resistance", TARGET_R),
                                (t_scores, "Tolerance", TARGET_T),
                                (c_scores, "Recovery", TARGET_C)]:
    print(f"\n{label} (target mean={target}):")
    print(f"  Mean: {np.mean(scores):.4f}  |  Std: {np.std(scores):.4f}")
    print(f"  Min:  {np.min(scores):.4f}  |  Max: {np.max(scores):.4f}")
    print(f"  25th: {np.percentile(scores,25):.4f}  |  50th: {np.percentile(scores,50):.4f}  "
          f"|  75th: {np.percentile(scores,75):.4f}  |  95th: {np.percentile(scores,95):.4f}  "
          f"|  99th: {np.percentile(scores,99):.4f}")

print(f"\nSurvival if infected:")
print(f"  Mean: {np.mean(surv_probs)*100:.3f}%  |  99th: {np.percentile(surv_probs,99)*100:.3f}%  "
      f"|  Max: {np.max(surv_probs)*100:.3f}%")

print(f"\nBreeding value:")
print(f"  Mean: {np.mean(breeding_value)*100:.2f}%  |  99th: {np.percentile(breeding_value,99)*100:.2f}%  "
      f"|  Max: {np.max(breeding_value)*100:.2f}%")

# Screening table
print(f"\n{'=' * 70}")
print(f"RESISTANCE SCREENING TABLE")
print(f"{'=' * 70}")
print(f"{'Threshold':>10} {'Fraction':>10} {'Ratio':>12} {'Screen(50%)':>14} {'Screen(95%)':>14}")
print(f"{'─'*10} {'─'*10} {'─'*12} {'─'*14} {'─'*14}")
for thresh in [0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]:
    p = np.mean(r_scores >= thresh)
    if p > 0:
        n95 = int(np.ceil(np.log(0.05) / np.log(1 - p)))
        n50 = int(np.ceil(np.log(0.50) / np.log(1 - p)))
        ratio = f"1 in {int(1/p):,}"
        print(f"  r≥{thresh:.2f}    {p*100:>7.2f}%   {ratio:>12}   {n50:>12,}   {n95:>12,}")
    else:
        print(f"  r≥{thresh:.2f}        0%       —          ∞              ∞")

# Expected counts by population size
print(f"\n{'=' * 70}")
print(f"EXPECTED COUNTS BY POPULATION SIZE")
print(f"{'=' * 70}")
header = f"{'Threshold':>10}" + "".join(f"{'N='+str(n):>12}" for n in POP_SIZES)
print(header)
print("─" * (10 + 12 * len(POP_SIZES)))
for thresh in [0.20, 0.25, 0.30, 0.35, 0.40, 0.50]:
    p = np.mean(r_scores >= thresh)
    row = f"  r≥{thresh:.2f}  "
    for N in POP_SIZES:
        count = p * N
        if count >= 10:
            row += f"  {int(count):>9,}"
        elif count >= 1:
            row += f"  {count:>9.1f}"
        else:
            row += f"  {count:>9.2f}"
    print(row)

# Joint probabilities
print(f"\nJoint trait probabilities:")
combos = [
    ("r≥0.25 & t≥0.20 & c≥0.04", (r_scores >= 0.25) & (t_scores >= 0.20) & (c_scores >= 0.04)),
    ("r≥0.30 & t≥0.15 & c≥0.03", (r_scores >= 0.30) & (t_scores >= 0.15) & (c_scores >= 0.03)),
    ("r≥0.30 (resistance only)", r_scores >= 0.30),
    ("r≥0.20 & t≥0.15", (r_scores >= 0.20) & (t_scores >= 0.15)),
]
for label, mask in combos:
    frac = np.mean(mask)
    if frac > 0:
        print(f"  {label}: {frac*100:.3f}% (1 in {int(1/frac):,})")
    else:
        print(f"  {label}: 0%")

print(f"\nTrait correlations: r-t={np.corrcoef(r_scores,t_scores)[0,1]:.4f}, "
      f"r-c={np.corrcoef(r_scores,c_scores)[0,1]:.4f}, "
      f"t-c={np.corrcoef(t_scores,c_scores)[0,1]:.4f}")

print(f"\nAll figures saved to: {outdir}")
