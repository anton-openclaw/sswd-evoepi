#!/usr/bin/env python3
"""Breeding & crossing analysis for Pycnopodia conservation genetics.

Simulates:
1. Post-epidemic trait distributions (using evolved means from validation)
2. Complementary crossing — pairing individuals with different protective loci
3. Multi-generation selective breeding to create "super resistant" offspring
4. How many founders needed and generations required

Usage:
    python3 scripts/breeding_analysis.py

Output: results/breeding_analysis/*.png
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import json

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from sswd_evoepi.genetics import (
    initialize_trait_effect_sizes,
    initialize_genotypes_three_trait,
    compute_trait_batch,
    trait_slices,
    N_RESISTANCE_DEFAULT,
    N_TOLERANCE_DEFAULT,
    N_RECOVERY_DEFAULT,
    N_LOCI,
)

# ── Parameters ───────────────────────────────────────────────────────
SEED = 42
N_LOCI_R = N_RESISTANCE_DEFAULT  # 17

# Post-epidemic trait means from validation (satellite SST, R→S, Prentice)
# Using Howe Sound as representative surviving population (largest survivor group)
POST_EPIDEMIC = {
    'Sitka':      {'N': 10,  'r': 0.1757, 't': 0.0890, 'c': 0.0108},
    'Howe Sound': {'N': 133, 'r': 0.1787, 't': 0.1366, 'c': 0.0278},
    'SJI':        {'N': 3,   'r': 0.2487, 't': 0.1266, 'c': 0.0000},
}
# Pre-epidemic baseline
PRE_EPIDEMIC = {'r': 0.15, 't': 0.10, 'c': 0.02}

# Breeding parameters
N_FOUNDERS = 10_000      # Large pool to characterize distributions
N_CROSSES = 5_000        # Crosses to simulate per generation
N_OFFSPRING_PER_CROSS = 10  # Offspring per cross
N_GENERATIONS = 8         # Breeding generations
TOP_K_SELECT = 50         # Keep top K individuals per generation

# Disease parameters
TAU_MAX = 0.85
RHO_REC = 0.05
MU_I2D_REF = 0.563
MU_I1I2_REF = 0.434
C_EARLY_THRESH = 0.5

# ── Initialize ───────────────────────────────────────────────────────
rng = np.random.default_rng(SEED)

effects_r = initialize_trait_effect_sizes(rng, N_RESISTANCE_DEFAULT, 1.0)
effects_t = initialize_trait_effect_sizes(rng, N_TOLERANCE_DEFAULT, 1.0)
effects_c = initialize_trait_effect_sizes(rng, N_RECOVERY_DEFAULT, 1.0)

res_s, tol_s, rec_s = trait_slices(N_RESISTANCE_DEFAULT, N_TOLERANCE_DEFAULT, N_RECOVERY_DEFAULT)


def make_population(n, target_r, target_t, target_c, rng_local):
    """Create a population with given trait means."""
    geno = initialize_genotypes_three_trait(
        n, effects_r, effects_t, effects_c, rng_local,
        target_mean_r=target_r, target_mean_t=target_t, target_mean_c=target_c,
    )
    alive = np.ones(n, dtype=bool)
    r = compute_trait_batch(geno, effects_r, alive, res_s)
    t = compute_trait_batch(geno, effects_t, alive, tol_s)
    c = compute_trait_batch(geno, effects_c, alive, rec_s)
    return geno, r, t, c


def cross(parent1_geno, parent2_geno, n_offspring, rng_local):
    """Mendelian crossing: for each locus, offspring gets one allele from each parent."""
    offspring = np.zeros((n_offspring, N_LOCI, 2), dtype=np.int8)
    for locus in range(N_LOCI):
        # From parent 1: randomly pick allele 0 or 1
        picks1 = rng_local.integers(0, 2, size=n_offspring)
        offspring[:, locus, 0] = parent1_geno[locus, :][picks1]
        # From parent 2: randomly pick allele 0 or 1
        picks2 = rng_local.integers(0, 2, size=n_offspring)
        offspring[:, locus, 1] = parent2_geno[locus, :][picks2]
    return offspring


def compute_r_from_geno(geno_batch):
    """Compute resistance scores for a batch of genotypes."""
    n = geno_batch.shape[0]
    alive = np.ones(n, dtype=bool)
    return compute_trait_batch(geno_batch, effects_r, alive, res_s)


def locus_profile(geno, trait_slice):
    """Get per-locus protective allele count (0, 1, or 2) for one individual."""
    return geno[trait_slice, :].sum(axis=1)


def complementarity_score(geno1, geno2, trait_slice):
    """How complementary are two individuals? 
    Score = number of loci where at least one parent has a protective allele
    but the other has it at a different locus."""
    p1 = locus_profile(geno1, trait_slice)
    p2 = locus_profile(geno2, trait_slice)
    # Union: loci where at least one parent has ≥1 protective allele
    union = np.sum((p1 > 0) | (p2 > 0))
    # Overlap: loci where both have protective alleles
    overlap = np.sum((p1 > 0) & (p2 > 0))
    return union, overlap, union - overlap  # unique contribution


def survival_probability(r, t, c):
    """P(survive if infected)."""
    eff_rate = MU_I2D_REF * (1 - t * TAU_MAX)
    eff_rate = max(eff_rate, MU_I2D_REF * 0.05)
    mean_i2 = 1 / eff_rate
    p_rec_daily = RHO_REC * c
    p_rec_i2 = 1 - (1 - p_rec_daily) ** mean_i2
    if c > C_EARLY_THRESH:
        p_early = RHO_REC * 2.0 * (c - C_EARLY_THRESH)
        p_rec_i1 = 1 - (1 - p_early) ** (1.0 / MU_I1I2_REF)
    else:
        p_rec_i1 = 0
    return p_rec_i1 + (1 - p_rec_i1) * p_rec_i2


outdir = Path(__file__).resolve().parent.parent / "results" / "breeding_analysis"
outdir.mkdir(parents=True, exist_ok=True)

# ══════════════════════════════════════════════════════════════════════
# PART 1: Pre- vs Post-epidemic distributions comparison
# ══════════════════════════════════════════════════════════════════════
print("Generating pre- and post-epidemic populations...")

rng_pre = np.random.default_rng(SEED)
geno_pre, r_pre, t_pre, c_pre = make_population(
    N_FOUNDERS, PRE_EPIDEMIC['r'], PRE_EPIDEMIC['t'], PRE_EPIDEMIC['c'], rng_pre)

# Post-epidemic: use Howe Sound means (largest surviving pop)
hs = POST_EPIDEMIC['Howe Sound']
rng_post = np.random.default_rng(SEED + 1)
geno_post, r_post, t_post, c_post = make_population(
    N_FOUNDERS, hs['r'], hs['t'], hs['c'], rng_post)

fig, axes = plt.subplots(1, 3, figsize=(16, 5.5))
for ax, (pre, post, label, color) in zip(axes, [
    (r_pre, r_post, "Resistance (r)", "#2196F3"),
    (t_pre, t_post, "Tolerance (t)", "#4CAF50"),
    (c_pre, c_post, "Recovery (c)", "#FF9800"),
]):
    bins = np.linspace(0, max(np.max(pre), np.max(post)) * 1.2, 60)
    ax.hist(pre, bins=bins, density=True, alpha=0.5, color='gray', label=f'Pre-epidemic (mean={np.mean(pre):.3f})')
    ax.hist(post, bins=bins, density=True, alpha=0.6, color=color, label=f'Post-epidemic (mean={np.mean(post):.3f})')
    ax.set_xlabel(label, fontsize=12)
    ax.set_ylabel("Density" if ax == axes[0] else "")
    ax.legend(fontsize=9)

fig.suptitle("Pre- vs Post-Epidemic Trait Distributions\n(Post = Howe Sound evolved means after 20yr)",
             fontsize=14, fontweight='bold')
plt.tight_layout()
fig.savefig(outdir / "fig1_pre_vs_post.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved fig1_pre_vs_post.png")

# ══════════════════════════════════════════════════════════════════════
# PART 2: Locus-level analysis — where are the protective alleles?
# ══════════════════════════════════════════════════════════════════════
print("Analyzing locus-level allele distributions...")

# Get top 1% most resistant in post-epidemic population
top_1pct_idx = np.argsort(r_post)[-int(N_FOUNDERS * 0.01):]
top_5pct_idx = np.argsort(r_post)[-int(N_FOUNDERS * 0.05):]

# Per-locus protective allele frequency in different groups
def allele_freq_per_locus(geno_batch, indices, trait_slice):
    """Frequency of protective allele at each locus."""
    subset = geno_batch[indices]
    return subset[:, trait_slice, :].mean(axis=(0, 2))

all_idx = np.arange(N_FOUNDERS)
freq_all = allele_freq_per_locus(geno_post, all_idx, res_s)
freq_top5 = allele_freq_per_locus(geno_post, top_5pct_idx, res_s)
freq_top1 = allele_freq_per_locus(geno_post, top_1pct_idx, res_s)

fig, axes = plt.subplots(2, 1, figsize=(14, 9))

# Top: allele frequencies by locus
ax = axes[0]
x = np.arange(N_LOCI_R)
width = 0.25
ax.bar(x - width, freq_all, width, label='All individuals', color='gray', alpha=0.6)
ax.bar(x, freq_top5, width, label='Top 5% resistant', color='#42A5F5', alpha=0.7)
ax.bar(x + width, freq_top1, width, label='Top 1% resistant', color='#1565C0', alpha=0.8)
ax.set_xlabel("Resistance Locus Index", fontsize=12)
ax.set_ylabel("Protective Allele Frequency", fontsize=12)
ax.set_title("Protective Allele Frequency by Locus", fontsize=13, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels([f'L{i}' for i in range(N_LOCI_R)], fontsize=8)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.2, axis='y')

# Bottom: effect sizes overlaid
ax2 = axes[1]
ax2.bar(x, effects_r, color='#FF7043', alpha=0.7, label='Effect size')
ax2.set_xlabel("Resistance Locus Index", fontsize=12)
ax2.set_ylabel("Effect Size (contribution to r)", fontsize=12)
ax2.set_title("Per-Locus Effect Sizes (sorted descending)", fontsize=13, fontweight='bold')
ax2.set_xticks(x)
ax2.set_xticklabels([f'L{i}' for i in range(N_LOCI_R)], fontsize=8)
ax2.legend(fontsize=10)
ax2.grid(True, alpha=0.2, axis='y')

# Annotate cumulative
cumsum = np.cumsum(effects_r)
for i in [0, 4, 9, 16]:
    ax2.annotate(f'{cumsum[i]:.1%} cum.', xy=(i, effects_r[i]),
                xytext=(i + 0.5, effects_r[i] + 0.01), fontsize=8, color='gray')

plt.tight_layout()
fig.savefig(outdir / "fig2_locus_analysis.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved fig2_locus_analysis.png")

# ══════════════════════════════════════════════════════════════════════
# PART 3: Complementarity analysis — finding good crossing pairs
# ══════════════════════════════════════════════════════════════════════
print("Analyzing complementarity of top individuals...")

# Take top 50 most resistant
top50_idx = np.argsort(r_post)[-50:]
top50_r = r_post[top50_idx]
top50_geno = geno_post[top50_idx]

# Compute complementarity matrix
n_top = len(top50_idx)
union_matrix = np.zeros((n_top, n_top))
unique_matrix = np.zeros((n_top, n_top))

for i in range(n_top):
    for j in range(i+1, n_top):
        union, overlap, unique = complementarity_score(
            top50_geno[i], top50_geno[j], res_s)
        union_matrix[i, j] = union
        union_matrix[j, i] = union
        unique_matrix[i, j] = unique
        unique_matrix[j, i] = unique

fig, axes = plt.subplots(1, 2, figsize=(15, 6))

ax = axes[0]
im = ax.imshow(union_matrix, cmap='YlOrRd', aspect='auto')
ax.set_title("Locus Union (total covered loci)", fontsize=12, fontweight='bold')
ax.set_xlabel("Parent index (ranked by r)")
ax.set_ylabel("Parent index (ranked by r)")
plt.colorbar(im, ax=ax, label="# loci with ≥1 protective allele")

ax = axes[1]
im = ax.imshow(unique_matrix, cmap='YlGnBu', aspect='auto')
ax.set_title("Unique Contribution (non-overlapping loci)", fontsize=12, fontweight='bold')
ax.set_xlabel("Parent index (ranked by r)")
ax.set_ylabel("Parent index (ranked by r)")
plt.colorbar(im, ax=ax, label="# loci unique to one parent")

fig.suptitle("Complementarity of Top 50 Resistant Individuals\n(Higher = better crossing pair)",
             fontsize=14, fontweight='bold')
plt.tight_layout()
fig.savefig(outdir / "fig3_complementarity.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved fig3_complementarity.png")

# ══════════════════════════════════════════════════════════════════════
# PART 4: Simulated crosses — what do offspring look like?
# ══════════════════════════════════════════════════════════════════════
print("Simulating crosses between top individuals...")

# Strategy 1: Cross two most resistant individuals
idx_best = np.argmax(r_post)
idx_2nd = np.argsort(r_post)[-2]

# Strategy 2: Cross most complementary pair from top 20
top20_idx = np.argsort(r_post)[-20:]
best_comp_score = -1
best_pair = (0, 0)
for i in range(20):
    for j in range(i+1, 20):
        _, _, unique = complementarity_score(
            geno_post[top20_idx[i]], geno_post[top20_idx[j]], res_s)
        if unique > best_comp_score:
            best_comp_score = unique
            best_pair = (top20_idx[i], top20_idx[j])

strategies = {
    'Top 2 by r': (idx_best, idx_2nd),
    'Most complementary\n(from top 20)': best_pair,
}

# Also: random pair from top 5%
top5_pool = np.argsort(r_post)[-int(N_FOUNDERS * 0.05):]

fig, axes = plt.subplots(1, 3, figsize=(16, 5.5))

# For each strategy, simulate many crosses
for ax, (strat_name, (p1, p2)) in zip(axes[:2], strategies.items()):
    offspring_geno = cross(geno_post[p1], geno_post[p2], 1000, rng)
    offspring_r = compute_r_from_geno(offspring_geno)
    
    parent_r1 = r_post[p1]
    parent_r2 = r_post[p2]
    
    ax.hist(offspring_r, bins=40, density=True, color='#7B1FA2', alpha=0.7, edgecolor='white')
    ax.axvline(parent_r1, color='red', linewidth=2, linestyle='--', label=f'Parent 1: r={parent_r1:.3f}')
    ax.axvline(parent_r2, color='blue', linewidth=2, linestyle='--', label=f'Parent 2: r={parent_r2:.3f}')
    ax.axvline(np.mean(offspring_r), color='black', linewidth=2, label=f'Offspring mean: r={np.mean(offspring_r):.3f}')
    ax.axvline(np.max(offspring_r), color='darkred', linewidth=1.5, linestyle=':', 
               label=f'Best offspring: r={np.max(offspring_r):.3f}')
    ax.set_xlabel("Resistance (r)", fontsize=12)
    ax.set_ylabel("Density", fontsize=12)
    ax.set_title(strat_name, fontsize=11, fontweight='bold')
    ax.legend(fontsize=8)

    # Complementarity info
    u, o, uniq = complementarity_score(geno_post[p1], geno_post[p2], res_s)
    ax.text(0.02, 0.98, f"Locus coverage: {u}/17\nUnique loci: {uniq}",
            transform=ax.transAxes, fontsize=9, va='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

# Third panel: random pairs from top 5% (distribution of offspring maxima)
ax = axes[2]
best_offspring_r = []
for _ in range(500):
    i, j = rng.choice(len(top5_pool), size=2, replace=False)
    offspring_geno = cross(geno_post[top5_pool[i]], geno_post[top5_pool[j]], 
                          N_OFFSPRING_PER_CROSS, rng)
    offspring_r = compute_r_from_geno(offspring_geno)
    best_offspring_r.append(np.max(offspring_r))

best_offspring_r = np.array(best_offspring_r)
ax.hist(best_offspring_r, bins=40, density=True, color='#00897B', alpha=0.7, edgecolor='white')
ax.axvline(np.mean(best_offspring_r), color='black', linewidth=2, 
           label=f'Mean best: r={np.mean(best_offspring_r):.3f}')
ax.axvline(np.percentile(best_offspring_r, 95), color='red', linestyle='--', linewidth=1.5,
           label=f'95th %ile: r={np.percentile(best_offspring_r, 95):.3f}')
ax.set_xlabel("Best Offspring Resistance (r)", fontsize=12)
ax.set_ylabel("Density", fontsize=12)
ax.set_title("Random pairs from top 5%\n(best of 10 offspring each)", fontsize=11, fontweight='bold')
ax.legend(fontsize=9)

fig.suptitle("Offspring Distributions from Single Crosses", fontsize=14, fontweight='bold')
plt.tight_layout()
fig.savefig(outdir / "fig4_cross_offspring.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved fig4_cross_offspring.png")

# ══════════════════════════════════════════════════════════════════════
# PART 5: Multi-generation selective breeding
# ══════════════════════════════════════════════════════════════════════
print("Simulating multi-generation selective breeding...")

def run_breeding_program(geno_pool, n_generations, n_select, n_crosses, 
                         n_offspring, rng_local, strategy='top_r'):
    """Simulate selective breeding over multiple generations.
    
    Each generation:
    1. Evaluate all individuals
    2. Select top n_select by resistance
    3. Make n_crosses random pairings from selected pool
    4. Generate n_offspring per cross
    5. New pool = all offspring
    
    Returns per-generation stats.
    """
    history = {'gen': [], 'mean_r': [], 'max_r': [], 'p99_r': [],
               'mean_t': [], 'mean_c': [], 'pool_size': [],
               'mean_loci_covered': []}
    
    current_pool = geno_pool.copy()
    
    for gen in range(n_generations + 1):
        alive = np.ones(current_pool.shape[0], dtype=bool)
        r_vals = compute_trait_batch(current_pool, effects_r, alive, res_s)
        t_vals = compute_trait_batch(current_pool, effects_t, alive, tol_s)
        c_vals = compute_trait_batch(current_pool, effects_c, alive, rec_s)
        
        # Per-individual locus coverage
        loci_covered = []
        for i in range(min(len(current_pool), 200)):  # Sample to save time
            profile = locus_profile(current_pool[i], res_s)
            loci_covered.append(np.sum(profile > 0))
        
        history['gen'].append(gen)
        history['mean_r'].append(float(np.mean(r_vals)))
        history['max_r'].append(float(np.max(r_vals)))
        history['p99_r'].append(float(np.percentile(r_vals, 99)))
        history['mean_t'].append(float(np.mean(t_vals)))
        history['mean_c'].append(float(np.mean(c_vals)))
        history['pool_size'].append(len(current_pool))
        history['mean_loci_covered'].append(float(np.mean(loci_covered)))
        
        if gen == n_generations:
            break
        
        # Selection
        if strategy == 'top_r':
            selected_idx = np.argsort(r_vals)[-n_select:]
        elif strategy == 'complementary':
            # First select top 2× by r, then pick diverse pairs
            candidates = np.argsort(r_vals)[-n_select * 2:]
            selected_idx = candidates[-n_select:]  # Still use top, but pair smartly
        
        selected_geno = current_pool[selected_idx]
        
        # Crossing
        all_offspring = []
        for _ in range(n_crosses):
            p1, p2 = rng_local.choice(len(selected_geno), size=2, replace=False)
            offspring = cross(selected_geno[p1], selected_geno[p2], n_offspring, rng_local)
            all_offspring.append(offspring)
        
        current_pool = np.concatenate(all_offspring, axis=0)
    
    return history


# Run breeding programs with different starting pools
print("  Running: post-epidemic founders...")
rng_breed = np.random.default_rng(SEED + 100)
hist_post = run_breeding_program(
    geno_post, N_GENERATIONS, TOP_K_SELECT, N_CROSSES, 
    N_OFFSPRING_PER_CROSS, rng_breed)

print("  Running: pre-epidemic founders...")
rng_breed2 = np.random.default_rng(SEED + 200)
hist_pre = run_breeding_program(
    geno_pre, N_GENERATIONS, TOP_K_SELECT, N_CROSSES,
    N_OFFSPRING_PER_CROSS, rng_breed2)

# Also run with fewer founders (realistic scenario: 200 wild-caught)
print("  Running: 200 wild-caught founders...")
rng_small = np.random.default_rng(SEED + 300)
small_pool_idx = rng_small.choice(N_FOUNDERS, size=200, replace=False)
geno_small = geno_post[small_pool_idx]
rng_breed3 = np.random.default_rng(SEED + 400)
hist_small = run_breeding_program(
    geno_small, N_GENERATIONS, min(TOP_K_SELECT, 30), 
    min(N_CROSSES, 500), N_OFFSPRING_PER_CROSS, rng_breed3)

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Top-left: Mean resistance over generations
ax = axes[0, 0]
for hist, label, color, ls in [
    (hist_post, f'Post-epidemic (N={N_FOUNDERS:,})', '#2196F3', '-'),
    (hist_pre, f'Pre-epidemic (N={N_FOUNDERS:,})', 'gray', '--'),
    (hist_small, 'Wild-caught (N=200)', '#FF9800', '-'),
]:
    ax.plot(hist['gen'], hist['mean_r'], f'{ls}o', color=color, linewidth=2, 
            markersize=6, label=label)
ax.set_xlabel("Generation", fontsize=12)
ax.set_ylabel("Mean Resistance (r)", fontsize=12)
ax.set_title("Mean Resistance Over Breeding Generations", fontsize=12, fontweight='bold')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Top-right: Max resistance over generations
ax = axes[0, 1]
for hist, label, color, ls in [
    (hist_post, 'Post-epidemic', '#2196F3', '-'),
    (hist_pre, 'Pre-epidemic', 'gray', '--'),
    (hist_small, 'Wild-caught (200)', '#FF9800', '-'),
]:
    ax.plot(hist['gen'], hist['max_r'], f'{ls}o', color=color, linewidth=2, 
            markersize=6, label=label)
    ax.plot(hist['gen'], hist['p99_r'], f'{ls}', color=color, linewidth=1, alpha=0.5)
ax.axhline(0.50, color='red', linestyle=':', alpha=0.5, label='r=0.50 (elite)')
ax.axhline(0.70, color='darkred', linestyle=':', alpha=0.5, label='r=0.70 (very high)')
ax.set_xlabel("Generation", fontsize=12)
ax.set_ylabel("Max Resistance (r)", fontsize=12)
ax.set_title("Best Individual Over Breeding Generations", fontsize=12, fontweight='bold')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Bottom-left: Loci coverage
ax = axes[1, 0]
for hist, label, color, ls in [
    (hist_post, 'Post-epidemic', '#2196F3', '-'),
    (hist_pre, 'Pre-epidemic', 'gray', '--'),
    (hist_small, 'Wild-caught (200)', '#FF9800', '-'),
]:
    ax.plot(hist['gen'], hist['mean_loci_covered'], f'{ls}o', color=color, 
            linewidth=2, markersize=6, label=label)
ax.axhline(17, color='green', linestyle=':', alpha=0.5, label='All 17 loci')
ax.set_xlabel("Generation", fontsize=12)
ax.set_ylabel("Mean # Resistance Loci Covered", fontsize=12)
ax.set_title("Locus Coverage Over Generations", fontsize=12, fontweight='bold')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_ylim(0, 18)

# Bottom-right: Breeding value (combined survival)
ax = axes[1, 1]
for hist, label, color, ls in [
    (hist_post, 'Post-epidemic', '#2196F3', '-'),
    (hist_pre, 'Pre-epidemic', 'gray', '--'),
    (hist_small, 'Wild-caught (200)', '#FF9800', '-'),
]:
    # Approximate breeding value from mean r (ignoring t, c evolution for clarity)
    bv = [r + (1 - r) * survival_probability(r, 0.10, 0.02) for r in hist['mean_r']]
    ax.plot(hist['gen'], [b * 100 for b in bv], f'{ls}o', color=color,
            linewidth=2, markersize=6, label=label)
ax.set_xlabel("Generation", fontsize=12)
ax.set_ylabel("Breeding Value (%)", fontsize=12)
ax.set_title("Expected Survival per Exposure", fontsize=12, fontweight='bold')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

fig.suptitle(f"Selective Breeding Program Simulation\n(Select top {TOP_K_SELECT}, "
             f"{N_CROSSES:,} crosses × {N_OFFSPRING_PER_CROSS} offspring/cross per gen)",
             fontsize=14, fontweight='bold')
plt.tight_layout()
fig.savefig(outdir / "fig5_breeding_program.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved fig5_breeding_program.png")

# ══════════════════════════════════════════════════════════════════════
# PART 6: Practical screening table for post-epidemic
# ══════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(12, 6))

# Compare screening requirements pre vs post
thresholds_r = np.arange(0.15, 0.55, 0.01)

for scores, label, color, ls in [
    (r_pre, 'Pre-epidemic', 'gray', '--'),
    (r_post, 'Post-epidemic (Howe Sound means)', '#2196F3', '-'),
]:
    p_exceed = np.array([np.mean(scores >= t) for t in thresholds_r])
    n_needed = np.where(p_exceed > 0, np.ceil(np.log(0.05) / np.log(1 - p_exceed)), np.inf)
    ax.semilogy(thresholds_r, n_needed, color=color, linewidth=2.5, linestyle=ls, label=label)

for n_ref, lbl in [(50, "50 (field trip)"), (200, "200 (survey)"),
                    (1000, "1,000 (major effort)")]:
    ax.axhline(n_ref, color='orange', linestyle=':', alpha=0.4)
    ax.text(0.155, n_ref * 1.15, lbl, fontsize=8, color='orange', alpha=0.7)

ax.set_xlabel("Minimum Resistance Threshold", fontsize=13)
ax.set_ylabel("Number to Screen (95% confidence)", fontsize=13)
ax.set_title("Screening Effort: Pre- vs Post-Epidemic", fontsize=14, fontweight='bold')
ax.legend(fontsize=12)
ax.set_ylim(1, 1e5)
ax.set_xlim(0.15, 0.50)
ax.grid(True, alpha=0.3)

plt.tight_layout()
fig.savefig(outdir / "fig6_screening_comparison.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved fig6_screening_comparison.png")

# ══════════════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("BREEDING ANALYSIS SUMMARY")
print("=" * 70)

print(f"\nPost-epidemic resistance (Howe Sound means, r̄={hs['r']:.4f}):")
print(f"  Mean: {np.mean(r_post):.4f}  |  Max: {np.max(r_post):.4f}  |  99th: {np.percentile(r_post,99):.4f}")

print(f"\nScreening comparison (95% confidence):")
print(f"  {'Threshold':>10} {'Pre-epidemic':>15} {'Post-epidemic':>15} {'Improvement':>15}")
print(f"  {'─'*10} {'─'*15} {'─'*15} {'─'*15}")
for thresh in [0.20, 0.25, 0.30, 0.35, 0.40]:
    p_pre = np.mean(r_pre >= thresh)
    p_post = np.mean(r_post >= thresh)
    n_pre = int(np.ceil(np.log(0.05) / np.log(1 - p_pre))) if p_pre > 0 else float('inf')
    n_post = int(np.ceil(np.log(0.05) / np.log(1 - p_post))) if p_post > 0 else float('inf')
    if n_pre < float('inf') and n_post < float('inf'):
        improvement = f"{n_pre/n_post:.1f}× easier"
    else:
        improvement = "—"
    print(f"  r≥{thresh:.2f}     {n_pre:>12,}   {n_post:>12,}   {improvement:>15}")

print(f"\nBreeding program results ({N_GENERATIONS} generations):")
for hist, label in [(hist_post, 'Post-epidemic 10K'), (hist_pre, 'Pre-epidemic 10K'), 
                     (hist_small, 'Wild-caught 200')]:
    print(f"  {label}:")
    print(f"    Gen 0: mean r={hist['mean_r'][0]:.3f}, max r={hist['max_r'][0]:.3f}")
    print(f"    Gen {N_GENERATIONS}: mean r={hist['mean_r'][-1]:.3f}, max r={hist['max_r'][-1]:.3f}")
    print(f"    Loci coverage: {hist['mean_loci_covered'][0]:.1f} → {hist['mean_loci_covered'][-1]:.1f} / 17")

# Complementarity summary
print(f"\nComplementarity (top 50 by r):")
print(f"  Max locus union: {np.max(union_matrix):.0f} / 17")
print(f"  Max unique contribution: {np.max(unique_matrix):.0f}")
print(f"  Mean union: {union_matrix[union_matrix > 0].mean():.1f}")

# Best cross info
p1, p2 = best_pair
print(f"\nBest complementary pair (from top 20):")
print(f"  Parent 1: r={r_post[p1]:.3f}")
print(f"  Parent 2: r={r_post[p2]:.3f}")
u, o, uniq = complementarity_score(geno_post[p1], geno_post[p2], res_s)
print(f"  Union: {u}/17 loci covered, {uniq} unique")
offspring_test = cross(geno_post[p1], geno_post[p2], 1000, rng)
offspring_r_test = compute_r_from_geno(offspring_test)
print(f"  Offspring: mean r={np.mean(offspring_r_test):.3f}, max r={np.max(offspring_r_test):.3f}")

print(f"\nAll figures saved to: {outdir}")
