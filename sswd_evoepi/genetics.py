"""Genetics and evolution module for SSWD-EvoEpi.

Implements the 51-locus diploid three-trait architecture:
  - 17 resistance loci (indices 0–16): immune exclusion (reduce infection probability)
  - 17 tolerance loci (indices 17–33): damage limitation (reduce mortality while infected)
  - 17 recovery loci (indices 34–50): pathogen clearance (clear infection, transition S→R)

Partition is configurable via GeneticsSection (n_resistance, n_tolerance, n_recovery).
Total always constrained to N_LOCI=51.

Core responsibilities:
  - Per-trait effect size initialization (Exp(λ), normalized, sorted descending)
  - Trait score computation (vectorized + single) for any trait block
  - Three-trait genotype initialization at Hardy-Weinberg equilibrium
  - Mutation (bidirectional, μ=10⁻⁸ per locus per generation)
  - Allele frequency tracking
  - Heterozygosity (observed + expected)
  - Additive genetic variance (V_A) per trait
  - Genetic diagnostics per node per generation
  - Genotype bank operations (Tier 2 compression/expansion)

NO cost of resistance (CE-1: Willem's decision).
NO EF1A overdominant locus (removed: Wares 2016 finding was Pisaster, not Pycnopodia).

References:
  - three-trait-genetic-architecture-spec.md §5.1–§5.10
  - CODE_ERRATA CE-1: cost_resistance removed
  - CODE_ERRATA CE-3: exponential effect sizes confirmed

Authors: Anton 🔬 & Willem Weertman
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Tuple

import numpy as np

from sswd_evoepi.types import (
    N_LOCI,
    N_RESISTANCE_DEFAULT,
    N_TOLERANCE_DEFAULT,
    N_RECOVERY_DEFAULT,
    allocate_genotypes,
    trait_slices,
)


# ═══════════════════════════════════════════════════════════════════════
# CONSTANTS
# ═══════════════════════════════════════════════════════════════════════

MU_PER_LOCUS: float = 1e-8                 # per allele per generation (Lynch 2010)

# Genotype bank size for Tier 2 nodes
N_BANK: int = 100

# Default trait slices (17/17/17 partition)
RESISTANCE_SLICE, TOLERANCE_SLICE, RECOVERY_SLICE = trait_slices(
    N_RESISTANCE_DEFAULT, N_TOLERANCE_DEFAULT, N_RECOVERY_DEFAULT
)


# ═══════════════════════════════════════════════════════════════════════
# EFFECT SIZE INITIALIZATION
# ═══════════════════════════════════════════════════════════════════════


def initialize_trait_effect_sizes(
    rng: np.random.Generator,
    n_loci: int,
    total_weight: float = 1.0,
) -> np.ndarray:
    """Draw and normalize effect sizes for one trait from Exp(λ).

    Effect sizes are sorted DESCENDING so the first locus within each
    trait block has the largest effect. Sum equals ``total_weight``
    exactly, ensuring trait score ∈ [0, total_weight] for a fully
    homozygous-derived individual.

    Args:
        rng: NumPy random Generator (for reproducibility).
        n_loci: Number of loci for this trait.
        total_weight: Maximum total trait weight (default 1.0).

    Returns:
        (n_loci,) float64 array, sorted descending. Sum = total_weight.
    """
    raw = rng.exponential(scale=1.0, size=n_loci)
    normalized = raw / raw.sum() * total_weight
    normalized.sort()
    return normalized[::-1].copy()


# Backward compatibility alias (used by model.py make_effect_sizes)
def initialize_effect_sizes(
    rng: np.random.Generator,
    n_additive: int = N_LOCI,
    w_add: float = 1.0,
) -> np.ndarray:
    """Backward-compatible wrapper around initialize_trait_effect_sizes.

    .. deprecated:: Use initialize_trait_effect_sizes directly.
    """
    return initialize_trait_effect_sizes(rng, n_additive, w_add)


# ═══════════════════════════════════════════════════════════════════════
# TRAIT SCORE COMPUTATION
# ═══════════════════════════════════════════════════════════════════════


def compute_trait_batch(
    genotypes: np.ndarray,
    effects: np.ndarray,
    alive_mask: np.ndarray,
    locus_slice: slice,
) -> np.ndarray:
    """Vectorized trait score computation for one trait.

    score_i = Σ e_l × (g_l0 + g_l1) / 2, summed over the trait's loci.

    Args:
        genotypes: (max_agents, N_LOCI, 2) int8.
        effects: (n_trait_loci,) float64 effect sizes.
        alive_mask: (max_agents,) bool.
        locus_slice: slice object for this trait's loci.

    Returns:
        (max_agents,) float32 — trait score for alive agents, 0.0 for dead.
    """
    max_agents = genotypes.shape[0]
    scores = np.zeros(max_agents, dtype=np.float32)

    alive_idx = np.where(alive_mask)[0]
    if len(alive_idx) == 0:
        return scores

    # (n_alive, n_trait_loci) allele sums → dot with effects
    allele_sums = genotypes[alive_idx, locus_slice, :].sum(axis=2)
    allele_means = allele_sums.astype(np.float64) * 0.5
    scores[alive_idx] = (allele_means @ effects).astype(np.float32)

    return scores


def compute_trait_single(
    genotype: np.ndarray,
    effects: np.ndarray,
    locus_slice: slice,
) -> float:
    """Compute trait score for a single individual.

    score = Σ e_l × (g_l0 + g_l1) / 2 over the trait's loci.

    Args:
        genotype: (N_LOCI, 2) int8 array.
        effects: (n_trait_loci,) float64 effect sizes.
        locus_slice: slice object for this trait's loci.

    Returns:
        Trait score (float).
    """
    allele_means = genotype[locus_slice, :].sum(axis=1).astype(np.float64) * 0.5
    return float(np.dot(allele_means, effects))


def compute_resistance_batch(
    genotypes: np.ndarray,
    effects: np.ndarray,
    alive_mask: np.ndarray,
    locus_slice: slice = RESISTANCE_SLICE,
) -> np.ndarray:
    """Backward-compatible resistance score computation.

    Thin wrapper around compute_trait_batch using RESISTANCE_SLICE.

    Args:
        genotypes: (max_agents, N_LOCI, 2) int8.
        effects: (n_resistance_loci,) float64 effect sizes.
        alive_mask: (max_agents,) bool.
        locus_slice: slice object (default: RESISTANCE_SLICE).

    Returns:
        (max_agents,) float32 — resistance scores.
    """
    return compute_trait_batch(genotypes, effects, alive_mask, locus_slice)


def compute_resistance_single(
    genotype: np.ndarray,
    effects: np.ndarray,
    locus_slice: slice = RESISTANCE_SLICE,
) -> float:
    """Backward-compatible single-agent resistance score.

    Thin wrapper around compute_trait_single using RESISTANCE_SLICE.
    """
    return compute_trait_single(genotype, effects, locus_slice)


# ═══════════════════════════════════════════════════════════════════════
# UPDATE ALL TRAIT SCORES
# ═══════════════════════════════════════════════════════════════════════


def update_all_trait_scores(
    agents: np.ndarray,
    genotypes: np.ndarray,
    effects_r: np.ndarray,
    effects_t: np.ndarray,
    effects_c: np.ndarray,
    res_slice: slice = RESISTANCE_SLICE,
    tol_slice: slice = TOLERANCE_SLICE,
    rec_slice: slice = RECOVERY_SLICE,
) -> None:
    """Recompute and write r_i, t_i, c_i for all alive agents in-place.

    Args:
        agents: Agent structured array (AGENT_DTYPE), modified in place.
        genotypes: (max_agents, N_LOCI, 2) int8.
        effects_r: Resistance effect sizes.
        effects_t: Tolerance effect sizes.
        effects_c: Recovery effect sizes.
        res_slice: Resistance locus slice.
        tol_slice: Tolerance locus slice.
        rec_slice: Recovery locus slice.
    """
    alive = agents['alive']
    agents['resistance'] = compute_trait_batch(genotypes, effects_r, alive, res_slice)
    agents['tolerance'] = compute_trait_batch(genotypes, effects_t, alive, tol_slice)
    agents['recovery_ability'] = compute_trait_batch(genotypes, effects_c, alive, rec_slice)


def update_resistance_scores(
    agents: np.ndarray,
    genotypes: np.ndarray,
    effects: np.ndarray,
) -> None:
    """Backward-compatible resistance-only update.

    Writes agents['resistance'] only. For use during migration period;
    new code should use update_all_trait_scores().

    Args:
        agents: Agent structured array, modified in place.
        genotypes: (max_agents, N_LOCI, 2) int8.
        effects: (n_resistance_loci,) float64 effect sizes.
    """
    alive_mask = agents['alive']
    r_vals = compute_trait_batch(genotypes, effects, alive_mask, RESISTANCE_SLICE)
    agents['resistance'] = r_vals


# ═══════════════════════════════════════════════════════════════════════
# GENOTYPE INITIALIZATION
# ═══════════════════════════════════════════════════════════════════════


def initialize_genotypes_three_trait(
    n_agents: int,
    effects_r: np.ndarray,
    effects_t: np.ndarray,
    effects_c: np.ndarray,
    rng: np.random.Generator,
    target_mean_r: float = 0.15,
    target_mean_t: float = 0.10,
    target_mean_c: float = 0.08,
    beta_a: float = 2.0,
    beta_b: float = 8.0,
    n_resistance: int = N_RESISTANCE_DEFAULT,
    n_tolerance: int = N_TOLERANCE_DEFAULT,
    n_recovery: int = N_RECOVERY_DEFAULT,
) -> np.ndarray:
    """Initialize 51-locus genotypes with per-trait target means.

    Each trait block gets Beta-distributed per-locus allele frequencies,
    scaled so the expected population-mean trait score ≈ target.

    Args:
        n_agents: Number of individuals.
        effects_r: Resistance effect sizes (n_resistance,).
        effects_t: Tolerance effect sizes (n_tolerance,).
        effects_c: Recovery effect sizes (n_recovery,).
        rng: Random generator.
        target_mean_r: Target population-mean resistance (default 0.15).
        target_mean_t: Target population-mean tolerance (default 0.10).
        target_mean_c: Target population-mean recovery (default 0.08).
        beta_a: Beta distribution shape a.
        beta_b: Beta distribution shape b.
        n_resistance: Number of resistance loci.
        n_tolerance: Number of tolerance loci.
        n_recovery: Number of recovery loci.

    Returns:
        (n_agents, N_LOCI, 2) int8 genotype array.
    """
    geno = np.zeros((n_agents, N_LOCI, 2), dtype=np.int8)

    res_s, tol_s, rec_s = trait_slices(n_resistance, n_tolerance, n_recovery)

    for locus_slice, effects, target in [
        (res_s, effects_r, target_mean_r),
        (tol_s, effects_t, target_mean_t),
        (rec_s, effects_c, target_mean_c),
    ]:
        n_trait = effects.shape[0]
        raw_q = rng.beta(beta_a, beta_b, size=n_trait)

        # Scale so E[trait_score] ≈ target
        current = np.dot(effects, raw_q)
        if current > 0:
            scale = target / current
            q_vals = np.clip(raw_q * scale, 0.001, 0.5)
        else:
            q_vals = np.full(n_trait, 0.01)

        start = locus_slice.start
        for i in range(n_trait):
            l_idx = start + i
            geno[:, l_idx, 0] = (rng.random(n_agents) < q_vals[i]).astype(np.int8)
            geno[:, l_idx, 1] = (rng.random(n_agents) < q_vals[i]).astype(np.int8)

    return geno


def initialize_genotypes(
    n_agents: int,
    effects: np.ndarray,
    rng: np.random.Generator,
    target_mean_r: float = 0.15,
    q_additive: Optional[np.ndarray] = None,
    **kwargs,
) -> np.ndarray:
    """Backward-compatible genotype initialization (resistance only).

    Initializes resistance loci only. Tolerance and recovery loci are
    left at zero (no derived alleles). For full three-trait initialization,
    use initialize_genotypes_three_trait().

    Args:
        n_agents: Number of individuals.
        effects: (n_resistance,) float64 effect sizes.
        rng: Random generator.
        target_mean_r: Target population-mean resistance.
        q_additive: Optional per-locus allele frequencies.
        **kwargs: Absorbs legacy arguments (ef1a_q, etc.).

    Returns:
        (n_agents, N_LOCI, 2) int8 genotype array.
    """
    geno = np.zeros((n_agents, N_LOCI, 2), dtype=np.int8)
    n_res = len(effects)

    if q_additive is not None:
        q_vals = np.clip(q_additive, 0.0, 1.0)
    else:
        # Uniform q across resistance loci to hit target mean
        q_uniform = target_mean_r / effects.sum() if effects.sum() > 0 else 0.01
        q_uniform = np.clip(q_uniform, 0.001, 0.5)
        q_vals = np.full(n_res, q_uniform)

    for l_idx in range(n_res):
        geno[:, l_idx, 0] = (rng.random(n_agents) < q_vals[l_idx]).astype(np.int8)
        geno[:, l_idx, 1] = (rng.random(n_agents) < q_vals[l_idx]).astype(np.int8)

    return geno


def initialize_genotypes_beta(
    n_agents: int,
    effects: np.ndarray,
    rng: np.random.Generator,
    target_mean_r: float = 0.15,
    beta_a: float = 2.0,
    beta_b: float = 8.0,
    **kwargs,
) -> np.ndarray:
    """Backward-compatible Beta-distributed genotype initialization.

    Initializes resistance loci only with Beta-distributed per-locus
    allele frequencies. For full three-trait initialization, use
    initialize_genotypes_three_trait().

    Args:
        n_agents: Number of individuals.
        effects: (n_resistance,) float64 effect sizes.
        rng: Random generator.
        target_mean_r: Target population-mean resistance.
        beta_a: Beta distribution shape a.
        beta_b: Beta distribution shape b.
        **kwargs: Absorbs legacy arguments (ef1a_q, etc.).

    Returns:
        (n_agents, N_LOCI, 2) int8 genotype array.
    """
    n_res = len(effects)
    raw_q = rng.beta(beta_a, beta_b, size=n_res)

    # Scale so E[resistance] ≈ target_mean_r
    current = np.dot(effects, raw_q)
    if current > 0:
        scale = target_mean_r / current
        q_vals = np.clip(raw_q * scale, 0.001, 0.5)
    else:
        q_vals = np.full(n_res, 0.01)

    return initialize_genotypes(
        n_agents, effects, rng,
        target_mean_r=target_mean_r,
        q_additive=q_vals,
    )


# ═══════════════════════════════════════════════════════════════════════
# MUTATION
# ═══════════════════════════════════════════════════════════════════════


def apply_mutations(
    offspring_geno: np.ndarray,
    rng: np.random.Generator,
    mu: float = MU_PER_LOCUS,
) -> int:
    """Apply bidirectional point mutations to offspring genotypes.

    Each allele copy independently mutates with probability μ.
    Mutation flips the allele: 0→1 or 1→0. Operates on the full
    genotype array — no distinction between trait blocks.

    Args:
        offspring_geno: (n_offspring, N_LOCI, 2) int8 — modified in place.
        rng: Random generator.
        mu: Per-allele per-generation mutation rate.

    Returns:
        Number of mutations applied (typically 0 at μ=10⁻⁸).
    """
    n_offspring, n_loci, ploidy = offspring_geno.shape
    n_alleles = n_offspring * n_loci * ploidy
    n_mutations = rng.poisson(mu * n_alleles)

    if n_mutations == 0:
        return 0

    # Limit to n_alleles (can't mutate more sites than exist)
    n_mutations = min(n_mutations, n_alleles)

    # Random flat indices into the allele array
    mut_flat_idx = rng.choice(n_alleles, size=n_mutations, replace=False)

    # Vectorized mutation application using unravel_index
    if n_mutations > 0:
        i_indices, l_indices, a_indices = np.unravel_index(
            mut_flat_idx, (n_offspring, n_loci, ploidy)
        )
        offspring_geno[i_indices, l_indices, a_indices] = (
            1 - offspring_geno[i_indices, l_indices, a_indices]
        )

    return n_mutations


# ═══════════════════════════════════════════════════════════════════════
# ALLELE FREQUENCY & HETEROZYGOSITY
# ═══════════════════════════════════════════════════════════════════════


def compute_allele_frequencies(
    genotypes: np.ndarray,
    alive_mask: np.ndarray,
) -> np.ndarray:
    """Compute allele frequency at each locus for alive individuals.

    q_l = mean allele value across all alive individuals and both copies.

    Args:
        genotypes: (max_agents, N_LOCI, 2) int8.
        alive_mask: (max_agents,) bool.

    Returns:
        (N_LOCI,) float64 — allele frequency q at each locus.
        Returns zeros if no alive individuals.
    """
    alive_idx = np.where(alive_mask)[0]
    n = len(alive_idx)
    if n == 0:
        return np.zeros(N_LOCI, dtype=np.float64)

    alive_geno = genotypes[alive_idx]  # (n, N_LOCI, 2)
    return alive_geno.sum(axis=(0, 2)).astype(np.float64) / (2.0 * n)


def compute_heterozygosity(
    genotypes: np.ndarray,
    alive_mask: np.ndarray,
) -> Tuple[float, float]:
    """Compute observed and expected heterozygosity across all loci.

    H_o = mean fraction of heterozygous individuals per locus.
    H_e = mean 2pq per locus (expected under HWE).

    Args:
        genotypes: (max_agents, N_LOCI, 2) int8.
        alive_mask: (max_agents,) bool.

    Returns:
        (H_o, H_e) averaged across all N_LOCI loci.
    """
    alive_idx = np.where(alive_mask)[0]
    n = len(alive_idx)
    if n < 2:
        return 0.0, 0.0

    alive_geno = genotypes[alive_idx]  # (n, N_LOCI, 2)

    # Observed heterozygosity: fraction heterozygous at each locus
    het_per_locus = np.mean(alive_geno[:, :, 0] != alive_geno[:, :, 1], axis=0)
    H_o = float(np.mean(het_per_locus))

    # Expected heterozygosity: 2pq at each locus
    q = alive_geno.sum(axis=2).astype(np.float64).mean(axis=0) / 2.0
    H_e_per_locus = 2.0 * q * (1.0 - q)
    H_e = float(np.mean(H_e_per_locus))

    return H_o, H_e


def compute_additive_variance(
    allele_freq: np.ndarray,
    effects: np.ndarray,
    locus_slice: Optional[slice] = None,
) -> float:
    """Compute additive genetic variance V_A for a trait.

    V_A = Σ 2 × e_l² × q_l × (1 − q_l)

    Args:
        allele_freq: (N_LOCI,) float64 allele frequencies (full array).
        effects: (n_trait_loci,) float64 effect sizes.
        locus_slice: Slice into allele_freq for this trait's loci.
            If None, uses first len(effects) loci (backward compat).

    Returns:
        V_A (float).
    """
    if locus_slice is not None:
        q = allele_freq[locus_slice]
    else:
        q = allele_freq[:len(effects)]
    return float(2.0 * np.sum(effects**2 * q * (1.0 - q)))


# ═══════════════════════════════════════════════════════════════════════
# GENETIC DIAGNOSTICS
# ═══════════════════════════════════════════════════════════════════════


@dataclass
class GeneticDiagnostics:
    """Genetic summary statistics for a single node in a single generation.

    Three-trait architecture: tracks mean/var/V_A for resistance, tolerance,
    and recovery independently.
    """

    # Allele frequencies (all 51 loci)
    allele_freq: np.ndarray = field(
        default_factory=lambda: np.zeros(N_LOCI, dtype=np.float64)
    )

    # Per-trait distribution stats
    mean_resistance: float = 0.0
    var_resistance: float = 0.0
    va_resistance: float = 0.0

    mean_tolerance: float = 0.0
    var_tolerance: float = 0.0
    va_tolerance: float = 0.0

    mean_recovery: float = 0.0
    var_recovery: float = 0.0
    va_recovery: float = 0.0

    # Diversity (across all 51 loci)
    heterozygosity_obs: float = 0.0
    heterozygosity_exp: float = 0.0

    # Allele frequency change (relative to previous generation)
    delta_q_top3: np.ndarray = field(
        default_factory=lambda: np.zeros(3, dtype=np.float64)
    )
    mean_abs_delta_q: float = 0.0

    # Population count
    n_alive: int = 0


def compute_genetic_diagnostics(
    agents: np.ndarray,
    genotypes: np.ndarray,
    effects_r: np.ndarray,
    effects_t: Optional[np.ndarray] = None,
    effects_c: Optional[np.ndarray] = None,
    prev_allele_freq: Optional[np.ndarray] = None,
    res_slice: slice = RESISTANCE_SLICE,
    tol_slice: slice = TOLERANCE_SLICE,
    rec_slice: slice = RECOVERY_SLICE,
) -> GeneticDiagnostics:
    """Compute all genetic summary statistics for a node.

    Args:
        agents: Agent structured array.
        genotypes: (max_agents, N_LOCI, 2) int8.
        effects_r: Resistance effect sizes.
        effects_t: Tolerance effect sizes (optional for backward compat).
        effects_c: Recovery effect sizes (optional for backward compat).
        prev_allele_freq: (N_LOCI,) float64 from previous generation,
            or None if first generation.
        res_slice: Resistance locus slice.
        tol_slice: Tolerance locus slice.
        rec_slice: Recovery locus slice.

    Returns:
        GeneticDiagnostics with all fields populated.
    """
    alive_mask = agents['alive']
    n_alive = int(alive_mask.sum())

    diag = GeneticDiagnostics(n_alive=n_alive)

    if n_alive < 2:
        return diag

    # Allele frequencies
    diag.allele_freq = compute_allele_frequencies(genotypes, alive_mask)

    # Resistance distribution
    r_vals = agents['resistance'][alive_mask]
    diag.mean_resistance = float(np.mean(r_vals))
    diag.var_resistance = float(np.var(r_vals))
    diag.va_resistance = compute_additive_variance(
        diag.allele_freq, effects_r, res_slice
    )

    # Tolerance distribution
    if effects_t is not None:
        t_vals = agents['tolerance'][alive_mask]
        diag.mean_tolerance = float(np.mean(t_vals))
        diag.var_tolerance = float(np.var(t_vals))
        diag.va_tolerance = compute_additive_variance(
            diag.allele_freq, effects_t, tol_slice
        )

    # Recovery distribution
    if effects_c is not None:
        c_vals = agents['recovery_ability'][alive_mask]
        diag.mean_recovery = float(np.mean(c_vals))
        diag.var_recovery = float(np.var(c_vals))
        diag.va_recovery = compute_additive_variance(
            diag.allele_freq, effects_c, rec_slice
        )

    # Heterozygosity (all 51 loci)
    diag.heterozygosity_obs, diag.heterozygosity_exp = compute_heterozygosity(
        genotypes, alive_mask
    )

    # Allele frequency change
    if prev_allele_freq is not None:
        delta_q = diag.allele_freq - prev_allele_freq
        abs_delta = np.abs(delta_q)
        top3_idx = np.argsort(abs_delta)[-3:][::-1]
        diag.delta_q_top3 = delta_q[top3_idx].copy()
        diag.mean_abs_delta_q = float(np.mean(abs_delta))

    return diag


# ═══════════════════════════════════════════════════════════════════════
# Ne/N COMPUTATION (wrapper around reproduction module)
# ═══════════════════════════════════════════════════════════════════════


def compute_ne_ratio_from_offspring(offspring_per_parent: np.ndarray) -> float:
    """Compute Ne/N from realized offspring distribution.

    Ne = (4N − 2) / (Vk + 2), where Vk = variance in offspring count.

    Args:
        offspring_per_parent: (N_parents,) int array.

    Returns:
        Ne/N ratio. Returns 1.0 if < 2 parents.
    """
    N = len(offspring_per_parent)
    if N < 2:
        return 1.0
    Vk = float(np.var(offspring_per_parent, ddof=1))
    Ne = (4.0 * N - 2.0) / (Vk + 2.0)
    return Ne / N


# ═══════════════════════════════════════════════════════════════════════
# GENOTYPE BANK — TIER 2 OPERATIONS
# ═══════════════════════════════════════════════════════════════════════


@dataclass
class GenotypeBankState:
    """Genetic state at a Tier 2 node.

    Instead of individual genotypes, stores a bank of N_BANK representative
    diploid genotypes with associated weights and summary statistics.
    Tracks all three trait scores.
    """
    bank: np.ndarray                   # (N_BANK, N_LOCI, 2) int8
    bank_resistance: np.ndarray        # (N_BANK,) float32
    bank_tolerance: np.ndarray         # (N_BANK,) float32
    bank_recovery: np.ndarray          # (N_BANK,) float32
    bank_weights: np.ndarray           # (N_BANK,) float64 — frequency weights

    # Population summary statistics
    mean_resistance: float = 0.0
    var_resistance: float = 0.0
    mean_tolerance: float = 0.0
    var_tolerance: float = 0.0
    mean_recovery: float = 0.0
    var_recovery: float = 0.0
    allele_freq: np.ndarray = field(
        default_factory=lambda: np.zeros(N_LOCI, dtype=np.float64)
    )
    heterozygosity: float = 0.0

    def update_summary(
        self,
        effects_r: np.ndarray,
        effects_t: Optional[np.ndarray] = None,
        effects_c: Optional[np.ndarray] = None,
    ) -> None:
        """Recompute summary statistics from the genotype bank.

        Args:
            effects_r: Resistance effect sizes.
            effects_t: Tolerance effect sizes (optional).
            effects_c: Recovery effect sizes (optional).
        """
        n_bank = len(self.bank_weights)

        # Weighted mean and variance of resistance
        weighted_r = float(np.average(self.bank_resistance, weights=self.bank_weights))
        self.mean_resistance = weighted_r
        self.var_resistance = float(
            np.average(
                (self.bank_resistance - weighted_r) ** 2,
                weights=self.bank_weights,
            )
        )

        # Weighted mean and variance of tolerance
        weighted_t = float(np.average(self.bank_tolerance, weights=self.bank_weights))
        self.mean_tolerance = weighted_t
        self.var_tolerance = float(
            np.average(
                (self.bank_tolerance - weighted_t) ** 2,
                weights=self.bank_weights,
            )
        )

        # Weighted mean and variance of recovery
        weighted_c = float(np.average(self.bank_recovery, weights=self.bank_weights))
        self.mean_recovery = weighted_c
        self.var_recovery = float(
            np.average(
                (self.bank_recovery - weighted_c) ** 2,
                weights=self.bank_weights,
            )
        )

        # Allele frequencies (weighted)
        for l_idx in range(N_LOCI):
            q_l = np.average(
                self.bank[:, l_idx, :].sum(axis=1).astype(np.float64) / 2.0,
                weights=self.bank_weights,
            )
            self.allele_freq[l_idx] = q_l

        # Observed heterozygosity (weighted)
        het_per_locus = np.zeros(N_LOCI, dtype=np.float64)
        for l_idx in range(N_LOCI):
            is_het = (self.bank[:, l_idx, 0] != self.bank[:, l_idx, 1]).astype(np.float64)
            het_per_locus[l_idx] = np.average(is_het, weights=self.bank_weights)
        self.heterozygosity = float(np.mean(het_per_locus))


def compress_to_genotype_bank(
    genotypes: np.ndarray,
    alive_mask: np.ndarray,
    effects_r: np.ndarray,
    rng: np.random.Generator,
    n_bank: int = N_BANK,
    effects_t: Optional[np.ndarray] = None,
    effects_c: Optional[np.ndarray] = None,
    res_slice: slice = RESISTANCE_SLICE,
    tol_slice: slice = TOLERANCE_SLICE,
    rec_slice: slice = RECOVERY_SLICE,
) -> GenotypeBankState:
    """Compress individual genotypes to a Tier 2 genotype bank.

    Randomly samples ``n_bank`` individuals from the alive population.
    If fewer than ``n_bank`` are alive, pads with duplicates.

    Args:
        genotypes: (max_agents, N_LOCI, 2) int8.
        alive_mask: (max_agents,) bool.
        effects_r: Resistance effect sizes.
        rng: Random generator.
        n_bank: Bank size (default 100).
        effects_t: Tolerance effect sizes (optional).
        effects_c: Recovery effect sizes (optional).
        res_slice: Resistance locus slice.
        tol_slice: Tolerance locus slice.
        rec_slice: Recovery locus slice.

    Returns:
        GenotypeBankState with summary statistics computed.
    """
    alive_idx = np.where(alive_mask)[0]
    n_alive = len(alive_idx)

    if n_alive == 0:
        bank = np.zeros((n_bank, N_LOCI, 2), dtype=np.int8)
        return GenotypeBankState(
            bank=bank,
            bank_resistance=np.zeros(n_bank, dtype=np.float32),
            bank_tolerance=np.zeros(n_bank, dtype=np.float32),
            bank_recovery=np.zeros(n_bank, dtype=np.float32),
            bank_weights=np.ones(n_bank, dtype=np.float64) / n_bank,
        )

    if n_alive <= n_bank:
        selected = alive_idx
        if n_alive < n_bank:
            extra = rng.choice(alive_idx, size=n_bank - n_alive, replace=True)
            selected = np.concatenate([selected, extra])
    else:
        selected = rng.choice(alive_idx, size=n_bank, replace=False)

    bank_geno = genotypes[selected].copy()

    # Compute trait scores for bank entries
    bank_r = np.array(
        [compute_trait_single(bank_geno[i], effects_r, res_slice) for i in range(n_bank)],
        dtype=np.float32,
    )
    bank_t = np.zeros(n_bank, dtype=np.float32)
    bank_c = np.zeros(n_bank, dtype=np.float32)
    if effects_t is not None:
        bank_t = np.array(
            [compute_trait_single(bank_geno[i], effects_t, tol_slice) for i in range(n_bank)],
            dtype=np.float32,
        )
    if effects_c is not None:
        bank_c = np.array(
            [compute_trait_single(bank_geno[i], effects_c, rec_slice) for i in range(n_bank)],
            dtype=np.float32,
        )

    weights = np.ones(n_bank, dtype=np.float64) / n_bank

    state = GenotypeBankState(
        bank=bank_geno,
        bank_resistance=bank_r,
        bank_tolerance=bank_t,
        bank_recovery=bank_c,
        bank_weights=weights,
    )
    state.update_summary(effects_r, effects_t, effects_c)
    return state


def expand_genotype_bank(
    bank: GenotypeBankState,
    n_agents: int,
    effects_r: np.ndarray,
    rng: np.random.Generator,
    alpha_srs: float = 1.35,
    effects_t: Optional[np.ndarray] = None,
    effects_c: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Expand Tier 2 genotype bank to individual agent genotypes (Tier 1).

    Uses SRS-weighted sampling (NOT uniform) to preserve the heavy-tailed
    genetic structure (CE-8 / spec §10.3).

    Args:
        bank: Tier 2 genotype bank.
        n_agents: Number of individual genotypes to produce.
        effects_r: Resistance effect sizes.
        rng: Random generator.
        alpha_srs: Pareto shape for SRS weighting.
        effects_t: Tolerance effect sizes (accepted for API consistency).
        effects_c: Recovery effect sizes (accepted for API consistency).

    Returns:
        (n_agents, N_LOCI, 2) int8 genotype array.
    """
    n_bank = len(bank.bank_weights)

    # SRS-weighted sampling from bank
    raw_weights = rng.pareto(alpha_srs, size=n_bank) + 1.0
    combined = raw_weights * bank.bank_weights
    prob = combined / combined.sum()

    parent_idx = rng.choice(n_bank, size=n_agents, p=prob)
    agent_genotypes = bank.bank[parent_idx].copy()

    # Apply one round of mutation for minor variation
    apply_mutations(agent_genotypes, rng)

    return agent_genotypes


# ═══════════════════════════════════════════════════════════════════════
# F_ST COMPUTATION (across multiple nodes)
# ═══════════════════════════════════════════════════════════════════════


def compute_fst(allele_freqs_per_node: list[np.ndarray]) -> float:
    """Compute Weir-Cockerham-style F_ST across nodes.

    Uses the standard formula: F_ST = Var(q) / (q̄ × (1 − q̄)),
    averaged across all N_LOCI loci.

    Args:
        allele_freqs_per_node: List of (N_LOCI,) float64 arrays,
            one per node.

    Returns:
        F_ST (float). Returns 0.0 if < 2 nodes or no variation.
    """
    if len(allele_freqs_per_node) < 2:
        return 0.0

    freqs = np.array(allele_freqs_per_node)  # (n_nodes, N_LOCI)
    q_bar = freqs.mean(axis=0)               # (N_LOCI,)
    var_q = freqs.var(axis=0)                 # (N_LOCI,)

    denom = q_bar * (1.0 - q_bar)
    valid = denom > 1e-10  # only polymorphic loci
    if not np.any(valid):
        return 0.0

    fst_per_locus = var_q[valid] / denom[valid]
    return float(np.mean(fst_per_locus))


# ═══════════════════════════════════════════════════════════════════════
# HARDY-WEINBERG EQUILIBRIUM TEST
# ═══════════════════════════════════════════════════════════════════════


def hardy_weinberg_test(
    genotypes: np.ndarray,
    alive_mask: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute observed vs expected genotype frequencies under HWE.

    For each of the N_LOCI=51 loci, returns counts of (0/0, 0/1, 1/1)
    observed and expected genotype frequencies.

    Args:
        genotypes: (max_agents, N_LOCI, 2) int8.
        alive_mask: (max_agents,) bool.

    Returns:
        Tuple of:
          obs_freq: (N_LOCI, 3) float64 — observed frequencies of (0/0, het, 1/1)
          exp_freq: (N_LOCI, 3) float64 — HWE expected frequencies
          chi2: (N_LOCI,) float64 — chi-squared statistic per locus
    """
    alive_idx = np.where(alive_mask)[0]
    n = len(alive_idx)

    obs_freq = np.zeros((N_LOCI, 3), dtype=np.float64)
    exp_freq = np.zeros((N_LOCI, 3), dtype=np.float64)
    chi2 = np.zeros(N_LOCI, dtype=np.float64)

    if n < 2:
        return obs_freq, exp_freq, chi2

    alive_geno = genotypes[alive_idx]

    for l_idx in range(N_LOCI):
        allele_sums = alive_geno[:, l_idx, :].sum(axis=1)  # 0, 1, or 2
        n_00 = np.sum(allele_sums == 0)
        n_01 = np.sum(allele_sums == 1)
        n_11 = np.sum(allele_sums == 2)

        obs_freq[l_idx] = [n_00 / n, n_01 / n, n_11 / n]

        q = (n_01 + 2 * n_11) / (2.0 * n)
        p = 1.0 - q
        exp_freq[l_idx] = [p * p, 2 * p * q, q * q]

        # Chi-squared
        for g_idx in range(3):
            expected = exp_freq[l_idx, g_idx] * n
            if expected > 0:
                chi2[l_idx] += (obs_freq[l_idx, g_idx] * n - expected) ** 2 / expected

    return obs_freq, exp_freq, chi2
