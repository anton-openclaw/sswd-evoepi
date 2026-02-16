"""Genetics and evolution module for SSWD-EvoEpi.

Implements the 52-locus diploid resistance architecture:
  - 51 additive loci (exponentially distributed effect sizes)
  - 1 overdominant locus (EF1A analog: heterozygote advantage, lethal homozygote)

Core responsibilities:
  - Effect size initialization (Exp(λ), normalized, sorted descending)
  - Resistance score r_i computation (vectorized + single)
  - Genotype initialization at Hardy-Weinberg equilibrium
  - Mutation (bidirectional, μ=10⁻⁸ per locus per generation)
  - EF1A lethal homozygote elimination
  - Allele frequency tracking
  - Heterozygosity (observed + expected)
  - Additive genetic variance (V_A)
  - Genetic diagnostics per node per generation
  - Genotype bank operations (Tier 2 compression/expansion)

NO cost of resistance (CE-1: Willem's decision).

References:
  - genetics-evolution-spec.md §1–§14
  - CODE_ERRATA CE-1: cost_resistance removed; fecundity_mod always 1.0
  - CODE_ERRATA CE-3: exponential effect sizes confirmed

Authors: Anton 🔬 & Willem Weertman
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Tuple

import numpy as np

from sswd_evoepi.types import (
    IDX_EF1A,
    N_ADDITIVE,
    N_LOCI,
    allocate_genotypes,
)


# ═══════════════════════════════════════════════════════════════════════
# CONSTANTS
# ═══════════════════════════════════════════════════════════════════════

S_HET: float = 0.19                        # EF1A heterozygote advantage (Wares 2016)
W_OD: float = S_HET / (1.0 + S_HET)       # ≈ 0.160 — overdominant weight
W_ADD: float = 1.0 - W_OD                  # ≈ 0.840 — max additive weight

MU_PER_LOCUS: float = 1e-8                 # per allele per generation (Lynch 2010)

# Genotype bank size for Tier 2 nodes
N_BANK: int = 100


# ═══════════════════════════════════════════════════════════════════════
# EFFECT SIZE INITIALIZATION
# ═══════════════════════════════════════════════════════════════════════


def initialize_effect_sizes(
    rng: np.random.Generator,
    n_additive: int = N_ADDITIVE,
    w_add: float = W_ADD,
) -> np.ndarray:
    """Draw and normalize additive effect sizes from Exp(λ).

    Effect sizes are sorted DESCENDING so locus 0 = largest effect.
    Sum equals ``w_add`` exactly, ensuring r_i ∈ [0, 1] when combined
    with the overdominant component.

    Args:
        rng: NumPy random Generator (for reproducibility).
        n_additive: Number of additive loci (default 51).
        w_add: Maximum total additive weight (default ≈ 0.840).

    Returns:
        (n_additive,) float64 array, sorted descending. Sum = w_add.
    """
    raw = rng.exponential(scale=1.0, size=n_additive)
    normalized = raw / raw.sum() * w_add
    normalized.sort()
    return normalized[::-1].copy()


# ═══════════════════════════════════════════════════════════════════════
# RESISTANCE SCORE COMPUTATION
# ═══════════════════════════════════════════════════════════════════════


def compute_resistance_single(
    genotype: np.ndarray,
    effects: np.ndarray,
    w_od: float = W_OD,
) -> float:
    """Compute individual resistance score r_i from a single genotype.

    r_i = Σ ẽ_l × (g_l0 + g_l1)/2  +  δ_het × w_od

    Args:
        genotype: (N_LOCI, 2) int8 array.
        effects: (N_ADDITIVE,) float64 effect sizes.
        w_od: Overdominant weight for EF1A heterozygotes.

    Returns:
        r_i ∈ [0, 1].
    """
    # Additive component
    allele_means = genotype[:N_ADDITIVE, :].sum(axis=1).astype(np.float64) * 0.5
    additive = float(np.dot(allele_means, effects))

    # Overdominant: bonus only for heterozygotes (sum == 1)
    ef1a_sum = int(genotype[IDX_EF1A, 0]) + int(genotype[IDX_EF1A, 1])
    od_bonus = w_od if ef1a_sum == 1 else 0.0

    return additive + od_bonus


def compute_resistance_batch(
    genotypes: np.ndarray,
    effects: np.ndarray,
    alive_mask: np.ndarray,
    w_od: float = W_OD,
) -> np.ndarray:
    """Vectorized r_i computation for all alive agents at a node.

    Args:
        genotypes: (max_agents, N_LOCI, 2) int8.
        effects: (N_ADDITIVE,) float64 effect sizes.
        alive_mask: (max_agents,) bool — which agents are alive.
        w_od: Overdominant weight.

    Returns:
        (max_agents,) float32 — r_i for alive agents, 0.0 for dead.
    """
    max_agents = genotypes.shape[0]
    resistance = np.zeros(max_agents, dtype=np.float32)

    alive_idx = np.where(alive_mask)[0]
    if len(alive_idx) == 0:
        return resistance

    # Additive component: (n_alive, N_ADDITIVE) allele means → dot with effects
    allele_sums = genotypes[alive_idx, :N_ADDITIVE, :].sum(axis=2)  # (n_alive, 51)
    allele_means = allele_sums.astype(np.float64) * 0.5
    additive = allele_means @ effects  # (n_alive,)

    # Overdominant component
    ef1a_sum = genotypes[alive_idx, IDX_EF1A, :].sum(axis=1)  # (n_alive,)
    od_bonus = np.where(ef1a_sum == 1, w_od, 0.0)

    resistance[alive_idx] = (additive + od_bonus).astype(np.float32)

    return resistance


def update_resistance_scores(
    agents: np.ndarray,
    genotypes: np.ndarray,
    effects: np.ndarray,
    w_od: float = W_OD,
) -> None:
    """Recompute and write r_i for all alive agents in-place.

    Writes to agents['resistance']. CE-1: fecundity_mod stays 1.0.

    Args:
        agents: Agent structured array (AGENT_DTYPE), modified in place.
        genotypes: (max_agents, N_LOCI, 2) int8.
        effects: (N_ADDITIVE,) float64 effect sizes.
        w_od: Overdominant weight.
    """
    alive_mask = agents['alive']
    r_vals = compute_resistance_batch(genotypes, effects, alive_mask, w_od)
    agents['resistance'] = r_vals
    # CE-1: fecundity_mod always 1.0 — no cost of resistance
    agents['fecundity_mod'][alive_mask] = 1.0


# ═══════════════════════════════════════════════════════════════════════
# GENOTYPE INITIALIZATION
# ═══════════════════════════════════════════════════════════════════════


def initialize_genotypes(
    n_agents: int,
    effects: np.ndarray,
    rng: np.random.Generator,
    target_mean_r: float = 0.08,
    ef1a_q: float = 0.24,
    q_additive: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Initialize genotypes at approximate pre-SSWD Hardy-Weinberg equilibrium.

    Strategy:
      1. If ``q_additive`` is provided, use those per-locus frequencies.
         Otherwise compute a uniform q across additive loci to hit target mean r.
      2. Draw alleles at each locus from Bernoulli(q_l), independently per copy.
      3. EF1A: draw from Bernoulli(ef1a_q), then replace lethal ins/ins with
         heterozygotes (ensures no lethals in initial population).

    Args:
        n_agents: Number of individuals.
        effects: (N_ADDITIVE,) float64 effect sizes.
        rng: Random generator.
        target_mean_r: Target population-mean resistance (default 0.08).
        ef1a_q: Initial EF1A insertion allele frequency (default 0.24).
        q_additive: Optional (N_ADDITIVE,) per-locus allele frequencies.
            If None, a uniform q is computed from target_mean_r.

    Returns:
        (n_agents, N_LOCI, 2) int8 genotype array.
    """
    geno = np.zeros((n_agents, N_LOCI, 2), dtype=np.int8)

    # --- EF1A locus ---
    ef1a_alleles = rng.random((n_agents, 2)) < ef1a_q
    geno[:, IDX_EF1A, :] = ef1a_alleles.astype(np.int8)

    # Eliminate lethals: if both copies are 1 (ins/ins), flip second to 0
    lethal_mask = geno[:, IDX_EF1A, :].sum(axis=1) == 2
    geno[lethal_mask, IDX_EF1A, 1] = 0

    # --- Additive loci ---
    if q_additive is not None:
        q_vals = np.clip(q_additive, 0.0, 1.0)
    else:
        # EF1A mean bonus
        # After lethal elimination, heterozygote frequency ≈ 2pq (for small q²)
        ef1a_het_contrib = 2.0 * ef1a_q * (1.0 - ef1a_q) * W_OD
        target_additive = max(0.0, target_mean_r - ef1a_het_contrib)
        # Sum of effects × q_uniform = target_additive → q_uniform = target / sum
        q_uniform = target_additive / effects.sum() if effects.sum() > 0 else 0.01
        q_uniform = np.clip(q_uniform, 0.001, 0.5)
        q_vals = np.full(N_ADDITIVE, q_uniform)

    # Draw alleles independently per locus per copy
    for l_idx in range(N_ADDITIVE):
        geno[:, l_idx, 0] = (rng.random(n_agents) < q_vals[l_idx]).astype(np.int8)
        geno[:, l_idx, 1] = (rng.random(n_agents) < q_vals[l_idx]).astype(np.int8)

    return geno


def initialize_genotypes_beta(
    n_agents: int,
    effects: np.ndarray,
    rng: np.random.Generator,
    target_mean_r: float = 0.08,
    ef1a_q: float = 0.24,
    beta_a: float = 1.0,
    beta_b: float = 20.0,
) -> np.ndarray:
    """Initialize genotypes with per-locus frequencies from a Beta distribution.

    Draws per-locus allele frequencies from Beta(a, b), then scales them
    so the expected population mean r matches ``target_mean_r``. This
    produces realistic variation in allele frequencies across loci.

    Args:
        n_agents: Number of individuals.
        effects: (N_ADDITIVE,) float64 effect sizes.
        rng: Random generator.
        target_mean_r: Target population-mean resistance.
        ef1a_q: Initial EF1A allele frequency.
        beta_a: Beta distribution shape a (default 1.0).
        beta_b: Beta distribution shape b (default 20.0 → low mean, right-skewed).

    Returns:
        (n_agents, N_LOCI, 2) int8 genotype array.
    """
    raw_q = rng.beta(beta_a, beta_b, size=N_ADDITIVE)
    # Scale so E[r_additive] = target_additive
    ef1a_het_contrib = 2.0 * ef1a_q * (1.0 - ef1a_q) * W_OD
    target_additive = max(0.0, target_mean_r - ef1a_het_contrib)

    # Current expected additive = Σ e_l × q_l
    current_additive = np.dot(effects, raw_q)
    if current_additive > 0:
        scale = target_additive / current_additive
        q_vals = np.clip(raw_q * scale, 0.001, 0.5)
    else:
        q_vals = np.full(N_ADDITIVE, 0.01)

    return initialize_genotypes(
        n_agents, effects, rng,
        target_mean_r=target_mean_r,
        ef1a_q=ef1a_q,
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
    Mutation flips the allele: 0→1 or 1→0.

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
        # Convert flat indices to 3D coordinates all at once
        i_indices, l_indices, a_indices = np.unravel_index(
            mut_flat_idx, (n_offspring, n_loci, ploidy)
        )
        # Vectorized bit flip using advanced indexing
        offspring_geno[i_indices, l_indices, a_indices] = (
            1 - offspring_geno[i_indices, l_indices, a_indices]
        )

    return n_mutations


# ═══════════════════════════════════════════════════════════════════════
# EF1A LETHAL ELIMINATION
# ═══════════════════════════════════════════════════════════════════════


def eliminate_ef1a_lethals(offspring_geno: np.ndarray) -> np.ndarray:
    """Return boolean mask of viable offspring (not homozygous ins/ins at EF1A).

    Args:
        offspring_geno: (n_offspring, N_LOCI, 2) int8.

    Returns:
        (n_offspring,) bool mask — True = viable (not lethal).
    """
    ef1a_sum = offspring_geno[:, IDX_EF1A, :].sum(axis=1)
    return ef1a_sum < 2


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
    # q_l = sum over individuals and copies / (2n)
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
        (H_o, H_e) averaged across all loci.
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
) -> float:
    """Compute additive genetic variance V_A at additive loci.

    V_A = Σ 2 × e_l² × q_l × (1 − q_l)

    Args:
        allele_freq: (N_LOCI,) or (N_ADDITIVE,) float64 allele frequencies.
            If N_LOCI, only first N_ADDITIVE entries are used.
        effects: (N_ADDITIVE,) float64 effect sizes.

    Returns:
        V_A (float).
    """
    q = allele_freq[:N_ADDITIVE]
    return float(2.0 * np.sum(effects**2 * q * (1.0 - q)))


# ═══════════════════════════════════════════════════════════════════════
# GENETIC DIAGNOSTICS
# ═══════════════════════════════════════════════════════════════════════


@dataclass
class GeneticDiagnostics:
    """Genetic summary statistics for a single node in a single generation."""

    # Allele frequencies
    allele_freq: np.ndarray = field(
        default_factory=lambda: np.zeros(N_LOCI, dtype=np.float64)
    )

    # Resistance distribution
    mean_resistance: float = 0.0
    var_resistance: float = 0.0
    va_additive: float = 0.0     # Additive genetic variance V_A

    # EF1A statistics
    ef1a_allele_freq: float = 0.0
    ef1a_het_freq: float = 0.0   # Observed heterozygote frequency at EF1A

    # Diversity
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
    effects: np.ndarray,
    prev_allele_freq: Optional[np.ndarray] = None,
) -> GeneticDiagnostics:
    """Compute all genetic summary statistics for a node.

    Args:
        agents: Agent structured array.
        genotypes: (max_agents, N_LOCI, 2) int8.
        effects: (N_ADDITIVE,) float64 effect sizes.
        prev_allele_freq: (N_LOCI,) float64 from previous generation,
            or None if first generation.

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

    # Additive genetic variance
    diag.va_additive = compute_additive_variance(diag.allele_freq, effects)

    # EF1A
    diag.ef1a_allele_freq = float(diag.allele_freq[IDX_EF1A])
    alive_idx = np.where(alive_mask)[0]
    ef1a_het = genotypes[alive_idx, IDX_EF1A, 0] != genotypes[alive_idx, IDX_EF1A, 1]
    diag.ef1a_het_freq = float(np.mean(ef1a_het))

    # Heterozygosity
    diag.heterozygosity_obs, diag.heterozygosity_exp = compute_heterozygosity(
        genotypes, alive_mask
    )

    # Allele frequency change
    if prev_allele_freq is not None:
        delta_q = diag.allele_freq - prev_allele_freq
        diag.delta_q_top3 = delta_q[:3].copy()
        diag.mean_abs_delta_q = float(np.mean(np.abs(delta_q)))

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
    """
    bank: np.ndarray                   # (N_BANK, N_LOCI, 2) int8
    bank_resistance: np.ndarray        # (N_BANK,) float32
    bank_weights: np.ndarray           # (N_BANK,) float64 — frequency weights

    # Population summary statistics
    mean_resistance: float = 0.0
    var_resistance: float = 0.0
    allele_freq: np.ndarray = field(
        default_factory=lambda: np.zeros(N_LOCI, dtype=np.float64)
    )
    heterozygosity: float = 0.0

    def update_summary(self, effects: np.ndarray) -> None:
        """Recompute summary statistics from the genotype bank.

        Args:
            effects: (N_ADDITIVE,) float64 effect sizes.
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
    effects: np.ndarray,
    rng: np.random.Generator,
    n_bank: int = N_BANK,
) -> GenotypeBankState:
    """Compress individual genotypes to a Tier 2 genotype bank.

    Randomly samples ``n_bank`` individuals from the alive population.
    If fewer than ``n_bank`` are alive, pads with duplicates.

    Args:
        genotypes: (max_agents, N_LOCI, 2) int8.
        alive_mask: (max_agents,) bool.
        effects: (N_ADDITIVE,) float64 effect sizes.
        rng: Random generator.
        n_bank: Bank size (default 100).

    Returns:
        GenotypeBankState with summary statistics computed.
    """
    alive_idx = np.where(alive_mask)[0]
    n_alive = len(alive_idx)

    if n_alive == 0:
        # Empty bank
        bank = np.zeros((n_bank, N_LOCI, 2), dtype=np.int8)
        return GenotypeBankState(
            bank=bank,
            bank_resistance=np.zeros(n_bank, dtype=np.float32),
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
    bank_r = np.array(
        [compute_resistance_single(bank_geno[i], effects) for i in range(n_bank)],
        dtype=np.float32,
    )
    weights = np.ones(n_bank, dtype=np.float64) / n_bank

    state = GenotypeBankState(
        bank=bank_geno,
        bank_resistance=bank_r,
        bank_weights=weights,
    )
    state.update_summary(effects)
    return state


def expand_genotype_bank(
    bank: GenotypeBankState,
    n_agents: int,
    effects: np.ndarray,
    rng: np.random.Generator,
    alpha_srs: float = 1.35,
) -> np.ndarray:
    """Expand Tier 2 genotype bank to individual agent genotypes (Tier 1).

    Uses SRS-weighted sampling (NOT uniform) to preserve the heavy-tailed
    genetic structure (CE-8 / spec §10.3).

    Args:
        bank: Tier 2 genotype bank.
        n_agents: Number of individual genotypes to produce.
        effects: (N_ADDITIVE,) float64 effect sizes.
        rng: Random generator.
        alpha_srs: Pareto shape for SRS weighting.

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
    averaged across loci.

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

    For each locus, returns counts of (0/0, 0/1, 1/1) observed and
    expected genotype frequencies.

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
