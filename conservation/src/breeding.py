"""Breeding program simulation for conservation genetics.

Implements crossing strategies and multi-generation breeding programs
from Section 5 of the conservation report:
- Mendelian inheritance at 51 diploid loci (Eq. 5.6–5.9)
- Selection schemes: truncation, assortative, complementary, OCS (Eq. 5.10)
- Multi-generation simulator with per-generation tracking
- Within-family selection (Eq. 5.14)
- Fitness weighting across three traits (Eq. 5.1–5.3)

Uses the actual model genetics code from sswd_evoepi.genetics.
"""

import numpy as np
from typing import Optional, Literal
from dataclasses import dataclass, field

from sswd_evoepi.types import N_LOCI, trait_slices
from sswd_evoepi.genetics import (
    compute_trait_batch,
    compute_trait_single,
    initialize_trait_effect_sizes,
    apply_mutations,
    compute_allele_frequencies,
    compute_heterozygosity,
    compute_additive_variance,
)


# ═══════════════════════════════════════════════════════════════════════
# FITNESS WEIGHTING (Eq. 5.1–5.3)
# ═══════════════════════════════════════════════════════════════════════


def selection_index(
    r: np.ndarray,
    t: np.ndarray,
    c: np.ndarray,
    w_r: float = 1.0,
    w_t: float = 0.0,
    w_c: float = 0.0,
) -> np.ndarray:
    """Compute selection index as weighted sum of trait values.

    w(r, t, c) = w_r·r + w_t·t + w_c·c  (linear index)

    Default weights reflect the finding (Section 5.1) that resistance
    is ~1000× more important than tolerance and ~12× more important
    than recovery at population-mean trait values.

    Args:
        r: (N,) resistance scores.
        t: (N,) tolerance scores.
        c: (N,) recovery scores.
        w_r: Weight on resistance (default 1.0).
        w_t: Weight on tolerance (default 0.0).
        w_c: Weight on recovery (default 0.0).

    Returns:
        (N,) selection index values.
    """
    return w_r * r + w_t * t + w_c * c


# ═══════════════════════════════════════════════════════════════════════
# MENDELIAN CROSSING (Eq. 5.6–5.9)
# ═══════════════════════════════════════════════════════════════════════


def mendelian_cross(
    parent1_geno: np.ndarray,
    parent2_geno: np.ndarray,
    rng: np.random.Generator,
    n_offspring: int = 1,
) -> np.ndarray:
    """Produce offspring via Mendelian segregation at all 51 loci.

    At each locus, each parent transmits one of their two alleles
    chosen uniformly at random to the offspring (Eq. 5.6).

    Args:
        parent1_geno: (N_LOCI, 2) int8 genotype of parent 1.
        parent2_geno: (N_LOCI, 2) int8 genotype of parent 2.
        rng: NumPy random Generator.
        n_offspring: Number of offspring to produce.

    Returns:
        (n_offspring, N_LOCI, 2) int8 offspring genotypes.
    """
    offspring = np.empty((n_offspring, N_LOCI, 2), dtype=np.int8)

    # For each locus, randomly pick allele 0 or 1 from each parent
    picks_p1 = rng.integers(0, 2, size=(n_offspring, N_LOCI))
    picks_p2 = rng.integers(0, 2, size=(n_offspring, N_LOCI))

    for i in range(n_offspring):
        offspring[i, :, 0] = parent1_geno[np.arange(N_LOCI), picks_p1[i]]
        offspring[i, :, 1] = parent2_geno[np.arange(N_LOCI), picks_p2[i]]

    return offspring


def batch_cross(
    parent_pairs: np.ndarray,
    genotypes: np.ndarray,
    rng: np.random.Generator,
    n_offspring_per_pair: int = 1,
) -> np.ndarray:
    """Perform Mendelian crosses for multiple parent pairs.

    Args:
        parent_pairs: (n_pairs, 2) int indices into genotypes.
        genotypes: (N_individuals, N_LOCI, 2) int8.
        rng: Random generator.
        n_offspring_per_pair: Offspring per pair.

    Returns:
        (n_pairs * n_offspring_per_pair, N_LOCI, 2) int8.
    """
    n_pairs = len(parent_pairs)
    total = n_pairs * n_offspring_per_pair
    offspring = np.empty((total, N_LOCI, 2), dtype=np.int8)

    for k, (i, j) in enumerate(parent_pairs):
        start = k * n_offspring_per_pair
        end = start + n_offspring_per_pair
        offspring[start:end] = mendelian_cross(
            genotypes[i], genotypes[j], rng, n_offspring_per_pair
        )

    return offspring


# ═══════════════════════════════════════════════════════════════════════
# EXPECTED OFFSPRING RESISTANCE (Eq. 5.8)
# ═══════════════════════════════════════════════════════════════════════


def expected_offspring_trait(
    geno_i: np.ndarray,
    geno_j: np.ndarray,
    effects: np.ndarray,
    locus_slice: slice,
) -> float:
    """Expected offspring trait value from crossing parents i and j.

    E[τ_offspring(i,j)] = Σ_ℓ α_ℓ × (g_i,ℓ + g_j,ℓ) / 4  (Eq. 5.8)

    Args:
        geno_i: (N_LOCI, 2) parent i genotype.
        geno_j: (N_LOCI, 2) parent j genotype.
        effects: (n_trait_loci,) effect sizes.
        locus_slice: Slice for this trait's loci.

    Returns:
        Expected offspring trait value.
    """
    g_i = geno_i[locus_slice].sum(axis=1)  # (n_loci,) 0/1/2
    g_j = geno_j[locus_slice].sum(axis=1)
    q_offspring = (g_i + g_j) / 4.0
    return float(np.dot(effects, q_offspring))


def segregation_variance(
    geno_i: np.ndarray,
    geno_j: np.ndarray,
    effects: np.ndarray,
    locus_slice: slice,
) -> float:
    """Mendelian segregation variance within a cross.

    V_seg = Σ_ℓ (α_ℓ²/4) × h_ℓ^(i) × h_ℓ^(j)  (Eq. 5.9)

    where h_ℓ^(p) = 1 if parent p is heterozygous at locus ℓ.
    Note: also includes variance from homozygous-different combinations.

    Args:
        geno_i: (N_LOCI, 2) parent i genotype.
        geno_j: (N_LOCI, 2) parent j genotype.
        effects: (n_trait_loci,) effect sizes.
        locus_slice: Slice for trait loci.

    Returns:
        Segregation variance.
    """
    # Per-locus variance contribution from each parent
    # Parent contributes variance = 0 if homo, 1/4 if het
    g_i = geno_i[locus_slice]  # (n_loci, 2)
    g_j = geno_j[locus_slice]

    het_i = (g_i[:, 0] != g_i[:, 1]).astype(np.float64)
    het_j = (g_j[:, 0] != g_j[:, 1]).astype(np.float64)

    # Full segregation variance from the cross
    # Each parent independently contributes allele variance of 1/4 when het
    var_per_locus = effects**2 * (het_i + het_j) / 4.0
    return float(np.sum(var_per_locus))


# ═══════════════════════════════════════════════════════════════════════
# PAIRING STRATEGIES (Section 5.2)
# ═══════════════════════════════════════════════════════════════════════

SelectionScheme = Literal[
    "truncation", "assortative", "complementary", "ocs"
]


def truncation_select(
    scores: np.ndarray,
    n_select: int,
) -> np.ndarray:
    """Select the top n_select individuals by score (truncation selection).

    Args:
        scores: (N,) selection index values.
        n_select: Number to select.

    Returns:
        (n_select,) indices of selected individuals.
    """
    n_select = min(n_select, len(scores))
    return np.argsort(scores)[-n_select:][::-1]


def pair_random(
    selected: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    """Random mating: shuffle selected and pair sequentially (Eq. 5.5).

    Args:
        selected: (n,) indices of selected individuals.
        rng: Random generator.

    Returns:
        (n//2, 2) parent pairs.
    """
    perm = rng.permutation(selected)
    n_pairs = len(perm) // 2
    return perm[:2 * n_pairs].reshape(n_pairs, 2)


def pair_assortative(
    selected: np.ndarray,
    scores: np.ndarray,
) -> np.ndarray:
    """Assortative mating: pair best with best (Section 5.2.2).

    Rank selected by score, pair 1st↔2nd, 3rd↔4th, etc.

    Args:
        selected: (n,) indices of selected individuals.
        scores: (N,) trait scores for all individuals.

    Returns:
        (n//2, 2) parent pairs, sorted by average score descending.
    """
    ranked = selected[np.argsort(scores[selected])[::-1]]
    n_pairs = len(ranked) // 2
    return ranked[:2 * n_pairs].reshape(n_pairs, 2)


def pair_complementary(
    selected: np.ndarray,
    genotypes: np.ndarray,
    effects: np.ndarray,
    locus_slice: slice,
) -> np.ndarray:
    """Complementary mating: pair parents to maximize locus coverage (Section 5.2.3).

    Greedy algorithm: for each unmatched individual, find the partner
    that maximizes the expected offspring trait value (Eq. 5.8).

    Args:
        selected: (n,) indices of selected individuals.
        genotypes: (N_total, N_LOCI, 2) int8.
        effects: (n_trait_loci,) effect sizes.
        locus_slice: Slice for the primary trait.

    Returns:
        (n//2, 2) parent pairs.
    """
    remaining = list(selected)
    pairs = []

    while len(remaining) >= 2:
        # Take first remaining individual
        focal = remaining.pop(0)
        best_partner = None
        best_score = -np.inf

        for k, candidate in enumerate(remaining):
            score = expected_offspring_trait(
                genotypes[focal], genotypes[candidate],
                effects, locus_slice,
            )
            if score > best_score:
                best_score = score
                best_partner = k

        pairs.append((focal, remaining.pop(best_partner)))

    return np.array(pairs, dtype=np.intp) if pairs else np.empty((0, 2), dtype=np.intp)


def pair_ocs(
    selected: np.ndarray,
    scores: np.ndarray,
    genotypes: np.ndarray,
    target_ne: float = 50.0,
    rng: Optional[np.random.Generator] = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Optimal Contribution Selection: maximize gain subject to ΔF constraint (Eq. 5.10).

    Simplified version: solve for contribution vector c that maximizes
    c^T τ subject to c^T A c ≤ 1/(2 Ne*), where A is the genomic
    relationship matrix.

    Uses Lagrangian relaxation with iterative reweighting.

    Args:
        selected: (n,) indices of selected individuals.
        scores: (N,) trait scores for all individuals.
        genotypes: (N_total, N_LOCI, 2) int8.
        target_ne: Target effective population size.
        rng: Random generator (for stochastic pairing from contributions).

    Returns:
        Tuple of:
          contributions: (n,) float contribution weights (sum to 1).
          pairs: (n//2, 2) parent pairs sampled proportional to contributions.
    """
    n = len(selected)
    if n < 2:
        c = np.ones(n, dtype=np.float64)
        return c, np.empty((0, 2), dtype=np.intp)

    # Genomic relationship matrix for selected individuals
    sel_geno = genotypes[selected]  # (n, N_LOCI, 2)
    allele_sums = sel_geno.sum(axis=2).astype(np.float64)  # (n, N_LOCI)

    # Center and scale
    q_mean = allele_sums.mean(axis=0) / 2.0  # (N_LOCI,)
    denom = 2.0 * q_mean * (1.0 - q_mean)
    denom = np.where(denom < 1e-10, 1.0, denom)
    Z = (allele_sums / 2.0 - q_mean[np.newaxis, :]) / np.sqrt(denom)[np.newaxis, :]
    A = Z @ Z.T / N_LOCI  # (n, n) genomic relationship matrix

    tau = scores[selected].astype(np.float64)

    # Solve via Lagrangian: c* = A^{-1}(τ - λ1) / normalization
    # Use pseudoinverse for numerical stability
    delta_f_target = 1.0 / (2.0 * target_ne)

    try:
        A_inv = np.linalg.pinv(A + np.eye(n) * 1e-6)
    except np.linalg.LinAlgError:
        # Fallback: equal contributions
        c = np.ones(n, dtype=np.float64) / n
        if rng is not None:
            pairs = _sample_pairs_from_contributions(selected, c, rng)
        else:
            pairs = pair_random(selected, np.random.default_rng(42))
        return c, pairs

    # Binary search for Lagrange multiplier λ
    lam_low, lam_high = 0.0, float(tau.max()) * 10
    for _ in range(50):
        lam = (lam_low + lam_high) / 2.0
        c_raw = A_inv @ (tau - lam)
        c_raw = np.maximum(c_raw, 0.0)
        c_sum = c_raw.sum()
        if c_sum < 1e-10:
            lam_high = lam
            continue
        c = c_raw / c_sum
        delta_f = float(c @ A @ c)
        if delta_f > delta_f_target:
            lam_low = lam
        else:
            lam_high = lam

    # Normalize final contributions
    c_raw = A_inv @ (tau - lam_high)
    c_raw = np.maximum(c_raw, 0.0)
    c_sum = c_raw.sum()
    if c_sum < 1e-10:
        c = np.ones(n, dtype=np.float64) / n
    else:
        c = c_raw / c_sum

    if rng is not None:
        pairs = _sample_pairs_from_contributions(selected, c, rng)
    else:
        pairs = pair_random(selected, np.random.default_rng(42))

    return c, pairs


def _sample_pairs_from_contributions(
    selected: np.ndarray,
    contributions: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    """Sample parent pairs proportional to contribution weights.

    Each individual appears at most once per pair. Sampling without
    replacement, weighted by contributions.

    Args:
        selected: (n,) individual indices.
        contributions: (n,) normalized weights.
        rng: Random generator.

    Returns:
        (n//2, 2) parent pairs.
    """
    n = len(selected)
    n_pairs = n // 2

    # Sample 2*n_pairs indices without replacement, weighted
    prob = contributions / contributions.sum()
    try:
        drawn = rng.choice(n, size=2 * n_pairs, replace=False, p=prob)
    except ValueError:
        # If contributions are too concentrated, fall back to top-weighted
        drawn = np.argsort(prob)[-2 * n_pairs:]
        rng.shuffle(drawn)

    pairs = np.array(
        [(selected[drawn[2 * k]], selected[drawn[2 * k + 1]]) for k in range(n_pairs)],
        dtype=np.intp,
    )
    return pairs


# ═══════════════════════════════════════════════════════════════════════
# WITHIN-FAMILY SELECTION (Eq. 5.14)
# ═══════════════════════════════════════════════════════════════════════


def within_family_select(
    parent_pairs: np.ndarray,
    genotypes: np.ndarray,
    effects: np.ndarray,
    locus_slice: slice,
    rng: np.random.Generator,
    n_offspring_per_pair: int = 100,
    n_keep_per_family: int = 1,
) -> np.ndarray:
    """Within-family selection: produce many offspring, keep the best.

    Leverages high fecundity of sea stars (Eq. 5.14):
    E[τ_(m)] = E[τ_offspring] + σ_seg × E[Z_(m)]

    Args:
        parent_pairs: (n_pairs, 2) parent indices.
        genotypes: (N, N_LOCI, 2) int8.
        effects: (n_trait_loci,) effect sizes for the selected trait.
        locus_slice: Slice for the trait loci.
        rng: Random generator.
        n_offspring_per_pair: Offspring produced per cross.
        n_keep_per_family: Best offspring to keep from each family.

    Returns:
        (n_pairs * n_keep_per_family, N_LOCI, 2) int8 selected offspring.
    """
    n_pairs = len(parent_pairs)
    kept = []

    for i, j in parent_pairs:
        offspring = mendelian_cross(
            genotypes[i], genotypes[j], rng, n_offspring_per_pair
        )
        # Score all offspring
        alive_mask = np.ones(n_offspring_per_pair, dtype=bool)
        scores = compute_trait_batch(offspring, effects, alive_mask, locus_slice)
        # Keep the best
        best_idx = np.argsort(scores)[-n_keep_per_family:][::-1]
        kept.append(offspring[best_idx])

    if kept:
        return np.concatenate(kept, axis=0)
    return np.empty((0, N_LOCI, 2), dtype=np.int8)


# ═══════════════════════════════════════════════════════════════════════
# PER-GENERATION STATISTICS
# ═══════════════════════════════════════════════════════════════════════


@dataclass
class GenerationStats:
    """Statistics tracked per generation of a breeding program."""
    generation: int = 0
    n_individuals: int = 0

    # Trait means
    mean_r: float = 0.0
    mean_t: float = 0.0
    mean_c: float = 0.0

    # Trait maxima
    max_r: float = 0.0
    max_t: float = 0.0
    max_c: float = 0.0

    # Additive genetic variance
    va_r: float = 0.0
    va_t: float = 0.0
    va_c: float = 0.0

    # Diversity
    he: float = 0.0  # Expected heterozygosity
    ho: float = 0.0  # Observed heterozygosity
    f_mean: float = 0.0  # Mean inbreeding coefficient

    # Loci fixed (for primary trait only)
    loci_fixed: int = 0

    # Selection index
    mean_index: float = 0.0
    max_index: float = 0.0


def compute_generation_stats(
    genotypes: np.ndarray,
    effects_r: np.ndarray,
    effects_t: np.ndarray,
    effects_c: np.ndarray,
    generation: int,
    res_slice: slice,
    tol_slice: slice,
    rec_slice: slice,
    w_r: float = 1.0,
    w_t: float = 0.0,
    w_c: float = 0.0,
) -> GenerationStats:
    """Compute all tracking statistics for one generation.

    Args:
        genotypes: (N, N_LOCI, 2) int8 — all individuals in this generation.
        effects_r: Resistance effect sizes.
        effects_t: Tolerance effect sizes.
        effects_c: Recovery effect sizes.
        generation: Generation number.
        res_slice: Resistance locus slice.
        tol_slice: Tolerance locus slice.
        rec_slice: Recovery locus slice.
        w_r: Selection index weight on resistance.
        w_t: Selection index weight on tolerance.
        w_c: Selection index weight on recovery.

    Returns:
        GenerationStats with all fields populated.
    """
    n = len(genotypes)
    if n == 0:
        return GenerationStats(generation=generation)

    alive = np.ones(n, dtype=bool)

    # Trait scores
    r = compute_trait_batch(genotypes, effects_r, alive, res_slice)
    t = compute_trait_batch(genotypes, effects_t, alive, tol_slice)
    c = compute_trait_batch(genotypes, effects_c, alive, rec_slice)

    # Allele frequencies
    freq = compute_allele_frequencies_from_array(genotypes)

    # Additive variance
    va_r = compute_additive_variance(freq, effects_r, res_slice)
    va_t = compute_additive_variance(freq, effects_t, tol_slice)
    va_c = compute_additive_variance(freq, effects_c, rec_slice)

    # Heterozygosity
    ho, he = _heterozygosity_from_array(genotypes)

    # Mean inbreeding from excess homozygosity (Eq. 6.2)
    if he > 0:
        f_mean = 1.0 - ho / he
    else:
        f_mean = 0.0

    # Loci fixed for resistance
    res_freq = freq[res_slice]
    loci_fixed = int(np.sum((res_freq < 0.001) | (res_freq > 0.999)))

    # Selection index
    idx = selection_index(r, t, c, w_r, w_t, w_c)

    return GenerationStats(
        generation=generation,
        n_individuals=n,
        mean_r=float(np.mean(r)),
        mean_t=float(np.mean(t)),
        mean_c=float(np.mean(c)),
        max_r=float(np.max(r)),
        max_t=float(np.max(t)),
        max_c=float(np.max(c)),
        va_r=va_r,
        va_t=va_t,
        va_c=va_c,
        he=he,
        ho=ho,
        f_mean=f_mean,
        loci_fixed=loci_fixed,
        mean_index=float(np.mean(idx)),
        max_index=float(np.max(idx)),
    )


# ═══════════════════════════════════════════════════════════════════════
# HELPERS (genotype arrays without agent structured arrays)
# ═══════════════════════════════════════════════════════════════════════


def compute_allele_frequencies_from_array(genotypes: np.ndarray) -> np.ndarray:
    """Compute allele frequencies from a plain genotype array.

    Args:
        genotypes: (N, N_LOCI, 2) int8.

    Returns:
        (N_LOCI,) float64 allele frequencies.
    """
    n = len(genotypes)
    if n == 0:
        return np.zeros(N_LOCI, dtype=np.float64)
    return genotypes.sum(axis=(0, 2)).astype(np.float64) / (2.0 * n)


def _heterozygosity_from_array(genotypes: np.ndarray) -> tuple[float, float]:
    """Compute Ho and He from a plain genotype array.

    Returns:
        (Ho, He) averaged across all loci.
    """
    n = len(genotypes)
    if n < 2:
        return 0.0, 0.0

    het_per_locus = np.mean(genotypes[:, :, 0] != genotypes[:, :, 1], axis=0)
    ho = float(np.mean(het_per_locus))

    q = genotypes.sum(axis=2).astype(np.float64).mean(axis=0) / 2.0
    he_per_locus = 2.0 * q * (1.0 - q)
    he = float(np.mean(he_per_locus))

    return ho, he


# ═══════════════════════════════════════════════════════════════════════
# MULTI-GENERATION BREEDING SIMULATOR (Section 5.3, Eq. 5.12)
# ═══════════════════════════════════════════════════════════════════════


@dataclass
class BreedingResult:
    """Complete results from a multi-generation breeding program."""
    stats: list  # List[GenerationStats]
    final_genotypes: np.ndarray  # (N, N_LOCI, 2) int8
    scheme: str
    n_founders: int
    n_selected: int
    n_offspring_per_pair: int
    n_generations: int


def run_breeding_program(
    founder_genotypes: np.ndarray,
    effects_r: np.ndarray,
    effects_t: np.ndarray,
    effects_c: np.ndarray,
    rng: np.random.Generator,
    n_generations: int = 5,
    n_selected: int = 20,
    n_offspring_per_pair: int = 100,
    n_keep_per_family: int = 2,
    scheme: SelectionScheme = "truncation",
    n_r: int = 17,
    n_t: int = 17,
    n_c: int = 17,
    w_r: float = 1.0,
    w_t: float = 0.0,
    w_c: float = 0.0,
    mu: float = 1e-8,
    target_ne: float = 50.0,
) -> BreedingResult:
    """Run a multi-generation breeding program.

    Each generation:
    1. Score all individuals by selection index
    2. Select top n_selected (truncation on index)
    3. Pair selected parents according to scheme
    4. Produce n_offspring_per_pair offspring per pair
    5. Optionally apply within-family selection
    6. Apply mutations
    7. Record statistics

    Args:
        founder_genotypes: (N_founders, N_LOCI, 2) int8.
        effects_r: Resistance effect sizes.
        effects_t: Tolerance effect sizes.
        effects_c: Recovery effect sizes.
        rng: Random generator.
        n_generations: Number of breeding generations.
        n_selected: Parents selected per generation.
        n_offspring_per_pair: Offspring per cross.
        n_keep_per_family: Best kept per family (within-family selection).
        scheme: Pairing strategy.
        n_r: Number of resistance loci.
        n_t: Number of tolerance loci.
        n_c: Number of recovery loci.
        w_r: Selection index weight on resistance.
        w_t: Selection index weight on tolerance.
        w_c: Selection index weight on recovery.
        mu: Mutation rate.
        target_ne: Target Ne for OCS.

    Returns:
        BreedingResult with per-generation stats and final genotypes.
    """
    res_s, tol_s, rec_s = trait_slices(n_r, n_t, n_c)
    current_geno = founder_genotypes.copy()
    stats_history = []

    # Record founder generation (gen 0)
    gen0 = compute_generation_stats(
        current_geno, effects_r, effects_t, effects_c,
        generation=0,
        res_slice=res_s, tol_slice=tol_s, rec_slice=rec_s,
        w_r=w_r, w_t=w_t, w_c=w_c,
    )
    stats_history.append(gen0)

    for g in range(1, n_generations + 1):
        n_pop = len(current_geno)
        alive = np.ones(n_pop, dtype=bool)

        # Score everyone
        r = compute_trait_batch(current_geno, effects_r, alive, res_s)
        t = compute_trait_batch(current_geno, effects_t, alive, tol_s)
        c = compute_trait_batch(current_geno, effects_c, alive, rec_s)
        idx = selection_index(r, t, c, w_r, w_t, w_c)

        # Select
        actual_n_sel = min(n_selected, n_pop)
        selected = truncation_select(idx, actual_n_sel)

        # Pair
        if scheme == "assortative":
            pairs = pair_assortative(selected, idx)
        elif scheme == "complementary":
            pairs = pair_complementary(
                selected, current_geno, effects_r, res_s,
            )
        elif scheme == "ocs":
            _, pairs = pair_ocs(
                selected, idx, current_geno,
                target_ne=target_ne, rng=rng,
            )
        else:  # truncation (random pairing of selected)
            pairs = pair_random(selected, rng)

        if len(pairs) == 0:
            break

        # Produce offspring with within-family selection
        if n_keep_per_family < n_offspring_per_pair:
            offspring = within_family_select(
                pairs, current_geno, effects_r, res_s, rng,
                n_offspring_per_pair=n_offspring_per_pair,
                n_keep_per_family=n_keep_per_family,
            )
        else:
            offspring = batch_cross(
                pairs, current_geno, rng,
                n_offspring_per_pair=n_offspring_per_pair,
            )

        # Mutations
        if mu > 0:
            apply_mutations(offspring, rng, mu)

        current_geno = offspring

        # Record stats
        gen_stats = compute_generation_stats(
            current_geno, effects_r, effects_t, effects_c,
            generation=g,
            res_slice=res_s, tol_slice=tol_s, rec_slice=rec_s,
            w_r=w_r, w_t=w_t, w_c=w_c,
        )
        stats_history.append(gen_stats)

    return BreedingResult(
        stats=stats_history,
        final_genotypes=current_geno,
        scheme=scheme,
        n_founders=len(founder_genotypes),
        n_selected=n_selected,
        n_offspring_per_pair=n_offspring_per_pair,
        n_generations=n_generations,
    )
