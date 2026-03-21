"""Inbreeding and genetic diversity metrics for conservation genetics.

Implements equations from Section 6 of the conservation report:
- Genomic inbreeding coefficient F_i (Eq. 6.2)
- Genomic relationship matrix A_ij (Section 6.4/OCS)
- Effective population size N_e (Eq. 6.5)
- Expected heterozygosity H_e (Eq. 6.8)
- Allelic richness (polymorphic loci count)
- Rate of inbreeding ΔF = 1/(2Ne) (Eq. 6.3)
- Projected F after g generations (Eq. 6.4)

Uses sswd_evoepi.types for N_LOCI and trait_slices.
"""

import numpy as np
from typing import Optional

from sswd_evoepi.types import N_LOCI


# ═══════════════════════════════════════════════════════════════════════
# GENOMIC INBREEDING COEFFICIENT (Eq. 6.2)
# ═══════════════════════════════════════════════════════════════════════


def genomic_inbreeding(
    genotypes: np.ndarray,
    ref_freqs: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Compute genomic inbreeding coefficient for each individual.

    F_i = 1 - H_obs,i / H_exp  (Eq. 6.2)

    where H_obs,i is individual i's observed heterozygosity fraction,
    and H_exp is the expected heterozygosity from the reference
    allele frequencies (population-level or supplied).

    Args:
        genotypes: (N, N_LOCI, 2) int8 genotype array.
        ref_freqs: (N_LOCI,) reference allele frequencies for H_exp.
            If None, computed from the input genotypes.

    Returns:
        (N,) float64 inbreeding coefficients. Positive = excess
        homozygosity; negative = excess heterozygosity (outbred).
    """
    n = len(genotypes)
    if n == 0:
        return np.array([], dtype=np.float64)

    # Per-individual observed heterozygosity (fraction of het loci)
    het_per_ind = np.mean(genotypes[:, :, 0] != genotypes[:, :, 1], axis=1)  # (N,)

    # Expected heterozygosity from reference frequencies
    if ref_freqs is None:
        ref_freqs = genotypes.sum(axis=(0, 2)).astype(np.float64) / (2.0 * n)

    he_per_locus = 2.0 * ref_freqs * (1.0 - ref_freqs)  # (N_LOCI,)
    h_exp = float(np.mean(he_per_locus))

    if h_exp < 1e-10:
        return np.zeros(n, dtype=np.float64)

    f = 1.0 - het_per_ind / h_exp
    return f.astype(np.float64)


def mean_inbreeding(
    genotypes: np.ndarray,
    ref_freqs: Optional[np.ndarray] = None,
) -> float:
    """Population mean inbreeding coefficient.

    Args:
        genotypes: (N, N_LOCI, 2) int8.
        ref_freqs: Optional reference allele frequencies.

    Returns:
        Mean F across all individuals.
    """
    f = genomic_inbreeding(genotypes, ref_freqs)
    if len(f) == 0:
        return 0.0
    return float(np.mean(f))


# ═══════════════════════════════════════════════════════════════════════
# GENOMIC RELATIONSHIP MATRIX (Section 5.2.4 / 6.4)
# ═══════════════════════════════════════════════════════════════════════


def genomic_relationship_matrix(
    genotypes: np.ndarray,
    ref_freqs: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Compute the genomic relationship matrix (GRM).

    A_ij = (1/L) Σ_ℓ (x_iℓ - 2q_ℓ)(x_jℓ - 2q_ℓ) / [2q_ℓ(1-q_ℓ)]

    where x_iℓ = a_{i,ℓ,1} + a_{i,ℓ,2} is the allele dosage,
    and q_ℓ is the reference allele frequency.

    This is the VanRaden (2008) Method 1 GRM, standard in genomic
    selection and OCS (Eq. 5.10 constraint).

    Args:
        genotypes: (N, N_LOCI, 2) int8.
        ref_freqs: (N_LOCI,) reference allele frequencies.
            If None, computed from the input genotypes.

    Returns:
        (N, N) float64 relationship matrix. Diagonal ≈ 1 + F_i.
    """
    n = len(genotypes)
    if n == 0:
        return np.empty((0, 0), dtype=np.float64)

    # Allele dosage matrix: (N, N_LOCI)
    X = genotypes.sum(axis=2).astype(np.float64)  # 0, 1, or 2

    if ref_freqs is None:
        ref_freqs = X.mean(axis=0) / 2.0

    # Scaling denominator per locus
    denom = 2.0 * ref_freqs * (1.0 - ref_freqs)
    # Mask monomorphic loci
    poly = denom > 1e-10
    n_poly = poly.sum()

    if n_poly == 0:
        return np.eye(n, dtype=np.float64)

    # Center by 2q
    Z = X[:, poly] - 2.0 * ref_freqs[poly][np.newaxis, :]
    # Scale
    Z = Z / np.sqrt(denom[poly])[np.newaxis, :]

    A = Z @ Z.T / n_poly
    return A


def pairwise_kinship(
    genotypes: np.ndarray,
    ref_freqs: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Pairwise kinship coefficients (half the GRM off-diagonal).

    f_ij = A_ij / 2  for i ≠ j.

    Args:
        genotypes: (N, N_LOCI, 2) int8.
        ref_freqs: Optional reference allele frequencies.

    Returns:
        (N, N) float64 kinship matrix.
    """
    A = genomic_relationship_matrix(genotypes, ref_freqs)
    return A / 2.0


# ═══════════════════════════════════════════════════════════════════════
# EFFECTIVE POPULATION SIZE (Eq. 6.5)
# ═══════════════════════════════════════════════════════════════════════


def ne_from_family_variance(
    offspring_counts: np.ndarray,
) -> float:
    """Estimate N_e from variance in family sizes.

    N_e = (4N - 4) / (σ²_k + 2)  (Eq. 6.5)

    where σ²_k is the variance in offspring number per parent.

    Args:
        offspring_counts: (N_parents,) int — number of offspring per parent.

    Returns:
        Effective population size estimate.
    """
    n = len(offspring_counts)
    if n < 2:
        return float(n)

    vk = float(np.var(offspring_counts, ddof=1))
    ne = (4.0 * n - 4.0) / (vk + 2.0)
    return max(ne, 1.0)


def ne_from_sex_ratio(
    n_males: int,
    n_females: int,
) -> float:
    """Effective size from unequal sex ratio.

    N_e = 4 N_m N_f / (N_m + N_f)  (Eq. 6.4)

    Args:
        n_males: Number of breeding males.
        n_females: Number of breeding females.

    Returns:
        Effective population size.
    """
    total = n_males + n_females
    if total == 0:
        return 0.0
    return 4.0 * n_males * n_females / total


def ne_from_grm(
    genotypes: np.ndarray,
    ref_freqs: Optional[np.ndarray] = None,
) -> float:
    """Estimate N_e from the genomic relationship matrix.

    Uses the mean off-diagonal of the GRM:
    N_e ≈ 1 / (2 × mean(A_ij for i≠j))

    This is a rough estimate; more suitable for large samples.

    Args:
        genotypes: (N, N_LOCI, 2) int8.
        ref_freqs: Optional reference allele frequencies.

    Returns:
        Effective population size estimate.
    """
    A = genomic_relationship_matrix(genotypes, ref_freqs)
    n = len(A)
    if n < 2:
        return float(n)

    # Mean off-diagonal
    mask = ~np.eye(n, dtype=bool)
    mean_aij = float(np.mean(A[mask]))

    if mean_aij <= 0:
        return float(n)  # No detectable relatedness

    return 1.0 / (2.0 * mean_aij)


# ═══════════════════════════════════════════════════════════════════════
# RATE OF INBREEDING (Eq. 6.3)
# ═══════════════════════════════════════════════════════════════════════


def delta_f(ne: float) -> float:
    """Rate of inbreeding per generation.

    ΔF = 1 / (2 N_e)  (Eq. 6.3)

    Args:
        ne: Effective population size.

    Returns:
        ΔF per generation.
    """
    if ne <= 0:
        return 1.0
    return 1.0 / (2.0 * ne)


def projected_f(
    ne: float,
    n_generations: int,
    f_init: float = 0.0,
) -> np.ndarray:
    """Project inbreeding coefficient over multiple generations.

    F_g = 1 - (1 - F_0)(1 - ΔF)^g  (Eq. 6.4)

    Args:
        ne: Effective population size.
        n_generations: Number of generations to project.
        f_init: Initial inbreeding coefficient.

    Returns:
        (n_generations+1,) float64 array of F values (gen 0..n_generations).
    """
    df = delta_f(ne)
    gens = np.arange(n_generations + 1)
    f = 1.0 - (1.0 - f_init) * (1.0 - df) ** gens
    return f


# ═══════════════════════════════════════════════════════════════════════
# DIVERSITY METRICS (Eq. 6.8, Section 6.3)
# ═══════════════════════════════════════════════════════════════════════


def expected_heterozygosity(
    genotypes: np.ndarray,
) -> float:
    """Expected heterozygosity averaged across all loci.

    H_e = (1/L) Σ_ℓ 2q_ℓ(1 - q_ℓ)  (Eq. 6.8)

    Args:
        genotypes: (N, N_LOCI, 2) int8.

    Returns:
        Mean expected heterozygosity.
    """
    n = len(genotypes)
    if n < 1:
        return 0.0

    q = genotypes.sum(axis=(0, 2)).astype(np.float64) / (2.0 * n)
    he_per_locus = 2.0 * q * (1.0 - q)
    return float(np.mean(he_per_locus))


def observed_heterozygosity(
    genotypes: np.ndarray,
) -> float:
    """Observed heterozygosity averaged across all loci.

    H_o = (1/L) Σ_ℓ (fraction heterozygous at ℓ)

    Args:
        genotypes: (N, N_LOCI, 2) int8.

    Returns:
        Mean observed heterozygosity.
    """
    n = len(genotypes)
    if n < 1:
        return 0.0

    het_per_locus = np.mean(genotypes[:, :, 0] != genotypes[:, :, 1], axis=0)
    return float(np.mean(het_per_locus))


def allelic_richness(
    genotypes: np.ndarray,
    threshold: float = 0.0,
) -> int:
    """Count of polymorphic loci (allelic richness for biallelic system).

    A locus is polymorphic if 0 < q < 1 (both alleles present).
    With threshold > 0, requires min(q, 1-q) > threshold.

    Args:
        genotypes: (N, N_LOCI, 2) int8.
        threshold: Minimum minor allele frequency to count as polymorphic.

    Returns:
        Number of polymorphic loci.
    """
    n = len(genotypes)
    if n < 1:
        return 0

    q = genotypes.sum(axis=(0, 2)).astype(np.float64) / (2.0 * n)
    maf = np.minimum(q, 1.0 - q)

    return int(np.sum(maf > threshold))


def per_locus_he(
    genotypes: np.ndarray,
) -> np.ndarray:
    """Per-locus expected heterozygosity.

    Args:
        genotypes: (N, N_LOCI, 2) int8.

    Returns:
        (N_LOCI,) float64 expected heterozygosity per locus.
    """
    n = len(genotypes)
    if n < 1:
        return np.zeros(N_LOCI, dtype=np.float64)

    q = genotypes.sum(axis=(0, 2)).astype(np.float64) / (2.0 * n)
    return 2.0 * q * (1.0 - q)


def heterozygosity_loss_per_gen(ne: float) -> float:
    """Expected proportional loss of heterozygosity per generation.

    ΔH/H = -1/(2Ne)

    Args:
        ne: Effective population size.

    Returns:
        Fractional loss per generation (positive number).
    """
    if ne <= 0:
        return 1.0
    return 1.0 / (2.0 * ne)


# ═══════════════════════════════════════════════════════════════════════
# INBREEDING DEPRESSION (Eq. 6.7)
# ═══════════════════════════════════════════════════════════════════════


def inbreeding_depression(
    f: float,
    b_lethal_equiv: float = 6.0,
    baseline_fitness: float = 1.0,
) -> float:
    """Expected fitness under inbreeding depression.

    w̄(F) = w̄(0) × exp(-B × F)  (Eq. 6.7)

    Args:
        f: Inbreeding coefficient.
        b_lethal_equiv: Lethal equivalents per diploid genome.
            Marine invertebrates: B ∈ [2, 12], default 6.
        baseline_fitness: Fitness at F=0.

    Returns:
        Expected fitness.
    """
    return baseline_fitness * np.exp(-b_lethal_equiv * f)


def generations_to_f_threshold(
    ne: float,
    f_threshold: float = 0.10,
    f_init: float = 0.0,
) -> int:
    """Number of generations to reach a given inbreeding threshold.

    Solves: F_g = 1 - (1-F_0)(1-ΔF)^g ≥ F_threshold

    Args:
        ne: Effective population size.
        f_threshold: Target inbreeding level (default 10%).
        f_init: Starting inbreeding.

    Returns:
        Number of generations (int). Returns 0 if already above threshold,
        -1 if threshold is unreachable (ne → ∞).
    """
    if f_init >= f_threshold:
        return 0

    df = delta_f(ne)
    if df < 1e-15:
        return -1  # Effectively infinite Ne

    # g = ln((1-F_threshold)/(1-F_init)) / ln(1-ΔF)
    ratio = (1.0 - f_threshold) / (1.0 - f_init)
    if ratio <= 0:
        return 0
    g = np.log(ratio) / np.log(1.0 - df)
    return int(np.ceil(g))
