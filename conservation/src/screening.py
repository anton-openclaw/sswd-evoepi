"""Screening theory for founder selection.

Implements equations from Section 4 of the conservation report:
- Required sample size for exceedance (Eq. 4.1)
- Empirical exceedance from genotype arrays
- Expected maximum from sample (order statistics, Eq. 4.6)
- Multi-site optimal allocation (Eq. 4.8)
- Complementarity scoring between individuals (Eq. 4.9–4.10)

Uses sswd_evoepi.genetics for trait score computation.
"""

import numpy as np
from scipy import stats, optimize
from typing import Optional

from sswd_evoepi.types import N_LOCI, trait_slices
from sswd_evoepi.genetics import (
    compute_trait_batch,
    compute_trait_single,
)


# ═══════════════════════════════════════════════════════════════════════
# REQUIRED SAMPLE SIZE (Eq. 4.1)
# ═══════════════════════════════════════════════════════════════════════


def required_sample_size(
    exceedance_prob: float,
    confidence: float = 0.95,
) -> int:
    """Minimum sample size to find ≥1 individual above threshold.

    n(γ, p) = ⌈ln(1 - γ) / ln(1 - p)⌉  (Eq. 4.1)

    Args:
        exceedance_prob: P(τ ≥ τ*) — probability a random individual
            exceeds the threshold.
        confidence: Desired probability of finding at least one (default 0.95).

    Returns:
        Required sample size (int ≥ 1).

    Raises:
        ValueError: If exceedance_prob ≤ 0 or ≥ 1, or confidence not in (0,1).
    """
    if exceedance_prob <= 0:
        raise ValueError(f"exceedance_prob must be > 0, got {exceedance_prob}")
    if exceedance_prob >= 1:
        return 1
    if not (0 < confidence < 1):
        raise ValueError(f"confidence must be in (0, 1), got {confidence}")

    return int(np.ceil(
        np.log(1 - confidence) / np.log(1 - exceedance_prob)
    ))


def required_n_for_k(
    exceedance_prob: float,
    k: int = 1,
    confidence: float = 0.95,
) -> int:
    """Sample size to find at least k individuals above threshold.

    Uses the inverse of the binomial CDF:
    Find smallest n such that P(X ≥ k | n, p) ≥ γ.

    Args:
        exceedance_prob: P(τ ≥ τ*).
        k: Minimum number of qualifying individuals.
        confidence: Desired probability.

    Returns:
        Required sample size.
    """
    if exceedance_prob <= 0:
        raise ValueError(f"exceedance_prob must be > 0, got {exceedance_prob}")
    if exceedance_prob >= 1:
        return k

    # Binary search for n
    n_low, n_high = k, max(k, int(10 * k / exceedance_prob))
    # Ensure upper bound is sufficient
    while 1 - stats.binom.cdf(k - 1, n_high, exceedance_prob) < confidence:
        n_high *= 2

    while n_low < n_high:
        n_mid = (n_low + n_high) // 2
        prob_at_least_k = 1 - stats.binom.cdf(k - 1, n_mid, exceedance_prob)
        if prob_at_least_k >= confidence:
            n_high = n_mid
        else:
            n_low = n_mid + 1

    return n_low


# ═══════════════════════════════════════════════════════════════════════
# EMPIRICAL EXCEEDANCE FROM GENOTYPES
# ═══════════════════════════════════════════════════════════════════════


def empirical_exceedance(
    genotypes: np.ndarray,
    effects: np.ndarray,
    locus_slice: slice,
    threshold: float,
) -> float:
    """Empirical exceedance probability from a genotype array.

    P̂(τ ≥ τ*) = #{i : τ_i ≥ τ*} / N

    Args:
        genotypes: (N, N_LOCI, 2) int8.
        effects: (n_trait_loci,) effect sizes.
        locus_slice: Slice for the trait's loci.
        threshold: Trait value threshold τ*.

    Returns:
        Empirical exceedance probability.
    """
    n = len(genotypes)
    if n == 0:
        return 0.0

    alive = np.ones(n, dtype=bool)
    scores = compute_trait_batch(genotypes, effects, alive, locus_slice)
    return float(np.sum(scores >= threshold)) / n


def exceedance_curve(
    genotypes: np.ndarray,
    effects: np.ndarray,
    locus_slice: slice,
    thresholds: Optional[np.ndarray] = None,
    n_points: int = 100,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute the full exceedance curve P(τ ≥ τ*) vs τ*.

    Args:
        genotypes: (N, N_LOCI, 2) int8.
        effects: (n_trait_loci,) effect sizes.
        locus_slice: Slice for the trait's loci.
        thresholds: Explicit threshold array. If None, auto-generated.
        n_points: Number of threshold points (if auto-generated).

    Returns:
        (thresholds, exceedance_probs) — both (n_points,) arrays.
    """
    n = len(genotypes)
    if n == 0:
        t = np.linspace(0, 1, n_points)
        return t, np.zeros(n_points)

    alive = np.ones(n, dtype=bool)
    scores = compute_trait_batch(genotypes, effects, alive, locus_slice)

    if thresholds is None:
        lo = float(np.min(scores)) - 0.01
        hi = float(np.max(scores)) + 0.01
        thresholds = np.linspace(lo, hi, n_points)

    probs = np.array([float(np.sum(scores >= t)) / n for t in thresholds])
    return thresholds, probs


# ═══════════════════════════════════════════════════════════════════════
# EXPECTED MAXIMUM FROM SAMPLE (Eq. 4.6)
# ═══════════════════════════════════════════════════════════════════════


def expected_max_normal(
    mu: float,
    sigma: float,
    n: int,
) -> float:
    """Expected maximum trait value from n samples (normal approximation).

    E[τ_(n)] ≈ μ + σ × Φ⁻¹(n/(n+1))  (Eq. 4.6)

    Args:
        mu: Population mean.
        sigma: Population standard deviation.
        n: Sample size.

    Returns:
        Expected maximum.
    """
    if sigma <= 0 or n <= 1:
        return mu
    return mu + sigma * stats.norm.ppf(n / (n + 1))


def expected_max_empirical(
    genotypes: np.ndarray,
    effects: np.ndarray,
    locus_slice: slice,
    n_sample: int,
    rng: np.random.Generator,
    n_replicates: int = 1000,
) -> float:
    """Expected maximum trait value via Monte Carlo sampling.

    More accurate than the normal approximation for skewed distributions.

    Args:
        genotypes: (N_pop, N_LOCI, 2) int8 — the source population.
        effects: Effect sizes.
        locus_slice: Trait locus slice.
        n_sample: Sample size per draw.
        rng: Random generator.
        n_replicates: Number of MC replicates.

    Returns:
        Expected maximum (averaged over replicates).
    """
    n_pop = len(genotypes)
    if n_pop == 0:
        return 0.0

    alive = np.ones(n_pop, dtype=bool)
    all_scores = compute_trait_batch(genotypes, effects, alive, locus_slice)

    maxima = np.empty(n_replicates)
    for rep in range(n_replicates):
        idx = rng.choice(n_pop, size=min(n_sample, n_pop), replace=False)
        maxima[rep] = np.max(all_scores[idx])

    return float(np.mean(maxima))


def screening_effort_curve(
    genotypes: np.ndarray,
    effects: np.ndarray,
    locus_slice: slice,
    sample_sizes: np.ndarray,
    rng: np.random.Generator,
    n_replicates: int = 500,
) -> np.ndarray:
    """Expected max vs sample size (the diminishing returns curve).

    Args:
        genotypes: (N_pop, N_LOCI, 2) int8.
        effects: Effect sizes.
        locus_slice: Trait locus slice.
        sample_sizes: (K,) array of sample sizes to evaluate.
        rng: Random generator.
        n_replicates: MC replicates per sample size.

    Returns:
        (K,) expected maximum at each sample size.
    """
    n_pop = len(genotypes)
    alive = np.ones(n_pop, dtype=bool)
    all_scores = compute_trait_batch(genotypes, effects, alive, locus_slice)

    expected_maxes = np.empty(len(sample_sizes))
    for i, n_s in enumerate(sample_sizes):
        n_s = min(int(n_s), n_pop)
        maxima = np.empty(n_replicates)
        for rep in range(n_replicates):
            idx = rng.choice(n_pop, size=n_s, replace=False)
            maxima[rep] = np.max(all_scores[idx])
        expected_maxes[i] = np.mean(maxima)

    return expected_maxes


# ═══════════════════════════════════════════════════════════════════════
# MULTI-SITE OPTIMAL ALLOCATION (Eq. 4.8)
# ═══════════════════════════════════════════════════════════════════════


def multisite_allocation(
    site_means: np.ndarray,
    site_stds: np.ndarray,
    total_budget: int,
    site_max_n: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Optimal screening allocation across multiple sites.

    Maximizes E[max_k τ_(n_k)^(k)] subject to Σ n_k = N (Eq. 4.8).

    Uses the normal approximation E[max] ≈ μ + σ × Φ⁻¹(n/(n+1)).
    Greedy algorithm: iteratively allocate one sample to the site
    where it produces the largest marginal improvement.

    Args:
        site_means: (K,) mean trait at each site.
        site_stds: (K,) trait SD at each site.
        total_budget: Total samples available.
        site_max_n: (K,) maximum samples per site (population limits).
            If None, no per-site limit.

    Returns:
        (K,) int allocation per site.
    """
    k = len(site_means)
    alloc = np.ones(k, dtype=int)  # Start with 1 per site
    remaining = total_budget - k

    if site_max_n is None:
        site_max_n = np.full(k, total_budget, dtype=int)

    if remaining < 0:
        # Can't even sample every site once — allocate to best sites
        alloc = np.zeros(k, dtype=int)
        ranked = np.argsort(site_means)[::-1]
        for i in range(min(total_budget, k)):
            alloc[ranked[i]] = 1
        return alloc

    def _expected_max(site_idx, n):
        if n <= 0:
            return -np.inf
        return expected_max_normal(
            site_means[site_idx], site_stds[site_idx], n
        )

    # Greedy marginal allocation
    for _ in range(remaining):
        best_site = -1
        best_gain = -np.inf

        for s in range(k):
            if alloc[s] >= site_max_n[s]:
                continue
            current = _expected_max(s, alloc[s])
            proposed = _expected_max(s, alloc[s] + 1)
            gain = proposed - current
            if gain > best_gain:
                best_gain = gain
                best_site = s

        if best_site < 0:
            break

        alloc[best_site] += 1

    return alloc


# ═══════════════════════════════════════════════════════════════════════
# COMPLEMENTARITY SCORING (Eq. 4.9–4.10)
# ═══════════════════════════════════════════════════════════════════════


def locus_union(
    geno_i: np.ndarray,
    geno_j: np.ndarray,
    locus_slice: slice,
) -> int:
    """Count of loci where at least one parent carries a protective allele.

    U(i, j) = Σ_ℓ 1[(g_i,ℓ > 0) ∨ (g_j,ℓ > 0)]  (Eq. 4.9)

    Args:
        geno_i: (N_LOCI, 2) int8 genotype of individual i.
        geno_j: (N_LOCI, 2) int8 genotype of individual j.
        locus_slice: Slice for the trait's loci.

    Returns:
        Number of covered loci (int).
    """
    g_i = geno_i[locus_slice].sum(axis=1)  # (n_loci,) 0/1/2
    g_j = geno_j[locus_slice].sum(axis=1)
    covered = (g_i > 0) | (g_j > 0)
    return int(np.sum(covered))


def locus_overlap(
    geno_i: np.ndarray,
    geno_j: np.ndarray,
    locus_slice: slice,
) -> int:
    """Count of loci where both parents carry protective alleles.

    O(i, j) = Σ_ℓ 1[(g_i,ℓ > 0) ∧ (g_j,ℓ > 0)]

    Args:
        geno_i: (N_LOCI, 2) int8 genotype of individual i.
        geno_j: (N_LOCI, 2) int8 genotype of individual j.
        locus_slice: Slice for the trait's loci.

    Returns:
        Number of overlapping loci (int).
    """
    g_i = geno_i[locus_slice].sum(axis=1)
    g_j = geno_j[locus_slice].sum(axis=1)
    overlap = (g_i > 0) & (g_j > 0)
    return int(np.sum(overlap))


def complementarity_score(
    geno_i: np.ndarray,
    geno_j: np.ndarray,
    locus_slice: slice,
) -> int:
    """Complementarity score between two individuals.

    C(i, j) = U(i, j) - O(i, j)  (Eq. 4.10)

    High C means the parents cover different loci.

    Args:
        geno_i: (N_LOCI, 2) int8 genotype of individual i.
        geno_j: (N_LOCI, 2) int8 genotype of individual j.
        locus_slice: Slice for the trait's loci.

    Returns:
        Complementarity score (int ≥ 0).
    """
    return locus_union(geno_i, geno_j, locus_slice) - locus_overlap(geno_i, geno_j, locus_slice)


def unique_contribution(
    geno_i: np.ndarray,
    geno_j: np.ndarray,
    locus_slice: slice,
) -> tuple[int, int]:
    """Unique locus contributions of each individual in a pair.

    Args:
        geno_i: (N_LOCI, 2) int8 genotype of individual i.
        geno_j: (N_LOCI, 2) int8 genotype of individual j.
        locus_slice: Slice for the trait's loci.

    Returns:
        (unique_i, unique_j) — loci covered by only i, only j.
    """
    g_i = geno_i[locus_slice].sum(axis=1)
    g_j = geno_j[locus_slice].sum(axis=1)

    has_i = g_i > 0
    has_j = g_j > 0

    unique_i = int(np.sum(has_i & ~has_j))
    unique_j = int(np.sum(~has_i & has_j))
    return unique_i, unique_j


def complementarity_matrix(
    genotypes: np.ndarray,
    locus_slice: slice,
    metric: str = "union",
) -> np.ndarray:
    """Pairwise complementarity/union matrix for a set of individuals.

    Args:
        genotypes: (N, N_LOCI, 2) int8.
        locus_slice: Trait locus slice.
        metric: "union" (Eq. 4.9), "complementarity" (Eq. 4.10),
            or "overlap".

    Returns:
        (N, N) int matrix.
    """
    n = len(genotypes)
    mat = np.zeros((n, n), dtype=int)

    # Precompute allele dosages
    dosage = genotypes[:, locus_slice, :].sum(axis=2)  # (N, n_loci)

    for i in range(n):
        for j in range(i, n):
            has_i = dosage[i] > 0
            has_j = dosage[j] > 0

            if metric == "union":
                val = int(np.sum(has_i | has_j))
            elif metric == "overlap":
                val = int(np.sum(has_i & has_j))
            elif metric == "complementarity":
                union = int(np.sum(has_i | has_j))
                overlap = int(np.sum(has_i & has_j))
                val = union - overlap
            else:
                raise ValueError(f"Unknown metric: {metric}")

            mat[i, j] = val
            mat[j, i] = val

    return mat


def select_complementary_founders(
    genotypes: np.ndarray,
    effects: np.ndarray,
    locus_slice: slice,
    n_founders: int,
    w_trait: float = 0.5,
    w_complement: float = 0.5,
) -> np.ndarray:
    """Greedy selection of founders maximizing trait + complementarity.

    At each step, adds the individual that maximizes:
    w_trait × normalized_trait + w_complement × marginal_coverage

    Args:
        genotypes: (N_candidates, N_LOCI, 2) int8.
        effects: Effect sizes for the focal trait.
        locus_slice: Trait locus slice.
        n_founders: Number of founders to select.
        w_trait: Weight on individual trait value.
        w_complement: Weight on marginal locus coverage.

    Returns:
        (n_founders,) int indices of selected founders.
    """
    n = len(genotypes)
    n_founders = min(n_founders, n)

    alive = np.ones(n, dtype=bool)
    scores = compute_trait_batch(genotypes, effects, alive, locus_slice)

    # Normalize scores to [0, 1]
    s_min, s_max = float(scores.min()), float(scores.max())
    if s_max > s_min:
        norm_scores = (scores - s_min) / (s_max - s_min)
    else:
        norm_scores = np.ones(n, dtype=np.float32)

    dosage = genotypes[:, locus_slice, :].sum(axis=2)  # (N, n_loci)
    n_loci = dosage.shape[1]

    selected = []
    covered = np.zeros(n_loci, dtype=bool)
    remaining = set(range(n))

    for _ in range(n_founders):
        best_idx = -1
        best_score = -np.inf

        for idx in remaining:
            has_allele = dosage[idx] > 0
            marginal = float(np.sum(has_allele & ~covered)) / n_loci
            combined = w_trait * norm_scores[idx] + w_complement * marginal
            if combined > best_score:
                best_score = combined
                best_idx = idx

        if best_idx < 0:
            break

        selected.append(best_idx)
        remaining.discard(best_idx)
        covered |= (dosage[best_idx] > 0)

    return np.array(selected, dtype=np.intp)
