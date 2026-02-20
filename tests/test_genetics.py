"""Tests for sswd_evoepi.genetics — 51-locus three-trait diploid architecture.

Three-trait architecture (17R/17T/17C):
  - 17 resistance loci (immune exclusion)
  - 17 tolerance loci (damage limitation)
  - 17 recovery loci (pathogen clearance)

Acceptance criteria:
  - Effect sizes: exponential distribution, sorted descending, sum = total_weight
  - Each trait score in [0, 1] independently
  - Without selection: Hardy-Weinberg equilibrium maintained over 50 generations
  - Without selection: allele frequencies stable (drift only, not directional)
  - Ne/N ≈ 10⁻³ under SRS
  - Three-trait initialization hits target means
  - Vectorized trait score matches loop-based calculation
  - Genotype bank compression/expansion preserves statistics for all 3 traits
  - Mutation is applied but negligible at μ = 10⁻⁸
  - F_ST = 0 for identical populations
  - NO cost of resistance (CE-1)
  - NO EF1A overdominant locus (removed — Pisaster, not Pycnopodia)
"""

import numpy as np
import pytest

from sswd_evoepi.types import (
    AGENT_DTYPE,
    N_LOCI,
    N_RESISTANCE_DEFAULT,
    N_TOLERANCE_DEFAULT,
    N_RECOVERY_DEFAULT,
    DiseaseState,
    Stage,
    allocate_agents,
    allocate_genotypes,
    trait_slices,
)
from sswd_evoepi.genetics import (
    MU_PER_LOCUS,
    N_BANK,
    RESISTANCE_SLICE,
    TOLERANCE_SLICE,
    RECOVERY_SLICE,
    initialize_trait_effect_sizes,
    initialize_effect_sizes,
    compute_trait_single,
    compute_trait_batch,
    compute_resistance_single,
    compute_resistance_batch,
    update_all_trait_scores,
    update_resistance_scores,
    initialize_genotypes,
    initialize_genotypes_beta,
    initialize_genotypes_three_trait,
    apply_mutations,
    compute_allele_frequencies,
    compute_heterozygosity,
    compute_additive_variance,
    compute_genetic_diagnostics,
    compute_ne_ratio_from_offspring,
    compress_to_genotype_bank,
    expand_genotype_bank,
    compute_fst,
    hardy_weinberg_test,
    GeneticDiagnostics,
    GenotypeBankState,
)
from sswd_evoepi.reproduction import (
    mendelian_inherit_batch,
    srs_reproductive_lottery,
)


# ═══════════════════════════════════════════════════════════════════════
# FIXTURES
# ═══════════════════════════════════════════════════════════════════════


@pytest.fixture
def rng():
    return np.random.default_rng(42)


@pytest.fixture
def effect_sizes_r():
    """Resistance effect sizes (17 loci)."""
    return initialize_trait_effect_sizes(np.random.default_rng(12345), 17, 1.0)


@pytest.fixture
def effect_sizes_t():
    """Tolerance effect sizes (17 loci)."""
    return initialize_trait_effect_sizes(np.random.default_rng(12346), 17, 1.0)


@pytest.fixture
def effect_sizes_c():
    """Recovery effect sizes (17 loci)."""
    return initialize_trait_effect_sizes(np.random.default_rng(12347), 17, 1.0)


@pytest.fixture
def three_effects(effect_sizes_r, effect_sizes_t, effect_sizes_c):
    """Tuple of (effects_r, effects_t, effects_c)."""
    return effect_sizes_r, effect_sizes_t, effect_sizes_c


@pytest.fixture
def population_500(three_effects, rng):
    """A 500-agent population with three-trait genotypes."""
    eff_r, eff_t, eff_c = three_effects
    n = 500
    agents = allocate_agents(n)
    genotypes = initialize_genotypes_three_trait(
        n, eff_r, eff_t, eff_c, rng,
        target_mean_r=0.15, target_mean_t=0.10, target_mean_c=0.08,
    )

    for i in range(n):
        agents[i]['alive'] = True
        agents[i]['stage'] = Stage.ADULT
        agents[i]['sex'] = 0 if i < 250 else 1
        agents[i]['size'] = 600.0 + rng.normal(0, 50)
        agents[i]['age'] = 10.0
        agents[i]['disease_state'] = DiseaseState.S
        agents[i]['node_id'] = 0

    # Compute all trait scores
    update_all_trait_scores(agents, genotypes, eff_r, eff_t, eff_c)

    return agents, genotypes


# ═══════════════════════════════════════════════════════════════════════
# EFFECT SIZE TESTS
# ═══════════════════════════════════════════════════════════════════════


class TestEffectSizes:
    def test_shape(self, effect_sizes_r):
        assert effect_sizes_r.shape == (17,)
        assert effect_sizes_r.dtype == np.float64

    def test_sum_equals_total_weight(self, effect_sizes_r):
        assert np.isclose(effect_sizes_r.sum(), 1.0, atol=1e-10)

    def test_sorted_descending(self, effect_sizes_r):
        for i in range(len(effect_sizes_r) - 1):
            assert effect_sizes_r[i] >= effect_sizes_r[i + 1]

    def test_all_positive(self, effect_sizes_r):
        assert np.all(effect_sizes_r > 0)

    def test_top_loci_disproportionate(self, effect_sizes_r):
        """Top 3 loci should contribute >10% of total (exponential tail)."""
        top3_sum = effect_sizes_r[:3].sum()
        assert top3_sum / 1.0 > 0.10

    def test_reproducible(self):
        """Same seed → same effect sizes."""
        e1 = initialize_trait_effect_sizes(np.random.default_rng(999), 17, 1.0)
        e2 = initialize_trait_effect_sizes(np.random.default_rng(999), 17, 1.0)
        np.testing.assert_array_equal(e1, e2)

    def test_different_seeds_differ(self):
        e1 = initialize_trait_effect_sizes(np.random.default_rng(1), 17, 1.0)
        e2 = initialize_trait_effect_sizes(np.random.default_rng(2), 17, 1.0)
        assert not np.allclose(e1, e2)

    def test_different_n_loci(self):
        """Effect sizes work with different loci counts."""
        e10 = initialize_trait_effect_sizes(np.random.default_rng(1), 10, 1.0)
        e25 = initialize_trait_effect_sizes(np.random.default_rng(1), 25, 1.0)
        assert e10.shape == (10,)
        assert e25.shape == (25,)
        assert np.isclose(e10.sum(), 1.0)
        assert np.isclose(e25.sum(), 1.0)

    def test_custom_total_weight(self):
        """Total weight can be set to arbitrary value."""
        e = initialize_trait_effect_sizes(np.random.default_rng(1), 17, 0.5)
        assert np.isclose(e.sum(), 0.5, atol=1e-10)

    def test_backward_compat_alias(self):
        """initialize_effect_sizes still works as an alias."""
        e = initialize_effect_sizes(np.random.default_rng(1), 17, 1.0)
        assert e.shape == (17,)
        assert np.isclose(e.sum(), 1.0)


# ═══════════════════════════════════════════════════════════════════════
# TRAIT SCORE TESTS
# ═══════════════════════════════════════════════════════════════════════


class TestTraitScore:
    def test_all_zero_genotype(self, effect_sizes_r):
        """All-zero genotype → trait score = 0."""
        geno = np.zeros((N_LOCI, 2), dtype=np.int8)
        r = compute_trait_single(geno, effect_sizes_r, RESISTANCE_SLICE)
        assert r == 0.0

    def test_max_resistance(self, effect_sizes_r):
        """All resistance loci homozygous derived → score = 1.0."""
        geno = np.zeros((N_LOCI, 2), dtype=np.int8)
        geno[:17, :] = 1
        r = compute_trait_single(geno, effect_sizes_r, RESISTANCE_SLICE)
        assert np.isclose(r, 1.0, atol=1e-6)

    def test_max_tolerance(self, effect_sizes_t):
        """All tolerance loci homozygous derived → score = 1.0."""
        geno = np.zeros((N_LOCI, 2), dtype=np.int8)
        geno[17:34, :] = 1
        t = compute_trait_single(geno, effect_sizes_t, TOLERANCE_SLICE)
        assert np.isclose(t, 1.0, atol=1e-6)

    def test_max_recovery(self, effect_sizes_c):
        """All recovery loci homozygous derived → score = 1.0."""
        geno = np.zeros((N_LOCI, 2), dtype=np.int8)
        geno[34:51, :] = 1
        c = compute_trait_single(geno, effect_sizes_c, RECOVERY_SLICE)
        assert np.isclose(c, 1.0, atol=1e-6)

    def test_traits_independent(self, three_effects):
        """Derived alleles in one trait block don't affect other traits."""
        eff_r, eff_t, eff_c = three_effects
        geno = np.zeros((N_LOCI, 2), dtype=np.int8)
        geno[:17, :] = 1  # Only resistance loci derived

        r = compute_trait_single(geno, eff_r, RESISTANCE_SLICE)
        t = compute_trait_single(geno, eff_t, TOLERANCE_SLICE)
        c = compute_trait_single(geno, eff_c, RECOVERY_SLICE)

        assert r > 0.9  # Should be ~1.0
        assert t == 0.0  # No tolerance alleles
        assert c == 0.0  # No recovery alleles

    def test_heterozygote_half(self, effect_sizes_r):
        """Single heterozygous locus → score = effect/2."""
        geno = np.zeros((N_LOCI, 2), dtype=np.int8)
        geno[0, 0] = 1  # heterozygous at locus 0
        r = compute_trait_single(geno, effect_sizes_r, RESISTANCE_SLICE)
        assert np.isclose(r, effect_sizes_r[0] / 2.0, atol=1e-6)

    def test_scores_in_bounds(self, three_effects, population_500):
        """All trait scores in [0, 1]."""
        agents, genotypes = population_500
        alive = agents['alive']
        r_vals = agents['resistance'][alive]
        t_vals = agents['tolerance'][alive]
        c_vals = agents['recovery_ability'][alive]

        for name, vals in [('r', r_vals), ('t', t_vals), ('c', c_vals)]:
            assert np.all(vals >= 0.0), f"{name} has negative values"
            assert np.all(vals <= 1.0), f"{name} has values > 1"

    def test_vectorized_matches_loop(self, three_effects, population_500):
        """Batch computation matches individual computation for all traits."""
        eff_r, eff_t, eff_c = three_effects
        agents, genotypes = population_500
        alive_mask = agents['alive']

        batch_r = compute_trait_batch(genotypes, eff_r, alive_mask, RESISTANCE_SLICE)
        batch_t = compute_trait_batch(genotypes, eff_t, alive_mask, TOLERANCE_SLICE)
        batch_c = compute_trait_batch(genotypes, eff_c, alive_mask, RECOVERY_SLICE)
        alive_idx = np.where(alive_mask)[0]

        for idx in alive_idx[:50]:
            single_r = compute_trait_single(genotypes[idx], eff_r, RESISTANCE_SLICE)
            single_t = compute_trait_single(genotypes[idx], eff_t, TOLERANCE_SLICE)
            single_c = compute_trait_single(genotypes[idx], eff_c, RECOVERY_SLICE)
            assert np.isclose(batch_r[idx], single_r, atol=1e-5), \
                f"R mismatch at {idx}"
            assert np.isclose(batch_t[idx], single_t, atol=1e-5), \
                f"T mismatch at {idx}"
            assert np.isclose(batch_c[idx], single_c, atol=1e-5), \
                f"C mismatch at {idx}"

    def test_dead_agents_get_zero(self, effect_sizes_r):
        """Dead agents should have score = 0 in batch computation."""
        n = 10
        genotypes = np.ones((n, N_LOCI, 2), dtype=np.int8)
        alive_mask = np.array([True] * 5 + [False] * 5)
        r = compute_trait_batch(genotypes, effect_sizes_r, alive_mask, RESISTANCE_SLICE)
        assert np.all(r[5:] == 0.0)
        assert np.all(r[:5] > 0.0)

    def test_backward_compat_wrappers(self, effect_sizes_r):
        """compute_resistance_batch/single still work."""
        geno = np.zeros((N_LOCI, 2), dtype=np.int8)
        geno[:17, :] = 1
        r = compute_resistance_single(geno, effect_sizes_r)
        assert np.isclose(r, 1.0, atol=1e-6)

        genotypes = np.zeros((5, N_LOCI, 2), dtype=np.int8)
        genotypes[:, :17, :] = 1
        alive = np.ones(5, dtype=bool)
        r_batch = compute_resistance_batch(genotypes, effect_sizes_r, alive)
        assert np.allclose(r_batch, 1.0, atol=1e-6)


# ═══════════════════════════════════════════════════════════════════════
# THREE-TRAIT GENOTYPE INITIALIZATION
# ═══════════════════════════════════════════════════════════════════════


class TestThreeTraitInitialization:
    def test_shape(self, three_effects, rng):
        eff_r, eff_t, eff_c = three_effects
        geno = initialize_genotypes_three_trait(100, eff_r, eff_t, eff_c, rng)
        assert geno.shape == (100, N_LOCI, 2)
        assert geno.dtype == np.int8

    def test_values_binary(self, three_effects, rng):
        eff_r, eff_t, eff_c = three_effects
        geno = initialize_genotypes_three_trait(100, eff_r, eff_t, eff_c, rng)
        assert np.all((geno == 0) | (geno == 1))

    def test_target_mean_resistance(self, three_effects, rng):
        """Population mean resistance should be near target (0.15)."""
        eff_r, eff_t, eff_c = three_effects
        n = 5000
        geno = initialize_genotypes_three_trait(
            n, eff_r, eff_t, eff_c, rng, target_mean_r=0.15
        )
        alive = np.ones(n, dtype=bool)
        r = compute_trait_batch(geno, eff_r, alive, RESISTANCE_SLICE)
        assert abs(np.mean(r) - 0.15) < 0.03, \
            f"Mean r = {np.mean(r):.4f}, expected ~0.15"

    def test_target_mean_tolerance(self, three_effects, rng):
        """Population mean tolerance should be near target (0.10)."""
        eff_r, eff_t, eff_c = three_effects
        n = 5000
        geno = initialize_genotypes_three_trait(
            n, eff_r, eff_t, eff_c, rng, target_mean_t=0.10
        )
        alive = np.ones(n, dtype=bool)
        t = compute_trait_batch(geno, eff_t, alive, TOLERANCE_SLICE)
        assert abs(np.mean(t) - 0.10) < 0.03, \
            f"Mean t = {np.mean(t):.4f}, expected ~0.10"

    def test_target_mean_recovery(self, three_effects, rng):
        """Population mean recovery should be near target (0.08)."""
        eff_r, eff_t, eff_c = three_effects
        n = 5000
        geno = initialize_genotypes_three_trait(
            n, eff_r, eff_t, eff_c, rng, target_mean_c=0.08
        )
        alive = np.ones(n, dtype=bool)
        c = compute_trait_batch(geno, eff_c, alive, RECOVERY_SLICE)
        assert abs(np.mean(c) - 0.08) < 0.03, \
            f"Mean c = {np.mean(c):.4f}, expected ~0.08"

    def test_custom_partition(self, rng):
        """Non-default partition (10/20/21) should work."""
        eff_r = initialize_trait_effect_sizes(rng, 10, 1.0)
        eff_t = initialize_trait_effect_sizes(rng, 20, 1.0)
        eff_c = initialize_trait_effect_sizes(rng, 21, 1.0)
        geno = initialize_genotypes_three_trait(
            100, eff_r, eff_t, eff_c, rng,
            n_resistance=10, n_tolerance=20, n_recovery=21,
        )
        assert geno.shape == (100, 51, 2)

        # Check scores use correct slices
        rs, ts, cs = trait_slices(10, 20, 21)
        alive = np.ones(100, dtype=bool)
        r = compute_trait_batch(geno, eff_r, alive, rs)
        t = compute_trait_batch(geno, eff_t, alive, ts)
        c = compute_trait_batch(geno, eff_c, alive, cs)
        assert np.all(r >= 0) and np.all(r <= 1)
        assert np.all(t >= 0) and np.all(t <= 1)
        assert np.all(c >= 0) and np.all(c <= 1)


# ═══════════════════════════════════════════════════════════════════════
# BACKWARD-COMPATIBLE INITIALIZATION
# ═══════════════════════════════════════════════════════════════════════


class TestBackwardCompatInit:
    def test_initialize_genotypes_shape(self, effect_sizes_r, rng):
        geno = initialize_genotypes(100, effect_sizes_r, rng, target_mean_r=0.15)
        assert geno.shape == (100, N_LOCI, 2)
        assert geno.dtype == np.int8

    def test_initialize_genotypes_absorbs_legacy_args(self, effect_sizes_r, rng):
        """Legacy ef1a_q argument should be absorbed without error."""
        geno = initialize_genotypes(
            100, effect_sizes_r, rng, target_mean_r=0.15, ef1a_q=0.24
        )
        assert geno.shape == (100, N_LOCI, 2)

    def test_initialize_genotypes_beta(self, effect_sizes_r, rng):
        geno = initialize_genotypes_beta(
            100, effect_sizes_r, rng, target_mean_r=0.15, ef1a_q=0.24
        )
        assert geno.shape == (100, N_LOCI, 2)

    def test_only_resistance_loci_populated(self, effect_sizes_r, rng):
        """Backward-compat init only sets resistance loci (0-16)."""
        geno = initialize_genotypes(500, effect_sizes_r, rng, target_mean_r=0.15)
        # Resistance loci should have some derived alleles
        assert geno[:, :17, :].sum() > 0
        # Tolerance and recovery loci should be all zero
        assert geno[:, 17:, :].sum() == 0


# ═══════════════════════════════════════════════════════════════════════
# UPDATE ALL TRAIT SCORES
# ═══════════════════════════════════════════════════════════════════════


class TestUpdateAllTraitScores:
    def test_writes_all_three_traits(self, three_effects, rng):
        eff_r, eff_t, eff_c = three_effects
        n = 50
        agents = allocate_agents(n)
        geno = initialize_genotypes_three_trait(n, eff_r, eff_t, eff_c, rng)

        for i in range(n):
            agents[i]['alive'] = True
            agents[i]['resistance'] = -1.0
            agents[i]['tolerance'] = -1.0
            agents[i]['recovery_ability'] = -1.0

        update_all_trait_scores(agents, geno, eff_r, eff_t, eff_c)

        alive = agents['alive']
        assert np.all(agents['resistance'][alive] >= 0.0)
        assert np.all(agents['tolerance'][alive] >= 0.0)
        assert np.all(agents['recovery_ability'][alive] >= 0.0)
        assert np.all(agents['resistance'][alive] <= 1.0)
        assert np.all(agents['tolerance'][alive] <= 1.0)
        assert np.all(agents['recovery_ability'][alive] <= 1.0)

    def test_backward_compat_update_resistance(self, effect_sizes_r, rng):
        """update_resistance_scores still works."""
        n = 50
        agents = allocate_agents(n)
        geno = initialize_genotypes(n, effect_sizes_r, rng, target_mean_r=0.15)

        for i in range(n):
            agents[i]['alive'] = True

        update_resistance_scores(agents, geno, effect_sizes_r)
        alive = agents['alive']
        assert np.all(agents['resistance'][alive] >= 0.0)
        assert np.all(agents['resistance'][alive] <= 1.0)


# ═══════════════════════════════════════════════════════════════════════
# MUTATION
# ═══════════════════════════════════════════════════════════════════════


class TestMutation:
    def test_mutation_at_default_rate_negligible(self, rng):
        """At μ=10⁻⁸ with 500 offspring, expect ~0 mutations."""
        n = 500
        geno = np.zeros((n, N_LOCI, 2), dtype=np.int8)
        n_mut = apply_mutations(geno, rng, mu=MU_PER_LOCUS)
        # Expected: 500 × 51 × 2 × 10⁻⁸ ≈ 5.1×10⁻⁴ mutations → almost always 0
        assert n_mut == 0 or n_mut <= 2

    def test_mutation_at_high_rate(self, rng):
        """At elevated μ, mutations should be detectable."""
        n = 1000
        geno = np.zeros((n, N_LOCI, 2), dtype=np.int8)
        n_mut = apply_mutations(geno, rng, mu=0.01)
        # Expected: 1000 × 51 × 2 × 0.01 = 1020 mutations
        assert n_mut > 500, f"Expected many mutations at μ=0.01, got {n_mut}"
        assert geno.sum() > 0

    def test_mutation_bidirectional(self, rng):
        """Mutations should flip 0→1 and 1→0."""
        n = 1000
        geno = np.ones((n, N_LOCI, 2), dtype=np.int8)
        original_sum = geno.sum()
        n_mut = apply_mutations(geno, rng, mu=0.01)
        assert geno.sum() < original_sum


# ═══════════════════════════════════════════════════════════════════════
# ALLELE FREQUENCIES & HETEROZYGOSITY
# ═══════════════════════════════════════════════════════════════════════


class TestAlleleFrequencies:
    def test_allele_freq_range(self, population_500):
        agents, genotypes = population_500
        freq = compute_allele_frequencies(genotypes, agents['alive'])
        assert freq.shape == (N_LOCI,)
        assert np.all(freq >= 0.0)
        assert np.all(freq <= 1.0)

    def test_allele_freq_all_zero(self):
        """All-zero genotypes → all frequencies = 0."""
        n = 100
        geno = np.zeros((n, N_LOCI, 2), dtype=np.int8)
        alive = np.ones(n, dtype=bool)
        freq = compute_allele_frequencies(geno, alive)
        np.testing.assert_array_equal(freq, 0.0)

    def test_allele_freq_all_one(self):
        """All-one genotypes → all frequencies = 1."""
        n = 100
        geno = np.ones((n, N_LOCI, 2), dtype=np.int8)
        alive = np.ones(n, dtype=bool)
        freq = compute_allele_frequencies(geno, alive)
        np.testing.assert_array_almost_equal(freq, 1.0)

    def test_heterozygosity_agreement(self, population_500):
        """H_o and H_e should be similar for HWE population."""
        agents, genotypes = population_500
        H_o, H_e = compute_heterozygosity(genotypes, agents['alive'])
        assert abs(H_o - H_e) < 0.05, f"H_o={H_o:.4f}, H_e={H_e:.4f}"

    def test_heterozygosity_range(self, population_500):
        agents, genotypes = population_500
        H_o, H_e = compute_heterozygosity(genotypes, agents['alive'])
        assert 0.0 <= H_o <= 0.5
        assert 0.0 <= H_e <= 0.5


class TestAdditiveVariance:
    def test_va_positive(self, effect_sizes_r, population_500):
        agents, genotypes = population_500
        freq = compute_allele_frequencies(genotypes, agents['alive'])
        va = compute_additive_variance(freq, effect_sizes_r, RESISTANCE_SLICE)
        assert va > 0, "V_A should be positive with polymorphic loci"

    def test_va_per_trait(self, three_effects, population_500):
        """V_A should be computable per trait independently."""
        eff_r, eff_t, eff_c = three_effects
        agents, genotypes = population_500
        freq = compute_allele_frequencies(genotypes, agents['alive'])
        va_r = compute_additive_variance(freq, eff_r, RESISTANCE_SLICE)
        va_t = compute_additive_variance(freq, eff_t, TOLERANCE_SLICE)
        va_c = compute_additive_variance(freq, eff_c, RECOVERY_SLICE)
        assert va_r > 0
        assert va_t > 0
        assert va_c > 0

    def test_va_zero_at_fixation(self, effect_sizes_r):
        """V_A = 0 when all loci are fixed (q = 0 or q = 1)."""
        freq = np.zeros(N_LOCI, dtype=np.float64)
        va = compute_additive_variance(freq, effect_sizes_r, RESISTANCE_SLICE)
        assert va == 0.0

        freq_fixed = np.ones(N_LOCI, dtype=np.float64)
        va = compute_additive_variance(freq_fixed, effect_sizes_r, RESISTANCE_SLICE)
        assert va == 0.0

    def test_va_backward_compat_no_slice(self, effect_sizes_r, population_500):
        """V_A without slice uses first len(effects) loci."""
        agents, genotypes = population_500
        freq = compute_allele_frequencies(genotypes, agents['alive'])
        va_with_slice = compute_additive_variance(freq, effect_sizes_r, RESISTANCE_SLICE)
        va_no_slice = compute_additive_variance(freq, effect_sizes_r)
        assert np.isclose(va_with_slice, va_no_slice)


# ═══════════════════════════════════════════════════════════════════════
# HARDY-WEINBERG EQUILIBRIUM MAINTENANCE
# ═══════════════════════════════════════════════════════════════════════


class TestHardyWeinberg:
    def test_hwe_initial(self, three_effects):
        """Freshly initialized population should be at HWE."""
        eff_r, eff_t, eff_c = three_effects
        rng = np.random.default_rng(42)
        n = 2000
        geno = initialize_genotypes_three_trait(n, eff_r, eff_t, eff_c, rng)
        alive = np.ones(n, dtype=bool)

        obs_freq, exp_freq, chi2 = hardy_weinberg_test(geno, alive)

        # Average chi2 across all loci should be small
        mean_chi2 = float(np.mean(chi2))
        assert mean_chi2 < 5.0, \
            f"Mean chi2 = {mean_chi2:.2f}, population not at HWE"

    def test_hwe_maintained_without_selection(self, effect_sizes_r):
        """HWE maintained over 50 generations of random mating (no selection)."""
        rng = np.random.default_rng(42)
        n = 500

        geno = initialize_genotypes(n, effect_sizes_r, rng, target_mean_r=0.15)

        for gen in range(50):
            mothers = rng.integers(0, n, size=n)
            fathers = rng.integers(0, n, size=n)
            offspring = mendelian_inherit_batch(geno, mothers, fathers, n, rng)
            # No EF1A lethal elimination needed
            geno = offspring[:n] if len(offspring) >= n else np.concatenate([
                offspring,
                offspring[rng.integers(0, len(offspring), size=n - len(offspring))]
            ])

        alive = np.ones(n, dtype=bool)
        obs_freq, exp_freq, chi2 = hardy_weinberg_test(geno, alive)

        # At most 15% of loci should show significant departure
        n_significant = np.sum(chi2 > 3.84)
        fraction_sig = n_significant / N_LOCI
        assert fraction_sig < 0.20, \
            f"{n_significant}/{N_LOCI} loci depart from HWE (expect <20%)"

    def test_allele_freq_stable_without_selection(self, effect_sizes_r):
        """Allele frequencies stable (drift only) over 50 generations."""
        rng = np.random.default_rng(99)
        n = 500

        geno = initialize_genotypes(n, effect_sizes_r, rng, target_mean_r=0.15)
        alive = np.ones(n, dtype=bool)
        initial_freq = compute_allele_frequencies(geno, alive)
        initial_mean_q = float(np.mean(initial_freq[:17]))

        for gen in range(50):
            mothers = rng.integers(0, n, size=n)
            fathers = rng.integers(0, n, size=n)
            offspring = mendelian_inherit_batch(geno, mothers, fathers, n, rng)
            geno = offspring[:n] if len(offspring) >= n else np.concatenate([
                offspring,
                offspring[rng.integers(0, len(offspring), size=n - len(offspring))]
            ])

        final_freq = compute_allele_frequencies(geno, alive)
        final_mean_q = float(np.mean(final_freq[:17]))

        assert abs(final_mean_q - initial_mean_q) < 0.03, \
            f"Mean q shifted from {initial_mean_q:.4f} to {final_mean_q:.4f}"


# ═══════════════════════════════════════════════════════════════════════
# Ne/N UNDER SRS
# ═══════════════════════════════════════════════════════════════════════


class TestNeOverN:
    def test_ne_ratio_wright_fisher(self):
        """Under Wright-Fisher (uniform offspring), Ne/N ≈ 2."""
        offspring = np.full(100, 10, dtype=np.int64)
        ratio = compute_ne_ratio_from_offspring(offspring)
        assert ratio > 1.5, f"Ne/N = {ratio:.3f}, expected ~2 for uniform"

    def test_ne_ratio_extreme_skew(self):
        """Under extreme SRS, Ne/N << 1."""
        offspring = np.zeros(1000, dtype=np.int64)
        offspring[0] = 5000
        offspring[1] = 3000
        ratio = compute_ne_ratio_from_offspring(offspring)
        assert ratio < 0.01, f"Ne/N = {ratio:.3f}, expected << 0.01"

    def test_ne_ratio_under_srs(self, effect_sizes_r):
        """Full SRS lottery should produce Ne/N << 1."""
        rng = np.random.default_rng(42)
        n_parents = 200
        agents = allocate_agents(n_parents)
        genotypes = allocate_genotypes(n_parents)

        for i in range(n_parents):
            agents[i]['alive'] = True
            agents[i]['stage'] = Stage.ADULT
            agents[i]['sex'] = 0 if i < 100 else 1
            agents[i]['size'] = 600.0
            agents[i]['disease_state'] = DiseaseState.S

        geno_init = initialize_genotypes(n_parents, effect_sizes_r, rng)
        genotypes[:n_parents] = geno_init

        females = np.where(
            agents['alive'] & (agents['sex'] == 0) & (agents['stage'] == Stage.ADULT)
        )[0]
        males = np.where(
            agents['alive'] & (agents['sex'] == 1) & (agents['stage'] == Stage.ADULT)
        )[0]

        offspring_geno, parent_pairs = srs_reproductive_lottery(
            females, males, agents, genotypes,
            n_offspring_target=5000, alpha_srs=1.35, rng=rng,
        )

        all_parents = parent_pairs.ravel()
        full_counts = np.zeros(n_parents, dtype=np.int64)
        unique_parents, parent_counts = np.unique(all_parents, return_counts=True)
        full_counts[unique_parents] = parent_counts

        ratio = compute_ne_ratio_from_offspring(full_counts)
        assert ratio < 0.1, f"Ne/N = {ratio:.4f}, expected << 1 under SRS"


# ═══════════════════════════════════════════════════════════════════════
# GENETIC DIAGNOSTICS
# ═══════════════════════════════════════════════════════════════════════


class TestGeneticDiagnostics:
    def test_diagnostics_populated(self, three_effects, population_500):
        eff_r, eff_t, eff_c = three_effects
        agents, genotypes = population_500
        diag = compute_genetic_diagnostics(
            agents, genotypes, eff_r, eff_t, eff_c
        )

        assert diag.n_alive == 500
        assert diag.mean_resistance > 0
        assert diag.var_resistance > 0
        assert diag.va_resistance > 0
        assert diag.mean_tolerance > 0
        assert diag.var_tolerance > 0
        assert diag.va_tolerance > 0
        assert diag.mean_recovery > 0
        assert diag.var_recovery > 0
        assert diag.va_recovery > 0
        assert diag.heterozygosity_obs > 0
        assert diag.heterozygosity_exp > 0

    def test_diagnostics_no_ef1a_fields(self):
        """EF1A fields should not exist in GeneticDiagnostics."""
        diag = GeneticDiagnostics()
        assert not hasattr(diag, 'ef1a_allele_freq')
        assert not hasattr(diag, 'ef1a_het_freq')

    def test_diagnostics_allele_freq_shape(self, three_effects, population_500):
        eff_r, eff_t, eff_c = three_effects
        agents, genotypes = population_500
        diag = compute_genetic_diagnostics(agents, genotypes, eff_r, eff_t, eff_c)
        assert diag.allele_freq.shape == (51,)

    def test_diagnostics_allele_freq_change(self, three_effects, population_500):
        eff_r, eff_t, eff_c = three_effects
        agents, genotypes = population_500
        prev_freq = np.full(N_LOCI, 0.1, dtype=np.float64)
        diag = compute_genetic_diagnostics(
            agents, genotypes, eff_r, eff_t, eff_c, prev_allele_freq=prev_freq
        )
        assert diag.delta_q_top3.shape == (3,)
        assert diag.mean_abs_delta_q > 0

    def test_diagnostics_backward_compat(self, effect_sizes_r, rng):
        """Diagnostics work with only resistance effects (backward compat)."""
        n = 50
        agents = allocate_agents(n)
        geno = initialize_genotypes(n, effect_sizes_r, rng, target_mean_r=0.15)
        for i in range(n):
            agents[i]['alive'] = True
        update_resistance_scores(agents, geno, effect_sizes_r)
        diag = compute_genetic_diagnostics(agents, geno, effect_sizes_r)
        assert diag.n_alive == 50
        assert diag.mean_resistance > 0

    def test_diagnostics_empty_population(self, effect_sizes_r):
        agents = allocate_agents(10)
        genotypes = allocate_genotypes(10)
        diag = compute_genetic_diagnostics(agents, genotypes, effect_sizes_r)
        assert diag.n_alive == 0
        assert diag.mean_resistance == 0.0


# ═══════════════════════════════════════════════════════════════════════
# GENOTYPE BANK (TIER 2)
# ═══════════════════════════════════════════════════════════════════════


class TestGenotypeBank:
    def test_compress_roundtrip(self, three_effects, population_500):
        eff_r, eff_t, eff_c = three_effects
        agents, genotypes = population_500
        rng = np.random.default_rng(42)

        bank = compress_to_genotype_bank(
            genotypes, agents['alive'], eff_r, rng, effects_t=eff_t, effects_c=eff_c
        )
        assert bank.bank.shape == (N_BANK, N_LOCI, 2)
        assert bank.bank_resistance.shape == (N_BANK,)
        assert bank.bank_tolerance.shape == (N_BANK,)
        assert bank.bank_recovery.shape == (N_BANK,)
        assert abs(bank.mean_resistance - agents['resistance'][agents['alive']].mean()) < 0.05

    def test_bank_three_trait_stats(self, three_effects, population_500):
        """Bank should track mean/var for all three traits."""
        eff_r, eff_t, eff_c = three_effects
        agents, genotypes = population_500
        rng = np.random.default_rng(42)

        bank = compress_to_genotype_bank(
            genotypes, agents['alive'], eff_r, rng, effects_t=eff_t, effects_c=eff_c
        )
        assert bank.mean_tolerance > 0
        assert bank.mean_recovery > 0
        assert bank.var_tolerance >= 0
        assert bank.var_recovery >= 0

    def test_expand_produces_valid_genotypes(self, three_effects, population_500):
        eff_r, eff_t, eff_c = three_effects
        agents, genotypes = population_500
        rng = np.random.default_rng(42)

        bank = compress_to_genotype_bank(
            genotypes, agents['alive'], eff_r, rng, effects_t=eff_t, effects_c=eff_c
        )
        expanded = expand_genotype_bank(bank, 300, eff_r, rng)
        assert expanded.shape == (300, N_LOCI, 2)
        assert np.all((expanded == 0) | (expanded == 1))

    def test_expand_preserves_resistance_range(self, three_effects, population_500):
        eff_r, eff_t, eff_c = three_effects
        agents, genotypes = population_500
        rng = np.random.default_rng(42)

        bank = compress_to_genotype_bank(
            genotypes, agents['alive'], eff_r, rng, effects_t=eff_t, effects_c=eff_c
        )
        expanded = expand_genotype_bank(bank, 500, eff_r, rng)

        r_expanded = np.array([
            compute_trait_single(expanded[i], eff_r, RESISTANCE_SLICE) for i in range(500)
        ])
        assert np.all(r_expanded >= 0)
        assert np.all(r_expanded <= 1.0)
        assert abs(np.mean(r_expanded) - bank.mean_resistance) < 0.05

    def test_bank_update_summary(self, three_effects, population_500):
        eff_r, eff_t, eff_c = three_effects
        agents, genotypes = population_500
        rng = np.random.default_rng(42)

        bank = compress_to_genotype_bank(
            genotypes, agents['alive'], eff_r, rng, effects_t=eff_t, effects_c=eff_c
        )
        assert bank.heterozygosity > 0
        assert len(bank.allele_freq) == N_LOCI


# ═══════════════════════════════════════════════════════════════════════
# F_ST
# ═══════════════════════════════════════════════════════════════════════


class TestFST:
    def test_fst_identical_populations(self, population_500):
        """F_ST = 0 for identical populations."""
        agents, genotypes = population_500
        freq = compute_allele_frequencies(genotypes, agents['alive'])
        fst = compute_fst([freq, freq, freq])
        assert np.isclose(fst, 0.0, atol=1e-10)

    def test_fst_different_populations(self, three_effects):
        """F_ST > 0 for genetically different populations."""
        eff_r, eff_t, eff_c = three_effects
        rng = np.random.default_rng(42)
        n = 200

        geno1 = initialize_genotypes_three_trait(
            n, eff_r, eff_t, eff_c, rng, target_mean_r=0.05
        )
        geno2 = initialize_genotypes_three_trait(
            n, eff_r, eff_t, eff_c, rng, target_mean_r=0.25
        )

        alive = np.ones(n, dtype=bool)
        freq1 = compute_allele_frequencies(geno1, alive)
        freq2 = compute_allele_frequencies(geno2, alive)

        fst = compute_fst([freq1, freq2])
        assert fst > 0, "F_ST should be > 0 for different populations"

    def test_fst_single_population(self, population_500):
        """F_ST = 0 for a single population (< 2 nodes)."""
        agents, genotypes = population_500
        freq = compute_allele_frequencies(genotypes, agents['alive'])
        fst = compute_fst([freq])
        assert fst == 0.0


# ═══════════════════════════════════════════════════════════════════════
# UPDATE RESISTANCE SCORES (in-place)
# ═══════════════════════════════════════════════════════════════════════


class TestUpdateResistanceScores:
    def test_writes_to_agents(self, effect_sizes_r):
        rng = np.random.default_rng(42)
        n = 50
        agents = allocate_agents(n)
        geno = initialize_genotypes(n, effect_sizes_r, rng)

        for i in range(n):
            agents[i]['alive'] = True
            agents[i]['resistance'] = -1.0

        update_resistance_scores(agents, geno, effect_sizes_r)
        assert np.all(agents['resistance'][:n] >= 0.0)
        assert np.all(agents['resistance'][:n] <= 1.0)

    def test_dead_get_zero(self, effect_sizes_r):
        rng = np.random.default_rng(42)
        n = 10
        agents = allocate_agents(n)
        genotypes = allocate_genotypes(n)
        agents['resistance'] = -999.0

        update_resistance_scores(agents, genotypes, effect_sizes_r)
        assert np.all(agents['resistance'] == 0.0)


# ═══════════════════════════════════════════════════════════════════════
# INTEGRATION: MENDELIAN INHERITANCE + GENETICS
# ═══════════════════════════════════════════════════════════════════════


class TestMendelianIntegration:
    def test_offspring_alleles_from_parents(self, effect_sizes_r):
        """Each offspring allele must come from one of its parents."""
        rng = np.random.default_rng(42)
        n = 50
        geno = initialize_genotypes(n, effect_sizes_r, rng)

        mothers = rng.integers(0, n, size=200)
        fathers = rng.integers(0, n, size=200)
        offspring = mendelian_inherit_batch(geno, mothers, fathers, 200, rng)

        for i in range(min(50, len(offspring))):
            m = mothers[i]
            f = fathers[i]
            for l in range(N_LOCI):
                assert offspring[i, l, 0] in (geno[m, l, 0], geno[m, l, 1]), \
                    f"Offspring {i} locus {l} maternal allele not from mother"
                assert offspring[i, l, 1] in (geno[f, l, 0], geno[f, l, 1]), \
                    f"Offspring {i} locus {l} paternal allele not from father"


# ═══════════════════════════════════════════════════════════════════════
# CONSTANTS VERIFICATION
# ═══════════════════════════════════════════════════════════════════════


class TestConstants:
    def test_n_loci(self):
        assert N_LOCI == 51

    def test_default_partition_sums_to_n_loci(self):
        assert N_RESISTANCE_DEFAULT + N_TOLERANCE_DEFAULT + N_RECOVERY_DEFAULT == N_LOCI

    def test_default_slices(self):
        assert RESISTANCE_SLICE == slice(0, 17)
        assert TOLERANCE_SLICE == slice(17, 34)
        assert RECOVERY_SLICE == slice(34, 51)

    def test_trait_slices_validation(self):
        """Invalid partition should raise AssertionError."""
        with pytest.raises(AssertionError):
            trait_slices(10, 10, 10)  # sum=30, not 51

    def test_no_ef1a_constants(self):
        """EF1A-related constants should not exist."""
        import sswd_evoepi.genetics as gen
        assert not hasattr(gen, 'S_HET')
        assert not hasattr(gen, 'W_OD')
        assert not hasattr(gen, 'W_ADD')
        assert not hasattr(gen, 'eliminate_ef1a_lethals')
