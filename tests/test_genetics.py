"""Tests for sswd_evoepi.genetics — 52-locus diploid resistance architecture.

Acceptance criteria (Phase 6):
  - Effect sizes: exponential distribution, sorted descending, sum = W_ADD
  - Without selection: Hardy-Weinberg equilibrium maintained over 50 generations
  - Without selection: allele frequencies stable (drift only, not directional)
  - Ne/N ≈ 10⁻³ under SRS
  - EF1A homozygous lethal: ins/ins individuals eliminated
  - EF1A maintained near equilibrium frequency
  - r_i distribution is bimodal (ERRATA E6)
  - Resistance scores in [0, 1]
  - Vectorized r_i matches loop-based calculation
  - Genotype bank compression/expansion preserves statistics
  - Mutation is applied but negligible at μ = 10⁻⁸
  - F_ST = 0 for identical populations
  - NO cost of resistance (CE-1)
"""

import numpy as np
import pytest

from sswd_evoepi.types import (
    AGENT_DTYPE,
    IDX_EF1A,
    N_ADDITIVE,
    N_LOCI,
    DiseaseState,
    Stage,
    allocate_agents,
    allocate_genotypes,
)
from sswd_evoepi.genetics import (
    W_ADD,
    W_OD,
    S_HET,
    MU_PER_LOCUS,
    N_BANK,
    initialize_effect_sizes,
    compute_resistance_single,
    compute_resistance_batch,
    update_resistance_scores,
    initialize_genotypes,
    initialize_genotypes_beta,
    apply_mutations,
    eliminate_ef1a_lethals,
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
def effect_sizes():
    """Canonical effect sizes for testing (deterministic seed)."""
    return initialize_effect_sizes(np.random.default_rng(12345))


@pytest.fixture
def population_500(effect_sizes, rng):
    """A 500-agent population at HWE for genetics tests."""
    n = 500
    agents = allocate_agents(n)
    genotypes = initialize_genotypes(n, effect_sizes, rng, target_mean_r=0.08)

    for i in range(n):
        agents[i]['alive'] = True
        agents[i]['stage'] = Stage.ADULT
        agents[i]['sex'] = 0 if i < 250 else 1
        agents[i]['size'] = 600.0 + rng.normal(0, 50)
        agents[i]['age'] = 10.0
        agents[i]['disease_state'] = DiseaseState.S
        agents[i]['fecundity_mod'] = 1.0
        agents[i]['node_id'] = 0

    # Compute initial resistance
    update_resistance_scores(agents, genotypes, effect_sizes)

    return agents, genotypes


# ═══════════════════════════════════════════════════════════════════════
# EFFECT SIZE TESTS
# ═══════════════════════════════════════════════════════════════════════


class TestEffectSizes:
    def test_shape(self, effect_sizes):
        assert effect_sizes.shape == (N_ADDITIVE,)
        assert effect_sizes.dtype == np.float64

    def test_sum_equals_w_add(self, effect_sizes):
        assert np.isclose(effect_sizes.sum(), W_ADD, atol=1e-10)

    def test_sorted_descending(self, effect_sizes):
        for i in range(len(effect_sizes) - 1):
            assert effect_sizes[i] >= effect_sizes[i + 1]

    def test_all_positive(self, effect_sizes):
        assert np.all(effect_sizes > 0)

    def test_top_loci_disproportionate(self, effect_sizes):
        """Top 3 loci should contribute >10% of total (exponential tail)."""
        top3_sum = effect_sizes[:3].sum()
        assert top3_sum / W_ADD > 0.10

    def test_reproducible(self):
        """Same seed → same effect sizes."""
        e1 = initialize_effect_sizes(np.random.default_rng(999))
        e2 = initialize_effect_sizes(np.random.default_rng(999))
        np.testing.assert_array_equal(e1, e2)

    def test_different_seeds_differ(self):
        e1 = initialize_effect_sizes(np.random.default_rng(1))
        e2 = initialize_effect_sizes(np.random.default_rng(2))
        assert not np.allclose(e1, e2)


# ═══════════════════════════════════════════════════════════════════════
# RESISTANCE SCORE TESTS
# ═══════════════════════════════════════════════════════════════════════


class TestResistanceScore:
    def test_all_zero_genotype(self, effect_sizes):
        """All-zero genotype → r_i = 0."""
        geno = np.zeros((N_LOCI, 2), dtype=np.int8)
        r = compute_resistance_single(geno, effect_sizes)
        assert r == 0.0

    def test_max_additive_only(self, effect_sizes):
        """All additive loci homozygous resistant, EF1A wild-type → r_i = W_ADD."""
        geno = np.zeros((N_LOCI, 2), dtype=np.int8)
        geno[:N_ADDITIVE, :] = 1
        # EF1A wild-type (0,0)
        r = compute_resistance_single(geno, effect_sizes)
        assert np.isclose(r, W_ADD, atol=1e-6)

    def test_max_with_ef1a_het(self, effect_sizes):
        """All additive resistant + EF1A heterozygote → r_i = 1.0."""
        geno = np.zeros((N_LOCI, 2), dtype=np.int8)
        geno[:N_ADDITIVE, :] = 1
        geno[IDX_EF1A, 0] = 0
        geno[IDX_EF1A, 1] = 1
        r = compute_resistance_single(geno, effect_sizes)
        assert np.isclose(r, 1.0, atol=1e-6)

    def test_ef1a_heterozygote_bonus(self, effect_sizes):
        """EF1A heterozygote adds W_OD to resistance."""
        geno_wt = np.zeros((N_LOCI, 2), dtype=np.int8)
        geno_het = np.zeros((N_LOCI, 2), dtype=np.int8)
        geno_het[IDX_EF1A, 0] = 1

        r_wt = compute_resistance_single(geno_wt, effect_sizes)
        r_het = compute_resistance_single(geno_het, effect_sizes)
        assert np.isclose(r_het - r_wt, W_OD, atol=1e-10)

    def test_ef1a_homozygous_ins_no_bonus(self, effect_sizes):
        """EF1A homozygous ins/ins (lethal) gets no overdominant bonus."""
        geno = np.zeros((N_LOCI, 2), dtype=np.int8)
        geno[IDX_EF1A, :] = 1
        r = compute_resistance_single(geno, effect_sizes)
        # No OD bonus (ef1a_sum == 2, not 1)
        assert np.isclose(r, 0.0, atol=1e-10)

    def test_resistance_in_bounds(self, effect_sizes, population_500):
        """All resistance scores in [0, 1]."""
        agents, genotypes = population_500
        r_vals = agents['resistance'][agents['alive']]
        assert np.all(r_vals >= 0.0)
        assert np.all(r_vals <= 1.0)

    def test_vectorized_matches_loop(self, effect_sizes, population_500):
        """Batch computation matches individual computation."""
        agents, genotypes = population_500
        alive_mask = agents['alive']

        batch_r = compute_resistance_batch(genotypes, effect_sizes, alive_mask)
        alive_idx = np.where(alive_mask)[0]

        for idx in alive_idx[:50]:  # test first 50
            single_r = compute_resistance_single(genotypes[idx], effect_sizes)
            assert np.isclose(
                batch_r[idx], single_r, atol=1e-5
            ), f"Mismatch at idx {idx}: batch={batch_r[idx]}, single={single_r}"

    def test_dead_agents_get_zero(self, effect_sizes):
        """Dead agents should have r_i = 0 in batch computation."""
        n = 10
        genotypes = np.ones((n, N_LOCI, 2), dtype=np.int8)
        alive_mask = np.array([True] * 5 + [False] * 5)
        r = compute_resistance_batch(genotypes, effect_sizes, alive_mask)
        assert np.all(r[5:] == 0.0)
        assert np.all(r[:5] > 0.0)


# ═══════════════════════════════════════════════════════════════════════
# BIMODAL RESISTANCE DISTRIBUTION (ERRATA E6)
# ═══════════════════════════════════════════════════════════════════════


class TestBimodalDistribution:
    def test_bimodal_r_distribution(self, effect_sizes):
        """Pre-SSWD r_i distribution should be bimodal due to EF1A.

        Mode 1: ~0.015 for EF1A +/+ (no overdominant bonus)
        Mode 2: ~0.175 for EF1A +/ins (with W_OD bonus)
        """
        rng = np.random.default_rng(42)
        n = 5000
        geno = initialize_genotypes(n, effect_sizes, rng, target_mean_r=0.08)

        r_vals = np.array([
            compute_resistance_single(geno[i], effect_sizes) for i in range(n)
        ])

        # Identify EF1A heterozygotes
        ef1a_het = geno[:, IDX_EF1A, :].sum(axis=1) == 1
        ef1a_wt = geno[:, IDX_EF1A, :].sum(axis=1) == 0

        r_wt = r_vals[ef1a_wt]
        r_het = r_vals[ef1a_het]

        # Both groups should exist
        assert len(r_wt) > 100, "Too few +/+ individuals"
        assert len(r_het) > 100, "Too few +/ins individuals"

        # Heterozygotes should have higher mean r
        assert np.mean(r_het) > np.mean(r_wt) + W_OD * 0.5, \
            "EF1A heterozygotes should have notably higher resistance"

        # The separation should be roughly W_OD
        gap = np.mean(r_het) - np.mean(r_wt)
        assert abs(gap - W_OD) < 0.05, \
            f"Gap between modes ({gap:.3f}) should be close to W_OD ({W_OD:.3f})"

    def test_population_mean_near_target(self, effect_sizes):
        """Population mean r_i should be near target (0.08)."""
        rng = np.random.default_rng(123)
        n = 5000
        geno = initialize_genotypes(n, effect_sizes, rng, target_mean_r=0.08)
        r_vals = np.array([
            compute_resistance_single(geno[i], effect_sizes) for i in range(n)
        ])
        assert abs(np.mean(r_vals) - 0.08) < 0.03, \
            f"Mean r_i = {np.mean(r_vals):.4f}, expected ~0.08"


# ═══════════════════════════════════════════════════════════════════════
# GENOTYPE INITIALIZATION
# ═══════════════════════════════════════════════════════════════════════


class TestGenotypeInitialization:
    def test_shape(self, effect_sizes, rng):
        n = 100
        geno = initialize_genotypes(n, effect_sizes, rng)
        assert geno.shape == (n, N_LOCI, 2)
        assert geno.dtype == np.int8

    def test_values_binary(self, effect_sizes, rng):
        geno = initialize_genotypes(100, effect_sizes, rng)
        assert np.all((geno == 0) | (geno == 1))

    def test_no_ef1a_lethals(self, effect_sizes, rng):
        """Initial population should contain no ins/ins at EF1A."""
        geno = initialize_genotypes(1000, effect_sizes, rng, ef1a_q=0.24)
        ef1a_sum = geno[:, IDX_EF1A, :].sum(axis=1)
        assert np.all(ef1a_sum < 2), "Found EF1A ins/ins in initial population"

    def test_ef1a_allele_frequency(self, effect_sizes, rng):
        """EF1A allele frequency should be near target."""
        n = 5000
        geno = initialize_genotypes(n, effect_sizes, rng, ef1a_q=0.24)
        q_ef1a = geno[:, IDX_EF1A, :].sum() / (2.0 * n)
        # After lethal elimination, q will be slightly reduced
        assert 0.15 < q_ef1a < 0.30, f"EF1A q = {q_ef1a:.3f}, expected ~0.20-0.24"

    def test_beta_initialization(self, effect_sizes, rng):
        """Beta initialization should produce valid genotypes."""
        geno = initialize_genotypes_beta(500, effect_sizes, rng)
        assert geno.shape == (500, N_LOCI, 2)
        ef1a_sum = geno[:, IDX_EF1A, :].sum(axis=1)
        assert np.all(ef1a_sum < 2)


# ═══════════════════════════════════════════════════════════════════════
# EF1A LETHALITY
# ═══════════════════════════════════════════════════════════════════════


class TestEF1ALethality:
    def test_eliminate_lethals(self):
        """ins/ins individuals at EF1A should be eliminated."""
        n = 100
        geno = np.zeros((n, N_LOCI, 2), dtype=np.int8)
        # Make 20 of them ins/ins at EF1A
        geno[:20, IDX_EF1A, :] = 1
        # Make 30 heterozygous
        geno[20:50, IDX_EF1A, 0] = 1
        geno[20:50, IDX_EF1A, 1] = 0

        viable = eliminate_ef1a_lethals(geno)
        assert viable.sum() == 80  # 100 - 20 lethals
        assert np.all(viable[:20] == False)
        assert np.all(viable[20:] == True)

    def test_lethal_fraction_near_q_squared(self, rng):
        """Fraction of lethals should be ≈ q² for random mating."""
        n_parents = 1000
        q = 0.24
        geno = np.zeros((n_parents, N_LOCI, 2), dtype=np.int8)
        geno[:, IDX_EF1A, 0] = (rng.random(n_parents) < q).astype(np.int8)
        geno[:, IDX_EF1A, 1] = (rng.random(n_parents) < q).astype(np.int8)

        viable = eliminate_ef1a_lethals(geno)
        lethal_fraction = 1.0 - viable.mean()
        expected = q * q
        assert abs(lethal_fraction - expected) < 0.05, \
            f"Lethal fraction {lethal_fraction:.3f}, expected ~{expected:.3f}"


# ═══════════════════════════════════════════════════════════════════════
# MUTATION
# ═══════════════════════════════════════════════════════════════════════


class TestMutation:
    def test_mutation_at_default_rate_negligible(self, rng):
        """At μ=10⁻⁸ with 500 offspring, expect ~0 mutations."""
        n = 500
        geno = np.zeros((n, N_LOCI, 2), dtype=np.int8)
        original = geno.copy()
        n_mut = apply_mutations(geno, rng, mu=MU_PER_LOCUS)
        # Expected: 500 × 52 × 2 × 10⁻⁸ ≈ 5.2×10⁻⁴ mutations → almost always 0
        assert n_mut == 0 or n_mut <= 2

    def test_mutation_at_high_rate(self, rng):
        """At elevated μ, mutations should be detectable."""
        n = 1000
        geno = np.zeros((n, N_LOCI, 2), dtype=np.int8)
        n_mut = apply_mutations(geno, rng, mu=0.01)
        # Expected: 1000 × 52 × 2 × 0.01 = 1040 mutations
        assert n_mut > 500, f"Expected many mutations at μ=0.01, got {n_mut}"
        # Check that some alleles flipped to 1
        assert geno.sum() > 0

    def test_mutation_bidirectional(self, rng):
        """Mutations should flip 0→1 and 1→0."""
        n = 1000
        geno = np.ones((n, N_LOCI, 2), dtype=np.int8)
        original_sum = geno.sum()
        n_mut = apply_mutations(geno, rng, mu=0.01)
        # Some alleles should have flipped 1→0
        assert geno.sum() < original_sum


# ═══════════════════════════════════════════════════════════════════════
# ALLELE FREQUENCIES & HETEROZYGOSITY
# ═══════════════════════════════════════════════════════════════════════


class TestAlleleFrequencies:
    def test_allele_freq_range(self, effect_sizes, population_500):
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

    def test_heterozygosity_agreement(self, effect_sizes, population_500):
        """H_o and H_e should be similar for HWE population."""
        agents, genotypes = population_500
        H_o, H_e = compute_heterozygosity(genotypes, agents['alive'])
        # For a population initialized at HWE, these should agree roughly
        assert abs(H_o - H_e) < 0.05, f"H_o={H_o:.4f}, H_e={H_e:.4f}"

    def test_heterozygosity_range(self, effect_sizes, population_500):
        agents, genotypes = population_500
        H_o, H_e = compute_heterozygosity(genotypes, agents['alive'])
        assert 0.0 <= H_o <= 0.5
        assert 0.0 <= H_e <= 0.5


class TestAdditiveVariance:
    def test_va_positive(self, effect_sizes, population_500):
        agents, genotypes = population_500
        freq = compute_allele_frequencies(genotypes, agents['alive'])
        va = compute_additive_variance(freq, effect_sizes)
        assert va > 0, "V_A should be positive with polymorphic loci"

    def test_va_zero_at_fixation(self, effect_sizes):
        """V_A = 0 when all loci are fixed (q = 0 or q = 1)."""
        freq = np.zeros(N_LOCI, dtype=np.float64)
        va = compute_additive_variance(freq, effect_sizes)
        assert va == 0.0

        freq_fixed = np.ones(N_LOCI, dtype=np.float64)
        va = compute_additive_variance(freq_fixed, effect_sizes)
        assert va == 0.0


# ═══════════════════════════════════════════════════════════════════════
# HARDY-WEINBERG EQUILIBRIUM MAINTENANCE
# ═══════════════════════════════════════════════════════════════════════


class TestHardyWeinberg:
    def test_hwe_initial(self, effect_sizes):
        """Freshly initialized population should be at HWE."""
        rng = np.random.default_rng(42)
        n = 2000
        geno = initialize_genotypes(n, effect_sizes, rng)
        alive = np.ones(n, dtype=bool)

        obs_freq, exp_freq, chi2 = hardy_weinberg_test(geno, alive)

        # Average chi2 across additive loci should be small
        additive_chi2 = chi2[:N_ADDITIVE]
        mean_chi2 = float(np.mean(additive_chi2))
        # With 1 df, chi2 < 3.84 is non-significant at 0.05
        assert mean_chi2 < 5.0, \
            f"Mean chi2 = {mean_chi2:.2f}, population not at HWE"

    def test_hwe_maintained_without_selection(self, effect_sizes):
        """HWE maintained over 50 generations of random mating (no selection).

        This is the critical acceptance test: random Mendelian inheritance
        preserves HWE genotype frequencies.
        """
        rng = np.random.default_rng(42)
        n = 500

        # Initialize at HWE
        geno = initialize_genotypes(n, effect_sizes, rng, target_mean_r=0.08)

        for gen in range(50):
            # Random mating: pick n parent pairs uniformly (no selection)
            mothers = rng.integers(0, n, size=n)
            fathers = rng.integers(0, n, size=n)

            # Pad genotypes into (max_agents, N_LOCI, 2) shape for API
            max_n = n
            geno_padded = geno  # already (n, N_LOCI, 2)

            offspring = mendelian_inherit_batch(
                geno_padded, mothers, fathers, n, rng
            )

            # Eliminate EF1A lethals and resample
            viable = eliminate_ef1a_lethals(offspring)
            offspring = offspring[viable]

            # Maintain population size
            if len(offspring) >= n:
                geno = offspring[:n]
            else:
                # Pad back to n by resampling
                extra_idx = rng.integers(0, len(offspring), size=n - len(offspring))
                geno = np.concatenate([offspring, offspring[extra_idx]])

        # Check HWE at additive loci after 50 generations
        alive = np.ones(n, dtype=bool)
        obs_freq, exp_freq, chi2 = hardy_weinberg_test(geno, alive)
        additive_chi2 = chi2[:N_ADDITIVE]

        # At most 10% of loci should show significant departure
        n_significant = np.sum(additive_chi2 > 3.84)
        fraction_sig = n_significant / N_ADDITIVE
        assert fraction_sig < 0.15, \
            f"{n_significant}/{N_ADDITIVE} loci depart from HWE (expect <15%)"

    def test_allele_freq_stable_without_selection(self, effect_sizes):
        """Allele frequencies should be stable (drift only) over 50 generations.

        Without selection, mean allele frequency across loci should not shift
        directionally. Individual loci will drift, but the ensemble mean
        should remain near the starting value.
        """
        rng = np.random.default_rng(99)
        n = 500

        geno = initialize_genotypes(n, effect_sizes, rng, target_mean_r=0.08)
        alive = np.ones(n, dtype=bool)
        initial_freq = compute_allele_frequencies(geno, alive)
        initial_mean_q = float(np.mean(initial_freq[:N_ADDITIVE]))

        for gen in range(50):
            mothers = rng.integers(0, n, size=n)
            fathers = rng.integers(0, n, size=n)
            offspring = mendelian_inherit_batch(geno, mothers, fathers, n, rng)
            viable = eliminate_ef1a_lethals(offspring)
            offspring = offspring[viable]
            if len(offspring) >= n:
                geno = offspring[:n]
            else:
                extra = rng.integers(0, len(offspring), size=n - len(offspring))
                geno = np.concatenate([offspring, offspring[extra]])

        final_freq = compute_allele_frequencies(geno, alive)
        final_mean_q = float(np.mean(final_freq[:N_ADDITIVE]))

        # Mean allele frequency should not have shifted much
        # (drift, not directional change)
        assert abs(final_mean_q - initial_mean_q) < 0.03, \
            f"Mean q shifted from {initial_mean_q:.4f} to {final_mean_q:.4f}"


# ═══════════════════════════════════════════════════════════════════════
# Ne/N UNDER SRS
# ═══════════════════════════════════════════════════════════════════════


class TestNeOverN:
    def test_ne_ratio_wright_fisher(self):
        """Under Wright-Fisher (uniform offspring), Ne/N ≈ 1."""
        # Equal offspring per parent
        offspring = np.full(100, 10, dtype=np.int64)
        ratio = compute_ne_ratio_from_offspring(offspring)
        # Vk = 0 → Ne = (4N-2)/2 → Ne/N ≈ 2 (exact: (4*100-2)/(0+2)/100 = 1.99)
        assert ratio > 1.5, f"Ne/N = {ratio:.3f}, expected ~2 for uniform"

    def test_ne_ratio_extreme_skew(self):
        """Under extreme SRS, Ne/N << 1."""
        offspring = np.zeros(1000, dtype=np.int64)
        offspring[0] = 5000  # one parent has 5000 offspring
        offspring[1] = 3000
        ratio = compute_ne_ratio_from_offspring(offspring)
        assert ratio < 0.01, f"Ne/N = {ratio:.3f}, expected << 0.01"

    def test_ne_ratio_under_srs(self, effect_sizes):
        """Full SRS lottery should produce Ne/N << 1 (target ~10⁻³).

        This test runs the actual SRS reproductive lottery and checks
        that the offspring distribution produces low Ne/N.
        """
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
            agents[i]['fecundity_mod'] = 1.0

        geno_init = initialize_genotypes(n_parents, effect_sizes, rng)
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

        # Count offspring per parent
        all_parents = parent_pairs.ravel()
        _, counts = np.unique(all_parents, return_counts=True)

        # Pad with zeros for parents that got no offspring
        full_counts = np.zeros(n_parents, dtype=np.int64)
        unique_parents, parent_counts = np.unique(all_parents, return_counts=True)
        full_counts[unique_parents] = parent_counts

        ratio = compute_ne_ratio_from_offspring(full_counts)
        # Should be very low — not exactly 10⁻³ but much less than 1
        assert ratio < 0.1, f"Ne/N = {ratio:.4f}, expected << 1 under SRS"


# ═══════════════════════════════════════════════════════════════════════
# EF1A EQUILIBRIUM
# ═══════════════════════════════════════════════════════════════════════


class TestEF1AEquilibrium:
    def test_ef1a_maintained_over_generations(self, effect_sizes):
        """EF1A allele frequency should decrease predictably under lethal purging.

        Without disease (no heterozygote advantage in survival), the lethal
        ins/ins homozygote purges the insertion allele. Each generation:
          q' ≈ q / (1 + q)
        Starting at q=0.24, after 10 generations: q ≈ 1/(1/0.24 + 10) ≈ 0.07.

        The key test: q should still be > 0 (allele not lost) and should
        decrease (consistent with lethal purging, not drift alone).
        """
        rng = np.random.default_rng(42)
        n = 1000
        q_target = 0.24
        geno = initialize_genotypes(n, effect_sizes, rng, ef1a_q=q_target)

        alive = np.ones(n, dtype=bool)
        q_initial = compute_allele_frequencies(geno, alive)[IDX_EF1A]

        for gen in range(10):
            mothers = rng.integers(0, n, size=int(n * 1.1))
            fathers = rng.integers(0, n, size=int(n * 1.1))
            offspring = mendelian_inherit_batch(geno, mothers, fathers, len(mothers), rng)
            viable = eliminate_ef1a_lethals(offspring)
            offspring = offspring[viable]
            if len(offspring) >= n:
                geno = offspring[:n]
            else:
                extra = rng.integers(0, len(offspring), size=n - len(offspring))
                geno = np.concatenate([offspring, offspring[extra]])

        q_final = compute_allele_frequencies(geno, alive)[IDX_EF1A]
        # q should have decreased (lethal purging)
        assert q_final < q_initial, \
            f"EF1A q should decrease: initial={q_initial:.3f}, final={q_final:.3f}"
        # q should not be lost entirely in just 10 generations at n=1000
        assert q_final > 0.01, \
            f"EF1A q = {q_final:.3f}, should not be lost in only 10 generations"
        # Theoretical: q ≈ 1/(1/q0 + t) ≈ 0.07 after 10 generations
        q_expected = 1.0 / (1.0 / q_initial + 10)
        assert abs(q_final - q_expected) < 0.05, \
            f"EF1A q = {q_final:.3f}, expected ~{q_expected:.3f} (lethal purging model)"


# ═══════════════════════════════════════════════════════════════════════
# GENETIC DIAGNOSTICS
# ═══════════════════════════════════════════════════════════════════════


class TestGeneticDiagnostics:
    def test_diagnostics_populated(self, effect_sizes, population_500):
        agents, genotypes = population_500
        diag = compute_genetic_diagnostics(agents, genotypes, effect_sizes)

        assert diag.n_alive == 500
        assert diag.mean_resistance > 0
        assert diag.var_resistance > 0
        assert diag.va_additive > 0
        assert diag.ef1a_allele_freq > 0
        assert diag.ef1a_het_freq > 0
        assert diag.heterozygosity_obs > 0
        assert diag.heterozygosity_exp > 0

    def test_diagnostics_allele_freq_change(self, effect_sizes, population_500):
        agents, genotypes = population_500
        prev_freq = np.full(N_LOCI, 0.1, dtype=np.float64)
        diag = compute_genetic_diagnostics(
            agents, genotypes, effect_sizes, prev_allele_freq=prev_freq
        )
        assert diag.delta_q_top3.shape == (3,)
        assert diag.mean_abs_delta_q > 0

    def test_diagnostics_empty_population(self, effect_sizes):
        agents = allocate_agents(10)
        genotypes = allocate_genotypes(10)
        diag = compute_genetic_diagnostics(agents, genotypes, effect_sizes)
        assert diag.n_alive == 0
        assert diag.mean_resistance == 0.0


# ═══════════════════════════════════════════════════════════════════════
# GENOTYPE BANK (TIER 2)
# ═══════════════════════════════════════════════════════════════════════


class TestGenotypeBank:
    def test_compress_roundtrip(self, effect_sizes, population_500):
        agents, genotypes = population_500
        rng = np.random.default_rng(42)

        bank = compress_to_genotype_bank(
            genotypes, agents['alive'], effect_sizes, rng
        )
        assert bank.bank.shape == (N_BANK, N_LOCI, 2)
        assert bank.bank_resistance.shape == (N_BANK,)
        assert abs(bank.mean_resistance - agents['resistance'][agents['alive']].mean()) < 0.05

    def test_expand_produces_valid_genotypes(self, effect_sizes, population_500):
        agents, genotypes = population_500
        rng = np.random.default_rng(42)

        bank = compress_to_genotype_bank(
            genotypes, agents['alive'], effect_sizes, rng
        )
        expanded = expand_genotype_bank(bank, 300, effect_sizes, rng)
        assert expanded.shape == (300, N_LOCI, 2)
        assert np.all((expanded == 0) | (expanded == 1))

    def test_expand_preserves_resistance_range(self, effect_sizes, population_500):
        agents, genotypes = population_500
        rng = np.random.default_rng(42)

        bank = compress_to_genotype_bank(
            genotypes, agents['alive'], effect_sizes, rng
        )
        expanded = expand_genotype_bank(bank, 500, effect_sizes, rng)

        r_expanded = np.array([
            compute_resistance_single(expanded[i], effect_sizes) for i in range(500)
        ])
        assert np.all(r_expanded >= 0)
        assert np.all(r_expanded <= 1.0)
        # Mean should be in same ballpark
        assert abs(np.mean(r_expanded) - bank.mean_resistance) < 0.05

    def test_bank_update_summary(self, effect_sizes, population_500):
        agents, genotypes = population_500
        rng = np.random.default_rng(42)

        bank = compress_to_genotype_bank(
            genotypes, agents['alive'], effect_sizes, rng
        )
        assert bank.heterozygosity > 0
        assert len(bank.allele_freq) == N_LOCI


# ═══════════════════════════════════════════════════════════════════════
# F_ST
# ═══════════════════════════════════════════════════════════════════════


class TestFST:
    def test_fst_identical_populations(self, effect_sizes, population_500):
        """F_ST = 0 for identical populations."""
        agents, genotypes = population_500
        freq = compute_allele_frequencies(genotypes, agents['alive'])
        fst = compute_fst([freq, freq, freq])
        assert np.isclose(fst, 0.0, atol=1e-10)

    def test_fst_different_populations(self, effect_sizes):
        """F_ST > 0 for genetically different populations."""
        rng = np.random.default_rng(42)
        n = 200

        geno1 = initialize_genotypes(n, effect_sizes, rng, target_mean_r=0.05)
        geno2 = initialize_genotypes(n, effect_sizes, rng, target_mean_r=0.15)

        alive = np.ones(n, dtype=bool)
        freq1 = compute_allele_frequencies(geno1, alive)
        freq2 = compute_allele_frequencies(geno2, alive)

        fst = compute_fst([freq1, freq2])
        assert fst > 0, "F_ST should be > 0 for different populations"

    def test_fst_single_population(self, effect_sizes, population_500):
        """F_ST = 0 for a single population (< 2 nodes)."""
        agents, genotypes = population_500
        freq = compute_allele_frequencies(genotypes, agents['alive'])
        fst = compute_fst([freq])
        assert fst == 0.0


# ═══════════════════════════════════════════════════════════════════════
# NO COST OF RESISTANCE (CE-1)
# ═══════════════════════════════════════════════════════════════════════


class TestNoCostOfResistance:
    def test_fecundity_mod_always_one(self, effect_sizes, population_500):
        """CE-1: fecundity_mod should always be 1.0 after update."""
        agents, genotypes = population_500
        update_resistance_scores(agents, genotypes, effect_sizes)
        alive = agents['alive']
        assert np.all(agents['fecundity_mod'][alive] == 1.0), \
            "fecundity_mod should be 1.0 — no cost of resistance (CE-1)"


# ═══════════════════════════════════════════════════════════════════════
# UPDATE RESISTANCE SCORES (in-place)
# ═══════════════════════════════════════════════════════════════════════


class TestUpdateResistanceScores:
    def test_writes_to_agents(self, effect_sizes):
        rng = np.random.default_rng(42)
        n = 50
        agents = allocate_agents(n)
        genotypes = allocate_genotypes(n)
        geno = initialize_genotypes(n, effect_sizes, rng)
        genotypes[:n] = geno

        for i in range(n):
            agents[i]['alive'] = True
            agents[i]['resistance'] = -1.0  # sentinel

        update_resistance_scores(agents, genotypes, effect_sizes)
        assert np.all(agents['resistance'][:n] >= 0.0)
        assert np.all(agents['resistance'][:n] <= 1.0)

    def test_dead_unchanged(self, effect_sizes):
        rng = np.random.default_rng(42)
        n = 10
        agents = allocate_agents(n)
        genotypes = allocate_genotypes(n)
        # Don't set alive → all dead
        agents['resistance'] = -999.0

        update_resistance_scores(agents, genotypes, effect_sizes)
        # Dead agents should get 0.0 from batch (overwriting sentinel)
        # because compute_resistance_batch returns 0 for dead
        assert np.all(agents['resistance'] == 0.0)


# ═══════════════════════════════════════════════════════════════════════
# INTEGRATION: MENDELIAN INHERITANCE + GENETICS
# ═══════════════════════════════════════════════════════════════════════


class TestMendelianIntegration:
    def test_offspring_alleles_from_parents(self, effect_sizes):
        """Each offspring allele must come from one of its parents."""
        rng = np.random.default_rng(42)
        n = 50
        geno = initialize_genotypes(n, effect_sizes, rng)

        mothers = rng.integers(0, n, size=200)
        fathers = rng.integers(0, n, size=200)
        offspring = mendelian_inherit_batch(geno, mothers, fathers, 200, rng)

        for i in range(min(50, len(offspring))):
            m = mothers[i]
            f = fathers[i]
            for l in range(N_LOCI):
                # Maternal allele (copy 0) must be one of mother's alleles
                assert offspring[i, l, 0] in (geno[m, l, 0], geno[m, l, 1]), \
                    f"Offspring {i} locus {l} maternal allele not from mother"
                # Paternal allele (copy 1) must be one of father's alleles
                assert offspring[i, l, 1] in (geno[f, l, 0], geno[f, l, 1]), \
                    f"Offspring {i} locus {l} paternal allele not from father"


# ═══════════════════════════════════════════════════════════════════════
# CONSTANTS VERIFICATION
# ═══════════════════════════════════════════════════════════════════════


class TestConstants:
    def test_w_add_plus_w_od_equals_one(self):
        assert np.isclose(W_ADD + W_OD, 1.0, atol=1e-10)

    def test_w_od_formula(self):
        expected = S_HET / (1.0 + S_HET)
        assert np.isclose(W_OD, expected, atol=1e-10)

    def test_loci_counts(self):
        assert N_ADDITIVE == 51
        assert IDX_EF1A == 51
        assert N_LOCI == 52
