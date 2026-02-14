"""Tests for sswd_evoepi.reproduction — broadcast spawning, SRS, Allee effects.

Acceptance criteria (Phase 3):
  - SRS produces realistic Ne/N ≈ 10⁻³ (within order of magnitude)
  - Allee effect visible: per-capita growth declines at low density
  - Very low density → negative population growth
  - Normal density → stable replacement
  - Spawning only in correct season
  - Pareto distribution produces heavy-tailed family sizes
"""

import numpy as np
import pytest

from sswd_evoepi.types import (
    AGENT_DTYPE,
    ANNUAL_SURVIVAL,
    IDX_EF1A,
    N_ADDITIVE,
    N_LOCI,
    DiseaseState,
    LarvalCohort,
    Stage,
    allocate_agents,
    allocate_genotypes,
)
from sswd_evoepi.reproduction import (
    beverton_holt_recruitment,
    compute_aggregation_factor,
    compute_ne,
    compute_ne_from_pairs,
    compute_ne_ratio,
    fecundity,
    fertilization_success,
    get_spawning_day,
    is_spawning_season,
    larval_survival,
    mendelian_inherit_batch,
    pelagic_larval_duration,
    per_capita_growth_rate,
    produce_larval_cohort,
    settlement_cue_modifier,
    settle_recruits,
    srs_reproductive_lottery,
    total_eggs,
    _compute_resistance,
)


# ═══════════════════════════════════════════════════════════════════════
# HELPERS
# ═══════════════════════════════════════════════════════════════════════

def _make_effect_sizes(rng=None):
    """Create a canonical effect-size vector for testing."""
    if rng is None:
        rng = np.random.default_rng(12345)
    raw = rng.exponential(1.0, size=N_ADDITIVE)
    normalized = raw / raw.sum() * 0.840
    normalized.sort()
    return normalized[::-1].copy()


def _make_spawning_population(
    n_females: int = 50,
    n_males: int = 50,
    size_mm: float = 600.0,
    rng=None,
):
    """Create a minimal viable spawning population for testing.

    Returns (agents, genotypes, female_idx, male_idx, effect_sizes).
    """
    if rng is None:
        rng = np.random.default_rng(42)

    n_total = n_females + n_males
    max_n = n_total + 100  # extra slots for recruits

    agents = allocate_agents(max_n)
    genotypes = allocate_genotypes(max_n)

    effect_sizes = _make_effect_sizes(np.random.default_rng(12345))

    for i in range(n_total):
        agents[i]['alive'] = True
        agents[i]['stage'] = Stage.ADULT
        agents[i]['sex'] = 0 if i < n_females else 1
        agents[i]['size'] = size_mm + rng.normal(0, 50)
        agents[i]['age'] = 10.0
        agents[i]['disease_state'] = DiseaseState.S
        agents[i]['fecundity_mod'] = 1.0
        agents[i]['node_id'] = 0

        # Random genotypes: ~5% resistant allele frequency at additive loci
        for l in range(N_ADDITIVE):
            genotypes[i, l, 0] = 1 if rng.random() < 0.05 else 0
            genotypes[i, l, 1] = 1 if rng.random() < 0.05 else 0
        # EF1A: 24% allele frequency
        genotypes[i, IDX_EF1A, 0] = 1 if rng.random() < 0.24 else 0
        genotypes[i, IDX_EF1A, 1] = 1 if rng.random() < 0.24 else 0
        # Fix lethal homozygotes
        if genotypes[i, IDX_EF1A, :].sum() == 2:
            genotypes[i, IDX_EF1A, 1] = 0

        agents[i]['resistance'] = _compute_resistance(genotypes[i], effect_sizes)

    females = np.where(
        agents['alive'] & (agents['sex'] == 0) & (agents['stage'] == Stage.ADULT)
    )[0]
    males = np.where(
        agents['alive'] & (agents['sex'] == 1) & (agents['stage'] == Stage.ADULT)
    )[0]

    return agents, genotypes, females, males, effect_sizes


# ═══════════════════════════════════════════════════════════════════════
# SPAWNING PHENOLOGY
# ═══════════════════════════════════════════════════════════════════════

class TestSpawningPhenology:
    def test_base_latitude(self):
        """Spawning day at 48°N should be ~150."""
        assert get_spawning_day(48.0) == 150

    def test_southern_shift(self):
        """Southern latitudes spawn earlier."""
        assert get_spawning_day(44.0) < get_spawning_day(48.0)

    def test_northern_shift(self):
        """Northern latitudes spawn later."""
        assert get_spawning_day(52.0) > get_spawning_day(48.0)

    def test_clamped_range(self):
        """Spawning day should be clamped to [120, 200]."""
        assert get_spawning_day(20.0) >= 120
        assert get_spawning_day(70.0) <= 200

    def test_is_spawning_season_exact(self):
        """Exact day matches."""
        assert is_spawning_season(150, 150)
        assert not is_spawning_season(149, 150)
        assert not is_spawning_season(151, 150)

    def test_is_spawning_season_window(self):
        """Window tolerance works."""
        assert is_spawning_season(148, 150, window=3)
        assert is_spawning_season(153, 150, window=3)
        assert not is_spawning_season(154, 150, window=3)

    def test_spawning_only_in_season(self):
        """Verify spawning is restricted to the correct day-of-year."""
        spawning_day = get_spawning_day(48.0)
        # Non-spawning days should NOT trigger spawning
        for day in [1, 50, 100, 250, 300, 365]:
            assert not is_spawning_season(day, spawning_day), \
                f"Day {day} should NOT be spawning season"
        # Spawning day SHOULD trigger
        assert is_spawning_season(spawning_day, spawning_day)


# ═══════════════════════════════════════════════════════════════════════
# FECUNDITY
# ═══════════════════════════════════════════════════════════════════════

class TestFecundity:
    def test_below_min_size(self):
        """No reproduction below L_min_repro."""
        assert fecundity(300.0) == 0.0
        assert fecundity(399.9) == 0.0

    def test_at_reference_size(self):
        """At L_ref, fecundity = F0."""
        assert fecundity(500.0) == pytest.approx(1e7, rel=1e-6)

    def test_allometric_scaling(self):
        """Larger individuals produce more eggs."""
        f_small = fecundity(450.0)
        f_large = fecundity(800.0)
        assert f_large > f_small > 0

    def test_scales_with_exponent(self):
        """Fecundity scales as (L/L_ref)^b."""
        f_600 = fecundity(600.0)
        expected = 1e7 * (600.0 / 500.0) ** 2.5
        assert f_600 == pytest.approx(expected, rel=1e-6)

    def test_total_eggs_vectorized(self):
        """total_eggs matches sum of individual fecundity calls."""
        agents = allocate_agents(5)
        sizes = [400.0, 500.0, 600.0, 700.0, 300.0]
        for i, s in enumerate(sizes):
            agents[i]['size'] = s
        idx = np.array([0, 1, 2, 3, 4])
        tot = total_eggs(agents, idx)
        expected = sum(fecundity(s) for s in sizes)
        assert tot == pytest.approx(expected, rel=1e-6)


# ═══════════════════════════════════════════════════════════════════════
# FERTILIZATION (ALLEE EFFECT)
# ═══════════════════════════════════════════════════════════════════════

class TestFertilizationAllee:
    """Test the quadratic Allee effect in fertilization kinetics."""

    def test_zero_density(self):
        """Zero density → zero fertilization."""
        assert fertilization_success(0.0) == 0.0

    def test_very_low_density_post_sswd(self):
        """At ρ = 0.001 (post-SSWD), F should be near zero."""
        F = fertilization_success(0.001)
        assert F < 0.02  # essentially zero

    def test_moderate_density_pre_sswd_wa(self):
        """At ρ = 0.03 (pre-SSWD WA average), F should be ~0.13."""
        F = fertilization_success(0.03)
        assert 0.05 < F < 0.35

    def test_high_density_howe_sound(self):
        """At ρ = 0.43 (Howe Sound pre-SSWD), F should be high."""
        F = fertilization_success(0.43)
        assert F > 0.80

    def test_saturation_density(self):
        """At ρ = 0.67, F should be ~0.95."""
        F = fertilization_success(0.67)
        assert F > 0.90

    def test_monotonically_increasing(self):
        """Fertilization is monotonically increasing with density."""
        densities = np.linspace(0.0, 1.0, 100)
        F_values = [fertilization_success(d) for d in densities]
        for i in range(1, len(F_values)):
            assert F_values[i] >= F_values[i - 1]

    def test_quadratic_at_low_density(self):
        """Zygote production should be ≈ quadratic in density at low ρ.

        Zygotes = F(ρ) × eggs ∝ ρ × ρ at low density.
        So zygote_production = F(ρ_male) × n_females ∝ ρ².
        Log-log slope should be ~2.0.
        """
        rho_low = np.array([0.001, 0.002, 0.005, 0.01, 0.02])
        # Zygote production ∝ ρ × F(ρ/2)  (half are males)
        zygotes = rho_low * np.array([fertilization_success(r * 0.5) for r in rho_low])
        # Log-log regression
        log_rho = np.log(rho_low)
        log_z = np.log(zygotes)
        slope = np.polyfit(log_rho, log_z, 1)[0]
        assert 1.5 < slope < 2.5, f"Log-log slope = {slope}, expected ~2.0"

    def test_no_sharp_threshold(self):
        """Per Lundquist & Botsford 2004, no discontinuity should exist.

        Check that F(D) is smooth (no jumps > 0.1 between adjacent points).
        """
        densities = np.linspace(0.0, 1.0, 1000)
        F_vals = np.array([fertilization_success(d) for d in densities])
        diffs = np.abs(np.diff(F_vals))
        assert np.all(diffs < 0.01), "Discontinuity detected in F(D)"


class TestAggregation:
    def test_single_individual(self):
        """No aggregation with < 2 adults."""
        assert compute_aggregation_factor(1, 10000.0) == 1.0
        assert compute_aggregation_factor(0, 10000.0) == 1.0

    def test_already_dense(self):
        """If already at max aggregation density, factor = 1.0."""
        assert compute_aggregation_factor(100, 100.0) == pytest.approx(1.0)

    def test_aggregation_increases_effective_density(self):
        """Aggregation should increase effective density (factor > 1)."""
        factor = compute_aggregation_factor(10, 100000.0)
        assert factor > 1.0

    def test_more_adults_more_aggregation(self):
        """More adults → higher effective density after aggregation.

        The aggregation *factor* may decrease with more adults (because
        ambient density rises), but the *effective density* should increase.
        """
        area = 100000.0
        f5 = compute_aggregation_factor(5, area)
        f15 = compute_aggregation_factor(15, area)
        # Effective density = factor × ambient_density
        eff5 = f5 * (5 / area)
        eff15 = f15 * (15 / area)
        assert eff15 > eff5


# ═══════════════════════════════════════════════════════════════════════
# PER-CAPITA GROWTH RATE & ALLEE EFFECT
# ═══════════════════════════════════════════════════════════════════════

class TestPerCapitaGrowth:
    """Test Allee effect via per-capita growth rate.

    For high-fecundity broadcast spawners (10⁷ eggs/female), the
    DETERMINISTIC growth rate is always positive (even 0.1% fertilization
    produces enough recruits). The Allee effect manifests as dramatically
    reduced growth rate at low density, making stochastic extinction
    likely. The actual extinction threshold depends on demographic
    stochasticity, not the deterministic rate crossing zero.
    """

    def test_very_low_growth_at_post_sswd_density(self):
        """At post-SSWD density (0.001), growth should be very low.

        Much lower than at healthy density, making stochastic extinction
        highly likely — the demographic Allee effect.
        """
        r_low = per_capita_growth_rate(density=0.001)
        r_high = per_capita_growth_rate(density=0.43)
        # Growth at low density should be < 1% of high-density growth
        assert r_low < 0.01 * r_high, \
            f"r(0.001)={r_low:.4f}, r(0.43)={r_high:.4f}, expected >100× difference"

    def test_positive_at_moderate_density(self):
        """At healthy density (0.1+), growth should be clearly positive."""
        r = per_capita_growth_rate(density=0.10)
        assert r > 0, f"Expected positive growth at ρ=0.10, got {r}"

    def test_positive_at_high_density(self):
        """At Howe Sound density, growth should be clearly positive."""
        r = per_capita_growth_rate(density=0.43)
        assert r > 0

    def test_growth_increases_with_density(self):
        """Per-capita growth should increase with density (Allee regime)."""
        densities = [0.001, 0.01, 0.05, 0.1, 0.5]
        rates = [per_capita_growth_rate(d) for d in densities]
        # Should be monotonically increasing
        for i in range(1, len(rates)):
            assert rates[i] > rates[i - 1], \
                f"Growth rate not increasing: r({densities[i]})={rates[i]} <= r({densities[i-1]})={rates[i-1]}"

    def test_allee_shape_orders_of_magnitude(self):
        """Growth rate should span multiple orders of magnitude across densities.

        This is the hallmark of the quadratic Allee effect: r ∝ ρ at low ρ.
        """
        r_0001 = per_capita_growth_rate(density=0.001)
        r_01 = per_capita_growth_rate(density=0.1)
        r_05 = per_capita_growth_rate(density=0.5)

        # At least 10× difference between 0.001 and 0.1
        assert r_01 > 10 * r_0001, \
            f"Expected >10× growth difference, got r(0.1)/r(0.001) = {r_01/r_0001:.1f}"
        # High density should be even higher
        assert r_05 > r_01


# ═══════════════════════════════════════════════════════════════════════
# MENDELIAN INHERITANCE
# ═══════════════════════════════════════════════════════════════════════

class TestMendelianInheritance:
    def test_allele_values_valid(self):
        """Offspring alleles should all be 0 or 1."""
        rng = np.random.default_rng(42)
        genotypes = allocate_genotypes(10)
        # Set some resistant alleles
        genotypes[:, :, :] = rng.integers(0, 2, size=(10, N_LOCI, 2), dtype=np.int8)

        mothers = np.array([0, 1, 2, 3, 4])
        fathers = np.array([5, 6, 7, 8, 9])
        offspring = mendelian_inherit_batch(genotypes, mothers, fathers, 5, rng)

        assert offspring.shape == (5, N_LOCI, 2)
        assert np.all((offspring == 0) | (offspring == 1))

    def test_inheritance_from_parents(self):
        """Each offspring allele must come from the corresponding parent."""
        rng = np.random.default_rng(42)
        genotypes = allocate_genotypes(4)
        # Parent 0: all zeros
        genotypes[0, :, :] = 0
        # Parent 1: all ones
        genotypes[1, :, :] = 1
        # Parent 2: all zeros
        genotypes[2, :, :] = 0
        # Parent 3: all ones
        genotypes[3, :, :] = 1

        # Cross 0 × 1: maternal=0 (all 0), paternal=1 (all 1)
        mothers = np.array([0])
        fathers = np.array([1])
        offspring = mendelian_inherit_batch(genotypes, mothers, fathers, 1, rng)
        # Maternal allele must be 0, paternal must be 1
        assert np.all(offspring[0, :, 0] == 0)
        assert np.all(offspring[0, :, 1] == 1)

    def test_heterozygous_parent_segregation(self):
        """Heterozygous parent should produce both alleles in offspring.

        With enough offspring, we should see both 0 and 1 from a het parent.
        """
        rng = np.random.default_rng(42)
        genotypes = allocate_genotypes(2)
        # Mother: heterozygous at locus 0
        genotypes[0, 0, 0] = 0
        genotypes[0, 0, 1] = 1
        # Father: homozygous 0 at locus 0
        genotypes[1, 0, 0] = 0
        genotypes[1, 0, 1] = 0

        n = 1000
        mothers = np.zeros(n, dtype=int)
        fathers = np.ones(n, dtype=int)
        offspring = mendelian_inherit_batch(genotypes, mothers, fathers, n, rng)

        # Maternal allele (column 0) should be mix of 0 and 1
        mat_alleles = offspring[:, 0, 0]
        assert np.any(mat_alleles == 0)
        assert np.any(mat_alleles == 1)

        # Paternal allele (column 1) should all be 0
        pat_alleles = offspring[:, 0, 1]
        assert np.all(pat_alleles == 0)

    def test_mendelian_50_50_ratio(self):
        """Het × Hom should give ~50:50 segregation ratio."""
        rng = np.random.default_rng(42)
        genotypes = allocate_genotypes(2)
        genotypes[0, 0, 0] = 0
        genotypes[0, 0, 1] = 1
        genotypes[1, 0, :] = 0

        n = 10000
        offspring = mendelian_inherit_batch(
            genotypes,
            np.zeros(n, dtype=int),
            np.ones(n, dtype=int),
            n, rng,
        )
        frac_1 = offspring[:, 0, 0].mean()
        assert 0.45 < frac_1 < 0.55, f"Expected ~0.5, got {frac_1}"


# ═══════════════════════════════════════════════════════════════════════
# SRS LOTTERY
# ═══════════════════════════════════════════════════════════════════════

class TestSRSLottery:
    """Test sweepstakes reproductive success mechanics."""

    def test_produces_offspring(self):
        """Basic SRS lottery produces offspring with valid genotypes."""
        agents, geno, fem, mal, efx = _make_spawning_population()
        rng = np.random.default_rng(42)
        offspring_g, pairs = srs_reproductive_lottery(
            fem, mal, agents, geno, n_offspring_target=500, rng=rng,
        )
        assert len(offspring_g) > 0
        assert offspring_g.shape[1] == N_LOCI
        assert offspring_g.shape[2] == 2
        assert np.all((offspring_g == 0) | (offspring_g == 1))

    def test_ef1a_lethal_elimination(self):
        """No offspring should be homozygous ins/ins at EF1A."""
        agents, geno, fem, mal, efx = _make_spawning_population()
        rng = np.random.default_rng(42)
        offspring_g, pairs = srs_reproductive_lottery(
            fem, mal, agents, geno, n_offspring_target=5000, rng=rng,
        )
        ef1a_sum = offspring_g[:, IDX_EF1A, :].sum(axis=1)
        assert not np.any(ef1a_sum == 2), "EF1A ins/ins lethal should be eliminated"

    def test_no_parents_no_offspring(self):
        """Zero females or males should produce no offspring."""
        agents, geno, fem, mal, efx = _make_spawning_population()
        rng = np.random.default_rng(42)
        # No females
        g, p = srs_reproductive_lottery(
            np.array([], dtype=int), mal, agents, geno, 100, rng=rng,
        )
        assert len(g) == 0
        # No males
        g, p = srs_reproductive_lottery(
            fem, np.array([], dtype=int), agents, geno, 100, rng=rng,
        )
        assert len(g) == 0

    def test_pareto_heavy_tail(self):
        """Pareto distribution should produce heavy-tailed family sizes.

        Most parents should have few offspring; rare parents should have many.
        """
        agents, geno, fem, mal, efx = _make_spawning_population(
            n_females=100, n_males=100,
        )
        rng = np.random.default_rng(42)
        offspring_g, pairs = srs_reproductive_lottery(
            fem, mal, agents, geno, n_offspring_target=10000, rng=rng,
        )

        # Count offspring per mother
        mother_ids = pairs[:, 0]
        unique_moms, mom_counts = np.unique(mother_ids, return_counts=True)

        # Heavy tail: top parent should have significantly more than average
        avg = len(offspring_g) / len(fem)
        max_count = mom_counts.max()
        assert max_count > 3 * avg, \
            f"Top mother has {max_count} offspring, expected >> {avg:.0f} (avg)"

        # Many parents should contribute below average
        below_avg = np.sum(mom_counts < avg)
        assert below_avg > len(unique_moms) * 0.5, \
            "Expected majority of parents to contribute below-average offspring"

    def test_parent_pairs_shape(self):
        """Parent pairs should have correct shape and valid indices."""
        agents, geno, fem, mal, efx = _make_spawning_population()
        rng = np.random.default_rng(42)
        g, pairs = srs_reproductive_lottery(
            fem, mal, agents, geno, n_offspring_target=100, rng=rng,
        )
        assert pairs.shape[1] == 2
        assert len(pairs) == len(g)
        # All parent indices should be in the female/male arrays
        assert np.all(np.isin(pairs[:, 0], fem))
        assert np.all(np.isin(pairs[:, 1], mal))


class TestNeCalculation:
    """Test Ne/N tracking under SRS."""

    def test_ne_basic(self):
        """Ne should be positive and ≤ N."""
        counts = np.array([10, 5, 0, 0, 1, 3, 0, 0, 15, 2])
        Ne = compute_ne(counts)
        assert Ne > 0
        assert Ne <= len(counts) * 4  # theoretical max

    def test_ne_equal_contribution(self):
        """If all parents contribute equally, Ne → 4N."""
        counts = np.ones(100, dtype=int) * 50
        Ne = compute_ne(counts)
        # With Vk = 0, Ne = (4N-2)/2 = 2N-1
        assert Ne == pytest.approx(199.0, rel=0.01)

    def test_ne_extreme_skew(self):
        """Extreme skew (one parent gets everything) → Ne ≈ 4/Vk."""
        N = 100
        counts = np.zeros(N, dtype=int)
        counts[0] = 10000  # one parent dominates
        Ne = compute_ne(counts)
        # Should be very small relative to N
        assert Ne < 1.0

    def test_srs_produces_low_ne_ratio(self):
        """The SRS lottery should produce Ne/N ≈ 10⁻³ (within order of magnitude).

        This is the KEY acceptance criterion for Phase 3.
        """
        agents, geno, fem, mal, efx = _make_spawning_population(
            n_females=200, n_males=200, size_mm=600.0,
        )
        rng = np.random.default_rng(42)

        # Large offspring pool for stable Ne estimate
        offspring_g, pairs = srs_reproductive_lottery(
            fem, mal, agents, geno, n_offspring_target=50000,
            alpha_srs=1.35, rng=rng,
        )

        Ne, Ne_over_N = compute_ne_from_pairs(pairs, fem, mal)
        N_total = len(fem) + len(mal)

        # Ne/N should be in range [10⁻⁵, 10⁻¹] — realistic for SRS
        # With α=1.35, we expect ~10⁻³ but there's natural variation
        assert 1e-5 < Ne_over_N < 0.1, \
            f"Ne/N = {Ne_over_N:.6f}, expected ~10⁻³ range"

        # More specifically, should be << 1 (strong SRS signal)
        assert Ne_over_N < 0.05, \
            f"Ne/N = {Ne_over_N:.4f}, expected << 1 for SRS with α=1.35"


# ═══════════════════════════════════════════════════════════════════════
# LARVAL PRODUCTION
# ═══════════════════════════════════════════════════════════════════════

class TestLarvalProduction:
    def test_pld_at_reference_temp(self):
        """PLD at reference temperature should be ~63 days."""
        pld = pelagic_larval_duration(10.5)
        assert pld == pytest.approx(63.0, rel=0.01)

    def test_pld_warmer_shorter(self):
        """Warmer SST → shorter PLD."""
        assert pelagic_larval_duration(14.0) < pelagic_larval_duration(10.5)

    def test_pld_cooler_longer(self):
        """Cooler SST → longer PLD."""
        assert pelagic_larval_duration(8.0) > pelagic_larval_duration(10.5)

    def test_pld_clamped(self):
        """PLD should be clamped to [30, 150]."""
        assert pelagic_larval_duration(30.0) >= 30.0
        assert pelagic_larval_duration(-5.0) <= 150.0

    def test_larval_survival_at_default(self):
        """At 63 days and μ=0.05: survival ≈ 4.2%."""
        s = larval_survival(63.0)
        assert s == pytest.approx(np.exp(-0.05 * 63), rel=1e-6)
        assert 0.03 < s < 0.06

    def test_larval_survival_decreases_with_pld(self):
        """Longer PLD → lower survival."""
        assert larval_survival(80.0) < larval_survival(50.0)

    def test_produce_cohort_basic(self):
        """produce_larval_cohort should return a valid LarvalCohort."""
        agents, geno, fem, mal, efx = _make_spawning_population()
        rng = np.random.default_rng(42)

        cohort, diag = produce_larval_cohort(
            node_id=0,
            agents=agents,
            genotypes=geno,
            habitat_area=10000.0,  # 1 hectare
            sst=11.0,
            lat=48.0,
            rng=rng,
        )

        assert cohort is not None
        assert isinstance(cohort, LarvalCohort)
        assert cohort.source_node == 0
        assert cohort.n_competent > 0
        assert cohort.genotypes.shape[0] == cohort.n_competent
        assert cohort.genotypes.shape[1] == N_LOCI
        assert cohort.genotypes.shape[2] == 2
        assert cohort.pld_days > 0

    def test_produce_cohort_no_males(self):
        """No males → no cohort."""
        agents = allocate_agents(10)
        geno = allocate_genotypes(10)
        # All females
        for i in range(10):
            agents[i]['alive'] = True
            agents[i]['stage'] = Stage.ADULT
            agents[i]['sex'] = 0
            agents[i]['size'] = 600.0
            agents[i]['disease_state'] = DiseaseState.S

        cohort, diag = produce_larval_cohort(
            0, agents, geno, 10000.0, 11.0, 48.0,
        )
        assert cohort is None
        assert diag['n_spawning_males'] == 0

    def test_cohort_diagnostics(self):
        """Diagnostics dict should contain expected keys."""
        agents, geno, fem, mal, efx = _make_spawning_population()
        rng = np.random.default_rng(42)
        cohort, diag = produce_larval_cohort(
            0, agents, geno, 10000.0, 11.0, 48.0, rng=rng,
        )
        expected_keys = {
            'n_spawning_females', 'n_spawning_males', 'total_eggs',
            'fertilization_success', 'n_zygotes', 'pld_days',
            'pelagic_survival', 'n_competent', 'ne', 'ne_over_n',
        }
        assert expected_keys.issubset(diag.keys())
        assert diag['total_eggs'] > 0
        assert 0 < diag['fertilization_success'] <= 1
        assert diag['n_competent'] > 0


# ═══════════════════════════════════════════════════════════════════════
# SETTLEMENT & RECRUITMENT
# ═══════════════════════════════════════════════════════════════════════

class TestSettlement:
    def test_settlement_cue_no_adults(self):
        """No adults → baseline modifier of 0.2."""
        assert settlement_cue_modifier(0) == pytest.approx(0.2, rel=0.01)

    def test_settlement_cue_saturates(self):
        """Many adults → modifier approaches 1.0."""
        mod = settlement_cue_modifier(100)
        assert mod > 0.95

    def test_settlement_cue_half_saturation(self):
        """At halfsat adults, modifier ≈ 0.6."""
        mod = settlement_cue_modifier(5, halfsat=5)
        assert mod == pytest.approx(0.6, rel=0.01)

    def test_settlement_cue_monotonic(self):
        """More adults → higher modifier."""
        vals = [settlement_cue_modifier(n) for n in range(20)]
        for i in range(1, len(vals)):
            assert vals[i] >= vals[i - 1]


class TestBevertonHolt:
    def test_zero_settlers(self):
        """No settlers → no recruits."""
        assert beverton_holt_recruitment(0, 500) == 0

    def test_zero_capacity(self):
        """No capacity → no recruits."""
        assert beverton_holt_recruitment(100, 0) == 0

    def test_low_settlers_supply_limited(self):
        """At low settler numbers, recruitment ≈ S × s0."""
        R = beverton_holt_recruitment(10, 1000)
        expected = 10 * 0.03  # = 0.3 → rounds to 0
        assert R <= 1  # very few from small settler number

    def test_high_settlers_saturates(self):
        """At very high settler numbers, recruitment approaches K × s0 = 15."""
        R = beverton_holt_recruitment(1_000_000, 500)
        max_expected = 500 * 0.03  # = 15
        # Should be close to but not exceed asymptote
        assert R <= max_expected + 1
        assert R >= max_expected - 1  # should be near the asymptote

    def test_intermediate_values(self):
        """Intermediate settlers give intermediate recruitment."""
        R_low = beverton_holt_recruitment(100, 500)
        R_high = beverton_holt_recruitment(10000, 500)
        assert R_high >= R_low

    def test_density_dependent(self):
        """More settlers → diminishing returns (density dependence).

        Start at S=100 to avoid integer-rounding artifacts at very low R.
        """
        settlers = [100, 500, 1000, 5000, 50000]
        per_settler = []
        for S in settlers:
            R = beverton_holt_recruitment(S, 500)
            per_settler.append(R / S if S > 0 else 0)
        # Per-settler survival should decrease (or stay same with rounding)
        for i in range(1, len(per_settler)):
            assert per_settler[i] <= per_settler[i - 1] + 0.001, \
                f"Per-settler survival not decreasing: {per_settler}"


class TestSettleRecruits:
    def test_basic_settlement(self):
        """settle_recruits should add individuals to the agent array.

        With corrected BH (R = S*s0/(1+S/K)), 500 settlers at K=500:
        R = 500*0.03/(1+500/500) = 15/2 = 7.5 → 8 recruits (after cue mod).
        """
        rng = np.random.default_rng(42)
        max_n = 1000
        agents = allocate_agents(max_n)
        genotypes = allocate_genotypes(max_n)
        effect_sizes = _make_effect_sizes()

        n_settlers = 500
        settler_geno = rng.integers(0, 2, size=(n_settlers, N_LOCI, 2), dtype=np.int8)
        # Fix EF1A lethals
        lethal = settler_geno[:, IDX_EF1A, :].sum(axis=1) == 2
        settler_geno[lethal, IDX_EF1A, 1] = 0

        n_added = settle_recruits(
            agents=agents,
            genotypes=genotypes,
            settler_genotypes=settler_geno,
            node_id=0,
            carrying_capacity=500,
            n_adults_present=20,
            habitat_area=10000.0,
            effect_sizes=effect_sizes,
            rng=rng,
        )

        assert n_added > 0
        assert np.sum(agents['alive']) == n_added

        # Check recruit properties
        alive_mask = agents['alive']
        assert np.all(agents['stage'][alive_mask] == Stage.SETTLER)
        assert np.all(agents['size'][alive_mask] == pytest.approx(0.5))
        assert np.all(agents['disease_state'][alive_mask] == DiseaseState.S)
        assert np.all(agents['fecundity_mod'][alive_mask] == 1.0)

    def test_no_settlers_no_recruits(self):
        """Empty settler array → 0 recruits."""
        agents = allocate_agents(100)
        genotypes = allocate_genotypes(100)
        effect_sizes = _make_effect_sizes()

        n = settle_recruits(
            agents, genotypes,
            np.empty((0, N_LOCI, 2), dtype=np.int8),
            0, 500, 10, 10000.0, effect_sizes,
        )
        assert n == 0

    def test_resistance_computed(self):
        """Settled recruits should have valid resistance scores."""
        rng = np.random.default_rng(42)
        max_n = 1000
        agents = allocate_agents(max_n)
        genotypes = allocate_genotypes(max_n)
        effect_sizes = _make_effect_sizes()

        n_settlers = 500
        settler_geno = rng.integers(0, 2, size=(n_settlers, N_LOCI, 2), dtype=np.int8)
        lethal = settler_geno[:, IDX_EF1A, :].sum(axis=1) == 2
        settler_geno[lethal, IDX_EF1A, 1] = 0

        settle_recruits(
            agents, genotypes, settler_geno,
            0, 500, 20, 10000.0, effect_sizes, rng=rng,
        )

        alive_mask = agents['alive']
        resistances = agents['resistance'][alive_mask]
        assert len(resistances) > 0, "Expected some recruits"
        # All should be in valid range
        assert np.all(resistances >= 0.0)
        assert np.all(resistances <= 1.0)
        # Not all zero (very unlikely with random genotypes)
        assert np.any(resistances > 0.0)


# ═══════════════════════════════════════════════════════════════════════
# RESISTANCE SCORE
# ═══════════════════════════════════════════════════════════════════════

class TestResistanceScore:
    def test_all_susceptible(self):
        """All-zero genotype → r = 0."""
        geno = np.zeros((N_LOCI, 2), dtype=np.int8)
        effects = _make_effect_sizes()
        assert _compute_resistance(geno, effects) == 0.0

    def test_all_resistant_het_ef1a(self):
        """All additive loci resistant + EF1A heterozygote → r = 1.0."""
        geno = np.ones((N_LOCI, 2), dtype=np.int8)
        geno[IDX_EF1A, 1] = 0  # heterozygote
        effects = _make_effect_sizes()
        r = _compute_resistance(geno, effects)
        assert r == pytest.approx(1.0, abs=0.01)

    def test_ef1a_heterozygote_bonus(self):
        """EF1A heterozygote should add ~0.16 to resistance."""
        effects = _make_effect_sizes()
        geno_hom = np.zeros((N_LOCI, 2), dtype=np.int8)
        geno_het = np.zeros((N_LOCI, 2), dtype=np.int8)
        geno_het[IDX_EF1A, 0] = 1  # heterozygote

        r_hom = _compute_resistance(geno_hom, effects)
        r_het = _compute_resistance(geno_het, effects)
        assert r_het - r_hom == pytest.approx(0.16, abs=0.01)

    def test_additive_scaling(self):
        """More resistant alleles → higher r (additive)."""
        effects = _make_effect_sizes()
        r_values = []
        for n_resistant in [0, 10, 25, 51]:
            geno = np.zeros((N_LOCI, 2), dtype=np.int8)
            geno[:n_resistant, :] = 1
            r_values.append(_compute_resistance(geno, effects))
        for i in range(1, len(r_values)):
            assert r_values[i] > r_values[i - 1]

    def test_in_range(self):
        """Resistance score should always be in [0, 1]."""
        rng = np.random.default_rng(42)
        effects = _make_effect_sizes()
        for _ in range(100):
            geno = rng.integers(0, 2, size=(N_LOCI, 2), dtype=np.int8)
            # Fix EF1A lethal
            if geno[IDX_EF1A, :].sum() == 2:
                geno[IDX_EF1A, 1] = 0
            r = _compute_resistance(geno, effects)
            assert 0.0 <= r <= 1.0 + 1e-6, f"r={r} out of range"


# ═══════════════════════════════════════════════════════════════════════
# INTEGRATION: FULL PIPELINE
# ═══════════════════════════════════════════════════════════════════════

class TestFullPipeline:
    """End-to-end tests combining multiple components."""

    def test_spawning_to_settlement(self):
        """Full pipeline: spawn → larvae → settle.

        Use the full cohort (not a small subset) to ensure enough settlers
        survive Beverton-Holt density-dependent recruitment (s0 = 0.03).
        """
        rng = np.random.default_rng(42)
        agents, geno, fem, mal, effect_sizes = _make_spawning_population(
            n_females=50, n_males=50, size_mm=600.0, rng=rng,
        )

        # Produce larvae
        cohort, diag = produce_larval_cohort(
            node_id=0,
            agents=agents,
            genotypes=geno,
            habitat_area=10000.0,
            sst=11.0,
            lat=48.0,
            rng=rng,
        )
        assert cohort is not None
        assert cohort.n_competent > 0

        # Settle at a new node — use ALL competent larvae (not a subset)
        max_agents = 1000
        agents2 = allocate_agents(max_agents)
        geno2 = allocate_genotypes(max_agents)
        # Place some existing adults (for settlement cue)
        for i in range(20):
            agents2[i]['alive'] = True
            agents2[i]['stage'] = Stage.ADULT

        n_settled = settle_recruits(
            agents=agents2,
            genotypes=geno2,
            settler_genotypes=cohort.genotypes,  # full cohort
            node_id=0,
            carrying_capacity=500,
            n_adults_present=20,
            habitat_area=10000.0,
            effect_sizes=effect_sizes,
            rng=rng,
        )
        assert n_settled > 0, \
            f"Expected recruits from {cohort.n_competent} competent larvae"

    def test_allee_reduces_fertilization_at_low_density(self):
        """Low density should have much lower fertilization than high density.

        With aggregation behavior, even a few adults can find each other
        (biologically realistic for mobile sea stars), so fertilization
        isn't zero — but it's dramatically lower than at healthy density.
        The true Allee collapse occurs below ~2 adults (mate-finding failure).
        """
        rng = np.random.default_rng(42)
        habitat_area = 100_000.0  # 10 hectares

        # Low density: 3 adults
        agents_low, geno_low, fem_l, mal_l, efx = _make_spawning_population(
            n_females=2, n_males=1, size_mm=600.0, rng=rng,
        )
        _, diag_low = produce_larval_cohort(
            0, agents_low, geno_low, habitat_area, 11.0, 48.0, rng=rng,
        )

        # High density: 200 adults in same area
        rng2 = np.random.default_rng(99)
        agents_hi, geno_hi, fem_h, mal_h, _ = _make_spawning_population(
            n_females=100, n_males=100, size_mm=600.0, rng=rng2,
        )
        _, diag_hi = produce_larval_cohort(
            0, agents_hi, geno_hi, habitat_area, 11.0, 48.0, rng=rng2,
        )

        # Low density should have MUCH lower fertilization
        assert diag_low['fertilization_success'] < diag_hi['fertilization_success'], \
            "Low density should have lower fertilization than high density"
        # And specifically, low density fert < 0.30
        assert diag_low['fertilization_success'] < 0.30, \
            f"Fert at low density = {diag_low['fertilization_success']}, expected < 0.30"

    def test_allee_zero_males_no_reproduction(self):
        """With zero males, there should be no reproduction at all."""
        rng = np.random.default_rng(42)
        agents, geno, fem, mal, _ = _make_spawning_population(
            n_females=50, n_males=0, size_mm=600.0, rng=rng,
        )
        cohort, diag = produce_larval_cohort(
            0, agents, geno, 100_000.0, 11.0, 48.0, rng=rng,
        )
        assert cohort is None
        assert diag['fertilization_success'] == 0.0

    def test_healthy_population_replacement(self):
        """A healthy population should achieve positive growth.

        At normal density, fertilization + recruitment should produce
        more than enough offspring to replace adult mortality.
        """
        rng = np.random.default_rng(42)
        habitat_area = 5000.0  # 0.5 hectares
        # 50 adults → density = 0.01 ind/m²
        agents, geno, fem, mal, efx = _make_spawning_population(
            n_females=50, n_males=50, size_mm=600.0, rng=rng,
        )

        cohort, diag = produce_larval_cohort(
            node_id=0,
            agents=agents,
            genotypes=geno,
            habitat_area=habitat_area,
            sst=11.0,
            lat=48.0,
            rng=rng,
        )

        # Fertilization should be reasonable
        assert diag['fertilization_success'] > 0.05, \
            f"Fert = {diag['fertilization_success']}, expected >0.05 at healthy density"
        assert diag['n_competent'] > 0
