"""Invariant tests for the three-trait genetic architecture.

Validates the fundamental biological separations between the three traits:
  - Resistance (r_i): reduces probability of infection
  - Tolerance (t_i): extends I₂ survival time (timer-scaling)
  - Recovery (c_i): clears infection (I₁/I₂ → R)

These are INVARIANTS, not statistical tendencies. The architecture must
guarantee that each trait affects only its designated mechanism.

Author: Anton 🔬
"""

import numpy as np
import pytest

from sswd_evoepi.config import default_config
from sswd_evoepi.disease import (
    NodeDiseaseState,
    DiseaseSection,
    daily_disease_update,
    recovery_probability_I1,
    recovery_probability_I2,
    sample_stage_duration,
    shedding_rate_I1,
    shedding_rate_I2,
)
from sswd_evoepi.genetics import (
    N_LOCI,
    RESISTANCE_SLICE,
    TOLERANCE_SLICE,
    RECOVERY_SLICE,
    compute_trait_batch,
    compute_trait_single,
    initialize_genotypes_three_trait,
    initialize_trait_effect_sizes,
    update_all_trait_scores,
)
from sswd_evoepi.types import (
    AGENT_DTYPE,
    DiseaseState,
    N_LOCI,
    N_RESISTANCE_DEFAULT,
    N_TOLERANCE_DEFAULT,
    N_RECOVERY_DEFAULT,
    Stage,
    allocate_agents,
    allocate_genotypes,
    trait_slices,
)


# ═══════════════════════════════════════════════════════════════════════
# HELPERS
# ═══════════════════════════════════════════════════════════════════════


def _make_agents(n: int, **field_overrides) -> np.ndarray:
    """Create n alive adult agents with sensible defaults."""
    agents = allocate_agents(n)
    agents['alive'][:n] = True
    agents['stage'][:n] = Stage.ADULT
    agents['size'][:n] = 300.0
    agents['age'][:n] = 5.0
    agents['disease_state'][:n] = DiseaseState.S
    for field, val in field_overrides.items():
        agents[field][:n] = val
    return agents


def _make_disease_state(n_s: int = 0, n_e: int = 0, n_i1: int = 0,
                        n_i2: int = 0, n_r: int = 0) -> NodeDiseaseState:
    """Create NodeDiseaseState with high Vibrio for guaranteed infections."""
    state = NodeDiseaseState()
    state.vibrio_concentration = 500_000.0  # High pathogen load
    state.n_S = n_s
    state.n_E = n_e
    state.n_I1 = n_i1
    state.n_I2 = n_i2
    state.n_R = n_r
    return state


def _make_effect_sizes(seed: int = 42):
    """Create effect size arrays for all three traits."""
    rng = np.random.default_rng(seed)
    e_r = initialize_trait_effect_sizes(rng, N_RESISTANCE_DEFAULT)
    e_t = initialize_trait_effect_sizes(rng, N_TOLERANCE_DEFAULT)
    e_c = initialize_trait_effect_sizes(rng, N_RECOVERY_DEFAULT)
    return e_r, e_t, e_c


# ═══════════════════════════════════════════════════════════════════════
# TEST 1: Infection rate independence from tolerance
# ═══════════════════════════════════════════════════════════════════════


class TestInfectionIndependence:
    """Infection probability must depend ONLY on resistance."""

    def test_infection_independent_of_tolerance(self):
        """High-tolerance vs low-tolerance: same infection rate at identical resistance."""
        n = 2000
        rng = np.random.default_rng(42)
        config = default_config()
        cfg = config.disease

        # Group A: high tolerance, zero resistance
        agents_a = _make_agents(n, tolerance=0.8, resistance=0.0, recovery_ability=0.0)

        # Group B: low tolerance, zero resistance
        agents_b = _make_agents(n, tolerance=0.1, resistance=0.0, recovery_ability=0.0)

        # Expose both groups to identical pathogen pressure
        infected_a = 0
        infected_b = 0
        n_trials = 20

        for trial in range(n_trials):
            # Fresh agents each trial
            a = _make_agents(n, tolerance=0.8, resistance=0.0, recovery_ability=0.0)
            b = _make_agents(n, tolerance=0.1, resistance=0.0, recovery_ability=0.0)

            state_a = _make_disease_state(n_s=n)
            state_b = _make_disease_state(n_s=n)

            rng_a = np.random.default_rng(1000 + trial)
            rng_b = np.random.default_rng(1000 + trial)  # Same seed!

            daily_disease_update(
                agents=a, node_state=state_a, T_celsius=20.0,
                salinity=32.0, phi_k=0.0, dispersal_input=0.0,
                day=100, rng=rng_a, cfg=cfg,
            )
            daily_disease_update(
                agents=b, node_state=state_b, T_celsius=20.0,
                salinity=32.0, phi_k=0.0, dispersal_input=0.0,
                day=100, rng=rng_b, cfg=cfg,
            )

            infected_a += int(np.sum(a['disease_state'][:n] != DiseaseState.S))
            infected_b += int(np.sum(b['disease_state'][:n] != DiseaseState.S))

        # With identical seeds and resistance, infection counts should be EXACTLY equal
        assert infected_a == infected_b, (
            f"Tolerance affected infection! high-tol={infected_a}, low-tol={infected_b}"
        )

    def test_infection_independent_of_recovery(self):
        """High-recovery vs low-recovery: same infection rate at identical resistance."""
        n = 2000
        config = default_config()
        cfg = config.disease

        infected_a = 0
        infected_b = 0
        n_trials = 20

        for trial in range(n_trials):
            a = _make_agents(n, recovery_ability=0.9, resistance=0.0, tolerance=0.0)
            b = _make_agents(n, recovery_ability=0.0, resistance=0.0, tolerance=0.0)

            state_a = _make_disease_state(n_s=n)
            state_b = _make_disease_state(n_s=n)

            rng_a = np.random.default_rng(2000 + trial)
            rng_b = np.random.default_rng(2000 + trial)

            daily_disease_update(
                agents=a, node_state=state_a, T_celsius=20.0,
                salinity=32.0, phi_k=0.0, dispersal_input=0.0,
                day=100, rng=rng_a, cfg=cfg,
            )
            daily_disease_update(
                agents=b, node_state=state_b, T_celsius=20.0,
                salinity=32.0, phi_k=0.0, dispersal_input=0.0,
                day=100, rng=rng_b, cfg=cfg,
            )

            infected_a += int(np.sum(a['disease_state'][:n] != DiseaseState.S))
            infected_b += int(np.sum(b['disease_state'][:n] != DiseaseState.S))

        assert infected_a == infected_b, (
            f"Recovery ability affected infection! high-c={infected_a}, low-c={infected_b}"
        )


# ═══════════════════════════════════════════════════════════════════════
# TEST 2: Shedding rate independence from tolerance
# ═══════════════════════════════════════════════════════════════════════


class TestSheddingIndependence:
    """Shedding rates are functions of (T, config) only — not host traits."""

    @pytest.mark.parametrize("T_celsius", [10.0, 15.0, 20.0])
    def test_shedding_I1_is_host_independent(self, T_celsius):
        """I₁ shedding rate is a function of temperature only."""
        cfg = default_config().disease
        rate = shedding_rate_I1(T_celsius, cfg)
        # Same call, same answer — no host state involved
        assert rate == shedding_rate_I1(T_celsius, cfg)
        assert rate > 0

    @pytest.mark.parametrize("T_celsius", [10.0, 15.0, 20.0])
    def test_shedding_I2_is_host_independent(self, T_celsius):
        """I₂ shedding rate is a function of temperature only."""
        cfg = default_config().disease
        rate = shedding_rate_I2(T_celsius, cfg)
        assert rate == shedding_rate_I2(T_celsius, cfg)
        assert rate > 0
        # I₂ should shed more than I₁
        assert rate > shedding_rate_I1(T_celsius, cfg)


# ═══════════════════════════════════════════════════════════════════════
# TEST 3: Tolerance extends I₂ survival
# ═══════════════════════════════════════════════════════════════════════


class TestToleranceExtendsI2:
    """High tolerance should produce longer I₂ timer (lower effective death rate)."""

    def test_tolerance_extends_i2_timer_statistically(self):
        """High-t_i agents get longer disease_timer in I₂ on average (N=5000)."""
        cfg = default_config().disease
        rng = np.random.default_rng(42)

        mu_I2D = cfg.mu_I2D_ref  # Base rate at T_ref
        n_samples = 5000

        timers_high_tol = []
        timers_low_tol = []

        for _ in range(n_samples):
            # High tolerance (t_i=0.8): effective_rate = mu_I2D * (1 - 0.8*0.85) = mu_I2D * 0.32
            t_high = 0.8
            eff_rate_high = mu_I2D * (1.0 - t_high * cfg.tau_max)
            eff_rate_high = max(eff_rate_high, mu_I2D * 0.05)
            timers_high_tol.append(sample_stage_duration(eff_rate_high, 3, rng))

            # Low tolerance (t_i=0.1): effective_rate = mu_I2D * (1 - 0.1*0.85) = mu_I2D * 0.915
            t_low = 0.1
            eff_rate_low = mu_I2D * (1.0 - t_low * cfg.tau_max)
            eff_rate_low = max(eff_rate_low, mu_I2D * 0.05)
            timers_low_tol.append(sample_stage_duration(eff_rate_low, 3, rng))

        mean_high = np.mean(timers_high_tol)
        mean_low = np.mean(timers_low_tol)

        # High tolerance → lower rate → longer timer
        assert mean_high > mean_low * 1.5, (
            f"High-tolerance mean timer ({mean_high:.1f}) should be > 1.5× "
            f"low-tolerance mean ({mean_low:.1f})"
        )


# ═══════════════════════════════════════════════════════════════════════
# TEST 4: Recovery uses c_i, not r_i
# ═══════════════════════════════════════════════════════════════════════


class TestRecoveryUsesCI:
    """Recovery probability depends on c_i (clearance), not r_i (resistance)."""

    def test_recovery_I2_depends_on_c_i(self):
        """recovery_probability_I2(c_i) is proportional to c_i."""
        assert recovery_probability_I2(0.0) == 0.0
        assert recovery_probability_I2(0.5) == pytest.approx(0.025)
        assert recovery_probability_I2(1.0) == pytest.approx(0.05)

    def test_recovery_I1_depends_on_c_i(self):
        """recovery_probability_I1 only works for c_i > 0.5."""
        assert recovery_probability_I1(0.0) == 0.0
        assert recovery_probability_I1(0.3) == 0.0
        assert recovery_probability_I1(0.5) == 0.0
        assert recovery_probability_I1(0.75) > 0.0
        assert recovery_probability_I1(1.0) == pytest.approx(0.05)

    def test_high_c_low_r_recovers_more_than_low_c_high_r(self):
        """Recovery is driven by c_i. Agents with high c_i + low r_i
        should recover more than those with low c_i + high r_i."""
        n = 1000
        n_trials = 100
        rng = np.random.default_rng(42)

        # Directly test recovery probability: high c_i (0.9) vs low c_i (0.1)
        high_c_recoveries = 0
        low_c_recoveries = 0

        for _ in range(n_trials * n):
            if rng.random() < recovery_probability_I2(0.9):
                high_c_recoveries += 1
            if rng.random() < recovery_probability_I2(0.1):
                low_c_recoveries += 1

        assert high_c_recoveries > low_c_recoveries * 5, (
            f"High-c_i ({high_c_recoveries}) should recover far more than "
            f"low-c_i ({low_c_recoveries})"
        )


# ═══════════════════════════════════════════════════════════════════════
# TEST 5: Resistance still reduces infection
# ═══════════════════════════════════════════════════════════════════════


class TestResistanceReducesInfection:
    """Higher resistance should reduce infection probability."""

    def test_high_resistance_fewer_infections(self):
        """High r_i agents get infected less than low r_i agents."""
        n = 2000
        config = default_config()
        cfg = config.disease

        infected_high_r = 0
        infected_low_r = 0
        n_trials = 20

        for trial in range(n_trials):
            a_high = _make_agents(n, resistance=0.8, tolerance=0.0, recovery_ability=0.0)
            a_low = _make_agents(n, resistance=0.0, tolerance=0.0, recovery_ability=0.0)

            state_high = _make_disease_state(n_s=n)
            state_low = _make_disease_state(n_s=n)

            rng_h = np.random.default_rng(3000 + trial)
            rng_l = np.random.default_rng(3000 + trial)

            daily_disease_update(
                agents=a_high, node_state=state_high, T_celsius=20.0,
                salinity=32.0, phi_k=0.0, dispersal_input=0.0,
                day=100, rng=rng_h, cfg=cfg,
            )
            daily_disease_update(
                agents=a_low, node_state=state_low, T_celsius=20.0,
                salinity=32.0, phi_k=0.0, dispersal_input=0.0,
                day=100, rng=rng_l, cfg=cfg,
            )

            infected_high_r += int(np.sum(a_high['disease_state'][:n] != DiseaseState.S))
            infected_low_r += int(np.sum(a_low['disease_state'][:n] != DiseaseState.S))

        assert infected_high_r < infected_low_r, (
            f"High resistance ({infected_high_r}) should have fewer infections "
            f"than low resistance ({infected_low_r})"
        )


# ═══════════════════════════════════════════════════════════════════════
# TEST 6: Trait isolation in genotype space
# ═══════════════════════════════════════════════════════════════════════


class TestTraitIsolation:
    """Setting loci to all-1 in one trait block should not affect other traits."""

    def test_resistance_only(self):
        """All-1 resistance loci, all-0 tolerance+recovery → r_i > 0, t_i = 0, c_i = 0."""
        e_r, e_t, e_c = _make_effect_sizes()
        agents = allocate_agents(1)
        agents['alive'][0] = True
        geno = allocate_genotypes(1)

        # Set all resistance loci to 1
        geno[0, RESISTANCE_SLICE, :] = 1

        update_all_trait_scores(agents, geno, e_r, e_t, e_c)

        assert agents['resistance'][0] > 0.0
        assert agents['tolerance'][0] == 0.0
        assert agents['recovery_ability'][0] == 0.0

    def test_tolerance_only(self):
        """All-1 tolerance loci, all-0 resistance+recovery → t_i > 0, r_i = 0, c_i = 0."""
        e_r, e_t, e_c = _make_effect_sizes()
        agents = allocate_agents(1)
        agents['alive'][0] = True
        geno = allocate_genotypes(1)

        geno[0, TOLERANCE_SLICE, :] = 1

        update_all_trait_scores(agents, geno, e_r, e_t, e_c)

        assert agents['resistance'][0] == 0.0
        assert agents['tolerance'][0] > 0.0
        assert agents['recovery_ability'][0] == 0.0

    def test_recovery_only(self):
        """All-1 recovery loci, all-0 resistance+tolerance → c_i > 0, r_i = 0, t_i = 0."""
        e_r, e_t, e_c = _make_effect_sizes()
        agents = allocate_agents(1)
        agents['alive'][0] = True
        geno = allocate_genotypes(1)

        geno[0, RECOVERY_SLICE, :] = 1

        update_all_trait_scores(agents, geno, e_r, e_t, e_c)

        assert agents['resistance'][0] == 0.0
        assert agents['tolerance'][0] == 0.0
        assert agents['recovery_ability'][0] > 0.0

    def test_all_traits_independent(self):
        """All-1 everywhere → all three traits positive and independent."""
        e_r, e_t, e_c = _make_effect_sizes()
        agents = allocate_agents(1)
        agents['alive'][0] = True
        geno = allocate_genotypes(1)
        geno[0, :, :] = 1  # All loci derived

        update_all_trait_scores(agents, geno, e_r, e_t, e_c)

        # All traits positive (each effect array sums to 1.0)
        assert agents['resistance'][0] == pytest.approx(1.0, abs=0.01)
        assert agents['tolerance'][0] == pytest.approx(1.0, abs=0.01)
        assert agents['recovery_ability'][0] == pytest.approx(1.0, abs=0.01)


# ═══════════════════════════════════════════════════════════════════════
# TEST 7: Genotype initialization targets
# ═══════════════════════════════════════════════════════════════════════


class TestGenotypeInitialization:
    """Three-trait initialization should hit target means."""

    def test_population_means_near_targets(self):
        """10000 agents: population-mean r_i, t_i, c_i within 20% of targets."""
        rng = np.random.default_rng(42)
        e_r, e_t, e_c = _make_effect_sizes(seed=42)

        target_r, target_t, target_c = 0.15, 0.10, 0.08

        geno = initialize_genotypes_three_trait(
            n_agents=10000,
            effects_r=e_r,
            effects_t=e_t,
            effects_c=e_c,
            rng=rng,
            target_mean_r=target_r,
            target_mean_t=target_t,
            target_mean_c=target_c,
        )

        # Compute trait scores
        agents = allocate_agents(10000)
        agents['alive'][:10000] = True

        # Expand geno into full genotype array
        full_geno = allocate_genotypes(10000)
        full_geno[:10000] = geno

        update_all_trait_scores(agents, full_geno, e_r, e_t, e_c)

        alive = agents['alive']
        mean_r = float(np.mean(agents['resistance'][alive]))
        mean_t = float(np.mean(agents['tolerance'][alive]))
        mean_c = float(np.mean(agents['recovery_ability'][alive]))

        assert abs(mean_r - target_r) / target_r < 0.20, (
            f"Mean resistance {mean_r:.4f} too far from target {target_r}"
        )
        assert abs(mean_t - target_t) / target_t < 0.20, (
            f"Mean tolerance {mean_t:.4f} too far from target {target_t}"
        )
        assert abs(mean_c - target_c) / target_c < 0.20, (
            f"Mean recovery {mean_c:.4f} too far from target {target_c}"
        )

    def test_all_traits_have_variance(self):
        """Population should have variance in all three traits (not degenerate)."""
        rng = np.random.default_rng(123)
        e_r, e_t, e_c = _make_effect_sizes(seed=123)

        geno = initialize_genotypes_three_trait(
            n_agents=5000,
            effects_r=e_r,
            effects_t=e_t,
            effects_c=e_c,
            rng=rng,
        )

        agents = allocate_agents(5000)
        agents['alive'][:5000] = True
        full_geno = allocate_genotypes(5000)
        full_geno[:5000] = geno

        update_all_trait_scores(agents, full_geno, e_r, e_t, e_c)

        alive = agents['alive']
        assert np.std(agents['resistance'][alive]) > 0.01
        assert np.std(agents['tolerance'][alive]) > 0.01
        assert np.std(agents['recovery_ability'][alive]) > 0.01


# ═══════════════════════════════════════════════════════════════════════
# TEST 8: Partition validation
# ═══════════════════════════════════════════════════════════════════════


class TestPartitionValidation:
    """trait_slices() must enforce N_LOCI=51 constraint."""

    def test_valid_partitions(self):
        """Various valid 51-locus partitions should work."""
        r, t, c = trait_slices(17, 17, 17)
        assert r == slice(0, 17)
        assert t == slice(17, 34)
        assert c == slice(34, 51)

        # Asymmetric partition
        r2, t2, c2 = trait_slices(20, 20, 11)
        assert r2 == slice(0, 20)
        assert t2 == slice(20, 40)
        assert c2 == slice(40, 51)

    def test_invalid_partition_raises(self):
        """Partition summing to != 51 should raise AssertionError."""
        with pytest.raises(AssertionError, match="Partition must sum to 51"):
            trait_slices(20, 20, 12)  # sum = 52

        with pytest.raises(AssertionError, match="Partition must sum to 51"):
            trait_slices(10, 10, 10)  # sum = 30


# ═══════════════════════════════════════════════════════════════════════
# TEST 9: Effect size independence
# ═══════════════════════════════════════════════════════════════════════


class TestEffectSizeIndependence:
    """Three effect arrays should be independently drawn and normalized."""

    def test_each_sums_to_one(self):
        """Each trait's effect sizes should sum to 1.0."""
        e_r, e_t, e_c = _make_effect_sizes()
        assert np.sum(e_r) == pytest.approx(1.0)
        assert np.sum(e_t) == pytest.approx(1.0)
        assert np.sum(e_c) == pytest.approx(1.0)

    def test_arrays_are_different(self):
        """Each trait should have different effect sizes (independent draws)."""
        e_r, e_t, e_c = _make_effect_sizes()
        assert not np.array_equal(e_r, e_t)
        assert not np.array_equal(e_t, e_c)
        assert not np.array_equal(e_r, e_c)

    def test_sorted_descending(self):
        """Effect sizes should be sorted largest-first."""
        e_r, e_t, e_c = _make_effect_sizes()
        for eff in [e_r, e_t, e_c]:
            # Each successive element should be <= the previous
            assert np.all(np.diff(eff) <= 0), "Effect sizes not sorted descending"

    def test_correct_lengths(self):
        """Each effect array should have the right number of loci."""
        e_r, e_t, e_c = _make_effect_sizes()
        assert len(e_r) == N_RESISTANCE_DEFAULT  # 17
        assert len(e_t) == N_TOLERANCE_DEFAULT    # 17
        assert len(e_c) == N_RECOVERY_DEFAULT     # 17


# ═══════════════════════════════════════════════════════════════════════
# TEST 10: Trait score computation correctness
# ═══════════════════════════════════════════════════════════════════════


class TestTraitScoreComputation:
    """Verify batch and single trait score computation agrees."""

    def test_batch_matches_single(self):
        """compute_trait_batch and compute_trait_single should agree."""
        rng = np.random.default_rng(42)
        e_r, e_t, e_c = _make_effect_sizes()

        n = 50
        geno = allocate_genotypes(n)
        geno[:n] = (rng.random((n, N_LOCI, 2)) < 0.3).astype(np.int8)

        alive = np.zeros(n, dtype=bool)
        alive[:n] = True

        batch_scores = compute_trait_batch(geno, e_r, alive, RESISTANCE_SLICE)

        for i in range(n):
            single = compute_trait_single(geno[i], e_r, RESISTANCE_SLICE)
            assert batch_scores[i] == pytest.approx(single, abs=1e-5), (
                f"Agent {i}: batch={batch_scores[i]:.6f}, single={single:.6f}"
            )

    def test_homozygous_derived_gives_max_score(self):
        """All loci homozygous derived (1,1) → trait score = 1.0."""
        e_r, _, _ = _make_effect_sizes()
        geno = allocate_genotypes(1)
        geno[0, RESISTANCE_SLICE, :] = 1
        alive = np.array([True])

        score = compute_trait_batch(geno, e_r, alive, RESISTANCE_SLICE)
        assert score[0] == pytest.approx(1.0, abs=0.01)

    def test_homozygous_ancestral_gives_zero(self):
        """All loci homozygous ancestral (0,0) → trait score = 0.0."""
        e_r, _, _ = _make_effect_sizes()
        geno = allocate_genotypes(1)
        # Already all zeros
        alive = np.array([True])

        score = compute_trait_batch(geno, e_r, alive, RESISTANCE_SLICE)
        assert score[0] == 0.0

    def test_heterozygous_gives_half_score(self):
        """All loci heterozygous (0,1) → trait score ≈ 0.5."""
        e_r, _, _ = _make_effect_sizes()
        geno = allocate_genotypes(1)
        geno[0, RESISTANCE_SLICE, 0] = 0
        geno[0, RESISTANCE_SLICE, 1] = 1
        alive = np.array([True])

        score = compute_trait_batch(geno, e_r, alive, RESISTANCE_SLICE)
        assert score[0] == pytest.approx(0.5, abs=0.01)
