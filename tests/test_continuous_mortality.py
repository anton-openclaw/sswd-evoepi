"""Tests for continuous (daily) mortality + growth functions (Phase 4).

Validates:
  1. Annual→daily mortality probability conversion accuracy
  2. Statistical equivalence between daily×365 and annual mortality
  3. No year-boundary discontinuity in coupled simulations
  4. Von Bertalanffy daily growth trajectory matches annual prediction
  5. Stage transitions via daily growth
  6. Senescence causes elevated mortality in old agents
  7. Cause-of-death stamping (NATURAL vs SENESCENCE)
  8. Coupled simulation runs end-to-end with daily demographics
  9. Spatial simulation runs end-to-end with daily demographics

References:
  - continuous-mortality-spec.md
  - population-dynamics-spec.md §2
"""

import numpy as np
import pytest

from sswd_evoepi.config import (
    PopulationSection,
    SimulationConfig,
    default_config,
)
from sswd_evoepi.model import (
    CoupledSimResult,
    annual_natural_mortality,
    daily_natural_mortality,
    daily_growth_and_aging,
    initialize_population,
    make_effect_sizes,
    run_coupled_simulation,
    run_spatial_simulation,
    von_bertalanffy,
)
from sswd_evoepi.spatial import (
    NodeDefinition,
    build_network,
)
from sswd_evoepi.types import (
    AGENT_DTYPE,
    ANNUAL_SURVIVAL,
    DeathCause,
    Stage,
    allocate_agents,
)


# ═══════════════════════════════════════════════════════════════════════
# FIXTURES
# ═══════════════════════════════════════════════════════════════════════

@pytest.fixture
def default_cfg() -> SimulationConfig:
    return default_config()


@pytest.fixture
def pop_cfg() -> PopulationSection:
    return PopulationSection()


@pytest.fixture
def effect_sizes() -> np.ndarray:
    return make_effect_sizes(12345)


def _make_alive_agents(n: int, stage: int, age: float, size: float = 100.0) -> np.ndarray:
    """Helper: create n alive agents with given stage, age, and size."""
    agents = allocate_agents(n)
    agents['alive'][:n] = True
    agents['stage'][:n] = stage
    agents['age'][:n] = age
    agents['size'][:n] = size
    agents['sex'][:n] = 0
    agents['disease_state'][:n] = 0  # S
    # fecundity_mod removed (three-trait architecture)
    return agents


def _make_3node_network() -> 'MetapopulationNetwork':
    """Minimal 3-node network for spatial sim tests."""
    nodes = [
        NodeDefinition(
            node_id=0, name="North", lat=57.0, lon=-135.0,
            subregion="AK-SE", habitat_area=20000.0,
            carrying_capacity=500, is_fjord=False,
            sill_depth=np.inf, flushing_rate=0.8,
            mean_sst=8.0, sst_amplitude=3.5, sst_trend=0.015,
            salinity=32.0, depth_range=(5.0, 60.0),
        ),
        NodeDefinition(
            node_id=1, name="Middle", lat=48.5, lon=-123.0,
            subregion="SS", habitat_area=20000.0,
            carrying_capacity=500, is_fjord=False,
            sill_depth=np.inf, flushing_rate=0.5,
            mean_sst=10.0, sst_amplitude=4.0, sst_trend=0.02,
            salinity=30.0, depth_range=(5.0, 40.0),
        ),
        NodeDefinition(
            node_id=2, name="South", lat=36.6, lon=-122.0,
            subregion="CA-BJ", habitat_area=20000.0,
            carrying_capacity=500, is_fjord=False,
            sill_depth=np.inf, flushing_rate=0.8,
            mean_sst=14.0, sst_amplitude=2.5, sst_trend=0.025,
            salinity=33.5, depth_range=(10.0, 60.0),
        ),
    ]
    return build_network(nodes, seed=42)


# ═══════════════════════════════════════════════════════════════════════
# TEST 1: Daily mortality probability conversion
# ═══════════════════════════════════════════════════════════════════════

class TestDailyMortalityProbConversion:
    """Verify annual→daily mortality conversion: 1-(1-p_daily)^365 ≈ p_annual."""

    def test_daily_mortality_probability_conversion(self, pop_cfg):
        """For each stage's annual_survival, daily→annual round-trip within 1%."""
        annual_surv = np.array(pop_cfg.annual_survival, dtype=np.float64)
        annual_mort = 1.0 - annual_surv

        for stage_idx, (s, m) in enumerate(zip(annual_surv, annual_mort)):
            # Convert to daily mortality
            daily_mort = 1.0 - s ** (1.0 / 365.0)
            # Convert back to annual mortality
            reconstructed_annual_mort = 1.0 - (1.0 - daily_mort) ** 365

            # Should match within 1%
            if m > 0.001:
                rel_err = abs(reconstructed_annual_mort - m) / m
                assert rel_err < 0.01, (
                    f"Stage {stage_idx}: annual_mort={m:.4f}, "
                    f"reconstructed={reconstructed_annual_mort:.4f}, "
                    f"rel_err={rel_err:.4f}"
                )
            else:
                # For very low mortality, absolute error check
                assert abs(reconstructed_annual_mort - m) < 0.001

    def test_conversion_all_stages(self, pop_cfg):
        """Verify the specific default values: [0.001, 0.03, 0.90, 0.95, 0.98]."""
        expected_survivals = [0.001, 0.03, 0.90, 0.95, 0.98]
        assert pop_cfg.annual_survival == expected_survivals

        for s in expected_survivals:
            daily_mort = 1.0 - s ** (1.0 / 365.0)
            # Daily mortality should be small positive number
            assert 0 < daily_mort < 1, f"Bad daily_mort={daily_mort} for S={s}"
            # And the compound should reconstruct the annual rate
            annual_mort_recon = 1.0 - (1.0 - daily_mort) ** 365
            assert annual_mort_recon == pytest.approx(1.0 - s, rel=0.01)


# ═══════════════════════════════════════════════════════════════════════
# TEST 2: Statistical equivalence (daily×365 ≈ annual)
# ═══════════════════════════════════════════════════════════════════════

class TestDailyMortalityStatisticalEquivalence:
    """daily_natural_mortality called 365× should kill ~same number as annual."""

    def test_daily_mortality_statistical_equivalence(self, pop_cfg, effect_sizes):
        """Mean deaths: daily×365 vs annual×1, within 10% over 50 seeds."""
        n_agents = 5000
        n_seeds = 50

        daily_deaths_all = []
        annual_deaths_all = []

        for seed in range(n_seeds):
            # --- Daily run ---
            agents_d, _ = initialize_population(
                n_individuals=n_agents,
                max_agents=n_agents + 100,
                habitat_area=50000.0,
                effect_sizes=effect_sizes,
                pop_cfg=pop_cfg,
                rng=np.random.default_rng(seed),
            )
            rng_d = np.random.default_rng(seed + 10000)
            dead_daily = 0
            for _ in range(365):
                n_k, _ = daily_natural_mortality(agents_d, pop_cfg, rng_d)
                dead_daily += n_k
            daily_deaths_all.append(dead_daily)

            # --- Annual run ---
            agents_a, _ = initialize_population(
                n_individuals=n_agents,
                max_agents=n_agents + 100,
                habitat_area=50000.0,
                effect_sizes=effect_sizes,
                pop_cfg=pop_cfg,
                rng=np.random.default_rng(seed),
            )
            rng_a = np.random.default_rng(seed + 10000)
            n_k_ann, _ = annual_natural_mortality(agents_a, pop_cfg, rng_a)
            annual_deaths_all.append(n_k_ann)

        mean_daily = np.mean(daily_deaths_all)
        mean_annual = np.mean(annual_deaths_all)

        # They should be within 10% of each other
        # (daily compounds slightly differently due to population depletion
        # during the year, but with mixed stages the effect should be modest)
        assert mean_daily > 0, "No daily deaths at all — something is broken"
        assert mean_annual > 0, "No annual deaths at all — something is broken"

        rel_diff = abs(mean_daily - mean_annual) / max(mean_daily, mean_annual)
        assert rel_diff < 0.10, (
            f"Daily mean={mean_daily:.1f}, Annual mean={mean_annual:.1f}, "
            f"rel_diff={rel_diff:.3f} > 10%"
        )


# ═══════════════════════════════════════════════════════════════════════
# TEST 3: No year-boundary discontinuity
# ═══════════════════════════════════════════════════════════════════════

class TestNoYearBoundaryDiscontinuity:
    """Coupled sim daily_pop should not show sudden drops at year boundaries."""

    def test_no_year_boundary_discontinuity(self):
        """5-year disease-free run: max single-day drop < 3% of pop."""
        result = run_coupled_simulation(
            n_individuals=500,
            carrying_capacity=500,
            habitat_area=50000.0,
            n_years=5,
            disease_year=999,  # no disease
            seed=42,
            record_daily=True,
        )

        daily_pop = result.daily_pop
        assert daily_pop is not None
        assert len(daily_pop) == 5 * 365

        # Check single-day drops
        for i in range(1, len(daily_pop)):
            if daily_pop[i - 1] > 0:
                drop_frac = (daily_pop[i - 1] - daily_pop[i]) / daily_pop[i - 1]
                assert drop_frac < 0.03, (
                    f"Day {i}: pop dropped from {daily_pop[i-1]} to "
                    f"{daily_pop[i]} ({drop_frac*100:.1f}%)"
                )

    def test_year_boundary_days_not_special(self):
        """Days 365, 730, etc. should not have outsized drops vs other days."""
        result = run_coupled_simulation(
            n_individuals=500,
            carrying_capacity=500,
            habitat_area=50000.0,
            n_years=4,
            disease_year=999,
            seed=123,
            record_daily=True,
        )
        daily_pop = result.daily_pop
        assert daily_pop is not None

        # Compute all daily deltas
        deltas = np.diff(daily_pop.astype(np.float64))
        # Year boundary days (0-indexed day 365, 730, 1095)
        boundary_days = [364, 729, 1094]
        non_boundary = [i for i in range(len(deltas)) if i not in boundary_days]

        # Mean drop at boundaries should not be dramatically worse than
        # the overall mean drop
        mean_all = np.mean(deltas[non_boundary])
        for bd in boundary_days:
            if bd < len(deltas):
                # Boundary drop shouldn't be more than 5× the average daily drop
                assert deltas[bd] > mean_all * 5 or abs(deltas[bd]) < 20, (
                    f"Day {bd}: delta={deltas[bd]}, mean_daily_delta={mean_all:.2f}"
                )


# ═══════════════════════════════════════════════════════════════════════
# TEST 4: VB daily growth trajectory
# ═══════════════════════════════════════════════════════════════════════

class TestDailyGrowthVBTrajectory:
    """Daily VB growth should track the analytical VB curve."""

    def test_daily_growth_vb_trajectory(self, pop_cfg):
        """Agent starting ON VB curve at age=2yr → after 365 steps ≈ VB(age=3)."""
        # Start on the VB curve so analytical prediction is valid
        start_size = von_bertalanffy(
            2.0, pop_cfg.L_inf, pop_cfg.k_growth, pop_cfg.t0_growth
        )
        agents = _make_alive_agents(1, Stage.JUVENILE, age=2.0, size=start_size)

        # Expected size at age 3 from VB curve
        expected_size_3yr = von_bertalanffy(
            3.0, pop_cfg.L_inf, pop_cfg.k_growth, pop_cfg.t0_growth
        )

        # Run daily growth 365 times with noise suppressed by averaging
        # across many seeds
        final_sizes = []
        for seed in range(100):
            test_agents = agents.copy()
            rng = np.random.default_rng(seed)
            for _ in range(365):
                daily_growth_and_aging(test_agents, pop_cfg, rng)
            final_sizes.append(float(test_agents['size'][0]))

        mean_final = np.mean(final_sizes)

        # Mean across seeds should match VB prediction within 5%
        rel_err = abs(mean_final - expected_size_3yr) / expected_size_3yr
        assert rel_err < 0.05, (
            f"Mean final size={mean_final:.1f}mm, "
            f"VB prediction={expected_size_3yr:.1f}mm, "
            f"rel_err={rel_err:.3f}"
        )

    def test_daily_growth_monotonic(self, pop_cfg):
        """Size should never decrease (no-shrinking invariant)."""
        agents = _make_alive_agents(1, Stage.JUVENILE, age=2.0, size=100.0)
        rng = np.random.default_rng(42)
        prev_size = float(agents['size'][0])
        for _ in range(365):
            daily_growth_and_aging(agents, pop_cfg, rng)
            curr_size = float(agents['size'][0])
            assert curr_size >= prev_size - 1e-6, (
                f"Size decreased: {prev_size:.4f} → {curr_size:.4f}"
            )
            prev_size = curr_size


# ═══════════════════════════════════════════════════════════════════════
# TEST 5: Stage transitions via daily growth
# ═══════════════════════════════════════════════════════════════════════

class TestStageTransitionsDaily:
    """Agents near stage boundaries should transition after enough growth."""

    def test_settler_to_juvenile(self, pop_cfg):
        """Agent at 8mm SETTLER → JUVENILE after enough daily growth."""
        agents = _make_alive_agents(1, Stage.SETTLER, age=0.5, size=8.0)
        rng = np.random.default_rng(42)
        # Grow for 365 days — at size 8mm near 10mm threshold, should cross
        for _ in range(365):
            daily_growth_and_aging(agents, pop_cfg, rng)
        assert agents['stage'][0] >= Stage.JUVENILE, (
            f"Expected JUVENILE+, got stage={agents['stage'][0]}, "
            f"size={agents['size'][0]:.1f}mm"
        )

    def test_juvenile_to_subadult(self, pop_cfg):
        """Agent at 140mm JUVENILE → SUBADULT after enough growth."""
        agents = _make_alive_agents(1, Stage.JUVENILE, age=4.0, size=140.0)
        rng = np.random.default_rng(42)
        # 150mm threshold — should cross within a year or two
        for _ in range(730):  # 2 years
            daily_growth_and_aging(agents, pop_cfg, rng)
        assert agents['stage'][0] >= Stage.SUBADULT, (
            f"Expected SUBADULT+, got stage={agents['stage'][0]}, "
            f"size={agents['size'][0]:.1f}mm"
        )

    def test_subadult_to_adult(self, pop_cfg):
        """Agent at 380mm SUBADULT → ADULT after enough growth."""
        agents = _make_alive_agents(1, Stage.SUBADULT, age=8.0, size=380.0)
        rng = np.random.default_rng(42)
        for _ in range(730):  # 2 years
            daily_growth_and_aging(agents, pop_cfg, rng)
        assert agents['stage'][0] == Stage.ADULT, (
            f"Expected ADULT, got stage={agents['stage'][0]}, "
            f"size={agents['size'][0]:.1f}mm"
        )

    def test_stages_never_demote(self, pop_cfg):
        """Stages should only increase, never decrease."""
        agents = _make_alive_agents(1, Stage.ADULT, age=20.0, size=600.0)
        rng = np.random.default_rng(42)
        for _ in range(365):
            daily_growth_and_aging(agents, pop_cfg, rng)
        assert agents['stage'][0] == Stage.ADULT


# ═══════════════════════════════════════════════════════════════════════
# TEST 6: Senescence causes elevated mortality
# ═══════════════════════════════════════════════════════════════════════

class TestSenescenceDaily:
    """Old agents (60yr) should die more than young agents (10yr)."""

    def test_senescence_daily(self, pop_cfg):
        """100 seeds: old agents survive less than young agents."""
        n_agents = 500
        n_seeds = 100
        n_days = 365

        young_survivors = []
        old_survivors = []

        for seed in range(n_seeds):
            # Young agents (10yr, ADULT)
            young = _make_alive_agents(n_agents, Stage.ADULT, age=10.0, size=600.0)
            rng_y = np.random.default_rng(seed)
            for _ in range(n_days):
                daily_natural_mortality(young, pop_cfg, rng_y)
            young_survivors.append(int(np.sum(young['alive'])))

            # Old agents (60yr, ADULT, well past senescence_age=50)
            old = _make_alive_agents(n_agents, Stage.ADULT, age=60.0, size=800.0)
            rng_o = np.random.default_rng(seed)
            for _ in range(n_days):
                daily_natural_mortality(old, pop_cfg, rng_o)
            old_survivors.append(int(np.sum(old['alive'])))

        mean_young_surv = np.mean(young_survivors)
        mean_old_surv = np.mean(old_survivors)

        # Old agents should survive significantly less
        assert mean_old_surv < mean_young_surv, (
            f"Old agents survived more ({mean_old_surv:.1f}) "
            f"than young ({mean_young_surv:.1f})"
        )
        # Quantitative: old agents (60yr) should have noticeable extra mortality
        # senescence_mortality=0.10, (60-50)/20 = 0.5 → extra 5% annual mort
        # vs adult baseline 2% annual mort → old should lose ~3.5× more
        assert mean_old_surv < mean_young_surv * 0.99, (
            "Senescence effect too small to detect"
        )


# ═══════════════════════════════════════════════════════════════════════
# TEST 7: Cause of death stamped correctly
# ═══════════════════════════════════════════════════════════════════════

class TestCauseOfDeathStamping:
    """Dead agents should have correct cause_of_death values."""

    def test_cause_of_death_stamped_correctly(self, pop_cfg):
        """Mixed-age population: dead agents get NATURAL or SENESCENCE."""
        n_young = 300
        n_old = 200
        total = n_young + n_old

        agents = allocate_agents(total)

        # Young adults (age 10)
        agents['alive'][:n_young] = True
        agents['stage'][:n_young] = Stage.ADULT
        agents['age'][:n_young] = 10.0
        agents['size'][:n_young] = 600.0
        # fecundity_mod removed (three-trait architecture)

        # Old adults (age 60, past senescence_age=50)
        agents['alive'][n_young:total] = True
        agents['stage'][n_young:total] = Stage.ADULT
        agents['age'][n_young:total] = 60.0
        agents['size'][n_young:total] = 800.0
        # fecundity_mod removed (three-trait architecture)

        rng = np.random.default_rng(42)

        # Run daily mortality for a year to get some deaths
        for _ in range(365):
            daily_natural_mortality(agents, pop_cfg, rng)

        dead_mask = ~agents['alive'][:total]

        # All dead agents should have cause_of_death != 0 (ALIVE)
        dead_causes = agents['cause_of_death'][:total][dead_mask]
        assert np.all(dead_causes != DeathCause.ALIVE), (
            "Some dead agents still have cause_of_death=ALIVE"
        )

        # Agents that were old (idx >= n_young) and died should have SENESCENCE
        old_dead = ~agents['alive'][n_young:total]
        if np.any(old_dead):
            old_causes = agents['cause_of_death'][n_young:total][old_dead]
            n_senescence = np.sum(old_causes == DeathCause.SENESCENCE)
            assert n_senescence > 0, (
                "No old agents died of SENESCENCE — all got NATURAL?"
            )

        # Young agents (idx < n_young) that died should have NATURAL
        young_dead = ~agents['alive'][:n_young]
        if np.any(young_dead):
            young_causes = agents['cause_of_death'][:n_young][young_dead]
            assert np.all(young_causes == DeathCause.NATURAL), (
                f"Young agents got non-NATURAL causes: "
                f"{np.unique(young_causes)}"
            )

    def test_alive_agents_have_alive_cause(self, pop_cfg):
        """Living agents should have cause_of_death = ALIVE (0)."""
        agents = _make_alive_agents(100, Stage.ADULT, age=10.0)
        rng = np.random.default_rng(42)
        daily_natural_mortality(agents, pop_cfg, rng)

        alive_mask = agents['alive'][:100]
        alive_causes = agents['cause_of_death'][:100][alive_mask]
        assert np.all(alive_causes == DeathCause.ALIVE)


# ═══════════════════════════════════════════════════════════════════════
# TEST 8: Coupled simulation with daily demographics
# ═══════════════════════════════════════════════════════════════════════

class TestCoupledSimDailyDemographics:
    """run_coupled_simulation should work end-to-end with daily mortality."""

    def test_coupled_sim_runs_with_daily_demographics(self):
        """10-year coupled sim completes; yearly_natural_deaths all > 0."""
        result = run_coupled_simulation(
            n_individuals=500,
            carrying_capacity=500,
            habitat_area=50000.0,
            n_years=10,
            disease_year=999,  # no disease for cleaner test
            seed=42,
        )

        assert isinstance(result, CoupledSimResult)
        assert result.n_years == 10
        assert result.yearly_natural_deaths is not None

        # Every year should have some natural deaths (population > 0)
        for yr in range(10):
            assert result.yearly_natural_deaths[yr] > 0, (
                f"Year {yr}: zero natural deaths"
            )

        # Final population should be reasonable (not zero, not exploded)
        assert result.final_pop > 0
        assert result.final_pop < 5000  # shouldn't 10× from K=500

    def test_coupled_sim_with_disease(self):
        """10-year sim with disease at year 3 completes without error."""
        result = run_coupled_simulation(
            n_individuals=500,
            carrying_capacity=500,
            habitat_area=50000.0,
            n_years=10,
            disease_year=3,
            initial_infected=5,
            seed=42,
        )
        assert result.n_years == 10
        assert result.total_disease_deaths >= 0
        assert result.total_natural_deaths > 0


# ═══════════════════════════════════════════════════════════════════════
# TEST 9: Spatial simulation with daily demographics
# ═══════════════════════════════════════════════════════════════════════

class TestSpatialSimDailyDemographics:
    """Spatial sim should work end-to-end with daily mortality/growth."""

    def test_spatial_sim_runs_with_daily_demographics(self):
        """3-node, 5-year spatial sim completes without errors."""
        network = _make_3node_network()
        config = default_config()

        result = run_spatial_simulation(
            network=network,
            n_years=5,
            disease_year=999,  # no disease
            seed=42,
            config=config,
        )

        assert result.n_years == 5
        assert result.n_nodes == 3
        assert result.yearly_natural_deaths is not None

        # Each node should have natural deaths each year
        for node_i in range(3):
            for yr in range(5):
                assert result.yearly_natural_deaths[node_i, yr] > 0, (
                    f"Node {node_i}, Year {yr}: zero natural deaths"
                )

        # Population should survive
        assert result.final_total_pop > 0
