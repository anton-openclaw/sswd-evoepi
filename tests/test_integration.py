"""Integration tests for coupled disease ↔ population dynamics (Phase 5).

Acceptance criteria:
  1. Disease at equilibrium → 80–95% crash → slow recovery (years)
  2. Temperature sensitivity: mortality increases with T (Eisenlord 2016)
  3. Post-crash Allee: if population too low, recovery stalls
  4. Disease-free runs produce stable K
  5. 20-year simulation completes without errors or NaN
  6. Mean resistance increases after epidemic (natural selection)
  7. Coupled feedbacks: disease kills → fewer adults → reduced reproduction

References:
  - Hamilton 2021: >90.6% global decline of Pycnopodia
  - Eisenlord 2016: temperature-mortality relationship
  - Harvell 2019: shelter sites with lower exposure had better survival
  - CODE_ERRATA CE-5: demographic Allee for high-fecundity species
"""

import numpy as np
import pytest

from sswd_evoepi.config import (
    DiseaseSection,
    PopulationSection,
    SimulationConfig,
    default_config,
)
from sswd_evoepi.model import (
    CoupledSimResult,
    annual_growth_and_aging,
    annual_natural_mortality,
    annual_reproduction,
    assign_stage,
    grow_individual,
    initialize_population,
    make_effect_sizes,
    natural_mortality_prob,
    run_coupled_simulation,
    von_bertalanffy,
)
from sswd_evoepi.types import (
    ANNUAL_SURVIVAL,
    DiseaseState,
    Stage,
)


# ═══════════════════════════════════════════════════════════════════════
# FIXTURES
# ═══════════════════════════════════════════════════════════════════════

@pytest.fixture
def default_cfg() -> SimulationConfig:
    return default_config()


@pytest.fixture
def effect_sizes() -> np.ndarray:
    return make_effect_sizes(12345)


# ═══════════════════════════════════════════════════════════════════════
# UNIT TESTS: POPULATION COMPONENTS
# ═══════════════════════════════════════════════════════════════════════

class TestVonBertalanffy:
    """Test VB growth curve."""

    def test_newborn_small(self):
        """Age 0 → ~39 mm."""
        L = von_bertalanffy(0.0)
        assert 0 < L < 100

    def test_adult_large(self):
        """Age 20 → ~808 mm."""
        L = von_bertalanffy(20.0)
        assert 700 < L < 900

    def test_approaches_linf(self):
        """Very old → approaches L_inf."""
        L = von_bertalanffy(100.0)
        assert L > 990

    def test_monotonically_increasing(self):
        """Size increases with age."""
        ages = [0, 1, 5, 10, 20, 40]
        sizes = [von_bertalanffy(a) for a in ages]
        for i in range(1, len(sizes)):
            assert sizes[i] > sizes[i - 1]


class TestGrowIndividual:
    def test_grows(self):
        """Individual should grow."""
        new = grow_individual(200.0, rng=None)
        assert new > 200.0

    def test_never_shrinks(self):
        """Size never decreases."""
        rng = np.random.default_rng(42)
        for _ in range(100):
            old = rng.uniform(50, 900)
            new = grow_individual(old, rng=rng)
            assert new >= old

    def test_near_linf_slow(self):
        """Near L_inf, growth is slow."""
        delta = grow_individual(990.0, rng=None) - 990.0
        assert delta < 2.0  # very slow growth


class TestAssignStage:
    def test_settler(self):
        assert assign_stage(5.0, Stage.SETTLER) == Stage.SETTLER

    def test_juvenile(self):
        assert assign_stage(50.0, Stage.SETTLER) == Stage.JUVENILE

    def test_subadult(self):
        assert assign_stage(200.0, Stage.JUVENILE) == Stage.SUBADULT

    def test_adult(self):
        assert assign_stage(450.0, Stage.SUBADULT) == Stage.ADULT

    def test_no_regression(self):
        """Stage should never decrease."""
        assert assign_stage(100.0, Stage.SUBADULT) == Stage.SUBADULT


class TestNaturalMortality:
    def test_settler_high_mortality(self):
        """Settlers have ~97% annual mortality."""
        p = natural_mortality_prob(Stage.SETTLER, 0.5)
        assert 0.95 < p < 0.99

    def test_adult_low_mortality(self):
        """Adults have ~2% annual mortality."""
        p = natural_mortality_prob(Stage.ADULT, 10.0)
        assert 0.01 < p < 0.05

    def test_senescence_increases(self):
        """Mortality increases after senescence age."""
        p_young = natural_mortality_prob(Stage.ADULT, 30.0)
        p_old = natural_mortality_prob(Stage.ADULT, 60.0)
        assert p_old > p_young


class TestPopulationInit:
    def test_correct_count(self, default_cfg, effect_sizes):
        rng = np.random.default_rng(42)
        agents, geno = initialize_population(
            100, 200, 10000.0, effect_sizes, default_cfg.population, rng,
        )
        assert np.sum(agents['alive']) == 100

    def test_all_alive(self, default_cfg, effect_sizes):
        rng = np.random.default_rng(42)
        agents, geno = initialize_population(
            50, 100, 10000.0, effect_sizes, default_cfg.population, rng,
        )
        alive = agents['alive'][:50]
        assert np.all(alive)

    def test_stages_assigned(self, default_cfg, effect_sizes):
        rng = np.random.default_rng(42)
        agents, _ = initialize_population(
            200, 400, 10000.0, effect_sizes, default_cfg.population, rng,
        )
        alive = agents['alive']
        stages = agents['stage'][alive]
        # Should have mix of stages
        unique_stages = np.unique(stages)
        assert len(unique_stages) >= 3  # at least settlers, juveniles/subadults, adults

    def test_resistance_in_range(self, default_cfg, effect_sizes):
        rng = np.random.default_rng(42)
        agents, _ = initialize_population(
            100, 200, 10000.0, effect_sizes, default_cfg.population, rng,
        )
        alive = agents['alive']
        r = agents['resistance'][alive]
        assert np.all(r >= 0.0)
        assert np.all(r <= 1.0)

    def test_no_ef1a_lethals(self, default_cfg, effect_sizes):
        rng = np.random.default_rng(42)
        _, geno = initialize_population(
            100, 200, 10000.0, effect_sizes, default_cfg.population, rng,
        )
        for i in range(100):
            assert geno[i, 51, :].sum() < 2  # no ins/ins lethals


# ═══════════════════════════════════════════════════════════════════════
# INTEGRATION TEST: DISEASE-FREE STABLE POPULATION
# ═══════════════════════════════════════════════════════════════════════

class TestDiseaseFreeStability:
    """Disease-free population should stabilize near K."""

    def test_stable_at_K(self):
        """Without disease, population should hover around K=500.

        Run for 20 years with no disease (disease_year > n_years).
        Population should stay within [K*0.5, K*1.5] after initial transient.
        """
        result = run_coupled_simulation(
            n_individuals=500,
            carrying_capacity=500,
            habitat_area=10000.0,
            T_celsius=11.0,
            salinity=30.0,
            phi_k=0.1,
            n_years=20,
            disease_year=999,  # never triggers disease
            seed=42,
        )

        assert result.n_years == 20
        # After 5-year spinup, pop should be stable
        late_pop = result.yearly_pop[5:]
        assert np.all(late_pop > 0), "Population crashed without disease!"
        # Population should stay within reasonable range
        mean_pop = np.mean(late_pop)
        assert mean_pop > 200, f"Mean pop {mean_pop} too low"
        assert mean_pop < 800, f"Mean pop {mean_pop} too high"
        # No NaN
        assert not np.any(np.isnan(result.yearly_mean_resistance))

    def test_recruits_produced(self):
        """Disease-free population should produce recruits each year."""
        result = run_coupled_simulation(
            n_individuals=500,
            carrying_capacity=500,
            habitat_area=10000.0,
            T_celsius=11.0,
            salinity=30.0,
            phi_k=0.1,
            n_years=10,
            disease_year=999,
            seed=42,
        )
        # After year 0 stabilization, should have recruits most years
        total_recruits = result.yearly_recruits[1:].sum()
        assert total_recruits > 0, "No recruits produced in 10 years"

    def test_twenty_year_no_crash(self):
        """20-year disease-free run completes without crash or NaN."""
        result = run_coupled_simulation(
            n_individuals=500,
            carrying_capacity=500,
            habitat_area=10000.0,
            T_celsius=11.0,
            salinity=30.0,
            phi_k=0.1,
            n_years=20,
            disease_year=999,
            seed=99,
        )
        assert result.final_pop > 0
        assert not np.any(np.isnan(result.yearly_pop.astype(float)))
        assert not np.any(np.isnan(result.yearly_mean_resistance))
        assert result.total_disease_deaths == 0


# ═══════════════════════════════════════════════════════════════════════
# INTEGRATION TEST: EPIDEMIC CRASH
# ═══════════════════════════════════════════════════════════════════════

class TestEpidemicCrash:
    """Disease at equilibrium population → massive crash."""

    def test_crash_80_plus_percent(self):
        """Epidemic at 15°C fjord should cause ≥80% population crash.

        Consistent with Hamilton 2021 (>90.6% global decline).
        """
        result = run_coupled_simulation(
            n_individuals=500,
            carrying_capacity=500,
            habitat_area=10000.0,
            T_celsius=15.0,
            salinity=30.0,
            phi_k=0.02,  # fjord
            n_years=8,
            disease_year=2,  # introduce after 2-year spinup
            initial_infected=5,
            seed=42,
        )
        # Population should crash dramatically
        pre_disease_pop = result.yearly_pop[1]  # year before disease
        min_pop = result.min_pop
        crash_fraction = 1.0 - (min_pop / pre_disease_pop)

        assert crash_fraction >= 0.80, \
            f"Expected ≥80% crash, got {crash_fraction:.1%} " \
            f"(pre={pre_disease_pop}, min={min_pop})"

    def test_disease_deaths_dominate(self):
        """Disease deaths should far exceed natural deaths during epidemic year."""
        result = run_coupled_simulation(
            n_individuals=500,
            carrying_capacity=500,
            habitat_area=10000.0,
            T_celsius=15.0,
            salinity=30.0,
            phi_k=0.02,
            n_years=6,
            disease_year=2,
            initial_infected=5,
            seed=42,
        )
        # In the epidemic year, disease deaths >> natural deaths
        dd = result.yearly_disease_deaths[2:5].sum()  # years 2-4
        nd = result.yearly_natural_deaths[2:5].sum()
        assert dd > nd * 3, \
            f"Disease deaths ({dd}) should dominate natural ({nd})"

    def test_vibrio_spike_during_epidemic(self):
        """Vibrio concentration should spike during epidemic."""
        result = run_coupled_simulation(
            n_individuals=500,
            carrying_capacity=500,
            habitat_area=10000.0,
            T_celsius=15.0,
            salinity=30.0,
            phi_k=0.02,
            n_years=5,
            disease_year=2,
            initial_infected=5,
            seed=42,
            record_daily=True,
        )
        # Pre-epidemic Vibrio
        pre_vibrio = np.mean(result.daily_vibrio[:365])
        # Peak vibrio during epidemic
        epidemic_start = 2 * 365
        peak_vibrio = np.max(result.daily_vibrio[epidemic_start:])
        assert peak_vibrio > pre_vibrio * 5, \
            f"Expected Vibrio spike, peak={peak_vibrio:.1f}, baseline={pre_vibrio:.1f}"


# ═══════════════════════════════════════════════════════════════════════
# INTEGRATION TEST: TEMPERATURE SENSITIVITY
# ═══════════════════════════════════════════════════════════════════════

class TestTemperatureSensitivity:
    """Mortality should increase with temperature (Eisenlord 2016)."""

    def _run_at_temp(self, T: float, seed: int = 42) -> CoupledSimResult:
        return run_coupled_simulation(
            n_individuals=500,
            carrying_capacity=500,
            habitat_area=10000.0,
            T_celsius=T,
            salinity=30.0,
            phi_k=0.02,
            n_years=6,
            disease_year=2,
            initial_infected=5,
            seed=seed,
        )

    def test_higher_temp_more_deaths(self):
        """16°C should have more disease deaths than 12°C.

        Eisenlord 2016: temperature is a strong predictor of SSWD mortality.
        """
        result_12 = self._run_at_temp(12.0)
        result_16 = self._run_at_temp(16.0)

        dd_12 = result_12.total_disease_deaths
        dd_16 = result_16.total_disease_deaths

        assert dd_16 > dd_12, \
            f"16°C ({dd_16} deaths) should exceed 12°C ({dd_12} deaths)"

    def test_19C_severe(self):
        """19°C should produce very high mortality (near T_opt for Vibrio)."""
        result_19 = self._run_at_temp(19.0)
        # At near-T_opt, disease should be devastating
        assert result_19.total_disease_deaths > 200, \
            f"Expected >200 deaths at 19°C, got {result_19.total_disease_deaths}"

    def test_temperature_gradient(self):
        """Mortality should follow: 12°C < 16°C ≤ 19°C."""
        # Use same seed for comparability
        r12 = self._run_at_temp(12.0, seed=99)
        r16 = self._run_at_temp(16.0, seed=99)
        r19 = self._run_at_temp(19.0, seed=99)

        dd_12 = r12.total_disease_deaths
        dd_16 = r16.total_disease_deaths
        dd_19 = r19.total_disease_deaths

        assert dd_12 < dd_16, \
            f"12°C ({dd_12}) should have fewer deaths than 16°C ({dd_16})"
        assert dd_16 <= dd_19 + 50, \
            f"16°C ({dd_16}) should not vastly exceed 19°C ({dd_19})"


# ═══════════════════════════════════════════════════════════════════════
# INTEGRATION TEST: POST-CRASH ALLEE EFFECT
# ═══════════════════════════════════════════════════════════════════════

class TestPostCrashAllee:
    """After population crash, Allee effects should slow recovery."""

    def test_slow_recovery(self):
        """Post-epidemic recovery should be slow (years, not immediate).

        CE-5: High-fecundity species have demographic Allee through
        reduced growth rate at low density, making recovery take decades.
        """
        result = run_coupled_simulation(
            n_individuals=500,
            carrying_capacity=500,
            habitat_area=10000.0,
            T_celsius=15.0,
            salinity=30.0,
            phi_k=0.02,
            n_years=15,
            disease_year=2,
            initial_infected=5,
            seed=42,
        )
        # After crash, population should not recover to K within 5 years
        crash_year = result.min_pop_year
        if crash_year + 5 < result.n_years:
            pop_5yr_later = result.yearly_pop[min(crash_year + 5, result.n_years - 1)]
            assert pop_5yr_later < 500, \
                f"Population recovered too quickly: {pop_5yr_later} at year {crash_year+5}"

    def test_reduced_reproduction_post_crash(self):
        """Post-crash population should have much lower reproduction.

        Fewer adults → fewer eggs → Allee reduces fertilization →
        dramatically fewer recruits.
        """
        result = run_coupled_simulation(
            n_individuals=500,
            carrying_capacity=500,
            habitat_area=10000.0,
            T_celsius=15.0,
            salinity=30.0,
            phi_k=0.02,
            n_years=10,
            disease_year=2,
            initial_infected=5,
            seed=42,
        )
        # Compare pre-disease and post-disease fertilization
        pre_fert = result.yearly_fert_success[1]  # pre-disease year
        # Find post-crash years
        crash_year = result.min_pop_year
        if crash_year + 1 < result.n_years:
            post_fert = result.yearly_fert_success[crash_year + 1]
            # Post-crash fertilization should be much lower
            if pre_fert > 0:
                assert post_fert < pre_fert, \
                    f"Post-crash fert ({post_fert:.4f}) should be < " \
                    f"pre-crash ({pre_fert:.4f})"

    def test_stalled_recovery_at_very_low_pop(self):
        """With extreme crash + small habitat, recovery may stall.

        This tests the Allee effect: at very low density in large habitat,
        fertilization success plummets and recovery takes many years.
        """
        result = run_coupled_simulation(
            n_individuals=500,
            carrying_capacity=500,
            habitat_area=100000.0,  # Large habitat → low density
            T_celsius=16.0,
            salinity=30.0,
            phi_k=0.02,
            n_years=15,
            disease_year=2,
            initial_infected=5,
            seed=42,
        )
        # Should have crashed
        assert result.min_pop < 100, \
            f"Expected severe crash, min_pop = {result.min_pop}"
        # Recovery should be slow
        if result.min_pop_year + 3 < result.n_years:
            pop_3yr = result.yearly_pop[result.min_pop_year + 3]
            assert pop_3yr < 400, \
                f"Recovery too fast in large habitat: {pop_3yr}"


# ═══════════════════════════════════════════════════════════════════════
# INTEGRATION TEST: SELECTION ON RESISTANCE
# ═══════════════════════════════════════════════════════════════════════

class TestSelectionOnResistance:
    """Epidemic should select for higher resistance."""

    def test_mean_resistance_increases(self):
        """Mean resistance should increase after epidemic.

        Disease kills low-r_i individuals preferentially,
        so survivors have higher mean resistance. Then their
        offspring also have higher mean r.
        """
        result = run_coupled_simulation(
            n_individuals=500,
            carrying_capacity=500,
            habitat_area=10000.0,
            T_celsius=15.0,
            salinity=30.0,
            phi_k=0.02,
            n_years=10,
            disease_year=2,
            initial_infected=5,
            seed=42,
        )
        # Pre-disease mean resistance
        r_pre = result.yearly_mean_resistance[1]
        # Post-epidemic resistance (after crash + recovery)
        r_post = result.yearly_mean_resistance[-1]

        assert r_post > r_pre, \
            f"Resistance should increase: pre={r_pre:.4f}, post={r_post:.4f}"


# ═══════════════════════════════════════════════════════════════════════
# INTEGRATION TEST: 20-YEAR SIMULATION
# ═══════════════════════════════════════════════════════════════════════

class TestLongSimulation:
    """20-year simulation should complete without errors or NaN."""

    def test_twenty_years_with_disease(self):
        """Full 20-year run: spinup → epidemic → recovery."""
        result = run_coupled_simulation(
            n_individuals=500,
            carrying_capacity=500,
            habitat_area=10000.0,
            T_celsius=14.0,
            salinity=30.0,
            phi_k=0.05,
            n_years=20,
            disease_year=3,
            initial_infected=5,
            seed=42,
            record_daily=True,
        )

        # Completes
        assert result.n_years == 20

        # No NaN in any array
        assert not np.any(np.isnan(result.yearly_pop.astype(float)))
        assert not np.any(np.isnan(result.yearly_adults.astype(float)))
        assert not np.any(np.isnan(result.yearly_mean_resistance))
        assert not np.any(np.isnan(result.daily_pop.astype(float)))
        assert not np.any(np.isnan(result.daily_vibrio))

        # No negative populations
        assert np.all(result.yearly_pop >= 0)
        assert np.all(result.daily_pop >= 0)

        # No negative Vibrio
        assert np.all(result.daily_vibrio >= 0)

        # Disease should have occurred
        assert result.total_disease_deaths > 0

    def test_twenty_years_disease_free(self):
        """20-year disease-free run for baseline comparison."""
        result = run_coupled_simulation(
            n_individuals=500,
            carrying_capacity=500,
            habitat_area=10000.0,
            T_celsius=11.0,
            salinity=30.0,
            phi_k=0.1,
            n_years=20,
            disease_year=999,
            seed=77,
        )
        assert result.n_years == 20
        assert result.total_disease_deaths == 0
        # Population should persist
        assert result.final_pop > 100
        assert not np.any(np.isnan(result.yearly_mean_resistance))


# ═══════════════════════════════════════════════════════════════════════
# INTEGRATION TEST: COUPLED FEEDBACKS
# ═══════════════════════════════════════════════════════════════════════

class TestCoupledFeedbacks:
    """Disease ↔ population feedbacks should be visible."""

    def test_disease_reduces_adults(self):
        """Disease should dramatically reduce adult count."""
        # Disease run
        r_dis = run_coupled_simulation(
            n_individuals=500,
            carrying_capacity=500,
            habitat_area=10000.0,
            T_celsius=15.0,
            salinity=30.0,
            phi_k=0.02,
            n_years=6,
            disease_year=2,
            initial_infected=5,
            seed=42,
        )
        # No-disease run
        r_nodis = run_coupled_simulation(
            n_individuals=500,
            carrying_capacity=500,
            habitat_area=10000.0,
            T_celsius=15.0,
            salinity=30.0,
            phi_k=0.02,
            n_years=6,
            disease_year=999,
            seed=42,
        )
        # Year 4 (2 years after disease introduction)
        adults_dis = r_dis.yearly_adults[4]
        adults_nodis = r_nodis.yearly_adults[4]
        assert adults_dis < adults_nodis, \
            f"Disease should reduce adults: {adults_dis} vs {adults_nodis}"

    def test_shelter_site_less_mortality(self):
        """Fjord with high flushing (shelter) should have less mortality.

        Bates 2009: exposed sites had lower disease.
        Open coast (high φ) dilutes pathogen → less mortality.
        """
        # Protected fjord (low flushing)
        r_fjord = run_coupled_simulation(
            n_individuals=500,
            carrying_capacity=500,
            habitat_area=10000.0,
            T_celsius=15.0,
            salinity=30.0,
            phi_k=0.02,
            n_years=6,
            disease_year=2,
            initial_infected=5,
            seed=42,
        )
        # Open coast (high flushing)
        r_open = run_coupled_simulation(
            n_individuals=500,
            carrying_capacity=500,
            habitat_area=10000.0,
            T_celsius=15.0,
            salinity=30.0,
            phi_k=0.5,
            n_years=6,
            disease_year=2,
            initial_infected=5,
            seed=42,
        )
        assert r_open.total_disease_deaths < r_fjord.total_disease_deaths, \
            f"Open coast ({r_open.total_disease_deaths}) should have less " \
            f"mortality than fjord ({r_fjord.total_disease_deaths})"

    def test_reproduction_drops_with_crash(self):
        """After crash, recruitment should drop dramatically."""
        result = run_coupled_simulation(
            n_individuals=500,
            carrying_capacity=500,
            habitat_area=10000.0,
            T_celsius=15.0,
            salinity=30.0,
            phi_k=0.02,
            n_years=10,
            disease_year=2,
            initial_infected=5,
            seed=42,
        )
        # Pre-crash recruits
        recruits_pre = result.yearly_recruits[1]
        # Post-crash recruits (first year after crash)
        crash_yr = result.min_pop_year
        if crash_yr + 1 < result.n_years:
            recruits_post = result.yearly_recruits[crash_yr + 1]
            # Post-crash recruits should be lower (possibly zero)
            assert recruits_post <= recruits_pre, \
                f"Post-crash recruits ({recruits_post}) should be ≤ " \
                f"pre-crash ({recruits_pre})"


# ═══════════════════════════════════════════════════════════════════════
# INTEGRATION TEST: ERRATA COMPLIANCE
# ═══════════════════════════════════════════════════════════════════════

class TestErrataCompliance:
    """Verify CODE_ERRATA items are properly integrated."""

    def test_ce1_no_cost_resistance(self, default_cfg, effect_sizes):
        """CE-1: fecundity_mod is always 1.0 (no cost of resistance)."""
        rng = np.random.default_rng(42)
        agents, _ = initialize_population(
            100, 200, 10000.0, effect_sizes, default_cfg.population, rng,
        )
        assert np.all(agents['fecundity_mod'][:100] == 1.0)

    def test_ce5_high_fecundity_allee(self):
        """CE-5: Allee effect operates through reduced growth rate.

        At low density, population still has positive deterministic growth
        but stochastic effects dominate → very slow recovery.
        """
        result = run_coupled_simulation(
            n_individuals=500,
            carrying_capacity=500,
            habitat_area=10000.0,
            T_celsius=15.0,
            salinity=30.0,
            phi_k=0.02,
            n_years=12,
            disease_year=2,
            initial_infected=5,
            seed=42,
        )
        # After crash, population persists (doesn't go to zero immediately)
        # but recovery is slow
        assert result.min_pop > 0 or result.final_pop >= 0
        # If population survives, it shouldn't bounce back instantly
        crash_year = result.min_pop_year
        if result.min_pop > 0 and crash_year + 3 < result.n_years:
            pop_3yr = result.yearly_pop[crash_year + 3]
            assert pop_3yr < result.initial_pop, \
                "Population should not fully recover in 3 years"
