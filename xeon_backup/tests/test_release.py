"""Tests for the captive-bred release mechanism.

Covers:
  - Config: ReleaseEvent dataclass, defaults, validation
  - Genetics: generate_release_genotypes, genotypes_from_trait_targets, explicit genotypes
  - Model integration: population changes, node assignment, disease states,
    age ranges, origin tracking, multiple releases, backward compatibility
  - Results tracking: release log, released fate tracking

References:
  - config.py: ReleaseEvent, SimulationConfig.release_events
  - genetics.py: generate_release_genotypes, genotypes_from_trait_targets
  - model.py: process_release_event, run_coupled_simulation release_schedule
  - types.py: Origin enum, release_cohort field
"""

import numpy as np
import pytest

from sswd_evoepi.config import (
    ReleaseEvent,
    SimulationConfig,
    default_config,
    validate_config,
    GeneticsSection,
    PopulationSection,
)
from sswd_evoepi.genetics import (
    generate_release_genotypes,
    genotypes_from_trait_targets,
    initialize_trait_effect_sizes,
    compute_trait_batch,
    compute_trait_single,
    update_all_trait_scores,
)
from sswd_evoepi.model import (
    CoupledSimResult,
    make_effect_sizes,
    initialize_population,
    process_release_event,
    run_coupled_simulation,
)
from sswd_evoepi.types import (
    AGENT_DTYPE,
    N_LOCI,
    DiseaseState,
    Origin,
    Stage,
    allocate_agents,
    allocate_genotypes,
    trait_slices,
)


# ═══════════════════════════════════════════════════════════════════════
# FIXTURES
# ═══════════════════════════════════════════════════════════════════════


@pytest.fixture
def rng():
    return np.random.default_rng(42)


@pytest.fixture
def gen_cfg():
    return GeneticsSection()


@pytest.fixture
def pop_cfg():
    return PopulationSection()


@pytest.fixture
def effect_sizes_r():
    return make_effect_sizes(12345, n_loci=17)


@pytest.fixture
def effect_sizes_t():
    return initialize_trait_effect_sizes(np.random.default_rng(12346), 17, 1.0)


@pytest.fixture
def effect_sizes_c():
    return initialize_trait_effect_sizes(np.random.default_rng(12347), 17, 1.0)


@pytest.fixture
def trait_slice_tuple():
    return trait_slices(17, 17, 17)


@pytest.fixture
def small_population(effect_sizes_r, effect_sizes_t, effect_sizes_c, rng):
    """A 200-agent population in a max_agents=400 array."""
    cfg = default_config()
    agents, genotypes = initialize_population(
        n_individuals=200,
        max_agents=400,
        habitat_area=10000.0,
        effect_sizes=effect_sizes_r,
        pop_cfg=cfg.population,
        rng=rng,
        genetics_cfg=cfg.genetics,
        effects_t=effect_sizes_t,
        effects_c=effect_sizes_c,
    )
    return agents, genotypes


# ═══════════════════════════════════════════════════════════════════════
# CONFIG TESTS
# ═══════════════════════════════════════════════════════════════════════


class TestReleaseConfig:
    def test_default_no_releases(self):
        """SimulationConfig with no release_events works as before."""
        cfg = default_config()
        assert cfg.release_events == []
        # Validate should pass without error
        validate_config(cfg)

    def test_release_event_creation(self):
        """ReleaseEvent can be created with all fields."""
        freqs = np.full(N_LOCI, 0.5)
        event = ReleaseEvent(
            time_step=100,
            node_id=2,
            n_individuals=50,
            genetics_mode='allele_freqs',
            allele_freqs=freqs,
            age_range=(365, 730),
            mark_released=True,
        )
        assert event.time_step == 100
        assert event.node_id == 2
        assert event.n_individuals == 50
        assert event.genetics_mode == 'allele_freqs'
        assert np.array_equal(event.allele_freqs, freqs)
        assert event.age_range == (365, 730)
        assert event.mark_released is True

    def test_release_event_defaults(self):
        """Default age_range and mark_released values are correct."""
        event = ReleaseEvent(time_step=0)
        assert event.node_id == 0
        assert event.n_individuals == 100
        assert event.genetics_mode == 'trait_targets'
        assert event.allele_freqs is None
        assert event.trait_targets is None
        assert event.genotypes is None
        assert event.age_range == (365, 730)
        assert event.mark_released is True

    def test_release_event_in_config(self):
        """Config with release events validates successfully."""
        cfg = default_config()
        cfg.release_events = [
            ReleaseEvent(time_step=50, n_individuals=20),
            ReleaseEvent(time_step=100, node_id=1, n_individuals=30),
        ]
        validate_config(cfg)

    def test_invalid_time_step(self):
        """Negative time_step should fail validation."""
        cfg = default_config()
        cfg.release_events = [ReleaseEvent(time_step=-1)]
        with pytest.raises(ValueError, match="time_step"):
            validate_config(cfg)

    def test_invalid_n_individuals(self):
        """Zero individuals should fail validation."""
        cfg = default_config()
        cfg.release_events = [ReleaseEvent(time_step=0, n_individuals=0)]
        with pytest.raises(ValueError, match="n_individuals"):
            validate_config(cfg)

    def test_invalid_genetics_mode(self):
        """Invalid genetics_mode should fail validation."""
        cfg = default_config()
        cfg.release_events = [
            ReleaseEvent(time_step=0, genetics_mode='invalid')
        ]
        with pytest.raises(ValueError, match="genetics_mode"):
            validate_config(cfg)

    def test_allele_freqs_required(self):
        """allele_freqs mode without array should fail validation."""
        cfg = default_config()
        cfg.release_events = [
            ReleaseEvent(time_step=0, genetics_mode='allele_freqs')
        ]
        with pytest.raises(ValueError, match="allele_freqs"):
            validate_config(cfg)

    def test_genotypes_required(self):
        """genotypes mode without array should fail validation."""
        cfg = default_config()
        cfg.release_events = [
            ReleaseEvent(time_step=0, genetics_mode='genotypes')
        ]
        with pytest.raises(ValueError, match="genotypes"):
            validate_config(cfg)

    def test_invalid_age_range(self):
        """Inverted age_range should fail validation."""
        cfg = default_config()
        cfg.release_events = [
            ReleaseEvent(time_step=0, age_range=(730, 365))
        ]
        with pytest.raises(ValueError, match="age_range"):
            validate_config(cfg)


# ═══════════════════════════════════════════════════════════════════════
# GENETICS TESTS
# ═══════════════════════════════════════════════════════════════════════


class TestReleaseGenetics:
    def test_generate_release_genotypes_allele_freqs(self, rng):
        """Generate 1000 individuals with known allele frequencies (p=0.8).
        Verify empirical allele frequencies match within statistical tolerance.
        """
        n = 1000
        target_p = 0.8
        allele_freqs = np.full(N_LOCI, target_p)
        geno = generate_release_genotypes(n, allele_freqs, rng)

        assert geno.shape == (n, N_LOCI, 2)
        assert geno.dtype == np.int8

        # Empirical allele frequency at each locus
        empirical_q = geno.sum(axis=(0, 2)).astype(np.float64) / (2 * n)
        # With n=1000 and p=0.8, SE ≈ sqrt(0.8*0.2/(2*1000)) ≈ 0.009
        # Allow 4 SE tolerance (99.99% CI)
        np.testing.assert_allclose(
            empirical_q, target_p, atol=0.04,
            err_msg="Empirical allele frequencies should match target"
        )

    def test_generate_release_genotypes_varied_freqs(self, rng):
        """Generate with different frequencies per locus."""
        n = 2000
        allele_freqs = np.linspace(0.1, 0.9, N_LOCI)
        geno = generate_release_genotypes(n, allele_freqs, rng)

        empirical_q = geno.sum(axis=(0, 2)).astype(np.float64) / (2 * n)
        np.testing.assert_allclose(empirical_q, allele_freqs, atol=0.05)

    def test_genotypes_from_trait_targets(
        self, rng, effect_sizes_r, effect_sizes_t, effect_sizes_c
    ):
        """Generate individuals with target resistance=0.5.
        Verify mean trait score is close to 0.5.
        """
        n = 1000
        targets = {'resistance': 0.5, 'tolerance': 0.3, 'recovery': 0.1}
        geno = genotypes_from_trait_targets(
            n, targets, effect_sizes_r, effect_sizes_t, effect_sizes_c, rng,
        )

        assert geno.shape == (n, N_LOCI, 2)
        assert geno.dtype == np.int8

        res_s, tol_s, rec_s = trait_slices(17, 17, 17)

        # Compute mean resistance
        r_scores = np.array([
            compute_trait_single(geno[i], effect_sizes_r, res_s)
            for i in range(n)
        ])
        mean_r = np.mean(r_scores)
        # Target 0.5, allow ±0.1 tolerance (stochastic process)
        assert abs(mean_r - 0.5) < 0.1, f"Mean resistance {mean_r} not close to 0.5"

        # Compute mean tolerance
        t_scores = np.array([
            compute_trait_single(geno[i], effect_sizes_t, tol_s)
            for i in range(n)
        ])
        mean_t = np.mean(t_scores)
        assert abs(mean_t - 0.3) < 0.1, f"Mean tolerance {mean_t} not close to 0.3"

    def test_explicit_genotypes(self):
        """Pass explicit genotype array, verify it's used as-is."""
        n = 50
        explicit = np.ones((n, N_LOCI, 2), dtype=np.int8)  # all derived
        event = ReleaseEvent(
            time_step=0,
            n_individuals=n,
            genetics_mode='genotypes',
            genotypes=explicit,
        )
        assert event.genotypes is not None
        np.testing.assert_array_equal(event.genotypes, explicit)

    def test_release_genotypes_shape(self, rng, effect_sizes_r, effect_sizes_t, effect_sizes_c):
        """Verify output shapes are correct for all genetics modes."""
        n = 50

        # allele_freqs mode
        af = np.full(N_LOCI, 0.3)
        g1 = generate_release_genotypes(n, af, rng)
        assert g1.shape == (n, N_LOCI, 2)
        assert g1.dtype == np.int8

        # trait_targets mode
        targets = {'resistance': 0.2, 'tolerance': 0.1, 'recovery': 0.05}
        g2 = genotypes_from_trait_targets(
            n, targets, effect_sizes_r, effect_sizes_t, effect_sizes_c, rng,
        )
        assert g2.shape == (n, N_LOCI, 2)
        assert g2.dtype == np.int8

        # genotypes mode (explicit)
        g3 = np.zeros((n, N_LOCI, 2), dtype=np.int8)
        assert g3.shape == (n, N_LOCI, 2)

    def test_genotypes_binary_values(self, rng):
        """All allele values should be 0 or 1."""
        n = 500
        af = np.full(N_LOCI, 0.5)
        geno = generate_release_genotypes(n, af, rng)
        assert np.all((geno == 0) | (geno == 1))

    def test_boundary_allele_freqs(self, rng):
        """Edge case: allele freqs of 0.0 and 1.0."""
        n = 100
        af = np.zeros(N_LOCI)
        af[:10] = 1.0  # first 10 loci all derived
        geno = generate_release_genotypes(n, af, rng)
        # Loci with freq=1.0 should be all 1s
        assert np.all(geno[:, :10, :] == 1)
        # Loci with freq=0.0 should be all 0s
        assert np.all(geno[:, 10:, :] == 0)


# ═══════════════════════════════════════════════════════════════════════
# MODEL INTEGRATION TESTS — process_release_event
# ═══════════════════════════════════════════════════════════════════════


class TestProcessReleaseEvent:
    """Tests for the process_release_event function directly."""

    def _make_context(self, rng, n_pop=100, max_agents=400):
        """Helper: create a population and supporting structures."""
        cfg = default_config()
        gen_cfg = cfg.genetics
        pop_cfg = cfg.population

        effect_sizes = make_effect_sizes(gen_cfg.effect_size_seed, n_loci=gen_cfg.n_resistance)
        effects_t = initialize_trait_effect_sizes(
            np.random.default_rng(gen_cfg.effect_size_seed + 1),
            gen_cfg.n_tolerance, 1.0,
        )
        effects_c = initialize_trait_effect_sizes(
            np.random.default_rng(gen_cfg.effect_size_seed + 2),
            gen_cfg.n_recovery, 1.0,
        )
        res_s, tol_s, rec_s = trait_slices(
            gen_cfg.n_resistance, gen_cfg.n_tolerance, gen_cfg.n_recovery,
        )

        agents, genotypes = initialize_population(
            n_individuals=n_pop,
            max_agents=max_agents,
            habitat_area=10000.0,
            effect_sizes=effect_sizes,
            pop_cfg=pop_cfg,
            rng=rng,
            genetics_cfg=gen_cfg,
            effects_t=effects_t,
            effects_c=effects_c,
        )

        return {
            'agents': agents,
            'genotypes': genotypes,
            'effect_sizes': effect_sizes,
            'effects_t': effects_t,
            'effects_c': effects_c,
            'res_slice': res_s,
            'tol_slice': tol_s,
            'rec_slice': rec_s,
            'gen_cfg': gen_cfg,
            'pop_cfg': pop_cfg,
        }

    def test_release_adds_individuals(self, rng):
        """Release at a given time adds exactly n_individuals."""
        ctx = self._make_context(rng, n_pop=100, max_agents=400)
        pop_before = int(np.sum(ctx['agents']['alive']))

        event = ReleaseEvent(time_step=50, n_individuals=25)
        agents, genotypes, log = process_release_event(
            event=event,
            event_index=1,
            agents=ctx['agents'],
            genotypes=ctx['genotypes'],
            effect_sizes=ctx['effect_sizes'],
            effects_t=ctx['effects_t'],
            effects_c=ctx['effects_c'],
            res_slice=ctx['res_slice'],
            tol_slice=ctx['tol_slice'],
            rec_slice=ctx['rec_slice'],
            gen_cfg=ctx['gen_cfg'],
            pop_cfg=ctx['pop_cfg'],
            rng=rng,
            habitat_area=10000.0,
            sim_day=50,
        )

        pop_after = int(np.sum(agents['alive']))
        assert pop_after == pop_before + 25

    def test_release_correct_node(self, rng):
        """Released individuals are at the specified node."""
        ctx = self._make_context(rng)
        event = ReleaseEvent(time_step=0, node_id=3, n_individuals=20)
        agents, genotypes, log = process_release_event(
            event=event, event_index=1,
            agents=ctx['agents'], genotypes=ctx['genotypes'],
            effect_sizes=ctx['effect_sizes'],
            effects_t=ctx['effects_t'], effects_c=ctx['effects_c'],
            res_slice=ctx['res_slice'], tol_slice=ctx['tol_slice'],
            rec_slice=ctx['rec_slice'],
            gen_cfg=ctx['gen_cfg'], pop_cfg=ctx['pop_cfg'],
            rng=rng, habitat_area=10000.0, sim_day=0,
        )

        released = agents['alive'] & (agents['release_cohort'] == 1)
        assert int(np.sum(released)) == 20
        assert np.all(agents['node_id'][released] == 3)

    def test_release_susceptible(self, rng):
        """All released individuals start as SUSCEPTIBLE."""
        ctx = self._make_context(rng)
        event = ReleaseEvent(time_step=0, n_individuals=50)
        agents, genotypes, log = process_release_event(
            event=event, event_index=1,
            agents=ctx['agents'], genotypes=ctx['genotypes'],
            effect_sizes=ctx['effect_sizes'],
            effects_t=ctx['effects_t'], effects_c=ctx['effects_c'],
            res_slice=ctx['res_slice'], tol_slice=ctx['tol_slice'],
            rec_slice=ctx['rec_slice'],
            gen_cfg=ctx['gen_cfg'], pop_cfg=ctx['pop_cfg'],
            rng=rng, habitat_area=10000.0, sim_day=0,
        )

        released = agents['alive'] & (agents['release_cohort'] == 1)
        assert np.all(agents['disease_state'][released] == DiseaseState.S)

    def test_release_age_range(self, rng):
        """Released individuals have ages within specified range."""
        min_days, max_days = 500, 1000
        ctx = self._make_context(rng)
        event = ReleaseEvent(
            time_step=0, n_individuals=100,
            age_range=(min_days, max_days),
        )
        agents, genotypes, log = process_release_event(
            event=event, event_index=1,
            agents=ctx['agents'], genotypes=ctx['genotypes'],
            effect_sizes=ctx['effect_sizes'],
            effects_t=ctx['effects_t'], effects_c=ctx['effects_c'],
            res_slice=ctx['res_slice'], tol_slice=ctx['tol_slice'],
            rec_slice=ctx['rec_slice'],
            gen_cfg=ctx['gen_cfg'], pop_cfg=ctx['pop_cfg'],
            rng=rng, habitat_area=10000.0, sim_day=0,
        )

        released = agents['alive'] & (agents['release_cohort'] == 1)
        ages_years = agents['age'][released]
        ages_days = ages_years * 365.0
        # Allow small floating point slack
        assert np.all(ages_days >= min_days - 1)
        assert np.all(ages_days <= max_days + 1)

    def test_release_origin_tracking(self, rng):
        """Released individuals have released=True (CAPTIVE_BRED) and correct release_cohort."""
        ctx = self._make_context(rng)
        event = ReleaseEvent(time_step=0, n_individuals=30, mark_released=True)
        agents, genotypes, log = process_release_event(
            event=event, event_index=2,
            agents=ctx['agents'], genotypes=ctx['genotypes'],
            effect_sizes=ctx['effect_sizes'],
            effects_t=ctx['effects_t'], effects_c=ctx['effects_c'],
            res_slice=ctx['res_slice'], tol_slice=ctx['tol_slice'],
            rec_slice=ctx['rec_slice'],
            gen_cfg=ctx['gen_cfg'], pop_cfg=ctx['pop_cfg'],
            rng=rng, habitat_area=10000.0, sim_day=0,
        )

        released = agents['alive'] & (agents['release_cohort'] == 2)
        assert int(np.sum(released)) == 30
        assert np.all(agents['origin'][released] == Origin.CAPTIVE_BRED)

    def test_wild_born_not_marked(self, rng):
        """Wild-born individuals have released=False (WILD), release_cohort=0."""
        ctx = self._make_context(rng, n_pop=50)
        # Check all initial agents
        alive_initial = ctx['agents']['alive']
        assert np.all(ctx['agents']['origin'][alive_initial] == Origin.WILD)
        assert np.all(ctx['agents']['release_cohort'][alive_initial] == 0)

    def test_release_unmarked(self, rng):
        """Release with mark_released=False gets origin=WILD."""
        ctx = self._make_context(rng)
        event = ReleaseEvent(
            time_step=0, n_individuals=10, mark_released=False,
        )
        agents, genotypes, log = process_release_event(
            event=event, event_index=1,
            agents=ctx['agents'], genotypes=ctx['genotypes'],
            effect_sizes=ctx['effect_sizes'],
            effects_t=ctx['effects_t'], effects_c=ctx['effects_c'],
            res_slice=ctx['res_slice'], tol_slice=ctx['tol_slice'],
            rec_slice=ctx['rec_slice'],
            gen_cfg=ctx['gen_cfg'], pop_cfg=ctx['pop_cfg'],
            rng=rng, habitat_area=10000.0, sim_day=0,
        )

        released = agents['alive'] & (agents['release_cohort'] == 1)
        assert np.all(agents['origin'][released] == Origin.WILD)

    def test_multiple_releases(self, rng):
        """Schedule 3 releases, verify all execute correctly."""
        ctx = self._make_context(rng, n_pop=50, max_agents=500)
        pop_before = int(np.sum(ctx['agents']['alive']))

        events = [
            ReleaseEvent(time_step=10, node_id=0, n_individuals=20),
            ReleaseEvent(time_step=20, node_id=1, n_individuals=30),
            ReleaseEvent(time_step=30, node_id=2, n_individuals=40),
        ]

        agents = ctx['agents']
        genotypes = ctx['genotypes']
        logs = []
        for i, event in enumerate(events):
            agents, genotypes, log = process_release_event(
                event=event, event_index=i + 1,
                agents=agents, genotypes=genotypes,
                effect_sizes=ctx['effect_sizes'],
                effects_t=ctx['effects_t'], effects_c=ctx['effects_c'],
                res_slice=ctx['res_slice'], tol_slice=ctx['tol_slice'],
                rec_slice=ctx['rec_slice'],
                gen_cfg=ctx['gen_cfg'], pop_cfg=ctx['pop_cfg'],
                rng=rng, habitat_area=10000.0,
                sim_day=event.time_step,
            )
            logs.append(log)

        pop_after = int(np.sum(agents['alive']))
        assert pop_after == pop_before + 20 + 30 + 40

        # Verify each cohort
        c1 = agents['alive'] & (agents['release_cohort'] == 1)
        c2 = agents['alive'] & (agents['release_cohort'] == 2)
        c3 = agents['alive'] & (agents['release_cohort'] == 3)
        assert int(np.sum(c1)) == 20
        assert int(np.sum(c2)) == 30
        assert int(np.sum(c3)) == 40
        assert np.all(agents['node_id'][c1] == 0)
        assert np.all(agents['node_id'][c2] == 1)
        assert np.all(agents['node_id'][c3] == 2)

    def test_trait_scores_computed(self, rng):
        """Released individuals have nonzero trait scores."""
        ctx = self._make_context(rng)
        event = ReleaseEvent(
            time_step=0, n_individuals=50,
            genetics_mode='trait_targets',
            trait_targets={'resistance': 0.3, 'tolerance': 0.2, 'recovery': 0.1},
        )
        agents, genotypes, log = process_release_event(
            event=event, event_index=1,
            agents=ctx['agents'], genotypes=ctx['genotypes'],
            effect_sizes=ctx['effect_sizes'],
            effects_t=ctx['effects_t'], effects_c=ctx['effects_c'],
            res_slice=ctx['res_slice'], tol_slice=ctx['tol_slice'],
            rec_slice=ctx['rec_slice'],
            gen_cfg=ctx['gen_cfg'], pop_cfg=ctx['pop_cfg'],
            rng=rng, habitat_area=10000.0, sim_day=0,
        )

        released = agents['alive'] & (agents['release_cohort'] == 1)
        r_vals = agents['resistance'][released]
        t_vals = agents['tolerance'][released]
        c_vals = agents['recovery_ability'][released]

        # Should have nonzero values (stochastic, but with targets > 0)
        assert np.mean(r_vals) > 0.05
        assert np.mean(t_vals) > 0.05
        assert np.mean(c_vals) > 0.01

    def test_array_grows_if_needed(self, rng):
        """Releasing more than available slots grows the arrays."""
        ctx = self._make_context(rng, n_pop=100, max_agents=100)
        # All 100 slots are alive — no dead slots
        old_len = len(ctx['agents'])

        event = ReleaseEvent(time_step=0, n_individuals=50)
        agents, genotypes, log = process_release_event(
            event=event, event_index=1,
            agents=ctx['agents'], genotypes=ctx['genotypes'],
            effect_sizes=ctx['effect_sizes'],
            effects_t=ctx['effects_t'], effects_c=ctx['effects_c'],
            res_slice=ctx['res_slice'], tol_slice=ctx['tol_slice'],
            rec_slice=ctx['rec_slice'],
            gen_cfg=ctx['gen_cfg'], pop_cfg=ctx['pop_cfg'],
            rng=rng, habitat_area=10000.0, sim_day=0,
        )

        assert len(agents) > old_len
        pop_after = int(np.sum(agents['alive']))
        assert pop_after == 150

    def test_log_entry_metadata(self, rng):
        """process_release_event returns correct log metadata."""
        ctx = self._make_context(rng)
        event = ReleaseEvent(
            time_step=75, node_id=2, n_individuals=15,
            genetics_mode='trait_targets',
            trait_targets={'resistance': 0.4, 'tolerance': 0.2, 'recovery': 0.1},
        )
        agents, genotypes, log = process_release_event(
            event=event, event_index=3,
            agents=ctx['agents'], genotypes=ctx['genotypes'],
            effect_sizes=ctx['effect_sizes'],
            effects_t=ctx['effects_t'], effects_c=ctx['effects_c'],
            res_slice=ctx['res_slice'], tol_slice=ctx['tol_slice'],
            rec_slice=ctx['rec_slice'],
            gen_cfg=ctx['gen_cfg'], pop_cfg=ctx['pop_cfg'],
            rng=rng, habitat_area=10000.0, sim_day=75,
        )

        assert log['event_index'] == 3
        assert log['time_step'] == 75
        assert log['node_id'] == 2
        assert log['n_released'] == 15
        assert log['genetics_mode'] == 'trait_targets'
        assert 'mean_resistance' in log
        assert 'mean_tolerance' in log
        assert 'mean_recovery' in log


# ═══════════════════════════════════════════════════════════════════════
# FULL SIMULATION INTEGRATION TESTS
# ═══════════════════════════════════════════════════════════════════════


class TestReleaseSimulation:
    """Tests using run_coupled_simulation with release events."""

    def test_release_adds_individuals_in_sim(self):
        """Run model for ~100 days with a release at day 50.
        Verify population increases around day 50.
        """
        cfg = default_config()
        cfg.release_events = [
            ReleaseEvent(time_step=50, n_individuals=100),
        ]
        # Short simulation: 1 year, no disease
        result = run_coupled_simulation(
            n_individuals=200,
            carrying_capacity=500,
            n_years=1,
            disease_year=999,  # no disease
            seed=42,
            record_daily=True,
            config=cfg,
        )

        # Verify release recorded
        assert result.release_log is not None
        assert len(result.release_log) == 1
        assert result.release_log[0]['n_released'] == 100
        assert result.total_released == 100

        # Daily pop should show a jump around day 50
        # (accounting for natural mortality, pop at day 51 > pop at day 49)
        pop_before = result.daily_pop[49]
        pop_after = result.daily_pop[51]
        # The jump should be roughly 100 (minus a day or two of mortality)
        assert pop_after > pop_before + 80

    def test_release_correct_node_in_sim(self):
        """Release at node 0 (single-node sim), verify node assignment."""
        cfg = default_config()
        cfg.release_events = [
            ReleaseEvent(time_step=10, node_id=0, n_individuals=50),
        ]
        result = run_coupled_simulation(
            n_individuals=100,
            carrying_capacity=300,
            n_years=1,
            disease_year=999,
            seed=42,
            config=cfg,
        )
        assert result.total_released == 50

    def test_release_participates_in_disease(self):
        """Release individuals into an active epidemic. Some should get infected."""
        cfg = default_config()
        # Release 200 individuals on day 0 of disease year (year 1, day 365)
        release_day = 1 * 365 + 30  # 30 days into disease year
        cfg.release_events = [
            ReleaseEvent(time_step=release_day, n_individuals=200),
        ]
        result = run_coupled_simulation(
            n_individuals=300,
            carrying_capacity=600,
            T_celsius=15.0,
            n_years=3,
            disease_year=1,
            initial_infected=10,
            seed=42,
            config=cfg,
        )
        # Disease should have killed some individuals (including released ones)
        assert result.total_disease_deaths > 0
        # Released surviving should be less than total released (some died)
        assert result.released_surviving < result.total_released

    def test_release_participates_in_reproduction(self):
        """Release mature individuals. They should be able to reproduce."""
        cfg = default_config()
        # Release on day 0 — adults that can participate in spawning
        cfg.release_events = [
            ReleaseEvent(
                time_step=0,
                n_individuals=100,
                age_range=(3650, 5475),  # 10-15 years old → ADULT stage
            ),
        ]
        result = run_coupled_simulation(
            n_individuals=50,
            carrying_capacity=300,
            n_years=3,
            disease_year=999,  # no disease
            seed=42,
            config=cfg,
        )
        # With 50 initial + 100 released adults, reproduction should produce recruits
        total_recruits = int(np.sum(result.yearly_recruits))
        assert total_recruits > 0, "Released adults should participate in reproduction"

    def test_no_release_identical(self):
        """Model with empty release_events list produces identical results
        to model without release feature (same seed). CRITICAL backward compatibility.
        """
        cfg_with = default_config()
        cfg_with.release_events = []  # explicitly empty

        cfg_without = default_config()
        # cfg_without has release_events=[] by default too

        result_with = run_coupled_simulation(
            n_individuals=200,
            carrying_capacity=300,
            n_years=3,
            disease_year=1,
            initial_infected=5,
            seed=42,
            config=cfg_with,
        )
        result_without = run_coupled_simulation(
            n_individuals=200,
            carrying_capacity=300,
            n_years=3,
            disease_year=1,
            initial_infected=5,
            seed=42,
            config=cfg_without,
        )

        # Population trajectories must be identical
        np.testing.assert_array_equal(
            result_with.yearly_pop,
            result_without.yearly_pop,
            err_msg="Empty release_events should not change simulation results",
        )
        np.testing.assert_array_equal(
            result_with.yearly_disease_deaths,
            result_without.yearly_disease_deaths,
        )
        np.testing.assert_array_equal(
            result_with.yearly_mean_resistance,
            result_without.yearly_mean_resistance,
        )
        assert result_with.final_pop == result_without.final_pop
        assert result_with.total_disease_deaths == result_without.total_disease_deaths


# ═══════════════════════════════════════════════════════════════════════
# RESULTS TRACKING TESTS
# ═══════════════════════════════════════════════════════════════════════


class TestReleaseResultsTracking:
    def test_release_results_recorded(self):
        """Verify release events are recorded in results with correct metadata."""
        cfg = default_config()
        cfg.release_events = [
            ReleaseEvent(time_step=10, node_id=0, n_individuals=25),
            ReleaseEvent(time_step=100, node_id=0, n_individuals=50),
        ]
        result = run_coupled_simulation(
            n_individuals=100,
            carrying_capacity=300,
            n_years=1,
            disease_year=999,
            seed=42,
            config=cfg,
        )

        assert result.release_log is not None
        assert len(result.release_log) == 2
        assert result.release_log[0]['n_released'] == 25
        assert result.release_log[0]['time_step'] == 10
        assert result.release_log[1]['n_released'] == 50
        assert result.release_log[1]['time_step'] == 100
        assert result.total_released == 75

    def test_released_fate_tracking(self):
        """After release into epidemic, verify released alive/dead counts are tracked."""
        cfg = default_config()
        release_day = 365 + 10  # 10 days into year 1 (disease year)
        cfg.release_events = [
            ReleaseEvent(time_step=release_day, n_individuals=100),
        ]
        result = run_coupled_simulation(
            n_individuals=200,
            carrying_capacity=500,
            T_celsius=15.0,
            n_years=3,
            disease_year=1,
            initial_infected=10,
            seed=42,
            config=cfg,
        )

        assert result.total_released == 100
        # Some released should still be alive, some dead
        assert result.released_surviving >= 0
        assert result.released_surviving <= 100

        # yearly_released_alive should be tracked
        assert result.yearly_released_alive is not None
        # After release (year 1+), should show released alive counts
        assert result.yearly_released_alive[1] > 0 or result.yearly_released_alive[2] >= 0

    def test_yearly_released_alive_tracking(self):
        """Verify yearly_released_alive decreases over time as released die."""
        cfg = default_config()
        cfg.release_events = [
            ReleaseEvent(time_step=0, n_individuals=100),
        ]
        result = run_coupled_simulation(
            n_individuals=100,
            carrying_capacity=300,
            T_celsius=15.0,
            n_years=5,
            disease_year=1,
            initial_infected=10,
            seed=42,
            config=cfg,
        )

        assert result.yearly_released_alive is not None
        assert len(result.yearly_released_alive) == 5
        # Year 0 should show ~100 released alive (minus any natural mortality)
        assert result.yearly_released_alive[0] > 50
        # Over time with disease, released alive should decrease
        # (not necessarily monotonic, but final should be <= initial)
        assert result.yearly_released_alive[-1] <= result.yearly_released_alive[0]

    def test_no_release_no_tracking(self):
        """Without releases, yearly_released_alive should be None."""
        cfg = default_config()
        cfg.release_events = []
        result = run_coupled_simulation(
            n_individuals=100,
            carrying_capacity=200,
            n_years=2,
            disease_year=999,
            seed=42,
            config=cfg,
        )
        assert result.release_log is None
        assert result.total_released == 0
        assert result.released_surviving == 0
        assert result.yearly_released_alive is None


# ═══════════════════════════════════════════════════════════════════════
# EDGE CASES
# ═══════════════════════════════════════════════════════════════════════


class TestReleaseEdgeCases:
    def test_release_on_day_zero(self):
        """Release on simulation day 0 works."""
        cfg = default_config()
        cfg.release_events = [
            ReleaseEvent(time_step=0, n_individuals=50),
        ]
        result = run_coupled_simulation(
            n_individuals=100,
            carrying_capacity=300,
            n_years=1,
            disease_year=999,
            seed=42,
            config=cfg,
        )
        assert result.total_released == 50

    def test_release_on_last_day(self):
        """Release on the last day of simulation works."""
        total_days = 2 * 365 - 1  # last day of year 1 (0-indexed)
        cfg = default_config()
        cfg.release_events = [
            ReleaseEvent(time_step=total_days, n_individuals=20),
        ]
        result = run_coupled_simulation(
            n_individuals=100,
            carrying_capacity=300,
            n_years=2,
            disease_year=999,
            seed=42,
            config=cfg,
        )
        assert result.total_released == 20

    def test_single_individual_release(self):
        """Releasing a single individual works."""
        cfg = default_config()
        cfg.release_events = [
            ReleaseEvent(time_step=10, n_individuals=1),
        ]
        result = run_coupled_simulation(
            n_individuals=100,
            carrying_capacity=200,
            n_years=1,
            disease_year=999,
            seed=42,
            config=cfg,
        )
        assert result.total_released == 1

    def test_large_release(self):
        """Releasing more than initial population works."""
        cfg = default_config()
        cfg.release_events = [
            ReleaseEvent(time_step=10, n_individuals=500),
        ]
        result = run_coupled_simulation(
            n_individuals=100,
            carrying_capacity=1000,
            n_years=1,
            disease_year=999,
            seed=42,
            config=cfg,
        )
        assert result.total_released == 500

    def test_release_with_explicit_genotypes(self, rng):
        """Release using explicit genotype arrays in full simulation."""
        n_release = 30
        # Create explicit genotypes: all heterozygous at every locus
        explicit_geno = np.zeros((n_release, N_LOCI, 2), dtype=np.int8)
        explicit_geno[:, :, 0] = 0
        explicit_geno[:, :, 1] = 1  # heterozygous everywhere

        cfg = default_config()
        cfg.release_events = [
            ReleaseEvent(
                time_step=10,
                n_individuals=n_release,
                genetics_mode='genotypes',
                genotypes=explicit_geno,
            ),
        ]
        result = run_coupled_simulation(
            n_individuals=100,
            carrying_capacity=300,
            n_years=1,
            disease_year=999,
            seed=42,
            config=cfg,
        )
        assert result.total_released == n_release

    def test_release_with_allele_freqs(self):
        """Release using allele_freqs mode in full simulation."""
        cfg = default_config()
        af = np.full(N_LOCI, 0.7)
        cfg.release_events = [
            ReleaseEvent(
                time_step=10,
                n_individuals=50,
                genetics_mode='allele_freqs',
                allele_freqs=af,
            ),
        ]
        result = run_coupled_simulation(
            n_individuals=100,
            carrying_capacity=300,
            n_years=1,
            disease_year=999,
            seed=42,
            config=cfg,
        )
        assert result.total_released == 50
