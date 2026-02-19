"""Tests for continuous settlement data structures and settle_daily_cohorts.

Phase 1: Verifies LarvalCohort spawn_day/sst_at_spawn fields, the standalone
settle_daily_cohorts function (isolation, no modifications to main sim loop),
and PLD temperature dependence matching Hodin 2021 data.
"""

import numpy as np
import pytest

from sswd_evoepi.types import (
    AGENT_DTYPE,
    LarvalCohort,
    N_LOCI,
    DiseaseState,
    Stage,
    allocate_agents,
    allocate_genotypes,
)
from sswd_evoepi.config import PopulationSection
from sswd_evoepi.model import make_effect_sizes, settle_daily_cohorts
from sswd_evoepi.reproduction import pelagic_larval_duration


# ═══════════════════════════════════════════════════════════════════════
# HELPERS
# ═══════════════════════════════════════════════════════════════════════

def _make_population(n_alive: int, n_total: int, rng: np.random.Generator = None):
    """Create a test population with n_alive alive agents out of n_total slots.

    Returns (agents, genotypes, effect_sizes).
    Alive agents are adults at node 0 with disease_state=S.
    """
    if rng is None:
        rng = np.random.default_rng(42)
    agents = allocate_agents(n_total)
    genotypes = allocate_genotypes(n_total)

    for i in range(n_alive):
        agents[i]['alive'] = True
        agents[i]['stage'] = Stage.ADULT
        agents[i]['disease_state'] = DiseaseState.S
        agents[i]['size'] = 600.0
        agents[i]['age'] = 5.0
        agents[i]['sex'] = rng.integers(0, 2)
        agents[i]['node_id'] = 0
        # Random genotypes
        genotypes[i] = rng.integers(0, 2, size=(N_LOCI, 2), dtype=np.int8)

    effect_sizes = make_effect_sizes()
    return agents, genotypes, effect_sizes


def _make_cohort(
    n_competent: int,
    source_node: int = 0,
    spawn_day: int = 0,
    sst_at_spawn: float = 10.5,
    rng: np.random.Generator = None,
) -> LarvalCohort:
    """Create a test LarvalCohort with random genotypes."""
    if rng is None:
        rng = np.random.default_rng(99)
    geno = rng.integers(0, 2, size=(n_competent, N_LOCI, 2), dtype=np.int8)
    parents = np.zeros((n_competent, 2), dtype=np.int32)
    pld = pelagic_larval_duration(sst_at_spawn)
    return LarvalCohort(
        source_node=source_node,
        n_competent=n_competent,
        genotypes=geno,
        parent_pairs=parents,
        pld_days=pld,
        spawn_day=spawn_day,
        sst_at_spawn=sst_at_spawn,
    )


# ═══════════════════════════════════════════════════════════════════════
# DATA STRUCTURE TESTS
# ═══════════════════════════════════════════════════════════════════════


class TestLarvalCohortSpawnDay:
    """Verify spawn_day and sst_at_spawn fields on LarvalCohort."""

    def test_larval_cohort_has_spawn_day(self):
        cohort = LarvalCohort(
            source_node=0,
            n_competent=10,
            genotypes=np.zeros((10, N_LOCI, 2), dtype=np.int8),
            parent_pairs=np.zeros((10, 2), dtype=np.int32),
            pld_days=63.0,
            spawn_day=100,
            sst_at_spawn=12.0,
        )
        assert cohort.spawn_day == 100
        assert cohort.sst_at_spawn == 12.0

    def test_default_spawn_day(self):
        """Default spawn_day=0 and sst_at_spawn=10.5 when not provided."""
        cohort = LarvalCohort(
            source_node=0,
            n_competent=5,
            genotypes=np.zeros((5, N_LOCI, 2), dtype=np.int8),
            parent_pairs=np.zeros((5, 2), dtype=np.int32),
            pld_days=63.0,
        )
        assert cohort.spawn_day == 0
        assert cohort.sst_at_spawn == 10.5


# ═══════════════════════════════════════════════════════════════════════
# SETTLEMENT FUNCTION TESTS
# ═══════════════════════════════════════════════════════════════════════


class TestSettleDailyCohorts:
    """Tests for the standalone settle_daily_cohorts function."""

    def test_settle_empty_cohorts(self):
        """No cohorts → 0 recruits, population unchanged."""
        agents, genotypes, effect_sizes = _make_population(50, 100)
        pop_cfg = PopulationSection()
        rng = np.random.default_rng(42)

        n_before = int(np.sum(agents['alive']))
        result = settle_daily_cohorts(
            cohorts=[],
            agents=agents,
            genotypes=genotypes,
            carrying_capacity=100,
            pop_cfg=pop_cfg,
            rng=rng,
            effect_sizes=effect_sizes,
        )
        assert result == 0
        assert int(np.sum(agents['alive'])) == n_before

    def test_settle_single_cohort(self):
        """One cohort with larvae settles, agents appear in population."""
        agents, genotypes, effect_sizes = _make_population(50, 200)
        pop_cfg = PopulationSection()
        rng = np.random.default_rng(42)

        cohort = _make_cohort(n_competent=500, spawn_day=10, sst_at_spawn=10.5, rng=rng)
        n_before = int(np.sum(agents['alive']))

        n_settled = settle_daily_cohorts(
            cohorts=[cohort],
            agents=agents,
            genotypes=genotypes,
            carrying_capacity=200,
            pop_cfg=pop_cfg,
            rng=np.random.default_rng(42),
            effect_sizes=effect_sizes,
        )

        n_after = int(np.sum(agents['alive']))
        assert n_settled > 0, "Should settle at least some larvae"
        assert n_after == n_before + n_settled
        assert n_after <= 200, "Should not exceed carrying capacity"

    def test_settle_respects_carrying_capacity(self):
        """Settlement stops when population reaches K."""
        # Start near capacity
        K = 100
        agents, genotypes, effect_sizes = _make_population(95, K)
        pop_cfg = PopulationSection()
        rng = np.random.default_rng(42)

        # Try to settle a huge cohort
        cohort = _make_cohort(n_competent=10000, rng=rng)
        n_settled = settle_daily_cohorts(
            cohorts=[cohort],
            agents=agents,
            genotypes=genotypes,
            carrying_capacity=K,
            pop_cfg=pop_cfg,
            rng=np.random.default_rng(42),
            effect_sizes=effect_sizes,
        )

        n_alive = int(np.sum(agents['alive']))
        assert n_alive <= K, f"Population {n_alive} exceeds K={K}"

    def test_settle_at_capacity_returns_zero(self):
        """If already at K, no settlement occurs."""
        K = 100
        agents, genotypes, effect_sizes = _make_population(K, K)
        pop_cfg = PopulationSection()
        rng = np.random.default_rng(42)

        cohort = _make_cohort(n_competent=100, rng=rng)
        n_settled = settle_daily_cohorts(
            cohorts=[cohort],
            agents=agents,
            genotypes=genotypes,
            carrying_capacity=K,
            pop_cfg=pop_cfg,
            rng=np.random.default_rng(42),
            effect_sizes=effect_sizes,
        )
        assert n_settled == 0

    def test_settle_beverton_holt(self):
        """Fewer recruits when close to K (slot-limited by available capacity)."""
        pop_cfg = PopulationSection()
        K = 500
        rng_seed = 42

        # Low density: 50/500 alive → 450 available slots
        agents_low, geno_low, eff = _make_population(50, K)
        cohort_low = _make_cohort(n_competent=2000, rng=np.random.default_rng(99))
        n_low = settle_daily_cohorts(
            cohorts=[cohort_low],
            agents=agents_low,
            genotypes=geno_low,
            carrying_capacity=K,
            pop_cfg=pop_cfg,
            rng=np.random.default_rng(rng_seed),
            effect_sizes=eff,
        )

        # Near capacity: 490/500 alive → only 10 available slots
        agents_high, geno_high, eff2 = _make_population(490, K)
        cohort_high = _make_cohort(n_competent=2000, rng=np.random.default_rng(99))
        n_high = settle_daily_cohorts(
            cohorts=[cohort_high],
            agents=agents_high,
            genotypes=geno_high,
            carrying_capacity=K,
            pop_cfg=pop_cfg,
            rng=np.random.default_rng(rng_seed),
            effect_sizes=eff2,
        )

        # Near capacity: slot limitation should constrain recruitment
        assert n_low > n_high, (
            f"Expected more recruits at low density ({n_low}) than near-capacity ({n_high})"
        )
        assert n_high <= 10, f"Near capacity should settle ≤10 but got {n_high}"

    def test_settlers_are_susceptible(self):
        """All newly settled agents have disease_state == S."""
        agents, genotypes, effect_sizes = _make_population(30, 200)
        pop_cfg = PopulationSection()
        rng = np.random.default_rng(42)

        cohort = _make_cohort(n_competent=500, rng=np.random.default_rng(99))

        # Track which slots are dead before settlement
        dead_before = set(np.where(~agents['alive'])[0])

        n_settled = settle_daily_cohorts(
            cohorts=[cohort],
            agents=agents,
            genotypes=genotypes,
            carrying_capacity=200,
            pop_cfg=pop_cfg,
            rng=rng,
            effect_sizes=effect_sizes,
        )

        assert n_settled > 0
        # Find newly alive agents (were dead, now alive)
        now_alive = set(np.where(agents['alive'])[0])
        newly_settled = now_alive & dead_before

        for idx in newly_settled:
            assert agents[idx]['disease_state'] == DiseaseState.S, (
                f"Settler at slot {idx} has disease_state={agents[idx]['disease_state']}, expected S"
            )
            assert agents[idx]['stage'] == Stage.SETTLER
            assert agents[idx]['age'] == 0.0
            assert agents[idx]['size'] == pytest.approx(0.5)

    def test_settlers_have_genotypes(self):
        """Settled agents have non-trivial genotypes properly copied."""
        agents, genotypes, effect_sizes = _make_population(20, 200)
        pop_cfg = PopulationSection()
        rng = np.random.default_rng(42)

        # Create cohort with known genotypes (all 1s)
        n_comp = 100
        geno_all_ones = np.ones((n_comp, N_LOCI, 2), dtype=np.int8)
        cohort = LarvalCohort(
            source_node=0,
            n_competent=n_comp,
            genotypes=geno_all_ones,
            parent_pairs=np.zeros((n_comp, 2), dtype=np.int32),
            pld_days=63.0,
            spawn_day=0,
            sst_at_spawn=10.5,
        )

        dead_before = set(np.where(~agents['alive'])[0])

        n_settled = settle_daily_cohorts(
            cohorts=[cohort],
            agents=agents,
            genotypes=genotypes,
            carrying_capacity=200,
            pop_cfg=pop_cfg,
            rng=rng,
            effect_sizes=effect_sizes,
        )

        assert n_settled > 0
        now_alive = set(np.where(agents['alive'])[0])
        newly_settled = sorted(now_alive & dead_before)

        for idx in newly_settled:
            # Genotypes should be all 1s (from our cohort)
            assert np.all(genotypes[idx] == 1), (
                f"Genotype at slot {idx} not properly copied"
            )
            # Resistance should be computed and positive
            assert agents[idx]['resistance'] > 0, (
                f"Resistance at slot {idx} is {agents[idx]['resistance']}"
            )

    def test_multiple_cohorts_settle_sequentially(self):
        """Multiple cohorts settle in order, respecting available slots."""
        K = 200
        agents, genotypes, effect_sizes = _make_population(50, K)
        pop_cfg = PopulationSection()
        rng = np.random.default_rng(42)

        cohorts = [
            _make_cohort(n_competent=300, spawn_day=10, rng=np.random.default_rng(i))
            for i in range(5)
        ]

        n_settled = settle_daily_cohorts(
            cohorts=cohorts,
            agents=agents,
            genotypes=genotypes,
            carrying_capacity=K,
            pop_cfg=pop_cfg,
            rng=rng,
            effect_sizes=effect_sizes,
        )

        n_alive = int(np.sum(agents['alive']))
        assert n_settled > 0
        assert n_alive <= K
        assert n_alive == 50 + n_settled

    def test_settle_with_default_effect_sizes(self):
        """settle_daily_cohorts works when effect_sizes=None (uses default)."""
        agents, genotypes, _ = _make_population(30, 100)
        pop_cfg = PopulationSection()
        rng = np.random.default_rng(42)

        cohort = _make_cohort(n_competent=200, rng=np.random.default_rng(99))
        n_settled = settle_daily_cohorts(
            cohorts=[cohort],
            agents=agents,
            genotypes=genotypes,
            carrying_capacity=100,
            pop_cfg=pop_cfg,
            rng=rng,
            effect_sizes=None,  # Should use default
        )
        assert n_settled > 0


# ═══════════════════════════════════════════════════════════════════════
# PLD TEMPERATURE DEPENDENCE TESTS
# ═══════════════════════════════════════════════════════════════════════


class TestPLDTemperatureDependence:
    """Verify PLD values match Hodin 2021 / spec expectations."""

    def test_pld_at_8C(self):
        """8°C → ~71 days (cold, Sitka AK)."""
        pld = pelagic_larval_duration(8.0)
        assert 68 <= pld <= 74, f"PLD at 8°C = {pld:.1f}, expected ~71"

    def test_pld_at_10_5C(self):
        """10.5°C → 63 days (reference, Hodin lab)."""
        pld = pelagic_larval_duration(10.5)
        assert pld == pytest.approx(63.0, abs=0.5), f"PLD at 10.5°C = {pld:.1f}, expected 63"

    def test_pld_at_14C(self):
        """14°C → ~53 days (warm, Monterey)."""
        pld = pelagic_larval_duration(14.0)
        assert 50 <= pld <= 56, f"PLD at 14°C = {pld:.1f}, expected ~53"

    def test_pld_decreases_with_temperature(self):
        """Warmer → shorter PLD (monotonic decrease)."""
        temps = [6.0, 8.0, 10.0, 12.0, 14.0, 16.0]
        plds = [pelagic_larval_duration(t) for t in temps]
        for i in range(len(plds) - 1):
            assert plds[i] > plds[i + 1], (
                f"PLD should decrease with temperature: "
                f"{temps[i]}°C={plds[i]:.1f}d, {temps[i+1]}°C={plds[i+1]:.1f}d"
            )

    def test_pld_clamped(self):
        """PLD clamped to [30, 150] range at extreme temperatures."""
        assert pelagic_larval_duration(-20.0) == 150.0  # Extreme cold → max
        assert pelagic_larval_duration(50.0) == 30.0    # Extreme hot → min
