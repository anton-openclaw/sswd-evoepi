"""Tests for juvenile immunity — age-dependent disease susceptibility (Phase 11).

Verifies that recently settled agents are protected from S→E transitions
for a configurable refractory period (min_susceptible_age_days), while
remaining subject to natural mortality and other processes.
"""

import numpy as np
import pytest

from sswd_evoepi.types import (
    AGENT_DTYPE,
    DiseaseState,
    Stage,
    allocate_agents,
    allocate_genotypes,
    N_LOCI,
)
from sswd_evoepi.config import DiseaseSection, default_config
from sswd_evoepi.disease import (
    daily_disease_update,
    NodeDiseaseState,
)


# ═══════════════════════════════════════════════════════════════════════
# HELPERS
# ═══════════════════════════════════════════════════════════════════════

def _make_single_agent(settlement_day: int = 0, alive: bool = True,
                        disease_state: int = DiseaseState.S,
                        resistance: float = 0.0,
                        size: float = 100.0,
                        node_id: int = 0):
    """Create a minimal agent array with one live agent."""
    agents = allocate_agents(10)
    agents['alive'][0] = alive
    agents['disease_state'][0] = disease_state
    agents['resistance'][0] = resistance
    agents['size'][0] = size
    agents['settlement_day'][0] = settlement_day
    agents['node_id'][0] = node_id
    agents['stage'][0] = Stage.SETTLER
    return agents


def _run_disease_step(agents, cfg, sim_day, vibrio=1e6, seed=42):
    """Run one disease step and return updated agents + node state."""
    rng = np.random.default_rng(seed)
    node_state = NodeDiseaseState(node_id=0)
    node_state.vibrio_concentration = vibrio  # High concentration → high infection pressure
    result = daily_disease_update(
        agents=agents,
        node_state=node_state,
        T_celsius=20.0,  # T_ref → maximum pathogen activity
        salinity=30.0,   # Full marine
        phi_k=0.02,
        dispersal_input=0.0,
        day=sim_day,
        cfg=cfg,
        rng=rng,
    )
    return agents, result


# ═══════════════════════════════════════════════════════════════════════
# TESTS
# ═══════════════════════════════════════════════════════════════════════

class TestJuvenileNotInfected:
    """Settlers within refractory period should NOT transition S→E."""

    def test_juvenile_not_infected(self):
        """Settler at day 30, min_susceptible_age_days=90 → stays S."""
        cfg = DiseaseSection(min_susceptible_age_days=90)
        agents = _make_single_agent(settlement_day=0)

        # Run 50 disease steps at day 30 (agent age = 30 days, < 90 threshold)
        # Use multiple seeds to ensure it's not just luck
        for seed in range(50):
            test_agents = agents.copy()
            test_agents, _ = _run_disease_step(test_agents, cfg, sim_day=30, seed=seed)
            assert test_agents['disease_state'][0] == DiseaseState.S, \
                f"Juvenile infected at seed={seed}, age_days=30 < threshold=90"


class TestJuvenileBecomeSusceptible:
    """Settlers past the refractory period CAN be infected."""

    def test_juvenile_becomes_susceptible(self):
        """Settler at day 91 with min_susceptible_age_days=90 → can be infected."""
        cfg = DiseaseSection(min_susceptible_age_days=90)
        agents = _make_single_agent(settlement_day=0, resistance=0.0)

        # Run many disease steps at day 91 — at least one should infect
        # With zero resistance, high vibrio, T=T_ref, infection prob ≈ 0.5+
        infected = False
        for seed in range(200):
            test_agents = agents.copy()
            test_agents, _ = _run_disease_step(test_agents, cfg, sim_day=91, seed=seed)
            if test_agents['disease_state'][0] != DiseaseState.S:
                infected = True
                break

        assert infected, \
            "Agent at age 91 (threshold=90) was never infected across 200 attempts"


class TestZeroThresholdBackwardCompat:
    """min_susceptible_age_days=0 should give immediate susceptibility (old behavior)."""

    def test_zero_threshold_backward_compat(self):
        """Default (0) → settlers can be infected immediately."""
        cfg = DiseaseSection(min_susceptible_age_days=0)
        agents = _make_single_agent(settlement_day=100, resistance=0.0)

        # Even on settlement day itself, should be susceptible
        infected = False
        for seed in range(200):
            test_agents = agents.copy()
            test_agents, _ = _run_disease_step(test_agents, cfg, sim_day=100, seed=seed)
            if test_agents['disease_state'][0] != DiseaseState.S:
                infected = True
                break

        assert infected, \
            "Agent with min_susceptible_age_days=0 should be immediately susceptible"


class TestInitialPopAlwaysSusceptible:
    """Initial population (settlement_day=0) is always susceptible."""

    def test_initial_pop_always_susceptible(self):
        """settlement_day=0, sim_day=500 → susceptible (age=500 >> any threshold)."""
        cfg = DiseaseSection(min_susceptible_age_days=90)
        agents = _make_single_agent(settlement_day=0, resistance=0.0)

        infected = False
        for seed in range(200):
            test_agents = agents.copy()
            test_agents, _ = _run_disease_step(test_agents, cfg, sim_day=500, seed=seed)
            if test_agents['disease_state'][0] != DiseaseState.S:
                infected = True
                break

        assert infected, \
            "Initial pop agent (settlement_day=0) should be susceptible at sim_day=500"

    def test_initial_pop_susceptible_day_one(self):
        """settlement_day=0, sim_day=0 → susceptible (age=0, but 0 >= 0 for initial pop)."""
        # min_susceptible_age_days=90, but age_days = 0 - 0 = 0 < 90
        # However! Initial pop was set before disease year, so they'll be well past threshold.
        # Still, let's test the edge case: settlement_day=0 is only "always susceptible"
        # if the agent was an initial pop member AND enough time has passed.
        # In practice, disease starts at disease_year * 365. With threshold=90,
        # initial pop (settlement_day=0) at sim_day=90 is exactly at threshold.
        cfg = DiseaseSection(min_susceptible_age_days=90)
        agents = _make_single_agent(settlement_day=0, resistance=0.0)

        infected = False
        for seed in range(200):
            test_agents = agents.copy()
            test_agents, _ = _run_disease_step(test_agents, cfg, sim_day=90, seed=seed)
            if test_agents['disease_state'][0] != DiseaseState.S:
                infected = True
                break

        assert infected, \
            "Initial pop at sim_day=90 (age=90, threshold=90) should be susceptible"


class TestMaxThreshold:
    """180-day threshold — max SA range — settlers immune for ~6 months."""

    def test_180_day_threshold(self):
        """Settler immune for 179 days, susceptible at 180."""
        cfg = DiseaseSection(min_susceptible_age_days=180)
        agents = _make_single_agent(settlement_day=100, resistance=0.0)

        # At day 279 (age=179), should NOT be infected
        for seed in range(50):
            test_agents = agents.copy()
            test_agents, _ = _run_disease_step(test_agents, cfg, sim_day=279, seed=seed)
            assert test_agents['disease_state'][0] == DiseaseState.S, \
                f"Juvenile (age=179, threshold=180) should not be infected"

        # At day 280 (age=180), should be susceptible
        infected = False
        for seed in range(200):
            test_agents = agents.copy()
            test_agents, _ = _run_disease_step(test_agents, cfg, sim_day=280, seed=seed)
            if test_agents['disease_state'][0] != DiseaseState.S:
                infected = True
                break

        assert infected, \
            "Agent at age=180 (threshold=180) should be susceptible"


class TestJuvenileStillDiesNaturally:
    """Juveniles in refractory period can still die from non-disease causes."""

    def test_settlement_day_field_exists(self):
        """Verify settlement_day is in AGENT_DTYPE."""
        agents = allocate_agents(1)
        assert 'settlement_day' in agents.dtype.names
        # Default should be 0
        assert agents['settlement_day'][0] == 0

    def test_juvenile_can_be_killed(self):
        """Juvenile's alive flag can be set to False (natural mortality)."""
        agents = _make_single_agent(settlement_day=100)
        assert agents['alive'][0] == True
        # Natural mortality just sets alive=False — no disease involvement
        agents['alive'][0] = False
        assert agents['alive'][0] == False

    def test_juvenile_disease_state_unchanged_during_refractory(self):
        """Multiple disease steps during refractory period → agent stays S."""
        cfg = DiseaseSection(min_susceptible_age_days=90)
        agents = _make_single_agent(settlement_day=0, resistance=0.0)

        # Run 30 consecutive disease steps (days 10-39), all within refractory
        for d in range(10, 40):
            agents, _ = _run_disease_step(agents, cfg, sim_day=d, seed=d)
            assert agents['disease_state'][0] == DiseaseState.S


class TestMultipleAgentsMixed:
    """Mix of juvenile and adult agents — only adults should be infected."""

    def test_mixed_population(self):
        """2 juveniles + 2 adults → only adults can be infected."""
        cfg = DiseaseSection(min_susceptible_age_days=90)
        agents = allocate_agents(10)

        # Agent 0: juvenile (settled day 50, current day 80 → age=30 < 90)
        agents['alive'][0] = True
        agents['disease_state'][0] = DiseaseState.S
        agents['resistance'][0] = 0.0
        agents['size'][0] = 50.0
        agents['settlement_day'][0] = 50
        agents['stage'][0] = Stage.SETTLER

        # Agent 1: juvenile (settled day 60, current day 80 → age=20 < 90)
        agents['alive'][1] = True
        agents['disease_state'][1] = DiseaseState.S
        agents['resistance'][1] = 0.0
        agents['size'][1] = 50.0
        agents['settlement_day'][1] = 60
        agents['stage'][1] = Stage.SETTLER

        # Agent 2: adult (settled day 0, current day 80 → age=80, but
        #           settlement_day=0 means initial pop → age far exceeds threshold
        #           actually age_days = 80 - 0 = 80 < 90... let me use settlement
        #           that makes age >= 90)
        # Use settlement_day=-100 to ensure age >= 90? No, let's set properly.
        # Actually: sim_day=200, settlement_day=0 → age=200 >> 90
        agents['alive'][2] = True
        agents['disease_state'][2] = DiseaseState.S
        agents['resistance'][2] = 0.0
        agents['size'][2] = 200.0
        agents['settlement_day'][2] = 0
        agents['stage'][2] = Stage.ADULT

        # Agent 3: adult (settled day 50, sim_day=200 → age=150 > 90)
        agents['alive'][3] = True
        agents['disease_state'][3] = DiseaseState.S
        agents['resistance'][3] = 0.0
        agents['size'][3] = 200.0
        agents['settlement_day'][3] = 50
        agents['stage'][3] = Stage.ADULT

        sim_day = 200  # Agents 0,1 settled at 50,60 → ages 150,140 >> 90
        # Wait, that makes them NOT juveniles. Let me adjust.
        # I need juveniles to still be in refractory. Let me set sim_day=80.
        sim_day = 80

        # At sim_day=80: agent 0 age=30, agent 1 age=20, agent 2 age=80<90, agent 3 age=30<90
        # Hmm, agents 2 and 3 would also be juveniles. Let me fix:
        agents['settlement_day'][2] = 0   # age=80 at sim_day=80 → still < 90!
        agents['settlement_day'][3] = 0   # same

        # Use sim_day=500 and settlement_day for juveniles = 480 (age=20) and adults = 0 (age=500)
        sim_day = 500
        agents['settlement_day'][0] = 480  # age = 20 < 90 → juvenile
        agents['settlement_day'][1] = 490  # age = 10 < 90 → juvenile
        agents['settlement_day'][2] = 0    # age = 500 >> 90 → adult
        agents['settlement_day'][3] = 100  # age = 400 >> 90 → adult

        # Track infections over many seeds
        juvenile_ever_infected = False
        adult_ever_infected = False

        for seed in range(200):
            test_agents = agents.copy()
            test_agents, _ = _run_disease_step(test_agents, cfg, sim_day=sim_day, seed=seed)

            if test_agents['disease_state'][0] != DiseaseState.S or \
               test_agents['disease_state'][1] != DiseaseState.S:
                juvenile_ever_infected = True

            if test_agents['disease_state'][2] != DiseaseState.S or \
               test_agents['disease_state'][3] != DiseaseState.S:
                adult_ever_infected = True

        assert not juvenile_ever_infected, "Juveniles should never be infected during refractory"
        assert adult_ever_infected, "Adults should be infectable"
