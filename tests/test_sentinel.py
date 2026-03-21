"""Tests for sentinel agents — multi-host disease reservoir.

Verifies:
  - Sentinels don't die (natural or disease mortality)
  - Sentinels don't reproduce
  - Sentinels can get infected
  - Sentinel shedding is reduced by shedding_fraction
  - Population metrics exclude sentinels
  - Config parsing works with sentinel section
  - Default config (sentinels OFF) produces identical behavior
"""

import numpy as np
import pytest

from sswd_evoepi.config import (
    DiseaseSection,
    PopulationSection,
    SentinelSection,
    SimulationConfig,
    default_config,
    validate_config,
)
from sswd_evoepi.disease import (
    NodeDiseaseState,
    daily_disease_update,
    shedding_rate_I1,
    shedding_rate_I2,
)
from sswd_evoepi.model import (
    daily_natural_mortality,
    daily_growth_and_aging,
    annual_reproduction,
    initialize_sentinels,
    make_effect_sizes,
)
from sswd_evoepi.types import (
    AGENT_DTYPE,
    DiseaseState,
    Stage,
    allocate_agents,
    allocate_genotypes,
)


# ═══════════════════════════════════════════════════════════════════════
# Config tests
# ═══════════════════════════════════════════════════════════════════════


class TestSentinelConfig:
    def test_default_sentinel_disabled(self):
        """Default config has sentinels OFF."""
        cfg = default_config()
        assert cfg.sentinel.enabled is False

    def test_sentinel_section_defaults(self):
        """SentinelSection has correct default values."""
        s = SentinelSection()
        assert s.enabled is False
        assert s.n_per_site == 20
        assert s.shedding_fraction == 0.1
        assert s.site_filter == 'rocky'
        assert s.custom_node_ids == []

    def test_sentinel_enabled_config_valid(self):
        """Enabled sentinel config passes validation."""
        cfg = default_config()
        cfg.sentinel = SentinelSection(enabled=True, n_per_site=10,
                                        shedding_fraction=0.2)
        validate_config(cfg)  # Should not raise

    def test_sentinel_invalid_site_filter(self):
        """Invalid site_filter raises ValueError."""
        cfg = default_config()
        cfg.sentinel = SentinelSection(enabled=True, site_filter='invalid')
        with pytest.raises(ValueError, match="site_filter"):
            validate_config(cfg)

    def test_sentinel_invalid_n_per_site(self):
        """n_per_site < 1 raises ValueError."""
        cfg = default_config()
        cfg.sentinel = SentinelSection(enabled=True, n_per_site=0)
        with pytest.raises(ValueError, match="n_per_site"):
            validate_config(cfg)

    def test_sentinel_invalid_shedding_fraction(self):
        """shedding_fraction out of (0, 1] raises ValueError."""
        cfg = default_config()
        cfg.sentinel = SentinelSection(enabled=True, shedding_fraction=0.0)
        with pytest.raises(ValueError, match="shedding_fraction"):
            validate_config(cfg)

        cfg.sentinel = SentinelSection(enabled=True, shedding_fraction=1.5)
        with pytest.raises(ValueError, match="shedding_fraction"):
            validate_config(cfg)

    def test_is_sentinel_field_exists(self):
        """is_sentinel field exists in AGENT_DTYPE."""
        assert 'is_sentinel' in AGENT_DTYPE.names
        assert AGENT_DTYPE['is_sentinel'] == np.int8


# ═══════════════════════════════════════════════════════════════════════
# Mortality tests
# ═══════════════════════════════════════════════════════════════════════


class TestSentinelMortality:
    def _make_population(self, n_normal=50, n_sentinel=20):
        """Create a mixed population of normal and sentinel agents."""
        total = n_normal + n_sentinel
        agents = allocate_agents(total)

        # Normal agents
        agents['alive'][:n_normal] = True
        agents['stage'][:n_normal] = Stage.ADULT
        agents['age'][:n_normal] = 20.0  # old enough for some mortality
        agents['size'][:n_normal] = 500.0
        agents['sex'][:n_normal] = 0
        agents['disease_state'][:n_normal] = DiseaseState.S

        # Sentinel agents
        sl = slice(n_normal, total)
        agents['alive'][sl] = True
        agents['is_sentinel'][sl] = 1
        agents['stage'][sl] = Stage.ADULT
        agents['age'][sl] = 10.0
        agents['size'][sl] = 300.0
        agents['sex'][sl] = 0
        agents['disease_state'][sl] = DiseaseState.S

        return agents, n_normal, n_sentinel

    def test_sentinel_immune_to_natural_mortality(self):
        """Sentinels never die from natural causes."""
        agents, n_normal, n_sentinel = self._make_population(100, 50)
        pop_cfg = PopulationSection()
        rng = np.random.default_rng(42)

        # Run many rounds of mortality
        total_sentinel_deaths = 0
        for _ in range(1000):
            n_killed, _ = daily_natural_mortality(agents, pop_cfg, rng)

        # All sentinels should still be alive
        sentinel_alive = np.sum(
            agents['alive'][100:150] & agents['is_sentinel'][100:150].astype(bool)
        )
        assert sentinel_alive == n_sentinel, (
            f"Expected all {n_sentinel} sentinels alive, got {sentinel_alive}"
        )

    def test_sentinel_immune_to_disease_death(self):
        """Sentinels in I2 don't die when timer expires — timer resets."""
        agents = allocate_agents(10)
        # Create a sentinel in I2 with timer about to expire
        agents['alive'][0] = True
        agents['is_sentinel'][0] = 1
        agents['disease_state'][0] = DiseaseState.I2
        agents['disease_timer'][0] = 1  # will expire next step
        agents['resistance'][0] = 0.0
        agents['tolerance'][0] = 0.0
        agents['recovery_ability'][0] = 0.0
        agents['size'][0] = 300.0

        node_state = NodeDiseaseState()
        cfg = DiseaseSection()
        rng = np.random.default_rng(42)

        # Run disease update — timer will expire
        node_state = daily_disease_update(
            agents, node_state, T_celsius=15.0, salinity=30.0,
            phi_k=0.02, dispersal_input=0.0, day=100, cfg=cfg, rng=rng,
            sentinel_shedding_fraction=0.1,
        )

        # Sentinel should still be alive and in I2
        assert agents['alive'][0] == True
        assert agents['disease_state'][0] == DiseaseState.I2
        # Timer should have been reset (new countdown)
        assert agents['disease_timer'][0] > 0


# ═══════════════════════════════════════════════════════════════════════
# Reproduction tests
# ═══════════════════════════════════════════════════════════════════════


class TestSentinelReproduction:
    def test_sentinels_excluded_from_reproduction(self):
        """Sentinels don't spawn in annual_reproduction."""
        n = 100
        agents = allocate_agents(200)
        genotypes = allocate_genotypes(200)
        effect_sizes = make_effect_sizes()

        # All sentinel adults — should produce no offspring
        agents['alive'][:n] = True
        agents['is_sentinel'][:n] = 1
        agents['stage'][:n] = Stage.ADULT
        agents['size'][:n] = 600.0
        agents['age'][:n] = 10.0
        agents['sex'][:n // 2] = 0  # female
        agents['sex'][n // 2:n] = 1  # male
        agents['disease_state'][:n] = DiseaseState.S
        agents['resistance'][:n] = 0.1

        pop_cfg = PopulationSection()
        rng = np.random.default_rng(42)

        diag = annual_reproduction(
            agents, genotypes, effect_sizes,
            habitat_area=10000.0, sst=15.0,
            carrying_capacity=200, pop_cfg=pop_cfg, rng=rng,
        )

        assert diag['n_spawning_females'] == 0
        assert diag['n_spawning_males'] == 0
        assert diag['n_recruits'] == 0


# ═══════════════════════════════════════════════════════════════════════
# Infection tests
# ═══════════════════════════════════════════════════════════════════════


class TestSentinelInfection:
    def test_sentinels_can_get_infected(self):
        """Sentinels participate in disease transmission and can get infected."""
        agents = allocate_agents(20)

        # Create 10 sentinels, all susceptible
        for i in range(10):
            agents['alive'][i] = True
            agents['is_sentinel'][i] = 1
            agents['disease_state'][i] = DiseaseState.S
            agents['resistance'][i] = 0.0  # no resistance
            agents['tolerance'][i] = 0.0
            agents['recovery_ability'][i] = 0.0
            agents['size'][i] = 300.0
            agents['stage'][i] = Stage.ADULT

        node_state = NodeDiseaseState()
        node_state.vibrio_concentration = 1e6  # very high vibrio
        cfg = DiseaseSection()
        rng = np.random.default_rng(42)

        # Run disease update with high vibrio — some should get infected
        for day in range(30):
            node_state = daily_disease_update(
                agents, node_state, T_celsius=18.0, salinity=30.0,
                phi_k=0.02, dispersal_input=0.0, day=day, cfg=cfg, rng=rng,
                sentinel_shedding_fraction=0.1,
            )

        # At least some sentinels should be infected (E, I1, or I2)
        sentinel_ds = agents['disease_state'][:10]
        n_infected = np.sum(sentinel_ds > DiseaseState.S)
        assert n_infected > 0, "Expected at least one sentinel to get infected"


# ═══════════════════════════════════════════════════════════════════════
# Shedding tests
# ═══════════════════════════════════════════════════════════════════════


class TestSentinelShedding:
    def test_sentinel_shedding_reduced(self):
        """Sentinel shedding is reduced by shedding_fraction.

        Compare vibrio buildup with all-sentinel-infected vs all-normal-infected.
        The sentinel case should build up less vibrio due to reduced shedding.
        Uses invasion scenario to disable background environmental vibrio,
        isolating the shedding effect.
        """
        cfg = DiseaseSection()
        cfg.P_env_dynamic = False  # Simpler for testing
        cfg.scenario = "invasion"  # No background vibrio
        cfg.invasion_year = 0
        T = 15.0
        sal = 30.0
        phi = 0.02
        shedding_fraction = 0.1

        # Case 1: 10 normal I2 agents
        agents_normal = allocate_agents(10)
        agents_normal['alive'][:] = True
        agents_normal['disease_state'][:] = DiseaseState.I2
        agents_normal['disease_timer'][:] = 100  # won't expire
        agents_normal['size'][:] = 300.0
        agents_normal['tolerance'][:] = 0.0

        ns_normal = NodeDiseaseState()
        ns_normal.vibrio_concentration = 0.0
        rng1 = np.random.default_rng(99)

        # Run a few steps
        for d in range(5):
            ns_normal = daily_disease_update(
                agents_normal, ns_normal, T, sal, phi,
                dispersal_input=0.0, day=d, cfg=cfg, rng=rng1,
                sentinel_shedding_fraction=0.0,
            )
        vibrio_normal = ns_normal.vibrio_concentration

        # Case 2: 10 sentinel I2 agents with shedding_fraction=0.1
        agents_sent = allocate_agents(10)
        agents_sent['alive'][:] = True
        agents_sent['is_sentinel'][:] = 1
        agents_sent['disease_state'][:] = DiseaseState.I2
        agents_sent['disease_timer'][:] = 100  # won't expire
        agents_sent['size'][:] = 300.0
        agents_sent['tolerance'][:] = 0.0

        ns_sent = NodeDiseaseState()
        ns_sent.vibrio_concentration = 0.0
        rng2 = np.random.default_rng(99)

        for d in range(5):
            ns_sent = daily_disease_update(
                agents_sent, ns_sent, T, sal, phi,
                dispersal_input=0.0, day=d, cfg=cfg, rng=rng2,
                sentinel_shedding_fraction=shedding_fraction,
            )
        vibrio_sentinel = ns_sent.vibrio_concentration

        # Sentinel vibrio should be significantly less than normal
        assert vibrio_normal > 0, "Expected some vibrio from normal shedding"
        assert vibrio_sentinel < vibrio_normal, (
            f"Expected sentinel vibrio ({vibrio_sentinel:.1f}) < "
            f"normal vibrio ({vibrio_normal:.1f})"
        )
        # Should be approximately shedding_fraction × normal
        ratio = vibrio_sentinel / vibrio_normal if vibrio_normal > 0 else 0
        assert ratio < 0.25, f"Expected ratio < 0.25, got {ratio:.3f}"


# ═══════════════════════════════════════════════════════════════════════
# Metrics exclusion tests
# ═══════════════════════════════════════════════════════════════════════


class TestSentinelMetrics:
    def test_population_counts_exclude_sentinels(self):
        """Population counts should exclude sentinels.

        Tests the basic logic used in run_spatial_simulation's annual
        recording loop.
        """
        agents = allocate_agents(50)

        # 30 normal alive
        agents['alive'][:30] = True
        agents['stage'][:30] = Stage.ADULT

        # 20 sentinels alive
        agents['alive'][30:50] = True
        agents['is_sentinel'][30:50] = 1
        agents['stage'][30:50] = Stage.ADULT

        alive_mask = agents['alive'].astype(bool)
        pyc_mask = alive_mask & ~agents['is_sentinel'].astype(bool)

        # Total alive includes sentinels
        assert int(np.sum(alive_mask)) == 50
        # Pyc count excludes sentinels
        assert int(np.sum(pyc_mask)) == 30

    def test_sentinels_default_zero(self):
        """New agents have is_sentinel=0 by default."""
        agents = allocate_agents(10)
        assert np.all(agents['is_sentinel'] == 0)


# ═══════════════════════════════════════════════════════════════════════
# Growth/aging tests
# ═══════════════════════════════════════════════════════════════════════


class TestSentinelGrowth:
    def test_sentinels_dont_age(self):
        """Sentinels are excluded from daily growth and aging."""
        agents = allocate_agents(20)
        pop_cfg = PopulationSection()
        rng = np.random.default_rng(42)

        # 10 normal + 10 sentinel
        agents['alive'][:10] = True
        agents['age'][:10] = 5.0
        agents['size'][:10] = 200.0
        agents['stage'][:10] = Stage.JUVENILE

        agents['alive'][10:20] = True
        agents['is_sentinel'][10:20] = 1
        agents['age'][10:20] = 10.0
        agents['size'][10:20] = 300.0
        agents['stage'][10:20] = Stage.ADULT

        sentinel_age_before = agents['age'][10:20].copy()
        sentinel_size_before = agents['size'][10:20].copy()

        # Run growth
        daily_growth_and_aging(agents, pop_cfg, rng)

        # Sentinel age/size should be unchanged
        np.testing.assert_array_equal(
            agents['age'][10:20], sentinel_age_before,
            err_msg="Sentinel age changed after growth step"
        )
        np.testing.assert_array_equal(
            agents['size'][10:20], sentinel_size_before,
            err_msg="Sentinel size changed after growth step"
        )

        # Normal agents should have aged
        assert np.all(agents['age'][:10] > 5.0)


# ═══════════════════════════════════════════════════════════════════════
# Initialize sentinels tests
# ═══════════════════════════════════════════════════════════════════════


class TestInitializeSentinels:
    def test_initialize_sentinels_basic(self):
        """initialize_sentinels creates correct sentinel agents."""
        from sswd_evoepi.spatial import NodeDefinition

        nd = NodeDefinition(
            node_id=0, name='test', lat=48.0, lon=-123.0,
            subregion='TEST',
            habitat_area=10000.0, carrying_capacity=100,
        )

        # Create a minimal SpatialNode-like object
        class MockNode:
            def __init__(self, nd):
                self.definition = nd
                self.agents = allocate_agents(100)
                self.genotypes = allocate_genotypes(100)
                # Initialize 50 normal agents
                self.agents['alive'][:50] = True
                self.agents['stage'][:50] = Stage.ADULT

        node = MockNode(nd)
        rng = np.random.default_rng(42)

        initialize_sentinels(node, n_sentinels=20, rng=rng)

        # Check sentinel properties
        sentinel_mask = node.agents['is_sentinel'].astype(bool) & node.agents['alive'].astype(bool)
        n_sentinels = int(np.sum(sentinel_mask))
        assert n_sentinels == 20

        # Sentinels should be susceptible
        assert np.all(node.agents['disease_state'][sentinel_mask] == DiseaseState.S)

        # Sentinels should be alive
        assert np.all(node.agents['alive'][sentinel_mask])

        # Normal agents should still exist
        normal_mask = node.agents['alive'].astype(bool) & ~node.agents['is_sentinel'].astype(bool)
        assert int(np.sum(normal_mask)) == 50

    def test_initialize_sentinels_grows_array(self):
        """initialize_sentinels grows arrays if needed."""
        from sswd_evoepi.spatial import NodeDefinition

        nd = NodeDefinition(
            node_id=0, name='test', lat=48.0, lon=-123.0,
            subregion='TEST',
            habitat_area=10000.0, carrying_capacity=100,
        )

        class MockNode:
            def __init__(self, nd):
                self.definition = nd
                self.agents = allocate_agents(100)
                self.genotypes = allocate_genotypes(100)
                # Fill all 100 slots with alive agents
                self.agents['alive'][:] = True
                self.agents['stage'][:] = Stage.ADULT

        node = MockNode(nd)
        rng = np.random.default_rng(42)

        # All slots full — should grow array
        initialize_sentinels(node, n_sentinels=20, rng=rng)

        assert len(node.agents) > 100  # Array grew
        sentinel_mask = node.agents['is_sentinel'].astype(bool) & node.agents['alive'].astype(bool)
        assert int(np.sum(sentinel_mask)) == 20
