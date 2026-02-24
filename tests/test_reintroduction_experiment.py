"""Tests for the Monterey reintroduction experiment runner.

Covers:
  - Scenario generation (count, names, baseline)
  - 907-node network loading and Monterey site discovery
  - Config construction for treatment and baseline scenarios
  - Smoke test: single scenario completes end-to-end (small K, short run)

References:
  - experiments/reintroduction_monterey.py
  - scripts/genetic_backgrounds.py

Authors: Anton 🔬 & Willem Weertman
"""

import json
from pathlib import Path

import numpy as np
import pytest

from sswd_evoepi.config import ReleaseEvent, SimulationConfig
from sswd_evoepi.model import run_spatial_simulation
from sswd_evoepi.spatial import build_network, NodeDefinition

from scripts.genetic_backgrounds import get_background_by_name

# Import experiment infrastructure
from experiments.reintroduction_monterey import (
    MONTEREY_SITE_NAME,
    REGIONS,
    RELEASE_DAY,
    RESTORATION_LEVELS,
    SITES_PATH,
    ScenarioConfig,
    build_config_for_scenario,
    build_scenarios,
    load_node_definitions,
)


PROJECT_ROOT = Path(__file__).resolve().parent.parent


# ═══════════════════════════════════════════════════════════════════════
# SCENARIO GENERATION
# ═══════════════════════════════════════════════════════════════════════


class TestBuildScenarios:
    """Tests for the factorial scenario builder."""

    def test_build_scenarios_count(self):
        """19 unique scenarios × n_replicates."""
        scenarios = build_scenarios(n_replicates=1)
        assert len(scenarios) == 19  # 1 baseline + 3 levels × 6 genetics

        scenarios_3 = build_scenarios(n_replicates=3)
        assert len(scenarios_3) == 19 * 3

    def test_build_scenarios_baseline(self):
        """Baseline has no intervention."""
        scenarios = build_scenarios(n_replicates=1)
        baselines = [s for s in scenarios if s.name == 'baseline']
        assert len(baselines) == 1
        bl = baselines[0]
        assert bl.n_released == 0
        assert bl.restoration_level is None
        assert bl.genetics_name is None
        assert bl.allele_freqs is None

    def test_build_scenarios_names(self):
        """All 19 unique scenario names present."""
        scenarios = build_scenarios(n_replicates=1)
        names = {s.name for s in scenarios}
        assert 'baseline' in names

        # 3 levels × 6 genetics = 18 treatment names
        genetics_names = [
            'pre_sswd', 'survivors_2019', 'bred_1gen',
            'bred_2gen', 'bred_5gen', 'optimal',
        ]
        for level in RESTORATION_LEVELS:
            for gn in genetics_names:
                expected = f'{level}_{gn}'
                assert expected in names, f"Missing scenario: {expected}"

        assert len(names) == 19


# ═══════════════════════════════════════════════════════════════════════
# NETWORK LOADING
# ═══════════════════════════════════════════════════════════════════════


class TestLoadNodeDefinitions:
    """Tests for 907-node site loading."""

    @pytest.fixture(scope='class')
    def loaded_nodes(self):
        """Load once per class — reading 907 nodes is moderately expensive."""
        return load_node_definitions(K=100)

    def test_load_node_definitions(self, loaded_nodes):
        """Loads 907 nodes with correct types."""
        node_defs, sst_mapping, monterey_idx = loaded_nodes
        assert len(node_defs) == 907
        assert len(sst_mapping) == 907
        assert isinstance(monterey_idx, int)
        assert 0 <= monterey_idx < 907

    def test_monterey_site_index(self, loaded_nodes):
        """CA-C-043 found at the correct index with right coordinates."""
        node_defs, _, monterey_idx = loaded_nodes
        nd = node_defs[monterey_idx]
        assert nd.name == MONTEREY_SITE_NAME
        assert nd.lat == pytest.approx(36.62, abs=0.01)
        assert nd.subregion == 'CA-C'


# ═══════════════════════════════════════════════════════════════════════
# CONFIG CONSTRUCTION
# ═══════════════════════════════════════════════════════════════════════


class TestBuildConfig:
    """Tests for SimulationConfig construction from scenarios."""

    def test_build_config_for_scenario(self):
        """Treatment scenario produces config with release event."""
        bg = get_background_by_name('bred_2gen')
        scenario = ScenarioConfig(
            name='medium_bred_2gen',
            restoration_level='medium',
            n_released=500,
            genetics_name='bred_2gen',
            allele_freqs=bg['allele_freqs'],
            expected_traits=bg['expected_traits'],
            seed=42,
            replicate=0,
        )
        config = build_config_for_scenario(
            scenario, sst_mapping={0: 'Monterey'}, monterey_idx=0,
        )
        assert isinstance(config, SimulationConfig)
        assert len(config.release_events) == 1

        re = config.release_events[0]
        assert re.time_step == RELEASE_DAY
        assert re.node_id == 0
        assert re.n_individuals == 500
        assert re.genetics_mode == 'allele_freqs'
        assert re.allele_freqs is not None
        assert re.mark_released is True

    def test_build_config_baseline_no_release(self):
        """Baseline scenario has no release events."""
        scenario = ScenarioConfig(
            name='baseline',
            restoration_level=None,
            n_released=0,
            genetics_name=None,
            allele_freqs=None,
            expected_traits=None,
            seed=42,
            replicate=0,
        )
        config = build_config_for_scenario(
            scenario, sst_mapping={}, monterey_idx=0,
        )
        assert isinstance(config, SimulationConfig)
        assert len(config.release_events) == 0


# ═══════════════════════════════════════════════════════════════════════
# SMOKE TEST (end-to-end with tiny network)
# ═══════════════════════════════════════════════════════════════════════


class TestSingleScenarioRuns:
    """End-to-end smoke test on a minimal 3-node network.

    Uses K=100 and 5 years to keep runtime under ~60s.
    This does NOT test the full 907-node network (too expensive for CI)
    but validates that the experiment infrastructure wires together
    correctly: backgrounds → config → release events → simulation.
    """

    @pytest.fixture(scope='class')
    def mini_network(self):
        """Build a 3-node mini network."""
        nodes = []
        for i in range(3):
            nd = NodeDefinition(
                node_id=i,
                name=f'TEST-{i:03d}',
                lat=36.0 + i,
                lon=-122.0,
                subregion='CA-C',
                habitat_area=25000.0,
                carrying_capacity=100,
                is_fjord=False,
                sill_depth=np.inf,
                flushing_rate=0.5,
                mean_sst=14.0,
                sst_amplitude=3.0,
                sst_trend=0.02,
                salinity=32.0,
                depth_range=(5.0, 60.0),
            )
            nodes.append(nd)
        return build_network(nodes, D_L=400.0, D_P=15.0, seed=42)

    def test_single_scenario_runs(self, mini_network):
        """Treatment scenario runs to completion with release event."""
        bg = get_background_by_name('bred_2gen')
        release = ReleaseEvent(
            time_step=365 * 3,  # Release in year 3
            node_id=0,
            n_individuals=50,
            genetics_mode='allele_freqs',
            allele_freqs=bg['allele_freqs'],
            age_range=(365, 730),
            mark_released=True,
        )
        config = SimulationConfig(release_events=[release])

        result = run_spatial_simulation(
            network=mini_network,
            n_years=5,
            disease_year=2,
            initial_infected_per_node=3,
            seed=42,
            config=config,
        )

        # Basic sanity checks
        assert result.n_years == 5
        assert result.n_nodes == 3
        assert result.initial_total_pop > 0
        assert result.final_total_pop >= 0
        assert result.total_released == 50
        # Some released individuals should be tracked
        assert result.yearly_pop is not None
        assert result.yearly_pop.shape == (3, 5)
        # Genetic tracking present
        assert result.yearly_mean_resistance is not None
        assert result.yearly_mean_resistance.shape == (3, 5)
