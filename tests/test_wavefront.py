"""Tests for wavefront disease spread mechanics.

Chunk 4 of the wavefront implementation:
- Backward compatibility (wavefront_enabled=False → identical behavior)
- Wavefront seeding (only origin nodes get initial infections)
- P_env gating (unreached nodes have 0 background vibrio)
- Waterborne activation (pathogen dispersal triggers node activation)
- disease_arrival_day tracking
- activation_threshold behavior
"""
import numpy as np
import pytest

from sswd_evoepi.config import (
    SimulationConfig,
    DiseaseSection,
    validate_config,
)
from sswd_evoepi.disease import (
    NodeDiseaseState,
    daily_disease_update,
    update_vibrio_concentration,
)
from sswd_evoepi.spatial import make_5node_network


# ═══════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════

def _default_disease_cfg(**overrides) -> DiseaseSection:
    """Create a DiseaseSection with sensible defaults for testing."""
    kwargs = dict(scenario="ubiquitous", P_env_max=500.0, T_vbnc=12.0)
    kwargs.update(overrides)
    return DiseaseSection(**kwargs)


# ═══════════════════════════════════════════════════════════════════════
# Config validation tests
# ═══════════════════════════════════════════════════════════════════════

class TestWavefrontConfig:
    """Test wavefront config fields and validation."""

    def test_default_wavefront_disabled(self):
        cfg = DiseaseSection()
        assert cfg.wavefront_enabled is False
        assert cfg.disease_origin_nodes is None
        assert cfg.activation_threshold == 1.0

    def test_wavefront_enabled_requires_origin_nodes(self):
        cfg = SimulationConfig(
            disease=DiseaseSection(
                wavefront_enabled=True,
                disease_origin_nodes=None,
            )
        )
        with pytest.raises(ValueError, match="disease_origin_nodes required"):
            validate_config(cfg)

    def test_wavefront_enabled_requires_positive_threshold(self):
        cfg = SimulationConfig(
            disease=DiseaseSection(
                wavefront_enabled=True,
                disease_origin_nodes=[0],
                activation_threshold=-1.0,
            )
        )
        with pytest.raises(ValueError, match="activation_threshold must be positive"):
            validate_config(cfg)

    def test_wavefront_valid_config(self):
        cfg = SimulationConfig(
            disease=DiseaseSection(
                wavefront_enabled=True,
                disease_origin_nodes=[0, 1],
                activation_threshold=5.0,
            )
        )
        validate_config(cfg)  # Should not raise


# ═══════════════════════════════════════════════════════════════════════
# NodeDiseaseState tests
# ═══════════════════════════════════════════════════════════════════════

class TestNodeDiseaseState:
    def test_default_disease_reached(self):
        nds = NodeDiseaseState(node_id=0)
        assert nds.disease_reached is True
        assert nds.disease_arrival_day == -1

    def test_unreached_state(self):
        nds = NodeDiseaseState(node_id=5)
        nds.disease_reached = False
        assert nds.disease_reached is False


# ═══════════════════════════════════════════════════════════════════════
# P_env gating tests (unit level)
# ═══════════════════════════════════════════════════════════════════════

class TestPenvGating:
    def test_reached_has_env(self):
        """disease_reached=True → environmental_vibrio contributes."""
        cfg = _default_disease_cfg()
        P = update_vibrio_concentration(
            P_k=0.0, n_I1=0, n_I2=0, n_D_fresh=0,
            T_celsius=15.0, salinity=33.0, phi_k=0.5,
            dispersal_input=0.0, cfg=cfg, disease_reached=True,
        )
        assert P > 0.0

    def test_unreached_no_env(self):
        """disease_reached=False → no environmental vibrio."""
        cfg = _default_disease_cfg()
        P = update_vibrio_concentration(
            P_k=0.0, n_I1=0, n_I2=0, n_D_fresh=0,
            T_celsius=15.0, salinity=33.0, phi_k=0.5,
            dispersal_input=0.0, cfg=cfg, disease_reached=False,
        )
        assert P == 0.0

    def test_unreached_still_receives_dispersal(self):
        """Unreached nodes still get pathogen from neighbors."""
        cfg = _default_disease_cfg()
        P = update_vibrio_concentration(
            P_k=0.0, n_I1=0, n_I2=0, n_D_fresh=0,
            T_celsius=10.0, salinity=33.0, phi_k=0.5,
            dispersal_input=50.0, cfg=cfg, disease_reached=False,
        )
        assert P > 0.0

    def test_unreached_still_has_shedding(self):
        """Shedding works even at unreached nodes."""
        cfg = _default_disease_cfg()
        P = update_vibrio_concentration(
            P_k=0.0, n_I1=5, n_I2=2, n_D_fresh=1,
            T_celsius=15.0, salinity=33.0, phi_k=0.5,
            dispersal_input=0.0, cfg=cfg, disease_reached=False,
        )
        assert P > 0.0

    def test_backward_compat_default(self):
        """Default (no disease_reached arg) = reached = True."""
        cfg = _default_disease_cfg()
        P_explicit = update_vibrio_concentration(
            P_k=0.0, n_I1=0, n_I2=0, n_D_fresh=0,
            T_celsius=15.0, salinity=33.0, phi_k=0.5,
            dispersal_input=0.0, cfg=cfg, disease_reached=True,
        )
        P_default = update_vibrio_concentration(
            P_k=0.0, n_I1=0, n_I2=0, n_D_fresh=0,
            T_celsius=15.0, salinity=33.0, phi_k=0.5,
            dispersal_input=0.0, cfg=cfg,
        )
        assert P_explicit == P_default


# ═══════════════════════════════════════════════════════════════════════
# Integration tests using run_spatial_simulation
# ═══════════════════════════════════════════════════════════════════════

class TestWavefrontSpatial:
    """Integration tests using the 5-node network."""

    def test_no_wavefront_no_arrival_day(self):
        """wavefront_enabled=False → disease_arrival_day is None."""
        from sswd_evoepi.model import run_spatial_simulation
        network = make_5node_network(seed=42)
        config = SimulationConfig(disease=DiseaseSection(scenario="ubiquitous"))

        result = run_spatial_simulation(
            network=network, config=config,
            n_years=3, disease_year=1,
            initial_infected_per_node=3, seed=42,
        )
        assert result.disease_arrival_day is None

    def test_wavefront_origin_seeded_others_not(self):
        """Only origin node gets initial infections."""
        from sswd_evoepi.model import run_spatial_simulation
        network = make_5node_network(seed=42)

        config = SimulationConfig(
            disease=DiseaseSection(
                scenario="ubiquitous",
                wavefront_enabled=True,
                disease_origin_nodes=[4],  # Monterey (southernmost)
                activation_threshold=1e9,  # impossibly high — no spread
            )
        )

        result = run_spatial_simulation(
            network=network, config=config,
            n_years=3, disease_year=1,
            initial_infected_per_node=5, seed=42,
        )

        assert result.disease_arrival_day is not None
        # Origin node (4) should be reached
        assert result.disease_arrival_day[4] == 1 * 365
        # With impossibly high threshold, non-origin nodes should NOT be reached
        for i in range(4):
            assert result.disease_arrival_day[i] == -1, (
                f"Node {i} should not be reached with threshold=1e9"
            )

    def test_wavefront_arrival_day_ordering(self):
        """With nonzero D matrix, disease spreads — farther nodes later."""
        from sswd_evoepi.model import run_spatial_simulation
        from sswd_evoepi.spatial import NodeDefinition, build_network

        # Create a 3-node linear chain with SHORT distances (within max_range=50km)
        # so the D matrix actually has nonzero entries
        node_defs = [
            NodeDefinition(node_id=0, name="North", lat=48.5, lon=-123.0,
                subregion="N", habitat_area=1000, carrying_capacity=200,
                mean_sst=9.0, sst_amplitude=3.0, sst_trend=0.0,
                salinity=33.0, depth_range=(5, 30)),
            NodeDefinition(node_id=1, name="Mid", lat=48.35, lon=-123.0,
                subregion="M", habitat_area=1000, carrying_capacity=200,
                mean_sst=11.0, sst_amplitude=3.0, sst_trend=0.0,
                salinity=33.0, depth_range=(5, 30)),
            NodeDefinition(node_id=2, name="South", lat=48.2, lon=-123.0,
                subregion="S", habitat_area=1000, carrying_capacity=200,
                mean_sst=13.0, sst_amplitude=3.0, sst_trend=0.0,
                salinity=33.0, depth_range=(5, 30)),
        ]
        # ~17km and ~17km between adjacent nodes — within max_range=50km
        network = build_network(node_defs, D_P=15.0, seed=42)

        # Verify D matrix is nonzero
        assert network.D.sum() > 0, "D matrix should have nonzero entries"

        config = SimulationConfig(
            disease=DiseaseSection(
                scenario="ubiquitous",
                wavefront_enabled=True,
                disease_origin_nodes=[2],  # South (warmest)
                activation_threshold=0.001,
            )
        )

        result = run_spatial_simulation(
            network=network, config=config,
            n_years=5, disease_year=1,
            initial_infected_per_node=5, seed=42,
        )

        assert result.disease_arrival_day is not None
        origin_day = result.disease_arrival_day[2]
        assert origin_day == 365

        # Disease should spread to at least the adjacent node
        reached = [i for i in range(3) if result.disease_arrival_day[i] >= 0]
        assert len(reached) >= 2, (
            f"Expected disease to spread to at least 2 nodes, got {reached}"
        )

        # Mid node should arrive after South
        if result.disease_arrival_day[1] >= 0:
            assert result.disease_arrival_day[1] > origin_day
        # North should arrive after Mid (if reached)
        if result.disease_arrival_day[0] >= 0 and result.disease_arrival_day[1] >= 0:
            assert result.disease_arrival_day[0] >= result.disease_arrival_day[1]

    def test_wavefront_no_initial_vibrio(self):
        """With wavefront enabled, non-origin nodes should start with 0 vibrio."""
        from sswd_evoepi.model import run_spatial_simulation
        network = make_5node_network(seed=42)

        config = SimulationConfig(
            disease=DiseaseSection(
                scenario="ubiquitous",
                wavefront_enabled=True,
                disease_origin_nodes=[4],
                activation_threshold=1e9,
            )
        )

        # Check initial vibrio — run for just 1 year (spinup, no disease yet)
        result = run_spatial_simulation(
            network=network, config=config,
            n_years=1, disease_year=None,  # no disease
            initial_infected_per_node=0, seed=42,
        )

        # All nodes should have 0 vibrio after spinup with wavefront enabled
        for node in network.nodes:
            assert node.vibrio_concentration == 0.0, (
                f"Node {node.name} should have 0 vibrio with wavefront enabled"
            )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
