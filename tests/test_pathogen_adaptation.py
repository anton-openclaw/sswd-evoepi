"""Tests for pathogen thermal adaptation (Chunk 1).

Tests cover:
  - adapt_pathogen_thermal() function behavior (pool-driven)
  - Config validation for new fields
  - environmental_vibrio() with T_vbnc_local override
  - Backward compatibility (pathogen_adaptation=False)
"""

import numpy as np
import pytest

from sswd_evoepi.config import DiseaseSection, SimulationConfig, validate_config
from sswd_evoepi.disease import (
    NodeDiseaseState,
    adapt_pathogen_thermal,
    environmental_vibrio,
)


def _make_cfg(**kwargs) -> DiseaseSection:
    """Create a DiseaseSection with pathogen adaptation enabled by default."""
    defaults = dict(
        pathogen_adaptation=True,
        T_vbnc_initial=12.0,
        T_vbnc_min=6.0,
        pathogen_adapt_rate=0.001,
        pathogen_revert_rate=0.0,
        P_adapt_half=500.0,
    )
    defaults.update(kwargs)
    return DiseaseSection(**defaults)


def _make_node_state(T_vbnc_local: float = 12.0) -> NodeDiseaseState:
    """Create a NodeDiseaseState with a given local threshold."""
    nds = NodeDiseaseState()
    nds.T_vbnc_local = T_vbnc_local
    return nds


# ═══════════════════════════════════════════════════════════════════════
# adapt_pathogen_thermal() tests
# ═══════════════════════════════════════════════════════════════════════


class TestAdaptationDisabled:
    """When pathogen_adaptation=False, T_vbnc_local should never change."""

    def test_adaptation_disabled_by_default(self):
        cfg = DiseaseSection()  # pathogen_adaptation=False by default
        nds = _make_node_state(12.0)
        adapt_pathogen_thermal(nds, T_celsius=8.0, P_env_pool=1000.0, cfg=cfg)
        assert nds.T_vbnc_local == 12.0

    def test_adaptation_disabled_explicit(self):
        cfg = _make_cfg(pathogen_adaptation=False)
        nds = _make_node_state(12.0)
        adapt_pathogen_thermal(nds, T_celsius=5.0, P_env_pool=1000.0, cfg=cfg)
        assert nds.T_vbnc_local == 12.0


class TestAdaptationReducesThreshold:
    """When temp < T_vbnc_local and P_env_pool > 0, threshold should decrease."""

    def test_adaptation_reduces_threshold(self):
        cfg = _make_cfg(pathogen_adapt_rate=0.01)
        nds = _make_node_state(12.0)
        # T=8°C < 12°C, P_env_pool=500 → pool_factor = 500/500 = 1.0
        adapt_pathogen_thermal(nds, T_celsius=8.0, P_env_pool=500.0, cfg=cfg)
        assert nds.T_vbnc_local < 12.0

    def test_threshold_decreases_with_low_pool(self):
        cfg = _make_cfg(pathogen_adapt_rate=0.01)
        nds = _make_node_state(12.0)
        # Small pool still adapts
        adapt_pathogen_thermal(nds, T_celsius=8.0, P_env_pool=50.0, cfg=cfg)
        assert nds.T_vbnc_local < 12.0

    def test_threshold_decreases_with_high_pool(self):
        cfg = _make_cfg(pathogen_adapt_rate=0.01)
        nds = _make_node_state(12.0)
        # Large pool (saturated)
        adapt_pathogen_thermal(nds, T_celsius=8.0, P_env_pool=5000.0, cfg=cfg)
        assert nds.T_vbnc_local < 12.0


class TestNoAdaptationAboveThreshold:
    """When temp >= T_vbnc_local, no adaptation should occur."""

    def test_no_adaptation_above_threshold(self):
        cfg = _make_cfg()
        nds = _make_node_state(12.0)
        adapt_pathogen_thermal(nds, T_celsius=12.0, P_env_pool=1000.0, cfg=cfg)
        assert nds.T_vbnc_local == 12.0

    def test_no_adaptation_well_above_threshold(self):
        cfg = _make_cfg()
        nds = _make_node_state(12.0)
        adapt_pathogen_thermal(nds, T_celsius=20.0, P_env_pool=5000.0, cfg=cfg)
        assert nds.T_vbnc_local == 12.0


class TestNoAdaptationZeroPool:
    """When P_env_pool == 0, adaptation branch should NOT activate."""

    def test_no_adaptation_zero_pool(self):
        cfg = _make_cfg()
        nds = _make_node_state(12.0)
        adapt_pathogen_thermal(nds, T_celsius=5.0, P_env_pool=0.0, cfg=cfg)
        assert nds.T_vbnc_local == 12.0

    def test_no_adaptation_zero_pool_even_cold(self):
        """Zero pool means no adaptation even with big temp gap."""
        cfg = _make_cfg(pathogen_adapt_rate=1.0)
        nds = _make_node_state(12.0)
        adapt_pathogen_thermal(nds, T_celsius=0.0, P_env_pool=0.0, cfg=cfg)
        assert nds.T_vbnc_local == 12.0


class TestAdaptationWithPoolNoInfectedHosts:
    """Adaptation occurs with P_env_pool > 0 even when there are zero infected hosts.

    This is the key behavior change: the bacterial community drives adaptation,
    not infected host prevalence.
    """

    def test_pool_drives_adaptation_without_hosts(self):
        cfg = _make_cfg(pathogen_adapt_rate=0.01, P_adapt_half=500.0)
        nds = _make_node_state(12.0)
        # P_env_pool=1000 but conceptually zero infected hosts
        # (we no longer pass host counts at all)
        adapt_pathogen_thermal(nds, T_celsius=8.0, P_env_pool=1000.0, cfg=cfg)
        assert nds.T_vbnc_local < 12.0

    def test_pool_adaptation_persists_after_crash(self):
        """Simulate: hosts crash to zero but pool persists → adaptation continues."""
        cfg = _make_cfg(pathogen_adapt_rate=0.01, P_adapt_half=500.0)
        nds = _make_node_state(12.0)
        # 100 days of adaptation with only pool (no hosts needed)
        for _ in range(100):
            adapt_pathogen_thermal(nds, T_celsius=8.0, P_env_pool=600.0, cfg=cfg)
        assert nds.T_vbnc_local < 11.0  # significant adaptation


class TestAdaptationRespectsMinimum:
    """T_vbnc_local should never go below T_vbnc_min."""

    def test_adaptation_respects_minimum(self):
        cfg = _make_cfg(
            pathogen_adapt_rate=10.0,  # Very aggressive rate
            T_vbnc_min=8.0,
        )
        nds = _make_node_state(12.0)
        # Extreme: high rate, big pool, big temp gap
        adapt_pathogen_thermal(nds, T_celsius=0.0, P_env_pool=5000.0, cfg=cfg)
        assert nds.T_vbnc_local == 8.0

    def test_adaptation_minimum_reached_incrementally(self):
        cfg = _make_cfg(
            pathogen_adapt_rate=1.0,
            T_vbnc_min=10.0,
        )
        nds = _make_node_state(12.0)
        for _ in range(1000):
            adapt_pathogen_thermal(nds, T_celsius=5.0, P_env_pool=1000.0, cfg=cfg)
        assert nds.T_vbnc_local == pytest.approx(10.0)


class TestReversionWithoutBacteria:
    """When P_env_pool == 0, T_vbnc_local should revert toward T_vbnc_initial."""

    def test_reversion_when_pool_zero(self):
        cfg = _make_cfg(pathogen_revert_rate=0.1)
        nds = _make_node_state(10.0)  # already adapted (below initial=12)
        adapt_pathogen_thermal(nds, T_celsius=8.0, P_env_pool=0.0, cfg=cfg)
        assert nds.T_vbnc_local == pytest.approx(10.1)

    def test_reversion_capped_at_initial(self):
        cfg = _make_cfg(pathogen_revert_rate=5.0, T_vbnc_initial=12.0)
        nds = _make_node_state(11.5)
        adapt_pathogen_thermal(nds, T_celsius=8.0, P_env_pool=0.0, cfg=cfg)
        assert nds.T_vbnc_local == 12.0

    def test_no_reversion_when_pool_present(self):
        """Reversion only happens when P_env_pool == 0."""
        cfg = _make_cfg(pathogen_revert_rate=0.1)
        nds = _make_node_state(10.0)
        # Pool > 0 but temp above threshold → no adaptation, but also no reversion
        adapt_pathogen_thermal(nds, T_celsius=15.0, P_env_pool=100.0, cfg=cfg)
        assert nds.T_vbnc_local == 10.0

    def test_no_reversion_at_initial(self):
        """No reversion when already at T_vbnc_initial."""
        cfg = _make_cfg(pathogen_revert_rate=0.1, T_vbnc_initial=12.0)
        nds = _make_node_state(12.0)
        adapt_pathogen_thermal(nds, T_celsius=8.0, P_env_pool=0.0, cfg=cfg)
        assert nds.T_vbnc_local == 12.0

    def test_revert_rate_zero_means_no_reversion(self):
        """Default revert_rate=0.0 means no reversion even when pool is zero."""
        cfg = _make_cfg(pathogen_revert_rate=0.0)
        nds = _make_node_state(10.0)  # adapted below initial=12
        for _ in range(1000):
            adapt_pathogen_thermal(nds, T_celsius=8.0, P_env_pool=0.0, cfg=cfg)
        assert nds.T_vbnc_local == 10.0  # unchanged — reversion disabled


class TestAdaptationScaling:
    """Adaptation rate should scale with pool_factor and temperature gap."""

    def test_adaptation_rate_scales_with_pool_half_saturation(self):
        """Half-saturation: pool=P_adapt_half → pool_factor=1.0 (capped)."""
        cfg = _make_cfg(pathogen_adapt_rate=0.01, P_adapt_half=500.0)

        # Pool = half of P_adapt_half → factor = 250/500 = 0.5
        nds_half = _make_node_state(12.0)
        adapt_pathogen_thermal(nds_half, T_celsius=8.0, P_env_pool=250.0, cfg=cfg)
        delta_half = 12.0 - nds_half.T_vbnc_local

        # Pool = P_adapt_half → factor = 500/500 = 1.0
        nds_full = _make_node_state(12.0)
        adapt_pathogen_thermal(nds_full, T_celsius=8.0, P_env_pool=500.0, cfg=cfg)
        delta_full = 12.0 - nds_full.T_vbnc_local

        assert delta_full > delta_half
        ratio = delta_full / delta_half
        assert ratio == pytest.approx(2.0, rel=0.01)

    def test_adaptation_saturates_above_half(self):
        """Pool >> P_adapt_half → pool_factor capped at 1.0."""
        cfg = _make_cfg(pathogen_adapt_rate=0.01, P_adapt_half=500.0)

        # Pool = P_adapt_half → factor = 1.0
        nds_at = _make_node_state(12.0)
        adapt_pathogen_thermal(nds_at, T_celsius=8.0, P_env_pool=500.0, cfg=cfg)
        delta_at = 12.0 - nds_at.T_vbnc_local

        # Pool = 10x P_adapt_half → factor = min(10, 1.0) = 1.0
        nds_over = _make_node_state(12.0)
        adapt_pathogen_thermal(nds_over, T_celsius=8.0, P_env_pool=5000.0, cfg=cfg)
        delta_over = 12.0 - nds_over.T_vbnc_local

        # Both should be equal (saturated)
        assert delta_over == pytest.approx(delta_at, rel=1e-9)

    def test_adaptation_rate_scales_with_temp_gap(self):
        cfg = _make_cfg(pathogen_adapt_rate=0.01, P_adapt_half=500.0)

        # Small gap (11°C vs 12°C threshold → gap=1)
        nds_small = _make_node_state(12.0)
        adapt_pathogen_thermal(nds_small, T_celsius=11.0, P_env_pool=500.0, cfg=cfg)
        delta_small = 12.0 - nds_small.T_vbnc_local

        # Large gap (4°C vs 12°C threshold → gap=8)
        nds_large = _make_node_state(12.0)
        adapt_pathogen_thermal(nds_large, T_celsius=4.0, P_env_pool=500.0, cfg=cfg)
        delta_large = 12.0 - nds_large.T_vbnc_local

        assert delta_large > delta_small
        ratio = delta_large / delta_small
        assert ratio == pytest.approx(8.0, rel=0.01)


# ═══════════════════════════════════════════════════════════════════════
# Config validation tests
# ═══════════════════════════════════════════════════════════════════════


class TestConfigValidation:
    """Validation should catch invalid pathogen adaptation params."""

    def test_config_validation_min_greater_than_initial(self):
        """T_vbnc_min > T_vbnc_initial should raise ValueError."""
        config = SimulationConfig()
        config.disease.pathogen_adaptation = True
        config.disease.T_vbnc_min = 15.0
        config.disease.T_vbnc_initial = 12.0
        with pytest.raises(ValueError, match="T_vbnc_min"):
            validate_config(config)

    def test_config_validation_negative_adapt_rate(self):
        config = SimulationConfig()
        config.disease.pathogen_adaptation = True
        config.disease.pathogen_adapt_rate = -0.001
        with pytest.raises(ValueError, match="pathogen_adapt_rate"):
            validate_config(config)

    def test_config_validation_negative_revert_rate(self):
        config = SimulationConfig()
        config.disease.pathogen_adaptation = True
        config.disease.pathogen_revert_rate = -0.001
        with pytest.raises(ValueError, match="pathogen_revert_rate"):
            validate_config(config)

    def test_config_validation_valid_params(self):
        """Valid params should not raise."""
        config = SimulationConfig()
        config.disease.pathogen_adaptation = True
        config.disease.T_vbnc_initial = 12.0
        config.disease.T_vbnc_min = 6.0
        config.disease.pathogen_adapt_rate = 0.001
        config.disease.pathogen_revert_rate = 0.0
        config.disease.P_adapt_half = 500.0
        validate_config(config)  # should not raise

    def test_config_validation_skipped_when_disabled(self):
        """When pathogen_adaptation=False, invalid T_vbnc_min is OK."""
        config = SimulationConfig()
        config.disease.pathogen_adaptation = False
        config.disease.T_vbnc_min = 15.0  # invalid but ignored
        config.disease.T_vbnc_initial = 12.0
        validate_config(config)  # should not raise

    def test_P_adapt_half_default(self):
        """P_adapt_half should default to 500.0."""
        cfg = DiseaseSection()
        assert cfg.P_adapt_half == 500.0

    def test_P_adapt_half_invalid_zero(self):
        """P_adapt_half must be > 0."""
        config = SimulationConfig()
        config.disease.pathogen_adaptation = True
        config.disease.P_adapt_half = 0.0
        with pytest.raises(ValueError, match="P_adapt_half"):
            validate_config(config)

    def test_P_adapt_half_invalid_negative(self):
        """P_adapt_half must be > 0."""
        config = SimulationConfig()
        config.disease.pathogen_adaptation = True
        config.disease.P_adapt_half = -100.0
        with pytest.raises(ValueError, match="P_adapt_half"):
            validate_config(config)

    def test_pathogen_revert_rate_default_zero(self):
        """pathogen_revert_rate should default to 0.0."""
        cfg = DiseaseSection()
        assert cfg.pathogen_revert_rate == 0.0


# ═══════════════════════════════════════════════════════════════════════
# environmental_vibrio() with local threshold
# ═══════════════════════════════════════════════════════════════════════


class TestEnvironmentalVibrioLocalThreshold:
    """Verify that environmental_vibrio uses T_vbnc_local when provided."""

    def test_environmental_vibrio_uses_local_threshold(self):
        """With a lower T_vbnc_local, VBNC activation should be higher at same T."""
        cfg = DiseaseSection(T_vbnc=12.0, k_vbnc=1.0)

        # Without local override (uses global T_vbnc=12.0)
        env_global = environmental_vibrio(10.0, 30.0, cfg)

        # With local override at 8.0 (shifted sigmoid to the left)
        env_local = environmental_vibrio(10.0, 30.0, cfg, T_vbnc_local=8.0)

        # Lower T_vbnc_local → higher activation at T=10°C
        assert env_local > env_global

    def test_environmental_vibrio_none_uses_global(self):
        """When T_vbnc_local=None, should be identical to no override."""
        cfg = DiseaseSection(T_vbnc=12.0)
        env_none = environmental_vibrio(15.0, 30.0, cfg, T_vbnc_local=None)
        env_default = environmental_vibrio(15.0, 30.0, cfg)
        assert env_none == env_default

    def test_environmental_vibrio_local_equals_global_same_value(self):
        """When T_vbnc_local equals cfg.T_vbnc, results should be identical."""
        cfg = DiseaseSection(T_vbnc=12.0)
        env_global = environmental_vibrio(15.0, 30.0, cfg)
        env_local = environmental_vibrio(15.0, 30.0, cfg, T_vbnc_local=12.0)
        assert env_local == pytest.approx(env_global)


# ═══════════════════════════════════════════════════════════════════════
# Exact calculation verification
# ═══════════════════════════════════════════════════════════════════════


class TestAdaptationExactCalculation:
    """Verify exact delta_T calculation with pool-based formula."""

    def test_exact_delta_saturated_pool(self):
        cfg = _make_cfg(pathogen_adapt_rate=0.01, P_adapt_half=500.0)
        nds = _make_node_state(12.0)

        T_celsius = 8.0
        P_env_pool = 1000.0  # > P_adapt_half → pool_factor = 1.0

        adapt_pathogen_thermal(nds, T_celsius, P_env_pool=P_env_pool, cfg=cfg)

        # Expected: pool_factor = min(1000/500, 1.0) = 1.0
        # gap = 12-8 = 4, delta_T = 0.01 * 1.0 * 4 = 0.04
        expected = 12.0 - 0.04
        assert nds.T_vbnc_local == pytest.approx(expected)

    def test_exact_delta_partial_pool(self):
        cfg = _make_cfg(pathogen_adapt_rate=0.01, P_adapt_half=500.0)
        nds = _make_node_state(12.0)

        T_celsius = 8.0
        P_env_pool = 250.0  # half of P_adapt_half → pool_factor = 0.5

        adapt_pathogen_thermal(nds, T_celsius, P_env_pool=P_env_pool, cfg=cfg)

        # Expected: pool_factor = 250/500 = 0.5
        # gap = 12-8 = 4, delta_T = 0.01 * 0.5 * 4 = 0.02
        expected = 12.0 - 0.02
        assert nds.T_vbnc_local == pytest.approx(expected)

    def test_exact_delta_at_half_saturation(self):
        cfg = _make_cfg(pathogen_adapt_rate=0.01, P_adapt_half=500.0)
        nds = _make_node_state(12.0)

        T_celsius = 8.0
        P_env_pool = 500.0  # exactly P_adapt_half → pool_factor = 1.0

        adapt_pathogen_thermal(nds, T_celsius, P_env_pool=P_env_pool, cfg=cfg)

        # pool_factor = min(500/500, 1.0) = 1.0
        # delta_T = 0.01 * 1.0 * 4 = 0.04
        expected = 12.0 - 0.04
        assert nds.T_vbnc_local == pytest.approx(expected)


# ══════════════════════════════════════════════════════════════════════
# Chunk 2: Wavefront inheritance tests
# ══════════════════════════════════════════════════════════════════════

class TestInheritTVbnc:
    """Tests for _inherit_pathogen_traits helper function."""

    def test_inherit_basic(self):
        """Target node inherits weighted-average T_vbnc from reached sources."""
        from scipy.sparse import csr_matrix
        from sswd_evoepi.model import _inherit_pathogen_traits

        # 3 nodes: 0 and 1 reached, 2 is target
        nds_list = [NodeDiseaseState(node_id=i) for i in range(3)]
        nds_list[0].T_vbnc_local = 9.0   # adapted source
        nds_list[1].T_vbnc_local = 10.0  # less adapted source
        nds_list[2].T_vbnc_local = 12.0  # target (default)

        disease_reached = [True, True, False]
        P = np.array([100.0, 50.0, 0.0])  # vibrio concentrations

        # D^T matrix: row 2 (target) has weights from nodes 0 and 1
        # D^T[2,0] = 0.6, D^T[2,1] = 0.3
        D_T = csr_matrix(np.array([
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.6, 0.3, 0.0],
        ]))

        cfg = DiseaseSection(pathogen_adaptation=True, T_vbnc_initial=12.0)
        _inherit_pathogen_traits(2, nds_list, disease_reached, D_T, P, cfg)

        # Expected: weighted avg = (0.6*100*9.0 + 0.3*50*10.0) / (0.6*100 + 0.3*50)
        #         = (540 + 150) / (60 + 15) = 690 / 75 = 9.2
        assert nds_list[2].T_vbnc_local == pytest.approx(9.2)

    def test_inherit_no_reached_sources(self):
        """When no reached sources contribute, T_vbnc stays at default."""
        from scipy.sparse import csr_matrix
        from sswd_evoepi.model import _inherit_pathogen_traits

        nds_list = [NodeDiseaseState(node_id=i) for i in range(2)]
        nds_list[0].T_vbnc_local = 9.0
        nds_list[1].T_vbnc_local = 12.0  # target

        disease_reached = [False, False]  # nothing reached
        P = np.array([100.0, 0.0])

        D_T = csr_matrix(np.array([
            [0.0, 0.0],
            [0.5, 0.0],
        ]))

        cfg = DiseaseSection(pathogen_adaptation=True, T_vbnc_initial=12.0)
        _inherit_pathogen_traits(1, nds_list, disease_reached, D_T, P, cfg)

        assert nds_list[1].T_vbnc_local == 12.0  # unchanged

    def test_inherit_single_source(self):
        """With a single source, target inherits that source's T_vbnc."""
        from scipy.sparse import csr_matrix
        from sswd_evoepi.model import _inherit_pathogen_traits

        nds_list = [NodeDiseaseState(node_id=i) for i in range(2)]
        nds_list[0].T_vbnc_local = 8.5
        nds_list[1].T_vbnc_local = 12.0

        disease_reached = [True, False]
        P = np.array([200.0, 0.0])

        D_T = csr_matrix(np.array([
            [0.0, 0.0],
            [0.4, 0.0],
        ]))

        cfg = DiseaseSection(pathogen_adaptation=True, T_vbnc_initial=12.0)
        _inherit_pathogen_traits(1, nds_list, disease_reached, D_T, P, cfg)

        assert nds_list[1].T_vbnc_local == pytest.approx(8.5)


# ══════════════════════════════════════════════════════════════════════
# Chunk 3: Recording tests
# ══════════════════════════════════════════════════════════════════════

class TestYearlyTVbncRecording:
    """Tests for yearly_T_vbnc recording in SpatialSimResult."""

    def test_result_has_yearly_T_vbnc_field(self):
        """SpatialSimResult has the yearly_T_vbnc field."""
        from sswd_evoepi.model import SpatialSimResult
        r = SpatialSimResult()
        assert hasattr(r, 'yearly_T_vbnc')
        assert r.yearly_T_vbnc is None  # default None when disabled
