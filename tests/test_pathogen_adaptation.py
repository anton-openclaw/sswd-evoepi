"""Tests for pathogen thermal adaptation (Chunk 1).

Tests cover:
  - adapt_pathogen_thermal() function behavior
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
        pathogen_revert_rate=0.0005,
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
        adapt_pathogen_thermal(nds, T_celsius=8.0, n_I1=10, n_I2=5,
                               n_E=3, n_total=100, cfg=cfg)
        assert nds.T_vbnc_local == 12.0

    def test_adaptation_disabled_explicit(self):
        cfg = _make_cfg(pathogen_adaptation=False)
        nds = _make_node_state(12.0)
        adapt_pathogen_thermal(nds, T_celsius=5.0, n_I1=50, n_I2=50,
                               n_E=10, n_total=100, cfg=cfg)
        assert nds.T_vbnc_local == 12.0


class TestAdaptationReducesThreshold:
    """When temp < T_vbnc_local and infected hosts present, threshold should decrease."""

    def test_adaptation_reduces_threshold(self):
        cfg = _make_cfg(pathogen_adapt_rate=0.01)
        nds = _make_node_state(12.0)
        # T=8°C < 12°C, 10 infected out of 100 → should adapt
        adapt_pathogen_thermal(nds, T_celsius=8.0, n_I1=5, n_I2=5,
                               n_E=2, n_total=100, cfg=cfg)
        assert nds.T_vbnc_local < 12.0

    def test_threshold_decreases_with_I1_only(self):
        cfg = _make_cfg(pathogen_adapt_rate=0.01)
        nds = _make_node_state(12.0)
        adapt_pathogen_thermal(nds, T_celsius=8.0, n_I1=10, n_I2=0,
                               n_E=0, n_total=100, cfg=cfg)
        assert nds.T_vbnc_local < 12.0

    def test_threshold_decreases_with_I2_only(self):
        cfg = _make_cfg(pathogen_adapt_rate=0.01)
        nds = _make_node_state(12.0)
        adapt_pathogen_thermal(nds, T_celsius=8.0, n_I1=0, n_I2=10,
                               n_E=0, n_total=100, cfg=cfg)
        assert nds.T_vbnc_local < 12.0


class TestNoAdaptationAboveThreshold:
    """When temp >= T_vbnc_local, no adaptation should occur."""

    def test_no_adaptation_above_threshold(self):
        cfg = _make_cfg()
        nds = _make_node_state(12.0)
        adapt_pathogen_thermal(nds, T_celsius=12.0, n_I1=10, n_I2=5,
                               n_E=3, n_total=100, cfg=cfg)
        assert nds.T_vbnc_local == 12.0

    def test_no_adaptation_well_above_threshold(self):
        cfg = _make_cfg()
        nds = _make_node_state(12.0)
        adapt_pathogen_thermal(nds, T_celsius=20.0, n_I1=50, n_I2=50,
                               n_E=10, n_total=100, cfg=cfg)
        assert nds.T_vbnc_local == 12.0


class TestNoAdaptationNoInfected:
    """When no infected hosts, adaptation branch should NOT activate."""

    def test_no_adaptation_no_infected(self):
        cfg = _make_cfg()
        nds = _make_node_state(12.0)
        adapt_pathogen_thermal(nds, T_celsius=5.0, n_I1=0, n_I2=0,
                               n_E=5, n_total=100, cfg=cfg)
        # No infected (I1+I2=0), even though temp is below threshold
        # and there are exposed — no adaptation
        assert nds.T_vbnc_local == 12.0

    def test_no_adaptation_zero_population(self):
        cfg = _make_cfg()
        nds = _make_node_state(12.0)
        adapt_pathogen_thermal(nds, T_celsius=5.0, n_I1=0, n_I2=0,
                               n_E=0, n_total=0, cfg=cfg)
        assert nds.T_vbnc_local == 12.0


class TestAdaptationRespectsMinimum:
    """T_vbnc_local should never go below T_vbnc_min."""

    def test_adaptation_respects_minimum(self):
        cfg = _make_cfg(
            pathogen_adapt_rate=10.0,  # Very aggressive rate
            T_vbnc_min=8.0,
        )
        nds = _make_node_state(12.0)
        # Extreme conditions: high rate, many infected, big temp gap
        adapt_pathogen_thermal(nds, T_celsius=0.0, n_I1=100, n_I2=100,
                               n_E=0, n_total=200, cfg=cfg)
        assert nds.T_vbnc_local == 8.0

    def test_adaptation_minimum_reached_incrementally(self):
        cfg = _make_cfg(
            pathogen_adapt_rate=1.0,
            T_vbnc_min=10.0,
        )
        nds = _make_node_state(12.0)
        for _ in range(1000):
            adapt_pathogen_thermal(nds, T_celsius=5.0, n_I1=50, n_I2=50,
                                   n_E=0, n_total=100, cfg=cfg)
        assert nds.T_vbnc_local == pytest.approx(10.0)


class TestReversionWithoutDisease:
    """When no infected/exposed hosts, T_vbnc_local should revert toward T_vbnc_initial."""

    def test_reversion_without_disease(self):
        cfg = _make_cfg(pathogen_revert_rate=0.1)
        nds = _make_node_state(10.0)  # already adapted (below initial=12)
        adapt_pathogen_thermal(nds, T_celsius=8.0, n_I1=0, n_I2=0,
                               n_E=0, n_total=100, cfg=cfg)
        assert nds.T_vbnc_local == pytest.approx(10.1)

    def test_reversion_capped_at_initial(self):
        cfg = _make_cfg(pathogen_revert_rate=5.0, T_vbnc_initial=12.0)
        nds = _make_node_state(11.5)
        adapt_pathogen_thermal(nds, T_celsius=8.0, n_I1=0, n_I2=0,
                               n_E=0, n_total=100, cfg=cfg)
        assert nds.T_vbnc_local == 12.0

    def test_no_reversion_when_exposed_present(self):
        """Reversion only happens when BOTH infected=0 AND exposed=0."""
        cfg = _make_cfg(pathogen_revert_rate=0.1)
        nds = _make_node_state(10.0)
        adapt_pathogen_thermal(nds, T_celsius=8.0, n_I1=0, n_I2=0,
                               n_E=5, n_total=100, cfg=cfg)
        # E>0, so no reversion
        assert nds.T_vbnc_local == 10.0

    def test_no_reversion_at_initial(self):
        """No reversion when already at T_vbnc_initial."""
        cfg = _make_cfg(pathogen_revert_rate=0.1, T_vbnc_initial=12.0)
        nds = _make_node_state(12.0)
        adapt_pathogen_thermal(nds, T_celsius=8.0, n_I1=0, n_I2=0,
                               n_E=0, n_total=100, cfg=cfg)
        assert nds.T_vbnc_local == 12.0


class TestAdaptationScaling:
    """Adaptation rate should scale with prevalence and temperature gap."""

    def test_adaptation_rate_scales_with_prevalence(self):
        cfg = _make_cfg(pathogen_adapt_rate=0.01)

        # Low prevalence
        nds_low = _make_node_state(12.0)
        adapt_pathogen_thermal(nds_low, T_celsius=8.0, n_I1=5, n_I2=5,
                               n_E=0, n_total=1000, cfg=cfg)
        delta_low = 12.0 - nds_low.T_vbnc_local

        # High prevalence
        nds_high = _make_node_state(12.0)
        adapt_pathogen_thermal(nds_high, T_celsius=8.0, n_I1=50, n_I2=50,
                               n_E=0, n_total=200, cfg=cfg)
        delta_high = 12.0 - nds_high.T_vbnc_local

        assert delta_high > delta_low
        # Verify proportionality: prevalence ratio = (100/200) / (10/1000) = 50x
        ratio = delta_high / delta_low
        assert ratio == pytest.approx(50.0, rel=0.01)

    def test_adaptation_rate_scales_with_temp_gap(self):
        cfg = _make_cfg(pathogen_adapt_rate=0.01)

        # Small gap (11°C vs 12°C threshold → gap=1)
        nds_small = _make_node_state(12.0)
        adapt_pathogen_thermal(nds_small, T_celsius=11.0, n_I1=50, n_I2=0,
                               n_E=0, n_total=100, cfg=cfg)
        delta_small = 12.0 - nds_small.T_vbnc_local

        # Large gap (4°C vs 12°C threshold → gap=8)
        nds_large = _make_node_state(12.0)
        adapt_pathogen_thermal(nds_large, T_celsius=4.0, n_I1=50, n_I2=0,
                               n_E=0, n_total=100, cfg=cfg)
        delta_large = 12.0 - nds_large.T_vbnc_local

        assert delta_large > delta_small
        # Verify proportionality: gap ratio = 8/1 = 8x
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
        config.disease.pathogen_revert_rate = 0.0005
        validate_config(config)  # should not raise

    def test_config_validation_skipped_when_disabled(self):
        """When pathogen_adaptation=False, invalid T_vbnc_min is OK."""
        config = SimulationConfig()
        config.disease.pathogen_adaptation = False
        config.disease.T_vbnc_min = 15.0  # invalid but ignored
        config.disease.T_vbnc_initial = 12.0
        validate_config(config)  # should not raise


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
    """Verify exact delta_T calculation."""

    def test_exact_delta(self):
        cfg = _make_cfg(pathogen_adapt_rate=0.01)
        nds = _make_node_state(12.0)

        T_celsius = 8.0
        n_I1, n_I2, n_total = 10, 10, 100

        adapt_pathogen_thermal(nds, T_celsius, n_I1, n_I2,
                               n_E=0, n_total=n_total, cfg=cfg)

        # Expected: prevalence = 20/100 = 0.2, gap = 12-8 = 4
        # delta_T = 0.01 * 0.2 * 4 = 0.008
        expected = 12.0 - 0.008
        assert nds.T_vbnc_local == pytest.approx(expected)
