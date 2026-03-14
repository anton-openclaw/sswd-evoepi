"""Tests for seasonal salinity model.

Tests backward compatibility (fw_strength=0), spatial asymmetry
(AK vs CA), seasonal pulse, and integration with model.py.
"""

import math
import numpy as np
import pytest

from sswd_evoepi.salinity import (
    ocean_baseline,
    freshwater_melt_pulse,
    latitude_melt_factor,
    compute_salinity_array,
    _monthly_to_daily,
    _load_woa23,
    _PEAK_DAY,
    DAYS_PER_YEAR,
)
from sswd_evoepi.spatial import NodeDefinition


def _make_node(lat=48.0, lon=-123.0, fjord_depth_norm=0.0, name="test"):
    """Create a minimal NodeDefinition for salinity tests."""
    return NodeDefinition(
        node_id=0, name=name, lat=lat, lon=lon,
        subregion="TEST", habitat_area=10000.0,
        carrying_capacity=1000, fjord_depth_norm=fjord_depth_norm,
    )


# ── ocean_baseline ─────────────────────────────────────────────────

class TestOceanBaseline:
    def test_reference_latitude(self):
        """At 50°N, baseline should be exactly 31.32 psu."""
        assert ocean_baseline(50.0) == pytest.approx(31.32)

    def test_increases_with_latitude(self):
        """Higher latitude → slightly higher baseline (positive slope)."""
        assert ocean_baseline(60.0) > ocean_baseline(50.0)
        assert ocean_baseline(50.0) > ocean_baseline(40.0)

    def test_slope(self):
        """Slope should be 0.054 psu/degree."""
        diff = ocean_baseline(51.0) - ocean_baseline(50.0)
        assert diff == pytest.approx(0.054)


# ── freshwater_melt_pulse ──────────────────────────────────────────

class TestFreshwaterMeltPulse:
    def test_peak_at_june15(self):
        """Peak pulse should be 1.0 at day 166 (~June 15)."""
        assert freshwater_melt_pulse(_PEAK_DAY) == pytest.approx(1.0)

    def test_zero_in_winter(self):
        """Pulse should be 0 during winter (e.g., January 1)."""
        assert freshwater_melt_pulse(0) == 0.0

    def test_zero_in_december(self):
        """Pulse should be 0 in December."""
        assert freshwater_melt_pulse(350) == 0.0

    def test_symmetric_around_peak(self):
        """Pulse should be symmetric around peak day."""
        offset = 30
        before = freshwater_melt_pulse(_PEAK_DAY - offset)
        after = freshwater_melt_pulse(_PEAK_DAY + offset)
        assert before == pytest.approx(after, abs=1e-10)

    def test_non_negative(self):
        """Pulse should never be negative."""
        for d in range(DAYS_PER_YEAR):
            assert freshwater_melt_pulse(d) >= 0.0

    def test_bounded_zero_one(self):
        """Pulse should be in [0, 1]."""
        for d in range(DAYS_PER_YEAR):
            v = freshwater_melt_pulse(d)
            assert 0.0 <= v <= 1.0


# ── latitude_melt_factor ──────────────────────────────────────────

class TestLatitudeMeltFactor:
    def test_zero_at_35N(self):
        """No glacial melt at 35°N."""
        assert latitude_melt_factor(35.0) == 0.0

    def test_one_at_60N(self):
        """Full glacial influence at 60°N."""
        assert latitude_melt_factor(60.0) == 1.0

    def test_half_at_47_5N(self):
        """Half-way at 47.5°N."""
        assert latitude_melt_factor(47.5) == pytest.approx(0.5)

    def test_clamped_below_35(self):
        """Below 35°N should clamp to 0."""
        assert latitude_melt_factor(30.0) == 0.0

    def test_clamped_above_60(self):
        """Above 60°N should clamp to 1."""
        assert latitude_melt_factor(65.0) == 1.0


# ── compute_salinity_array ─────────────────────────────────────────

class TestComputeSalinityArray:
    def test_parametric_fallback_fw_zero(self):
        """fw_strength=0 with no WOA23 should return parametric baseline."""
        nodes = [_make_node(lat=50.0, fjord_depth_norm=0.5)]
        sal = compute_salinity_array(nodes, fw_strength=0.0,
                                     woa23_path='/nonexistent/path.npz')
        expected = ocean_baseline(50.0)
        assert sal.shape == (1, 365)
        np.testing.assert_allclose(sal[0, :], expected)

    def test_woa23_baseline_has_seasonal_variation(self):
        """fw_strength=0 with WOA23 should have seasonal variation."""
        woa = _load_woa23()
        if woa is None:
            pytest.skip("WOA23 data not available")
        nodes = [_make_node(lat=50.0, lon=-125.0, fjord_depth_norm=0.5)]
        sal = compute_salinity_array(nodes, fw_strength=0.0)
        assert sal.shape == (1, 365)
        sal_range = float(sal[0, :].max() - sal[0, :].min())
        # WOA23 should show some seasonal variation (>0.1 psu)
        assert sal_range > 0.1
        # But baseline should be in reasonable range
        assert sal[0, :].min() > 20.0
        assert sal[0, :].max() < 36.0

    def test_shape(self):
        """Output should be (N, 365)."""
        nodes = [_make_node(lat=lat) for lat in [45.0, 55.0, 60.0]]
        sal = compute_salinity_array(nodes, fw_strength=10.0)
        assert sal.shape == (3, 365)

    def test_dtype(self):
        """Output should be float32."""
        nodes = [_make_node()]
        sal = compute_salinity_array(nodes, fw_strength=0.0)
        assert sal.dtype == np.float32

    def test_no_depression_open_coast(self):
        """Open-coast node (fjord_depth_norm=0) should have no fjord depression."""
        nodes = [_make_node(lat=57.0, lon=-135.0, fjord_depth_norm=0.0)]
        sal_fw0 = compute_salinity_array(nodes, fw_strength=0.0)
        sal_fw20 = compute_salinity_array(nodes, fw_strength=20.0)
        # With fd_norm=0, fw_strength shouldn't matter
        np.testing.assert_allclose(sal_fw0[0, :], sal_fw20[0, :], atol=1e-5)

    def test_depression_at_fjord_peak(self):
        """Deep fjord at high latitude should show depression from baseline in June."""
        nodes = [_make_node(lat=58.0, lon=-135.0, fjord_depth_norm=0.9, name="AK-fjord")]
        sal_baseline = compute_salinity_array(nodes, fw_strength=0.0)
        sal_depressed = compute_salinity_array(nodes, fw_strength=15.0)
        june_base = float(sal_baseline[0, _PEAK_DAY])
        june_dep = float(sal_depressed[0, _PEAK_DAY])
        # Should be significantly depressed from baseline
        assert june_dep < june_base - 5.0
        # Winter should be near baseline (pulse ≈ 0)
        jan_base = float(sal_baseline[0, 0])
        jan_dep = float(sal_depressed[0, 0])
        assert jan_dep == pytest.approx(jan_base, abs=0.1)

    def test_no_depression_at_low_latitude(self):
        """CA node (36°N, fjord_depth_norm=0.5) should have near-zero depression."""
        nodes = [_make_node(lat=36.0, lon=-122.0, fjord_depth_norm=0.5, name="CA-site")]
        sal_baseline = compute_salinity_array(nodes, fw_strength=0.0)
        sal_depressed = compute_salinity_array(nodes, fw_strength=15.0)
        # latitude_melt_factor(36) = (36-35)/25 = 0.04, very small
        june_base = float(sal_baseline[0, _PEAK_DAY])
        june_dep = float(sal_depressed[0, _PEAK_DAY])
        assert june_base - june_dep < 1.0  # < 1 psu depression

    def test_ak_ca_asymmetry(self):
        """AK fjord should have much more depression than CA site in June."""
        ak = _make_node(lat=58.0, lon=-135.0, fjord_depth_norm=0.8, name="AK-fjord")
        ca = _make_node(lat=36.0, lon=-122.0, fjord_depth_norm=0.3, name="CA-coast")
        sal_base = compute_salinity_array([ak, ca], fw_strength=0.0)
        sal_dep = compute_salinity_array([ak, ca], fw_strength=15.0)

        ak_depression = float(sal_base[0, _PEAK_DAY] - sal_dep[0, _PEAK_DAY])
        ca_depression = float(sal_base[1, _PEAK_DAY] - sal_dep[1, _PEAK_DAY])

        # AK depression should be >>10× CA depression
        assert ak_depression > 5.0
        assert ca_depression < 1.0
        assert ak_depression > 10 * ca_depression

    def test_floor_enforced(self):
        """Salinity should never drop below s_floor."""
        nodes = [_make_node(lat=60.0, lon=-140.0, fjord_depth_norm=1.0)]
        sal = compute_salinity_array(nodes, fw_strength=50.0, s_floor=5.0)
        assert sal.min() >= 5.0

    def test_fw_depth_exp_sqrt(self):
        """fw_depth_exp=0.5 should amplify moderate fjord_depth_norm values."""
        node_moderate = _make_node(lat=55.0, lon=-130.0, fjord_depth_norm=0.25, name="moderate")
        sal_linear = compute_salinity_array([node_moderate], fw_strength=15.0, fw_depth_exp=1.0)
        sal_sqrt = compute_salinity_array([node_moderate], fw_strength=15.0, fw_depth_exp=0.5)

        # sqrt(0.25) = 0.5 > 0.25, so sqrt should give MORE depression
        june_linear = float(sal_linear[0, _PEAK_DAY])
        june_sqrt = float(sal_sqrt[0, _PEAK_DAY])
        assert june_sqrt < june_linear  # more depression with sqrt


class TestMonthlyToDaily:
    def test_shape(self):
        """Output should be (365,)."""
        monthly = np.array([30.0 + i * 0.1 for i in range(12)])
        daily = _monthly_to_daily(monthly)
        assert daily.shape == (365,)

    def test_preserves_mean(self):
        """Daily mean should be close to monthly mean."""
        monthly = np.array([30.0 + i * 0.2 for i in range(12)])
        daily = _monthly_to_daily(monthly)
        assert abs(float(daily.mean()) - monthly.mean()) < 0.15

    def test_smooth_interpolation(self):
        """Daily values should not have discontinuities."""
        monthly = np.array([32.0, 31.8, 31.5, 31.0, 30.0, 29.0,
                           28.5, 29.0, 30.0, 31.0, 31.5, 31.8])
        daily = _monthly_to_daily(monthly)
        # Max day-to-day change should be small (< 0.2 psu)
        diffs = np.abs(np.diff(daily))
        assert diffs.max() < 0.2


# ── Integration tests ─────────────────────────────────────────────

class TestSalinityIntegration:
    def test_config_fw_strength_default_zero(self):
        """Default config should have fw_strength=0 (mechanism OFF)."""
        from sswd_evoepi.config import default_config
        cfg = default_config()
        assert cfg.disease.fw_strength == 0.0

    def test_config_fw_depth_exp_default(self):
        """Default config should have fw_depth_exp=1.0."""
        from sswd_evoepi.config import default_config
        cfg = default_config()
        assert cfg.disease.fw_depth_exp == 1.0

    def test_model_imports_salinity(self):
        """model.py should import compute_salinity_array."""
        from sswd_evoepi.model import compute_salinity_array  # noqa: F401

    def test_node_definition_has_fjord_depth_norm(self):
        """NodeDefinition should have fjord_depth_norm field."""
        nd = _make_node(fjord_depth_norm=0.42)
        assert nd.fjord_depth_norm == 0.42

    def test_5node_network_fjord_depth_norm(self):
        """5-node test network should have fjord_depth_norm set."""
        from sswd_evoepi.spatial import get_5node_definitions
        defs = get_5node_definitions()
        # Howe Sound (node 1) should have non-zero fjord_depth_norm
        assert defs[1].fjord_depth_norm == 0.5
        # Sitka (node 0) should be 0 (open coast)
        assert defs[0].fjord_depth_norm == 0.0

    def test_calibration_runner_fw_strength_override(self):
        """calibration_runner apply_param_overrides should handle disease.fw_strength."""
        from sswd_evoepi.config import default_config
        from experiments.calibration_runner import apply_param_overrides
        cfg = default_config()
        apply_param_overrides(cfg, {'disease.fw_strength': 12.0})
        assert cfg.disease.fw_strength == 12.0

    def test_calibration_runner_fw_depth_exp_override(self):
        """calibration_runner apply_param_overrides should handle disease.fw_depth_exp."""
        from sswd_evoepi.config import default_config
        from experiments.calibration_runner import apply_param_overrides
        cfg = default_config()
        apply_param_overrides(cfg, {'disease.fw_depth_exp': 0.5})
        assert cfg.disease.fw_depth_exp == 0.5


# ── Quantitative validation (against design document) ─────────────

class TestSalinityPredictions:
    """Test predictions from salinity implementation sketch analysis."""

    def test_ak_fn_june_depression(self):
        """AK-FN-like site at fw_strength=15: significant depression in June."""
        # AK-FN median fjord_depth_norm ≈ 0.77, lat ≈ 57.5
        node = _make_node(lat=57.5, lon=-135.0, fjord_depth_norm=0.77, name="AK-FN")
        sal_base = compute_salinity_array([node], fw_strength=0.0)
        sal_dep = compute_salinity_array([node], fw_strength=15.0)
        june_base = float(sal_base[0, _PEAK_DAY])
        june_dep = float(sal_dep[0, _PEAK_DAY])
        depression = june_base - june_dep
        # Expected: fw_strength(15) × fd_norm(0.77) × f_melt(0.9) × pulse(1.0) ≈ 10.4
        assert depression > 8.0
        assert depression < 13.0

    def test_or_no_depression(self):
        """OR coast (lat=44, fjord_depth_norm≈0.0): near-zero depression."""
        node = _make_node(lat=44.0, lon=-124.0, fjord_depth_norm=0.0, name="OR-coast")
        sal_base = compute_salinity_array([node], fw_strength=0.0)
        sal_dep = compute_salinity_array([node], fw_strength=15.0)
        june_base = float(sal_base[0, _PEAK_DAY])
        june_dep = float(sal_dep[0, _PEAK_DAY])
        assert abs(june_base - june_dep) < 0.01
