"""Tests for the seasonal salinity model (sswd_evoepi.salinity)."""

import math
from types import SimpleNamespace

import numpy as np
import pytest

from sswd_evoepi.salinity import (
    compute_salinity_array,
    freshwater_melt_pulse,
    latitude_melt_factor,
    ocean_baseline,
    DAYS_PER_YEAR,
)


# ─── ocean_baseline ─────────────────────────────────────────────────

class TestOceanBaseline:
    def test_at_reference_lat(self):
        """At 50°N the offset term is zero → exactly 31.32."""
        assert ocean_baseline(50.0) == pytest.approx(31.32)

    def test_higher_lat(self):
        """At 60°N: 31.32 + 0.054×10 = 31.86."""
        assert ocean_baseline(60.0) == pytest.approx(31.86)

    def test_lower_lat(self):
        """At 40°N: 31.32 + 0.054×(-10) = 30.78."""
        assert ocean_baseline(40.0) == pytest.approx(30.78)


# ─── freshwater_melt_pulse ──────────────────────────────────────────

class TestFreshwaterMeltPulse:
    def test_peak_day(self):
        """Pulse is 1.0 on peak day (166 ≈ June 15)."""
        assert freshwater_melt_pulse(166) == pytest.approx(1.0)

    def test_zero_in_winter(self):
        """Pulse is 0 near day 350 (mid-December)."""
        assert freshwater_melt_pulse(350) == pytest.approx(0.0)

    def test_non_negative(self):
        """Pulse is never negative for any day of year."""
        for d in range(DAYS_PER_YEAR):
            assert freshwater_melt_pulse(d) >= 0.0

    def test_symmetry_around_peak(self):
        """Symmetric around peak day within the positive lobe."""
        val_before = freshwater_melt_pulse(166 - 30)
        val_after = freshwater_melt_pulse(166 + 30)
        assert val_before == pytest.approx(val_after, abs=1e-10)


# ─── latitude_melt_factor ───────────────────────────────────────────

class TestLatitudeMeltFactor:
    def test_zero_at_35(self):
        assert latitude_melt_factor(35.0) == pytest.approx(0.0)

    def test_one_at_60(self):
        assert latitude_melt_factor(60.0) == pytest.approx(1.0)

    def test_half_at_47_5(self):
        assert latitude_melt_factor(47.5) == pytest.approx(0.5)

    def test_clipped_below(self):
        """Below 35°N, factor should be 0."""
        assert latitude_melt_factor(30.0) == 0.0

    def test_clipped_above(self):
        """Above 60°N, factor should be 1."""
        assert latitude_melt_factor(70.0) == 1.0


# ─── compute_salinity_array ─────────────────────────────────────────

def _make_node(lat: float, fjord_depth_norm: float = 0.0):
    """Create a minimal node-like object for testing."""
    return SimpleNamespace(lat=lat, fjord_depth_norm=fjord_depth_norm)


class TestComputeSalinityArray:
    def test_fw_strength_zero_gives_constant_baseline(self):
        """With fw_strength=0, salinity is ocean baseline for every day."""
        nodes = [_make_node(50.0, 0.5), _make_node(45.0, 0.8)]
        result = compute_salinity_array(nodes, fw_strength=0.0)
        assert result.shape == (2, 365)
        # Every day should be the ocean baseline
        for i, nd in enumerate(nodes):
            expected = ocean_baseline(nd.lat)
            np.testing.assert_allclose(result[i, :], expected, atol=1e-5)

    def test_shape(self):
        nodes = [_make_node(50.0)]
        result = compute_salinity_array(nodes, fw_strength=5.0)
        assert result.shape == (1, 365)
        assert result.dtype == np.float32

    def test_summer_depression_high_lat_fjord(self):
        """High-lat fjord with fw_strength=15 should show summer depression."""
        nodes = [_make_node(lat=58.0, fjord_depth_norm=0.9)]
        result = compute_salinity_array(nodes, fw_strength=15.0)
        baseline = ocean_baseline(58.0)
        # Day 166 is peak depression
        summer_val = float(result[0, 166])
        winter_val = float(result[0, 0])
        # Summer should be depressed below baseline
        assert summer_val < baseline - 1.0
        # Winter should be at baseline (pulse is zero)
        assert winter_val == pytest.approx(baseline, abs=0.1)

    def test_open_coast_no_depression(self):
        """Open coast (fjord_depth_norm=0) should show no depression."""
        nodes = [_make_node(lat=58.0, fjord_depth_norm=0.0)]
        result = compute_salinity_array(nodes, fw_strength=15.0)
        baseline = ocean_baseline(58.0)
        np.testing.assert_allclose(result[0, :], baseline, atol=1e-5)

    def test_low_lat_no_depression(self):
        """Low latitude (35°N) has zero melt factor → no depression."""
        nodes = [_make_node(lat=35.0, fjord_depth_norm=1.0)]
        result = compute_salinity_array(nodes, fw_strength=15.0)
        baseline = ocean_baseline(35.0)
        np.testing.assert_allclose(result[0, :], baseline, atol=1e-5)

    def test_salinity_floor(self):
        """Extreme parameters should not push salinity below s_floor."""
        nodes = [_make_node(lat=60.0, fjord_depth_norm=1.0)]
        result = compute_salinity_array(nodes, fw_strength=100.0, s_floor=5.0)
        assert result.min() >= 5.0


# ─── Backward compatibility ─────────────────────────────────────────

class TestBackwardCompatibility:
    def test_fw_strength_zero_matches_static_salinity(self):
        """fw_strength=0 → precomputed salinity equals ocean_baseline for each node.

        This ensures the default config produces EXACTLY the same salinity
        as the old hardcoded nd.salinity path (which used constant values).
        """
        nodes = [
            _make_node(lat=57.06, fjord_depth_norm=0.0),
            _make_node(lat=49.52, fjord_depth_norm=0.5),
            _make_node(lat=36.62, fjord_depth_norm=0.0),
        ]
        result = compute_salinity_array(nodes, fw_strength=0.0)
        for i, nd in enumerate(nodes):
            expected = ocean_baseline(nd.lat)
            # Every single day must be exactly the baseline
            assert np.all(result[i, :] == pytest.approx(expected, abs=1e-5))
