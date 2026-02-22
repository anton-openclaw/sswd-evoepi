"""Tests for satellite SST integration.

Tests cover:
  - Loading satellite SST climatology from CSV files
  - Name normalization for node lookups
  - Time series generation from climatology
  - Warming trend application
  - Config validation for sst_source field
  - Backward compatibility (default = sinusoidal)
  - Comparison of satellite vs sinusoidal SST ranges

Phase 6 of PLAN-feb21-fixes.md.
"""

from __future__ import annotations

import csv
import os
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pytest

from sswd_evoepi.environment import (
    _normalize_node_name,
    generate_satellite_sst_series,
    load_sst_climatology,
    make_sst_timeseries,
)
from sswd_evoepi.config import (
    SimulationConfig,
    SimulationSection,
    default_config,
    validate_config,
)


# ── Fixtures ──────────────────────────────────────────────────────────

@pytest.fixture
def tmp_sst_dir(tmp_path):
    """Create a temporary directory with a valid climatology CSV."""
    # Create a realistic-ish sinusoidal climatology for a fake node
    days = np.arange(1, 366)
    sst_mean = 10.0 + 3.5 * np.cos(2 * np.pi * (days - 227) / 365)
    sst_std = np.full(365, 0.8)
    n_years = np.full(365, 24, dtype=int)

    csv_path = tmp_path / "Test_Node_climatology.csv"
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['day_of_year', 'sst_mean', 'sst_std', 'n_years'])
        for i in range(365):
            writer.writerow([days[i], f"{sst_mean[i]:.4f}", f"{sst_std[i]:.4f}", n_years[i]])

    return tmp_path


@pytest.fixture
def real_sst_dir():
    """Path to the real SST data directory, if it exists."""
    repo_root = Path(__file__).parent.parent
    sst_dir = repo_root / "data" / "sst"
    if not sst_dir.is_dir():
        pytest.skip("data/sst/ not found — satellite SST data not available")
    return sst_dir


# ── Node name normalization ──────────────────────────────────────────

class TestNodeNameNormalization:
    def test_spaces_to_underscores(self):
        assert _normalize_node_name("Howe Sound") == "Howe_Sound"

    def test_no_change_needed(self):
        assert _normalize_node_name("SJI") == "SJI"

    def test_strips_whitespace(self):
        assert _normalize_node_name("  Fort Bragg  ") == "Fort_Bragg"

    def test_multiple_spaces(self):
        assert _normalize_node_name("Howe  Sound") == "Howe_Sound"

    def test_already_underscored(self):
        assert _normalize_node_name("Howe_Sound") == "Howe_Sound"

    def test_preserves_case(self):
        assert _normalize_node_name("howe sound") == "howe_sound"
        assert _normalize_node_name("HOWE SOUND") == "HOWE_SOUND"


# ── load_sst_climatology ─────────────────────────────────────────────

class TestLoadSSTClimatology:
    def test_load_valid_climatology(self, tmp_sst_dir):
        """Load a valid 365-row climatology CSV and verify shape."""
        clim = load_sst_climatology("Test Node", data_dir=str(tmp_sst_dir))
        assert clim.shape == (365,)
        assert clim.dtype == np.float64

    def test_load_climatology_values_reasonable(self, tmp_sst_dir):
        """SST values should be in a reasonable range."""
        clim = load_sst_climatology("Test Node", data_dir=str(tmp_sst_dir))
        assert np.all(clim > -5.0), "SST values too low"
        assert np.all(clim < 35.0), "SST values too high"

    def test_missing_node_raises_file_not_found(self, tmp_sst_dir):
        """Nonexistent node should raise FileNotFoundError."""
        with pytest.raises(FileNotFoundError, match="SST climatology file not found"):
            load_sst_climatology("Nonexistent_Node", data_dir=str(tmp_sst_dir))

    def test_name_normalization_spaces(self, tmp_sst_dir):
        """'Test Node' and 'Test_Node' should both resolve."""
        clim_space = load_sst_climatology("Test Node", data_dir=str(tmp_sst_dir))
        clim_underscore = load_sst_climatology("Test_Node", data_dir=str(tmp_sst_dir))
        np.testing.assert_array_equal(clim_space, clim_underscore)

    def test_wrong_row_count_raises_value_error(self, tmp_path):
        """A CSV with != 365 rows should raise ValueError."""
        csv_path = tmp_path / "Bad_Node_climatology.csv"
        with open(csv_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['day_of_year', 'sst_mean', 'sst_std', 'n_years'])
            for i in range(100):  # Only 100 rows
                writer.writerow([i + 1, 10.0, 0.5, 20])

        with pytest.raises(ValueError, match="100 rows.*expected 365"):
            load_sst_climatology("Bad Node", data_dir=str(tmp_path))

    def test_load_real_sitka(self, real_sst_dir):
        """Load actual Sitka climatology from data/sst/."""
        clim = load_sst_climatology("Sitka", data_dir=str(real_sst_dir))
        assert clim.shape == (365,)
        # Sitka should be cold — roughly 5-14°C
        assert np.min(clim) > 3.0
        assert np.max(clim) < 18.0

    def test_load_all_11_nodes(self, real_sst_dir):
        """All 11 stepping-stone nodes should load successfully."""
        nodes = [
            "Sitka", "Ketchikan", "Haida Gwaii", "Bella Bella",
            "Howe Sound", "SJI", "Westport", "Newport",
            "Crescent City", "Fort Bragg", "Monterey",
        ]
        for node in nodes:
            clim = load_sst_climatology(node, data_dir=str(real_sst_dir))
            assert clim.shape == (365,), f"Wrong shape for {node}"
            assert np.all(np.isfinite(clim)), f"NaN/Inf in {node} climatology"


# ── generate_satellite_sst_series ─────────────────────────────────────

class TestGenerateSatelliteSSTSeries:
    def test_output_shape(self, tmp_sst_dir):
        """Output shape should be (n_years * 365,)."""
        n_years = 10
        sst = generate_satellite_sst_series(
            n_years=n_years, start_year=2013,
            node_name="Test Node", data_dir=str(tmp_sst_dir),
        )
        assert sst.shape == (n_years * 365,)

    def test_output_shape_single_year(self, tmp_sst_dir):
        """Single year should produce exactly 365 values."""
        sst = generate_satellite_sst_series(
            n_years=1, start_year=2020,
            node_name="Test_Node", data_dir=str(tmp_sst_dir),
        )
        assert sst.shape == (365,)

    def test_no_trend_tiles_identically(self, tmp_sst_dir):
        """With zero trend, each year should be identical to the climatology."""
        clim = load_sst_climatology("Test Node", data_dir=str(tmp_sst_dir))
        sst = generate_satellite_sst_series(
            n_years=5, start_year=2015,
            node_name="Test Node", trend_per_year=0.0,
            reference_year=2015, data_dir=str(tmp_sst_dir),
        )
        for yr in range(5):
            yr_slice = sst[yr * 365:(yr + 1) * 365]
            np.testing.assert_array_almost_equal(yr_slice, clim)

    def test_warming_trend_applied(self, tmp_sst_dir):
        """With a warming trend, later years should be warmer."""
        trend = 0.05  # °C/year
        sst = generate_satellite_sst_series(
            n_years=20, start_year=2000,
            node_name="Test Node", trend_per_year=trend,
            reference_year=2000, data_dir=str(tmp_sst_dir),
        )
        # Mean of first year vs last year
        mean_first = np.mean(sst[:365])
        mean_last = np.mean(sst[-365:])
        expected_diff = trend * 19  # 19 years apart
        actual_diff = mean_last - mean_first
        np.testing.assert_almost_equal(actual_diff, expected_diff, decimal=5)

    def test_trend_direction(self, tmp_sst_dir):
        """Positive trend → later years warmer; negative → cooler."""
        sst_warm = generate_satellite_sst_series(
            n_years=10, start_year=2010, node_name="Test Node",
            trend_per_year=0.1, reference_year=2010,
            data_dir=str(tmp_sst_dir),
        )
        sst_cool = generate_satellite_sst_series(
            n_years=10, start_year=2010, node_name="Test Node",
            trend_per_year=-0.1, reference_year=2010,
            data_dir=str(tmp_sst_dir),
        )
        assert np.mean(sst_warm[-365:]) > np.mean(sst_warm[:365])
        assert np.mean(sst_cool[-365:]) < np.mean(sst_cool[:365])

    def test_matches_sinusoidal_shape(self, tmp_sst_dir):
        """Satellite SST series should have the same shape as sinusoidal."""
        n_years = 5
        sat_sst = generate_satellite_sst_series(
            n_years=n_years, start_year=2013,
            node_name="Test Node", data_dir=str(tmp_sst_dir),
        )
        sin_sst = make_sst_timeseries(
            n_years=n_years, start_year=2013,
            mean_sst=10.0, amplitude=3.5,
        )
        assert sat_sst.shape == sin_sst.shape


# ── Satellite vs sinusoidal range comparison ──────────────────────────

class TestSatelliteVsSinusoidalRange:
    def test_satellite_sst_range_reasonable(self, real_sst_dir):
        """Satellite SST for each node should produce reasonable ranges."""
        nodes = [
            "Sitka", "Ketchikan", "Haida Gwaii", "Bella Bella",
            "Howe Sound", "SJI", "Westport", "Newport",
            "Crescent City", "Fort Bragg", "Monterey",
        ]
        for node in nodes:
            sst = generate_satellite_sst_series(
                n_years=5, start_year=2013, node_name=node,
                trend_per_year=0.0, data_dir=str(real_sst_dir),
            )
            sst_range = np.max(sst) - np.min(sst)
            assert sst_range > 1.0, f"{node}: SST range too narrow ({sst_range:.2f}°C)"
            assert sst_range < 20.0, f"{node}: SST range too wide ({sst_range:.2f}°C)"
            assert np.mean(sst) > 5.0, f"{node}: mean SST too low ({np.mean(sst):.2f}°C)"
            assert np.mean(sst) < 20.0, f"{node}: mean SST too high ({np.mean(sst):.2f}°C)"

    def test_sinusoidal_range_similar_to_satellite(self, real_sst_dir):
        """Sinusoidal and satellite SST should produce similar magnitude ranges.

        They won't match exactly (satellite has daily variability, sinusoidal
        is smooth), but both should be in the same ballpark (within ~5°C).
        """
        clim = load_sst_climatology("SJI", data_dir=str(real_sst_dir))
        sat_range = np.max(clim) - np.min(clim)

        # SJI sinusoidal: mean ~10.3°C, amplitude ~3°C → range ~6°C
        sin_sst = make_sst_timeseries(
            n_years=1, start_year=2013, mean_sst=10.3, amplitude=3.0,
        )
        sin_range = np.max(sin_sst) - np.min(sin_sst)

        # Both ranges should be within 5°C of each other
        assert abs(sat_range - sin_range) < 5.0, (
            f"Satellite range ({sat_range:.2f}) and sinusoidal range "
            f"({sin_range:.2f}) differ by more than 5°C"
        )


# ── Config SST source validation ─────────────────────────────────────

class TestConfigSSTSource:
    def test_default_is_sinusoidal(self):
        """Default config should use sinusoidal SST (backward compatible)."""
        config = default_config()
        assert config.simulation.sst_source == "sinusoidal"

    def test_sinusoidal_accepted(self):
        """'sinusoidal' should pass validation."""
        config = default_config()
        config.simulation.sst_source = "sinusoidal"
        validate_config(config)  # Should not raise

    def test_satellite_accepted(self):
        """'satellite' should pass validation."""
        config = default_config()
        config.simulation.sst_source = "satellite"
        # May warn about missing dir, but should not raise ValueError
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            validate_config(config)

    def test_bogus_source_rejected(self):
        """Invalid sst_source should raise ValueError."""
        config = default_config()
        config.simulation.sst_source = "bogus"
        with pytest.raises(ValueError, match="sst_source"):
            validate_config(config)

    def test_satellite_missing_dir_warns(self):
        """satellite source with missing data dir should warn."""
        config = default_config()
        config.simulation.sst_source = "satellite"
        config.simulation.sst_data_dir = "/nonexistent/path/to/sst"
        with pytest.warns(UserWarning, match="does not exist"):
            validate_config(config)

    def test_default_sst_data_dir(self):
        """Default sst_data_dir should be 'data/sst'."""
        config = default_config()
        assert config.simulation.sst_data_dir == "data/sst"
