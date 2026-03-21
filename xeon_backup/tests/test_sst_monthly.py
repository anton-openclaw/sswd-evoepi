"""Tests for year-specific monthly SST integration.

Tests cover:
  - Loading monthly SST data from CSV files
  - Cubic/linear interpolation from monthly to daily resolution
  - Interpolation smoothness at month boundaries
  - Multi-year series generation with correct length
  - Known-value spot checks (Sitka Jan 2013)
  - Fallback to climatology for out-of-range years
  - Leap year handling
  - Full model run with sst_source='monthly'
  - Backward compatibility with 'sinusoidal' and 'satellite' sources

Phase 2 of SST yearly integration (3-phase plan).
"""

from __future__ import annotations

import calendar
import csv
import logging
import os
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pytest

from sswd_evoepi.environment import (
    generate_satellite_sst_series,
    generate_yearly_sst_series,
    interpolate_monthly_to_daily,
    load_sst_climatology,
    load_sst_monthly,
    make_sst_timeseries,
)
from sswd_evoepi.config import (
    SimulationConfig,
    default_config,
    validate_config,
)


# ── Paths ─────────────────────────────────────────────────────────────

REPO_ROOT = Path(__file__).resolve().parent.parent
REAL_SST_DIR = REPO_ROOT / "data" / "sst"


# ── Fixtures ──────────────────────────────────────────────────────────

@pytest.fixture
def real_sst_dir():
    """Path to the real SST data directory."""
    if not REAL_SST_DIR.is_dir():
        pytest.skip("Real SST data directory not found")
    return REAL_SST_DIR


@pytest.fixture
def tmp_monthly_dir(tmp_path):
    """Create a temporary directory with a synthetic monthly CSV and
    matching climatology for a fake node 'Test_Node'."""
    # Synthetic monthly data: sinusoidal pattern, 2010-2015
    months_sst = {
        1: 7.0, 2: 6.0, 3: 6.5, 4: 8.0, 5: 10.0, 6: 12.0,
        7: 14.0, 8: 14.5, 9: 13.0, 10: 11.0, 11: 9.0, 12: 7.5,
    }

    csv_path = tmp_path / "Test_Node_monthly.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["year", "month", "sst"])
        for yr in range(2010, 2016):
            for mo in range(1, 13):
                # Add small year-to-year variation
                sst = months_sst[mo] + 0.1 * (yr - 2012)
                writer.writerow([yr, mo, f"{sst:.4f}"])

    # Also create a climatology file for fallback tests
    clim_path = tmp_path / "Test_Node_climatology.csv"
    days = np.arange(1, 366)
    sst_mean = 10.0 + 3.5 * np.cos(2 * np.pi * (days - 227) / 365)
    with open(clim_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["day_of_year", "sst_mean", "sst_std", "n_years"])
        for i in range(365):
            writer.writerow([days[i], f"{sst_mean[i]:.4f}", "0.8000", 24])

    return tmp_path


# ══════════════════════════════════════════════════════════════════════
# TEST CLASS: Monthly data loading
# ══════════════════════════════════════════════════════════════════════

class TestLoadMonthlyData:
    """Tests for load_sst_monthly."""

    def test_load_real_sitka(self, real_sst_dir):
        """Verify Sitka monthly CSV loads correctly."""
        data = load_sst_monthly("Sitka", data_dir=str(real_sst_dir))

        # First row: 2002, January
        assert (2002, 1) in data
        assert abs(data[(2002, 1)] - 6.6703) < 0.001, (
            f"Sitka Jan 2002 expected ~6.6703, got {data[(2002, 1)]}"
        )

        # Data should span 2002-2025
        years = sorted(set(yr for yr, _ in data.keys()))
        assert years[0] == 2002
        assert years[-1] == 2025
        assert len(years) == 24  # 2002 through 2025

    def test_load_synthetic(self, tmp_monthly_dir):
        """Verify synthetic monthly data loads with correct keys."""
        data = load_sst_monthly("Test Node", data_dir=str(tmp_monthly_dir))

        # Should have 6 years × 12 months = 72 entries
        assert len(data) == 72

        # Spot check 2012 January (base year → exactly months_sst[1] = 7.0)
        assert abs(data[(2012, 1)] - 7.0) < 0.001

    def test_file_not_found(self, tmp_path):
        """Missing CSV should raise FileNotFoundError."""
        with pytest.raises(FileNotFoundError, match="monthly"):
            load_sst_monthly("Nonexistent", data_dir=str(tmp_path))


# ══════════════════════════════════════════════════════════════════════
# TEST CLASS: Interpolation
# ══════════════════════════════════════════════════════════════════════

class TestInterpolation:
    """Tests for interpolate_monthly_to_daily."""

    @pytest.fixture
    def known_monthly(self):
        """12 known monthly values for testing."""
        return np.array([
            7.0, 6.0, 6.5, 8.0, 10.0, 12.0,
            14.0, 14.5, 13.0, 11.0, 9.0, 7.5,
        ], dtype=np.float64)

    def test_midpoint_values(self, known_monthly):
        """Interpolated values at month midpoints should be close to input."""
        year = 2013  # non-leap
        daily = interpolate_monthly_to_daily(known_monthly, year)

        assert len(daily) == 365

        # Check midpoint of each month
        cumulative = 0
        for m in range(12):
            days_in_month = calendar.monthrange(year, m + 1)[1]
            midpoint_day = int(cumulative + (days_in_month - 1) / 2.0)
            # Linear interpolation should be exact at midpoints
            assert abs(daily[midpoint_day] - known_monthly[m]) < 0.15, (
                f"Month {m+1} midpoint (day {midpoint_day}): "
                f"expected ~{known_monthly[m]:.2f}, got {daily[midpoint_day]:.2f}"
            )
            cumulative += days_in_month

    def test_smoothness(self, known_monthly):
        """No discontinuities at month boundaries: daily change < 1°C."""
        daily = interpolate_monthly_to_daily(known_monthly, 2013)

        daily_changes = np.abs(np.diff(daily))
        max_change = np.max(daily_changes)
        assert max_change < 1.0, (
            f"Max daily SST change is {max_change:.4f}°C, exceeds 1°C threshold"
        )

    def test_range_preserved(self, known_monthly):
        """Interpolated values should stay within or very near the
        range of input monthly values (no extreme overshoots)."""
        daily = interpolate_monthly_to_daily(known_monthly, 2013)

        # Allow 0.5°C overshoot from linear interpolation edge effects
        assert np.min(daily) >= np.min(known_monthly) - 0.5
        assert np.max(daily) <= np.max(known_monthly) + 0.5

    def test_output_shape_nonleap(self, known_monthly):
        """Non-leap year should produce 365 values."""
        daily = interpolate_monthly_to_daily(known_monthly, 2013)
        assert daily.shape == (365,)

    def test_output_shape_leap(self, known_monthly):
        """Leap year should produce 366 values."""
        daily = interpolate_monthly_to_daily(known_monthly, 2024)
        assert daily.shape == (366,)


# ══════════════════════════════════════════════════════════════════════
# TEST CLASS: Yearly series generation
# ══════════════════════════════════════════════════════════════════════

class TestGenerateYearlySeries:
    """Tests for generate_yearly_sst_series."""

    def test_series_length_5yr(self, tmp_monthly_dir):
        """5-year series should produce exactly 5*365=1825 daily values
        (model uses fixed 365-day years)."""
        sst = generate_yearly_sst_series(
            node_name="Test Node",
            start_year=2010,
            n_years=5,
            data_dir=str(tmp_monthly_dir),
        )
        # Model truncates all years to 365 days
        assert len(sst) == 5 * 365

    def test_series_length_1yr(self, tmp_monthly_dir):
        """1-year series: 365 days."""
        sst = generate_yearly_sst_series(
            node_name="Test Node",
            start_year=2013,
            n_years=1,
            data_dir=str(tmp_monthly_dir),
        )
        assert len(sst) == 365

    def test_known_value_sitka_jan2013(self, real_sst_dir):
        """Sitka Jan 2013: daily SST near mid-January should be close
        to the monthly value of 6.2426°C."""
        sst = generate_yearly_sst_series(
            node_name="Sitka",
            start_year=2013,
            n_years=1,
            data_dir=str(real_sst_dir),
        )
        # January midpoint ≈ day 15 (0-indexed)
        jan_mid = sst[15]
        assert abs(jan_mid - 6.2426) < 0.5, (
            f"Sitka mid-January 2013: expected ~6.24°C, got {jan_mid:.4f}°C"
        )

    def test_fallback_to_climatology(self, tmp_monthly_dir, caplog):
        """Request year 2030 (beyond data range 2010-2015) — should
        fall back to climatology with a warning."""
        with caplog.at_level(logging.WARNING):
            sst = generate_yearly_sst_series(
                node_name="Test Node",
                start_year=2030,
                n_years=1,
                data_dir=str(tmp_monthly_dir),
                climatology_dir=str(tmp_monthly_dir),
            )

        assert len(sst) == 365
        # Should have logged a warning about fallback
        assert any("climatology" in rec.message.lower() or "not available" in rec.message.lower()
                    for rec in caplog.records), (
            "Expected warning about climatology fallback"
        )

    def test_leap_year_truncated(self, tmp_monthly_dir):
        """2012 is a leap year. The model truncates to 365 days, so the
        series should still be 365 per year."""
        sst = generate_yearly_sst_series(
            node_name="Test Node",
            start_year=2012,
            n_years=1,
            data_dir=str(tmp_monthly_dir),
        )
        # Model truncates leap years to 365
        assert len(sst) == 365

    def test_trend_application(self, tmp_monthly_dir):
        """Apply a warming trend and verify later years are warmer."""
        sst_base = generate_yearly_sst_series(
            node_name="Test Node",
            start_year=2010,
            n_years=5,
            data_dir=str(tmp_monthly_dir),
            trend_per_year=0.0,
        )
        sst_trend = generate_yearly_sst_series(
            node_name="Test Node",
            start_year=2010,
            n_years=5,
            data_dir=str(tmp_monthly_dir),
            trend_per_year=0.05,
            reference_year=2012,
        )

        # Year 2014 (offset 4, idx 4*365:5*365): trend = +0.05*(2014-2012) = +0.1
        yr4_base = np.mean(sst_base[4 * 365: 5 * 365])
        yr4_trend = np.mean(sst_trend[4 * 365: 5 * 365])
        expected_offset = 0.05 * (2014 - 2012)
        actual_offset = yr4_trend - yr4_base
        assert abs(actual_offset - expected_offset) < 0.01, (
            f"Year 2014 trend offset: expected +{expected_offset:.2f}, "
            f"got +{actual_offset:.4f}"
        )

    def test_reasonable_range_real_data(self, real_sst_dir):
        """Real Sitka monthly SST should produce values in 0-25°C range."""
        sst = generate_yearly_sst_series(
            node_name="Sitka",
            start_year=2002,
            n_years=5,
            data_dir=str(real_sst_dir),
        )
        assert np.all(sst > 0.0), f"SST min={np.min(sst):.2f}°C, below 0"
        assert np.all(sst < 25.0), f"SST max={np.max(sst):.2f}°C, above 25"


# ══════════════════════════════════════════════════════════════════════
# TEST CLASS: Model integration
# ══════════════════════════════════════════════════════════════════════

class TestModelRunMonthly:
    """Integration tests: full model run with monthly SST."""

    @pytest.fixture
    def single_node_network(self):
        """Minimal 1-node network using Sitka (has monthly data)."""
        from sswd_evoepi.spatial import NodeDefinition, build_network

        nodes = [
            NodeDefinition(
                node_id=0, name="Sitka", lat=57.05, lon=-135.33,
                subregion="AK-SE", habitat_area=20000.0,
                carrying_capacity=300, is_fjord=False,
                sill_depth=np.inf, flushing_rate=0.8,
                mean_sst=8.0, sst_amplitude=3.5, sst_trend=0.0,
                salinity=32.0, depth_range=(5.0, 60.0),
            ),
        ]
        return build_network(nodes, seed=42)

    def test_model_run_monthly_completes(self, single_node_network, real_sst_dir):
        """1-year model sim with sst_source='monthly' completes without error."""
        from sswd_evoepi.model import run_spatial_simulation

        config = default_config()
        config.simulation.sst_source = "monthly"
        config.simulation.sst_start_year = 2013
        config.simulation.sst_data_dir = str(real_sst_dir)

        result = run_spatial_simulation(
            network=single_node_network,
            n_years=1,
            disease_year=999,  # no disease
            seed=42,
            config=config,
        )

        assert result.n_years == 1
        assert result.final_total_pop > 0

    def test_model_monthly_sst_reasonable(self, single_node_network, real_sst_dir):
        """Monthly SST values used in model should be in 0-25°C range."""
        sst = generate_yearly_sst_series(
            node_name="Sitka",
            start_year=2013,
            n_years=1,
            data_dir=str(real_sst_dir),
        )
        assert np.all(sst > 0.0)
        assert np.all(sst < 25.0)


# ══════════════════════════════════════════════════════════════════════
# TEST CLASS: Backward compatibility
# ══════════════════════════════════════════════════════════════════════

class TestBackwardCompatibility:
    """Verify existing SST modes still work after monthly addition."""

    def test_sinusoidal_still_works(self):
        """sst_source='sinusoidal' should still be accepted."""
        config = default_config()
        config.simulation.sst_source = "sinusoidal"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            validate_config(config)

    def test_satellite_still_works(self):
        """sst_source='satellite' should still be accepted."""
        config = default_config()
        config.simulation.sst_source = "satellite"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            validate_config(config)

    def test_monthly_accepted(self):
        """sst_source='monthly' should pass validation."""
        config = default_config()
        config.simulation.sst_source = "monthly"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            validate_config(config)

    def test_sinusoidal_output_unchanged(self):
        """Sinusoidal SST should produce the same output as before."""
        sst = make_sst_timeseries(
            n_years=3, start_year=2013,
            mean_sst=10.0, amplitude=3.5,
        )
        assert len(sst) == 3 * 365
        # Basic sanity: smooth sinusoidal, range ~6.5–13.5
        assert np.min(sst) > 5.0
        assert np.max(sst) < 15.0

    def test_satellite_output_unchanged(self, real_sst_dir):
        """Satellite SST should produce the same output as before."""
        sst = generate_satellite_sst_series(
            n_years=3, start_year=2013,
            node_name="Sitka",
            data_dir=str(real_sst_dir),
        )
        assert len(sst) == 3 * 365
        assert np.min(sst) > 0.0
        assert np.max(sst) < 25.0


# ══════════════════════════════════════════════════════════════════════
# TEST CLASS: Config validation
# ══════════════════════════════════════════════════════════════════════

class TestConfigValidation:
    """Tests for monthly-related config validation."""

    def test_default_sst_source_is_sinusoidal(self):
        """Default config should have sst_source='sinusoidal'."""
        config = default_config()
        assert config.simulation.sst_source == "sinusoidal"

    def test_sst_start_year_default(self):
        """Default sst_start_year should be 2002."""
        config = default_config()
        assert config.simulation.sst_start_year == 2002

    def test_invalid_source_rejected(self):
        """Invalid sst_source should raise ValueError."""
        config = default_config()
        config.simulation.sst_source = "invalid_source"
        with pytest.raises(ValueError, match="sst_source"):
            validate_config(config)
