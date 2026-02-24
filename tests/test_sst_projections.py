"""Tests for SST projection (CMIP6 scenario) integration.

Tests cover:
  - Config: sst_scenario field validation
  - Loading: load_sst_projection_monthly
  - Integration: generate_yearly_sst_series with scenario parameter
  - Continuity: smooth OISST→projection transition
  - Backward compat: observed_only default doesn't break anything
"""

import csv
import os
import tempfile

import numpy as np
import pytest

from sswd_evoepi.config import (
    SimulationConfig,
    SimulationSection,
    default_config,
    validate_config,
)
from sswd_evoepi.environment import (
    generate_yearly_sst_series,
    load_sst_projection_monthly,
)


# ─── Fixtures ─────────────────────────────────────────────────────────

@pytest.fixture
def projection_dir(tmp_path):
    """Create a temporary projection directory with synthetic data."""
    proj_dir = tmp_path / "projections"
    proj_dir.mkdir()

    # Create synthetic SSP245 monthly CSV for TestNode
    ssp245_path = proj_dir / "TestNode_ssp245_monthly.csv"
    with open(ssp245_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['year', 'month', 'sst'])
        for yr in range(2026, 2101):
            for mo in range(1, 13):
                # Simple sinusoidal + warming trend
                sst = 10.0 + 5.0 * np.cos(2 * np.pi * (mo - 8) / 12)
                sst += 0.03 * (yr - 2026)  # 0.03°C/yr warming
                writer.writerow([yr, mo, f"{sst:.4f}"])

    # Create placeholder SSP585 CSV
    ssp585_path = proj_dir / "TestNode_ssp585_placeholder_monthly.csv"
    with open(ssp585_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['year', 'month', 'sst'])
        for yr in range(2026, 2101):
            for mo in range(1, 13):
                sst = 10.0 + 5.0 * np.cos(2 * np.pi * (mo - 8) / 12)
                sst += 0.08 * (yr - 2026)  # faster warming
                writer.writerow([yr, mo, f"{sst:.4f}"])

    return str(proj_dir)


@pytest.fixture
def obs_dir(tmp_path):
    """Create a temporary observation directory with synthetic OISST data."""
    obs_dir = tmp_path / "obs"
    obs_dir.mkdir()

    # Create synthetic monthly OISST for TestNode (2002-2025)
    monthly_path = obs_dir / "TestNode_monthly.csv"
    with open(monthly_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['year', 'month', 'sst'])
        for yr in range(2002, 2026):
            for mo in range(1, 13):
                sst = 10.0 + 5.0 * np.cos(2 * np.pi * (mo - 8) / 12)
                sst += 0.02 * (yr - 2015)  # slight observed warming
                writer.writerow([yr, mo, f"{sst:.4f}"])

    # Create climatology for TestNode (for fallback tests)
    clim_path = obs_dir / "TestNode_climatology.csv"
    with open(clim_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['day_of_year', 'sst_mean', 'sst_std', 'n_years'])
        for doy in range(1, 366):
            sst = 10.0 + 5.0 * np.cos(2 * np.pi * (doy - 227) / 365)
            writer.writerow([doy, f"{sst:.4f}", "0.5000", "20"])

    return str(obs_dir)


# ─── Config Tests ─────────────────────────────────────────────────────

class TestConfigSSTScenario:
    """Test sst_scenario configuration field."""

    def test_default_is_observed_only(self):
        config = default_config()
        assert config.simulation.sst_scenario == 'observed_only'

    def test_default_projection_dir(self):
        config = default_config()
        assert config.simulation.sst_projection_dir == 'data/sst/projections'

    def test_default_obs_end_year(self):
        config = default_config()
        assert config.simulation.sst_obs_end_year == 2025

    def test_valid_scenarios_accepted(self):
        for scenario in ['observed_only', 'ssp126', 'ssp245', 'ssp370', 'ssp585']:
            config = SimulationConfig()
            config.simulation.sst_scenario = scenario
            # Should not raise (projection dir warning is fine)
            try:
                validate_config(config)
            except ValueError:
                pytest.fail(f"sst_scenario='{scenario}' should be valid")

    def test_invalid_scenario_rejected(self):
        config = SimulationConfig()
        config.simulation.sst_scenario = 'rcp85'
        with pytest.raises(ValueError, match="sst_scenario"):
            validate_config(config)


# ─── Load Projection Tests ───────────────────────────────────────────

class TestLoadProjectionMonthly:
    """Test load_sst_projection_monthly."""

    def test_load_exact_match(self, projection_dir):
        data = load_sst_projection_monthly('TestNode', 'ssp245', projection_dir)
        assert len(data) == 75 * 12  # 2026-2100
        assert (2026, 1) in data
        assert (2100, 12) in data

    def test_load_placeholder_fallback(self, projection_dir):
        """Placeholder file is found when exact match doesn't exist."""
        data = load_sst_projection_monthly('TestNode', 'ssp585', projection_dir)
        assert len(data) == 75 * 12
        assert (2026, 1) in data

    def test_missing_scenario_raises(self, projection_dir):
        with pytest.raises(FileNotFoundError, match="ssp370"):
            load_sst_projection_monthly('TestNode', 'ssp370', projection_dir)

    def test_missing_node_raises(self, projection_dir):
        with pytest.raises(FileNotFoundError, match="Nonexistent"):
            load_sst_projection_monthly('Nonexistent', 'ssp245', projection_dir)

    def test_values_reasonable(self, projection_dir):
        data = load_sst_projection_monthly('TestNode', 'ssp245', projection_dir)
        ssts = list(data.values())
        assert min(ssts) > -5.0
        assert max(ssts) < 35.0


# ─── Yearly Series with Scenario Tests ───────────────────────────────

class TestGenerateYearlySeriesWithScenario:
    """Test generate_yearly_sst_series with scenario parameter."""

    def test_observed_only_default_backward_compat(self, obs_dir):
        """Default scenario='observed_only' works like before."""
        series = generate_yearly_sst_series(
            node_name='TestNode',
            start_year=2010,
            n_years=10,
            data_dir=obs_dir,
            climatology_dir=obs_dir,
        )
        assert series.shape == (10 * 365,)
        assert np.all(np.isfinite(series))

    def test_ssp245_uses_projections_after_obs(self, obs_dir, projection_dir):
        """SSP245 scenario uses projections for years > 2025."""
        series = generate_yearly_sst_series(
            node_name='TestNode',
            start_year=2024,
            n_years=4,  # 2024, 2025, 2026, 2027
            data_dir=obs_dir,
            climatology_dir=obs_dir,
            scenario='ssp245',
            projection_dir=projection_dir,
            obs_end_year=2025,
        )
        assert series.shape == (4 * 365,)
        assert np.all(np.isfinite(series))

    def test_ssp585_uses_placeholder(self, obs_dir, projection_dir):
        """SSP585 falls back to placeholder file."""
        series = generate_yearly_sst_series(
            node_name='TestNode',
            start_year=2024,
            n_years=4,
            data_dir=obs_dir,
            climatology_dir=obs_dir,
            scenario='ssp585',
            projection_dir=projection_dir,
            obs_end_year=2025,
        )
        assert series.shape == (4 * 365,)
        assert np.all(np.isfinite(series))

    def test_projection_years_warmer_than_obs(self, obs_dir, projection_dir):
        """SSP245 projections for 2090s should be warmer than 2020s observations."""
        obs_series = generate_yearly_sst_series(
            node_name='TestNode',
            start_year=2020,
            n_years=5,
            data_dir=obs_dir,
            climatology_dir=obs_dir,
        )
        proj_series = generate_yearly_sst_series(
            node_name='TestNode',
            start_year=2096,
            n_years=5,
            data_dir=obs_dir,
            climatology_dir=obs_dir,
            scenario='ssp245',
            projection_dir=projection_dir,
            obs_end_year=2025,
        )
        assert np.mean(proj_series) > np.mean(obs_series)

    def test_ssp585_warmer_than_ssp245(self, obs_dir, projection_dir):
        """SSP585 should be warmer than SSP245 by end of century."""
        ssp245 = generate_yearly_sst_series(
            node_name='TestNode',
            start_year=2090,
            n_years=5,
            data_dir=obs_dir,
            climatology_dir=obs_dir,
            scenario='ssp245',
            projection_dir=projection_dir,
            obs_end_year=2025,
        )
        ssp585 = generate_yearly_sst_series(
            node_name='TestNode',
            start_year=2090,
            n_years=5,
            data_dir=obs_dir,
            climatology_dir=obs_dir,
            scenario='ssp585',
            projection_dir=projection_dir,
            obs_end_year=2025,
        )
        assert np.mean(ssp585) > np.mean(ssp245)

    def test_output_length_long_run(self, obs_dir, projection_dir):
        """Full 2002-2100 run produces correct output length."""
        series = generate_yearly_sst_series(
            node_name='TestNode',
            start_year=2002,
            n_years=99,  # 2002-2100
            data_dir=obs_dir,
            climatology_dir=obs_dir,
            scenario='ssp245',
            projection_dir=projection_dir,
            obs_end_year=2025,
        )
        assert series.shape == (99 * 365,)
        assert np.all(np.isfinite(series))


# ─── Continuity Tests ────────────────────────────────────────────────

class TestTransitionContinuity:
    """Test OISST→projection transition is within reasonable bounds."""

    def test_no_extreme_jumps(self, obs_dir, projection_dir):
        """Day-to-day SST change at 2025/2026 boundary should be reasonable."""
        series = generate_yearly_sst_series(
            node_name='TestNode',
            start_year=2024,
            n_years=4,
            data_dir=obs_dir,
            climatology_dir=obs_dir,
            scenario='ssp245',
            projection_dir=projection_dir,
            obs_end_year=2025,
        )
        # Get the boundary: last day of year 2 (=2025) and first day of year 3 (=2026)
        boundary_idx = 2 * 365 - 1  # last day of 2025
        jump = abs(series[boundary_idx + 1] - series[boundary_idx])
        # Daily SST change should be < 5°C (generous; typical < 1°C)
        assert jump < 5.0, f"Transition jump of {jump:.2f}°C is too large"


# ─── Real Data Tests ─────────────────────────────────────────────────

class TestRealDataProjections:
    """Tests using actual CMIP6 projection data (skip if not available)."""

    @pytest.fixture
    def real_proj_dir(self):
        proj_dir = 'data/sst/projections'
        if not os.path.isdir(proj_dir):
            pytest.skip("Projection data not available")
        return proj_dir

    @pytest.fixture
    def real_obs_dir(self):
        obs_dir = 'data/sst'
        if not os.path.isfile(os.path.join(obs_dir, 'Sitka_monthly.csv')):
            pytest.skip("OISST data not available")
        return obs_dir

    def test_load_real_sitka_ssp245(self, real_proj_dir):
        data = load_sst_projection_monthly('Sitka', 'ssp245', real_proj_dir)
        assert len(data) > 0
        assert (2026, 1) in data
        ssts = list(data.values())
        # Sitka SST should be between -2 and 25°C
        assert min(ssts) > -2.0
        assert max(ssts) < 25.0

    def test_real_full_series_sitka(self, real_obs_dir, real_proj_dir):
        """Generate a full 2002-2100 SST series for Sitka with SSP245."""
        series = generate_yearly_sst_series(
            node_name='Sitka',
            start_year=2002,
            n_years=99,
            data_dir=real_obs_dir,
            climatology_dir=real_obs_dir,
            scenario='ssp245',
            projection_dir=real_proj_dir,
            obs_end_year=2025,
        )
        assert series.shape == (99 * 365,)
        assert np.all(np.isfinite(series))

        # Check warming trend: 2090s should be warmer than 2010s
        mean_2010s = np.mean(series[8*365:18*365])   # years 2010-2019
        mean_2090s = np.mean(series[88*365:98*365])   # years 2090-2099
        assert mean_2090s > mean_2010s, (
            f"Expected 2090s warmer than 2010s: {mean_2090s:.2f} vs {mean_2010s:.2f}"
        )

    def test_all_11_nodes_have_ssp245(self, real_proj_dir):
        """All 11 stepping-stone nodes have SSP245 projections."""
        nodes = [
            'Sitka', 'Ketchikan', 'Haida_Gwaii', 'Bella_Bella', 'Howe_Sound',
            'SJI', 'Westport', 'Newport', 'Crescent_City', 'Fort_Bragg', 'Monterey',
        ]
        for node in nodes:
            data = load_sst_projection_monthly(node, 'ssp245', real_proj_dir)
            assert len(data) > 0, f"No SSP245 data for {node}"
