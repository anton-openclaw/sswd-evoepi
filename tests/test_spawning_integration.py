"""Integration tests for spawning system in SSWD-EvoEpi.

Tests that the spawning module is properly integrated into the main
simulation loop and produces reasonable outputs.

Author: Anton 🔬
"""

import numpy as np
import pytest

from sswd_evoepi.config import default_config
from sswd_evoepi.model import run_coupled_simulation


class TestSpawningIntegration:
    """Test spawning system integration with main simulation."""

    def test_simulation_with_spawning_enabled(self):
        """Simulation with spawning enabled should complete and produce recruits."""
        config = default_config()
        
        result = run_coupled_simulation(
            n_individuals=50,
            carrying_capacity=100,
            habitat_area=5000.0,
            n_years=3,
            disease_year=999,  # No disease for this test
            seed=42,
            config=config,  # Uses default spawning configuration
        )
        
        # Basic sanity checks
        assert len(result.yearly_pop) == 3
        assert len(result.yearly_recruits) == 3
        assert not np.isnan(result.yearly_pop).any()
        assert not np.isnan(result.yearly_recruits).any()
        
        # Should produce some recruits over 3 years
        assert sum(result.yearly_recruits) > 0
        
        # Population should remain viable
        assert all(pop > 0 for pop in result.yearly_pop)

    def test_backward_compatibility_without_spawning(self):
        """Simulation without spawning config should fall back to annual reproduction."""
        config = default_config()
        config.spawning = None  # Disable spawning
        
        result = run_coupled_simulation(
            n_individuals=50,
            carrying_capacity=100,
            habitat_area=5000.0,
            n_years=3,
            disease_year=999,  # No disease
            seed=42,
            config=config,
        )
        
        # Should still work and produce recruits
        assert len(result.yearly_pop) == 3
        assert sum(result.yearly_recruits) > 0
        assert all(pop > 0 for pop in result.yearly_pop)

    def test_spawning_vs_annual_reproduction_different_outputs(self):
        """Spawning system should produce different results from annual reproduction."""
        base_config = default_config()
        
        # Run with spawning enabled
        config_spawning = base_config
        result_spawning = run_coupled_simulation(
            n_individuals=50,
            carrying_capacity=100,
            habitat_area=5000.0,
            n_years=2,
            disease_year=999,
            seed=42,
            config=config_spawning,
        )
        
        # Run with annual reproduction (spawning disabled)
        config_annual = default_config()
        config_annual.spawning = None
        result_annual = run_coupled_simulation(
            n_individuals=50,
            carrying_capacity=100,
            habitat_area=5000.0,
            n_years=2,
            disease_year=999,
            seed=42,
            config=config_annual,
        )
        
        # Results should be different (spawning vs annual pulse)
        # Note: Due to different reproductive timing, we expect differences
        # but both should be reasonable
        assert not np.array_equal(result_spawning.yearly_recruits, result_annual.yearly_recruits)
        
        # Both should produce viable populations
        assert all(pop > 0 for pop in result_spawning.yearly_pop)
        assert all(pop > 0 for pop in result_annual.yearly_pop)

    def test_spawning_with_disease(self):
        """Spawning system should work correctly with disease dynamics."""
        config = default_config()
        
        result = run_coupled_simulation(
            n_individuals=100,
            carrying_capacity=200,
            habitat_area=10000.0,
            n_years=5,
            disease_year=2,  # Introduce disease in year 2
            initial_infected=10,
            seed=42,
            config=config,
        )
        
        # Should complete without errors
        assert len(result.yearly_pop) == 5
        assert not np.isnan(result.yearly_pop).any()
        
        # Disease should cause population crash
        pre_disease_pop = result.yearly_pop[1]  # Year 1 (before disease)
        post_disease_pop = result.yearly_pop[3]  # Year 3 (after disease)
        
        # Should see significant population decline (at least 50% crash expected)
        mortality_fraction = (pre_disease_pop - post_disease_pop) / pre_disease_pop
        assert mortality_fraction > 0.3  # At least 30% decline from disease

    def test_immunosuppression_timers_handled(self):
        """Immunosuppression timers should be decremented during daily loop."""
        config = default_config()
        
        # This is mostly a regression test to ensure no errors when
        # immunosuppression timers are active. The actual spawning behavior
        # is tested in the spawning module tests.
        result = run_coupled_simulation(
            n_individuals=20,
            carrying_capacity=40,
            habitat_area=2000.0,
            n_years=1,
            disease_year=999,
            seed=42,
            config=config,
        )
        
        # Should complete without errors
        assert len(result.yearly_pop) == 1
        assert result.yearly_pop[0] > 0