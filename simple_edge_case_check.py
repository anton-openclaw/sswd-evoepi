#!/usr/bin/env python3
"""
Simple edge case verification for spawning system
Focusing on basic functionality without complex API interactions.
"""

import numpy as np
import sys
sys.path.append('.')

from sswd_evoepi.model import run_coupled_simulation
from sswd_evoepi.config import load_config
from sswd_evoepi.spawning import in_spawning_season, reset_spawning_season
from sswd_evoepi.types import AGENT_DTYPE

def test_basic_simulation():
    """Test basic simulation runs without crashing"""
    print("Testing basic simulation with spawning...")
    
    try:
        result = run_coupled_simulation(
            n_individuals=50,
            n_years=2,
            disease_year=5,  # No disease
            seed=42,
            record_daily=False
        )
        
        # Check basic properties
        pop_data = result.yearly_pop
        assert len(pop_data) == 2, f"Expected 2 years of data, got {len(pop_data)}"
        assert np.all(pop_data >= 0), f"Negative population: {np.min(pop_data)}"
        assert not np.any(np.isnan(pop_data)), "NaN found in population data"
        
        print(f"  Final population: {pop_data[-1]}")
        print("  ✓ Basic simulation test passed")
        
    except Exception as e:
        print(f"  ✗ Basic simulation test failed: {e}")
        raise

def test_season_detection():
    """Test spawning season detection logic"""
    print("Testing spawning season detection...")
    
    # Test known season boundaries
    test_cases = [
        (120, True, "April - should be in season"),  # April
        (320, True, "November - should be in season"),  # November 
        (250, False, "September - should be out of season"),  # September
        (1, True, "January - should be in season"),  # January
        (200, False, "July end - should be out of season"),  # July end
    ]
    
    for doy, expected, description in test_cases:
        result = in_spawning_season(doy, 305, 196)
        assert result == expected, f"{description}: got {result}, expected {expected}"
        
    print("  ✓ Season detection test passed")

def test_small_population():
    """Test simulation with very small population"""
    print("Testing small population (10 individuals)...")
    
    try:
        result = run_coupled_simulation(
            n_individuals=10,
            n_years=2,
            disease_year=5,  # No disease
            seed=42
        )
        
        # Should handle small populations gracefully
        final_pop = result.yearly_pop[-1]
        print(f"  Final population: {final_pop} (started with 10)")
        
        # Population might decline but shouldn't crash
        assert final_pop >= 0, f"Population went negative: {final_pop}"
        print("  ✓ Small population test passed")
        
    except Exception as e:
        print(f"  ✗ Small population test failed: {e}")
        raise

def test_with_disease():
    """Test simulation with disease and spawning"""
    print("Testing simulation with disease...")
    
    try:
        result = run_coupled_simulation(
            n_individuals=100,
            n_years=5,
            disease_year=2,  # Introduce disease in year 2
            initial_infected=3,
            seed=42
        )
        
        # Check that disease affects population
        pop_data = result.yearly_pop
        disease_impact = pop_data[1] - pop_data[2]  # Population change during disease year
        
        print(f"  Population before disease: {pop_data[1]}")
        print(f"  Population after disease: {pop_data[2]}")
        print(f"  Disease impact: {disease_impact}")
        
        # Should have some population remaining
        assert pop_data[-1] > 0, f"Population extinct: {pop_data[-1]}"
        print("  ✓ Disease simulation test passed")
        
    except Exception as e:
        print(f"  ✗ Disease simulation test failed: {e}")
        raise

def test_agent_structure():
    """Test basic agent data structure"""
    print("Testing agent data structure...")
    
    try:
        # Create small agent array to check structure
        agents = np.zeros(5, dtype=AGENT_DTYPE)
        
        # Check that basic fields exist
        assert 'alive' in agents.dtype.names, "Missing 'alive' field"
        assert 'sex' in agents.dtype.names, "Missing 'sex' field"
        assert 'age_days' in agents.dtype.names, "Missing 'age_days' field"
        
        # Check spawning-related fields
        spawning_fields = [f for f in agents.dtype.names if 'spawning' in f.lower()]
        print(f"  Spawning-related fields: {spawning_fields}")
        
        # Test reset function doesn't crash
        reset_spawning_season(agents)
        
        print("  ✓ Agent structure test passed")
        
    except Exception as e:
        print(f"  ✗ Agent structure test failed: {e}")
        raise

if __name__ == "__main__":
    print("Running simplified edge case verification for spawning system...")
    print("=" * 70)
    
    try:
        test_basic_simulation()
        test_season_detection()
        test_small_population()
        test_with_disease()
        test_agent_structure()
        
        print("=" * 70)
        print("✓ ALL SIMPLIFIED EDGE CASE TESTS PASSED!")
        
    except Exception as e:
        print("=" * 70)
        print(f"✗ EDGE CASE TEST FAILED: {e}")
        sys.exit(1)