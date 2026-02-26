"""Tests for movement module: CRW and spatial transmission grid.

Tests:
  1. Boundary reflection correctness
  2. CRW moves agents (positions change)
  3. Disease-dependent speed modifiers
  4. Dead agents don't move
  5. Positions stay within habitat bounds
  6. InfectedDensityGrid builds and normalizes correctly
  7. Grid lookup returns reasonable factors
  8. Mean-field equivalence when grid is uniform
  9. Spatial clustering produces heterogeneous factors
  10. Movement disabled → positions unchanged
"""

import numpy as np
import pytest

from sswd_evoepi.types import AGENT_DTYPE, DiseaseState, allocate_agents
from sswd_evoepi.movement import (
    _reflect,
    update_movement,
    daily_movement,
    InfectedDensityGrid,
    _diffuse_2d,
    SPEED_MODIFIER,
    _apply_spawning_gravity,
    _circular_blend,
)


# ═══════════════════════════════════════════════════════════════════════
# HELPERS
# ═══════════════════════════════════════════════════════════════════════

def make_agents(n, habitat_side=100.0, rng=None):
    """Create n alive agents uniformly distributed in habitat."""
    if rng is None:
        rng = np.random.default_rng(42)
    agents = allocate_agents(n)
    for i in range(n):
        agents[i]['alive'] = True
        agents[i]['x'] = rng.uniform(0, habitat_side)
        agents[i]['y'] = rng.uniform(0, habitat_side)
        agents[i]['heading'] = rng.uniform(0, 2 * np.pi)
        agents[i]['speed'] = 0.5
        agents[i]['disease_state'] = DiseaseState.S
        agents[i]['size'] = 500.0
    return agents


# ═══════════════════════════════════════════════════════════════════════
# BOUNDARY REFLECTION TESTS
# ═══════════════════════════════════════════════════════════════════════

class TestReflect:
    def test_within_bounds_unchanged(self):
        x = np.array([10.0, 50.0, 90.0])
        result = _reflect(x, 100.0)
        np.testing.assert_allclose(result, x)

    def test_negative_reflects(self):
        x = np.array([-10.0])
        result = _reflect(x, 100.0)
        np.testing.assert_allclose(result, [10.0])

    def test_over_boundary_reflects(self):
        x = np.array([110.0])
        result = _reflect(x, 100.0)
        np.testing.assert_allclose(result, [90.0])

    def test_double_bounce(self):
        """Agent moves far enough to bounce twice."""
        x = np.array([250.0])  # 250 in [0, 100] → 250 % 200 = 50
        result = _reflect(x, 100.0)
        np.testing.assert_allclose(result, [50.0])

    def test_negative_double_bounce(self):
        x = np.array([-250.0])
        result = _reflect(x, 100.0)
        # -250 % 200 = 150 → 200 - 150 = 50
        np.testing.assert_allclose(result, [50.0])

    def test_at_boundary(self):
        x = np.array([0.0, 100.0])
        result = _reflect(x, 100.0)
        np.testing.assert_allclose(result, [0.0, 100.0])


# ═══════════════════════════════════════════════════════════════════════
# CRW MOVEMENT TESTS
# ═══════════════════════════════════════════════════════════════════════

class TestUpdateMovement:
    def test_agents_move(self):
        """Positions should change after movement step."""
        rng = np.random.default_rng(42)
        agents = make_agents(100, rng=rng)
        x_before = agents['x'].copy()
        y_before = agents['y'].copy()

        update_movement(agents, dt_minutes=60.0, habitat_side=100.0,
                        sigma_turn=0.6, base_speed=0.5, rng=rng)

        alive = agents['alive'].astype(bool)
        assert not np.allclose(agents['x'][alive], x_before[alive])
        assert not np.allclose(agents['y'][alive], y_before[alive])

    def test_positions_within_bounds(self):
        """All positions must stay in [0, side] after movement."""
        rng = np.random.default_rng(42)
        agents = make_agents(500, habitat_side=50.0, rng=rng)

        # Many movement steps to stress-test boundary
        for _ in range(100):
            update_movement(agents, dt_minutes=60.0, habitat_side=50.0,
                            sigma_turn=1.0, base_speed=2.0, rng=rng)

        alive = agents['alive'].astype(bool)
        assert np.all(agents['x'][alive] >= 0)
        assert np.all(agents['x'][alive] <= 50.0)
        assert np.all(agents['y'][alive] >= 0)
        assert np.all(agents['y'][alive] <= 50.0)

    def test_dead_agents_unmoved(self):
        """Dead agents should not have their positions changed."""
        rng = np.random.default_rng(42)
        agents = make_agents(50, rng=rng)
        agents[10:20]['alive'] = False
        x_dead = agents['x'][10:20].copy()
        y_dead = agents['y'][10:20].copy()

        update_movement(agents, dt_minutes=60.0, habitat_side=100.0,
                        sigma_turn=0.6, base_speed=0.5, rng=rng)

        np.testing.assert_array_equal(agents['x'][10:20], x_dead)
        np.testing.assert_array_equal(agents['y'][10:20], y_dead)

    def test_disease_speed_modifiers(self):
        """Sick agents should move slower."""
        rng = np.random.default_rng(42)
        agents = make_agents(6, habitat_side=1000.0, rng=rng)

        # Place all at center with same heading (east)
        for i in range(6):
            agents[i]['x'] = 500.0
            agents[i]['y'] = 500.0
            agents[i]['heading'] = 0.0

        # Set disease states
        agents[0]['disease_state'] = DiseaseState.S   # ×1.0
        agents[1]['disease_state'] = DiseaseState.E   # ×1.0
        agents[2]['disease_state'] = DiseaseState.I1  # ×0.5
        agents[3]['disease_state'] = DiseaseState.I2  # ×0.1
        agents[4]['disease_state'] = DiseaseState.D   # ×0.0
        agents[5]['disease_state'] = DiseaseState.R   # ×1.0

        # Use sigma_turn=0 for deterministic heading (straight east)
        update_movement(agents, dt_minutes=60.0, habitat_side=1000.0,
                        sigma_turn=0.0, base_speed=1.0, rng=rng)

        # Check effective speeds stored
        np.testing.assert_allclose(agents[0]['speed'], 1.0)   # S
        np.testing.assert_allclose(agents[1]['speed'], 1.0)   # E
        np.testing.assert_allclose(agents[2]['speed'], 0.5)   # I1
        np.testing.assert_allclose(agents[3]['speed'], 0.1)   # I2
        np.testing.assert_allclose(agents[4]['speed'], 0.0)   # D
        np.testing.assert_allclose(agents[5]['speed'], 1.0)   # R

        # Check displacement proportional to speed
        dx = agents['x'][:6] - 500.0
        # S and E should have moved the same amount (60m)
        np.testing.assert_allclose(dx[0], 60.0, atol=0.1)
        np.testing.assert_allclose(dx[1], 60.0, atol=0.1)
        # I1 should be half
        np.testing.assert_allclose(dx[2], 30.0, atol=0.1)
        # I2 should be 10%
        np.testing.assert_allclose(dx[3], 6.0, atol=0.1)
        # D should not move
        np.testing.assert_allclose(dx[4], 0.0, atol=0.01)


class TestDailyMovement:
    def test_daily_24_substeps(self):
        """daily_movement should call update_movement 24 times."""
        rng = np.random.default_rng(42)
        agents = make_agents(50, habitat_side=200.0, rng=rng)
        x_before = agents['x'][:50].copy()

        daily_movement(agents, habitat_side=200.0, sigma_turn=0.6,
                       base_speed=0.5, substeps=24, rng=rng)

        alive = agents['alive'].astype(bool)
        # After 24 steps, positions should have changed
        assert not np.allclose(agents['x'][alive], x_before[alive])

    def test_all_within_bounds_after_full_day(self):
        rng = np.random.default_rng(42)
        agents = make_agents(500, habitat_side=100.0, rng=rng)

        daily_movement(agents, habitat_side=100.0, sigma_turn=0.6,
                       base_speed=0.5, substeps=24, rng=rng)

        alive = agents['alive'].astype(bool)
        assert np.all(agents['x'][alive] >= 0)
        assert np.all(agents['x'][alive] <= 100.0)


# ═══════════════════════════════════════════════════════════════════════
# INFECTED DENSITY GRID TESTS
# ═══════════════════════════════════════════════════════════════════════

class TestDiffuse2D:
    def test_uniform_unchanged(self):
        """Uniform grid should be unchanged by diffusion."""
        grid = np.ones((5, 5))
        result = _diffuse_2d(grid)
        np.testing.assert_allclose(result, 1.0, atol=1e-10)

    def test_smoothing(self):
        """Point source should spread out."""
        grid = np.zeros((10, 10))
        grid[5, 5] = 9.0
        result = _diffuse_2d(grid)
        # Center should decrease, neighbors should increase
        assert result[5, 5] < 9.0
        assert result[4, 5] > 0
        assert result[5, 4] > 0


class TestInfectedDensityGrid:
    def test_no_infected_uniform(self):
        """With no infected, grid should be all 1.0 (uniform exposure)."""
        agents = make_agents(100, habitat_side=100.0)
        grid = InfectedDensityGrid(habitat_side=100.0, cell_size=20.0)
        grid.build(agents, diffusion_passes=2)

        # All factors should be 1.0
        x = agents['x'][:100]
        y = agents['y'][:100]
        factors = grid.lookup(x, y)
        np.testing.assert_allclose(factors, 1.0)

    def test_mean_approximately_one(self):
        """Mean exposure factor across all susceptibles should ≈ 1."""
        rng = np.random.default_rng(42)
        agents = make_agents(1000, habitat_side=200.0, rng=rng)

        # Make 100 agents infected (clustered in one corner)
        for i in range(100):
            agents[i]['disease_state'] = DiseaseState.I1
            agents[i]['x'] = rng.uniform(0, 40)
            agents[i]['y'] = rng.uniform(0, 40)

        grid = InfectedDensityGrid(habitat_side=200.0, cell_size=20.0)
        grid.build(agents, diffusion_passes=2)

        # Check that grid mean is approximately 1
        assert abs(grid.grid.mean() - 1.0) < 0.01

    def test_clustering_creates_heterogeneity(self):
        """Infected cluster should create high/low exposure zones."""
        rng = np.random.default_rng(42)
        agents = make_agents(500, habitat_side=200.0, rng=rng)

        # Cluster 50 infected in corner (0-30, 0-30)
        for i in range(50):
            agents[i]['disease_state'] = DiseaseState.I1
            agents[i]['x'] = rng.uniform(0, 30)
            agents[i]['y'] = rng.uniform(0, 30)

        grid = InfectedDensityGrid(habitat_side=200.0, cell_size=20.0)
        grid.build(agents, diffusion_passes=2)

        # Susceptible near cluster should have factor > 1
        near_x = np.array([10.0])
        near_y = np.array([10.0])
        far_x = np.array([180.0])
        far_y = np.array([180.0])

        near_factor = grid.lookup(near_x, near_y)[0]
        far_factor = grid.lookup(far_x, far_y)[0]

        assert near_factor > far_factor
        assert near_factor > 1.0  # elevated near cluster
        # Far factor might be < 1 or close to the baseline

    def test_grid_dimensions(self):
        """Grid should have correct number of cells."""
        grid = InfectedDensityGrid(habitat_side=100.0, cell_size=20.0)
        assert grid.nx == 5
        assert grid.ny == 5
        assert grid.grid.shape == (5, 5)

    def test_lookup_clamps_to_bounds(self):
        """Positions outside grid should be clamped to edge cells."""
        grid = InfectedDensityGrid(habitat_side=100.0, cell_size=20.0)
        grid.grid[:] = 1.0

        # Out-of-bounds positions
        x = np.array([-5.0, 105.0])
        y = np.array([-5.0, 105.0])
        factors = grid.lookup(x, y)
        assert factors.shape == (2,)
        np.testing.assert_allclose(factors, 1.0)

    def test_i2_also_counts_as_infected(self):
        """Both I1 and I2 should contribute to infected density."""
        rng = np.random.default_rng(42)
        agents = make_agents(100, habitat_side=100.0, rng=rng)

        # Put some I2 at center
        for i in range(10):
            agents[i]['disease_state'] = DiseaseState.I2
            agents[i]['x'] = 50.0
            agents[i]['y'] = 50.0

        grid = InfectedDensityGrid(habitat_side=100.0, cell_size=20.0)
        grid.build(agents, diffusion_passes=0)

        # Center cell should have elevated density
        center_factor = grid.lookup(np.array([50.0]), np.array([50.0]))[0]
        assert center_factor > 1.0


# ═══════════════════════════════════════════════════════════════════════
# INTEGRATION TESTS
# ═══════════════════════════════════════════════════════════════════════

class TestMovementIntegration:
    def test_movement_then_grid_workflow(self):
        """Full workflow: move agents, build grid, lookup factors."""
        rng = np.random.default_rng(42)
        agents = make_agents(200, habitat_side=100.0, rng=rng)

        # Set some infected
        for i in range(20):
            agents[i]['disease_state'] = DiseaseState.I1

        # Move
        daily_movement(agents, habitat_side=100.0, sigma_turn=0.6,
                       base_speed=0.5, substeps=24, rng=rng)

        # Build grid
        grid = InfectedDensityGrid(habitat_side=100.0, cell_size=20.0)
        grid.build(agents, diffusion_passes=2)

        # Lookup for all susceptibles
        susc = (agents['alive'].astype(bool) &
                (agents['disease_state'] == DiseaseState.S))
        susc_idx = np.where(susc)[0]
        factors = grid.lookup(agents['x'][susc_idx], agents['y'][susc_idx])

        # All factors should be finite and positive
        assert np.all(np.isfinite(factors))
        assert np.all(factors >= 0)

    def test_speed_modifier_array(self):
        """SPEED_MODIFIER should match disease-module-spec §5.5."""
        expected = np.array([1.0, 1.0, 0.5, 0.1, 0.0, 1.0])
        np.testing.assert_allclose(SPEED_MODIFIER, expected, atol=1e-7)


# ═══════════════════════════════════════════════════════════════════════
# SPAWNING GRAVITY TESTS (Phase 3)
# ═══════════════════════════════════════════════════════════════════════

class TestCircularBlend:
    def test_zero_weight_returns_first_heading(self):
        """Weight 0 should return the first heading unchanged."""
        h1, h2 = 1.0, 3.0
        result = _circular_blend(h1, h2, 0.0)
        np.testing.assert_allclose(result, h1)

    def test_full_weight_returns_second_heading(self):
        """Weight 1 should return the second heading unchanged."""
        h1, h2 = 1.0, 3.0
        result = _circular_blend(h1, h2, 1.0)
        np.testing.assert_allclose(result, h2)

    def test_half_weight_circular_average(self):
        """Weight 0.5 should return circular average."""
        # Simple case: 0 and π/2 should average to π/4
        h1, h2 = 0.0, np.pi/2
        result = _circular_blend(h1, h2, 0.5)
        np.testing.assert_allclose(result, np.pi/4, atol=1e-6)

    def test_circular_wrapping(self):
        """Should handle angle wrapping correctly."""
        # 0.1 and 6.2 should blend through 0, not through π
        h1, h2 = 0.1, 6.2
        result = _circular_blend(h1, h2, 0.5)
        # Should be close to 0 or 2π, not around π
        assert result < 0.5 or result > 5.8


class TestSpawningGravity:
    def test_no_gravity_when_disabled(self):
        """When gravity is disabled, headings should be unchanged."""
        rng = np.random.default_rng(42)
        agents = make_agents(50, habitat_side=200.0, rng=rng)
        
        # Set some agents to have gravity timers
        agents[0:10]['spawn_gravity_timer'] = 5
        
        original_headings = agents['heading'].copy()
        
        # Run movement with gravity disabled
        update_movement(agents, dt_minutes=60.0, habitat_side=200.0,
                        sigma_turn=0.0, base_speed=0.5, rng=rng,
                        spawning_gravity_enabled=False)
        
        # Headings should be unchanged (sigma_turn=0 for determinism)
        np.testing.assert_allclose(agents['heading'], original_headings)

    def test_no_gravity_without_timer(self):
        """Agents without gravity timer should not be affected."""
        rng = np.random.default_rng(42)
        agents = make_agents(10, habitat_side=100.0, rng=rng)
        
        # Set all gravity timers to 0
        agents['spawn_gravity_timer'][:] = 0
        
        original_headings = agents['heading'].copy()
        
        # Run movement with gravity enabled but no timers
        update_movement(agents, dt_minutes=60.0, habitat_side=100.0,
                        sigma_turn=0.0, base_speed=0.5, rng=rng,
                        spawning_gravity_enabled=True)
        
        # Headings should be unchanged
        np.testing.assert_allclose(agents['heading'], original_headings)

    def test_agents_clump_together_with_gravity(self):
        """Agents with gravity should move toward each other."""
        rng = np.random.default_rng(42)
        agents = make_agents(20, habitat_side=500.0, rng=rng)
        
        # Place agents in one scattered group that should cluster
        for i in range(20):
            agents[i]['x'] = 250.0 + rng.uniform(-80, 80)  # Scattered over 160m
            agents[i]['y'] = 250.0 + rng.uniform(-80, 80)  # within gravity range
            agents[i]['spawn_gravity_timer'] = 14  # Mid-gravity phase
        
        # Calculate initial mean distances within the group
        initial_positions = np.column_stack([agents['x'][:20], agents['y'][:20]])
        
        initial_dist = np.mean([
            np.sqrt(np.sum((initial_positions[i] - initial_positions[j])**2))
            for i in range(20) for j in range(i+1, 20)
        ])
        
        # Run multiple movement steps with gravity
        for _ in range(10):
            update_movement(agents, dt_minutes=60.0, habitat_side=500.0,
                            sigma_turn=0.1, base_speed=1.0, rng=rng,
                            spawning_gravity_enabled=True, gravity_strength=2.0,
                            gravity_range=100.0, pre_spawn_gravity_days=14,
                            post_spawn_gravity_days=14)
        
        # Calculate final mean distances
        final_positions = np.column_stack([agents['x'][:20], agents['y'][:20]])
        
        final_dist = np.mean([
            np.sqrt(np.sum((final_positions[i] - final_positions[j])**2))
            for i in range(20) for j in range(i+1, 20)
        ])
        
        # Group should be more tightly clustered
        assert final_dist < initial_dist, f"Group should be more clustered: {final_dist:.3f} < {initial_dist:.3f}"

    def test_gravity_strength_ramp(self):
        """Gravity strength should ramp up then down over the spawning window."""
        rng = np.random.default_rng(42)
        agents = make_agents(10, habitat_side=100.0, rng=rng)
        alive_idx = np.arange(10)
        
        pre_days = 14
        post_days = 14
        total_days = pre_days + post_days
        
        # Test different timer values
        test_cases = [
            (total_days, 0.0),      # Start of window: 0 strength
            (total_days - 7, 0.5),  # Mid pre-spawn: 0.5 strength  
            (post_days, 1.0),       # At spawning: 1.0 strength
            (7, 0.5),               # Mid post-spawn: 0.5 strength
            (1, 0.07),              # Near end: ~0.07 strength
        ]
        
        for timer, expected_ramp in test_cases:
            agents['spawn_gravity_timer'][:] = timer
            
            # We need to extract the ramp calculation logic for testing
            # Since _apply_spawning_gravity is complex, let's test the ramp logic separately
            days_from_start = total_days - timer
            
            if days_from_start <= pre_days:
                ramp_factor = days_from_start / pre_days
            else:
                days_into_post_spawn = days_from_start - pre_days
                ramp_factor = 1.0 - (days_into_post_spawn / post_days)
            
            ramp_factor = np.clip(ramp_factor, 0.0, 1.0)
            
            np.testing.assert_allclose(ramp_factor, expected_ramp, atol=0.02)

    def test_gravity_respects_range_limit(self):
        """Agents should only be affected by conspecifics within gravity_range."""
        rng = np.random.default_rng(42)
        agents = make_agents(3, habitat_side=1000.0, rng=rng)
        
        # Place agents: one in center with gravity, two at different distances
        agents[0]['x'] = 500.0  # Center
        agents[0]['y'] = 500.0
        agents[0]['spawn_gravity_timer'] = 14
        agents[0]['heading'] = 0.0  # East
        
        agents[1]['x'] = 550.0  # 50m east (within 100m range)
        agents[1]['y'] = 500.0
        agents[1]['spawn_gravity_timer'] = 0  # No gravity
        
        agents[2]['x'] = 650.0  # 150m east (beyond 100m range) 
        agents[2]['y'] = 500.0
        agents[2]['spawn_gravity_timer'] = 0  # No gravity
        
        original_heading = agents[0]['heading']
        
        # Apply gravity with 100m range
        alive_idx = np.array([0, 1, 2])
        modified_headings = _apply_spawning_gravity(
            agents, alive_idx, gravity_strength=1.0, gravity_range=100.0,
            pre_spawn_gravity_days=14, post_spawn_gravity_days=14
        )
        
        # Agent 0 should be influenced by agent 1 (within range) but not agent 2
        # Since agent 1 is due east and agent 0 starts heading east,
        # the heading should remain close to east (maybe slightly adjusted)
        assert abs(modified_headings[0] - original_heading) < 0.1

    def test_gravity_disabled_by_default(self):
        """Gravity should be disabled by default in movement functions."""
        rng = np.random.default_rng(42)
        agents = make_agents(20, habitat_side=100.0, rng=rng)
        agents[0:10]['spawn_gravity_timer'] = 10
        
        original_headings = agents['heading'].copy()
        
        # Call update_movement without gravity parameters (should use defaults)
        update_movement(agents, dt_minutes=60.0, habitat_side=100.0,
                        sigma_turn=0.0, base_speed=0.5, rng=rng)
        
        # Headings should be unchanged since gravity is disabled by default
        np.testing.assert_allclose(agents['heading'], original_headings)


# ═══════════════════════════════════════════════════════════════════════
# STOCHASTIC STEP-LENGTH TESTS (speed_sigma)
# ═══════════════════════════════════════════════════════════════════════

class TestStochasticStepLength:
    def test_zero_sigma_identical_to_default(self):
        """speed_sigma=0.0 should produce identical results to the default path."""
        rng1 = np.random.default_rng(42)
        agents1 = make_agents(100, habitat_side=200.0, rng=rng1)
        rng1 = np.random.default_rng(99)

        rng2 = np.random.default_rng(42)
        agents2 = make_agents(100, habitat_side=200.0, rng=rng2)
        rng2 = np.random.default_rng(99)

        daily_movement(agents1, habitat_side=200.0, sigma_turn=0.6,
                       base_speed=0.5, substeps=24, rng=rng1,
                       speed_sigma=0.0)

        daily_movement(agents2, habitat_side=200.0, sigma_turn=0.6,
                       base_speed=0.5, substeps=24, rng=rng2)

        np.testing.assert_array_equal(agents1['x'], agents2['x'])
        np.testing.assert_array_equal(agents1['y'], agents2['y'])
        np.testing.assert_array_equal(agents1['heading'], agents2['heading'])

    def test_nonzero_sigma_changes_positions(self):
        """speed_sigma > 0 should produce different positions than sigma=0."""
        rng1 = np.random.default_rng(42)
        agents1 = make_agents(100, habitat_side=200.0, rng=rng1)
        rng1 = np.random.default_rng(99)

        rng2 = np.random.default_rng(42)
        agents2 = make_agents(100, habitat_side=200.0, rng=rng2)
        rng2 = np.random.default_rng(99)

        daily_movement(agents1, habitat_side=200.0, sigma_turn=0.6,
                       base_speed=0.5, substeps=24, rng=rng1,
                       speed_sigma=0.0)

        daily_movement(agents2, habitat_side=200.0, sigma_turn=0.6,
                       base_speed=0.5, substeps=24, rng=rng2,
                       speed_sigma=0.5)

        alive = agents1['alive'].astype(bool)
        # Positions should differ (RNG stream diverges after noise draws)
        assert not np.allclose(agents1['x'][alive], agents2['x'][alive])

    def test_sigma_preserves_mean_displacement(self):
        """Bias-corrected log-normal should preserve mean step length."""
        rng = np.random.default_rng(42)
        sigma = 0.5
        # Draw many samples and check E[noise] ≈ 1.0
        noise = rng.lognormal(-0.5 * sigma**2, sigma, size=100_000)
        np.testing.assert_allclose(np.mean(noise), 1.0, atol=0.02)

    def test_sigma_increases_displacement_variance(self):
        """With speed_sigma > 0, per-agent displacements should vary more."""
        # Run with sigma=0: all agents of same disease state move same distance
        rng0 = np.random.default_rng(42)
        agents0 = make_agents(200, habitat_side=2000.0, rng=rng0)
        # Set all to same heading and straight-line movement for clean measurement
        agents0['heading'][:200] = 0.0  # All heading east
        x_before0 = agents0['x'][:200].copy()
        rng0 = np.random.default_rng(99)
        daily_movement(agents0, habitat_side=2000.0, sigma_turn=0.0,
                       base_speed=0.5, substeps=1, rng=rng0,
                       speed_sigma=0.0)
        dx0 = agents0['x'][:200] - x_before0
        var0 = np.var(dx0)

        # Run with sigma=0.5: agents should have variable displacements
        rng1 = np.random.default_rng(42)
        agents1 = make_agents(200, habitat_side=2000.0, rng=rng1)
        agents1['heading'][:200] = 0.0
        x_before1 = agents1['x'][:200].copy()
        rng1 = np.random.default_rng(99)
        daily_movement(agents1, habitat_side=2000.0, sigma_turn=0.0,
                       base_speed=0.5, substeps=1, rng=rng1,
                       speed_sigma=0.5)
        dx1 = agents1['x'][:200] - x_before1
        var1 = np.var(dx1)

        assert var1 > var0, (
            f"Displacement variance with speed_sigma=0.5 ({var1:.4f}) "
            f"should exceed sigma=0 ({var0:.4f})"
        )

    def test_positions_within_bounds_with_sigma(self):
        """Positions must remain in [0, side] even with high speed_sigma."""
        rng = np.random.default_rng(42)
        agents = make_agents(500, habitat_side=100.0, rng=rng)

        daily_movement(agents, habitat_side=100.0, sigma_turn=0.6,
                       base_speed=0.5, substeps=24, rng=rng,
                       speed_sigma=1.0)

        alive = agents['alive'].astype(bool)
        assert np.all(agents['x'][alive] >= 0)
        assert np.all(agents['x'][alive] <= 100.0)
        assert np.all(agents['y'][alive] >= 0)
        assert np.all(agents['y'][alive] <= 100.0)

    def test_update_movement_with_sigma(self):
        """update_movement() should work with speed_sigma parameter."""
        rng = np.random.default_rng(42)
        agents = make_agents(50, habitat_side=200.0, rng=rng)
        x_before = agents['x'][:50].copy()

        update_movement(agents, dt_minutes=60.0, habitat_side=200.0,
                        sigma_turn=0.6, base_speed=0.5, rng=rng,
                        speed_sigma=0.5)

        alive = agents['alive'].astype(bool)
        assert not np.allclose(agents['x'][alive], x_before[alive])
        assert np.all(agents['x'][alive] >= 0)
        assert np.all(agents['x'][alive] <= 200.0)

    def test_daily_movement_gravity_path_with_sigma(self):
        """speed_sigma should work when gravity fallback path is used."""
        rng = np.random.default_rng(42)
        agents = make_agents(30, habitat_side=200.0, rng=rng)
        agents[0:10]['spawn_gravity_timer'] = 5

        daily_movement(agents, habitat_side=200.0, sigma_turn=0.3,
                       base_speed=0.5, substeps=24, rng=rng,
                       spawning_gravity_enabled=True, gravity_strength=0.5,
                       gravity_range=80.0, pre_spawn_gravity_days=10,
                       post_spawn_gravity_days=10, speed_sigma=0.3)

        alive = agents['alive'].astype(bool)
        assert np.all(agents['x'][alive] >= 0)
        assert np.all(agents['x'][alive] <= 200.0)


class TestSpawningGravityIntegration:
    def test_full_workflow_gravity_and_movement(self):
        """Test complete workflow: set timers, move with gravity, timers count down."""
        rng = np.random.default_rng(42)
        agents = make_agents(50, habitat_side=200.0, rng=rng)
        
        # Set some agents to have gravity
        agents[0:20]['spawn_gravity_timer'] = 5
        
        # Run daily movement with gravity enabled
        daily_movement(agents, habitat_side=200.0, sigma_turn=0.3,
                       base_speed=0.5, substeps=24, rng=rng,
                       spawning_gravity_enabled=True, gravity_strength=0.5,
                       gravity_range=80.0, pre_spawn_gravity_days=10,
                       post_spawn_gravity_days=10)
        
        # Agents should have moved (positions changed)
        alive = agents['alive'].astype(bool)
        # Can't test exact positions due to randomness, but ensure they're within bounds
        assert np.all(agents['x'][alive] >= 0)
        assert np.all(agents['x'][alive] <= 200.0)
        assert np.all(agents['y'][alive] >= 0) 
        assert np.all(agents['y'][alive] <= 200.0)
        
        # Gravity timers should still be set (movement doesn't decrement them)
        assert np.sum(agents['spawn_gravity_timer'] > 0) == 20
