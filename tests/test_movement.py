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
