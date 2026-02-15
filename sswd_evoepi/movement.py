"""Agent movement via correlated random walk (CRW).

Pycnopodia helianthoides movement:
  - Speed: 0.1–2.0 m/min (Kay & Emlet 2002)
  - Pattern: correlated random walk with directional persistence
  - Disease-dependent: I1 ×0.5, I2 ×0.1, D ×0.0 of base speed

The CRW updates heading with wrapped-normal turning noise:
    heading += Normal(0, sigma_turn)
    x += speed × cos(heading) × dt
    y += speed × sin(heading) × dt

Boundary: elastic reflection off square habitat [0, side] × [0, side].

Spatial transmission: optional grid-based local Vibrio exposure.
Instead of every susceptible seeing the node-average P_k, each sees
P_k scaled by the local density of infected agents relative to the
node mean.  This creates emergent spatial clustering of disease
without changing the Vibrio ODE.

References:
  - modeling-approach.md §2.4 (1-hour substeps, 24/day)
  - disease-module-spec §5.5 (speed modifiers)
  - integration-architecture-spec §2.2 (originally "disabled in v1")
"""

from __future__ import annotations

import numpy as np
from typing import Optional, Tuple

from .types import DiseaseState

# ═══════════════════════════════════════════════════════════════════════
# CONSTANTS
# ═══════════════════════════════════════════════════════════════════════

# Speed modifiers by disease state (disease-module-spec §5.5)
SPEED_MODIFIER = np.array([1.0, 1.0, 0.5, 0.1, 0.0, 1.0], dtype=np.float32)
#                          S     E    I1    I2    D    R

TWO_PI = 2.0 * np.pi


# ═══════════════════════════════════════════════════════════════════════
# BOUNDARY REFLECTION
# ═══════════════════════════════════════════════════════════════════════

def _reflect(x: np.ndarray, side: float) -> np.ndarray:
    """Elastic reflection on [0, side].

    Uses modular arithmetic to handle arbitrary displacements
    (multiple bounces in one step, though unlikely at normal speeds).
    """
    period = 2.0 * side
    # Map to [0, period)
    x = x % period
    # Reflect the upper half
    x = np.where(x > side, period - x, x)
    return x


# ═══════════════════════════════════════════════════════════════════════
# CORRELATED RANDOM WALK
# ═══════════════════════════════════════════════════════════════════════

def update_movement(
    agents: np.ndarray,
    dt_minutes: float,
    habitat_side: float,
    sigma_turn: float,
    base_speed: float,
    rng: np.random.Generator,
    spawning_gravity_enabled: bool = False,
    gravity_strength: float = 0.3,
    gravity_range: float = 100.0,
    pre_spawn_gravity_days: int = 14,
    post_spawn_gravity_days: int = 14,
) -> None:
    """Move all alive agents one CRW sub-step (in-place).

    Args:
        agents: Structured array with AGENT_DTYPE fields.
        dt_minutes: Timestep in minutes (typically 60 for hourly).
        habitat_side: Side length of square habitat (m).
        sigma_turn: Turning angle standard deviation (radians).
            0 = straight lines, π = uniform random turns.
        base_speed: Base undisturbed movement speed (m/min).
        rng: NumPy random generator.
        spawning_gravity_enabled: Enable spawning gravity aggregation.
        gravity_strength: Maximum speed bias toward conspecifics (m/min).
        gravity_range: Sensory detection range for conspecifics (m).
        pre_spawn_gravity_days: Days before spawning that gravity activates.
        post_spawn_gravity_days: Days after spawning that gravity persists.
    """
    alive_mask = agents['alive'].astype(bool)
    n_alive = int(alive_mask.sum())
    if n_alive == 0:
        return

    alive_idx = np.where(alive_mask)[0]

    # 1. Update headings: add wrapped-normal turning noise
    turn_angles = rng.normal(0.0, sigma_turn, size=n_alive).astype(np.float32)
    agents['heading'][alive_idx] += turn_angles
    agents['heading'][alive_idx] %= TWO_PI

    # 1.5. Apply spawning gravity bias (Phase 3)
    if spawning_gravity_enabled:
        agents['heading'][alive_idx] = _apply_spawning_gravity(
            agents, alive_idx, gravity_strength, gravity_range, 
            pre_spawn_gravity_days, post_spawn_gravity_days
        )

    # 2. Compute effective speed: base × disease modifier
    disease_states = agents['disease_state'][alive_idx]
    speed_mods = SPEED_MODIFIER[disease_states]
    effective_speed = base_speed * speed_mods
    agents['speed'][alive_idx] = effective_speed

    # 3. Compute displacement
    headings = agents['heading'][alive_idx]
    dx = effective_speed * np.cos(headings) * dt_minutes
    dy = effective_speed * np.sin(headings) * dt_minutes

    # 4. Update positions with boundary reflection
    new_x = agents['x'][alive_idx] + dx
    new_y = agents['y'][alive_idx] + dy

    agents['x'][alive_idx] = _reflect(new_x, habitat_side)
    agents['y'][alive_idx] = _reflect(new_y, habitat_side)


def daily_movement(
    agents: np.ndarray,
    habitat_side: float,
    sigma_turn: float,
    base_speed: float,
    substeps: int,
    rng: np.random.Generator,
    spawning_gravity_enabled: bool = False,
    gravity_strength: float = 0.3,
    gravity_range: float = 100.0,
    pre_spawn_gravity_days: int = 14,
    post_spawn_gravity_days: int = 14,
) -> None:
    """Run one full day of CRW movement (multiple sub-steps).

    Args:
        agents: Structured array with AGENT_DTYPE fields.
        habitat_side: Side length of square habitat (m).
        sigma_turn: Turning angle std dev (radians).
        base_speed: Base movement speed (m/min).
        substeps: Number of sub-steps per day (24 = hourly).
        rng: NumPy random generator.
        spawning_gravity_enabled: Enable spawning gravity aggregation.
        gravity_strength: Maximum speed bias toward conspecifics (m/min).
        gravity_range: Sensory detection range for conspecifics (m).
        pre_spawn_gravity_days: Days before spawning that gravity activates.
        post_spawn_gravity_days: Days after spawning that gravity persists.
    """
    dt_minutes = (24.0 * 60.0) / substeps  # minutes per sub-step
    for _ in range(substeps):
        update_movement(agents, dt_minutes, habitat_side, sigma_turn,
                        base_speed, rng, spawning_gravity_enabled, 
                        gravity_strength, gravity_range, 
                        pre_spawn_gravity_days, post_spawn_gravity_days)


# ═══════════════════════════════════════════════════════════════════════
# SPAWNING GRAVITY (Phase 3)
# ═══════════════════════════════════════════════════════════════════════

def _apply_spawning_gravity(
    agents: np.ndarray,
    alive_idx: np.ndarray,
    gravity_strength: float,
    gravity_range: float,
    pre_spawn_gravity_days: int,
    post_spawn_gravity_days: int,
) -> np.ndarray:
    """Apply spawning gravity bias to agent headings.
    
    For agents with spawn_gravity_timer > 0, add a bias vector toward 
    the center-of-mass of conspecifics within gravity_range.
    Gravity strength ramps linearly over the spawning window.
    
    Args:
        agents: Full agent array.
        alive_idx: Indices of alive agents.
        gravity_strength: Maximum speed bias toward conspecifics (m/min).
        gravity_range: Sensory detection range for conspecifics (m).
        pre_spawn_gravity_days: Days before spawning that gravity activates.
        post_spawn_gravity_days: Days after spawning that gravity persists.
        
    Returns:
        Modified headings array for alive agents.
    """
    headings = agents['heading'][alive_idx].copy()
    
    # Get agents with active gravity (spawn_gravity_timer > 0)
    gravity_timers = agents['spawn_gravity_timer'][alive_idx]
    gravity_mask = gravity_timers > 0
    
    if not np.any(gravity_mask):
        return headings  # No agents with active gravity
    
    # Get positions of all alive agents and agents with gravity
    all_positions = np.column_stack([agents['x'][alive_idx], agents['y'][alive_idx]])
    gravity_indices_local = np.where(gravity_mask)[0]  # Indices within alive_idx
    
    total_gravity_days = pre_spawn_gravity_days + post_spawn_gravity_days
    
    # Process each agent with active gravity
    for i in gravity_indices_local:
        agent_pos = all_positions[i]
        timer = gravity_timers[i]
        
        # Calculate gravity strength ramp (linear)
        # Timer counts down: total_days → 0
        # Days from start: 0 → total_days  
        days_from_start = total_gravity_days - timer
        
        if days_from_start <= pre_spawn_gravity_days:
            # Pre-spawn: ramp up from 0 to max
            ramp_factor = days_from_start / pre_spawn_gravity_days
        else:
            # Post-spawn: ramp down from max to 0
            days_into_post_spawn = days_from_start - pre_spawn_gravity_days
            ramp_factor = 1.0 - (days_into_post_spawn / post_spawn_gravity_days)
        
        ramp_factor = np.clip(ramp_factor, 0.0, 1.0)
        effective_gravity = gravity_strength * ramp_factor
        
        if effective_gravity <= 0:
            continue
        
        # Find neighbors within gravity range
        distances = np.sqrt(np.sum((all_positions - agent_pos)**2, axis=1))
        neighbor_mask = (distances <= gravity_range) & (distances > 0)  # Exclude self
        
        if not np.any(neighbor_mask):
            continue  # No neighbors within range
        
        # Calculate center of mass of neighbors
        neighbor_positions = all_positions[neighbor_mask]
        com = np.mean(neighbor_positions, axis=0)
        
        # Calculate angle to center of mass
        dx = com[0] - agent_pos[0]
        dy = com[1] - agent_pos[1]
        
        if dx == 0 and dy == 0:
            continue  # Agent at center of mass
            
        angle_to_com = np.arctan2(dy, dx)
        
        # Blend gravity direction with current heading
        current_heading = headings[i]
        blended_heading = _circular_blend(current_heading, angle_to_com, effective_gravity)
        headings[i] = blended_heading
    
    return headings


def _circular_blend(heading1: float, heading2: float, weight: float) -> float:
    """Blend two circular angles (headings) with given weight.
    
    Properly handles the circular nature of angles (e.g., blending 0.1 and 6.2 
    should go through 0, not through π).
    
    Args:
        heading1: First heading (radians).
        heading2: Second heading (radians).
        weight: Blend weight for heading2 [0, 1]. 0 = pure heading1, 1 = pure heading2.
        
    Returns:
        Blended heading (radians).
    """
    # Convert to complex numbers on unit circle
    z1 = np.exp(1j * heading1)
    z2 = np.exp(1j * heading2)
    
    # Linear interpolation in complex plane
    z_blend = (1 - weight) * z1 + weight * z2
    
    # Convert back to angle
    blended = np.angle(z_blend)
    
    # Ensure positive angle [0, 2π)
    if blended < 0:
        blended += TWO_PI
        
    return blended


# ═══════════════════════════════════════════════════════════════════════
# SPATIAL TRANSMISSION GRID
# ═══════════════════════════════════════════════════════════════════════

class InfectedDensityGrid:
    """Grid-based local infected density for spatial transmission.

    Divides the square habitat into cells and computes infected density
    per cell.  Used to scale P_k so that susceptibles near clusters of
    infected see higher local Vibrio exposure.

    The local exposure factor for susceptible i is:
        f_i = local_infected_density / mean_infected_density

    When f_i > 1: more infected nearby than average → higher exposure.
    When f_i < 1: fewer infected nearby → lower exposure.
    Mean(f_i) ≈ 1 across the node, preserving the average FOI.

    Grid is rebuilt each day (before transmission step).
    """

    def __init__(self, habitat_side: float, cell_size: float = 20.0):
        """Initialize grid.

        Args:
            habitat_side: Side length of square habitat (m).
            cell_size: Grid cell size (m).  Smaller = finer spatial
                resolution but noisier.  20m ≈ 10 arm-span of adult
                Pycnopodia, reasonable contact scale.
        """
        self.habitat_side = habitat_side
        self.cell_size = cell_size
        self.nx = max(1, int(np.ceil(habitat_side / cell_size)))
        self.ny = self.nx
        self.grid = np.zeros((self.nx, self.ny), dtype=np.float64)

    def build(
        self,
        agents: np.ndarray,
        diffusion_passes: int = 2,
    ) -> None:
        """Build infected density grid from current agent positions.

        Args:
            agents: Structured array with AGENT_DTYPE.
            diffusion_passes: Number of 3×3 averaging passes.
                0 = no smoothing (very noisy).
                2 = moderate smoothing (~40m effective radius).
        """
        self.grid[:] = 0.0

        alive = agents['alive'].astype(bool)
        ds = agents['disease_state']

        # Infected = I1 or I2 (actively shedding Vibrio)
        infected_mask = alive & ((ds == DiseaseState.I1) | (ds == DiseaseState.I2))
        n_infected = int(infected_mask.sum())

        if n_infected == 0:
            # No infected → uniform (factor = 1 everywhere)
            self.grid[:] = 1.0
            return

        # Bin infected agents into grid cells
        ix = agents['x'][infected_mask]
        iy = agents['y'][infected_mask]
        cx = np.clip((ix / self.cell_size).astype(int), 0, self.nx - 1)
        cy = np.clip((iy / self.cell_size).astype(int), 0, self.ny - 1)

        np.add.at(self.grid, (cx, cy), 1.0)

        # Smooth with diffusion (3×3 averaging)
        for _ in range(diffusion_passes):
            self.grid = _diffuse_2d(self.grid)

        # Normalize: mean grid value = 1.0
        mean_val = self.grid.mean()
        if mean_val > 0:
            self.grid /= mean_val
        else:
            self.grid[:] = 1.0

    def lookup(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        """Look up local exposure factor for agent positions.

        Args:
            x: Agent x-positions (m), array of shape (n,).
            y: Agent y-positions (m), array of shape (n,).

        Returns:
            Array of shape (n,) with local exposure factors.
            Mean ≈ 1.0 across all agents.
        """
        cx = np.clip((x / self.cell_size).astype(int), 0, self.nx - 1)
        cy = np.clip((y / self.cell_size).astype(int), 0, self.ny - 1)
        return self.grid[cx, cy]


def _diffuse_2d(grid: np.ndarray) -> np.ndarray:
    """One step of 3×3 averaging (diffusion on 2D grid).

    Uses reflect-padding at boundaries.  Pure NumPy, no scipy.
    """
    nx, ny = grid.shape
    # Pad with reflection
    padded = np.pad(grid, 1, mode='reflect')
    result = np.zeros_like(grid)
    for di in range(3):
        for dj in range(3):
            result += padded[di:di + nx, dj:dj + ny]
    return result / 9.0
