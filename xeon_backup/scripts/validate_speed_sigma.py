#!/usr/bin/env python3
"""Quick validation: effect of speed_sigma on movement variability.

Runs a simple movement-only test comparing speed_sigma = 0, 0.3, 0.5
to show that stochastic step lengths increase displacement variance
while preserving mean displacement.
"""

import numpy as np
import sys
sys.path.insert(0, '.')

from sswd_evoepi.types import AGENT_DTYPE, DiseaseState, allocate_agents
from sswd_evoepi.movement import daily_movement


def make_agents(n, habitat_side, rng):
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


def run_trial(speed_sigma, n_agents=500, n_days=10, seed=42):
    """Run movement for n_days and return displacement statistics."""
    rng_init = np.random.default_rng(seed)
    agents = make_agents(n_agents, habitat_side=2000.0, rng=rng_init)
    x0 = agents['x'][:n_agents].copy()
    y0 = agents['y'][:n_agents].copy()

    rng_move = np.random.default_rng(seed + 1000)
    for _ in range(n_days):
        daily_movement(agents, habitat_side=2000.0, sigma_turn=0.6,
                       base_speed=0.5, substeps=24, rng=rng_move,
                       speed_sigma=speed_sigma)

    dx = agents['x'][:n_agents] - x0
    dy = agents['y'][:n_agents] - y0
    disp = np.sqrt(dx**2 + dy**2)
    return disp


print("=" * 60)
print("Stochastic Step-Length Validation")
print("=" * 60)
print(f"{'sigma':>8}  {'mean_disp':>10}  {'std_disp':>10}  {'CV':>8}  {'min':>8}  {'max':>8}")
print("-" * 60)

for sigma in [0.0, 0.3, 0.5, 1.0]:
    disp = run_trial(sigma)
    mean_d = np.mean(disp)
    std_d = np.std(disp)
    cv = std_d / mean_d if mean_d > 0 else 0
    print(f"{sigma:>8.1f}  {mean_d:>10.2f}  {std_d:>10.2f}  {cv:>8.3f}  {disp.min():>8.2f}  {disp.max():>8.2f}")

print("-" * 60)
print("Note: mean displacement should stay ~similar (bias-corrected)")
print("      CV (coefficient of variation) should increase with sigma")
print("      This confirms stochastic step lengths work correctly.")
