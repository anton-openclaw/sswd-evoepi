"""Timing comparison script for spawning vectorization - AFTER optimization.

Tests performance of spawning system after vectorization.
"""

import time
import numpy as np
from sswd_evoepi.spawning import spawning_step
from sswd_evoepi.types import AGENT_DTYPE, Stage
from sswd_evoepi.config import SpawningSection, DiseaseSection

def create_test_agents(n_agents: int, node_latitude: float, rng: np.random.Generator) -> tuple:
    """Create test agent population for timing."""
    agents = np.zeros(n_agents, dtype=AGENT_DTYPE)
    genotypes = np.random.randint(0, 2, size=(n_agents, 52, 2), dtype=np.int8)
    
    # Set up test agents
    agents['alive'] = True
    agents['stage'] = Stage.ADULT
    agents['node_id'] = 1
    agents['x'] = rng.uniform(-1000, 1000, n_agents)  # Spread over 2km x 2km area
    agents['y'] = rng.uniform(-1000, 1000, n_agents)
    agents['sex'] = rng.choice([0, 1], n_agents)
    agents['size'] = rng.uniform(450, 600, n_agents)  # All reproductive size
    agents['spawning_ready'] = rng.choice([0, 1], n_agents, p=[0.3, 0.7])  # 70% ready
    agents['has_spawned'] = 0
    agents['spawn_refractory'] = 0
    agents['spawn_gravity_timer'] = 0
    agents['immunosuppression_timer'] = 0
    agents['last_spawn_day'] = rng.choice([0, 100, 105, 110], n_agents, p=[0.4, 0.2, 0.2, 0.2])  # Some recent spawners
    
    return agents, genotypes

def run_spawning_timing(n_agents: int, n_days: int, description: str):
    """Run spawning simulation and time it."""
    print(f"\n=== {description} ===")
    print(f"Agents: {n_agents}, Days: {n_days}")
    
    # Setup
    rng = np.random.Generator(np.random.PCG64(42))  # Fixed seed
    agents, genotypes = create_test_agents(n_agents, 48.5, rng)
    
    # Config
    spawning_config = SpawningSection()
    spawning_config.season_start_doy = 305
    spawning_config.season_end_doy = 196
    spawning_config.p_spontaneous_female = 0.05
    spawning_config.p_spontaneous_male = 0.05
    spawning_config.cascade_window = 14
    spawning_config.cascade_radius = 200.0
    spawning_config.induction_female_to_male = 0.80
    spawning_config.induction_male_to_female = 0.30
    spawning_config.male_max_bouts = 3
    spawning_config.male_refractory_days = 7
    spawning_config.peak_doy = 105
    spawning_config.peak_width_days = 45
    spawning_config.lat_shift_per_deg = 3.0
    spawning_config.gravity_enabled = False  # Simplify for timing
    
    disease_config = DiseaseSection()
    disease_config.immunosuppression_enabled = False  # Simplify for timing
    
    # Time the spawning steps
    start_time = time.time()
    total_cohorts = 0
    spawning_times = []
    
    for day in range(100, 100 + n_days):  # Peak spawning season
        step_start = time.time()
        cohorts = spawning_step(
            agents, genotypes, day, 48.5, spawning_config, disease_config, rng
        )
        step_time = time.time() - step_start
        spawning_times.append(step_time)
        total_cohorts += len(cohorts)
        
        # Progress indicator
        if day % 10 == 0:
            print(f"  Day {day}: {step_time:.4f}s, {len(cohorts)} cohorts")
    
    total_time = time.time() - start_time
    avg_step_time = np.mean(spawning_times)
    
    print(f"Total time: {total_time:.3f}s")
    print(f"Average step time: {avg_step_time:.4f}s")
    print(f"Total cohorts: {total_cohorts}")
    
    return total_time, avg_step_time

if __name__ == "__main__":
    print("AFTER VECTORIZATION - Timing Test")
    time_100_20 = run_spawning_timing(100, 20, "100 agents, 20 days")
    time_200_10 = run_spawning_timing(200, 10, "200 agents, 10 days")
    
    print(f"\n=== AFTER VECTORIZATION SUMMARY ===")
    print(f"100 agents, 20 days: {time_100_20[0]:.3f}s total, {time_100_20[1]:.4f}s/step")
    print(f"200 agents, 10 days: {time_200_10[0]:.3f}s total, {time_200_10[1]:.4f}s/step")
    
    # Comparison with baseline
    baseline_100_20 = (18.669, 0.9335)
    baseline_200_10 = (16.703, 1.6703)
    
    speedup_100_total = baseline_100_20[0] / time_100_20[0]
    speedup_100_step = baseline_100_20[1] / time_100_20[1]
    speedup_200_total = baseline_200_10[0] / time_200_10[0]
    speedup_200_step = baseline_200_10[1] / time_200_10[1]
    
    print(f"\n=== SPEEDUP ANALYSIS ===")
    print(f"100 agents, 20 days: {speedup_100_total:.2f}× total speedup, {speedup_100_step:.2f}× step speedup")
    print(f"200 agents, 10 days: {speedup_200_total:.2f}× total speedup, {speedup_200_step:.2f}× step speedup")
    print(f"Average speedup: {(speedup_100_step + speedup_200_step)/2:.2f}× per step")