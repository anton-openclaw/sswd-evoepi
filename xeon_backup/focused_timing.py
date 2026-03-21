"""Focused timing test for vectorized spawning functions."""

import time
import numpy as np
from sswd_evoepi.spawning import _get_recent_spawners_mask, _check_cascade_induction, spawning_step
from sswd_evoepi.types import AGENT_DTYPE, Stage
from sswd_evoepi.config import SpawningSection, DiseaseSection

def time_recent_spawners_mask(n_agents: int, n_iterations: int = 1000):
    """Time the recent spawners mask function."""
    print(f"\nTiming _get_recent_spawners_mask with {n_agents} agents, {n_iterations} iterations...")
    
    # Create test data
    agents = np.zeros(n_agents, dtype=AGENT_DTYPE)
    agents['last_spawn_day'] = np.random.choice([0, 100, 105, 110, 350], n_agents)
    
    current_doy = 108
    cascade_window = 14
    
    # Time the function
    start_time = time.time()
    for _ in range(n_iterations):
        mask = _get_recent_spawners_mask(agents, current_doy, cascade_window)
    end_time = time.time()
    
    total_time = end_time - start_time
    avg_time = total_time / n_iterations
    
    print(f"  Total time: {total_time:.4f}s")
    print(f"  Average per call: {avg_time*1000:.4f}ms")
    print(f"  Recent spawners found: {np.sum(mask)}")
    
    return avg_time

def time_cascade_induction(n_targets: int, n_inducers: int, n_iterations: int = 100):
    """Time the cascade induction function."""
    print(f"\nTiming _check_cascade_induction with {n_targets} targets, {n_inducers} inducers, {n_iterations} iterations...")
    
    # Create test agents spread over 2km x 2km area
    n_agents = n_targets + n_inducers + 100  # Ensure enough agents
    agents = np.zeros(n_agents, dtype=AGENT_DTYPE)
    agents['x'] = np.random.uniform(-1000, 1000, n_agents)
    agents['y'] = np.random.uniform(-1000, 1000, n_agents)
    
    target_indices = np.arange(n_targets)
    inducer_indices = np.arange(n_targets, n_targets + n_inducers)
    cascade_radius = 200.0
    induction_probability = 0.8
    rng = np.random.Generator(np.random.PCG64(42))
    
    # Time the function
    start_time = time.time()
    total_induced = 0
    for _ in range(n_iterations):
        induced = _check_cascade_induction(
            agents, target_indices, inducer_indices,
            cascade_radius, induction_probability, rng
        )
        total_induced += len(induced)
    end_time = time.time()
    
    total_time = end_time - start_time
    avg_time = total_time / n_iterations
    
    print(f"  Total time: {total_time:.4f}s")
    print(f"  Average per call: {avg_time*1000:.4f}ms")
    print(f"  Average induced per call: {total_induced/n_iterations:.1f}")
    print(f"  T×I complexity: {n_targets * n_inducers:,}")
    
    return avg_time

def time_spawning_step_simple(n_agents: int, n_steps: int = 10):
    """Time a simplified spawning step."""
    print(f"\nTiming spawning_step with {n_agents} agents, {n_steps} steps...")
    
    # Setup
    rng = np.random.Generator(np.random.PCG64(42))
    agents = np.zeros(n_agents, dtype=AGENT_DTYPE)
    genotypes = np.random.randint(0, 2, size=(n_agents, 52, 2), dtype=np.int8)
    
    # Set up test agents
    agents['alive'] = True
    agents['stage'] = Stage.ADULT
    agents['node_id'] = 1
    agents['x'] = rng.uniform(-1000, 1000, n_agents)
    agents['y'] = rng.uniform(-1000, 1000, n_agents)
    agents['sex'] = rng.choice([0, 1], n_agents)
    agents['size'] = rng.uniform(450, 600, n_agents)
    agents['spawning_ready'] = rng.choice([0, 1], n_agents, p=[0.3, 0.7])
    agents['has_spawned'] = 0
    agents['spawn_refractory'] = 0
    agents['spawn_gravity_timer'] = 0
    agents['immunosuppression_timer'] = 0
    agents['last_spawn_day'] = rng.choice([0, 100, 105, 110], n_agents, p=[0.4, 0.2, 0.2, 0.2])
    
    # Config
    spawning_config = SpawningSection()
    spawning_config.season_start_doy = 305
    spawning_config.season_end_doy = 196
    spawning_config.p_spontaneous_female = 0.02  # Lower to reduce larval generation overhead
    spawning_config.p_spontaneous_male = 0.02
    spawning_config.cascade_window = 14
    spawning_config.cascade_radius = 200.0
    spawning_config.induction_female_to_male = 0.80
    spawning_config.induction_male_to_female = 0.30
    spawning_config.male_max_bouts = 3
    spawning_config.male_refractory_days = 7
    spawning_config.peak_doy = 105
    spawning_config.peak_width_days = 45
    spawning_config.lat_shift_per_deg = 3.0
    spawning_config.gravity_enabled = False
    
    disease_config = DiseaseSection()
    disease_config.immunosuppression_enabled = False
    
    # Time the spawning steps
    start_time = time.time()
    total_cohorts = 0
    
    for day in range(105, 105 + n_steps):  # Peak spawning season
        cohorts = spawning_step(
            agents, genotypes, day, 48.5, spawning_config, disease_config, rng
        )
        total_cohorts += len(cohorts)
    
    total_time = time.time() - start_time
    avg_time = total_time / n_steps
    
    print(f"  Total time: {total_time:.4f}s")
    print(f"  Average per step: {avg_time*1000:.2f}ms")
    print(f"  Total cohorts: {total_cohorts}")
    
    return avg_time

if __name__ == "__main__":
    print("Focused Timing Test for Vectorized Spawning Functions")
    print("=" * 60)
    
    # Test individual function performance
    mask_time_100 = time_recent_spawners_mask(100, 1000)
    mask_time_1000 = time_recent_spawners_mask(1000, 1000)
    
    induction_time_small = time_cascade_induction(50, 50, 100)  # 2,500 pairs
    induction_time_medium = time_cascade_induction(100, 100, 50)  # 10,000 pairs
    induction_time_large = time_cascade_induction(200, 200, 20)  # 40,000 pairs
    
    # Test full spawning step performance (focused)
    step_time_100 = time_spawning_step_simple(100, 5)
    step_time_300 = time_spawning_step_simple(300, 5)
    
    print("\n" + "=" * 60)
    print("PERFORMANCE SUMMARY:")
    print(f"_get_recent_spawners_mask: {mask_time_100*1000:.4f}ms (100 agents) -> {mask_time_1000*1000:.4f}ms (1000 agents)")
    print(f"_check_cascade_induction: {induction_time_small*1000:.2f}ms (50×50) -> {induction_time_large*1000:.2f}ms (200×200)")
    print(f"spawning_step: {step_time_100*1000:.2f}ms (100 agents) -> {step_time_300*1000:.2f}ms (300 agents)")
    
    # Estimate speedup based on bottleneck elimination
    print(f"\nESTIMATED IMPROVEMENT:")
    print(f"- Recent spawner mask is now O(N) vectorized instead of O(N) Python loop")
    print(f"- Cascade induction is now O(T×I) vectorized instead of O(T×I) nested Python loops")
    print(f"- Distance calculations skip sqrt and use broadcasting")
    print(f"- Expected speedup: 2-5× for cascade-heavy scenarios")