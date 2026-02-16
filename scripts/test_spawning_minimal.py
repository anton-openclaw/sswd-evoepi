#!/usr/bin/env python3
"""
Minimal test of spawning functionality to debug issues.
"""

import numpy as np
from sswd_evoepi.types import AGENT_DTYPE, Stage, N_LOCI, allocate_agents, allocate_genotypes
from sswd_evoepi.config import SpawningSection, DiseaseSection
from sswd_evoepi.spawning import spawning_step
from sswd_evoepi.genetics import initialize_genotypes, initialize_effect_sizes

def test_minimal():
    print("Creating minimal test population...")
    
    rng = np.random.default_rng(42)
    n_agents = 10
    
    # Create agents
    agents = allocate_agents(n_agents)
    genotypes = allocate_genotypes(n_agents)
    
    # Initialize genetics
    effect_sizes = initialize_effect_sizes(rng)
    geno_init = initialize_genotypes(n_agents, effect_sizes, rng)
    genotypes[:n_agents] = geno_init
    
    print(f"Created {n_agents} agents and genotypes")
    
    # Set up agents
    for i in range(n_agents):
        agents[i]['node_id'] = 1
        agents[i]['x'] = rng.uniform(0, 100)
        agents[i]['y'] = rng.uniform(0, 100)
        agents[i]['stage'] = Stage.ADULT
        agents[i]['age'] = 1000
        agents[i]['size'] = 450.0
        agents[i]['sex'] = i % 2  # Alternating sex
        agents[i]['alive'] = True
        agents[i]['spawning_ready'] = False
        agents[i]['has_spawned'] = 0
        agents[i]['spawn_refractory'] = 0
        agents[i]['last_spawn_day'] = -999
    
    print(f"Set up agents: {np.sum(agents['alive'])} alive, {np.sum(agents['sex'] == 0)} female, {np.sum(agents['sex'] == 1)} male")
    
    # Create configs
    spawning_config = SpawningSection()
    spawning_config.p_spontaneous_female = 0.1  # High probability for testing
    spawning_config.p_spontaneous_male = 0.1
    
    disease_config = DiseaseSection()
    
    print("Running spawning step...")
    
    # Try one spawning step
    try:
        cohorts = spawning_step(
            agents, genotypes, 105,  # Peak spawning day
            48.5,  # Latitude
            spawning_config, disease_config, rng
        )
        
        print(f"Spawning step successful!")
        print(f"  Cohorts produced: {len(cohorts)}")
        print(f"  Total larvae: {sum(cohort.n_competent for cohort in cohorts)}")
        print(f"  Spawning events: {np.sum(agents['has_spawned'])}")
        print(f"  Ready agents: {np.sum(agents['spawning_ready'])}")
        
        # Print agent states
        for i in range(n_agents):
            print(f"  Agent {i}: sex={agents[i]['sex']}, ready={agents[i]['spawning_ready']}, spawned={agents[i]['has_spawned']}")
        
        return True
        
    except Exception as e:
        print(f"ERROR in spawning step: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = test_minimal()
    print(f"Test {'PASSED' if success else 'FAILED'}")