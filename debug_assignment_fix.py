#!/usr/bin/env python3
import sys; sys.path.insert(0, '.')
import numpy as np
from sswd_evoepi.types import allocate_agents, allocate_genotypes, Stage, DiseaseState

# Test different ways to fix the batch assignment
agents = allocate_agents(100)
rng = np.random.default_rng(42)

dead_slots = np.where(~agents['alive'])[0]
n_slots = 12
slots = dead_slots[:n_slots]

print("Testing different assignment methods:")

# Method 1: Original (fails)
agents1 = agents.copy()
agents1[slots]['alive'] = True
success1 = np.sum(agents1[slots]['alive'])
print(f"Method 1 (original): {success1}/{n_slots} successful")

# Method 2: Using a loop (should work)
agents2 = agents.copy()
for slot in slots:
    agents2[slot]['alive'] = True
success2 = np.sum(agents2[slots]['alive'])
print(f"Method 2 (loop): {success2}/{n_slots} successful")

# Method 3: Assign to field directly
agents3 = agents.copy()
agents3['alive'][slots] = True
success3 = np.sum(agents3[slots]['alive'])
print(f"Method 3 (field direct): {success3}/{n_slots} successful")

# Method 4: Using np.put
agents4 = agents.copy()
np.put(agents4['alive'], slots, True)
success4 = np.sum(agents4[slots]['alive'])
print(f"Method 4 (np.put): {success4}/{n_slots} successful")

print(f"\nRecommended fix: Use agents['alive'][slots] = True instead of agents[slots]['alive'] = True")