#!/usr/bin/env python3
import sys; sys.path.insert(0, '.')
import numpy as np
from sswd_evoepi.types import allocate_agents, allocate_genotypes, Stage, DiseaseState
from sswd_evoepi.genetics import initialize_effect_sizes

# Test the agent assignment logic directly
agents = allocate_agents(100)
genotypes = allocate_genotypes(100)
rng = np.random.default_rng(42)

print("Before assignment:")
print(f"  alive agents: {np.sum(agents['alive'])}")
print(f"  first 5 alive values: {agents['alive'][:5]}")

# Find dead slots (should be all of them initially)
dead_slots = np.where(~agents['alive'])[0]
print(f"  dead_slots length: {len(dead_slots)}")
print(f"  first 5 dead slots: {dead_slots[:5]}")

# Select some slots to assign
n_slots = 12
slots = dead_slots[:n_slots]
print(f"  slots to assign: {slots}")

# Test individual assignment
print("\nTesting individual assignment:")
agents[slots[0]]['alive'] = True
print(f"  After setting slot {slots[0]}: alive = {agents[slots[0]]['alive']}")
print(f"  Total alive: {np.sum(agents['alive'])}")

# Reset for batch test
agents[slots[0]]['alive'] = False

# Test batch assignment
print("\nTesting batch assignment:")
print(f"  Before batch: alive sum = {np.sum(agents['alive'])}")
agents[slots]['alive'] = True
print(f"  After batch: alive sum = {np.sum(agents['alive'])}")
print(f"  First few alive values: {agents['alive'][:15]}")
print(f"  Alive values at assigned slots: {agents[slots]['alive']}")

# Test if the assignment actually worked
successful_assignments = np.sum(agents[slots]['alive'])
print(f"  Successful assignments: {successful_assignments} out of {n_slots}")

if successful_assignments != n_slots:
    print("  ASSIGNMENT FAILED!")
    # Debug each slot individually
    for i, slot in enumerate(slots):
        print(f"    slot {slot}: alive = {agents[slot]['alive']}")
else:
    print("  Assignment successful!")