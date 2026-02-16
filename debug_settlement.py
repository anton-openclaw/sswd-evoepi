#!/usr/bin/env python3
import sys; sys.path.insert(0, '.')
import numpy as np
from sswd_evoepi.types import allocate_agents, allocate_genotypes, Stage, DiseaseState
from sswd_evoepi.reproduction import settle_recruits
from sswd_evoepi.genetics import initialize_effect_sizes

# Recreate the failing test setup
agents = allocate_agents(100)
genotypes = allocate_genotypes(100)
rng = np.random.default_rng(42)

# Arriving settlers (12 competent larvae)
n_settlers = 12
settler_geno = rng.integers(0, 2, (n_settlers, 52, 2), dtype=np.int8)

# Remove lethal EF1A homozygotes
IDX_EF1A = 51  # Last locus
ef1a_sums = settler_geno[:, IDX_EF1A, :].sum(axis=1)
lethal = ef1a_sums == 2
print(f"Lethal EF1A homozygotes: {np.sum(lethal)} out of {n_settlers}")
settler_geno[lethal, IDX_EF1A, 1] = 0  # Fix them

effect_sizes_arr = initialize_effect_sizes(rng)

print(f"Before settlement: {np.sum(agents['alive'])} alive")
print(f"Dead slots available: {len(np.where(~agents['alive'])[0])}")

n_added = settle_recruits(
    agents=agents,
    genotypes=genotypes,
    settler_genotypes=settler_geno,
    node_id=0,
    carrying_capacity=500,
    n_adults_present=20,
    habitat_area=10000.0,
    effect_sizes=effect_sizes_arr,
    rng=rng,
)

print(f"After settlement:")
print(f"  n_added returned: {n_added}")
print(f"  agents alive: {np.sum(agents['alive'])}")
print(f"  stages: {np.unique(agents['stage'][agents['alive']], return_counts=True)}")

# Debug the specific agents that should be settlers
if n_added > 0:
    alive_indices = np.where(agents['alive'])[0]
    print(f"  alive indices: {alive_indices}")
    print(f"  their stages: {agents['stage'][alive_indices]}")