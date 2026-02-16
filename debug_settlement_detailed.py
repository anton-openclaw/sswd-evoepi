#!/usr/bin/env python3
import sys; sys.path.insert(0, '.')
import numpy as np
from sswd_evoepi.types import allocate_agents, allocate_genotypes, Stage, DiseaseState
from sswd_evoepi.reproduction import settlement_cue_modifier, beverton_holt_recruitment
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
settler_geno[lethal, IDX_EF1A, 1] = 0  # Fix them

# Manually step through the settlement process
n_arriving = len(settler_geno)
n_adults_present = 20
carrying_capacity = 500

print(f"Step 1: n_arriving = {n_arriving}")

# Settlement cue Allee effect
cue_mod = settlement_cue_modifier(n_adults_present)
effective_settlers = max(0, int(n_arriving * cue_mod))
print(f"Step 2: cue_mod = {cue_mod:.3f}, effective_settlers = {effective_settlers}")

# Density-dependent recruitment (use full K; slot availability limits placement)
current_alive = int(np.sum(agents['alive']))
n_recruits = beverton_holt_recruitment(effective_settlers, carrying_capacity, 0.03)
available_slots = max(0, carrying_capacity - current_alive)
print(f"Step 3: current_alive = {current_alive}, n_recruits = {n_recruits}, available_slots = {available_slots}")

n_recruits_final = min(n_recruits, available_slots, effective_settlers, n_arriving)
print(f"Step 4: n_recruits_final = {n_recruits_final}")

print(f"Settlement should succeed: {n_recruits_final > 0}")

# Now try the actual function
from sswd_evoepi.reproduction import settle_recruits
effect_sizes_arr = initialize_effect_sizes(rng)

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

print(f"\nActual result: n_added = {n_added}")
print(f"Agents alive: {np.sum(agents['alive'])}")