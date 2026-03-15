#!/usr/bin/env python3
"""CDT Timing Analysis: Investigate disease onset timing at origin and nearby nodes.

Run a single-seed simulation with W257 params at K=5000,
track daily vibrio concentrations, infection counts, and cumulative dose
at origin nodes + their nearest neighbors.
"""

import sys
sys.path.insert(0, '/home/starbot/.openclaw/workspace/sswd-evoepi')

import numpy as np
from sswd_evoepi.config import SimulationConfig, DiseaseSection
from sswd_evoepi.disease import (
    force_of_infection, infection_probability,
    update_vibrio_concentration, vibrio_decay_rate,
    shedding_rate_I1, shedding_rate_I2, environmental_vibrio,
)

# ── W257 parameters ──
cfg = DiseaseSection()
# Overrides for W257
cfg.K_half = 1_200_000.0
cfg.cumulative_dose_threshold = 1000.0
cfg.seed_vibrio = 2000.0

# Channel Islands origin nodes
ORIGIN_NODES = [322, 319, 632, 633, 634]

# ── Quick analytical model for origin node ──
# Simulate vibrio concentration and secondary infections at a single origin node
# Key question: how quickly do secondary infections appear with seed_vibrio=2000?

print("=" * 80)
print("CDT TIMING ANALYSIS")
print("=" * 80)

# Parameters
T_sst = 16.0  # Southern CA SST ~ 16°C in summer
salinity = 33.5  # Typical ocean salinity
phi = 0.08  # Moderate flushing
N_pop = 100  # Population at node (K=5000 spread over 896 nodes... but origin nodes would have ~5-10)
# Actually K=5000 per node. Let's use that.
N_pop = 5000
n_initial_E = 5
r_mean = 0.15  # Initial mean resistance

print(f"\nNode parameters:")
print(f"  SST = {T_sst}°C, Salinity = {salinity} psu, Flushing = {phi}")
print(f"  Population = {N_pop}, Initial E = {n_initial_E}")
print(f"  K_half = {cfg.K_half:,.0f}, seed_vibrio = {cfg.seed_vibrio}")
print(f"  CDT = {cfg.cumulative_dose_threshold}")

# Shedding rates
shed_I1 = shedding_rate_I1(T_sst, cfg)
shed_I2 = shedding_rate_I2(T_sst, cfg)
print(f"\nShedding rates at {T_sst}°C:")
print(f"  I1: {shed_I1:.1f} bact/mL/day/individual")
print(f"  I2: {shed_I2:.1f} bact/mL/day/individual")

# Decay rate
xi = vibrio_decay_rate(T_sst)
print(f"  Decay: {xi:.4f} /day (half-life = {0.693/xi:.1f} days)")

# Environmental vibrio (P_env background)
P_env = environmental_vibrio(T_sst, salinity, cfg)
print(f"  Environmental P_env: {P_env:.1f} bact/mL")

print(f"\n--- Force of infection analysis ---")
# At different P_k levels, what's the daily infection probability?
for P_k in [100, 500, 1000, 2000, 5000, 10000, 50000, 100000, 500000, 1000000, 2000000]:
    lam = force_of_infection(P_k, r_mean, salinity, 150.0, cfg)
    p = infection_probability(lam)
    expected = p * (N_pop - n_initial_E)
    print(f"  P_k={P_k:>10,}: λ={lam:.6f}, p={p:.6f}, E[new infections]={expected:.1f}/day")

# Simulate first 200 days
print(f"\n--- Day-by-day simulation at origin node ---")
print(f"{'Day':>4} {'P_k':>10} {'n_S':>6} {'n_E':>5} {'n_I1':>5} {'n_I2':>5} {'n_D':>5} {'new_E':>6} {'cum_inf':>8}")

P_k = cfg.seed_vibrio  # Start with seed vibrio
n_S = N_pop - n_initial_E
n_E = n_initial_E
n_I1 = 0
n_I2 = 0
n_D = 0
n_D_fresh = 0
cum_infections = n_initial_E

rng = np.random.default_rng(137)

# Disease stage durations (rough)
# E→I1: ~3-7 days, I1→I2: ~5-10 days, I2→death: ~5-14 days
# Using Arrhenius-adjusted means at 16°C
from sswd_evoepi.disease import arrhenius, sample_stage_duration, K_SHAPE_E, K_SHAPE_I1, K_SHAPE_I2
mu_EI1 = arrhenius(cfg.mu_EI1_ref, cfg.Ea_EI1, T_sst)
mu_I1I2 = arrhenius(cfg.mu_I1I2_ref, cfg.Ea_I1I2, T_sst)
mu_I2D = arrhenius(cfg.mu_I2D_ref, cfg.Ea_I2D, T_sst)
print(f"\nStage durations at {T_sst}°C:")
print(f"  E→I1: mean={mu_EI1:.1f} days")
print(f"  I1→I2: mean={mu_I1I2:.1f} days")
print(f"  I2→D: mean={mu_I2D:.1f} days")

# Track individuals as simple lists of (state, timer)
individuals = []
for _ in range(n_initial_E):
    individuals.append(['E', sample_stage_duration(mu_EI1, K_SHAPE_E, rng)])
for _ in range(N_pop - n_initial_E):
    individuals.append(['S', 0])

first_secondary = None
pct_10 = None
pct_50 = None

for day in range(730):  # 2 years
    # Count states
    n_S = sum(1 for s, _ in individuals if s == 'S')
    n_E = sum(1 for s, _ in individuals if s == 'E')
    n_I1 = sum(1 for s, _ in individuals if s == 'I1')
    n_I2 = sum(1 for s, _ in individuals if s == 'I2')
    n_D = sum(1 for s, _ in individuals if s == 'D')
    n_D_fresh = sum(1 for s, t in individuals if s == 'D' and t < 3)

    # Update vibrio (no dispersal input for simplicity)
    P_k = update_vibrio_concentration(
        P_k, n_I1, n_I2, n_D_fresh, T_sst, salinity, phi,
        0.0, cfg, disease_reached=True,
        P_env_pool=0.0 if not cfg.P_env_dynamic else 0.0,
    )

    # New infections
    new_E = 0
    for idx in range(len(individuals)):
        state, timer = individuals[idx]
        if state == 'S':
            r_i = max(0, rng.normal(r_mean, 0.05))
            lam = force_of_infection(P_k, r_i, salinity, 150.0, cfg)
            p = infection_probability(lam)
            if rng.random() < p:
                individuals[idx] = ['E', sample_stage_duration(mu_EI1, K_SHAPE_E, rng)]
                new_E += 1
                cum_infections += 1

    # State transitions
    for idx in range(len(individuals)):
        state, timer = individuals[idx]
        if state == 'E':
            timer -= 1
            if timer <= 0:
                individuals[idx] = ['I1', sample_stage_duration(mu_I1I2, K_SHAPE_I1, rng)]
            else:
                individuals[idx] = ['E', timer]
        elif state == 'I1':
            timer -= 1
            if timer <= 0:
                individuals[idx] = ['I2', sample_stage_duration(mu_I2D, K_SHAPE_I2, rng)]
            else:
                individuals[idx] = ['I1', timer]
        elif state == 'I2':
            timer -= 1
            if timer <= 0:
                individuals[idx] = ['D', 0]
            else:
                individuals[idx] = ['I2', timer]
        elif state == 'D':
            individuals[idx] = ['D', timer + 1]  # Increment days dead

    if first_secondary is None and new_E > 0 and day > 0:
        first_secondary = day
    if pct_10 is None and cum_infections >= N_pop * 0.10:
        pct_10 = day
    if pct_50 is None and cum_infections >= N_pop * 0.50:
        pct_50 = day

    if day < 60 or day % 30 == 0 or new_E > 0 and day < 180:
        print(f"{day:4d} {P_k:10.1f} {n_S:6d} {n_E:5d} {n_I1:5d} {n_I2:5d} {n_D:5d} {new_E:6d} {cum_infections:8d}")

print(f"\n--- Timing Summary ---")
print(f"  First secondary infection: day {first_secondary}")
print(f"  10% infected: day {pct_10}")
print(f"  50% infected: day {pct_50}")
print(f"  Final cumulative infections: {cum_infections}")

# Now analyze CDT timing for non-origin nodes
print(f"\n{'='*80}")
print("CDT ANALYSIS FOR NON-ORIGIN (WAVEFRONT) NODES")
print(f"{'='*80}")

# At what rate does dispersal arrive at neighboring nodes?
# Dispersal = D_ij * P_j where D_ij is the dispersal kernel coefficient
# For nearby nodes (~50-100 km), D_ij is small

# Let's estimate: peak P at origin ~ 5000-20000 bact/mL
# D_ij for a node 50km away: using exponential kernel D_P=50km
from sswd_evoepi.config import SpatialSection
sp = SpatialSection()
print(f"\nDispersal kernel: D_P={sp.D_P} km, D_P_max_range={sp.D_P_max_range} km")

# d_ij = r_total * exp(-dist/D_P) / sum normalization
# For simplicity, estimate raw kernel value
for dist in [10, 25, 50, 100, 150, 200]:
    raw_kernel = np.exp(-dist / sp.D_P)
    print(f"  dist={dist:4d}km: raw_kernel={raw_kernel:.6f}")

# With CDT=1000 and different daily dose inputs, how many days to activate?
print(f"\nDays to CDT activation (CDT={cfg.cumulative_dose_threshold}, no decay):")
for daily_dose in [0.1, 0.5, 1, 5, 10, 20, 50, 100, 200]:
    days_needed = cfg.cumulative_dose_threshold / daily_dose
    months = days_needed / 30
    print(f"  daily_dose={daily_dose:7.1f} bact/mL/day → {days_needed:8.1f} days ({months:5.1f} months)")

# With decay
print(f"\nWith dose_decay_rate=0.02 (2%/day):")
for daily_dose in [0.1, 0.5, 1, 5, 10, 20, 50, 100, 200]:
    # Steady state: dose = daily_dose / decay_rate
    ss = daily_dose / 0.02
    if ss > cfg.cumulative_dose_threshold:
        # Find exact day by simulation
        cum = 0
        for d in range(10000):
            cum += daily_dose
            cum *= (1.0 - 0.02)
            if cum > cfg.cumulative_dose_threshold:
                print(f"  daily_dose={daily_dose:7.1f}: activates day {d+1} ({(d+1)/30:.1f} months), ss={ss:.0f}")
                break
        else:
            print(f"  daily_dose={daily_dose:7.1f}: NEVER (steady state {ss:.0f} < CDT)")
    else:
        print(f"  daily_dose={daily_dose:7.1f}: NEVER (steady state {ss:.0f} < CDT)")

# Sweep CDT values
print(f"\n{'='*80}")
print("CDT SWEEP: Days to activation at different CDT values")
print(f"{'='*80}")
print(f"\nAssuming daily dispersal dose = 20 bact/mL/day (nearby node), no decay:")
for cdt in [100, 250, 500, 750, 1000, 1500, 2000]:
    days = cdt / 20
    print(f"  CDT={cdt:5d}: {days:6.1f} days ({days/30:.1f} months)")

print(f"\nAssuming daily dispersal dose = 5 bact/mL/day (medium-distance node), no decay:")
for cdt in [100, 250, 500, 750, 1000, 1500, 2000]:
    days = cdt / 5
    print(f"  CDT={cdt:5d}: {days:6.1f} days ({days/30:.1f} months)")
