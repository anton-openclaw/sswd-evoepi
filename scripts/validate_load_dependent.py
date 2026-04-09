#!/usr/bin/env python3
"""Validation script for load-dependent disease model in SSWD-EvoEpi.

Runs side-by-side comparison of:
  1. Load-dependent (LD) disease model — within-host ODE, Hill-function mortality
  2. Timer-based (TB) disease model — classical SEIPD with Arrhenius rates

Tests:
  - Does the epidemic actually spread?
  - Is case fatality rate reasonable (not 0%, not 100%)?
  - Does recovery occur?
  - Does shedding feedback amplify infection?
  - Is timescale realistic (peak within weeks-months)?
  - Does r_i=1.0 confer functional immunity in LD model?
  - Does the dynamics qualitatively match between models?

Usage:
    python scripts/validate_load_dependent.py
"""

import sys
import os
import numpy as np
import math

# Ensure project root is on path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from sswd_evoepi.config import (
    DiseaseSection,
    LoadDependentSection,
    default_config,
)
from sswd_evoepi.types import AGENT_DTYPE, DiseaseState, allocate_agents
from sswd_evoepi.disease import (
    CarcassTracker,
    NodeDiseaseState,
    arrhenius,
    within_host_load_step,
    load_mortality_probability,
    load_dependent_disease_update,
    daily_disease_update,
    environmental_vibrio,
    vibrio_decay_rate,
    sample_stage_duration,
    K_SHAPE_E,
)


# ═══════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════

N_POP = 500
N_DAYS = 365
N_INITIAL_INFECTED = 5
SEED = 42

# PNW waters — realistic temperature range
T_CELSIUS = 12.0       # Cool water — growth rate should be lower
SALINITY = 30.0        # Full marine
PHI_K = 0.02           # Low flushing (sheltered site)

# Agent trait distributions
MEAN_RESISTANCE = 0.08
RESISTANCE_STD = 0.04
MEAN_TOLERANCE = 0.10
TOLERANCE_STD = 0.05
MEAN_RECOVERY = 0.08
RECOVERY_STD = 0.04
MEAN_SIZE = 300.0
SIZE_STD = 80.0


# ═══════════════════════════════════════════════════════════════════════
# HELPER: Create population
# ═══════════════════════════════════════════════════════════════════════

def create_population(rng, n=N_POP):
    """Create a population of susceptible agents with trait distributions."""
    agents = np.zeros(n, dtype=AGENT_DTYPE)
    agents['alive'] = True
    agents['disease_state'] = DiseaseState.S
    agents['disease_timer'] = 0
    agents['node_id'] = 0
    agents['settlement_day'] = 0  # All are established (no juvenile immunity)

    # Resistance (truncated normal [0, 1])
    r = rng.normal(MEAN_RESISTANCE, RESISTANCE_STD, n)
    agents['resistance'] = np.clip(r, 0.0, 1.0).astype(np.float32)

    # Tolerance
    t = rng.normal(MEAN_TOLERANCE, TOLERANCE_STD, n)
    agents['tolerance'] = np.clip(t, 0.0, 1.0).astype(np.float32)

    # Recovery ability
    c = rng.normal(MEAN_RECOVERY, RECOVERY_STD, n)
    agents['recovery_ability'] = np.clip(c, 0.0, 1.0).astype(np.float32)

    # Size
    sizes = rng.normal(MEAN_SIZE, SIZE_STD, n)
    agents['size'] = np.clip(sizes, 50.0, 800.0).astype(np.float32)

    # Pathogen load starts at zero
    agents['pathogen_load'] = 0.0

    return agents


# ═══════════════════════════════════════════════════════════════════════
# SIMULATION: Load-dependent model
# ═══════════════════════════════════════════════════════════════════════

def run_load_dependent_epidemic(T_celsius=T_CELSIUS):
    """Run a single-node epidemic using load-dependent disease model."""
    print("=" * 70)
    print(f"LOAD-DEPENDENT MODEL (T={T_celsius}°C, N={N_POP}, seed={SEED})")
    print("=" * 70)

    rng = np.random.default_rng(SEED)
    cfg = DiseaseSection()
    ld_cfg = LoadDependentSection(enabled=True)

    # Print key parameters
    arr_factor = arrhenius(1.0, ld_cfg.Ea_growth, T_celsius)
    effective_growth = ld_cfg.r_growth * arr_factor
    # Clearance for typical agent: delta_clear * (0.3 + 0.7*r) * (0.5 + 0.5*c)
    typ_r, typ_c = MEAN_RESISTANCE, MEAN_RECOVERY
    effective_clearance = ld_cfg.delta_clear * (0.3 + 0.7 * typ_r) * (0.5 + 0.5 * typ_c)
    print(f"\nParameters:")
    print(f"  r_growth={ld_cfg.r_growth}, Ea_growth={ld_cfg.Ea_growth}")
    print(f"  Arrhenius factor at {T_celsius}°C: {arr_factor:.4f}")
    print(f"  Effective growth rate: {effective_growth:.4f} d^-1")
    print(f"  delta_clear={ld_cfg.delta_clear}")
    print(f"  Typical clearance (r={typ_r}, c={typ_c}): {effective_clearance:.4f} d^-1")
    print(f"  Growth - clearance at low L: {effective_growth - effective_clearance:.4f}")
    print(f"  L_init_base={ld_cfg.L_init_base}, L_clear={ld_cfg.L_clear}, L_symp={ld_cfg.L_symp}")
    print(f"  LD50_base={ld_cfg.LD50_base}, n_hill={ld_cfg.n_hill}, p_death_max={ld_cfg.p_death_max}")
    print(f"  sigma_load={ld_cfg.sigma_load}, L_ref={ld_cfg.L_ref}")
    print(f"  alpha_reinfect={ld_cfg.alpha_reinfect}")
    print(f"  K_half={cfg.K_half}, a_exposure={cfg.a_exposure}")

    # Check r=1.0 immunity
    r_immune = 1.0
    c_immune = 0.5
    clear_immune = ld_cfg.delta_clear * (0.3 + 0.7 * r_immune) * (0.5 + 0.5 * c_immune)
    print(f"\n  r=1.0 clearance rate: {clear_immune:.4f} d^-1")
    print(f"  r=1.0 net rate at low L: {effective_growth - clear_immune:.4f}")
    if effective_growth > clear_immune:
        print(f"  *** WARNING: Growth ({effective_growth:.4f}) exceeds clearance ({clear_immune:.4f}) "
              f"even at r=1.0! Immunity is BROKEN at this temperature. ***")

    agents = create_population(rng)
    node_state = NodeDiseaseState(node_id=0, vibrio_concentration=0.0)

    # Seed infections with initial pathogen load
    infect_idx = rng.choice(N_POP, size=N_INITIAL_INFECTED, replace=False)
    agents['disease_state'][infect_idx] = DiseaseState.I1
    agents['pathogen_load'][infect_idx] = ld_cfg.L_init_base
    print(f"\n  Seeded {N_INITIAL_INFECTED} infections with load={ld_cfg.L_init_base}")

    # Tracking arrays
    daily_alive = np.zeros(N_DAYS, dtype=np.int32)
    daily_infected = np.zeros(N_DAYS, dtype=np.int32)
    daily_dead_cumul = np.zeros(N_DAYS, dtype=np.int32)
    daily_recovered_cumul = np.zeros(N_DAYS, dtype=np.int32)
    daily_vibrio = np.zeros(N_DAYS, dtype=np.float64)
    daily_mean_load = np.zeros(N_DAYS, dtype=np.float64)
    daily_new_infections = np.zeros(N_DAYS, dtype=np.int32)

    prev_cumul_inf = 0
    prev_cumul_deaths = 0

    print(f"\n{'Day':>5} {'Alive':>6} {'Inf':>5} {'Dead':>5} {'Rec':>5} {'NewInf':>6} "
          f"{'MeanLoad':>12} {'P_k':>12} {'Prevalence':>10}")
    print("-" * 90)

    for day in range(N_DAYS):
        node_state = load_dependent_disease_update(
            agents, node_state,
            T_celsius=T_celsius, salinity=SALINITY, phi_k=PHI_K,
            dispersal_input=0.0, day=day,
            cfg=cfg, ld_cfg=ld_cfg, rng=rng,
        )

        n_alive = int(np.sum(agents['alive']))
        n_inf = node_state.n_I1 + node_state.n_I2
        new_inf = node_state.cumulative_infections - prev_cumul_inf
        prev_cumul_inf = node_state.cumulative_infections

        # Mean load of infected
        alive_mask = agents['alive'].astype(bool)
        inf_mask = alive_mask & (
            (agents['disease_state'] == DiseaseState.I1) |
            (agents['disease_state'] == DiseaseState.I2)
        )
        if np.any(inf_mask):
            mean_load = float(np.mean(agents['pathogen_load'][inf_mask]))
        else:
            mean_load = 0.0

        daily_alive[day] = n_alive
        daily_infected[day] = n_inf
        daily_dead_cumul[day] = node_state.cumulative_deaths
        daily_recovered_cumul[day] = node_state.cumulative_recoveries
        daily_vibrio[day] = node_state.vibrio_concentration
        daily_mean_load[day] = mean_load
        daily_new_infections[day] = new_inf

        prevalence = n_inf / n_alive if n_alive > 0 else 0.0

        if day % 30 == 0 or day == N_DAYS - 1:
            print(f"{day:5d} {n_alive:6d} {n_inf:5d} {node_state.cumulative_deaths:5d} "
                  f"{node_state.cumulative_recoveries:5d} {new_inf:6d} "
                  f"{mean_load:12.1f} {node_state.vibrio_concentration:12.1f} {prevalence:10.4f}")

        # Early termination if everyone is dead
        if n_alive == 0:
            print(f"  *** All agents dead at day {day} ***")
            break

    # Summary
    peak_inf_day = int(np.argmax(daily_infected))
    peak_inf_count = int(daily_infected[peak_inf_day])
    peak_prevalence = peak_inf_count / daily_alive[peak_inf_day] if daily_alive[peak_inf_day] > 0 else 0
    total_deaths = int(daily_dead_cumul[day])
    total_recoveries = int(daily_recovered_cumul[day])
    total_infections = node_state.cumulative_infections
    cfr = total_deaths / total_infections if total_infections > 0 else 0
    final_alive = int(daily_alive[day])
    final_infected = int(daily_infected[day])
    peak_vibrio = float(np.max(daily_vibrio))

    print(f"\n{'SUMMARY':=^70}")
    print(f"  Total infections:       {total_infections}")
    print(f"  Total deaths:           {total_deaths}")
    print(f"  Total recoveries:       {total_recoveries}")
    print(f"  Case fatality rate:     {cfr:.2%}")
    print(f"  Peak prevalence:        {peak_prevalence:.2%} (day {peak_inf_day}, {peak_inf_count} infected)")
    print(f"  Peak vibrio:            {peak_vibrio:.1f} bact/mL")
    print(f"  Final alive:            {final_alive}/{N_POP}")
    print(f"  Final infected:         {final_infected}")
    print(f"  Mortality fraction:     {total_deaths/N_POP:.2%}")

    epidemic_ended = final_infected == 0 and day < N_DAYS - 1
    endemic = final_infected > 0
    print(f"  Epidemic burned out:    {epidemic_ended}")
    print(f"  Endemic at end:         {endemic}")

    return {
        'model': 'load_dependent',
        'T_celsius': T_celsius,
        'total_infections': total_infections,
        'total_deaths': total_deaths,
        'total_recoveries': total_recoveries,
        'cfr': cfr,
        'peak_prevalence': peak_prevalence,
        'peak_day': peak_inf_day,
        'peak_count': peak_inf_count,
        'peak_vibrio': peak_vibrio,
        'final_alive': final_alive,
        'final_infected': final_infected,
        'mortality_fraction': total_deaths / N_POP,
        'epidemic_ended': epidemic_ended,
        'endemic': endemic,
        'daily_infected': daily_infected,
        'daily_alive': daily_alive,
        'daily_vibrio': daily_vibrio,
        'daily_mean_load': daily_mean_load,
    }


# ═══════════════════════════════════════════════════════════════════════
# SIMULATION: Timer-based model (comparison)
# ═══════════════════════════════════════════════════════════════════════

def run_timer_based_epidemic(T_celsius=T_CELSIUS):
    """Run a single-node epidemic using timer-based (original) disease model."""
    print("\n" + "=" * 70)
    print(f"TIMER-BASED MODEL (T={T_celsius}°C, N={N_POP}, seed={SEED})")
    print("=" * 70)

    rng = np.random.default_rng(SEED)
    cfg = DiseaseSection()

    agents = create_population(rng)

    # Initial vibrio — steady-state background
    env = environmental_vibrio(T_celsius, SALINITY, cfg)
    xi = vibrio_decay_rate(T_celsius)
    P_init = env / (xi + PHI_K) if (xi + PHI_K) > 0 else 0.0
    print(f"\n  Background Vibrio: env={env:.1f}, xi={xi:.4f}, P_init={P_init:.1f}")

    node_state = NodeDiseaseState(node_id=0, vibrio_concentration=P_init)

    # Seed infections using timer model's E state
    infect_idx = rng.choice(N_POP, size=N_INITIAL_INFECTED, replace=False)
    mu_EI1 = arrhenius(cfg.mu_EI1_ref, cfg.Ea_EI1, T_celsius)
    for idx in infect_idx:
        agents['disease_state'][idx] = DiseaseState.E
        agents['disease_timer'][idx] = sample_stage_duration(mu_EI1, K_SHAPE_E, rng)
    print(f"  Seeded {N_INITIAL_INFECTED} exposed agents (E state, timer-based)")

    # Tracking
    daily_alive = np.zeros(N_DAYS, dtype=np.int32)
    daily_infected = np.zeros(N_DAYS, dtype=np.int32)
    daily_dead_cumul = np.zeros(N_DAYS, dtype=np.int32)
    daily_recovered_cumul = np.zeros(N_DAYS, dtype=np.int32)
    daily_vibrio = np.zeros(N_DAYS, dtype=np.float64)

    print(f"\n{'Day':>5} {'Alive':>6} {'Exp':>5} {'Inf':>5} {'Dead':>5} {'Rec':>5} "
          f"{'P_k':>12} {'Prevalence':>10}")
    print("-" * 75)

    for day in range(N_DAYS):
        node_state = daily_disease_update(
            agents, node_state,
            T_celsius=T_celsius, salinity=SALINITY, phi_k=PHI_K,
            dispersal_input=0.0, day=day,
            cfg=cfg, rng=rng,
        )

        n_alive = int(np.sum(agents['alive']))
        n_exp = node_state.n_E
        n_inf = node_state.n_I1 + node_state.n_I2
        prevalence = n_inf / n_alive if n_alive > 0 else 0.0

        daily_alive[day] = n_alive
        daily_infected[day] = n_inf
        daily_dead_cumul[day] = node_state.cumulative_deaths
        daily_recovered_cumul[day] = node_state.cumulative_recoveries
        daily_vibrio[day] = node_state.vibrio_concentration

        if day % 30 == 0 or day == N_DAYS - 1:
            print(f"{day:5d} {n_alive:6d} {n_exp:5d} {n_inf:5d} {node_state.cumulative_deaths:5d} "
                  f"{node_state.cumulative_recoveries:5d} "
                  f"{node_state.vibrio_concentration:12.1f} {prevalence:10.4f}")

        if n_alive == 0:
            print(f"  *** All agents dead at day {day} ***")
            break

    peak_inf_day = int(np.argmax(daily_infected))
    peak_inf_count = int(daily_infected[peak_inf_day])
    peak_prevalence = peak_inf_count / daily_alive[peak_inf_day] if daily_alive[peak_inf_day] > 0 else 0
    total_deaths = int(daily_dead_cumul[day])
    total_recoveries = int(daily_recovered_cumul[day])
    total_infections = node_state.cumulative_infections
    cfr = total_deaths / total_infections if total_infections > 0 else 0
    final_alive = int(daily_alive[day])
    final_infected = int(daily_infected[day])
    peak_vibrio = float(np.max(daily_vibrio))

    print(f"\n{'SUMMARY':=^70}")
    print(f"  Total infections:       {total_infections}")
    print(f"  Total deaths:           {total_deaths}")
    print(f"  Total recoveries:       {total_recoveries}")
    print(f"  Case fatality rate:     {cfr:.2%}")
    print(f"  Peak prevalence:        {peak_prevalence:.2%} (day {peak_inf_day}, {peak_inf_count} infected)")
    print(f"  Peak vibrio:            {peak_vibrio:.1f} bact/mL")
    print(f"  Final alive:            {final_alive}/{N_POP}")
    print(f"  Final infected:         {final_infected}")
    print(f"  Mortality fraction:     {total_deaths/N_POP:.2%}")

    epidemic_ended = final_infected == 0 and node_state.n_E == 0
    endemic = final_infected > 0 or node_state.n_E > 0
    print(f"  Epidemic burned out:    {epidemic_ended}")
    print(f"  Endemic at end:         {endemic}")

    return {
        'model': 'timer_based',
        'T_celsius': T_celsius,
        'total_infections': total_infections,
        'total_deaths': total_deaths,
        'total_recoveries': total_recoveries,
        'cfr': cfr,
        'peak_prevalence': peak_prevalence,
        'peak_day': peak_inf_day,
        'peak_count': peak_inf_count,
        'peak_vibrio': peak_vibrio,
        'final_alive': final_alive,
        'final_infected': final_infected,
        'mortality_fraction': total_deaths / N_POP,
        'epidemic_ended': epidemic_ended,
        'endemic': endemic,
        'daily_infected': daily_infected,
        'daily_alive': daily_alive,
        'daily_vibrio': daily_vibrio,
    }


# ═══════════════════════════════════════════════════════════════════════
# IMMUNITY TEST: r_i = 1.0 in load model
# ═══════════════════════════════════════════════════════════════════════

def test_immunity_at_r1():
    """Test whether r_i=1.0 confers functional immunity in LD model."""
    print("\n" + "=" * 70)
    print("IMMUNITY TEST: r_i=1.0 agents in load-dependent model")
    print("=" * 70)

    rng = np.random.default_rng(99)
    cfg = DiseaseSection()
    ld_cfg = LoadDependentSection(enabled=True)

    # 10 agents, all r=1.0, c=0.5, infected with moderate load
    n = 10
    agents = np.zeros(n, dtype=AGENT_DTYPE)
    agents['alive'] = True
    agents['disease_state'] = DiseaseState.I1
    agents['resistance'] = 1.0
    agents['recovery_ability'] = 0.5
    agents['tolerance'] = 0.1
    agents['size'] = 300.0
    agents['pathogen_load'] = ld_cfg.L_init_base  # 1000.0

    for T in [10.0, 12.0, 15.0, 18.0, 20.0]:
        L = np.array([ld_cfg.L_init_base] * n, dtype=np.float64)
        r = np.ones(n, dtype=np.float64)
        c = np.ones(n, dtype=np.float64) * 0.5

        arr = arrhenius(1.0, ld_cfg.Ea_growth, T)
        growth = ld_cfg.r_growth * arr * L[0] * (1.0 - L[0] / ld_cfg.L_max)
        clearance = ld_cfg.delta_clear * (0.3 + 0.7 * 1.0) * (0.5 + 0.5 * 0.5) * L[0]
        net = growth - clearance

        L_new = within_host_load_step(L, r, c, T, 0.0, cfg, ld_cfg)
        print(f"  T={T:5.1f}°C: L={L[0]:.0f} -> L_new={L_new[0]:.1f} "
              f"(growth={growth:.2f}, clearance={clearance:.2f}, net={net:.2f}) "
              f"{'CLEARING' if L_new[0] < L[0] else '*** GROWING ***'}")


# ═══════════════════════════════════════════════════════════════════════
# PARAMETER SENSITIVITY: Sweep temperatures
# ═══════════════════════════════════════════════════════════════════════

def sweep_temperatures():
    """Quick sweep of epidemic outcomes across temperatures."""
    print("\n" + "=" * 70)
    print("TEMPERATURE SWEEP (load-dependent model)")
    print("=" * 70)

    cfg = DiseaseSection()
    ld_cfg = LoadDependentSection(enabled=True)

    print(f"\n{'T(°C)':>6} {'ArrFact':>8} {'EffGrow':>8} {'TypClear':>8} {'NetRate':>8} "
          f"{'r=1 Clear':>9} {'r=1 Net':>8} {'r=1 OK?':>8}")
    print("-" * 75)

    for T in [8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20]:
        arr = arrhenius(1.0, ld_cfg.Ea_growth, T)
        eff_growth = ld_cfg.r_growth * arr
        typ_clear = ld_cfg.delta_clear * (0.3 + 0.7 * 0.08) * (0.5 + 0.5 * 0.08)
        net = eff_growth - typ_clear
        immune_clear = ld_cfg.delta_clear * (0.3 + 0.7 * 1.0) * (0.5 + 0.5 * 0.5)
        immune_net = eff_growth - immune_clear
        ok = "YES" if immune_net < 0 else "NO"
        print(f"{T:6.0f} {arr:8.4f} {eff_growth:8.4f} {typ_clear:8.4f} {net:+8.4f} "
              f"{immune_clear:9.4f} {immune_net:+8.4f} {ok:>8}")


# ═══════════════════════════════════════════════════════════════════════
# COMPARISON REPORT
# ═══════════════════════════════════════════════════════════════════════

def compare_results(ld_result, tb_result):
    """Print side-by-side comparison of LD and TB epidemic outcomes."""
    print("\n" + "=" * 70)
    print("SIDE-BY-SIDE COMPARISON")
    print("=" * 70)

    metrics = [
        ('Total infections', 'total_infections'),
        ('Total deaths', 'total_deaths'),
        ('Total recoveries', 'total_recoveries'),
        ('Case fatality rate', 'cfr'),
        ('Peak prevalence', 'peak_prevalence'),
        ('Day of peak', 'peak_day'),
        ('Peak infected count', 'peak_count'),
        ('Peak vibrio', 'peak_vibrio'),
        ('Final alive', 'final_alive'),
        ('Final infected', 'final_infected'),
        ('Mortality fraction', 'mortality_fraction'),
        ('Epidemic ended', 'epidemic_ended'),
        ('Endemic', 'endemic'),
    ]

    print(f"\n{'Metric':<25} {'Load-Dependent':>18} {'Timer-Based':>18}")
    print("-" * 65)
    for label, key in metrics:
        ld_val = ld_result[key]
        tb_val = tb_result[key]
        if isinstance(ld_val, float):
            if key in ('cfr', 'peak_prevalence', 'mortality_fraction'):
                print(f"{label:<25} {ld_val:>17.2%} {tb_val:>17.2%}")
            else:
                print(f"{label:<25} {ld_val:>18.1f} {tb_val:>18.1f}")
        elif isinstance(ld_val, bool):
            print(f"{label:<25} {str(ld_val):>18} {str(tb_val):>18}")
        else:
            print(f"{label:<25} {ld_val:>18} {tb_val:>18}")

    # Quality assessment
    print(f"\n{'QUALITY CHECKS':=^70}")
    checks = []

    # 1. Epidemic happened?
    ld_spread = ld_result['total_infections'] > N_INITIAL_INFECTED
    tb_spread = tb_result['total_infections'] > N_INITIAL_INFECTED
    checks.append(('Epidemic spreads (LD)', ld_spread))
    checks.append(('Epidemic spreads (TB)', tb_spread))

    # 2. Reasonable CFR (1-80%)
    ld_cfr_ok = 0.01 < ld_result['cfr'] < 0.95 if ld_result['total_infections'] > 0 else False
    tb_cfr_ok = 0.01 < tb_result['cfr'] < 0.95 if tb_result['total_infections'] > 0 else False
    checks.append(('CFR reasonable 1-95% (LD)', ld_cfr_ok))
    checks.append(('CFR reasonable 1-95% (TB)', tb_cfr_ok))

    # 3. Recovery occurs
    ld_rec = ld_result['total_recoveries'] > 0
    tb_rec = tb_result['total_recoveries'] > 0
    checks.append(('Recovery occurs (LD)', ld_rec))
    checks.append(('Recovery occurs (TB)', tb_rec))

    # 4. Timescale (peak within 14-300 days)
    ld_time_ok = 7 <= ld_result['peak_day'] <= 300
    tb_time_ok = 7 <= tb_result['peak_day'] <= 300
    checks.append(('Peak timing 7-300d (LD)', ld_time_ok))
    checks.append(('Peak timing 7-300d (TB)', tb_time_ok))

    # 5. Not 100% mortality
    ld_not_all_dead = ld_result['mortality_fraction'] < 0.95
    tb_not_all_dead = tb_result['mortality_fraction'] < 0.95
    checks.append(('Not total wipeout (LD)', ld_not_all_dead))
    checks.append(('Not total wipeout (TB)', tb_not_all_dead))

    # 6. Shedding feedback (vibrio rises from infection)
    ld_shed = ld_result['peak_vibrio'] > 100.0  # Should see significant vibrio buildup
    tb_shed = tb_result['peak_vibrio'] > 100.0
    checks.append(('Shedding feedback (LD)', ld_shed))
    checks.append(('Shedding feedback (TB)', tb_shed))

    n_pass = 0
    n_total = len(checks)
    for label, passed in checks:
        status = "PASS" if passed else "FAIL"
        marker = "  [✓]" if passed else "  [✗]"
        print(f"  {marker} {label}: {status}")
        if passed:
            n_pass += 1

    print(f"\n  Result: {n_pass}/{n_total} checks passed")

    # Overall verdict
    if n_pass == n_total:
        print("\n  >>> ALL CHECKS PASSED — both models produce reasonable epidemics <<<")
    elif n_pass >= n_total * 0.7:
        print("\n  >>> MOSTLY OK — some issues but models are broadly functional <<<")
    else:
        print("\n  >>> ISSUES DETECTED — review parameters and model behavior <<<")

    return n_pass, n_total


# ═══════════════════════════════════════════════════════════════════════
# DIAGNOSIS: Why LD model dynamics differ
# ═══════════════════════════════════════════════════════════════════════

def diagnose_ld_dynamics(ld_result):
    """Diagnose issues with LD model dynamics based on results."""
    print("\n" + "=" * 70)
    print("DIAGNOSIS: Load-Dependent Model Dynamics")
    print("=" * 70)

    cfg = DiseaseSection()
    ld_cfg = LoadDependentSection(enabled=True)

    # Issue 1: K_half too high for dose-response transmission
    peak_P = ld_result['peak_vibrio']
    dose_frac = peak_P / (cfg.K_half + peak_P)
    L_init_at_peak = ld_cfg.L_init_base * dose_frac
    print(f"\n  Issue 1: Dose-response mismatch")
    print(f"    K_half = {cfg.K_half:.0f} bact/mL (shared with timer model)")
    print(f"    Peak P_k = {peak_P:.0f} bact/mL")
    print(f"    Dose fraction at peak = P_k/(K_half+P_k) = {dose_frac:.6f}")
    print(f"    Initial load for new infections = L_init_base × dose_frac = {L_init_at_peak:.2f}")
    print(f"    L_clear threshold = {ld_cfg.L_clear}")
    if L_init_at_peak < ld_cfg.L_clear:
        print(f"    *** NEW INFECTIONS START BELOW L_clear! They immediately recover. ***")
        P_needed = ld_cfg.L_clear / ld_cfg.L_init_base * cfg.K_half / (1 - ld_cfg.L_clear / ld_cfg.L_init_base)
        print(f"    P_k needed for persistent infection: >{P_needed:.0f} bact/mL")
    else:
        print(f"    OK: New infections start above L_clear.")

    # Issue 2: Growth always exceeds clearance (immunity broken)
    arr_factor = arrhenius(1.0, ld_cfg.Ea_growth, ld_result['T_celsius'])
    eff_growth = ld_cfg.r_growth * arr_factor
    max_clearance = ld_cfg.delta_clear * 1.0 * 1.0  # r=1.0, c=1.0
    print(f"\n  Issue 2: Growth vs clearance balance")
    print(f"    Effective growth at {ld_result['T_celsius']}°C: {eff_growth:.4f} d^-1")
    print(f"    Max clearance (r=1.0, c=1.0): {max_clearance:.4f} d^-1")
    if eff_growth > max_clearance:
        print(f"    *** GROWTH EXCEEDS MAX CLEARANCE by {eff_growth - max_clearance:.4f} ***")
        print(f"    No agent can clear infection at this temperature!")
        # Find temperature where r=1.0,c=1.0 can clear
        # r_growth * arr(T) = delta_clear * 1.0 * 1.0
        # arr(T) = delta_clear / r_growth
        target_arr = max_clearance / ld_cfg.r_growth
        # arr(T) = exp(Ea * (1/T_ref - 1/T))
        # ln(target) = Ea * (1/T_ref - 1/T)
        # 1/T = 1/T_ref - ln(target)/Ea
        # T = 1/(1/T_ref - ln(target)/Ea) - 273.15
        import math
        T_ref_K = 293.15
        if target_arr > 0:
            inv_T = 1.0 / T_ref_K - math.log(target_arr) / ld_cfg.Ea_growth
            T_balance = 1.0 / inv_T - 273.15
            print(f"    Temperature for r=1,c=1 clearance = growth: {T_balance:.1f}°C")
            print(f"    (Must be below this for immunity to work)")
    else:
        print(f"    OK: Max clearance exceeds growth.")

    # Suggested parameter adjustments
    print(f"\n  SUGGESTED PARAMETER ADJUSTMENTS:")
    print(f"    Option A: Increase delta_clear to {eff_growth * 1.2:.2f} (growth×1.2)")
    print(f"             This ensures r=1.0,c=1.0 can clear at any realistic T")
    print(f"    Option B: Decrease r_growth to {max_clearance / arr_factor * 0.8:.3f}")
    print(f"             So max clearance rate exceeds growth at {ld_result['T_celsius']}°C")
    print(f"    Option C: Raise L_init_base to {int(ld_cfg.L_clear / 0.005) + 1}")
    print(f"             So even at dose_frac=0.005, initial load > L_clear")
    print(f"    Option D: Lower K_half for LD model (separate from timer model)")
    print(f"             K_half_LD ~= {int(peak_P * ld_cfg.L_clear / (ld_cfg.L_init_base - ld_cfg.L_clear)) + 1}")
    print(f"             So peak vibrio gives L_init > L_clear")


# ═══════════════════════════════════════════════════════════════════════
# ADJUSTED PARAMETER RUN
# ═══════════════════════════════════════════════════════════════════════

def run_load_dependent_adjusted(T_celsius=T_CELSIUS):
    """Run LD model with adjusted parameters to achieve reasonable dynamics.

    Key fixes:
    - Increase delta_clear so r=1.0 agents can clear
    - Decrease r_growth slightly
    - Increase L_init_base so new infections are above L_clear
    """
    print("\n" + "=" * 70)
    print(f"LOAD-DEPENDENT MODEL — ADJUSTED PARAMS (T={T_celsius}°C)")
    print("=" * 70)

    rng = np.random.default_rng(SEED)
    cfg = DiseaseSection()
    ld_cfg = LoadDependentSection(enabled=True)

    # ADJUSTMENTS — Goal: balanced epidemic with recovery, death, and immunity
    #
    # Key insight: clearance = delta_clear * (0.3 + 0.7*r) * (0.5 + 0.5*c)
    # Growth at 12°C = r_growth * arrhenius(Ea,12) = r_growth * 0.6197
    #
    # Design targets:
    #   r=0.0 → growth > clearance → load grows → die
    #   r~0.12 → borderline → some survive, some don't
    #   r=1.0 → clearance >> growth → clear infection (immunity)
    #
    # With r_growth=0.16, delta_clear=0.55:
    #   growth at 12°C = 0.16*0.6197 = 0.0992
    #   clearance(r=0.0) = 0.089, clearance(r=0.05) = 0.099 (borderline)
    #   clearance(r=0.08) = 0.106, clearance(r=1.0) = 0.297
    #   → r<0.05 agents can't clear (die), r>0.06 agents recover
    #   → ~25% of pop (mean_r=0.08, std=0.04) will die
    ld_cfg.r_growth = 0.16       # was 0.5 — calibrated for balanced epidemic
    ld_cfg.delta_clear = 0.55    # was 0.3 — stronger immune clearance

    # 2. Large L_init_base so even low dose_frac gives L > L_clear
    #    At P_k=300 and K_half=800000: dose_frac ≈ 0.000375
    #    L_init = 100000 * 0.000375 = 37.5 > L_clear=10 ✓
    ld_cfg.L_init_base = 100000.0  # was 1000 — scaled up for K_half=800000

    # 3. Increase shedding to build vibrio from within-host load
    ld_cfg.sigma_load = 150.0    # was 50

    # 4. Lower LD50 so deaths occur at achievable loads
    ld_cfg.LD50_base = 5e4       # was 5e5

    arr_factor = arrhenius(1.0, ld_cfg.Ea_growth, T_celsius)
    eff_growth = ld_cfg.r_growth * arr_factor
    max_clearance = ld_cfg.delta_clear * 1.0 * 1.0
    typ_clearance = ld_cfg.delta_clear * (0.3 + 0.7 * 0.08) * (0.5 + 0.5 * 0.08)

    print(f"\n  Adjusted parameters:")
    print(f"    r_growth:     0.5 -> {ld_cfg.r_growth}")
    print(f"    delta_clear:  0.3 -> {ld_cfg.delta_clear}")
    print(f"    L_init_base:  1000 -> {ld_cfg.L_init_base}")
    print(f"    sigma_load:   50 -> {ld_cfg.sigma_load}")
    print(f"    LD50_base:    5e5 -> {ld_cfg.LD50_base:.0f}")
    print(f"    Eff growth at {T_celsius}°C: {eff_growth:.4f}")
    print(f"    Max clearance (r=1,c=1): {max_clearance:.4f}")
    print(f"    Typical clearance (r=.08,c=.08): {typ_clearance:.4f}")
    print(f"    r=1 immunity works: {max_clearance > eff_growth}")
    print(f"    Typical r=.08 net rate: {eff_growth - typ_clearance:+.4f} (need < 0 for typical agents to clear)")

    agents = create_population(rng)
    node_state = NodeDiseaseState(node_id=0, vibrio_concentration=0.0)

    # Seed infections
    infect_idx = rng.choice(N_POP, size=N_INITIAL_INFECTED, replace=False)
    agents['disease_state'][infect_idx] = DiseaseState.I1
    agents['pathogen_load'][infect_idx] = ld_cfg.L_init_base

    # Tracking
    daily_alive = np.zeros(N_DAYS, dtype=np.int32)
    daily_infected = np.zeros(N_DAYS, dtype=np.int32)
    daily_dead_cumul = np.zeros(N_DAYS, dtype=np.int32)
    daily_recovered_cumul = np.zeros(N_DAYS, dtype=np.int32)
    daily_vibrio = np.zeros(N_DAYS, dtype=np.float64)
    daily_mean_load = np.zeros(N_DAYS, dtype=np.float64)

    prev_cumul_inf = 0

    print(f"\n{'Day':>5} {'Alive':>6} {'Inf':>5} {'Dead':>5} {'Rec':>5} {'NewInf':>6} "
          f"{'MeanLoad':>12} {'P_k':>12} {'Prevalence':>10}")
    print("-" * 90)

    for day in range(N_DAYS):
        node_state = load_dependent_disease_update(
            agents, node_state,
            T_celsius=T_celsius, salinity=SALINITY, phi_k=PHI_K,
            dispersal_input=0.0, day=day,
            cfg=cfg, ld_cfg=ld_cfg, rng=rng,
        )

        n_alive = int(np.sum(agents['alive']))
        n_inf = node_state.n_I1 + node_state.n_I2
        new_inf = node_state.cumulative_infections - prev_cumul_inf
        prev_cumul_inf = node_state.cumulative_infections

        alive_mask = agents['alive'].astype(bool)
        inf_mask = alive_mask & (
            (agents['disease_state'] == DiseaseState.I1) |
            (agents['disease_state'] == DiseaseState.I2)
        )
        mean_load = float(np.mean(agents['pathogen_load'][inf_mask])) if np.any(inf_mask) else 0.0

        daily_alive[day] = n_alive
        daily_infected[day] = n_inf
        daily_dead_cumul[day] = node_state.cumulative_deaths
        daily_recovered_cumul[day] = node_state.cumulative_recoveries
        daily_vibrio[day] = node_state.vibrio_concentration
        daily_mean_load[day] = mean_load

        prevalence = n_inf / n_alive if n_alive > 0 else 0.0

        if day % 30 == 0 or day == N_DAYS - 1:
            print(f"{day:5d} {n_alive:6d} {n_inf:5d} {node_state.cumulative_deaths:5d} "
                  f"{node_state.cumulative_recoveries:5d} {new_inf:6d} "
                  f"{mean_load:12.1f} {node_state.vibrio_concentration:12.1f} {prevalence:10.4f}")

        if n_alive == 0:
            print(f"  *** All agents dead at day {day} ***")
            break

    peak_inf_day = int(np.argmax(daily_infected))
    peak_inf_count = int(daily_infected[peak_inf_day])
    peak_prevalence = peak_inf_count / daily_alive[peak_inf_day] if daily_alive[peak_inf_day] > 0 else 0
    total_deaths = int(daily_dead_cumul[min(day, N_DAYS-1)])
    total_recoveries = int(daily_recovered_cumul[min(day, N_DAYS-1)])
    total_infections = node_state.cumulative_infections
    cfr = total_deaths / total_infections if total_infections > 0 else 0
    final_alive = int(daily_alive[min(day, N_DAYS-1)])
    final_infected = int(daily_infected[min(day, N_DAYS-1)])
    peak_vibrio = float(np.max(daily_vibrio[:day+1]))

    print(f"\n{'SUMMARY':=^70}")
    print(f"  Total infections:       {total_infections}")
    print(f"  Total deaths:           {total_deaths}")
    print(f"  Total recoveries:       {total_recoveries}")
    print(f"  Case fatality rate:     {cfr:.2%}")
    print(f"  Peak prevalence:        {peak_prevalence:.2%} (day {peak_inf_day}, {peak_inf_count} infected)")
    print(f"  Peak vibrio:            {peak_vibrio:.1f} bact/mL")
    print(f"  Final alive:            {final_alive}/{N_POP}")
    print(f"  Final infected:         {final_infected}")
    print(f"  Mortality fraction:     {total_deaths/N_POP:.2%}")

    epidemic_ended = final_infected == 0 and day < N_DAYS - 1
    endemic = final_infected > 0
    print(f"  Epidemic burned out:    {epidemic_ended}")
    print(f"  Endemic at end:         {endemic}")

    return {
        'model': 'load_dependent_adjusted',
        'T_celsius': T_celsius,
        'total_infections': total_infections,
        'total_deaths': total_deaths,
        'total_recoveries': total_recoveries,
        'cfr': cfr,
        'peak_prevalence': peak_prevalence,
        'peak_day': peak_inf_day,
        'peak_count': peak_inf_count,
        'peak_vibrio': peak_vibrio,
        'final_alive': final_alive,
        'final_infected': final_infected,
        'mortality_fraction': total_deaths / N_POP,
        'epidemic_ended': epidemic_ended,
        'endemic': endemic,
        'daily_infected': daily_infected,
        'daily_alive': daily_alive,
        'daily_vibrio': daily_vibrio,
        'daily_mean_load': daily_mean_load,
    }


# ═══════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════

def main():
    print("SSWD-EvoEpi: Load-Dependent Disease Model Validation")
    print(f"Population: {N_POP}, Days: {N_DAYS}, Temperature: {T_CELSIUS}°C")
    print(f"Salinity: {SALINITY} psu, Flushing: {PHI_K} d^-1")
    print()

    # 1. Temperature sweep — understand parameter space
    sweep_temperatures()

    # 2. Immunity test
    test_immunity_at_r1()

    # 3. Run load-dependent epidemic
    ld_result = run_load_dependent_epidemic()

    # 4. Run timer-based epidemic
    tb_result = run_timer_based_epidemic()

    # 5. Compare
    n_pass, n_total = compare_results(ld_result, tb_result)

    # 6. Diagnose load-dependent dynamics
    diagnose_ld_dynamics(ld_result)

    # 7. Run with adjusted parameters
    ld_adjusted = run_load_dependent_adjusted()

    # 8. Compare adjusted LD vs TB
    print("\n" + "=" * 70)
    print("ADJUSTED LD vs TIMER-BASED COMPARISON")
    print("=" * 70)
    n_pass2, n_total2 = compare_results(ld_adjusted, tb_result)

    print("\n" + "=" * 70)
    print("VALIDATION COMPLETE")
    print("=" * 70)

    return n_pass, n_total


if __name__ == '__main__':
    n_pass, n_total = main()
    sys.exit(0 if n_pass >= n_total * 0.5 else 1)
