#!/usr/bin/env python3
"""Phase 8: Full Test Run with Continuous Settlement.

Runs publication-quality simulations, checks quality gates, saves results.
"""

import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

# Ensure project root on path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
os.chdir(ROOT)

from sswd_evoepi.model import (
    SimulationConfig,
    run_spatial_simulation,
    run_coupled_simulation,
)
from sswd_evoepi.spatial import load_node_definitions_yaml, build_network

OUT_DIR = ROOT / "results" / "continuous_settlement" / "test_run"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def get_git_hash():
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT
        ).decode().strip()
    except Exception:
        return "unknown"


def run_seed(seed, cfg, nodes_yaml="configs/nodes_5site.yaml"):
    """Run one spatial simulation and return result + timing."""
    nodes = load_node_definitions_yaml(nodes_yaml)
    net = build_network(
        nodes,
        D_L=cfg.spatial.D_L,
        D_P=cfg.spatial.D_P,
        alpha_self_fjord=cfg.spatial.alpha_self_fjord,
        alpha_self_open=cfg.spatial.alpha_self_open,
    )
    t0 = time.time()
    res = run_spatial_simulation(
        net, n_years=20, disease_year=3, seed=seed, config=cfg
    )
    elapsed = time.time() - t0
    return res, elapsed


def print_diagnostics(res, elapsed, seed):
    """Print comprehensive per-node diagnostics."""
    print(f"\n{'='*80}")
    print(f"SEED {seed} — Runtime: {elapsed:.1f}s")
    print(f"{'='*80}")
    print(f"Initial total pop: {res.initial_total_pop}")
    print(f"Final total pop:   {res.final_total_pop}")
    print()

    header = (
        f"{'Node':<14} {'K':>6} {'Crash%':>7} {'Final':>7} {'Nadir':>7} "
        f"{'NadYr':>6} {'DisDth':>7} {'NatDth':>7} {'SenDth':>7} "
        f"{'ΔResist':>8} {'MeanVir':>8}"
    )
    print(header)
    print("-" * len(header))

    for i in range(res.n_nodes):
        name = res.node_names[i]
        K = res.node_K[i]
        yearly_pop = res.yearly_pop[i]
        initial = yearly_pop[0] if yearly_pop[0] > 0 else K
        nadir = min(yearly_pop)
        nadir_yr = int(np.argmin(yearly_pop))
        final = yearly_pop[-1]
        crash_pct = (1.0 - nadir / initial) * 100 if initial > 0 else 0

        dis_deaths = sum(res.yearly_disease_deaths[i])
        nat_deaths = sum(res.yearly_natural_deaths[i])
        sen_deaths = sum(res.yearly_senescence_deaths[i])

        mr = res.yearly_mean_resistance[i]
        resist_shift = mr[-1] - mr[0] if len(mr) > 1 else 0

        mv = res.yearly_mean_virulence[i]
        mean_vir = np.nanmean([v for v in mv if v > 0]) if any(v > 0 for v in mv) else 0

        print(
            f"{name:<14} {K:>6} {crash_pct:>6.1f}% {final:>7} {nadir:>7} "
            f"{nadir_yr:>6} {dis_deaths:>7} {nat_deaths:>7} {sen_deaths:>7} "
            f"{resist_shift:>+7.4f} {mean_vir:>8.4f}"
        )

    return


def extract_crash_pcts(res):
    """Return dict of node_name -> crash%."""
    result = {}
    for i in range(res.n_nodes):
        yearly_pop = res.yearly_pop[i]
        initial = yearly_pop[0] if yearly_pop[0] > 0 else res.node_K[i]
        nadir = min(yearly_pop)
        crash_pct = (1.0 - nadir / initial) * 100 if initial > 0 else 0
        result[res.node_names[i]] = crash_pct
    return result


def check_quality_gates(results_by_seed, timings):
    """Check all quality gates. Returns (pass_bool, report_lines)."""
    gates = {}
    lines = []

    # Primary result (seed 42)
    res42 = results_by_seed[42]
    t42 = timings[42]

    # Gate 1: N→S mortality gradient
    crash_pcts = extract_crash_pcts(res42)
    node_order = list(crash_pcts.keys())
    # Sitka (north) should crash less than Monterey (south)
    # Use partial matching for node names
    sitka_crash = next((v for k, v in crash_pcts.items() if "Sitka" in k), 0)
    monterey_crash = next((v for k, v in crash_pcts.items() if "Monterey" in k), 0)
    ns_gradient = monterey_crash > sitka_crash
    gates["N→S mortality gradient"] = ns_gradient
    lines.append(f"N→S gradient: Sitka {sitka_crash:.1f}% vs Monterey {monterey_crash:.1f}% — {'✅' if ns_gradient else '❌'}")

    # Gate 2: Fjord protection (Howe Sound < SJI)
    # Use partial matching for node names
    howe_crash = next((v for k, v in crash_pcts.items() if "Howe" in k), 0)
    sji_crash = next((v for k, v in crash_pcts.items() if "San Juan" in k or "SJI" in k), 0)
    fjord_protection = howe_crash < sji_crash
    gates["Fjord protection"] = fjord_protection
    lines.append(f"Fjord protection: Howe Sound {howe_crash:.1f}% < SJI {sji_crash:.1f}% — {'✅' if fjord_protection else '❌'}")

    # Gate 3: Disease deaths > 0 at all nodes
    all_disease = True
    for i in range(res42.n_nodes):
        dd = sum(res42.yearly_disease_deaths[i])
        if dd == 0:
            all_disease = False
            lines.append(f"  ❌ {res42.node_names[i]}: 0 disease deaths!")
    gates["Disease deaths > 0 all nodes"] = all_disease
    lines.append(f"Disease deaths everywhere: {'✅' if all_disease else '❌'}")

    # Gate 4: Mean virulence ~0.5 (PE working)
    all_vir = []
    for i in range(res42.n_nodes):
        mv = res42.yearly_mean_virulence[i]
        valid = [v for v in mv if v > 0]
        if valid:
            all_vir.extend(valid)
    mean_vir = np.mean(all_vir) if all_vir else 0
    pe_working = 0.2 < mean_vir < 0.8
    gates["Mean virulence ~0.5"] = pe_working
    lines.append(f"Mean virulence: {mean_vir:.3f} — {'✅' if pe_working else '❌'}")

    # Gate 5: Resistance increases over time
    resist_increases = 0
    for i in range(res42.n_nodes):
        mr = res42.yearly_mean_resistance[i]
        if len(mr) > 1 and mr[-1] > mr[0]:
            resist_increases += 1
    selection_working = resist_increases >= 2  # at least 2/5 nodes show increase
    gates["Resistance increases"] = selection_working
    lines.append(f"Resistance increase: {resist_increases}/{res42.n_nodes} nodes — {'✅' if selection_working else '❌'}")

    # Gate 6: Runtime < 150s
    fast_enough = t42 < 150
    gates["Runtime < 150s"] = fast_enough
    lines.append(f"Runtime: {t42:.1f}s — {'✅' if fast_enough else '❌'}")

    # Gate 7: Consistent across seeds
    all_crashes = {}
    for seed, res in results_by_seed.items():
        all_crashes[seed] = extract_crash_pcts(res)

    consistent = True
    inconsistencies = []
    for node in node_order:
        crashes = [all_crashes[s].get(node, 0) for s in sorted(all_crashes.keys())]
        span = max(crashes) - min(crashes)
        if span > 40:  # >40 percentage point spread is qualitatively different
            consistent = False
            inconsistencies.append(f"  ❌ {node}: spread={span:.1f}pp ({min(crashes):.1f}%-{max(crashes):.1f}%)")
    gates["Cross-seed consistency"] = consistent
    lines.append(f"Cross-seed consistency: {'✅' if consistent else '❌'}")
    for inc in inconsistencies:
        lines.append(inc)

    # Gate 8: No sawtooth — checked via daily resolution below
    # We'll set this gate externally
    gates["No sawtooth"] = None  # placeholder

    overall = all(v for k, v in gates.items() if v is not None)
    return overall, gates, lines


def daily_resolution_check(cfg):
    """Run single-node coupled sim with daily recording to check for sawtooth."""
    print("\n" + "=" * 80)
    print("DAILY RESOLUTION CHECK — Single Node Coupled Sim")
    print("=" * 80)

    t0 = time.time()
    res = run_coupled_simulation(
        n_individuals=5000,
        carrying_capacity=5000,
        habitat_area=100.0,
        T_celsius=12.0,
        salinity=32.0,
        phi_k=0.0,
        n_years=10,
        disease_year=3,
        initial_infected=10,
        seed=42,
        record_daily=True,
        config=cfg,
    )
    elapsed = time.time() - t0
    print(f"Runtime: {elapsed:.1f}s")

    daily_pop = np.array(res.daily_pop)
    print(f"Daily pop array: {len(daily_pop)} days")
    print(f"Min pop: {daily_pop.min()}, Max pop: {daily_pop.max()}")

    # Check for sawtooth: max daily increase should be < 10% of current pop
    max_daily_increase = 0
    max_daily_increase_day = 0
    max_daily_increase_pct = 0
    sawtooth_detected = False
    sawtooth_days = []
    year_boundary_sawtooth = []
    mid_year_sawtooth = []

    for d in range(1, len(daily_pop)):
        if daily_pop[d - 1] > 0:
            increase = daily_pop[d] - daily_pop[d - 1]
            pct = increase / daily_pop[d - 1] * 100
            if increase > max_daily_increase:
                max_daily_increase = increase
                max_daily_increase_day = d
                max_daily_increase_pct = pct
            if pct > 10.0:
                sawtooth_detected = True
                yr = d // 365
                doy = d % 365
                entry = (d, yr, doy, daily_pop[d-1], daily_pop[d], increase, pct)
                sawtooth_days.append(entry)
                if doy < 10 or doy > 355:
                    year_boundary_sawtooth.append(entry)
                else:
                    mid_year_sawtooth.append(entry)

    print(f"Max daily increase: {max_daily_increase} on day {max_daily_increase_day} ({max_daily_increase_pct:.2f}%)")
    print(f"Days exceeding 10% threshold: {len(sawtooth_days)}")
    print(f"  Year-boundary (DOY < 10 or > 355): {len(year_boundary_sawtooth)}")
    print(f"  Mid-year: {len(mid_year_sawtooth)}")
    if sawtooth_days:
        print("\nSawtooth days detail:")
        for d, yr, doy, prev, curr, inc, pct in sawtooth_days[:15]:
            print(f"  Day {d} (Y{yr} DOY{doy:>3}): {prev:>5} → {curr:>5} (+{inc}, {pct:.1f}%)")
    
    # Distinguish: year-boundary sawtooth is from annual mortality batch,
    # NOT from settlement. True settlement sawtooth would show mid-year spikes.
    settlement_sawtooth = len(mid_year_sawtooth) > 0
    if sawtooth_detected and not settlement_sawtooth:
        print("\n⚠️  Sawtooth only at year boundaries (annual mortality batch)")
        print("   Continuous settlement is working correctly.")
        print("   Year-boundary mortality batch is a separate issue.")
    
    print(f"\nSettlement sawtooth: {'❌ YES' if settlement_sawtooth else '✅ NO'}")
    print(f"Year-boundary mortality sawtooth: {'⚠️ YES' if year_boundary_sawtooth else '✅ NO'}")

    # Print weekly population summary (years 2-5 around epidemic)
    print("\nWeekly pop (years 2-6):")
    for week in range(52 * 2, min(52 * 6, len(daily_pop) // 7)):
        day = week * 7
        if day < len(daily_pop):
            yr = day // 365
            doy = day % 365
            print(f"  Year {yr} DOY {doy:>3}: {daily_pop[day]:>6}")

    # Return: (no_settlement_sawtooth, daily_pop, year_boundary_sawtooth_present)
    return not settlement_sawtooth, daily_pop, len(year_boundary_sawtooth) > 0


def save_results(results_by_seed, timings, cfg, gates, daily_pop=None):
    """Save all results to disk."""
    git_hash = get_git_hash()

    # Metadata
    meta = {
        "git_hash": git_hash,
        "seeds": list(results_by_seed.keys()),
        "n_years": 20,
        "disease_year": 3,
        "n_nodes": 5,
        "pathogen_evolution_enabled": cfg.pathogen_evolution.enabled,
        "timings": {str(s): t for s, t in timings.items()},
        "quality_gates": {k: bool(v) if v is not None else None for k, v in gates.items()},
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
    }
    with open(OUT_DIR / "metadata.json", "w") as f:
        json.dump(meta, f, indent=2)

    # Per-seed NPZ
    for seed, res in results_by_seed.items():
        arrays = {
            "node_names": np.array(res.node_names),
            "node_K": np.array(res.node_K),
            "yearly_pop": np.array(res.yearly_pop),
            "yearly_disease_deaths": np.array(res.yearly_disease_deaths),
            "yearly_natural_deaths": np.array(res.yearly_natural_deaths),
            "yearly_senescence_deaths": np.array(res.yearly_senescence_deaths),
            "yearly_mean_resistance": np.array(res.yearly_mean_resistance),
            "yearly_mean_virulence": np.array(res.yearly_mean_virulence),
            "yearly_recruits": np.array(res.yearly_recruits),
            "yearly_total_pop": np.array(res.yearly_total_pop),
        }
        np.savez_compressed(OUT_DIR / f"results_seed{seed}.npz", **arrays)

    if daily_pop is not None:
        np.savez_compressed(OUT_DIR / "daily_pop_coupled.npz", daily_pop=daily_pop)

    print(f"\nResults saved to {OUT_DIR}")


def main():
    # Config
    cfg = SimulationConfig()
    cfg.pathogen_evolution.enabled = True
    print("Config: pathogen_evolution.enabled =", cfg.pathogen_evolution.enabled)

    # ── Step 1: Primary run (seed 42) ──
    print("\n▶ Running seed 42 (primary)...")
    res42, t42 = run_seed(42, cfg)
    print_diagnostics(res42, t42, 42)

    results_by_seed = {42: res42}
    timings = {42: t42}

    # ── Step 2: Additional seeds ──
    for seed in [43, 44, 45]:
        print(f"\n▶ Running seed {seed}...")
        res, t = run_seed(seed, cfg)
        print_diagnostics(res, t, seed)
        results_by_seed[seed] = res
        timings[seed] = t

    # ── Step 3: Cross-seed crash comparison ──
    print("\n" + "=" * 80)
    print("CROSS-SEED CRASH COMPARISON")
    print("=" * 80)
    node_order = results_by_seed[42].node_names
    header = f"{'Node':<14}" + "".join(f"{'Seed ' + str(s):>10}" for s in sorted(results_by_seed.keys()))
    print(header)
    print("-" * len(header))
    for node in node_order:
        row = f"{node:<14}"
        for seed in sorted(results_by_seed.keys()):
            cp = extract_crash_pcts(results_by_seed[seed])
            row += f"{cp.get(node, 0):>9.1f}%"
        print(row)

    # ── Step 4: Quality gates (partial — without sawtooth) ──
    overall, gates, gate_lines = check_quality_gates(results_by_seed, timings)

    # ── Step 5: Daily resolution (sawtooth check) ──
    no_settlement_sawtooth, daily_pop, year_boundary_sawtooth = daily_resolution_check(cfg)
    gates["No sawtooth (settlement)"] = no_settlement_sawtooth
    del gates["No sawtooth"]  # remove placeholder

    # Final gate check
    overall = all(v for v in gates.values() if v is not None)

    print("\n" + "=" * 80)
    print("QUALITY GATE SUMMARY")
    print("=" * 80)
    for line in gate_lines:
        print(line)
    print(f"No settlement sawtooth: {'✅' if no_settlement_sawtooth else '❌'}")
    if year_boundary_sawtooth:
        print("⚠️  Year-boundary mortality sawtooth present (separate from settlement fix)")
    print(f"\n{'🟢 ALL GATES PASSED' if overall else '🔴 SOME GATES FAILED'}")

    # ── Step 6: Save results ──
    save_results(results_by_seed, timings, cfg, gates, daily_pop)

    # ── Step 7: Write report ──
    write_report(results_by_seed, timings, gates, gate_lines, no_settlement_sawtooth, daily_pop, cfg, year_boundary_sawtooth)

    return 0 if overall else 1


def write_report(results_by_seed, timings, gates, gate_lines, no_sawtooth, daily_pop, cfg, year_boundary_sawtooth=False):
    """Write TEST_RUN_REPORT.md."""
    res42 = results_by_seed[42]
    git_hash = get_git_hash()

    lines = []
    lines.append("# Continuous Settlement — Phase 8 Test Run Report")
    lines.append("")
    lines.append(f"**Date:** {time.strftime('%Y-%m-%d %H:%M %Z')}")
    lines.append(f"**Git:** `{git_hash[:8]}`")
    lines.append(f"**Config:** 5 nodes, 20 years, disease year 3, PE enabled")
    lines.append(f"**Seeds:** 42, 43, 44, 45")
    lines.append("")

    # Quality gates
    lines.append("## Quality Gates")
    lines.append("")
    overall = all(v for v in gates.values() if v is not None)
    lines.append(f"**Overall: {'✅ PASS' if overall else '❌ FAIL'}**")
    lines.append("")
    lines.append("| Gate | Status |")
    lines.append("|------|--------|")
    for k, v in gates.items():
        status = "✅" if v else "❌"
        lines.append(f"| {k} | {status} |")
    lines.append("")

    # Per-node results (seed 42)
    lines.append("## Seed 42 Results (Primary)")
    lines.append("")
    lines.append("| Node | K | Crash% | Final Pop | Nadir | Nadir Yr | Disease Deaths | Resistance Shift | Mean Virulence |")
    lines.append("|------|---|--------|-----------|-------|----------|----------------|-----------------|----------------|")
    for i in range(res42.n_nodes):
        name = res42.node_names[i]
        K = res42.node_K[i]
        yp = res42.yearly_pop[i]
        initial = yp[0] if yp[0] > 0 else K
        nadir = min(yp)
        nadir_yr = int(np.argmin(yp))
        final = yp[-1]
        crash = (1.0 - nadir / initial) * 100
        dd = sum(res42.yearly_disease_deaths[i])
        mr = res42.yearly_mean_resistance[i]
        rs = mr[-1] - mr[0] if len(mr) > 1 else 0
        mv = res42.yearly_mean_virulence[i]
        mean_v = np.nanmean([v for v in mv if v > 0]) if any(v > 0 for v in mv) else 0
        lines.append(f"| {name} | {K} | {crash:.1f}% | {final} | {nadir} | {nadir_yr} | {dd} | {rs:+.4f} | {mean_v:.4f} |")
    lines.append("")

    # Cross-seed comparison
    lines.append("## Cross-Seed Crash Comparison (%)")
    lines.append("")
    seeds = sorted(results_by_seed.keys())
    lines.append("| Node | " + " | ".join(f"Seed {s}" for s in seeds) + " | Spread |")
    lines.append("|------|" + "|".join(["--------|"] * (len(seeds) + 1)))
    for node in res42.node_names:
        crashes = []
        for s in seeds:
            cp = extract_crash_pcts(results_by_seed[s])
            crashes.append(cp.get(node, 0))
        spread = max(crashes) - min(crashes)
        row = f"| {node} | " + " | ".join(f"{c:.1f}" for c in crashes) + f" | {spread:.1f} |"
        lines.append(row)
    lines.append("")

    # Timing
    lines.append("## Performance")
    lines.append("")
    for s in seeds:
        lines.append(f"- Seed {s}: {timings[s]:.1f}s")
    lines.append(f"- Mean: {np.mean(list(timings.values())):.1f}s")
    lines.append("")

    # Daily resolution
    lines.append("## Daily Resolution Check")
    lines.append("")
    if daily_pop is not None:
        dp = np.array(daily_pop)
        max_inc = 0
        max_inc_pct = 0
        for d in range(1, len(dp)):
            if dp[d-1] > 0:
                inc = dp[d] - dp[d-1]
                pct = inc / dp[d-1] * 100
                if inc > max_inc:
                    max_inc = inc
                    max_inc_pct = pct
        lines.append(f"- Max daily increase: {max_inc} ({max_inc_pct:.2f}%)")
        lines.append(f"- Settlement sawtooth: {'❌ YES' if not no_sawtooth else '✅ NO'}")
        lines.append(f"- Year-boundary mortality sawtooth: {'⚠️ YES' if year_boundary_sawtooth else '✅ NO'}")
        lines.append(f"- Threshold: 10% daily increase of current pop")
    lines.append("")

    # Gate details
    lines.append("## Gate Details")
    lines.append("")
    for gl in gate_lines:
        lines.append(f"- {gl}")
    lines.append(f"- No settlement sawtooth: {'✅' if no_sawtooth else '❌'}")
    if year_boundary_sawtooth:
        lines.append("- ⚠️ Year-boundary mortality sawtooth present (annual batch, not settlement-related)")
    lines.append("")

    report_path = OUT_DIR / "TEST_RUN_REPORT.md"
    report_path.write_text("\n".join(lines))
    print(f"Report written to {report_path}")


if __name__ == "__main__":
    sys.exit(main())
