#!/usr/bin/env python3
"""BUILD PHASE 10: Full 5-node integrated 20-year epidemic simulation.

Runs a 20-year simulation:
  - Years 0-4: disease-free spinup (verify demographic equilibrium)
  - Year 5: introduce disease (ubiquitous scenario: seeded everywhere)
  - Years 5-20: epidemic + recovery dynamics

Records full output, runs verification checks, generates text summary.

Usage:
    python3 run_epidemic_20yr.py
"""

import json
import os
import sys
import time
from pathlib import Path

import numpy as np

# Add project root to path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

from sswd_evoepi.config import default_config
from sswd_evoepi.model import (
    SpatialSimResult,
    make_effect_sizes,
    run_spatial_simulation,
)
from sswd_evoepi.spatial import make_5node_network, get_5node_definitions
from sswd_evoepi.types import N_LOCI, IDX_EF1A, N_ADDITIVE


# ═══════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════

N_YEARS = 20
DISEASE_YEAR = 5           # Introduce disease at year 5
INITIAL_INFECTED = 5       # Per node
SEED = 42
OUTPUT_DIR = project_root / "results" / "5node_epidemic_20yr"


# ═══════════════════════════════════════════════════════════════════════
# RUN SIMULATION
# ═══════════════════════════════════════════════════════════════════════

def run_simulation():
    """Run the full 20-year epidemic simulation."""
    print("=" * 72)
    print("SSWD-EvoEpi: 20-Year Epidemic Simulation (5-node network)")
    print("=" * 72)
    print()

    config = default_config()
    network = make_5node_network(seed=SEED)

    print(f"Network: {network.n_nodes} nodes")
    for node in network.nodes:
        nd = node.definition
        print(f"  [{nd.node_id}] {nd.name}: K={nd.carrying_capacity}, "
              f"SST={nd.mean_sst}°C, φ={nd.flushing_rate:.3f}, "
              f"sal={nd.salinity:.0f}psu, {'fjord' if nd.is_fjord else 'coast'}")
    print()
    print(f"Simulation: {N_YEARS} years, disease at year {DISEASE_YEAR}")
    print(f"Seed: {SEED}")
    print()

    t0 = time.time()

    def progress(year, total):
        elapsed = time.time() - t0
        pct = (year + 1) / total * 100
        rate = elapsed / (year + 1) if year > 0 else 0
        eta = rate * (total - year - 1)
        print(f"  Year {year:3d}/{total} ({pct:5.1f}%)  "
              f"elapsed={elapsed:.0f}s  eta={eta:.0f}s", flush=True)

    result = run_spatial_simulation(
        network=network,
        n_years=N_YEARS,
        disease_year=DISEASE_YEAR,
        initial_infected_per_node=INITIAL_INFECTED,
        seed=SEED,
        config=config,
        progress_callback=progress,
    )

    elapsed = time.time() - t0
    print()
    print(f"Simulation complete in {elapsed:.1f}s")
    print()

    return result, config


# ═══════════════════════════════════════════════════════════════════════
# VERIFICATION CHECKS
# ═══════════════════════════════════════════════════════════════════════

def verify_results(result: SpatialSimResult):
    """Run verification checks on the simulation results.

    Returns (passed, failed, warnings) lists of check descriptions.
    """
    passed = []
    failed = []
    warnings = []

    N = result.n_nodes
    DY = DISEASE_YEAR
    names = result.node_names
    K_vals = result.node_K

    # ── 1. Pre-epidemic: all nodes near K ────────────────────────────
    print("─── Verification Checks ───")
    print()

    # Check year before disease (year DY-1)
    pre_year = DY - 1
    for i in range(N):
        pop = result.yearly_pop[i, pre_year]
        K = K_vals[i]
        ratio = pop / K if K > 0 else 0
        if 0.4 < ratio < 1.5:
            passed.append(f"[Node {i} {names[i]}] Pre-epidemic pop={pop} "
                         f"({ratio:.1%} of K={K}) ✓")
        else:
            failed.append(f"[Node {i} {names[i]}] Pre-epidemic pop={pop} "
                         f"({ratio:.1%} of K={K}) — expected 40-150%")

    # ── 2. Epidemic onset: disease spreads, significant mortality ────
    # Check years DY+1 to DY+3 for disease deaths
    for i in range(N):
        total_dd = sum(result.yearly_disease_deaths[i, DY:min(DY+4, N_YEARS)])
        pre_pop = result.yearly_pop[i, DY]
        mort_frac = total_dd / pre_pop if pre_pop > 0 else 0
        if mort_frac > 0.10:
            passed.append(f"[Node {i} {names[i]}] Epidemic mortality={mort_frac:.1%} "
                         f"({total_dd} deaths in 3yr) ✓")
        else:
            warnings.append(f"[Node {i} {names[i]}] Low epidemic mortality={mort_frac:.1%} "
                           f"({total_dd} deaths) — may need parameter tuning")

    # ── 3. North-south temperature gradient → Monterey worst ─────────
    # Node 4 (Monterey, 14°C) should have highest mortality
    # Node 0 (Sitka, 8°C) should have lowest
    total_dd = np.zeros(N)
    for i in range(N):
        total_dd[i] = result.yearly_disease_deaths[i, DY:].sum()

    worst_node = int(np.argmax(total_dd))
    least_node = int(np.argmin(total_dd))

    if worst_node == 4:  # Monterey
        passed.append(f"Monterey (warmest) has highest disease mortality ✓")
    else:
        warnings.append(f"Expected Monterey as worst-affected, got {names[worst_node]}")

    if least_node == 0:  # Sitka
        passed.append(f"Sitka (coldest) has lowest disease mortality ✓")
    elif least_node == 1:  # Howe Sound (fjord protection)
        passed.append(f"Howe Sound (fjord) has lowest disease mortality ✓")
    else:
        warnings.append(f"Expected Sitka/Howe Sound as least-affected, got {names[least_node]}")

    # ── 4. Fjord protection: Howe Sound lower mortality ──────────────
    howe_dd = total_dd[1]
    sji_dd = total_dd[2]  # Same latitude, non-fjord
    if howe_dd <= sji_dd:
        passed.append(f"Howe Sound ({howe_dd:.0f} deaths) ≤ SJI ({sji_dd:.0f}) — "
                     f"fjord protection ✓")
    else:
        warnings.append(f"Howe Sound ({howe_dd:.0f} deaths) > SJI ({sji_dd:.0f}) — "
                       f"fjord protection not evident")

    # ── 5. Genetic selection: mean resistance increases post-epidemic ─
    for i in range(N):
        r_before = result.yearly_mean_resistance[i, DY] if DY < N_YEARS else 0
        last_yr = min(N_YEARS - 1, DY + 10)
        r_after = result.yearly_mean_resistance[i, last_yr]
        delta_r = r_after - r_before
        pop_after = result.yearly_pop[i, last_yr]
        if pop_after < 5:
            warnings.append(f"[Node {i} {names[i]}] Pop too low ({pop_after}) "
                           f"for meaningful resistance check")
        elif delta_r > 0:
            passed.append(f"[Node {i} {names[i]}] Mean resistance increased: "
                         f"{r_before:.4f} → {r_after:.4f} (Δ={delta_r:+.4f}) ✓")
        else:
            warnings.append(f"[Node {i} {names[i]}] Resistance did not increase: "
                           f"{r_before:.4f} → {r_after:.4f} (Δ={delta_r:+.4f})")

    # ── 6. Population minimum during epidemic ────────────────────────
    for i in range(N):
        min_pop = int(np.min(result.yearly_pop[i, DY:]))
        pre_pop = result.yearly_pop[i, DY]
        crash_frac = 1.0 - min_pop / pre_pop if pre_pop > 0 else 0
        if crash_frac > 0.30:
            passed.append(f"[Node {i} {names[i]}] Population crash: "
                         f"{crash_frac:.0%} decline (min={min_pop}) ✓")
        else:
            warnings.append(f"[Node {i} {names[i]}] Modest crash: "
                           f"{crash_frac:.0%} (min={min_pop})")

    # ── 7. No node goes to exact zero (extinction check) ────────────
    for i in range(N):
        min_pop = int(np.min(result.yearly_pop[i, :]))
        if min_pop == 0:
            warnings.append(f"[Node {i} {names[i]}] Went EXTINCT (pop=0)")
        else:
            passed.append(f"[Node {i} {names[i]}] Never fully extinct (min={min_pop}) ✓")

    # ── 8. Peak disease prevalence ───────────────────────────────────
    if result.peak_disease_prevalence is not None:
        for i in range(N):
            peak = result.peak_disease_prevalence[i]
            if peak > 0.01:
                passed.append(f"[Node {i} {names[i]}] Peak disease prevalence: "
                             f"{peak:.1%} ✓")
            else:
                warnings.append(f"[Node {i} {names[i]}] Very low peak prevalence: "
                               f"{peak:.1%}")

    # ── 9. Allele frequency shifts at top loci ───────────────────────
    if (result.pre_epidemic_allele_freq is not None
            and result.post_epidemic_allele_freq is not None):
        for i in range(N):
            pre_q = result.pre_epidemic_allele_freq[i, :3]
            post_q = result.post_epidemic_allele_freq[i, :3]
            delta = post_q - pre_q
            if np.any(delta > 0.001):
                passed.append(f"[Node {i} {names[i]}] Top-locus Δq: "
                             f"{delta[0]:+.4f}, {delta[1]:+.4f}, {delta[2]:+.4f} ✓")
            else:
                warnings.append(f"[Node {i} {names[i]}] No positive allele shift at top loci")

    # ── 10. No NaN or negative values ────────────────────────────────
    has_nan = False
    for arr_name, arr in [
        ('yearly_pop', result.yearly_pop),
        ('yearly_mean_resistance', result.yearly_mean_resistance),
        ('yearly_vibrio_max', result.yearly_vibrio_max),
    ]:
        if np.any(np.isnan(arr)):
            failed.append(f"NaN found in {arr_name}")
            has_nan = True
    if not has_nan:
        passed.append("No NaN values in output ✓")

    if np.any(result.yearly_pop < 0):
        failed.append("Negative population values found")
    else:
        passed.append("No negative population values ✓")

    return passed, failed, warnings


# ═══════════════════════════════════════════════════════════════════════
# SAVE RESULTS
# ═══════════════════════════════════════════════════════════════════════

def save_results(result: SpatialSimResult, output_dir: Path):
    """Save all simulation output to files."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save as npz (compressed NumPy archive)
    arrays = {}
    for attr in [
        'yearly_pop', 'yearly_adults', 'yearly_recruits',
        'yearly_natural_deaths', 'yearly_disease_deaths',
        'yearly_mean_resistance', 'yearly_vibrio_max',
        'yearly_ef1a_freq', 'yearly_va', 'yearly_allele_freq_top3',
        'yearly_ne_ratio', 'yearly_total_pop',
        'yearly_total_larvae_dispersed', 'peak_disease_prevalence',
    ]:
        val = getattr(result, attr, None)
        if val is not None:
            arrays[attr] = val

    if result.pre_epidemic_allele_freq is not None:
        arrays['pre_epidemic_allele_freq'] = result.pre_epidemic_allele_freq
    if result.post_epidemic_allele_freq is not None:
        arrays['post_epidemic_allele_freq'] = result.post_epidemic_allele_freq
    if result.node_K is not None:
        arrays['node_K'] = result.node_K

    np.savez_compressed(output_dir / "simulation_data.npz", **arrays)

    # Save metadata as JSON
    metadata = {
        'n_years': result.n_years,
        'n_nodes': result.n_nodes,
        'node_names': result.node_names,
        'disease_year': result.disease_year,
        'seed': result.seed,
        'initial_total_pop': result.initial_total_pop,
        'final_total_pop': result.final_total_pop,
    }
    with open(output_dir / "metadata.json", 'w') as f:
        json.dump(metadata, f, indent=2)

    print(f"Results saved to {output_dir}/")


# ═══════════════════════════════════════════════════════════════════════
# TEXT SUMMARY
# ═══════════════════════════════════════════════════════════════════════

def generate_summary(result: SpatialSimResult, passed, failed, warnings) -> str:
    """Generate human-readable text summary."""
    lines = []
    lines.append("=" * 72)
    lines.append("SSWD-EvoEpi: 20-Year Epidemic Simulation Results")
    lines.append("=" * 72)
    lines.append("")
    lines.append(f"Simulation: {result.n_years} years, {result.n_nodes} nodes")
    lines.append(f"Disease introduction: year {result.disease_year}")
    lines.append(f"Seed: {result.seed}")
    lines.append(f"Total initial population: {result.initial_total_pop}")
    lines.append(f"Total final population: {result.final_total_pop}")
    lines.append("")

    # Per-node summary table
    lines.append("─── Per-Node Summary ───")
    lines.append(f"{'Node':<25s} {'K':>5s} {'Pre-N':>6s} {'Min-N':>6s} "
                 f"{'Final-N':>7s} {'Peak Prev':>10s} {'Tot Deaths':>10s} "
                 f"{'Crash%':>7s}")
    lines.append("─" * 85)

    DY = result.disease_year or 0
    effect_sizes = make_effect_sizes(12345)

    for i in range(result.n_nodes):
        name = result.node_names[i]
        K = result.node_K[i]
        pre_n = result.yearly_pop[i, DY]
        min_n = int(np.min(result.yearly_pop[i, DY:]))
        final_n = result.yearly_pop[i, -1]
        peak_prev = result.peak_disease_prevalence[i] if result.peak_disease_prevalence is not None else 0
        tot_dd = int(result.yearly_disease_deaths[i, DY:].sum())
        crash = 1.0 - min_n / pre_n if pre_n > 0 else 0

        lines.append(f"{name:<25s} {K:>5d} {pre_n:>6d} {min_n:>6d} "
                     f"{final_n:>7d} {peak_prev:>9.1%} {tot_dd:>10d} "
                     f"{crash:>6.1%}")

    lines.append("")

    # Genetics summary
    lines.append("─── Genetics Summary ───")
    lines.append(f"{'Node':<25s} {'r̄ pre':>7s} {'r̄ post':>7s} {'r̄ yr20':>7s} "
                 f"{'Δr̄':>7s} {'EF1A pre':>9s} {'EF1A yr20':>10s}")
    lines.append("─" * 78)

    for i in range(result.n_nodes):
        name = result.node_names[i]
        r_pre = result.yearly_mean_resistance[i, DY]
        post_yr = min(DY + 3, result.n_years - 1)
        r_post = result.yearly_mean_resistance[i, post_yr]
        r_final = result.yearly_mean_resistance[i, -1]
        delta_r = r_final - r_pre
        ef1a_pre = result.yearly_ef1a_freq[i, DY]
        ef1a_final = result.yearly_ef1a_freq[i, -1]

        lines.append(f"{name:<25s} {r_pre:>7.4f} {r_post:>7.4f} {r_final:>7.4f} "
                     f"{delta_r:>+7.4f} {ef1a_pre:>9.4f} {ef1a_final:>10.4f}")

    lines.append("")

    # Allele frequency shifts at top 3 loci
    if (result.pre_epidemic_allele_freq is not None
            and result.post_epidemic_allele_freq is not None):
        lines.append("─── Allele Frequency Shifts (top 3 effect-size loci) ───")
        lines.append(f"{'Node':<25s} {'Δq₁':>8s} {'Δq₂':>8s} {'Δq₃':>8s}")
        lines.append("─" * 55)
        for i in range(result.n_nodes):
            name = result.node_names[i]
            pre_q = result.pre_epidemic_allele_freq[i, :3]
            post_q = result.post_epidemic_allele_freq[i, :3]
            delta = post_q - pre_q
            lines.append(f"{name:<25s} {delta[0]:>+8.4f} {delta[1]:>+8.4f} {delta[2]:>+8.4f}")
        lines.append("")

    # Functional extinction check
    lines.append("─── Functional Extinction Check ───")
    for i in range(result.n_nodes):
        name = result.node_names[i]
        min_pop = int(np.min(result.yearly_pop[i, :]))
        min_yr = int(np.argmin(result.yearly_pop[i, :]))
        K = result.node_K[i]
        if min_pop == 0:
            lines.append(f"  {name}: EXTINCT at year {min_yr}")
        elif min_pop < K * 0.01:
            lines.append(f"  {name}: FUNCTIONALLY EXTINCT (min={min_pop}, "
                        f"{min_pop/K:.1%} of K) at year {min_yr}")
        else:
            lines.append(f"  {name}: Surviving (min={min_pop}, "
                        f"{min_pop/K:.1%} of K) at year {min_yr}")
    lines.append("")

    # Year-by-year population trajectory
    lines.append("─── Year-by-Year Total Population ───")
    for yr in range(result.n_years):
        marker = " ← disease" if yr == DY else ""
        lines.append(f"  Year {yr:3d}: {result.yearly_total_pop[yr]:6d}{marker}")
    lines.append("")

    # Verification results
    lines.append("─── Verification Results ───")
    lines.append(f"  Passed:   {len(passed)}")
    lines.append(f"  Failed:   {len(failed)}")
    lines.append(f"  Warnings: {len(warnings)}")
    lines.append("")

    if failed:
        lines.append("FAILURES:")
        for f in failed:
            lines.append(f"  ✗ {f}")
        lines.append("")

    if warnings:
        lines.append("WARNINGS:")
        for w in warnings:
            lines.append(f"  ⚠ {w}")
        lines.append("")

    if passed:
        lines.append("PASSED:")
        for p in passed:
            lines.append(f"  ✓ {p}")
        lines.append("")

    return "\n".join(lines)


# ═══════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════

def main():
    # Run simulation
    result, config = run_simulation()

    # Verification
    passed, failed, warnings = verify_results(result)

    # Generate summary
    summary = generate_summary(result, passed, failed, warnings)
    print(summary)

    # Save results
    save_results(result, OUTPUT_DIR)

    # Save text summary
    summary_path = OUTPUT_DIR / "summary.txt"
    with open(summary_path, 'w') as f:
        f.write(summary)
    print(f"Summary saved to {summary_path}")

    # Exit code
    if failed:
        print(f"\n⚠ {len(failed)} FAILURES — see above")
        return 1
    else:
        print(f"\n✓ All checks passed ({len(passed)} passed, {len(warnings)} warnings)")
        return 0


if __name__ == "__main__":
    sys.exit(main())
