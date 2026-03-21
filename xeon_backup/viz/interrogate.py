"""CLI drill-down tool for SSWD-EvoEpi simulation results.

Inspect specific node/year combinations with detailed output.

Usage:
    python -m viz.interrogate results/5node_epidemic_20yr/ --node 2 --year 7
    python -m viz.interrogate results/5node_epidemic_20yr/ --node 2 --year 7 --individuals
    python -m viz.interrogate results/5node_epidemic_20yr/ --list-nodes
    python -m viz.interrogate results/5node_epidemic_20yr/ --summary
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from viz.plot_utils import SimulationData


# ═══════════════════════════════════════════════════════════════════════
# DISPLAY HELPERS
# ═══════════════════════════════════════════════════════════════════════

def _divider(title: str = "", width: int = 70) -> str:
    if title:
        pad = width - len(title) - 4
        return f"\n{'═' * 2} {title} {'═' * max(pad, 2)}"
    return "═" * width


def _bar(value: float, max_val: float, width: int = 30) -> str:
    """Simple text bar chart."""
    if max_val <= 0:
        return "░" * width
    filled = int(value / max_val * width)
    filled = min(filled, width)
    return "█" * filled + "░" * (width - filled)


# ═══════════════════════════════════════════════════════════════════════
# INTERROGATION FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════

def list_nodes(data: SimulationData) -> None:
    """Print available nodes."""
    print(_divider("Available Nodes"))
    for i, name in enumerate(data.node_names):
        K = int(data.node_K[i])
        final_pop = int(data.yearly_pop[i, -1])
        print(f"  [{i}] {name:30s}  K={K:5d}  Final N={final_pop:5d}")
    print(f"\n  Years: 0 — {data.n_years - 1}")
    print(f"  Disease year: {data.disease_year}")


def summary_all(data: SimulationData) -> None:
    """Print quick summary across all nodes."""
    print(_divider("Simulation Summary"))
    print(f"  Years:        {data.n_years}")
    print(f"  Nodes:        {data.n_nodes}")
    print(f"  Disease year: {data.disease_year}")
    print(f"  Seed:         {data.seed}")
    print(f"  Initial pop:  {data.initial_total_pop:,}")
    print(f"  Final pop:    {data.final_total_pop:,}")
    print()

    # Per-node one-line summary
    print(f"  {'Node':<30s} {'K':>5s} {'Final':>6s} {'Crash%':>7s} {'DisDth':>7s} {'r̄ fin':>7s} {'EF1A':>6s}")
    print(f"  {'─'*30} {'─'*5} {'─'*6} {'─'*7} {'─'*7} {'─'*7} {'─'*6}")

    for i in range(data.n_nodes):
        name = data.node_names[i]
        K = int(data.node_K[i])
        final_pop = int(data.yearly_pop[i, -1])
        pre_pop = max(int(data.yearly_pop[i, data.disease_year]), 1)
        min_pop = int(data.yearly_pop[i, data.disease_year:].min())
        crash = (1 - min_pop / pre_pop) * 100
        total_dis = int(data.yearly_disease_deaths[i].sum())
        r_final = float(data.yearly_mean_resistance[i, -1])
        ef1a = float(data.yearly_ef1a_freq[i, -1])

        print(f"  {name:<30s} {K:5d} {final_pop:6d} {crash:6.1f}% {total_dis:7d} {r_final:7.4f} {ef1a:6.3f}")


def interrogate_node_year(data: SimulationData, node_id: int, year: int) -> None:
    """Print detailed state for a specific node at a specific year."""
    if node_id < 0 or node_id >= data.n_nodes:
        print(f"Error: node_id {node_id} out of range [0, {data.n_nodes - 1}]")
        sys.exit(1)
    if year < 0 or year >= data.n_years:
        print(f"Error: year {year} out of range [0, {data.n_years - 1}]")
        sys.exit(1)

    name = data.node_names[node_id]
    K = int(data.node_K[node_id])

    print(_divider(f"{name} — Year {year}"))
    phase = "Pre-epidemic" if year < data.disease_year else "Epidemic"
    print(f"  Phase: {phase}")
    print()

    # Population
    pop = int(data.yearly_pop[node_id, year])
    adults = int(data.yearly_adults[node_id, year])
    recruits = int(data.yearly_recruits[node_id, year])
    non_adult = pop - adults
    pop_frac = pop / max(K, 1)

    print(_divider("Population"))
    print(f"  Total population:  {pop:6d}  ({pop_frac:.1%} of K={K})")
    print(f"  Adults:            {adults:6d}  ({adults/max(pop,1):.1%})")
    print(f"  Juv + Subadults:   {non_adult:6d}  ({non_adult/max(pop,1):.1%})")
    print(f"  Recruits (year):   {recruits:6d}")
    print(f"  Pop bar: {_bar(pop, K)}  {pop}/{K}")

    # Demographics
    nat_d = int(data.yearly_natural_deaths[node_id, year])
    dis_d = int(data.yearly_disease_deaths[node_id, year])
    total_d = nat_d + dis_d

    print(_divider("Mortality"))
    print(f"  Natural deaths:    {nat_d:6d}")
    print(f"  Disease deaths:    {dis_d:6d}")
    print(f"  Total deaths:      {total_d:6d}")
    if total_d > 0:
        print(f"  Disease fraction:  {dis_d / total_d:.1%}")

    # Disease
    vibrio = float(data.yearly_vibrio_max[node_id, year])
    peak_prev = float(data.peak_disease_prevalence[node_id])

    print(_divider("Disease State"))
    print(f"  Peak Vibrio:       {vibrio:.0f} bact/mL")
    print(f"  Peak prevalence:   {peak_prev:.1%} (lifetime max)")

    # Genetics
    r_mean = float(data.yearly_mean_resistance[node_id, year])
    ef1a = float(data.yearly_ef1a_freq[node_id, year])
    va = float(data.yearly_va[node_id, year])
    ne_ratio = float(data.yearly_ne_ratio[node_id, year])

    print(_divider("Genetics"))
    print(f"  Mean resistance (r̄): {r_mean:.6f}")
    print(f"  EF1A frequency:      {ef1a:.4f}")
    print(f"  Additive variance:   {va:.6f}")
    print(f"  Ne / N ratio:        {ne_ratio:.4f}")

    # Top-3 allele frequencies
    top3 = data.yearly_allele_freq_top3[node_id, year]  # (3,)
    print(f"\n  Top-3 locus allele frequencies:")
    labels = ["Locus 1 (largest)", "Locus 2", "Locus 3"]
    for j in range(3):
        print(f"    {labels[j]:<20s}: {top3[j]:.4f}  {_bar(top3[j], 1.0, 20)}")

    # Full 52-locus snapshots (if available for this year)
    if year == data.disease_year:
        print(_divider("52-Locus Snapshot (pre-epidemic)"))
        freq = data.pre_epidemic_allele_freq[node_id]
        _print_locus_table(freq)
    elif year == data.disease_year + 2:
        print(_divider("52-Locus Snapshot (post-epidemic)"))
        freq = data.post_epidemic_allele_freq[node_id]
        _print_locus_table(freq)
    elif year > data.disease_year + 2:
        # Show delta
        print(_divider("Δq (post - pre epidemic)"))
        delta = data.post_epidemic_allele_freq[node_id] - data.pre_epidemic_allele_freq[node_id]
        _print_delta_table(delta)

    # Trajectory context: show surrounding years
    print(_divider("Trajectory Context (±3 years)"))
    yr_start = max(0, year - 3)
    yr_end = min(data.n_years, year + 4)
    print(f"  {'Year':>4s} {'Pop':>6s} {'Adults':>7s} {'Recr':>6s} {'NatD':>5s} {'DisD':>5s} {'r̄':>8s} {'EF1A':>6s}")
    for yr in range(yr_start, yr_end):
        marker = " ◄" if yr == year else ""
        print(f"  {yr:4d} {int(data.yearly_pop[node_id, yr]):6d} "
              f"{int(data.yearly_adults[node_id, yr]):7d} "
              f"{int(data.yearly_recruits[node_id, yr]):6d} "
              f"{int(data.yearly_natural_deaths[node_id, yr]):5d} "
              f"{int(data.yearly_disease_deaths[node_id, yr]):5d} "
              f"{float(data.yearly_mean_resistance[node_id, yr]):8.4f} "
              f"{float(data.yearly_ef1a_freq[node_id, yr]):6.3f}{marker}")
    print()


def _print_locus_table(freq: np.ndarray) -> None:
    """Print 52-locus allele frequencies in a compact table."""
    n_loci = len(freq)
    # Sort by frequency descending for readability
    order = np.argsort(freq)[::-1]

    print(f"  {'Rank':>4s} {'Locus':>5s} {'Freq':>7s}  Bar")
    for rank, idx in enumerate(order[:10]):
        locus_name = "EF1A" if idx == 51 else f"L{idx}"
        print(f"  {rank+1:4d} {locus_name:>5s} {freq[idx]:7.4f}  {_bar(freq[idx], 1.0, 25)}")
    print(f"  ... ({n_loci - 10} more loci)")

    # Stats
    print(f"\n  Min freq:  {freq.min():.4f}")
    print(f"  Max freq:  {freq.max():.4f}")
    print(f"  Mean freq: {freq.mean():.4f}")
    print(f"  Std freq:  {freq.std():.4f}")


def _print_delta_table(delta: np.ndarray) -> None:
    """Print allele frequency changes, sorted by magnitude."""
    n_loci = len(delta)
    order = np.argsort(np.abs(delta))[::-1]

    print(f"  {'Rank':>4s} {'Locus':>5s} {'Δq':>8s}  Direction")
    for rank, idx in enumerate(order[:10]):
        locus_name = "EF1A" if idx == 51 else f"L{idx}"
        direction = "↑" if delta[idx] > 0 else "↓" if delta[idx] < 0 else "—"
        print(f"  {rank+1:4d} {locus_name:>5s} {delta[idx]:+8.4f}  {direction}")
    print(f"  ... ({n_loci - 10} more loci)")

    # Stats
    pos = (delta > 0).sum()
    neg = (delta < 0).sum()
    print(f"\n  Increased: {pos} loci, Decreased: {neg} loci")
    print(f"  Mean Δq:   {delta.mean():+.4f}")
    print(f"  Max |Δq|:  {np.abs(delta).max():.4f}")


def dump_individuals(data: SimulationData, node_id: int, year: int) -> None:
    """Dump individual-level data for a node/year.

    Since we only have aggregate yearly data (not per-individual snapshots),
    we reconstruct an approximate individual-level view based on:
    - Population count and stage distribution
    - Mean resistance and Va (to generate representative individuals)
    - EF1A frequency

    For true individual dumps, the simulation would need to save per-agent
    snapshots at each timestep (future feature).
    """
    if node_id < 0 or node_id >= data.n_nodes:
        print(f"Error: node_id {node_id} out of range")
        sys.exit(1)
    if year < 0 or year >= data.n_years:
        print(f"Error: year {year} out of range")
        sys.exit(1)

    name = data.node_names[node_id]
    pop = int(data.yearly_pop[node_id, year])
    adults = int(data.yearly_adults[node_id, year])
    recruits = int(data.yearly_recruits[node_id, year])
    non_adult = pop - adults
    r_mean = float(data.yearly_mean_resistance[node_id, year])
    va = float(data.yearly_va[node_id, year])
    ef1a_q = float(data.yearly_ef1a_freq[node_id, year])

    print(_divider(f"Individual-Level Reconstruction: {name}, Year {year}"))
    print(f"  ⚠ Note: This is a RECONSTRUCTED view from aggregate statistics.")
    print(f"  True individual snapshots require per-timestep agent dumps (future feature).\n")

    if pop == 0:
        print("  Population is EXTINCT — no individuals to display.")
        return

    # Reconstruct representative individuals
    rng = np.random.RandomState(data.seed + node_id * 1000 + year)

    # Generate resistance values from N(r_mean, sqrt(Va))
    r_std = np.sqrt(max(va, 1e-10))
    r_values = np.clip(rng.normal(r_mean, r_std, pop), 0, 1)

    # Assign stages
    stages = []
    recruit_count = min(recruits, pop)
    adult_count = min(adults, pop - recruit_count)
    juvsub_count = pop - recruit_count - adult_count

    stages = (["recruit"] * recruit_count +
              ["juvenile/subadult"] * juvsub_count +
              ["adult"] * adult_count)

    # Generate EF1A genotypes
    ef1a_genotypes = []
    for _ in range(pop):
        a1 = rng.random() < ef1a_q
        a2 = rng.random() < ef1a_q
        if a1 and a2:
            ef1a_genotypes.append("ins/ins")
        elif a1 or a2:
            ef1a_genotypes.append("wt/ins")
        else:
            ef1a_genotypes.append("wt/wt")

    # Print header
    print(f"  {'ID':>5s} {'Stage':<18s} {'r_i':>8s} {'EF1A':>8s} {'Disease':>10s}")
    print(f"  {'─'*5} {'─'*18} {'─'*8} {'─'*8} {'─'*10}")

    # Show first 50 (or all if small population)
    n_show = min(pop, 50)
    for idx in range(n_show):
        # Disease state: randomly assign based on year's disease fraction
        dis_d = data.yearly_disease_deaths[node_id, year]
        dis_prob = min(dis_d / max(pop, 1), 0.5)
        disease = "infected" if rng.random() < dis_prob else "healthy"

        print(f"  {idx:5d} {stages[idx]:<18s} {r_values[idx]:8.4f} {ef1a_genotypes[idx]:>8s} {disease:>10s}")

    if pop > n_show:
        print(f"  ... ({pop - n_show} more individuals)")

    # Summary stats
    print(f"\n  Population: {pop}")
    print(f"  Stages: {recruit_count} recruits, {juvsub_count} juv/sub, {adult_count} adults")
    print(f"  Resistance: mean={r_values.mean():.4f}, std={r_values.std():.4f}, "
          f"min={r_values.min():.4f}, max={r_values.max():.4f}")

    # EF1A genotype counts
    from collections import Counter
    geno_counts = Counter(ef1a_genotypes)
    total = len(ef1a_genotypes)
    print(f"  EF1A genotypes: "
          f"wt/wt={geno_counts.get('wt/wt', 0)} ({geno_counts.get('wt/wt', 0)/total:.1%}), "
          f"wt/ins={geno_counts.get('wt/ins', 0)} ({geno_counts.get('wt/ins', 0)/total:.1%}), "
          f"ins/ins={geno_counts.get('ins/ins', 0)} ({geno_counts.get('ins/ins', 0)/total:.1%})")
    print()


# ═══════════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Interrogate SSWD-EvoEpi simulation results at node/year level.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python -m viz.interrogate results/5node_epidemic_20yr/ --list-nodes
  python -m viz.interrogate results/5node_epidemic_20yr/ --summary
  python -m viz.interrogate results/5node_epidemic_20yr/ --node 2 --year 7
  python -m viz.interrogate results/5node_epidemic_20yr/ --node 2 --year 7 --individuals
        """,
    )
    parser.add_argument("result_dir", help="Directory with simulation_data.npz + metadata.json")
    parser.add_argument("--node", "-n", type=int, default=None,
                        help="Node index (0-based)")
    parser.add_argument("--year", "-y", type=int, default=None,
                        help="Year index")
    parser.add_argument("--list-nodes", action="store_true",
                        help="List available nodes and exit")
    parser.add_argument("--summary", action="store_true",
                        help="Print quick summary and exit")
    parser.add_argument("--individuals", action="store_true",
                        help="Dump reconstructed individual-level data")

    args = parser.parse_args()
    result_dir = Path(args.result_dir)
    if not result_dir.exists():
        print(f"Error: {result_dir} does not exist", file=sys.stderr)
        sys.exit(1)

    data = SimulationData.load(result_dir)

    if args.list_nodes:
        list_nodes(data)
        return

    if args.summary:
        summary_all(data)
        return

    if args.node is None or args.year is None:
        if args.node is None and args.year is None:
            summary_all(data)
            return
        parser.error("Both --node and --year are required for interrogation")

    if args.individuals:
        dump_individuals(data, args.node, args.year)
    else:
        interrogate_node_year(data, args.node, args.year)


if __name__ == "__main__":
    main()
