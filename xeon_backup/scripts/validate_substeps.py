#!/usr/bin/env python3
"""Validate substeps=6 vs substeps=24: multi-seed distributional comparison.

Single-seed comparisons are meaningless because changing substeps changes
the number of RNG draws per day, shifting the entire random sequence.
Instead, we run many seeds and compare the *distributions* of outcomes.

Key question: Do substeps=6 and substeps=24 produce statistically
equivalent distributions of population-level outcomes?

Output:
  - Console: progress + summary statistics
  - audit/substep_validation.md: formal results report
"""

import sys
import os
import time
import json
from pathlib import Path
from dataclasses import dataclass
from typing import List, Dict, Tuple

import numpy as np
from scipy import stats

# Add project root to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from sswd_evoepi.config import SimulationConfig, default_config
from sswd_evoepi.spatial import (
    make_5node_network,
    get_5node_definitions,
    build_network,
    NodeDefinition,
)
from sswd_evoepi.model import run_spatial_simulation


# ═══════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════

N_SEEDS = 10           # Seeds 0-9
SUBSTEP_VALUES = [6, 24]
N_YEARS = 5
DISEASE_YEAR = 1
K_SCALE = 2000 / 3500  # Scale from default total K (~3500) to ~2000


# ═══════════════════════════════════════════════════════════════════════
# HELPERS
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class RunResult:
    """Metrics from a single simulation run."""
    seed: int
    substeps: int
    final_total_pop: int
    final_pop_per_node: List[int]
    total_disease_deaths: int
    crash_pct: float  # (1 - min_pop/initial_pop) * 100
    initial_total_pop: int
    min_pop: int
    runtime_sec: float


def make_scaled_network(seed: int) -> Tuple:
    """Create a 5-node network with K scaled down to ~2000 total."""
    node_defs_orig = get_5node_definitions()
    node_defs = []
    for nd in node_defs_orig:
        new_K = max(50, int(nd.carrying_capacity * K_SCALE))
        # Scale habitat area proportionally
        new_area = nd.habitat_area * K_SCALE
        node_defs.append(NodeDefinition(
            node_id=nd.node_id,
            name=nd.name,
            lat=nd.lat,
            lon=nd.lon,
            subregion=nd.subregion,
            habitat_area=new_area,
            carrying_capacity=new_K,
            is_fjord=nd.is_fjord,
            sill_depth=nd.sill_depth,
            flushing_rate=nd.flushing_rate,
            mean_sst=nd.mean_sst,
            sst_amplitude=nd.sst_amplitude,
            sst_trend=nd.sst_trend,
            salinity=nd.salinity,
            depth_range=nd.depth_range,
        ))
    total_K = sum(nd.carrying_capacity for nd in node_defs)
    network = build_network(node_defs, seed=seed)
    return network, total_K, node_defs


def run_single(seed: int, substeps: int) -> RunResult:
    """Run one simulation with given seed and substep count."""
    config = default_config()
    config.movement.substeps_per_day = substeps

    # Build network with this seed (network construction uses the seed)
    network, total_K, node_defs = make_scaled_network(seed=seed)

    t0 = time.time()
    result = run_spatial_simulation(
        network=network,
        n_years=N_YEARS,
        disease_year=DISEASE_YEAR,
        initial_infected_per_node=3,
        seed=seed,
        config=config,
    )
    runtime = time.time() - t0

    # Extract per-node final populations
    final_pop_per_node = []
    for i in range(result.n_nodes):
        final_pop_per_node.append(int(result.yearly_pop[i, -1]))

    final_total = int(result.yearly_total_pop[-1])
    initial_total = result.initial_total_pop

    # Crash percentage: how much did population drop at minimum?
    min_pop = int(np.min(result.yearly_total_pop))
    crash_pct = (1.0 - min_pop / initial_total) * 100.0 if initial_total > 0 else 0.0

    total_dd = int(np.sum(result.yearly_disease_deaths)) if result.yearly_disease_deaths is not None else 0

    return RunResult(
        seed=seed,
        substeps=substeps,
        final_total_pop=final_total,
        final_pop_per_node=final_pop_per_node,
        total_disease_deaths=total_dd,
        crash_pct=crash_pct,
        initial_total_pop=initial_total,
        min_pop=min_pop,
        runtime_sec=runtime,
    )


def cohen_d(group1: np.ndarray, group2: np.ndarray) -> float:
    """Compute Cohen's d effect size."""
    n1, n2 = len(group1), len(group2)
    var1, var2 = np.var(group1, ddof=1), np.var(group2, ddof=1)
    pooled_std = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
    if pooled_std == 0:
        return 0.0
    return (np.mean(group1) - np.mean(group2)) / pooled_std


# ═══════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════

def main():
    print("=" * 70)
    print("SUBSTEP VALIDATION: substeps=6 vs substeps=24")
    print("Multi-seed distributional comparison (10 seeds × 2 settings)")
    print("=" * 70)
    print(f"\nConfig: {N_SEEDS} seeds, {N_YEARS} years, disease_year={DISEASE_YEAR}")
    print(f"Network: 5 nodes, K scaled to ~{int(3500 * K_SCALE)} total")
    print(f"Expected runtime: ~20-40 minutes\n")

    results: Dict[int, List[RunResult]] = {s: [] for s in SUBSTEP_VALUES}
    total_runs = N_SEEDS * len(SUBSTEP_VALUES)
    run_count = 0

    for substeps in SUBSTEP_VALUES:
        print(f"\n{'─' * 50}")
        print(f"Running substeps={substeps} ({N_SEEDS} seeds)")
        print(f"{'─' * 50}")

        for seed in range(N_SEEDS):
            run_count += 1
            print(f"  [{run_count}/{total_runs}] seed={seed}, substeps={substeps} ...", end=" ", flush=True)
            try:
                r = run_single(seed, substeps)
                results[substeps].append(r)
                print(f"done ({r.runtime_sec:.1f}s) | pop={r.final_total_pop}, "
                      f"crash={r.crash_pct:.1f}%, deaths={r.total_disease_deaths}")
            except Exception as e:
                print(f"FAILED: {e}")
                import traceback
                traceback.print_exc()

    # ── Statistical Analysis ─────────────────────────────────────────
    print("\n" + "=" * 70)
    print("STATISTICAL ANALYSIS")
    print("=" * 70)

    metrics = {}
    for substeps in SUBSTEP_VALUES:
        runs = results[substeps]
        if not runs:
            print(f"  No results for substeps={substeps}!")
            continue
        metrics[substeps] = {
            'final_pop': np.array([r.final_total_pop for r in runs]),
            'crash_pct': np.array([r.crash_pct for r in runs]),
            'disease_deaths': np.array([r.total_disease_deaths for r in runs]),
            'initial_pop': np.array([r.initial_total_pop for r in runs]),
            'min_pop': np.array([r.min_pop for r in runs]),
            'runtime': np.array([r.runtime_sec for r in runs]),
        }

    if len(metrics) < 2:
        print("ERROR: Not enough data for comparison.")
        return

    analysis_lines = []  # for the markdown report

    for metric_name in ['final_pop', 'crash_pct', 'disease_deaths', 'min_pop']:
        g6 = metrics[6][metric_name]
        g24 = metrics[24][metric_name]

        mean6, std6 = np.mean(g6), np.std(g6, ddof=1)
        mean24, std24 = np.mean(g24), np.std(g24, ddof=1)

        # Welch's t-test
        t_stat, p_welch = stats.ttest_ind(g6, g24, equal_var=False)

        # Mann-Whitney U test (non-parametric)
        u_stat, p_mw = stats.mannwhitneyu(g6, g24, alternative='two-sided')

        # Effect size (Cohen's d)
        d = cohen_d(g6, g24)

        print(f"\n  {metric_name}:")
        print(f"    substeps=6:  mean={mean6:.1f}, std={std6:.1f}, range=[{g6.min():.1f}, {g6.max():.1f}]")
        print(f"    substeps=24: mean={mean24:.1f}, std={std24:.1f}, range=[{g24.min():.1f}, {g24.max():.1f}]")
        print(f"    Welch's t-test: t={t_stat:.3f}, p={p_welch:.4f}")
        print(f"    Mann-Whitney U: U={u_stat:.1f}, p={p_mw:.4f}")
        print(f"    Cohen's d: {d:.3f} ({'negligible' if abs(d) < 0.2 else 'small' if abs(d) < 0.5 else 'medium' if abs(d) < 0.8 else 'large'})")

        analysis_lines.append({
            'metric': metric_name,
            'mean_6': mean6, 'std_6': std6,
            'mean_24': mean24, 'std_24': std24,
            'values_6': g6.tolist(), 'values_24': g24.tolist(),
            't_stat': t_stat, 'p_welch': p_welch,
            'u_stat': u_stat, 'p_mw': p_mw,
            'cohens_d': d,
        })

    # Runtime comparison
    rt6 = metrics[6]['runtime']
    rt24 = metrics[24]['runtime']
    speedup = np.mean(rt24) / np.mean(rt6) if np.mean(rt6) > 0 else 0
    print(f"\n  Runtime:")
    print(f"    substeps=6:  mean={np.mean(rt6):.1f}s")
    print(f"    substeps=24: mean={np.mean(rt24):.1f}s")
    print(f"    Speedup: {speedup:.2f}×")

    # ── Summary Verdict ──────────────────────────────────────────────
    any_significant = any(a['p_welch'] < 0.05 for a in analysis_lines)
    any_large_effect = any(abs(a['cohens_d']) >= 0.8 for a in analysis_lines)

    print("\n" + "=" * 70)
    print("VERDICT")
    print("=" * 70)
    if not any_significant and not any_large_effect:
        verdict = "YES — substeps=6 is statistically equivalent to substeps=24"
        print(f"  ✅ {verdict}")
        print("  No significant differences (p > 0.05) and no large effect sizes.")
        print(f"  Recommended: Use substeps=6 for ~{speedup:.1f}× speedup.")
    elif any_significant and any_large_effect:
        verdict = "NO — substeps=6 produces significantly different dynamics"
        print(f"  ❌ {verdict}")
        print("  Significant differences detected with large effect sizes.")
        print("  Recommended: Keep substeps=24.")
    else:
        verdict = "MARGINAL — some differences detected but effect sizes are small"
        print(f"  ⚠️  {verdict}")
        sig_metrics = [a['metric'] for a in analysis_lines if a['p_welch'] < 0.05]
        print(f"  Significant metrics: {sig_metrics}")
        print(f"  Consider substeps=6 if speedup ({speedup:.1f}×) is worth the tradeoff.")

    # ── Write audit report ───────────────────────────────────────────
    audit_dir = Path(__file__).resolve().parent.parent / "audit"
    audit_dir.mkdir(exist_ok=True)
    report_path = audit_dir / "substep_validation.md"

    with open(report_path, "w") as f:
        f.write("# Substep Validation: substeps=6 vs substeps=24\n\n")
        f.write(f"**Date:** {time.strftime('%Y-%m-%d %H:%M')}\n\n")
        f.write("## Methodology\n\n")
        f.write("Single-seed comparisons are meaningless because changing substeps\n")
        f.write("changes the number of RNG draws per day, shifting the entire random\n")
        f.write("sequence. Instead, we compare **distributions** across many seeds.\n\n")
        f.write(f"- **Seeds:** 0–{N_SEEDS-1} ({N_SEEDS} per setting)\n")
        f.write(f"- **Network:** 5 nodes, K scaled to ~{int(3500 * K_SCALE)} total\n")
        f.write(f"- **Duration:** {N_YEARS} years, disease introduced at year {DISEASE_YEAR}\n")
        f.write(f"- **Total runs:** {total_runs}\n\n")

        f.write("## Results Summary\n\n")
        f.write("| Metric | substeps=6 (mean±std) | substeps=24 (mean±std) | Welch p | M-W p | Cohen's d | Significant? |\n")
        f.write("|--------|-----------------------|------------------------|---------|-------|-----------|-------------|\n")
        for a in analysis_lines:
            sig = "Yes" if a['p_welch'] < 0.05 else "No"
            d_label = f"{a['cohens_d']:.3f}"
            f.write(f"| {a['metric']} | {a['mean_6']:.1f} ± {a['std_6']:.1f} | "
                    f"{a['mean_24']:.1f} ± {a['std_24']:.1f} | "
                    f"{a['p_welch']:.4f} | {a['p_mw']:.4f} | {d_label} | {sig} |\n")

        f.write(f"\n## Runtime\n\n")
        f.write(f"- **substeps=6:** {np.mean(rt6):.1f}s mean\n")
        f.write(f"- **substeps=24:** {np.mean(rt24):.1f}s mean\n")
        f.write(f"- **Speedup:** {speedup:.2f}×\n\n")

        f.write("## Raw Data\n\n")
        for substeps in SUBSTEP_VALUES:
            f.write(f"### substeps={substeps}\n\n")
            f.write("| Seed | Final Pop | Crash % | Disease Deaths | Min Pop | Runtime (s) |\n")
            f.write("|------|-----------|---------|----------------|---------|-------------|\n")
            for r in results[substeps]:
                f.write(f"| {r.seed} | {r.final_total_pop} | {r.crash_pct:.1f}% | "
                        f"{r.total_disease_deaths} | {r.min_pop} | {r.runtime_sec:.1f} |\n")
            f.write("\n")

        f.write("## Per-Seed Comparison (node-level final populations)\n\n")
        for seed in range(N_SEEDS):
            r6 = [r for r in results[6] if r.seed == seed]
            r24 = [r for r in results[24] if r.seed == seed]
            if r6 and r24:
                f.write(f"**Seed {seed}:**\n")
                f.write(f"- substeps=6:  {r6[0].final_pop_per_node} (total={r6[0].final_total_pop})\n")
                f.write(f"- substeps=24: {r24[0].final_pop_per_node} (total={r24[0].final_total_pop})\n\n")

        f.write("## Verdict\n\n")
        f.write(f"**{verdict}**\n\n")
        if not any_significant and not any_large_effect:
            f.write("No statistically significant differences were found between substeps=6 and\n")
            f.write("substeps=24 across any population-level metric (all Welch p > 0.05, all\n")
            f.write(f"Cohen's d < 0.8). The {speedup:.1f}× runtime speedup makes substeps=6\n")
            f.write("a safe optimization.\n")
        elif any_significant and any_large_effect:
            f.write("Statistically significant differences with large effect sizes were found.\n")
            f.write("substeps=24 should be retained for accurate dynamics.\n")
        else:
            f.write("Some statistical differences were detected, but effect sizes are small.\n")
            f.write(f"The {speedup:.1f}× speedup may justify using substeps=6 depending on\n")
            f.write("the acceptable tolerance for population-level accuracy.\n")

    print(f"\n📄 Report written to: {report_path}")


if __name__ == "__main__":
    main()
