#!/usr/bin/env python3
"""Analysis 4: Reintroduction Scenarios.

Predicts outcomes of different reintroduction strategies using the full
spatial SSWD-EvoEpi model. Defines a factorial grid of scenarios varying
release size, location, timing, and genetics, then runs each scenario as
an ensemble and produces publication-quality figures.

This analysis uses the release mechanism (model.process_release_event)
which supports three genetics modes: allele_freqs, trait_targets, and
explicit genotypes.

Usage:
    python 04_reintroduction_scenarios.py [--dry-run] [--output-dir DIR]
    python 04_reintroduction_scenarios.py --max-scenarios 10 --seeds 5

See: conservation report Section 9.4
"""

import sys
import json
import time
import warnings
import argparse
from pathlib import Path
from dataclasses import dataclass, field
from itertools import product
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import yaml

# ── Project imports ──────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from sswd_evoepi.config import (
    ReleaseEvent,
    SimulationConfig,
    default_config,
    load_config,
    GeneticsSection,
)
from sswd_evoepi.model import (
    run_spatial_simulation,
    SpatialSimResult,
)
from sswd_evoepi.spatial import (
    MetapopulationNetwork,
    make_5node_network,
)
from sswd_evoepi.types import N_LOCI, Origin, trait_slices


# ═══════════════════════════════════════════════════════════════════════
# SCENARIO DEFINITIONS
# ═══════════════════════════════════════════════════════════════════════

# Release sizes to test (individuals per event)
RELEASE_SIZES = [20, 50, 100, 200, 500]

# Timing: years after SSWD onset (year 0 = disease introduction)
RELEASE_TIMINGS = {
    "immediate": 3,       # Year 3 = immediate response
    "5yr_post": 5,        # Year 5
    "10yr_post": 10,      # Year 10
}

# Genetics presets: (resistance, tolerance, recovery) trait targets
GENETICS_PRESETS = {
    "wild_caught": {
        "mode": "trait_targets",
        "resistance": 0.15,
        "tolerance": 0.10,
        "recovery": 0.02,
        "description": "Wild-type genetics (no selection)",
    },
    "moderate_selection": {
        "mode": "trait_targets",
        "resistance": 0.30,
        "tolerance": 0.10,
        "recovery": 0.02,
        "description": "Moderate resistance from ~3 gen breeding",
    },
    "high_selection": {
        "mode": "allele_freqs",
        "resistance_freq": 0.85,
        "tolerance_freq": 0.12,
        "recovery_freq": 0.03,
        "description": "High resistance from 8+ gen breeding",
    },
}

# 5-node network: name → node_id
NODE_MAP = {
    "Sitka": 0,
    "Howe Sound": 1,
    "SJI": 2,
    "Newport": 3,
    "Monterey": 4,
}


@dataclass
class ScenarioDef:
    """Definition for a single reintroduction scenario."""
    name: str
    release_size: int
    node_ids: List[int]
    node_names: List[str]
    release_year: int          # Year in simulation to release
    genetics_preset: str       # Key into GENETICS_PRESETS
    description: str = ""


def build_scenario_grid() -> List[ScenarioDef]:
    """Build the factorial scenario grid.

    Grid dimensions:
      - Release sizes: 20, 50, 100, 200, 500
      - Locations: all 5 nodes (single-site), plus multi-site
      - Timing: immediate (yr3), 5yr post, 10yr post
      - Genetics: wild-caught, moderate, high selection

    Returns:
        List of ScenarioDef objects.
    """
    scenarios = []
    scenario_id = 0

    # Single-site releases at each node
    for size in RELEASE_SIZES:
        for timing_name, timing_year in RELEASE_TIMINGS.items():
            for gen_name in GENETICS_PRESETS:
                for node_name, node_id in NODE_MAP.items():
                    s = ScenarioDef(
                        name=f"S{scenario_id:04d}_{node_name}_n{size}_{timing_name}_{gen_name}",
                        release_size=size,
                        node_ids=[node_id],
                        node_names=[node_name],
                        release_year=timing_year,
                        genetics_preset=gen_name,
                        description=(
                            f"{size} individuals at {node_name}, "
                            f"{timing_name}, {gen_name}"
                        ),
                    )
                    scenarios.append(s)
                    scenario_id += 1

    # Multi-site releases (north + mid + south)
    multi_site = {
        "spread_3site": (["Sitka", "SJI", "Monterey"], [0, 2, 4]),
        "all_5site": (
            list(NODE_MAP.keys()),
            list(NODE_MAP.values()),
        ),
    }
    for ms_name, (ms_node_names, ms_node_ids) in multi_site.items():
        for size in RELEASE_SIZES:
            for timing_name, timing_year in RELEASE_TIMINGS.items():
                for gen_name in GENETICS_PRESETS:
                    s = ScenarioDef(
                        name=f"S{scenario_id:04d}_{ms_name}_n{size}_{timing_name}_{gen_name}",
                        release_size=size,
                        node_ids=ms_node_ids,
                        node_names=ms_node_names,
                        release_year=timing_year,
                        genetics_preset=gen_name,
                        description=(
                            f"{size}/site at {ms_name} ({len(ms_node_ids)} sites), "
                            f"{timing_name}, {gen_name}"
                        ),
                    )
                    scenarios.append(s)
                    scenario_id += 1

    return scenarios


def scenario_to_release_events(scenario: ScenarioDef) -> List[ReleaseEvent]:
    """Convert a ScenarioDef to a list of ReleaseEvent objects.

    Args:
        scenario: Scenario definition.

    Returns:
        List of ReleaseEvent objects (one per node).
    """
    gen_preset = GENETICS_PRESETS[scenario.genetics_preset]
    release_day = scenario.release_year * 365  # Convert year to day

    events = []
    for node_id in scenario.node_ids:
        # Build genetics kwargs
        mode = gen_preset["mode"]
        if mode == "allele_freqs":
            freqs = np.zeros(N_LOCI)
            freqs[:17] = gen_preset["resistance_freq"]
            freqs[17:34] = gen_preset["tolerance_freq"]
            freqs[34:] = gen_preset["recovery_freq"]
            event = ReleaseEvent(
                time_step=release_day,
                node_id=node_id,
                n_individuals=scenario.release_size,
                genetics_mode="allele_freqs",
                allele_freqs=freqs,
                age_range=(365, 1095),
                mark_released=True,
            )
        else:
            targets = {
                "resistance": gen_preset.get("resistance", 0.15),
                "tolerance": gen_preset.get("tolerance", 0.10),
                "recovery": gen_preset.get("recovery", 0.02),
            }
            event = ReleaseEvent(
                time_step=release_day,
                node_id=node_id,
                n_individuals=scenario.release_size,
                genetics_mode="trait_targets",
                trait_targets=targets,
                age_range=(365, 1095),
                mark_released=True,
            )
        events.append(event)

    return events


# ═══════════════════════════════════════════════════════════════════════
# SIMULATION RUNNER
# ═══════════════════════════════════════════════════════════════════════

def run_scenario(
    scenario: ScenarioDef,
    seed: int,
    config: SimulationConfig,
    n_years: int = 20,
    disease_year: int = 0,
) -> Dict[str, Any]:
    """Run a single scenario replicate.

    Args:
        scenario: Scenario definition.
        seed: Random seed.
        config: Base SimulationConfig (release_events will be overwritten).
        n_years: Years to simulate.
        disease_year: Year to introduce disease.

    Returns:
        Dict with outcome metrics.
    """
    # Convert scenario to release events
    release_events = scenario_to_release_events(scenario)
    config.release_events = release_events

    # Build network
    network = make_5node_network(seed=seed)

    # Run simulation
    result = run_spatial_simulation(
        network=network,
        n_years=n_years,
        disease_year=disease_year,
        seed=seed,
        config=config,
    )

    # Extract metrics
    n_nodes = result.n_nodes
    node_names = result.node_names or [f"Node_{i}" for i in range(n_nodes)]

    metrics = {}
    if result.yearly_pop is not None:
        pop_total = np.sum(result.yearly_pop, axis=0)
        metrics["initial_pop"] = int(pop_total[0])
        metrics["final_pop"] = int(pop_total[-1])
        metrics["min_pop"] = int(np.min(pop_total))
        metrics["min_pop_year"] = int(np.argmin(pop_total))
        metrics["crash_pct"] = float(
            (1.0 - pop_total[-1] / max(pop_total[0], 1)) * 100
        )
        metrics["persistence"] = bool(pop_total[-1] > 0)
        metrics["n_extinct_nodes"] = int(np.sum(result.yearly_pop[:, -1] == 0))
        metrics["yearly_pop_total"] = pop_total.tolist()

        # Per-node final populations
        metrics["final_pop_per_node"] = {
            node_names[i]: int(result.yearly_pop[i, -1])
            for i in range(n_nodes)
        }

    # Trait evolution
    if result.yearly_mean_resistance is not None:
        metrics["final_mean_resistance_per_node"] = {
            node_names[i]: float(result.yearly_mean_resistance[i, -1])
            for i in range(n_nodes)
            if result.yearly_pop[i, -1] > 0
        }
    if result.yearly_mean_tolerance is not None:
        metrics["final_mean_tolerance_per_node"] = {
            node_names[i]: float(result.yearly_mean_tolerance[i, -1])
            for i in range(n_nodes)
            if result.yearly_pop[i, -1] > 0
        }

    # Released individual fate
    metrics["total_released"] = result.total_released
    metrics["released_surviving"] = result.released_surviving
    if result.total_released > 0:
        metrics["released_survival_pct"] = float(
            result.released_surviving / result.total_released * 100
        )
        if result.yearly_released_alive is not None:
            metrics["yearly_released_alive_total"] = (
                np.sum(result.yearly_released_alive, axis=0).tolist()
            )

    return {
        "scenario": scenario.name,
        "seed": seed,
        "metrics": metrics,
    }


# ═══════════════════════════════════════════════════════════════════════
# ANALYSIS & FIGURES
# ═══════════════════════════════════════════════════════════════════════

def generate_figures(
    all_results: List[dict],
    scenarios: List[ScenarioDef],
    output_dir: Path,
) -> List[str]:
    """Generate publication-quality figures from scenario results.

    Figures:
      a. Population recovery curves by release size
      b. Released vs wild-born survival comparison
      c. Genetic introgression: trait values over time
      d. Site × release_size heatmap of recovery probability
      e. Minimum viable release size by site

    Args:
        all_results: List of {scenario_name: {seed: metrics}} results.
        scenarios: List of ScenarioDef objects.
        output_dir: Directory for figure files.

    Returns:
        List of generated figure paths.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import Normalize
        from matplotlib.cm import ScalarMappable
    except ImportError:
        warnings.warn("matplotlib not available — skipping figure generation")
        return []

    fig_dir = output_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    generated = []

    # ── Helper: group results by scenario attributes ─────────────────
    def group_by_attribute(results, scenarios, attr_fn):
        """Group scenario results by an attribute function."""
        groups = {}
        for res, scen in zip(results, scenarios):
            key = attr_fn(scen)
            if key not in groups:
                groups[key] = []
            groups[key].append(res)
        return groups

    # ── Figure (a): Population recovery curves by release size ───────
    try:
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        # Filter: SJI node, 10yr timing, high selection
        for size in RELEASE_SIZES:
            trajectories = []
            for res in all_results:
                scen_name = res.get("scenario", "")
                if (f"_SJI_n{size}_10yr_post_high_selection" in scen_name
                        and "yearly_pop_total" in res.get("metrics", {})):
                    trajectories.append(res["metrics"]["yearly_pop_total"])

            if trajectories:
                arr = np.array(trajectories)
                mean_traj = np.mean(arr, axis=0)
                std_traj = np.std(arr, axis=0)
                years = np.arange(len(mean_traj))
                ax.plot(years, mean_traj, label=f"n={size}", linewidth=2)
                ax.fill_between(
                    years, mean_traj - std_traj, mean_traj + std_traj, alpha=0.2
                )

        ax.set_xlabel("Year")
        ax.set_ylabel("Total population")
        ax.set_title("Population recovery by release size\n(SJI, 10yr post-SSWD, high selection)")
        ax.legend()
        ax.set_yscale("symlog", linthresh=10)
        fig.tight_layout()
        path = fig_dir / "a_recovery_by_release_size.png"
        fig.savefig(path, dpi=150)
        plt.close(fig)
        generated.append(str(path))
    except Exception as e:
        warnings.warn(f"Figure (a) failed: {e}")

    # ── Figure (b): Released vs wild-born survival ───────────────────
    try:
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        sizes_for_plot = []
        released_survival = []
        wild_final_frac = []

        for res in all_results:
            m = res.get("metrics", {})
            scen_name = res.get("scenario", "")
            if ("_SJI_" in scen_name and "_10yr_post_high_selection" in scen_name
                    and "released_survival_pct" in m):
                # Extract release size from scenario name
                for sz in RELEASE_SIZES:
                    if f"_n{sz}_" in scen_name:
                        sizes_for_plot.append(sz)
                        released_survival.append(m["released_survival_pct"])
                        total = m.get("final_pop", 0)
                        rel_alive = m.get("released_surviving", 0)
                        wild = total - rel_alive
                        wild_frac = wild / max(total, 1) * 100
                        wild_final_frac.append(wild_frac)
                        break

        if sizes_for_plot:
            ax.scatter(sizes_for_plot, released_survival, label="Released survival %",
                       s=60, alpha=0.7)
            ax.set_xlabel("Release size")
            ax.set_ylabel("Survival / composition (%)")
            ax.set_title("Released individual survival by release size")
            ax.legend()
            fig.tight_layout()
            path = fig_dir / "b_released_vs_wild_survival.png"
            fig.savefig(path, dpi=150)
            generated.append(str(path))
        plt.close(fig)
    except Exception as e:
        warnings.warn(f"Figure (b) failed: {e}")

    # ── Figure (c): Genetic introgression (trait evolution) ──────────
    try:
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        traits = ["resistance", "tolerance"]
        trait_labels = ["Mean resistance", "Mean tolerance"]

        for ax, trait, label in zip(axes[:2], traits, trait_labels):
            for gen_name in ["wild_caught", "moderate_selection", "high_selection"]:
                trait_key = f"final_mean_{trait}_per_node"
                vals = []
                for res in all_results:
                    m = res.get("metrics", {})
                    scen_name = res.get("scenario", "")
                    if (f"_SJI_n200_10yr_post_{gen_name}" in scen_name
                            and trait_key in m and "SJI" in m.get(trait_key, {})):
                        # Use SJI if available
                        node_key = [k for k in m[trait_key] if "San Juan" in k or "SJI" in k]
                        if not node_key:
                            node_key = list(m[trait_key].keys())
                        if node_key:
                            vals.append(m[trait_key][node_key[0]])

                if vals:
                    desc = GENETICS_PRESETS[gen_name]["description"]
                    ax.bar(desc[:15], np.mean(vals), yerr=np.std(vals) if len(vals) > 1 else 0,
                           capsize=5, alpha=0.7)

            ax.set_ylabel(label)
            ax.set_title(f"Final {trait} at SJI")
            ax.tick_params(axis='x', rotation=30)

        # Third panel: released alive over time
        ax = axes[2]
        for gen_name in ["wild_caught", "high_selection"]:
            for res in all_results:
                scen_name = res.get("scenario", "")
                m = res.get("metrics", {})
                if (f"_SJI_n200_10yr_post_{gen_name}" in scen_name
                        and "yearly_released_alive_total" in m):
                    years = np.arange(len(m["yearly_released_alive_total"]))
                    ax.plot(years, m["yearly_released_alive_total"],
                            label=gen_name.replace("_", " "), alpha=0.7)
                    break  # One line per preset

        ax.set_xlabel("Year")
        ax.set_ylabel("Released individuals alive")
        ax.set_title("Released cohort decay")
        ax.legend(fontsize=8)

        fig.tight_layout()
        path = fig_dir / "c_genetic_introgression.png"
        fig.savefig(path, dpi=150)
        plt.close(fig)
        generated.append(str(path))
    except Exception as e:
        warnings.warn(f"Figure (c) failed: {e}")

    # ── Figure (d): Site × release_size heatmap ──────────────────────
    try:
        fig, ax = plt.subplots(1, 1, figsize=(8, 6))
        node_names_ordered = list(NODE_MAP.keys())
        heatmap = np.full((len(node_names_ordered), len(RELEASE_SIZES)), np.nan)

        for i, node_name in enumerate(node_names_ordered):
            for j, size in enumerate(RELEASE_SIZES):
                persistence_vals = []
                for res in all_results:
                    scen_name = res.get("scenario", "")
                    m = res.get("metrics", {})
                    if (f"_{node_name}_n{size}_10yr_post_high_selection" in scen_name
                            and "persistence" in m):
                        persistence_vals.append(float(m["persistence"]))

                if persistence_vals:
                    heatmap[i, j] = np.mean(persistence_vals)

        im = ax.imshow(heatmap, aspect="auto", cmap="RdYlGn", vmin=0, vmax=1)
        ax.set_xticks(range(len(RELEASE_SIZES)))
        ax.set_xticklabels(RELEASE_SIZES)
        ax.set_yticks(range(len(node_names_ordered)))
        ax.set_yticklabels(node_names_ordered)
        ax.set_xlabel("Release size")
        ax.set_ylabel("Release site")
        ax.set_title("Recovery probability\n(10yr post-SSWD, high selection)")
        fig.colorbar(im, label="P(persistence)")

        # Add text annotations
        for i in range(heatmap.shape[0]):
            for j in range(heatmap.shape[1]):
                if not np.isnan(heatmap[i, j]):
                    ax.text(j, i, f"{heatmap[i, j]:.2f}",
                            ha="center", va="center", fontsize=9,
                            color="white" if heatmap[i, j] < 0.5 else "black")

        fig.tight_layout()
        path = fig_dir / "d_site_size_heatmap.png"
        fig.savefig(path, dpi=150)
        plt.close(fig)
        generated.append(str(path))
    except Exception as e:
        warnings.warn(f"Figure (d) failed: {e}")

    # ── Figure (e): Minimum viable release size by site ──────────────
    try:
        fig, ax = plt.subplots(1, 1, figsize=(8, 5))
        min_viable = {}

        for node_name in NODE_MAP:
            for size in RELEASE_SIZES:
                persistence_vals = []
                for res in all_results:
                    scen_name = res.get("scenario", "")
                    m = res.get("metrics", {})
                    if (f"_{node_name}_n{size}_10yr_post_high_selection" in scen_name
                            and "persistence" in m):
                        persistence_vals.append(float(m["persistence"]))

                if persistence_vals and np.mean(persistence_vals) >= 0.5:
                    if node_name not in min_viable:
                        min_viable[node_name] = size
                    break

        if min_viable:
            sites = list(min_viable.keys())
            sizes_v = [min_viable[s] for s in sites]
            bars = ax.barh(sites, sizes_v, color="steelblue", alpha=0.8)
            ax.set_xlabel("Minimum release size for ≥50% persistence")
            ax.set_title("Minimum viable release size by site\n(high selection, 10yr post-SSWD)")

            for bar, val in zip(bars, sizes_v):
                ax.text(bar.get_width() + 5, bar.get_y() + bar.get_height() / 2,
                        str(val), va="center")

        fig.tight_layout()
        path = fig_dir / "e_minimum_viable_release.png"
        fig.savefig(path, dpi=150)
        plt.close(fig)
        generated.append(str(path))
    except Exception as e:
        warnings.warn(f"Figure (e) failed: {e}")

    return generated


# ═══════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Analysis 4: Reintroduction scenarios"
    )
    parser.add_argument("--dry-run", action="store_true",
                        help="Print scenario grid without running")
    parser.add_argument("--output-dir", type=str, default=None)
    parser.add_argument("--params", type=str, default=None)
    parser.add_argument("--seeds", type=int, default=None,
                        help="Seeds per scenario (default: from params)")
    parser.add_argument("--max-scenarios", type=int, default=None,
                        help="Limit scenarios (for testing)")
    parser.add_argument("--calibrated", type=str, default=None,
                        help="Path to calibrated params JSON (from ABC-SMC)")
    parser.add_argument("--n-years", type=int, default=20,
                        help="Simulation years (default: 20)")
    args = parser.parse_args()

    print("=" * 60)
    print("Analysis 4: Reintroduction Scenarios")
    print("=" * 60)

    # Load analysis params
    params_path = args.params or str(Path(__file__).parent / "params.yaml")
    with open(params_path) as f:
        params = yaml.safe_load(f)

    # Load calibrated params if available
    calibrated_path = args.calibrated or params.get("calibrated_params_file")
    if calibrated_path and Path(calibrated_path).exists():
        print(f"\n  Loading calibrated params: {calibrated_path}")
        with open(calibrated_path) as f:
            calibrated = json.load(f)
    else:
        calibrated = None
        warnings.warn(
            "⚠️  No calibrated parameters found. Running with DEFAULT params.\n"
            "  Results are illustrative only — rerun after ABC-SMC calibration.",
            UserWarning,
            stacklevel=2,
        )

    # Build config
    config = default_config()
    # TODO: Apply calibrated params to config when available
    # if calibrated:
    #     config = apply_calibrated_params(config, calibrated)

    # Build scenario grid
    scenarios = build_scenario_grid()
    n_seeds = args.seeds or params.get("reintroduction", {}).get("n_seeds", 5)

    print(f"\n  Total scenarios: {len(scenarios)}")
    print(f"  Seeds per scenario: {n_seeds}")
    print(f"  Total runs: {len(scenarios) * n_seeds}")

    if args.dry_run:
        print("\n── SCENARIO GRID (dry run) ──")
        # Summarize by category
        by_size = {}
        by_timing = {}
        by_genetics = {}
        for s in scenarios:
            by_size[s.release_size] = by_size.get(s.release_size, 0) + 1
            yr_label = [k for k, v in RELEASE_TIMINGS.items() if v == s.release_year]
            yr_label = yr_label[0] if yr_label else str(s.release_year)
            by_timing[yr_label] = by_timing.get(yr_label, 0) + 1
            by_genetics[s.genetics_preset] = by_genetics.get(s.genetics_preset, 0) + 1

        print(f"\n  By release size: {dict(sorted(by_size.items()))}")
        print(f"  By timing: {by_timing}")
        print(f"  By genetics: {by_genetics}")

        print("\n  First 20 scenarios:")
        for s in scenarios[:20]:
            print(f"    {s.name}: {s.description}")
        if len(scenarios) > 20:
            print(f"    ... and {len(scenarios) - 20} more")

        est_hours = len(scenarios) * n_seeds * 0.04  # ~2.4 min per run
        print(f"\n  Estimated compute time (single core):")
        print(f"    {est_hours:.0f} core-hours")
        print(f"    (~{est_hours/48:.1f} wall-hours on Xeon 48-core)")
        return

    # Output directory
    output_dir = Path(args.output_dir) if args.output_dir else \
        Path(__file__).parent / ".." / params.get("output", {}).get("analysis4", "../results/04_reintroduction")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Limit scenarios for testing
    if args.max_scenarios:
        scenarios = scenarios[:args.max_scenarios]
        print(f"  ⚠️  Limited to {args.max_scenarios} scenarios")

    # Run all scenarios
    all_results = []
    total_runs = len(scenarios) * n_seeds
    run_count = 0
    t_start = time.time()

    for i, scenario in enumerate(scenarios):
        scenario_results = []

        for seed_idx in range(n_seeds):
            seed = 42 + seed_idx * 1000
            run_count += 1

            elapsed = time.time() - t_start
            rate = run_count / max(elapsed, 1)
            eta = (total_runs - run_count) / max(rate, 0.001)

            print(f"  [{run_count}/{total_runs}] {scenario.name} seed={seed} "
                  f"(ETA {eta/60:.0f}min)...", end="", flush=True)

            t0 = time.time()
            try:
                result = run_scenario(
                    scenario, seed, config,
                    n_years=args.n_years,
                    disease_year=0,
                )
                dt = time.time() - t0
                m = result.get("metrics", {})
                print(f" {dt:.1f}s — pop={m.get('final_pop', '?')}")
                scenario_results.append(result)
                all_results.append(result)
            except Exception as e:
                dt = time.time() - t0
                print(f" FAILED ({dt:.1f}s): {e}")
                continue

        # Save per-scenario ensemble
        if scenario_results:
            scen_dir = output_dir / scenario.name
            scen_dir.mkdir(parents=True, exist_ok=True)
            with open(scen_dir / "results.json", "w") as f:
                json.dump(scenario_results, f, indent=2, default=str)

    # Save all results
    with open(output_dir / "all_results.json", "w") as f:
        json.dump(all_results, f, indent=2, default=str)

    total_time = time.time() - t_start
    print(f"\n  Completed {len(all_results)}/{total_runs} runs in {total_time/60:.1f} min")

    # Generate figures
    print("\n  Generating figures...")
    fig_paths = generate_figures(all_results, scenarios, output_dir)
    for fp in fig_paths:
        print(f"    ✅ {fp}")

    # Save master summary
    summary = {
        "n_scenarios": len(scenarios),
        "n_seeds": n_seeds,
        "total_runs": len(all_results),
        "total_time_min": total_time / 60,
        "calibrated": calibrated_path is not None,
        "figures": fig_paths,
        "todos": [
            "Calibrate model parameters (ABC-SMC) before production runs",
            "Extend to 11-node network after calibration",
            "Run full grid on Xeon (48 cores)",
        ] if not calibrated else [],
    }
    with open(output_dir / "analysis4_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\n✅ Analysis 4 complete.")
    print(f"   Output: {output_dir}")
    if not calibrated:
        print("   ⚠️  Results use default params — rerun after calibration.")


if __name__ == "__main__":
    main()
