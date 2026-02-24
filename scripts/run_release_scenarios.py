#!/usr/bin/env python3
"""Run reintroduction scenarios from YAML configuration files.

Parses YAML scenario files, resolves node names to IDs, converts
genetics config to ReleaseEvent objects, runs the spatial model with
specified releases, and saves per-replicate results as JSON.

Usage:
    python scripts/run_release_scenarios.py scenarios/small_pilot.yaml
    python scripts/run_release_scenarios.py scenarios/*.yaml --cores 8
    python scripts/run_release_scenarios.py scenarios/large_program.yaml --dry-run

See scenarios/ directory for example YAML configurations.

References:
    - sswd_evoepi/config.py: ReleaseEvent, SimulationConfig
    - sswd_evoepi/model.py: run_spatial_simulation, process_release_event
    - sswd_evoepi/spatial.py: MetapopulationNetwork
"""

import argparse
import json
import sys
import time
import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import yaml

# ── Project imports ──────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from sswd_evoepi.config import (
    ReleaseEvent,
    SimulationConfig,
    default_config,
    load_config,
)
from sswd_evoepi.model import (
    run_coupled_simulation,
    run_spatial_simulation,
    CoupledSimResult,
    SpatialSimResult,
)
from sswd_evoepi.spatial import (
    MetapopulationNetwork,
    NodeDefinition,
    SpatialNode,
    make_5node_network,
)
from sswd_evoepi.types import N_LOCI, Origin, trait_slices


# ═══════════════════════════════════════════════════════════════════════
# NODE NAME RESOLUTION
# ═══════════════════════════════════════════════════════════════════════

# Canonical 11-node stepping-stone network (from SA R4)
# Maps short names → (node_id, lat, lon, mean_sst, sst_amplitude,
#                      salinity, is_fjord, flushing_rate)
NODE_CATALOG = {
    "Sitka":         (0,  57.06, -135.34, 8.0,  3.5, 32.0, False, 0.80),
    "Ketchikan":     (1,  55.34, -131.64, 8.5,  3.5, 31.0, False, 0.50),
    "Haida Gwaii":   (2,  53.25, -132.07, 9.0,  3.0, 31.5, False, 0.60),
    "Bella Bella":   (3,  52.16, -128.15, 9.5,  3.5, 28.0, False, 0.40),
    "Howe Sound":    (4,  49.52, -123.25, 10.0, 4.0, 22.0, True,  0.03),
    "SJI":           (5,  48.53, -123.02, 10.5, 4.0, 30.0, False, 0.30),
    "Westport":      (6,  46.89, -124.10, 11.0, 3.5, 32.0, False, 0.50),
    "Newport":       (7,  44.63, -124.05, 11.5, 3.0, 33.0, False, 0.60),
    "Crescent City": (8,  41.76, -124.20, 12.0, 2.5, 33.0, False, 0.50),
    "Fort Bragg":    (9,  39.45, -123.80, 12.5, 2.5, 33.5, False, 0.50),
    "Monterey":      (10, 36.62, -121.90, 13.0, 2.5, 33.5, False, 0.40),
}

# Aliases for convenience
NODE_ALIASES = {
    "San Juan Islands": "SJI",
    "San Juan": "SJI",
    "FHL": "SJI",
    "Friday Harbor": "SJI",
}

# 5-node canonical network mapping (for smaller runs)
FIVE_NODE_MAP = {
    "Sitka": 0,
    "Howe Sound": 1,
    "SJI": 2,
    "San Juan Islands": 2,
    "Newport": 3,
    "Monterey": 4,
}


def resolve_node_name(name: str, network_size: int = 11) -> int:
    """Resolve a node name to its integer ID.

    Args:
        name: Node name (e.g. "SJI", "Monterey", "San Juan Islands").
        network_size: Number of nodes in the network (5 or 11).

    Returns:
        Integer node ID.

    Raises:
        ValueError: If name cannot be resolved.
    """
    # Check aliases first
    canonical = NODE_ALIASES.get(name, name)

    if network_size <= 5:
        if canonical in FIVE_NODE_MAP:
            return FIVE_NODE_MAP[canonical]
        raise ValueError(
            f"Unknown node '{name}' for 5-node network. "
            f"Valid: {list(FIVE_NODE_MAP.keys())}"
        )

    if canonical in NODE_CATALOG:
        return NODE_CATALOG[canonical][0]

    # Try case-insensitive match
    for cat_name, spec in NODE_CATALOG.items():
        if cat_name.lower() == canonical.lower():
            return spec[0]

    raise ValueError(
        f"Unknown node '{name}'. Valid: {list(NODE_CATALOG.keys())} "
        f"+ aliases {list(NODE_ALIASES.keys())}"
    )


# ═══════════════════════════════════════════════════════════════════════
# SCENARIO PARSING
# ═══════════════════════════════════════════════════════════════════════

def parse_genetics_config(gen_cfg: dict) -> dict:
    """Convert YAML genetics config to ReleaseEvent kwargs.

    Supports three modes:
      - allele_freqs: Per-locus protective allele frequencies
      - trait_targets: Target mean trait values
      - explicit: Direct genotype arrays (not from YAML)

    Args:
        gen_cfg: Genetics section from YAML.

    Returns:
        Dict with genetics_mode + appropriate arrays/dicts.
    """
    mode = gen_cfg.get("mode", "trait_targets")

    if mode == "allele_freqs":
        r_freq = gen_cfg.get("resistance_freq", 0.15)
        t_freq = gen_cfg.get("tolerance_freq", 0.10)
        c_freq = gen_cfg.get("recovery_freq", 0.02)

        # Build per-locus frequency array: 17R + 17T + 17C = 51 loci
        freqs = np.zeros(N_LOCI)
        freqs[:17] = r_freq
        freqs[17:34] = t_freq
        freqs[34:] = c_freq

        return {
            "genetics_mode": "allele_freqs",
            "allele_freqs": freqs,
        }

    elif mode == "trait_targets":
        targets = {
            "resistance": gen_cfg.get("resistance", 0.15),
            "tolerance": gen_cfg.get("tolerance", 0.10),
            "recovery": gen_cfg.get("recovery", 0.02),
        }
        return {
            "genetics_mode": "trait_targets",
            "trait_targets": targets,
        }

    else:
        raise ValueError(f"Unknown genetics mode: {mode}. Use 'allele_freqs' or 'trait_targets'.")


def parse_scenario(yaml_path: str, network_size: int = 5) -> Tuple[dict, List[ReleaseEvent], dict]:
    """Parse a scenario YAML file into config, release events, and output spec.

    Args:
        yaml_path: Path to YAML scenario file.
        network_size: Network size (5 or 11) for node resolution.

    Returns:
        (scenario_meta, release_events, output_config)
    """
    with open(yaml_path) as f:
        spec = yaml.safe_load(f)

    scenario_meta = {
        "name": spec.get("name", Path(yaml_path).stem),
        "description": spec.get("description", ""),
        "base_config": spec.get("base_config", {}),
        "source_file": str(yaml_path),
    }

    # Parse releases
    release_events = []
    for i, rel in enumerate(spec.get("releases", [])):
        # Resolve node name to ID
        node_name = rel.get("node", "SJI")
        node_id = resolve_node_name(node_name, network_size=network_size)

        # Parse genetics
        gen_cfg = rel.get("genetics", {"mode": "trait_targets"})
        gen_kwargs = parse_genetics_config(gen_cfg)

        # Parse age range
        age_range = tuple(rel.get("age_range", [365, 730]))

        event = ReleaseEvent(
            time_step=rel["time_step"],
            node_id=node_id,
            n_individuals=rel.get("n_individuals", 100),
            age_range=age_range,
            mark_released=True,
            **gen_kwargs,
        )
        release_events.append(event)

    output_config = spec.get("output", {})
    if "dir" not in output_config:
        output_config["dir"] = f"results/release_scenarios/{Path(yaml_path).stem}/"

    return scenario_meta, release_events, output_config


# ═══════════════════════════════════════════════════════════════════════
# NETWORK BUILDING
# ═══════════════════════════════════════════════════════════════════════

def build_scenario_network(
    base_cfg: dict,
    seed: int = 42,
    k_per_node: int = 5000,
) -> MetapopulationNetwork:
    """Build the simulation network for a scenario.

    Currently uses the 5-node canonical network. Will be extended to
    support 11-node and custom networks after calibration.

    Args:
        base_cfg: Base configuration dict from scenario YAML.
        seed: Random seed for network construction.
        k_per_node: Carrying capacity per node.

    Returns:
        MetapopulationNetwork ready for simulation.
    """
    # For now, always use 5-node canonical network
    # TODO: Support 11-node and custom networks via base_cfg
    return make_5node_network(seed=seed)


# ═══════════════════════════════════════════════════════════════════════
# RESULT EXTRACTION
# ═══════════════════════════════════════════════════════════════════════

def extract_replicate_results(
    result: SpatialSimResult,
    scenario_meta: dict,
    seed: int,
    track_released: bool = True,
) -> dict:
    """Extract structured results from a spatial simulation run.

    Args:
        result: SpatialSimResult from run_spatial_simulation.
        scenario_meta: Scenario metadata dict.
        seed: Random seed used for this replicate.
        track_released: Whether to extract released individual tracking.

    Returns:
        Dict with all result metrics in JSON-serializable format.
    """
    n_nodes = result.n_nodes
    n_years = result.yearly_pop.shape[1] if result.yearly_pop is not None else 0
    node_names = result.node_names or [f"Node_{i}" for i in range(n_nodes)]

    out = {
        "scenario": scenario_meta["name"],
        "seed": seed,
        "n_nodes": n_nodes,
        "n_years": n_years,
        "node_names": node_names,
    }

    # Population trajectories
    if result.yearly_pop is not None:
        out["yearly_pop_per_node"] = {
            node_names[i]: result.yearly_pop[i].tolist()
            for i in range(n_nodes)
        }
        out["yearly_pop_total"] = result.yearly_total_pop.tolist() if result.yearly_total_pop is not None else None

    # Disease metrics
    if result.yearly_infected is not None:
        out["yearly_infected_per_node"] = {
            node_names[i]: result.yearly_infected[i].tolist()
            for i in range(n_nodes)
        }
    if result.yearly_disease_deaths is not None:
        out["yearly_disease_deaths_per_node"] = {
            node_names[i]: result.yearly_disease_deaths[i].tolist()
            for i in range(n_nodes)
        }

    # Genetics / trait evolution
    if result.yearly_mean_resistance is not None:
        out["yearly_mean_resistance_per_node"] = {
            node_names[i]: result.yearly_mean_resistance[i].tolist()
            for i in range(n_nodes)
        }
    if result.yearly_mean_tolerance is not None:
        out["yearly_mean_tolerance_per_node"] = {
            node_names[i]: result.yearly_mean_tolerance[i].tolist()
            for i in range(n_nodes)
        }
    if result.yearly_mean_recovery is not None:
        out["yearly_mean_recovery_per_node"] = {
            node_names[i]: result.yearly_mean_recovery[i].tolist()
            for i in range(n_nodes)
        }

    # Release tracking
    if track_released and result.release_log:
        out["release_log"] = result.release_log
        out["total_released"] = result.total_released
        out["released_surviving"] = result.released_surviving
        if result.yearly_released_alive is not None:
            out["yearly_released_alive_per_node"] = {
                node_names[i]: result.yearly_released_alive[i].tolist()
                for i in range(n_nodes)
            }

    # Summary metrics
    if result.yearly_pop is not None and n_years > 0:
        initial_total = int(np.sum(result.yearly_pop[:, 0]))
        final_total = int(np.sum(result.yearly_pop[:, -1]))
        pop_trajectory = np.sum(result.yearly_pop, axis=0)
        min_pop = int(np.min(pop_trajectory))
        min_pop_year = int(np.argmin(pop_trajectory))

        out["metrics"] = {
            "initial_total_pop": initial_total,
            "final_total_pop": final_total,
            "crash_pct": float(1.0 - final_total / max(initial_total, 1)) * 100,
            "min_pop": min_pop,
            "min_pop_year": min_pop_year,
            "n_extinct_nodes": int(np.sum(result.yearly_pop[:, -1] == 0)),
            "persistence": final_total > 0,
        }

        # Released individual metrics
        if track_released and result.total_released > 0:
            out["metrics"]["total_released"] = result.total_released
            out["metrics"]["released_surviving"] = result.released_surviving
            out["metrics"]["released_survival_pct"] = (
                float(result.released_surviving / result.total_released * 100)
            )
    else:
        out["metrics"] = {}

    return out


def compute_ensemble_summary(replicate_results: List[dict]) -> dict:
    """Compute ensemble summary statistics across replicates.

    Args:
        replicate_results: List of per-replicate result dicts.

    Returns:
        Dict with mean/std/min/max for key metrics.
    """
    if not replicate_results:
        return {}

    metrics_keys = [
        "final_total_pop", "crash_pct", "min_pop", "n_extinct_nodes",
        "released_surviving", "released_survival_pct",
    ]

    summary = {
        "n_replicates": len(replicate_results),
        "scenario": replicate_results[0].get("scenario", "unknown"),
    }

    for key in metrics_keys:
        vals = [
            r["metrics"][key]
            for r in replicate_results
            if key in r.get("metrics", {}) and r["metrics"][key] is not None
        ]
        if vals:
            arr = np.array(vals, dtype=float)
            summary[key] = {
                "mean": float(np.mean(arr)),
                "std": float(np.std(arr)),
                "min": float(np.min(arr)),
                "max": float(np.max(arr)),
                "median": float(np.median(arr)),
            }

    # Persistence probability
    persistence_vals = [
        r["metrics"].get("persistence", False)
        for r in replicate_results
        if "metrics" in r
    ]
    if persistence_vals:
        summary["persistence_probability"] = float(np.mean(persistence_vals))

    return summary


# ═══════════════════════════════════════════════════════════════════════
# MAIN RUNNER
# ═══════════════════════════════════════════════════════════════════════

def run_scenario_file(
    yaml_path: str,
    output_dir: Optional[str] = None,
    dry_run: bool = False,
    network_size: int = 5,
    base_config_path: Optional[str] = None,
    verbose: bool = True,
) -> Optional[dict]:
    """Run all replicates for a single scenario YAML file.

    Args:
        yaml_path: Path to scenario YAML.
        output_dir: Override output directory (default: from YAML).
        dry_run: If True, parse and print config without running.
        network_size: Network size (5 or 11) for node resolution.
        base_config_path: Optional path to base config YAML for model params.
        verbose: Print progress messages.

    Returns:
        Ensemble summary dict, or None if dry_run.
    """
    scenario_meta, release_events, output_config = parse_scenario(
        yaml_path, network_size=network_size
    )

    out_dir = Path(output_dir or output_config.get("dir", "results/release_scenarios/"))
    out_dir.mkdir(parents=True, exist_ok=True)

    seeds = output_config.get("seeds", [42])
    n_replicates = output_config.get("n_replicates", len(seeds))
    # Ensure we have enough seeds
    while len(seeds) < n_replicates:
        seeds.append(seeds[-1] + 1000)
    seeds = seeds[:n_replicates]

    track_released = output_config.get("track_released", True)
    base_cfg = scenario_meta.get("base_config", {})

    if verbose:
        print(f"\n{'=' * 60}")
        print(f"Scenario: {scenario_meta['name']}")
        print(f"  Source: {yaml_path}")
        print(f"  Releases: {len(release_events)}")
        for i, ev in enumerate(release_events):
            mode = ev.genetics_mode
            print(f"    [{i+1}] day {ev.time_step}, node {ev.node_id}, "
                  f"n={ev.n_individuals}, mode={mode}")
        print(f"  Replicates: {n_replicates}")
        print(f"  Seeds: {seeds}")
        print(f"  Output: {out_dir}")

    if dry_run:
        print("\n  [DRY RUN] — would run simulation here")
        return None

    # Build config
    if base_config_path and Path(base_config_path).exists():
        config = load_config(base_config_path)
    else:
        config = default_config()

    # Override from scenario base_config
    n_years = base_cfg.get("n_years", 20)
    disease_year = base_cfg.get("disease_year", 3)

    # Attach release events to config
    config.release_events = release_events

    # Run replicates
    replicate_results = []
    for rep_idx, seed in enumerate(seeds):
        if verbose:
            print(f"\n  Replicate {rep_idx + 1}/{n_replicates} (seed={seed})...",
                  end="", flush=True)

        t0 = time.time()

        # Build fresh network per replicate
        network = build_scenario_network(base_cfg, seed=seed)

        # Run spatial simulation
        result = run_spatial_simulation(
            network=network,
            n_years=n_years,
            disease_year=disease_year,
            seed=seed,
            config=config,
        )

        elapsed = time.time() - t0

        # Extract results
        rep_result = extract_replicate_results(
            result, scenario_meta, seed, track_released=track_released,
        )
        rep_result["elapsed_seconds"] = elapsed
        replicate_results.append(rep_result)

        # Save individual replicate
        rep_file = out_dir / f"replicate_seed{seed}.json"
        with open(rep_file, "w") as f:
            json.dump(rep_result, f, indent=2, default=str)

        if verbose:
            m = rep_result.get("metrics", {})
            print(f" done ({elapsed:.1f}s) — "
                  f"final_pop={m.get('final_total_pop', '?')}, "
                  f"crash={m.get('crash_pct', '?'):.1f}%")

    # Compute ensemble summary
    summary = compute_ensemble_summary(replicate_results)
    summary["scenario_meta"] = scenario_meta
    summary["output_config"] = {
        k: v for k, v in output_config.items()
        if not isinstance(v, (np.ndarray,))  # JSON-safe
    }

    # Save summary
    summary_file = out_dir / "ensemble_summary.json"
    with open(summary_file, "w") as f:
        json.dump(summary, f, indent=2, default=str)

    if verbose:
        print(f"\n  Ensemble summary:")
        if "final_total_pop" in summary:
            s = summary["final_total_pop"]
            print(f"    Final pop: {s['mean']:.0f} ± {s['std']:.0f}")
        if "persistence_probability" in summary:
            print(f"    Persistence: {summary['persistence_probability']:.0%}")
        if "released_survival_pct" in summary:
            s = summary["released_survival_pct"]
            print(f"    Released survival: {s['mean']:.1f}% ± {s['std']:.1f}%")
        print(f"  Saved: {summary_file}")

    return summary


# ═══════════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Run reintroduction scenarios from YAML config files.",
        epilog="Example: python scripts/run_release_scenarios.py scenarios/small_pilot.yaml",
    )
    parser.add_argument(
        "scenarios", nargs="+",
        help="YAML scenario file(s) to run",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Parse and display config without running simulations",
    )
    parser.add_argument(
        "--output-dir", type=str, default=None,
        help="Override output directory (default: from YAML)",
    )
    parser.add_argument(
        "--network-size", type=int, default=5, choices=[5, 11],
        help="Network size for node resolution (default: 5)",
    )
    parser.add_argument(
        "--base-config", type=str, default=None,
        help="Base config YAML for model parameters",
    )
    parser.add_argument(
        "--quiet", action="store_true",
        help="Suppress progress output",
    )
    args = parser.parse_args()

    print("=" * 60)
    print("SSWD-EvoEpi Release Scenario Runner")
    print("=" * 60)

    all_summaries = []
    for scenario_path in args.scenarios:
        if not Path(scenario_path).exists():
            print(f"\n  ⚠️  File not found: {scenario_path}")
            continue

        summary = run_scenario_file(
            yaml_path=scenario_path,
            output_dir=args.output_dir,
            dry_run=args.dry_run,
            network_size=args.network_size,
            base_config_path=args.base_config,
            verbose=not args.quiet,
        )
        if summary is not None:
            all_summaries.append(summary)

    if all_summaries and len(all_summaries) > 1:
        print(f"\n{'=' * 60}")
        print(f"Ran {len(all_summaries)} scenarios successfully.")

    print("\n✅ Done.")


if __name__ == "__main__":
    main()
