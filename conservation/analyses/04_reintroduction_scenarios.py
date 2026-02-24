#!/usr/bin/env python3
"""Analysis 4: Reintroduction Scenarios (SKELETON).

Predicts outcomes of different reintroduction strategies using the
full spatial SSWD-EvoEpi model. This is the most computationally
expensive analysis (~2-5 days on Xeon).

STATUS: Skeleton template — requires model extensions before execution.

TODOs before this analysis can run:
  1. [ ] Release mechanism in model (inject individuals at node/time)
  2. [ ] Calibrated model parameters (ABC-SMC)
  3. [ ] Origin tracking in model output (wild vs captive-bred)
  4. [ ] Long-run stability validation (50-year runs)
  5. [ ] Computational budget allocation (Xeon scheduling)

Usage:
    python 04_reintroduction_scenarios.py [--dry-run] [--output-dir DIR]

See: conservation report Section 9.4
"""

import sys
import json
import warnings
import argparse
from pathlib import Path
from dataclasses import dataclass
from itertools import product

import numpy as np
import yaml

# ── Project imports ──────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from sswd_evoepi.types import N_LOCI, trait_slices


# ═══════════════════════════════════════════════════════════════════════
# SCENARIO DEFINITIONS
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class ReintroductionScenario:
    """Definition of a single reintroduction scenario."""
    name: str
    release_size: int           # Number of individuals per release
    release_node_ids: list      # Target node(s) for release
    release_year: int           # Year of first release
    release_month: int          # Month of release (1-12)
    release_frequency: str      # "one_time", "annual", "biennial"
    n_releases: int             # Total number of release events
    breeding_generations: int   # Generations of captive breeding
    breeding_strategy: str      # "truncation", "complementary", etc.
    forward_years: int          # Years to simulate after first release
    description: str = ""


def build_scenario_grid(params: dict) -> list:
    """Build the full factorial scenario grid.

    Returns:
        List of ReintroductionScenario objects.
    """
    reintro = params["reintroduction"]
    sites = params["sites"]

    # Release location strategies
    location_strategies = {
        "single_north": [0],             # Sitka only
        "single_south": [10],            # Monterey only
        "single_mid": [5],               # SJI only
        "multi_spread": [0, 5, 10],      # North + mid + south
        "stepping_stone": [2, 5, 8],     # Every ~3 nodes
    }

    scenarios = []
    scenario_id = 0

    for release_size in reintro["release_sizes"]:
        for breed_gen in reintro["release_generations"]:
            for freq in reintro["release_frequency"]:
                for loc_name, node_ids in location_strategies.items():
                    for fwd_years in reintro["forward_years"]:
                        # Spring release
                        n_releases = 1 if freq == "one_time" else \
                            min(5, fwd_years // (2 if freq == "biennial" else 1))

                        scenario = ReintroductionScenario(
                            name=f"S{scenario_id:04d}_{loc_name}_n{release_size}_g{breed_gen}_{freq}",
                            release_size=release_size,
                            release_node_ids=node_ids,
                            release_year=2028,  # Assumed program start
                            release_month=4,    # Spring release
                            release_frequency=freq,
                            n_releases=n_releases,
                            breeding_generations=breed_gen,
                            breeding_strategy="complementary",
                            forward_years=fwd_years,
                            description=f"{release_size} individuals at {loc_name}, "
                                        f"{breed_gen} gen breeding, {freq}, "
                                        f"{fwd_years}yr forward",
                        )
                        scenarios.append(scenario)
                        scenario_id += 1

    return scenarios


# ═══════════════════════════════════════════════════════════════════════
# MODEL RUN LOOP (SKELETON)
# ═══════════════════════════════════════════════════════════════════════

def run_scenario(
    scenario: ReintroductionScenario,
    params: dict,
    calibrated: dict,
    seed: int,
) -> dict:
    """Run a single reintroduction scenario.

    TODO: Implement once model extensions are ready.

    Steps:
    1. Load calibrated model
    2. Run 2013-2026 to establish post-epidemic baseline
    3. Generate captive-bred population (using Analysis 3 breeding)
    4. Schedule release events per scenario definition
    5. Run forward simulation
    6. Track outcomes

    Args:
        scenario: Scenario definition.
        params: Analysis parameters.
        calibrated: Calibrated model parameters.
        seed: Random seed.

    Returns:
        Dict with outcome metrics.
    """
    # TODO: Replace with actual model run
    warnings.warn(
        "⚠️  Reintroduction simulation not yet implemented. "
        "Returning placeholder results.",
        UserWarning,
        stacklevel=2,
    )

    rng = np.random.default_rng(seed)

    # Placeholder: generate fake outcome metrics
    # These mimic the expected output format
    n_years = scenario.forward_years
    years = np.arange(2028, 2028 + n_years + 1)

    # Fake population trajectory (logistic recovery with noise)
    k_total = params["population"]["K_per_node"] * len(scenario.release_node_ids)
    n_init = scenario.release_size * len(scenario.release_node_ids) * 0.01  # Crash
    pop = np.zeros(len(years))
    pop[0] = n_init
    for i in range(1, len(years)):
        growth = 0.1 * pop[i-1] * (1 - pop[i-1] / k_total)
        pop[i] = max(0, pop[i-1] + growth + rng.normal(0, pop[i-1] * 0.1))

    return {
        "scenario": scenario.name,
        "seed": seed,
        "years": years.tolist(),
        "total_population": pop.tolist(),
        "per_node_population": {},  # TODO: per-node trajectories
        "metrics": {
            "persistence_probability": float(pop[-1] > 10),
            "final_population": float(pop[-1]),
            "population_at_10yr": float(pop[min(10, n_years)]),
            "population_at_20yr": float(pop[min(20, n_years)]) if n_years >= 20 else None,
            "mean_resistance_final": rng.uniform(0.15, 0.35),  # TODO
            "allele_spread_nodes": 0,  # TODO: how many nodes have introduced alleles
            "time_to_recovery": None,  # TODO: years to reach 50% K
        },
        "is_placeholder": True,
    }


# ═══════════════════════════════════════════════════════════════════════
# OUTPUT FORMAT
# ═══════════════════════════════════════════════════════════════════════

EXPECTED_OUTPUT_FORMAT = """
Expected output format for each scenario × seed:
{
    "scenario": "S0001_single_north_n1000_g5_annual",
    "seed": 42,
    "years": [2028, 2029, ..., 2048],
    "total_population": [1000, 1234, ...],
    "per_node_population": {
        "Sitka": [1000, 900, ...],
        "Ketchikan": [0, 12, ...],  // Spread via larvae
        ...
    },
    "allele_frequencies": {
        "year_2030": {"Sitka": [0.2, 0.3, ...], ...},
        "year_2040": {...},
    },
    "mean_resistance": {
        "Sitka": [0.15, 0.16, ...],
        ...
    },
    "metrics": {
        "persistence_probability": 0.85,
        "final_population": 23456,
        "time_to_recovery": 15,  // years to 50% K
        "allele_spread_nodes": 5,
        "mean_resistance_final": 0.28,
    },
}

Ensemble summary (across seeds):
{
    "scenario": "S0001_...",
    "n_seeds": 50,
    "persistence_probability": 0.82 ± 0.05,
    "time_to_recovery": {"mean": 15.3, "std": 4.2, ...},
    ...
}
"""


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
    args = parser.parse_args()

    print("=" * 60)
    print("Analysis 4: Reintroduction Scenarios")
    print("=" * 60)

    params_path = args.params or str(Path(__file__).parent / "params.yaml")
    with open(params_path) as f:
        params = yaml.safe_load(f)

    # Build scenario grid
    scenarios = build_scenario_grid(params)
    n_seeds = args.seeds or params["reintroduction"]["n_seeds"]

    print(f"\n  Total scenarios: {len(scenarios)}")
    print(f"  Seeds per scenario: {n_seeds}")
    print(f"  Total runs: {len(scenarios) * n_seeds}")

    if args.dry_run:
        print("\n── SCENARIO GRID (dry run) ──")
        for s in scenarios[:20]:
            print(f"  {s.name}: {s.description}")
        if len(scenarios) > 20:
            print(f"  ... and {len(scenarios) - 20} more")
        print("\n  Estimated compute time (Xeon, 48 cores):")
        est_hours = len(scenarios) * n_seeds * 0.02  # ~1.2 min per run
        print(f"    {est_hours:.0f} core-hours "
              f"(~{est_hours/48:.1f} wall-hours)")
        return

    # Output directory
    output_dir = Path(args.output_dir) if args.output_dir else \
        Path(__file__).parent / ".." / params["output"]["analysis4"]
    output_dir.mkdir(parents=True, exist_ok=True)

    # Limit scenarios for testing
    if args.max_scenarios:
        scenarios = scenarios[:args.max_scenarios]
        print(f"  Limited to {args.max_scenarios} scenarios")

    print("\n⚠️  WARNING: Reintroduction simulation requires model extensions.")
    print("  Running with PLACEHOLDER results.\n")

    # TODO: Load calibrated params
    calibrated = None

    # Run scenarios
    all_results = []
    for i, scenario in enumerate(scenarios):
        print(f"  [{i+1}/{len(scenarios)}] {scenario.name}...", end="", flush=True)

        scenario_results = []
        for seed in range(n_seeds):
            result = run_scenario(scenario, params, calibrated, seed)
            scenario_results.append(result)

        # Ensemble summary
        persistence = np.mean([r["metrics"]["persistence_probability"]
                              for r in scenario_results])
        final_pop = np.mean([r["metrics"]["final_population"]
                            for r in scenario_results])

        all_results.append({
            "scenario": scenario.name,
            "config": {
                "release_size": scenario.release_size,
                "release_nodes": scenario.release_node_ids,
                "breeding_generations": scenario.breeding_generations,
                "frequency": scenario.release_frequency,
                "forward_years": scenario.forward_years,
            },
            "ensemble": {
                "persistence_probability": persistence,
                "final_population_mean": final_pop,
            },
            "n_seeds": n_seeds,
        })
        print(f" persistence={persistence:.2f}")

    # Save
    with open(output_dir / "reintroduction_results.json", "w") as f:
        json.dump({
            "results": all_results,
            "n_scenarios": len(scenarios),
            "n_seeds": n_seeds,
            "is_placeholder": True,
            "todos": [
                "Implement release mechanism in model",
                "Calibrate model parameters (ABC-SMC)",
                "Add origin tracking for released individuals",
                "Validate long-run stability (50yr)",
                "Run on Xeon with full scenario grid",
            ],
        }, f, indent=2)

    print(f"\n  Saved: {output_dir / 'reintroduction_results.json'}")
    print("\n✅ Analysis 4 complete (placeholder mode).")
    print("   Rerun after implementing model extensions.")


if __name__ == "__main__":
    main()
