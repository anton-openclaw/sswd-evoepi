#!/usr/bin/env python3
"""Analysis 3: Breeding Program Optimization.

Compares breeding strategies and determines optimal program design
for enhancing SSWD resistance in captive Pycnopodia populations.

Strategies compared:
  - Random mating of selected parents (baseline)
  - Assortative mating by resistance
  - Complementary mating (maximize locus coverage)
  - Optimal contribution selection (gain constrained by ΔF)

Tracks per generation: mean/max resistance, V_A, H_e, F, N_e.
Produces Pareto frontier of genetic gain vs diversity retention.

Usage:
    python 03_breeding_optimization.py [--input-dir DIR] [--output-dir DIR]

Dependencies: Analysis 1 results (for founder genotypes).
See: conservation report Section 9.3, Eq. 5.1-5.14
"""

import sys
import json
import warnings
import argparse
from pathlib import Path
from dataclasses import asdict

import numpy as np
import yaml

# ── Project imports ──────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from sswd_evoepi.types import N_LOCI, trait_slices
from sswd_evoepi.genetics import initialize_trait_effect_sizes

from conservation.src.breeding import (
    run_breeding_program,
    GenerationStats,
    BreedingResult,
)
from conservation.src.screening import select_complementary_founders
from conservation.src.inbreeding import (
    projected_f,
    generations_to_f_threshold,
)
from conservation.src.viz import (
    plot_breeding_trajectory,
    plot_breeding_comparison,
    plot_pareto_frontier,
)


# ═══════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════

def load_params(params_path: str = None) -> dict:
    if params_path is None:
        params_path = Path(__file__).parent / "params.yaml"
    with open(params_path) as f:
        return yaml.safe_load(f)


def load_founder_genotypes(input_dir: Path, site_name: str = None) -> tuple:
    """Load founder genotypes from Analysis 1.

    If site_name is given, loads from that site. Otherwise picks
    the site with highest mean resistance (best founders).

    Returns:
        (genotypes, effects_r, effects_t, effects_c, metadata)
    """
    summary_file = input_dir / "ensemble_summary.json"
    if not summary_file.exists():
        raise FileNotFoundError(
            f"Analysis 1 results not found: {summary_file}\n"
            "Run 01_current_genetic_state.py first."
        )

    with open(summary_file) as f:
        ensemble = json.load(f)

    # Find best site if not specified
    if site_name is None:
        summary = ensemble["summary"]
        site_name = max(
            summary.keys(),
            key=lambda s: summary[s]["mean_r"]["mean"]
        )
        print(f"  Auto-selected founder site: {site_name}")

    # Load site data
    safe_name = site_name.lower().replace(" ", "_")
    site_file = input_dir / f"{safe_name}.json"
    if not site_file.exists():
        raise FileNotFoundError(f"Site file not found: {site_file}")

    with open(site_file) as f:
        data = json.load(f)

    seed0 = data["seed_0_full"]
    genotypes = np.array(seed0["genotypes"], dtype=np.int8)
    effects_r = np.array(seed0["effects"]["resistance"])
    effects_t = np.array(seed0["effects"]["tolerance"])
    effects_c = np.array(seed0["effects"]["recovery"])

    return genotypes, effects_r, effects_t, effects_c, {
        "site": site_name,
        "n_available": seed0["n_survivors"],
        "is_placeholder": data.get("is_placeholder", True),
    }


# ═══════════════════════════════════════════════════════════════════════
# BREEDING EXPERIMENTS
# ═══════════════════════════════════════════════════════════════════════

def run_strategy_comparison(
    genotypes: np.ndarray,
    effects_r: np.ndarray,
    effects_t: np.ndarray,
    effects_c: np.ndarray,
    params: dict,
    n_founders: int,
    n_select: int,
    n_generations: int,
    n_replicates: int,
) -> dict:
    """Run breeding programs under all strategies.

    Args:
        genotypes: Full founder pool genotypes.
        effects_r/t/c: Effect size arrays.
        params: Analysis parameters.
        n_founders: Number of founders to draw.
        n_select: Parents selected per generation.
        n_generations: Breeding generations.
        n_replicates: Replicate seeds per strategy.

    Returns:
        Dict mapping strategy → list of BreedingResults.
    """
    gen = params["genetics"]
    breed = params["breeding"]
    n_r, n_t, n_c = gen["n_resistance"], gen["n_tolerance"], gen["n_recovery"]
    res_s, tol_s, rec_s = trait_slices(n_r, n_t, n_c)

    strategies = breed["strategies"]
    weights = breed["selection_weights"]["resistance_only"]
    results = {}

    for strategy in strategies:
        results[strategy] = []
        print(f"    {strategy}:", end="", flush=True)

        for rep in range(n_replicates):
            rng = np.random.default_rng(rep)

            # Draw founders (with replacement if needed)
            n_avail = len(genotypes)
            if n_founders <= n_avail:
                idx = rng.choice(n_avail, size=n_founders, replace=False)
            else:
                idx = rng.choice(n_avail, size=n_founders, replace=True)
            founders = genotypes[idx]

            # Run breeding program
            result = run_breeding_program(
                founder_genotypes=founders,
                effects_r=effects_r,
                effects_t=effects_t,
                effects_c=effects_c,
                rng=rng,
                n_generations=n_generations,
                n_selected=min(n_select, n_founders),
                n_offspring_per_pair=breed["n_offspring_per_pair"],
                n_keep_per_family=breed["n_keep_per_family"],
                scheme=strategy,
                n_r=n_r, n_t=n_t, n_c=n_c,
                w_r=weights["w_r"],
                w_t=weights["w_t"],
                w_c=weights["w_c"],
                target_ne=breed["target_ne"],
            )
            results[strategy].append(result)

            if (rep + 1) % 10 == 0:
                print(f" {rep+1}", end="", flush=True)

        print(" done.")

    return results


def summarize_strategies(results: dict) -> dict:
    """Compute summary statistics across replicates.

    Returns:
        Dict mapping strategy → {generations, mean_r, max_r, ...}
        with mean and std across replicates.
    """
    summary = {}

    for strategy, replicate_results in results.items():
        n_reps = len(replicate_results)
        n_gens = len(replicate_results[0].stats)

        # Collect per-gen metrics across replicates
        arrays = {
            "mean_r": np.zeros((n_reps, n_gens)),
            "max_r": np.zeros((n_reps, n_gens)),
            "va_r": np.zeros((n_reps, n_gens)),
            "he": np.zeros((n_reps, n_gens)),
            "f_mean": np.zeros((n_reps, n_gens)),
            "loci_fixed": np.zeros((n_reps, n_gens)),
        }

        for rep, result in enumerate(replicate_results):
            for g, stats in enumerate(result.stats):
                arrays["mean_r"][rep, g] = stats.mean_r
                arrays["max_r"][rep, g] = stats.max_r
                arrays["va_r"][rep, g] = stats.va_r
                arrays["he"][rep, g] = stats.he
                arrays["f_mean"][rep, g] = stats.f_mean
                arrays["loci_fixed"][rep, g] = stats.loci_fixed

        generations = np.arange(n_gens)

        summary[strategy] = {
            "generations": generations.tolist(),
            "n_replicates": n_reps,
        }
        for key, arr in arrays.items():
            summary[strategy][key] = {
                "mean": np.mean(arr, axis=0).tolist(),
                "std": np.std(arr, axis=0).tolist(),
                "min": np.min(arr, axis=0).tolist(),
                "max": np.max(arr, axis=0).tolist(),
            }

        # Endpoint gains
        delta_r = arr[:, -1] - arr[:, 0] if "mean_r" in key else \
            arrays["mean_r"][:, -1] - arrays["mean_r"][:, 0]
        summary[strategy]["endpoint"] = {
            "delta_r_mean": float(np.mean(delta_r)),
            "delta_r_std": float(np.std(delta_r)),
            "final_he_mean": float(np.mean(arrays["he"][:, -1])),
            "final_he_std": float(np.std(arrays["he"][:, -1])),
            "final_f_mean": float(np.mean(arrays["f_mean"][:, -1])),
            "final_f_std": float(np.std(arrays["f_mean"][:, -1])),
        }

    return summary


def compute_pareto_frontier(summary: dict) -> dict:
    """Compute gain–diversity Pareto frontier across strategies.

    Uses endpoint mean resistance gain vs final He retention.

    Returns:
        Dict with gain, diversity, labels, and pareto indices.
    """
    strategies = list(summary.keys())
    gains = np.array([summary[s]["endpoint"]["delta_r_mean"] for s in strategies])
    diversity = np.array([summary[s]["endpoint"]["final_he_mean"] for s in strategies])

    # Find Pareto-optimal
    n = len(strategies)
    is_pareto = np.ones(n, dtype=bool)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if gains[j] >= gains[i] and diversity[j] >= diversity[i]:
                if gains[j] > gains[i] or diversity[j] > diversity[i]:
                    is_pareto[i] = False
                    break

    return {
        "strategies": strategies,
        "gain": gains.tolist(),
        "diversity": diversity.tolist(),
        "pareto_optimal": np.where(is_pareto)[0].tolist(),
    }


# ═══════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Analysis 3: Breeding program optimization"
    )
    parser.add_argument("--input-dir", type=str, default=None,
                        help="Analysis 1 output directory")
    parser.add_argument("--output-dir", type=str, default=None,
                        help="Output directory")
    parser.add_argument("--params", type=str, default=None,
                        help="Path to params.yaml")
    parser.add_argument("--founders", type=int, default=100,
                        help="Number of founders (default: 100)")
    parser.add_argument("--generations", type=int, default=10,
                        help="Breeding generations (default: 10)")
    parser.add_argument("--replicates", type=int, default=None,
                        help="Seeds per strategy (default: from params)")
    parser.add_argument("--site", type=str, default=None,
                        help="Source site for founders (default: best)")
    args = parser.parse_args()

    print("=" * 60)
    print("Analysis 3: Breeding Program Optimization")
    print("=" * 60)

    params = load_params(args.params)
    breed = params["breeding"]

    # Directories
    input_dir = Path(args.input_dir) if args.input_dir else \
        Path(__file__).parent / ".." / params["output"]["analysis1"]
    output_dir = Path(args.output_dir) if args.output_dir else \
        Path(__file__).parent / ".." / params["output"]["analysis3"]
    output_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = output_dir / "figures"
    fig_dir.mkdir(exist_ok=True)

    n_reps = args.replicates or breed["n_replicates"]

    # Load founders
    print("\nLoading founder genotypes...")
    genotypes, effects_r, effects_t, effects_c, meta = \
        load_founder_genotypes(input_dir, args.site)
    print(f"  Source: {meta['site']} ({meta['n_available']} individuals)")
    if meta["is_placeholder"]:
        warnings.warn(
            "⚠️  Founders from PLACEHOLDER data. "
            "Results are illustrative only.",
            UserWarning,
        )

    n_select = min(20, args.founders // 2)

    # Run strategy comparison
    print(f"\nRunning breeding programs "
          f"({args.founders} founders, {args.generations} gen, "
          f"{n_reps} reps, {n_select} selected/gen)...")
    results = run_strategy_comparison(
        genotypes, effects_r, effects_t, effects_c,
        params, args.founders, n_select,
        args.generations, n_reps,
    )

    # Summarize
    print("\nSummarizing results...")
    summary = summarize_strategies(results)

    # Pareto frontier
    pareto = compute_pareto_frontier(summary)

    # ── Figures ──────────────────────────────────────────────────────
    print("\nGenerating figures...")
    import matplotlib.pyplot as plt

    # Comparison plot: mean resistance
    comparison_data = {}
    for strategy in summary:
        comparison_data[strategy] = {
            "generations": np.array(summary[strategy]["generations"]),
            "mean_r": np.array(summary[strategy]["mean_r"]["mean"]),
            "max_r": np.array(summary[strategy]["max_r"]["mean"]),
            "va_r": np.array(summary[strategy]["va_r"]["mean"]),
            "he": np.array(summary[strategy]["he"]["mean"]),
            "f_mean": np.array(summary[strategy]["f_mean"]["mean"]),
        }

    fig = plot_breeding_comparison(comparison_data, metric="mean_r")
    fig.savefig(fig_dir / "strategy_comparison_mean_r.png", dpi=150,
                bbox_inches="tight")
    plt.close(fig)

    fig = plot_breeding_comparison(comparison_data, metric="he")
    fig.savefig(fig_dir / "strategy_comparison_he.png", dpi=150,
                bbox_inches="tight")
    plt.close(fig)

    fig = plot_breeding_comparison(comparison_data, metric="f_mean")
    fig.savefig(fig_dir / "strategy_comparison_f.png", dpi=150,
                bbox_inches="tight")
    plt.close(fig)

    # Per-strategy trajectory (best strategy)
    best_strategy = max(summary.keys(),
                        key=lambda s: summary[s]["endpoint"]["delta_r_mean"])
    best = comparison_data[best_strategy]
    fig = plot_breeding_trajectory(
        best["generations"],
        best["mean_r"], best["max_r"], best["va_r"],
        best["he"], best["f_mean"],
        loci_fixed=np.array(summary[best_strategy]["loci_fixed"]["mean"]),
    )
    fig.suptitle(f"Best Strategy: {best_strategy}", fontsize=14,
                 fontweight="bold", y=1.01)
    fig.savefig(fig_dir / f"trajectory_{best_strategy}.png", dpi=150,
                bbox_inches="tight")
    plt.close(fig)

    # Pareto frontier
    fig = plot_pareto_frontier(
        np.array(pareto["gain"]),
        np.array(pareto["diversity"]),
        labels=pareto["strategies"],
        pareto_idx=np.array(pareto["pareto_optimal"]),
    )
    fig.savefig(fig_dir / "pareto_gain_diversity.png", dpi=150,
                bbox_inches="tight")
    plt.close(fig)

    # ── Save results ─────────────────────────────────────────────────
    print("\nSaving results...")
    output = {
        "summary": summary,
        "pareto": pareto,
        "config": {
            "n_founders": args.founders,
            "n_select": n_select,
            "n_generations": args.generations,
            "n_replicates": n_reps,
            "source_site": meta["site"],
        },
        "is_placeholder": meta["is_placeholder"],
    }
    with open(output_dir / "breeding_results.json", "w") as f:
        json.dump(output, f, indent=2)
    print(f"  Saved: {output_dir / 'breeding_results.json'}")

    # Print summary
    print("\n" + "=" * 70)
    print("BREEDING STRATEGY COMPARISON")
    print("=" * 70)
    print(f"{'Strategy':<18} {'Δr̄':>8} {'Final Hₑ':>10} {'Final F':>10} "
          f"{'Pareto':>8}")
    print("-" * 70)
    for i, strategy in enumerate(pareto["strategies"]):
        ep = summary[strategy]["endpoint"]
        is_p = "✓" if i in pareto["pareto_optimal"] else ""
        print(f"{strategy:<18} {ep['delta_r_mean']:>8.4f} "
              f"{ep['final_he_mean']:>10.4f} "
              f"{ep['final_f_mean']:>10.4f} "
              f"{is_p:>8}")

    print(f"\nBest gain: {best_strategy} (Δr̄ = "
          f"{summary[best_strategy]['endpoint']['delta_r_mean']:.4f})")
    print("\n✅ Analysis 3 complete.")


if __name__ == "__main__":
    main()
