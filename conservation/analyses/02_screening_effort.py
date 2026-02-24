#!/usr/bin/env python3
"""Analysis 2: Screening Effort by Site.

For each site's predicted 2026 population (from Analysis 1), computes:
  - Required sample sizes for various resistance thresholds
  - Expected best-of-n curves (diminishing returns)
  - Multi-site optimal screening allocation
  - Complementarity analysis of top founder candidates

Produces:
  - Per-site screening effort tables
  - Cross-site comparison figures
  - Optimal allocation bar chart

Usage:
    python 02_screening_effort.py [--input-dir DIR] [--output-dir DIR]

Dependencies: Analysis 1 results must exist.
See: conservation report Section 9.2, Eq. 4.1-4.10
"""

import sys
import json
import warnings
import argparse
from pathlib import Path

import numpy as np
import yaml

# ── Project imports ──────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from sswd_evoepi.types import N_LOCI, trait_slices
from sswd_evoepi.genetics import compute_trait_batch

from conservation.src.screening import (
    required_sample_size,
    required_n_for_k,
    empirical_exceedance,
    exceedance_curve,
    expected_max_normal,
    screening_effort_curve,
    multisite_allocation,
    select_complementary_founders,
)
from conservation.src.viz import (
    plot_screening_effort,
    plot_multisite_allocation,
    plot_exceedance,
    plot_complementarity_heatmap,
)


# ═══════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════

def load_params(params_path: str = None) -> dict:
    if params_path is None:
        params_path = Path(__file__).parent / "params.yaml"
    with open(params_path) as f:
        return yaml.safe_load(f)


def load_analysis1_results(input_dir: Path) -> dict:
    """Load Analysis 1 per-site results.

    Returns:
        Dict mapping site_name → full result dict (seed 0).
    """
    results = {}
    summary_file = input_dir / "ensemble_summary.json"

    if not summary_file.exists():
        raise FileNotFoundError(
            f"Analysis 1 results not found: {summary_file}\n"
            "Run 01_current_genetic_state.py first."
        )

    with open(summary_file) as f:
        ensemble = json.load(f)

    # Load per-site data (seed 0 for genotypes)
    for site_file in sorted(input_dir.glob("*.json")):
        if site_file.name == "ensemble_summary.json":
            continue
        with open(site_file) as f:
            data = json.load(f)
        site_name = data["site"]
        results[site_name] = data

    if not results:
        raise FileNotFoundError(
            f"No per-site JSON files found in {input_dir}"
        )

    return results, ensemble


# ═══════════════════════════════════════════════════════════════════════
# SCREENING ANALYSIS
# ═══════════════════════════════════════════════════════════════════════

def compute_screening_table(
    genotypes: np.ndarray,
    effects_r: np.ndarray,
    res_slice: slice,
    thresholds: list,
    confidence_levels: list,
) -> dict:
    """Compute screening effort table for one site.

    For each (threshold, confidence) pair, returns:
      - Empirical exceedance probability
      - Required sample size for ≥1 individual
      - Required sample size for ≥5 individuals

    Args:
        genotypes: (N, N_LOCI, 2) int8.
        effects_r: Resistance effect sizes.
        res_slice: Resistance locus slice.
        thresholds: Resistance thresholds to evaluate.
        confidence_levels: Confidence levels for sample size.

    Returns:
        Nested dict: {threshold: {confidence: {n_1, n_5, exceedance}}}.
    """
    table = {}

    for thresh in thresholds:
        # Empirical exceedance
        p_exceed = empirical_exceedance(
            genotypes, effects_r, res_slice, thresh
        )

        table[thresh] = {
            "exceedance_prob": p_exceed,
            "sample_sizes": {},
        }

        for conf in confidence_levels:
            if p_exceed > 0:
                n_1 = required_sample_size(p_exceed, conf)
                n_5 = required_n_for_k(p_exceed, k=5, confidence=conf)
            else:
                n_1 = float("inf")
                n_5 = float("inf")

            table[thresh]["sample_sizes"][conf] = {
                "n_for_1": n_1,
                "n_for_5": n_5,
            }

    return table


def compute_effort_curves(
    genotypes: np.ndarray,
    effects_r: np.ndarray,
    res_slice: slice,
    sample_sizes: list,
    rng: np.random.Generator,
    mc_reps: int = 500,
) -> dict:
    """Compute expected maximum resistance vs sample size.

    Args:
        genotypes: (N, N_LOCI, 2) int8.
        effects_r: Resistance effect sizes.
        res_slice: Resistance locus slice.
        sample_sizes: Sample sizes to evaluate.
        rng: Random generator.
        mc_reps: MC replicates.

    Returns:
        Dict with sample_sizes and expected_maxes arrays.
    """
    sizes = np.array(sample_sizes)
    expected_maxes = screening_effort_curve(
        genotypes, effects_r, res_slice, sizes, rng, mc_reps
    )

    return {
        "sample_sizes": sizes.tolist(),
        "expected_maxes": expected_maxes.tolist(),
    }


def compute_optimal_allocation(
    site_results: dict,
    params: dict,
) -> dict:
    """Compute optimal multi-site screening allocation.

    Args:
        site_results: Dict mapping site_name → screening table.
        params: Analysis parameters.

    Returns:
        Dict with allocation per site and expected returns.
    """
    screening_params = params["screening"]
    total_budget = screening_params["total_budget"]

    site_names = sorted(site_results.keys())
    n_sites = len(site_names)

    # Use mean and std of resistance at each site
    site_means = np.zeros(n_sites)
    site_stds = np.zeros(n_sites)
    site_max_n = np.zeros(n_sites, dtype=int)

    for i, name in enumerate(site_names):
        r = site_results[name]
        site_means[i] = r["traits"]["resistance"]["mean"]
        site_stds[i] = r["traits"]["resistance"]["std"]
        site_max_n[i] = r["n_survivors"]

    # Optimal allocation
    alloc = multisite_allocation(
        site_means, site_stds, total_budget, site_max_n
    )

    # Expected max at each site given allocation
    expected_maxes = np.array([
        expected_max_normal(site_means[i], site_stds[i], alloc[i])
        for i in range(n_sites)
    ])

    return {
        "site_names": site_names,
        "allocations": alloc.tolist(),
        "site_means": site_means.tolist(),
        "site_stds": site_stds.tolist(),
        "expected_maxes": expected_maxes.tolist(),
        "total_budget": total_budget,
        "best_site": site_names[int(np.argmax(expected_maxes))],
        "best_expected_max": float(np.max(expected_maxes)),
    }


# ═══════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Analysis 2: Screening effort by site"
    )
    parser.add_argument("--input-dir", type=str, default=None,
                        help="Analysis 1 output directory")
    parser.add_argument("--output-dir", type=str, default=None,
                        help="Output directory")
    parser.add_argument("--params", type=str, default=None,
                        help="Path to params.yaml")
    args = parser.parse_args()

    print("=" * 60)
    print("Analysis 2: Screening Effort by Site")
    print("=" * 60)

    params = load_params(args.params)
    gen = params["genetics"]
    scr = params["screening"]
    n_r, n_t, n_c = gen["n_resistance"], gen["n_tolerance"], gen["n_recovery"]
    res_s, tol_s, rec_s = trait_slices(n_r, n_t, n_c)

    # Directories
    input_dir = Path(args.input_dir) if args.input_dir else \
        Path(__file__).parent / ".." / params["output"]["analysis1"]
    output_dir = Path(args.output_dir) if args.output_dir else \
        Path(__file__).parent / ".." / params["output"]["analysis2"]
    output_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = output_dir / "figures"
    fig_dir.mkdir(exist_ok=True)

    # Load Analysis 1 results
    print("\nLoading Analysis 1 results...")
    site_data, ensemble = load_analysis1_results(input_dir)
    print(f"  Loaded {len(site_data)} sites")

    if ensemble.get("is_placeholder", True):
        warnings.warn(
            "⚠️  Analysis 1 used PLACEHOLDER data. "
            "Screening results are illustrative only.",
            UserWarning,
        )

    rng = np.random.default_rng(42)
    all_screening = {}

    # Per-site screening analysis
    print("\nComputing per-site screening tables...")
    for site_name, data in sorted(site_data.items()):
        print(f"  {site_name}...", end="", flush=True)

        seed0 = data["seed_0_full"]
        genotypes = np.array(seed0["genotypes"], dtype=np.int8)
        effects_r = np.array(seed0["effects"]["resistance"])

        # Screening table
        table = compute_screening_table(
            genotypes, effects_r, res_s,
            scr["resistance_thresholds"],
            scr["confidence_levels"],
        )

        # Effort curves
        effort = compute_effort_curves(
            genotypes, effects_r, res_s,
            scr["sample_sizes"], rng, scr["mc_replicates"],
        )

        # Exceedance curve
        thresholds, exc_probs = exceedance_curve(
            genotypes, effects_r, res_s, n_points=200
        )

        all_screening[site_name] = {
            "table": table,
            "effort": effort,
            "exceedance": {
                "thresholds": thresholds.tolist(),
                "probs": exc_probs.tolist(),
            },
            "n_available": seed0["n_survivors"],
            "traits": seed0["traits"],
        }

        # Per-site figures
        fig = plot_screening_effort(
            np.array(effort["sample_sizes"]),
            np.array(effort["expected_maxes"]),
            trait_name="resistance",
            target_value=0.30,
        )
        safe_name = site_name.lower().replace(" ", "_")
        fig.savefig(fig_dir / f"effort_{safe_name}.png", dpi=150,
                    bbox_inches="tight")
        import matplotlib.pyplot as plt
        plt.close(fig)

        fig = plot_exceedance(thresholds, exc_probs, target_prob=0.05)
        fig.savefig(fig_dir / f"exceedance_{safe_name}.png", dpi=150,
                    bbox_inches="tight")
        plt.close(fig)

        print(" done.")

    # Multi-site optimal allocation
    print("\nComputing optimal multi-site allocation...")
    seed0_results = {
        name: data["seed_0_full"] for name, data in site_data.items()
    }
    allocation = compute_optimal_allocation(seed0_results, params)

    # Allocation figure
    fig = plot_multisite_allocation(
        allocation["site_names"],
        np.array(allocation["allocations"]),
        np.array(allocation["site_means"]),
    )
    fig.savefig(fig_dir / "optimal_allocation.png", dpi=150,
                bbox_inches="tight")
    plt.close(fig)
    print(f"  Best site: {allocation['best_site']} "
          f"(E[max r] = {allocation['best_expected_max']:.3f})")

    # Save all results
    print("\nSaving results...")
    results = {
        "screening": all_screening,
        "allocation": allocation,
        "is_placeholder": ensemble.get("is_placeholder", True),
        "params": {
            "thresholds": scr["resistance_thresholds"],
            "confidence_levels": scr["confidence_levels"],
            "total_budget": scr["total_budget"],
        },
    }

    with open(output_dir / "screening_results.json", "w") as f:
        json.dump(results, f, indent=2)

    # Print summary table
    print("\n" + "=" * 90)
    print("SCREENING EFFORT SUMMARY (95% confidence, k≥1)")
    print("=" * 90)
    print(f"{'Site':<16} {'N_avail':>8} {'r̄':>7} "
          + "".join(f"{'r*='+str(t):>10}" for t in scr["resistance_thresholds"]))
    print("-" * 90)
    for site_name in sorted(all_screening.keys(),
                            key=lambda s: site_data[s]["seed_0_full"]["lat"],
                            reverse=True):
        s = all_screening[site_name]
        row = f"{site_name:<16} {s['n_available']:>8} "
        row += f"{s['traits']['resistance']['mean']:>7.3f} "
        for t in scr["resistance_thresholds"]:
            n_req = s["table"][t]["sample_sizes"][0.95]["n_for_1"]
            if n_req == float("inf"):
                row += f"{'∞':>10}"
            else:
                row += f"{n_req:>10}"
        print(row)

    print(f"\n{'Optimal allocation':>16} (budget={scr['total_budget']}):")
    for name, alloc in zip(allocation["site_names"], allocation["allocations"]):
        print(f"  {name:<16} → {alloc:>4} samples")

    print("\n✅ Analysis 2 complete.")


if __name__ == "__main__":
    main()
