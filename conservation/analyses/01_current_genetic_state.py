#!/usr/bin/env python3
"""Analysis 1: Current Genetic State of Wild Pycnopodia (2026).

Runs the full SSWD-EvoEpi model from 2013 (pre-epidemic) to 2026
for each of the 11 stepping-stone sites. Extracts endpoint genotypes
and computes trait distributions, genetic diversity metrics, and
allele frequencies.

Produces:
  - Per-site JSON results with trait distributions and diversity
  - Latitude × trait heatmap figure
  - Ensemble summary statistics (across replicate seeds)

Designed to run with default parameters (with warnings) until
calibrated params are available.

Usage:
    python 01_current_genetic_state.py [--seeds N] [--output-dir DIR]

See: conservation report Section 9.1, Eq. 4.2-4.3
"""

import sys
import os
import json
import warnings
import argparse
from pathlib import Path
from datetime import datetime

import numpy as np
import yaml

# ── Project imports ──────────────────────────────────────────────────
# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from sswd_evoepi.types import N_LOCI, trait_slices
from sswd_evoepi.genetics import (
    compute_trait_batch,
    initialize_trait_effect_sizes,
    compute_allele_frequencies,
    compute_heterozygosity,
    compute_additive_variance,
)

# Conservation module
from conservation.src.trait_math import (
    trait_mean,
    trait_variance,
    exceedance_probability,
)
from conservation.src.inbreeding import (
    genomic_inbreeding,
    mean_inbreeding,
    expected_heterozygosity,
    observed_heterozygosity,
    allelic_richness,
    ne_from_grm,
)


# ═══════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════

def load_params(params_path: str = None) -> dict:
    """Load analysis parameters from YAML config."""
    if params_path is None:
        params_path = Path(__file__).parent / "params.yaml"
    with open(params_path) as f:
        return yaml.safe_load(f)


def load_calibrated_params(params: dict) -> dict:
    """Load calibrated model parameters if available.

    Returns:
        Calibrated parameter dict, or None with a warning.
    """
    cal_file = params.get("calibrated_params_file")
    if cal_file is None:
        warnings.warn(
            "⚠️  No calibrated parameters available. "
            "Using default model parameters. "
            "Results are PLACEHOLDER — do not cite.\n"
            "Set 'calibrated_params_file' in params.yaml after ABC-SMC calibration.",
            UserWarning,
            stacklevel=2,
        )
        return None

    cal_path = Path(__file__).parent / cal_file
    if not cal_path.exists():
        warnings.warn(
            f"⚠️  Calibrated params file not found: {cal_path}\n"
            "Falling back to default parameters.",
            UserWarning,
            stacklevel=2,
        )
        return None

    with open(cal_path) as f:
        return json.load(f)


# ═══════════════════════════════════════════════════════════════════════
# MODEL SIMULATION
# ═══════════════════════════════════════════════════════════════════════

def run_single_site(
    site_info: dict,
    params: dict,
    calibrated: dict,
    seed: int,
) -> dict:
    """Run the model for one site from 2013 to 2026.

    TODO: Plug in full model run once calibrated.
    Currently generates synthetic placeholder data using the
    genetic architecture to demonstrate the analysis pipeline.

    Args:
        site_info: Site definition from params.yaml.
        params: Full parameter dictionary.
        calibrated: Calibrated params (or None for defaults).
        seed: Random seed.

    Returns:
        Dict with genotypes, trait scores, and metadata.
    """
    rng = np.random.default_rng(seed)

    # ── Genetic parameters ───────────────────────────────────────────
    gen = params["genetics"]
    n_r = gen["n_resistance"]
    n_t = gen["n_tolerance"]
    n_c = gen["n_recovery"]
    res_s, tol_s, rec_s = trait_slices(n_r, n_t, n_c)

    # TODO: Replace with calibrated effect sizes when available
    effects_r = initialize_trait_effect_sizes(n_r, gen["target_mean_r"], rng)
    effects_t = initialize_trait_effect_sizes(n_t, gen["target_mean_t"], rng)
    effects_c = initialize_trait_effect_sizes(n_c, gen["target_mean_c"], rng)

    # ── Generate placeholder population ──────────────────────────────
    # TODO: Replace with actual model run:
    #   1. Initialize pre-epidemic population
    #   2. config = build_config(params, calibrated, site_info)
    #   3. model = SpatialModel(config)
    #   4. model.run(start=2013, end=2026, seed=seed)
    #   5. Extract endpoint genotypes from model.agents

    # Synthetic: crash to ~2% survival, shift allele frequencies
    K = params["population"]["K_per_node"]
    n_survivors = max(10, int(K * rng.uniform(0.005, 0.030)))

    # Temperature-dependent crash severity (warmer → worse)
    sst = site_info.get("mean_sst", 10.0)
    crash_severity = 1.0 - np.clip((sst - 7.0) / 8.0, 0, 1) * 0.02

    genotypes = np.zeros((n_survivors, N_LOCI, 2), dtype=np.int8)

    # Initialize with allele frequencies shifted by selection
    # Higher resistance frequencies at cold sites (stronger selection signal)
    base_q_r = gen["target_mean_r"] + rng.normal(0, 0.02)
    selection_shift = rng.uniform(0.02, 0.08)  # Post-epidemic selection
    q_r = np.clip(base_q_r + selection_shift, 0.01, 0.99)

    for locus in range(N_LOCI):
        if locus < n_r:
            q = q_r + rng.normal(0, 0.03)  # Resistance loci
        elif locus < n_r + n_t:
            q = gen["target_mean_t"] + rng.normal(0, 0.02)  # Tolerance
        else:
            q = gen["target_mean_c"] + rng.normal(0, 0.01)  # Recovery

        q = np.clip(q, 0.01, 0.99)
        for i in range(n_survivors):
            genotypes[i, locus, 0] = rng.binomial(1, q)
            genotypes[i, locus, 1] = rng.binomial(1, q)

    # ── Compute traits ───────────────────────────────────────────────
    alive = np.ones(n_survivors, dtype=bool)
    r_scores = compute_trait_batch(genotypes, effects_r, alive, res_s)
    t_scores = compute_trait_batch(genotypes, effects_t, alive, tol_s)
    c_scores = compute_trait_batch(genotypes, effects_c, alive, rec_s)

    # ── Diversity metrics ────────────────────────────────────────────
    he = expected_heterozygosity(genotypes)
    ho = observed_heterozygosity(genotypes)
    f_mean = mean_inbreeding(genotypes)
    n_poly = allelic_richness(genotypes)
    ne_est = ne_from_grm(genotypes) if n_survivors >= 5 else float(n_survivors)

    # ── Allele frequencies ───────────────────────────────────────────
    allele_freqs = genotypes.sum(axis=(0, 2)).astype(np.float64) / (2.0 * n_survivors)

    return {
        "site": site_info["name"],
        "site_id": site_info["id"],
        "lat": site_info["lat"],
        "lon": site_info["lon"],
        "seed": seed,
        "n_survivors": n_survivors,
        "is_placeholder": calibrated is None,
        "traits": {
            "resistance": {
                "mean": float(np.mean(r_scores)),
                "std": float(np.std(r_scores)),
                "min": float(np.min(r_scores)),
                "max": float(np.max(r_scores)),
                "median": float(np.median(r_scores)),
                "scores": r_scores.tolist(),
            },
            "tolerance": {
                "mean": float(np.mean(t_scores)),
                "std": float(np.std(t_scores)),
                "min": float(np.min(t_scores)),
                "max": float(np.max(t_scores)),
                "median": float(np.median(t_scores)),
                "scores": t_scores.tolist(),
            },
            "recovery": {
                "mean": float(np.mean(c_scores)),
                "std": float(np.std(c_scores)),
                "min": float(np.min(c_scores)),
                "max": float(np.max(c_scores)),
                "median": float(np.median(c_scores)),
                "scores": c_scores.tolist(),
            },
        },
        "diversity": {
            "He": he,
            "Ho": ho,
            "F_mean": f_mean,
            "n_polymorphic": n_poly,
            "Ne_estimate": ne_est,
        },
        "allele_frequencies": allele_freqs.tolist(),
        "genotypes_shape": list(genotypes.shape),
        # Store genotypes as nested list for JSON serialization
        "genotypes": genotypes.tolist(),
        "effects": {
            "resistance": effects_r.tolist(),
            "tolerance": effects_t.tolist(),
            "recovery": effects_c.tolist(),
        },
    }


# ═══════════════════════════════════════════════════════════════════════
# ENSEMBLE ANALYSIS
# ═══════════════════════════════════════════════════════════════════════

def run_ensemble(params: dict, calibrated: dict, n_seeds: int = None) -> dict:
    """Run model for all sites across multiple seeds.

    Args:
        params: Analysis parameters.
        calibrated: Calibrated params (or None).
        n_seeds: Number of replicate seeds (overrides params).

    Returns:
        Nested dict: {site_name: {seed: result_dict}}.
    """
    if n_seeds is None:
        n_seeds = params["simulation"]["n_seeds"]

    sites = params["sites"]
    results = {}

    for site in sites:
        site_name = site["name"]
        results[site_name] = {}
        print(f"  Running {site_name} ({n_seeds} seeds)...", end="", flush=True)

        for seed in range(n_seeds):
            results[site_name][seed] = run_single_site(
                site, params, calibrated, seed=seed
            )

        print(" done.")

    return results


def summarize_ensemble(results: dict) -> dict:
    """Compute ensemble statistics across seeds.

    Args:
        results: Output from run_ensemble.

    Returns:
        Dict with per-site summaries (means ± std across seeds).
    """
    summary = {}

    for site_name, seed_results in results.items():
        seeds = sorted(seed_results.keys())
        n = len(seeds)

        # Collect per-seed metrics
        metrics = {
            "n_survivors": [],
            "mean_r": [], "mean_t": [], "mean_c": [],
            "max_r": [], "max_t": [], "max_c": [],
            "He": [], "Ho": [], "F_mean": [], "Ne": [],
        }

        for s in seeds:
            r = seed_results[s]
            metrics["n_survivors"].append(r["n_survivors"])
            metrics["mean_r"].append(r["traits"]["resistance"]["mean"])
            metrics["mean_t"].append(r["traits"]["tolerance"]["mean"])
            metrics["mean_c"].append(r["traits"]["recovery"]["mean"])
            metrics["max_r"].append(r["traits"]["resistance"]["max"])
            metrics["max_t"].append(r["traits"]["tolerance"]["max"])
            metrics["max_c"].append(r["traits"]["recovery"]["max"])
            metrics["He"].append(r["diversity"]["He"])
            metrics["Ho"].append(r["diversity"]["Ho"])
            metrics["F_mean"].append(r["diversity"]["F_mean"])
            metrics["Ne"].append(r["diversity"]["Ne_estimate"])

        summary[site_name] = {
            "lat": seed_results[seeds[0]]["lat"],
            "lon": seed_results[seeds[0]]["lon"],
            "n_seeds": n,
        }
        for key, vals in metrics.items():
            arr = np.array(vals)
            summary[site_name][key] = {
                "mean": float(np.mean(arr)),
                "std": float(np.std(arr)),
                "min": float(np.min(arr)),
                "max": float(np.max(arr)),
                "median": float(np.median(arr)),
            }

    return summary


# ═══════════════════════════════════════════════════════════════════════
# VISUALIZATION
# ═══════════════════════════════════════════════════════════════════════

def plot_latitude_heatmap(summary: dict, output_dir: Path):
    """Generate latitude × trait heatmap.

    Shows mean trait values at each site, ordered by latitude.
    The key figure for visualizing the latitudinal gradient in
    post-epidemic genetic state.
    """
    import matplotlib.pyplot as plt

    sites = sorted(summary.keys(), key=lambda s: summary[s]["lat"], reverse=True)
    lats = [summary[s]["lat"] for s in sites]

    traits = ["mean_r", "mean_t", "mean_c"]
    trait_labels = ["Resistance", "Tolerance", "Recovery"]

    data = np.zeros((len(sites), len(traits)))
    for i, site in enumerate(sites):
        for j, trait in enumerate(traits):
            data[i, j] = summary[site][trait]["mean"]

    fig, ax = plt.subplots(figsize=(8, 10))
    im = ax.imshow(data, aspect="auto", cmap="YlOrRd")

    ax.set_yticks(range(len(sites)))
    ax.set_yticklabels([f"{s} ({lats[i]:.1f}°N)" for i, s in enumerate(sites)],
                       fontsize=10)
    ax.set_xticks(range(len(traits)))
    ax.set_xticklabels(trait_labels, fontsize=11, fontweight="bold")

    # Annotate cells
    for i in range(len(sites)):
        for j in range(len(traits)):
            val = data[i, j]
            text_color = "white" if val > data.mean() else "black"
            ax.text(j, i, f"{val:.3f}", ha="center", va="center",
                    fontsize=9, color=text_color)

    cbar = fig.colorbar(im, ax=ax, shrink=0.6, label="Trait Mean")
    ax.set_title("Post-Epidemic Trait Values by Site (2026)\n"
                 "11-Node Stepping-Stone Network",
                 fontsize=13, fontweight="bold")

    # Add placeholder warning if applicable
    is_placeholder = any(
        summary[s].get("is_placeholder", True)
        for s in sites
    )
    if is_placeholder:
        fig.text(0.5, 0.02,
                 "⚠ PLACEHOLDER — Default parameters, not calibrated",
                 ha="center", fontsize=10, color="red", style="italic")

    fig.tight_layout()
    fig.savefig(output_dir / "latitude_trait_heatmap.png", dpi=150,
                bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_dir / 'latitude_trait_heatmap.png'}")


def plot_diversity_gradient(summary: dict, output_dir: Path):
    """Plot genetic diversity metrics along the latitudinal gradient."""
    import matplotlib.pyplot as plt

    sites = sorted(summary.keys(), key=lambda s: summary[s]["lat"])
    lats = [summary[s]["lat"] for s in sites]

    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)

    # He
    he_means = [summary[s]["He"]["mean"] for s in sites]
    he_stds = [summary[s]["He"]["std"] for s in sites]
    axes[0].errorbar(lats, he_means, yerr=he_stds, fmt="o-",
                     color="#2196F3", capsize=3)
    axes[0].set_ylabel("Expected Heterozygosity (Hₑ)", fontsize=11)
    axes[0].set_title("Genetic Diversity Gradient", fontsize=13,
                      fontweight="bold")

    # Ne
    ne_means = [summary[s]["Ne"]["mean"] for s in sites]
    ne_stds = [summary[s]["Ne"]["std"] for s in sites]
    axes[1].errorbar(lats, ne_means, yerr=ne_stds, fmt="s-",
                     color="#4CAF50", capsize=3)
    axes[1].set_ylabel("Effective Pop. Size (Nₑ)", fontsize=11)

    # N survivors
    n_means = [summary[s]["n_survivors"]["mean"] for s in sites]
    n_stds = [summary[s]["n_survivors"]["std"] for s in sites]
    axes[2].errorbar(lats, n_means, yerr=n_stds, fmt="^-",
                     color="#FF9800", capsize=3)
    axes[2].set_ylabel("N Survivors", fontsize=11)
    axes[2].set_xlabel("Latitude (°N)", fontsize=11)

    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.savefig(output_dir / "diversity_gradient.png", dpi=150,
                bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_dir / 'diversity_gradient.png'}")


# ═══════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Analysis 1: Current genetic state of wild Pycnopodia"
    )
    parser.add_argument("--seeds", type=int, default=None,
                        help="Number of replicate seeds (default: from params.yaml)")
    parser.add_argument("--output-dir", type=str, default=None,
                        help="Output directory (default: from params.yaml)")
    parser.add_argument("--params", type=str, default=None,
                        help="Path to params.yaml")
    args = parser.parse_args()

    print("=" * 60)
    print("Analysis 1: Current Genetic State (2013 → 2026)")
    print("=" * 60)

    # Load configuration
    params = load_params(args.params)
    calibrated = load_calibrated_params(params)

    # Output directory
    if args.output_dir:
        output_dir = Path(args.output_dir)
    else:
        output_dir = Path(__file__).parent / ".." / params["output"]["analysis1"]
    output_dir.mkdir(parents=True, exist_ok=True)

    n_seeds = args.seeds or params["simulation"]["n_seeds"]
    print(f"  Sites: {len(params['sites'])}")
    print(f"  Seeds: {n_seeds}")
    print(f"  Output: {output_dir}")
    if calibrated is None:
        print("  ⚠️  PLACEHOLDER MODE — using default parameters")
    print()

    # Run ensemble
    print("Running ensemble simulations...")
    results = run_ensemble(params, calibrated, n_seeds=n_seeds)

    # Summarize
    print("\nComputing ensemble statistics...")
    summary = summarize_ensemble(results)

    # Save results
    print("\nSaving results...")

    # Per-site full results (one file per site, seed=0 only for genotypes)
    for site_name, seed_results in results.items():
        site_file = output_dir / f"{site_name.lower().replace(' ', '_')}.json"
        # Save seed 0 with genotypes, others as summaries only
        save_data = {
            "site": site_name,
            "summary": summary[site_name],
            "seed_0_full": seed_results[0],
            "n_seeds": n_seeds,
            "timestamp": datetime.now().isoformat(),
            "is_placeholder": calibrated is None,
        }
        with open(site_file, "w") as f:
            json.dump(save_data, f, indent=2)

    # Ensemble summary
    summary_file = output_dir / "ensemble_summary.json"
    with open(summary_file, "w") as f:
        json.dump({
            "summary": summary,
            "n_seeds": n_seeds,
            "n_sites": len(params["sites"]),
            "timestamp": datetime.now().isoformat(),
            "is_placeholder": calibrated is None,
            "params_hash": "TODO",  # TODO: hash of params for reproducibility
        }, f, indent=2)
    print(f"  Saved: {summary_file}")

    # Figures
    print("\nGenerating figures...")
    fig_dir = output_dir / "figures"
    fig_dir.mkdir(exist_ok=True)
    plot_latitude_heatmap(summary, fig_dir)
    plot_diversity_gradient(summary, fig_dir)

    # Print summary table
    print("\n" + "=" * 80)
    print("SITE SUMMARY (ensemble means)")
    print("=" * 80)
    print(f"{'Site':<16} {'Lat':>5} {'N':>6} {'r̄':>7} {'t̄':>7} {'c̄':>7} "
          f"{'Hₑ':>6} {'Nₑ':>7}")
    print("-" * 80)
    for site in sorted(summary.keys(), key=lambda s: summary[s]["lat"], reverse=True):
        s = summary[site]
        print(f"{site:<16} {s['lat']:>5.1f} "
              f"{s['n_survivors']['mean']:>6.0f} "
              f"{s['mean_r']['mean']:>7.3f} "
              f"{s['mean_t']['mean']:>7.3f} "
              f"{s['mean_c']['mean']:>7.3f} "
              f"{s['He']['mean']:>6.3f} "
              f"{s['Ne']['mean']:>7.1f}")

    print("\n✅ Analysis 1 complete.")
    if calibrated is None:
        print("⚠️  Results are PLACEHOLDERS — rerun after ABC-SMC calibration.")


if __name__ == "__main__":
    main()
