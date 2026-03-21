#!/usr/bin/env python3
"""Analysis 5: Integrated Recommendations.

Synthesizes results from Analyses 1-4 into actionable conservation
recommendations. Generates summary tables and publication figures.

Produces:
  - Summary recommendation tables (LaTeX + markdown)
  - "Money figures" for the paper:
    1. Latitude × trait heatmap (from Analysis 1)
    2. Screening effort curves with site overlay
    3. Breeding strategy trajectories
    4. Reintroduction scenario comparison
  - Timeline from program initiation to recovery targets
  - Key uncertainties summary

Usage:
    python 05_recommendations.py [--output-dir DIR]

Dependencies: Results from Analyses 1-4.
See: conservation report Section 9.5
"""

import sys
import json
import warnings
import argparse
from pathlib import Path
from datetime import datetime

import numpy as np
import yaml

# ── Project imports ──────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))


# ═══════════════════════════════════════════════════════════════════════
# LOAD RESULTS
# ═══════════════════════════════════════════════════════════════════════

def load_params(params_path: str = None) -> dict:
    if params_path is None:
        params_path = Path(__file__).parent / "params.yaml"
    with open(params_path) as f:
        return yaml.safe_load(f)


def load_all_results(params: dict) -> dict:
    """Load results from all analyses.

    Returns dict with keys: analysis1, analysis2, analysis3, analysis4.
    Missing analyses return None with a warning.
    """
    base = Path(__file__).parent / ".."
    results = {}

    for key, subdir in [
        ("analysis1", params["output"]["analysis1"]),
        ("analysis2", params["output"]["analysis2"]),
        ("analysis3", params["output"]["analysis3"]),
        ("analysis4", params["output"]["analysis4"]),
    ]:
        result_dir = base / subdir
        # Look for the main result file
        candidates = [
            result_dir / "ensemble_summary.json",
            result_dir / "screening_results.json",
            result_dir / "breeding_results.json",
            result_dir / "reintroduction_results.json",
        ]

        found = None
        for c in candidates:
            if c.exists():
                with open(c) as f:
                    found = json.load(f)
                break

        if found is None:
            warnings.warn(
                f"⚠️  {key} results not found in {result_dir}",
                UserWarning,
            )
        results[key] = found

    return results


# ═══════════════════════════════════════════════════════════════════════
# SUMMARY TABLES
# ═══════════════════════════════════════════════════════════════════════

def generate_site_summary_table(a1: dict) -> str:
    """Generate markdown summary table of site-level results."""
    if a1 is None:
        return "⚠️ Analysis 1 results not available.\n"

    summary = a1.get("summary", {})
    lines = [
        "## Site-Level Genetic State Summary (2026)",
        "",
        "| Site | Lat (°N) | N survivors | Mean r | Mean t | Mean c | Hₑ | Nₑ |",
        "|------|----------|-------------|--------|--------|--------|----|----|",
    ]

    for site in sorted(summary.keys(),
                       key=lambda s: summary[s]["lat"], reverse=True):
        s = summary[site]
        lines.append(
            f"| {site} | {s['lat']:.1f} | "
            f"{s['n_survivors']['mean']:.0f} ± {s['n_survivors']['std']:.0f} | "
            f"{s['mean_r']['mean']:.3f} | "
            f"{s['mean_t']['mean']:.3f} | "
            f"{s['mean_c']['mean']:.3f} | "
            f"{s['He']['mean']:.3f} | "
            f"{s['Ne']['mean']:.0f} |"
        )

    return "\n".join(lines)


def generate_screening_summary(a2: dict) -> str:
    """Generate markdown screening recommendation table."""
    if a2 is None:
        return "⚠️ Analysis 2 results not available.\n"

    alloc = a2.get("allocation", {})
    lines = [
        "## Screening Allocation Recommendations",
        "",
        f"Total budget: {alloc.get('total_budget', '?')} individuals",
        f"Best site: {alloc.get('best_site', '?')} "
        f"(E[max r] = {alloc.get('best_expected_max', 0):.3f})",
        "",
        "| Site | Allocation | Expected Max r |",
        "|------|------------|----------------|",
    ]

    site_names = alloc.get("site_names", [])
    allocations = alloc.get("allocations", [])
    expected = alloc.get("expected_maxes", [])

    for name, n, em in zip(site_names, allocations, expected):
        lines.append(f"| {name} | {n} | {em:.3f} |")

    return "\n".join(lines)


def generate_breeding_summary(a3: dict) -> str:
    """Generate markdown breeding strategy comparison."""
    if a3 is None:
        return "⚠️ Analysis 3 results not available.\n"

    summary = a3.get("summary", {})
    lines = [
        "## Breeding Strategy Comparison",
        "",
        "| Strategy | Δr̄ | Final Hₑ | Final F | Pareto |",
        "|----------|-----|----------|---------|--------|",
    ]

    pareto = a3.get("pareto", {})
    pareto_idx = set(pareto.get("pareto_optimal", []))

    for i, (strategy, data) in enumerate(summary.items()):
        ep = data.get("endpoint", {})
        is_p = "✓" if i in pareto_idx else ""
        lines.append(
            f"| {strategy} | "
            f"{ep.get('delta_r_mean', 0):.4f} ± {ep.get('delta_r_std', 0):.4f} | "
            f"{ep.get('final_he_mean', 0):.4f} | "
            f"{ep.get('final_f_mean', 0):.4f} | "
            f"{is_p} |"
        )

    return "\n".join(lines)


def generate_reintroduction_summary(a4: dict) -> str:
    """Generate markdown reintroduction scenario comparison."""
    if a4 is None:
        return "⚠️ Analysis 4 results not available.\n"

    results_list = a4.get("results", [])
    if not results_list:
        return "⚠️ No reintroduction scenarios completed.\n"

    lines = [
        "## Reintroduction Scenario Summary",
        "",
        "| Scenario | Release Size | Nodes | Persistence | Final Pop |",
        "|----------|-------------|-------|-------------|-----------|",
    ]

    # Show top 10 by persistence
    sorted_results = sorted(
        results_list,
        key=lambda r: r["ensemble"]["persistence_probability"],
        reverse=True,
    )

    for r in sorted_results[:10]:
        cfg = r["config"]
        ens = r["ensemble"]
        lines.append(
            f"| {r['scenario'][:30]} | "
            f"{cfg['release_size']} | "
            f"{len(cfg['release_nodes'])} | "
            f"{ens['persistence_probability']:.2f} | "
            f"{ens['final_population_mean']:.0f} |"
        )

    if len(sorted_results) > 10:
        lines.append(f"\n*... and {len(sorted_results) - 10} more scenarios*")

    return "\n".join(lines)


# ═══════════════════════════════════════════════════════════════════════
# MONEY FIGURES
# ═══════════════════════════════════════════════════════════════════════

def generate_money_figures(results: dict, output_dir: Path):
    """Generate the key publication figures.

    These are the figures that tell the story:
    1. Latitudinal gradient in post-epidemic resistance
    2. Screening effort comparison across sites
    3. Breeding trajectories (best vs baseline)
    4. Reintroduction persistence by strategy
    """
    import matplotlib.pyplot as plt

    fig_dir = output_dir / "figures"
    fig_dir.mkdir(exist_ok=True)

    a1 = results.get("analysis1")
    a2 = results.get("analysis2")
    a3 = results.get("analysis3")

    # ── Figure 1: Integrated site × trait panel ──────────────────────
    if a1 is not None:
        summary = a1.get("summary", {})
        sites = sorted(summary.keys(),
                       key=lambda s: summary[s]["lat"], reverse=True)

        fig, axes = plt.subplots(1, 3, figsize=(14, 8), sharey=True)
        traits = [("mean_r", "Resistance", "#2196F3"),
                  ("mean_t", "Tolerance", "#FF9800"),
                  ("mean_c", "Recovery", "#4CAF50")]

        for ax, (trait_key, label, color) in zip(axes, traits):
            means = [summary[s][trait_key]["mean"] for s in sites]
            stds = [summary[s][trait_key]["std"] for s in sites]
            lats = [summary[s]["lat"] for s in sites]

            ax.barh(range(len(sites)), means, xerr=stds,
                    color=color, alpha=0.7, edgecolor="white")
            ax.set_yticks(range(len(sites)))
            ax.set_yticklabels([f"{s} ({lats[i]:.0f}°N)"
                               for i, s in enumerate(sites)],
                              fontsize=9)
            ax.set_xlabel(f"Mean {label}", fontsize=11)
            ax.set_title(label, fontsize=12, fontweight="bold")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

        fig.suptitle("Post-Epidemic Trait Values by Site (2026)",
                     fontsize=14, fontweight="bold")
        fig.tight_layout()
        fig.savefig(fig_dir / "fig1_site_traits.png", dpi=300,
                    bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: fig1_site_traits.png")

    # ── Figure 2: Screening effort overlay ───────────────────────────
    if a2 is not None:
        screening = a2.get("screening", {})
        if screening:
            fig, ax = plt.subplots(figsize=(10, 6))

            colors = plt.cm.viridis(np.linspace(0, 1, len(screening)))
            for (site_name, data), color in zip(
                sorted(screening.items()), colors
            ):
                effort = data.get("effort", {})
                if effort:
                    ax.plot(effort["sample_sizes"],
                            effort["expected_maxes"],
                            "o-", color=color, linewidth=1.5,
                            markersize=3, label=site_name)

            ax.axhline(0.30, color="red", linestyle="--", alpha=0.5,
                       label="Target r = 0.30")
            ax.set_xscale("log")
            ax.set_xlabel("Sample Size", fontsize=11)
            ax.set_ylabel("Expected Best Resistance", fontsize=11)
            ax.set_title("Screening Effort: Diminishing Returns by Site",
                         fontsize=13, fontweight="bold")
            ax.legend(fontsize=8, ncol=2, frameon=False)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

            fig.tight_layout()
            fig.savefig(fig_dir / "fig2_screening_curves.png", dpi=300,
                        bbox_inches="tight")
            plt.close(fig)
            print(f"  Saved: fig2_screening_curves.png")

    # ── Figure 3: Breeding trajectories ──────────────────────────────
    if a3 is not None:
        summary = a3.get("summary", {})
        if summary:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

            colors = {"truncation": "#2196F3", "assortative": "#FF9800",
                      "complementary": "#4CAF50", "ocs": "#9C27B0"}

            for strategy, data in summary.items():
                gens = data["generations"]
                color = colors.get(strategy, "gray")

                ax1.plot(gens, data["mean_r"]["mean"],
                         "o-", color=color, linewidth=2,
                         markersize=4, label=strategy)
                ax1.fill_between(
                    gens,
                    np.array(data["mean_r"]["mean"]) - np.array(data["mean_r"]["std"]),
                    np.array(data["mean_r"]["mean"]) + np.array(data["mean_r"]["std"]),
                    alpha=0.15, color=color,
                )

                ax2.plot(gens, data["he"]["mean"],
                         "o-", color=color, linewidth=2,
                         markersize=4, label=strategy)

            ax1.set_xlabel("Generation", fontsize=11)
            ax1.set_ylabel("Mean Resistance", fontsize=11)
            ax1.set_title("Genetic Gain", fontsize=12, fontweight="bold")
            ax1.legend(fontsize=10, frameon=False)
            ax1.spines["top"].set_visible(False)
            ax1.spines["right"].set_visible(False)

            ax2.set_xlabel("Generation", fontsize=11)
            ax2.set_ylabel("Expected Heterozygosity (Hₑ)", fontsize=11)
            ax2.set_title("Diversity Retention", fontsize=12,
                          fontweight="bold")
            ax2.legend(fontsize=10, frameon=False)
            ax2.spines["top"].set_visible(False)
            ax2.spines["right"].set_visible(False)

            fig.suptitle("Breeding Strategy Comparison",
                         fontsize=14, fontweight="bold")
            fig.tight_layout()
            fig.savefig(fig_dir / "fig3_breeding_trajectories.png",
                        dpi=300, bbox_inches="tight")
            plt.close(fig)
            print(f"  Saved: fig3_breeding_trajectories.png")

    print("  Figures complete.")


# ═══════════════════════════════════════════════════════════════════════
# RECOMMENDATIONS
# ═══════════════════════════════════════════════════════════════════════

def generate_recommendations(results: dict) -> str:
    """Synthesize actionable recommendations.

    TODO: Update with actual insights after calibration.
    """
    lines = [
        "## Conservation Recommendations",
        "",
        "### 1. Sampling Protocol",
        "- **Priority sites**: Sites with highest predicted resistance (northern)",
        "- **Sample size**: Per screening analysis (Analysis 2)",
        "- **Genotyping**: Full 51-locus panel (Schiebelhut GWAS loci)",
        "- **Timeline**: Begin sampling ASAP — standing variation is eroding",
        "",
        "### 2. Breeding Program Design",
        "- **Strategy**: [TBD — from Pareto frontier analysis]",
        "- **Founders**: [TBD — from founder sensitivity analysis]",
        "- **Generations**: [TBD — from gain trajectory analysis]",
        "- **Selection target**: Resistance (primary), with diversity constraint",
        "",
        "### 3. Release Strategy",
        "- **Locations**: [TBD — from reintroduction scenario analysis]",
        "- **Size**: [TBD — critical threshold from scenario analysis]",
        "- **Timing**: Spring release (pre-disease season)",
        "- **Frequency**: [TBD — one-time vs. annual]",
        "",
        "### 4. Timeline",
        "- Year 0: Wild sampling + genotyping",
        "- Year 0-1: Founder selection + captive breeding setup",
        "- Year 1-N: Breeding program (N = generations to target)",
        "- Year N+1: First release",
        "- Year N+1 onwards: Monitoring + adaptive management",
        "",
        "### 5. Key Uncertainties",
        "- Calibrated disease parameters (ABC-SMC pending)",
        "- Inbreeding depression severity in Pycnopodia",
        "- Captive breeding feasibility at scale",
        "- Climate change effects on SST and disease dynamics",
        "- Multi-species reservoir dynamics",
        "",
        "### 6. Adaptive Management Decision Points",
        "- After genotyping: Are resistance alleles present at sufficient frequency?",
        "- After 3 breeding generations: Is gain on track?",
        "- After first release: Survival and disease exposure data",
        "- After 5 years: Population trajectory vs predictions",
    ]

    return "\n".join(lines)


# ═══════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Analysis 5: Integrated recommendations"
    )
    parser.add_argument("--output-dir", type=str, default=None)
    parser.add_argument("--params", type=str, default=None)
    args = parser.parse_args()

    print("=" * 60)
    print("Analysis 5: Integrated Recommendations")
    print("=" * 60)

    params = load_params(args.params)
    output_dir = Path(args.output_dir) if args.output_dir else \
        Path(__file__).parent / ".." / params["output"]["analysis5"]
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load all results
    print("\nLoading results from Analyses 1-4...")
    results = load_all_results(params)

    n_available = sum(1 for v in results.values() if v is not None)
    print(f"  Available: {n_available}/4 analyses")

    is_placeholder = any(
        (r or {}).get("is_placeholder", True) for r in results.values()
    )
    if is_placeholder:
        print("  ⚠️  Some results are PLACEHOLDERS")

    # Generate summary tables
    print("\nGenerating summary tables...")
    sections = [
        ("# Conservation Genetics: Integrated Analysis Report", ""),
        ("", f"*Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}*"),
        ("", ""),
        generate_site_summary_table(results["analysis1"]),
        ("", ""),
        generate_screening_summary(results["analysis2"]),
        ("", ""),
        generate_breeding_summary(results["analysis3"]),
        ("", ""),
        generate_reintroduction_summary(results["analysis4"]),
        ("", ""),
        generate_recommendations(results),
    ]

    # Write markdown report
    report_lines = []
    for section in sections:
        if isinstance(section, str):
            report_lines.append(section)
        elif isinstance(section, tuple):
            report_lines.extend(section)

    report_text = "\n".join(report_lines)
    report_path = output_dir / "recommendations_report.md"
    with open(report_path, "w") as f:
        f.write(report_text)
    print(f"  Saved: {report_path}")

    # Generate figures
    print("\nGenerating publication figures...")
    generate_money_figures(results, output_dir)

    # Save machine-readable summary
    summary = {
        "analysis_status": {
            "analysis1": results["analysis1"] is not None,
            "analysis2": results["analysis2"] is not None,
            "analysis3": results["analysis3"] is not None,
            "analysis4": results["analysis4"] is not None,
        },
        "is_placeholder": is_placeholder,
        "timestamp": datetime.now().isoformat(),
    }
    with open(output_dir / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("\n" + "=" * 60)
    print("REPORT SUMMARY")
    print("=" * 60)
    print(f"  Report: {report_path}")
    print(f"  Figures: {output_dir / 'figures'}")
    if is_placeholder:
        print("\n  ⚠️  PLACEHOLDER RESULTS — rerun after calibration")
    print("\n✅ Analysis 5 complete.")


if __name__ == "__main__":
    main()
