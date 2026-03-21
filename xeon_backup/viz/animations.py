"""Diagnostic animations for SSWD-EvoEpi simulation results.

Three animations:
  1. Epidemic Spread — geographic network, node size=pop, color=prevalence
  2. Allele Frequency Evolution — 52 loci over time, per-node colored lines
  3. Population Pyramid — age/size distribution collapsing and rebuilding

Usage:
    python -m viz.animations results/5node_epidemic_20yr/
    python -m viz.animations results/5node_epidemic_20yr/ --only spread
    python -m viz.animations results/5node_epidemic_20yr/ --fps 4

All output saved as GIF to <result_dir>/animations/
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np

from viz.plot_utils import SimulationData, NODE_COLORS


# ═══════════════════════════════════════════════════════════════════════
# GEOGRAPHIC LAYOUT
# ═══════════════════════════════════════════════════════════════════════

# Approximate (lon, lat) — we plot lon on x, lat on y
NODE_GEO = {
    "Sitka, AK":              (-135.34, 57.06),
    "Howe Sound, BC":         (-123.25, 49.52),
    "San Juan Islands, WA":   (-123.02, 48.53),
    "Newport, OR":            (-124.05, 44.63),
    "Monterey, CA":           (-121.90, 36.62),
}

# Edges: fully connected, but we only draw adjacent coastal pairs
# (aesthetics — full mesh is too cluttered)
COASTAL_EDGES = [
    (0, 1),  # Sitka → Howe Sound
    (1, 2),  # Howe Sound → SJI
    (2, 3),  # SJI → Newport
    (3, 4),  # Newport → Monterey
]


def _node_positions(names: list[str]) -> np.ndarray:
    """Return (n_nodes, 2) array of [lon, lat]."""
    pos = np.zeros((len(names), 2))
    for i, name in enumerate(names):
        if name in NODE_GEO:
            pos[i] = NODE_GEO[name]
        else:
            # Fallback: spread vertically
            pos[i] = [-125 + i * 2, 40 + i * 4]
    return pos


# ═══════════════════════════════════════════════════════════════════════
# ANIMATION 1: Epidemic Spread Network
# ═══════════════════════════════════════════════════════════════════════

def animate_epidemic_spread(
    data: SimulationData,
    output_path: Path,
    fps: int = 2,
    dpi: int = 120,
) -> Path:
    """Animate epidemic spreading through the node network.

    Each frame = one year.
    - Node size ∝ population (shrinks during epidemic)
    - Node color = disease prevalence (green → yellow → red)
    - Edge thickness = proportional to larval flow (constant here; could be
      made dynamic with per-year dispersal data)
    - Text annotation shows N and mean resistance

    Args:
        data: Loaded simulation data.
        output_path: GIF output path.
        fps: Frames per second.
        dpi: Resolution.

    Returns:
        Path to saved GIF.
    """
    pos = _node_positions(data.node_names)
    n_nodes = data.n_nodes
    n_years = data.n_years

    # Prevalence colormap: green (0) → yellow (0.5) → red (1)
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "epidemic", ["#2ca02c", "#ffdd57", "#e41a1c"]
    )
    norm = mcolors.Normalize(vmin=0, vmax=1.0)

    # Max population for consistent node sizing
    max_pop = max(data.yearly_pop.max(), 1)
    min_marker = 80
    max_marker = 800

    fig, ax = plt.subplots(figsize=(9, 11))

    # Static elements: coastline approximation (simple line)
    coast_lon = [-136, -135, -131, -128, -126, -124, -123.5, -124, -124.5, -122, -121]
    coast_lat = [60, 57, 54, 51, 49, 48, 46, 44, 42, 38, 36]
    ax.plot(coast_lon, coast_lat, color="#aaddff", linewidth=3, alpha=0.4, zorder=0)

    # Pre-create artists
    edge_lines = []
    for i_from, i_to in COASTAL_EDGES:
        line, = ax.plot(
            [pos[i_from, 0], pos[i_to, 0]],
            [pos[i_from, 1], pos[i_to, 1]],
            color="#999999", linewidth=1.5, alpha=0.5, zorder=1,
        )
        edge_lines.append(line)

    scatter = ax.scatter(
        pos[:, 0], pos[:, 1],
        s=np.full(n_nodes, 400), c=["#2ca02c"] * n_nodes,
        edgecolors="black", linewidth=1.5, zorder=3,
    )

    labels = []
    for i in range(n_nodes):
        txt = ax.text(
            pos[i, 0] + 0.8, pos[i, 1] + 0.5,
            "", fontsize=7, fontweight="bold", zorder=4,
            bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.85, ec="none"),
        )
        labels.append(txt)

    title = ax.set_title("", fontsize=14, fontweight="bold")

    ax.set_xlim(-138, -119)
    ax.set_ylim(34, 60)
    ax.set_xlabel("Longitude", fontsize=10)
    ax.set_ylabel("Latitude", fontsize=10)
    ax.set_aspect(1.3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Legend
    legend_patches = [
        mpatches.Patch(color="#2ca02c", label="Healthy (prev=0)"),
        mpatches.Patch(color="#ffdd57", label="Moderate (prev≈0.5)"),
        mpatches.Patch(color="#e41a1c", label="Severe (prev≈1.0)"),
    ]
    ax.legend(handles=legend_patches, loc="lower left", fontsize=8,
              framealpha=0.9, title="Disease Prevalence")

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.3, pad=0.02, label="Disease Prevalence")

    def update(year: int):
        pops = data.yearly_pop[:, year]
        dis_d = data.yearly_disease_deaths[:, year]
        nat_d = data.yearly_natural_deaths[:, year]
        total_d = dis_d + nat_d

        # Disease mortality fraction as prevalence proxy
        # For pre-disease years, this is 0; during epidemic, high
        dis_frac = np.where(total_d > 0, dis_d / np.maximum(total_d, 1), 0.0)

        # Node sizes (proportional to pop, with min)
        sizes = np.clip(pops / max_pop * max_marker, min_marker, max_marker)
        scatter.set_sizes(sizes)

        # Node colors
        colors = [cmap(norm(f)) for f in dis_frac]
        scatter.set_facecolors(colors)

        # Labels
        for i in range(n_nodes):
            name_short = data.node_names[i].split(",")[0]
            r_mean = data.yearly_mean_resistance[i, year]
            labels[i].set_text(
                f"{name_short}\nN={pops[i]}  r̄={r_mean:.3f}"
            )

        # Edge thickness: thicker during epidemic years (larval flow implied)
        for idx, (i_from, i_to) in enumerate(COASTAL_EDGES):
            base_width = 1.5
            # During epidemic years, edges get thinner (reduced connectivity)
            if year >= data.disease_year and (pops[i_from] < 50 or pops[i_to] < 50):
                edge_lines[idx].set_linewidth(0.5)
                edge_lines[idx].set_alpha(0.2)
            else:
                edge_lines[idx].set_linewidth(base_width)
                edge_lines[idx].set_alpha(0.5)

        phase = "Pre-epidemic" if year < data.disease_year else "Epidemic"
        title.set_text(
            f"Epidemic Spread — Year {year} ({phase})\n"
            f"Total Pop: {pops.sum():,}"
        )

        return [scatter] + labels + edge_lines + [title]

    anim = FuncAnimation(
        fig, update, frames=n_years, interval=1000 // fps, blit=False,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    anim.save(str(output_path), writer=PillowWriter(fps=fps), dpi=dpi)
    plt.close(fig)
    print(f"  ✓ Epidemic spread animation: {output_path}")
    return output_path


# ═══════════════════════════════════════════════════════════════════════
# ANIMATION 2: Allele Frequency Evolution
# ═══════════════════════════════════════════════════════════════════════

def animate_allele_evolution(
    data: SimulationData,
    output_path: Path,
    fps: int = 2,
    dpi: int = 120,
) -> Path:
    """Animate allele frequency changes across 52 loci over time.

    X-axis = 52 loci (0-50 additive + 51 EF1A)
    Y-axis = allele frequency
    Per-node colored lines.
    Vertical dashed line at epidemic start.
    Watch resistance alleles rise during/after epidemic.

    Uses pre/post epidemic snapshots + interpolation for intermediate years
    (we only have full 52-locus snapshots at two timepoints, plus top-3 yearly).

    For the full animation, we show:
    - Pre-epidemic snapshot for years < disease_year
    - Linear interpolation from pre→post for years disease_year → disease_year+2
    - Post-epidemic snapshot for years > disease_year+2

    Args:
        data: Loaded simulation data.
        output_path: GIF output path.
        fps: Frames per second.
        dpi: Resolution.

    Returns:
        Path to saved GIF.
    """
    n_nodes = data.n_nodes
    n_years = data.n_years
    n_loci = data.pre_epidemic_allele_freq.shape[1]  # 52
    disease_yr = data.disease_year

    fig, ax = plt.subplots(figsize=(14, 6))

    loci_x = np.arange(n_loci)

    # Create line artists per node
    lines = []
    for i in range(n_nodes):
        line, = ax.plot(
            loci_x, np.zeros(n_loci),
            color=NODE_COLORS[i], linewidth=1.8, alpha=0.8,
            label=data.node_names[i],
            marker=".", markersize=3,
        )
        lines.append(line)

    # EF1A marker
    ef1a_line = ax.axvline(51, color="gold", linewidth=2, linestyle="-", alpha=0.7)
    ax.text(51, 0.97, "EF1A", fontsize=8, color="goldenrod", ha="center",
            fontweight="bold", transform=ax.get_xaxis_transform())

    # Epidemic onset marker (will be toggled)
    onset_text = ax.text(
        0.02, 0.95, "", transform=ax.transAxes,
        fontsize=11, fontweight="bold", color="#e41a1c",
        bbox=dict(boxstyle="round", fc="white", ec="#e41a1c", alpha=0.9),
    )

    title = ax.set_title("", fontsize=13, fontweight="bold")
    ax.set_xlabel("Locus Index (0-50: additive, 51: EF1A)", fontsize=10)
    ax.set_ylabel("Allele Frequency", fontsize=10)
    ax.set_xlim(-1, n_loci)
    ax.set_ylim(0, 1.0)
    ax.legend(loc="upper right", fontsize=8, framealpha=0.9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", alpha=0.3)

    # Pre-compute interpolated allele freq snapshots per year per node
    # We have: pre_epidemic (at disease_year), post_epidemic (at disease_year+2)
    # For years in between, linearly interpolate
    # For other years, use the nearest snapshot
    def get_freq_snapshot(node_id: int, year: int) -> np.ndarray:
        pre = data.pre_epidemic_allele_freq[node_id]  # (52,)
        post = data.post_epidemic_allele_freq[node_id]  # (52,)

        if year <= disease_yr:
            return pre.copy()
        elif year >= disease_yr + 2:
            # After post-epidemic, apply gentle drift (random walk seeded)
            rng = np.random.RandomState(data.seed + node_id * 100 + year)
            drift_years = year - (disease_yr + 2)
            # Small drift: σ per year ≈ 0.003 (SRS-like)
            drift = rng.normal(0, 0.003, n_loci) * np.sqrt(drift_years)
            freq = np.clip(post + drift, 0.01, 0.99)
            return freq
        else:
            # Linear interpolation
            t = (year - disease_yr) / 2.0
            return pre * (1 - t) + post * t

    def update(year: int):
        for i in range(n_nodes):
            freq = get_freq_snapshot(i, year)
            lines[i].set_ydata(freq)

        phase = "Pre-epidemic" if year < disease_yr else "Epidemic"
        title.set_text(
            f"Allele Frequency Evolution — Year {year} ({phase})"
        )

        if year == disease_yr:
            onset_text.set_text("!! DISEASE INTRODUCED")
        elif year == disease_yr + 1:
            onset_text.set_text(">> Epidemic wave")
        elif year == disease_yr + 2:
            onset_text.set_text(">> Post-epidemic snapshot")
        elif year > disease_yr + 2:
            onset_text.set_text(f">> Drift + selection (yr {year - disease_yr} post)")
        else:
            onset_text.set_text("")

        return lines + [title, onset_text]

    anim = FuncAnimation(
        fig, update, frames=n_years, interval=1000 // fps, blit=False,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    anim.save(str(output_path), writer=PillowWriter(fps=fps), dpi=dpi)
    plt.close(fig)
    print(f"  ✓ Allele evolution animation: {output_path}")
    return output_path


# ═══════════════════════════════════════════════════════════════════════
# ANIMATION 3: Population Pyramid Over Time
# ═══════════════════════════════════════════════════════════════════════

def animate_population_pyramid(
    data: SimulationData,
    output_path: Path,
    fps: int = 2,
    dpi: int = 120,
) -> Path:
    """Animate population structure (age/stage distribution) over time.

    Shows stacked bar chart of population structure across all nodes:
    - Recruits (age 0)
    - Juveniles (age 1-2)
    - Subadults (age 3-4)
    - Adults (age 5+)

    We approximate the stage distribution from available data:
    - adults from yearly_adults
    - recruits from yearly_recruits
    - non-adults = pop - adults (includes recruits from previous year, juveniles, subadults)

    Color-coded by disease state proxy (disease deaths / pop that year).

    Args:
        data: Loaded simulation data.
        output_path: GIF output path.
        fps: Frames per second.
        dpi: Resolution.

    Returns:
        Path to saved GIF.
    """
    n_nodes = data.n_nodes
    n_years = data.n_years

    fig, axes = plt.subplots(1, n_nodes, figsize=(3.2 * n_nodes, 7), sharey=True)
    if n_nodes == 1:
        axes = [axes]

    fig.suptitle("", fontsize=14, fontweight="bold", y=0.98)

    # Stage categories
    categories = ["Recruits", "Juveniles\n+ Subadults", "Adults"]
    cat_colors_healthy = ["#66c2a5", "#fc8d62", "#8da0cb"]
    cat_colors_diseased = ["#a6d854", "#e78ac3", "#e41a1c"]
    y_positions = np.arange(len(categories))

    # Pre-create bar containers (one set per node)
    bars_list = []
    label_texts = []
    for idx, ax in enumerate(axes):
        bars = ax.barh(y_positions, [0, 0, 0], height=0.6,
                       color=cat_colors_healthy, edgecolor="white", linewidth=0.5)
        bars_list.append(bars)

        ax.set_yticks(y_positions)
        ax.set_yticklabels(categories if idx == 0 else [], fontsize=9)
        short_name = data.node_names[idx].split(",")[0]
        ax.set_xlabel(short_name, fontsize=9, fontweight="bold")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        # N text in corner
        txt = ax.text(0.95, 0.95, "", transform=ax.transAxes,
                      fontsize=8, ha="right", va="top",
                      bbox=dict(boxstyle="round", fc="white", alpha=0.8))
        label_texts.append(txt)

    # Set consistent x-axis limit
    max_K = int(data.node_K.max())
    for ax in axes:
        ax.set_xlim(0, max_K * 1.1)

    def update(year: int):
        for i in range(n_nodes):
            pop = int(data.yearly_pop[i, year])
            adults = int(data.yearly_adults[i, year])
            recruits = int(data.yearly_recruits[i, year])
            non_adult = max(0, pop - adults)  # juveniles + subadults (includes prev-year recruits)

            # Split non_adult into recruits (current year) and juveniles+subadults
            juv_sub = max(0, non_adult - recruits)

            values = [recruits, juv_sub, adults]

            # Disease intensity coloring
            dis_d = data.yearly_disease_deaths[i, year]
            dis_frac = dis_d / max(pop, 1)

            for j, bar in enumerate(bars_list[i]):
                bar.set_width(values[j])
                # Blend color: healthy → diseased based on disease fraction
                h_color = np.array(mcolors.to_rgba(cat_colors_healthy[j]))
                d_color = np.array(mcolors.to_rgba(cat_colors_diseased[j]))
                t = min(max(dis_frac, 0.0), 1.0)
                blended = np.clip(h_color * (1 - t) + d_color * t, 0.0, 1.0)
                bar.set_facecolor(tuple(blended))

            label_texts[i].set_text(f"N={pop}")

        phase = "Pre-epidemic" if year < data.disease_year else "Epidemic"
        total = int(data.yearly_total_pop[year])
        fig.suptitle(
            f"Population Structure — Year {year} ({phase})  |  Total: {total:,}",
            fontsize=13, fontweight="bold",
        )

        return [b for bars in bars_list for b in bars] + label_texts

    anim = FuncAnimation(
        fig, update, frames=n_years, interval=1000 // fps, blit=False,
    )

    fig.tight_layout(rect=[0, 0, 1, 0.95])

    output_path.parent.mkdir(parents=True, exist_ok=True)
    anim.save(str(output_path), writer=PillowWriter(fps=fps), dpi=dpi)
    plt.close(fig)
    print(f"  ✓ Population pyramid animation: {output_path}")
    return output_path


# ═══════════════════════════════════════════════════════════════════════
# GENERATE ALL
# ═══════════════════════════════════════════════════════════════════════

def generate_all_animations(
    result_dir: str | Path,
    output_dir: Optional[str | Path] = None,
    fps: int = 2,
    dpi: int = 120,
    only: Optional[str] = None,
) -> list[Path]:
    """Generate all animations from a simulation result directory.

    Args:
        result_dir: Directory with simulation_data.npz + metadata.json.
        output_dir: Output directory (default: <result_dir>/animations/).
        fps: Frames per second for GIFs.
        dpi: Resolution.
        only: Generate only one animation ("spread", "allele", "pyramid").

    Returns:
        List of generated file paths.
    """
    result_dir = Path(result_dir)
    if output_dir is None:
        output_dir = result_dir / "animations"
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading simulation data from {result_dir}...")
    data = SimulationData.load(result_dir)
    print(f"  {data.n_nodes} nodes, {data.n_years} years, disease at year {data.disease_year}")

    outputs = []

    if only is None or only == "spread":
        print("\n🗺️  Animation 1: Epidemic Spread Network...")
        path = animate_epidemic_spread(
            data, output_dir / "epidemic_spread.gif", fps=fps, dpi=dpi
        )
        outputs.append(path)

    if only is None or only == "allele":
        print("\n🧬 Animation 2: Allele Frequency Evolution...")
        path = animate_allele_evolution(
            data, output_dir / "allele_evolution.gif", fps=fps, dpi=dpi
        )
        outputs.append(path)

    if only is None or only == "pyramid":
        print("\n📊 Animation 3: Population Pyramid...")
        path = animate_population_pyramid(
            data, output_dir / "population_pyramid.gif", fps=fps, dpi=dpi
        )
        outputs.append(path)

    print(f"\n✅ {len(outputs)} animations generated in {output_dir}/")
    return outputs


# ═══════════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Generate diagnostic animations from SSWD-EvoEpi simulation results.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python -m viz.animations results/5node_epidemic_20yr/
  python -m viz.animations results/5node_epidemic_20yr/ --only spread --fps 4
  python -m viz.animations results/5node_epidemic_20yr/ --output-dir my_animations/
        """,
    )
    parser.add_argument("result_dir", help="Directory with simulation_data.npz + metadata.json")
    parser.add_argument("--output-dir", "-o", default=None,
                        help="Output directory (default: <result_dir>/animations/)")
    parser.add_argument("--fps", type=int, default=2, help="Frames per second (default: 2)")
    parser.add_argument("--dpi", type=int, default=120, help="Resolution (default: 120)")
    parser.add_argument("--only", choices=["spread", "allele", "pyramid"],
                        help="Generate only one animation")

    args = parser.parse_args()
    result_dir = Path(args.result_dir)
    if not result_dir.exists():
        print(f"Error: {result_dir} does not exist", file=sys.stderr)
        sys.exit(1)

    generate_all_animations(result_dir, args.output_dir, args.fps, args.dpi, args.only)


if __name__ == "__main__":
    main()
