"""Reusable matplotlib plotting functions for SSWD-EvoEpi.

Each function takes loaded simulation data (from NPZ + metadata) and returns
a matplotlib Figure. Designed for static PNG export and embedding in reports.

Usage:
    data = SimulationData.load("results/5node_epidemic_20yr/")
    fig = plot_population_trajectory(data, node_id=2)
    fig.savefig("pop_sji.png", dpi=150, bbox_inches="tight")
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

import matplotlib
matplotlib.use("Agg")  # Non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np


# ═══════════════════════════════════════════════════════════════════════
# DATA LOADER
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class SimulationData:
    """Container for loaded simulation output."""

    # Metadata
    n_years: int
    n_nodes: int
    node_names: List[str]
    disease_year: int
    seed: int
    initial_total_pop: int
    final_total_pop: int

    # Per-node yearly arrays (n_nodes, n_years)
    yearly_pop: np.ndarray
    yearly_adults: np.ndarray
    yearly_recruits: np.ndarray
    yearly_natural_deaths: np.ndarray
    yearly_disease_deaths: np.ndarray
    yearly_mean_resistance: np.ndarray
    yearly_vibrio_max: np.ndarray
    yearly_ef1a_freq: np.ndarray
    yearly_va: np.ndarray
    yearly_ne_ratio: np.ndarray

    # Per-node yearly allele freqs at top 3 loci (n_nodes, n_years, 3)
    yearly_allele_freq_top3: np.ndarray

    # Totals (n_years,)
    yearly_total_pop: np.ndarray
    yearly_total_larvae: np.ndarray

    # Snapshots
    peak_disease_prevalence: np.ndarray   # (n_nodes,)
    pre_epidemic_allele_freq: np.ndarray  # (n_nodes, n_loci)
    post_epidemic_allele_freq: np.ndarray # (n_nodes, n_loci)
    node_K: np.ndarray                    # (n_nodes,)

    @classmethod
    def load(cls, result_dir: str | Path) -> "SimulationData":
        """Load simulation results from a directory."""
        result_dir = Path(result_dir)

        with open(result_dir / "metadata.json") as f:
            meta = json.load(f)

        d = np.load(result_dir / "simulation_data.npz", allow_pickle=True)

        return cls(
            n_years=meta["n_years"],
            n_nodes=meta["n_nodes"],
            node_names=meta["node_names"],
            disease_year=meta["disease_year"],
            seed=meta["seed"],
            initial_total_pop=meta["initial_total_pop"],
            final_total_pop=meta["final_total_pop"],
            yearly_pop=d["yearly_pop"],
            yearly_adults=d["yearly_adults"],
            yearly_recruits=d["yearly_recruits"],
            yearly_natural_deaths=d["yearly_natural_deaths"],
            yearly_disease_deaths=d["yearly_disease_deaths"],
            yearly_mean_resistance=d["yearly_mean_resistance"],
            yearly_vibrio_max=d["yearly_vibrio_max"],
            yearly_ef1a_freq=d["yearly_ef1a_freq"],
            yearly_va=d["yearly_va"],
            yearly_ne_ratio=d["yearly_ne_ratio"],
            yearly_allele_freq_top3=d["yearly_allele_freq_top3"],
            yearly_total_pop=d["yearly_total_pop"],
            yearly_total_larvae=d["yearly_total_larvae_dispersed"],
            peak_disease_prevalence=d["peak_disease_prevalence"],
            pre_epidemic_allele_freq=d["pre_epidemic_allele_freq"],
            post_epidemic_allele_freq=d["post_epidemic_allele_freq"],
            node_K=d["node_K"],
        )

    @property
    def years(self) -> np.ndarray:
        """Year indices as array."""
        return np.arange(self.n_years)


# ═══════════════════════════════════════════════════════════════════════
# STYLE CONSTANTS
# ═══════════════════════════════════════════════════════════════════════

# Node colors (consistent across all plots)
NODE_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]
# Compartment/stage colors
STAGE_COLORS = {"recruits": "#66c2a5", "non_adult": "#fc8d62", "adults": "#8da0cb"}
DISEASE_COLORS = {"natural": "#888888", "disease": "#e41a1c"}


def _style_axis(ax, xlabel: str = "", ylabel: str = "", title: str = ""):
    """Apply consistent styling to an axis."""
    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)
    if title:
        ax.set_title(title, fontsize=12, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=9)


def _add_disease_vline(ax, disease_year: int, label: bool = True):
    """Add a vertical line at disease introduction."""
    ax.axvline(disease_year, color="#e41a1c", linewidth=1.5, linestyle="--",
               alpha=0.7, zorder=5)
    if label:
        ax.text(disease_year + 0.2, ax.get_ylim()[1] * 0.95, "Disease →",
                fontsize=8, color="#e41a1c", va="top")


# ═══════════════════════════════════════════════════════════════════════
# PLOT FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════

def plot_population_trajectory(data: SimulationData, node_id: int) -> plt.Figure:
    """Population trajectory for a single node.

    Shows total population as area, adults overlaid, carrying capacity K as
    dashed line, and recruits as bars at the bottom.

    Args:
        data: Loaded simulation data.
        node_id: Node index (0-based).

    Returns:
        matplotlib Figure.
    """
    fig, ax = plt.subplots(figsize=(10, 4))
    years = data.years
    name = data.node_names[node_id]
    K = data.node_K[node_id]

    pop = data.yearly_pop[node_id]
    adults = data.yearly_adults[node_id]
    recruits = data.yearly_recruits[node_id]
    non_adult = pop - adults

    # Stacked area: non-adults (bottom) + adults (top)
    ax.fill_between(years, 0, non_adult, alpha=0.4,
                    color=STAGE_COLORS["non_adult"], label="Juveniles + Subadults")
    ax.fill_between(years, non_adult, pop, alpha=0.6,
                    color=STAGE_COLORS["adults"], label="Adults")

    # Recruits as thin bars
    ax.bar(years, recruits, width=0.4, color=STAGE_COLORS["recruits"],
           alpha=0.7, label="Annual Recruits", zorder=3)

    # K line
    ax.axhline(K, color="black", linestyle=":", linewidth=1, alpha=0.5, label=f"K = {K}")

    # Disease line
    _add_disease_vline(ax, data.disease_year)

    _style_axis(ax, "Year", "Count", f"Population — {name}")
    ax.legend(loc="upper right", fontsize=8, framealpha=0.9)
    ax.set_xlim(0, data.n_years - 1)
    ax.set_ylim(0, max(K * 1.1, pop.max() * 1.1))

    fig.tight_layout()
    return fig


def plot_population_all_nodes(data: SimulationData) -> plt.Figure:
    """Population trajectories for all nodes in subplots."""
    fig, axes = plt.subplots(data.n_nodes, 1, figsize=(12, 3 * data.n_nodes),
                             sharex=True)
    if data.n_nodes == 1:
        axes = [axes]

    years = data.years
    for i, ax in enumerate(axes):
        name = data.node_names[i]
        K = data.node_K[i]
        pop = data.yearly_pop[i]
        adults = data.yearly_adults[i]
        non_adult = pop - adults
        recruits = data.yearly_recruits[i]

        ax.fill_between(years, 0, non_adult, alpha=0.4,
                        color=STAGE_COLORS["non_adult"])
        ax.fill_between(years, non_adult, pop, alpha=0.6,
                        color=STAGE_COLORS["adults"])
        ax.bar(years, recruits, width=0.4, color=STAGE_COLORS["recruits"],
               alpha=0.7, zorder=3)
        ax.axhline(K, color="black", linestyle=":", linewidth=1, alpha=0.5)
        ax.axvline(data.disease_year, color="#e41a1c", linewidth=1.5,
                   linestyle="--", alpha=0.7)

        ax.set_ylabel("Count", fontsize=9)
        ax.set_title(f"{name} (K={K})", fontsize=10, fontweight="bold")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_ylim(0, K * 1.15)

    axes[-1].set_xlabel("Year", fontsize=10)
    axes[0].legend(["Juv+Sub", "Adults", "Recruits", "K", "Disease"],
                   loc="upper right", fontsize=7, ncol=5)

    fig.suptitle("Population Dynamics by Node", fontsize=14, fontweight="bold", y=1.01)
    fig.tight_layout()
    return fig


def plot_disease_trajectory(data: SimulationData, node_id: int) -> plt.Figure:
    """Disease impact trajectory for a single node.

    Shows disease deaths vs natural deaths (stacked bars), Vibrio concentration
    (line on secondary axis).

    Args:
        data: Loaded simulation data.
        node_id: Node index (0-based).

    Returns:
        matplotlib Figure.
    """
    fig, ax1 = plt.subplots(figsize=(10, 4))
    years = data.years
    name = data.node_names[node_id]

    nat_d = data.yearly_natural_deaths[node_id]
    dis_d = data.yearly_disease_deaths[node_id]

    # Stacked bars: natural (bottom) + disease (top)
    ax1.bar(years, nat_d, width=0.6, color=DISEASE_COLORS["natural"],
            alpha=0.7, label="Natural Deaths")
    ax1.bar(years, dis_d, width=0.6, bottom=nat_d,
            color=DISEASE_COLORS["disease"], alpha=0.8, label="Disease Deaths")

    _style_axis(ax1, "Year", "Deaths/Year", f"Disease Impact — {name}")
    ax1.legend(loc="upper left", fontsize=8)

    # Secondary axis: Vibrio
    ax2 = ax1.twinx()
    vibrio = data.yearly_vibrio_max[node_id]
    ax2.plot(years, vibrio, color="#ff7f0e", linewidth=2, alpha=0.8,
             label="Peak Vibrio (bact/mL)")
    ax2.set_ylabel("Peak Vibrio (bact/mL)", fontsize=10, color="#ff7f0e")
    ax2.tick_params(axis="y", labelcolor="#ff7f0e", labelsize=9)
    ax2.spines["top"].set_visible(False)

    # Disease line
    ax1.axvline(data.disease_year, color="#e41a1c", linewidth=1.5,
                linestyle="--", alpha=0.7)
    ax1.set_xlim(-0.5, data.n_years - 0.5)

    fig.tight_layout()
    return fig


def plot_disease_all_nodes(data: SimulationData) -> plt.Figure:
    """Disease impact for all nodes."""
    fig, axes = plt.subplots(data.n_nodes, 1, figsize=(12, 3 * data.n_nodes),
                             sharex=True)
    if data.n_nodes == 1:
        axes = [axes]

    years = data.years
    for i, ax in enumerate(axes):
        nat_d = data.yearly_natural_deaths[i]
        dis_d = data.yearly_disease_deaths[i]

        ax.bar(years, nat_d, width=0.6, color=DISEASE_COLORS["natural"], alpha=0.7)
        ax.bar(years, dis_d, width=0.6, bottom=nat_d,
               color=DISEASE_COLORS["disease"], alpha=0.8)
        ax.axvline(data.disease_year, color="#e41a1c", linewidth=1.5,
                   linestyle="--", alpha=0.7)

        ax.set_ylabel("Deaths", fontsize=9)
        ax.set_title(f"{data.node_names[i]}", fontsize=10, fontweight="bold")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    axes[-1].set_xlabel("Year", fontsize=10)
    axes[0].legend(["Natural", "Disease"], loc="upper right", fontsize=8)

    fig.suptitle("Disease Mortality by Node", fontsize=14, fontweight="bold", y=1.01)
    fig.tight_layout()
    return fig


def plot_allele_frequencies(data: SimulationData, node_id: int,
                            loci: Optional[List[int]] = None) -> plt.Figure:
    """Allele frequency trajectories for top loci.

    Args:
        data: Loaded simulation data.
        node_id: Node index.
        loci: Which loci indices to plot (default: top 3 stored in data).

    Returns:
        matplotlib Figure.
    """
    fig, ax = plt.subplots(figsize=(10, 4))
    years = data.years
    name = data.node_names[node_id]

    # Top 3 loci from yearly tracking
    freq = data.yearly_allele_freq_top3[node_id]  # (n_years, 3)
    labels = ["Locus 1 (largest effect)", "Locus 2", "Locus 3"]
    colors = ["#e41a1c", "#377eb8", "#4daf4a"]

    for j in range(3):
        ax.plot(years, freq[:, j], color=colors[j], linewidth=2,
                marker="o", markersize=3, label=labels[j])

    _add_disease_vline(ax, data.disease_year, label=False)
    _style_axis(ax, "Year", "Allele Frequency",
                f"Top-3 Locus Allele Frequencies — {name}")
    ax.legend(loc="best", fontsize=8)
    ax.set_xlim(0, data.n_years - 1)
    ax.set_ylim(0, 1)

    fig.tight_layout()
    return fig


def plot_network_state(data: SimulationData, year: int) -> plt.Figure:
    """Network overview at a specific year.

    Nodes positioned by approximate lat/lon. Size = population, color = disease
    mortality (that year), edges indicate connectivity.

    Args:
        data: Loaded simulation data.
        year: Year index to display.

    Returns:
        matplotlib Figure.
    """
    # Approximate positions for the 5-node test network
    node_positions = {
        "Sitka, AK": (-135.34, 57.06),
        "Howe Sound, BC": (-123.25, 49.52),
        "San Juan Islands, WA": (-123.02, 48.53),
        "Newport, OR": (-124.05, 44.63),
        "Monterey, CA": (-121.90, 36.62),
    }

    fig, ax = plt.subplots(figsize=(8, 10))

    pops = data.yearly_pop[:, year]
    dis_d = data.yearly_disease_deaths[:, year]
    total_d = dis_d + data.yearly_natural_deaths[:, year]

    # Disease mortality fraction (0 = no disease, 1 = all disease)
    dis_frac = np.where(total_d > 0, dis_d / total_d, 0.0)

    # Colormap: green (healthy) → red (diseased)
    cmap = plt.cm.RdYlGn_r
    norm = mcolors.Normalize(vmin=0, vmax=1)

    # Draw connectivity edges (fully connected, width ~ inverse distance)
    for i in range(data.n_nodes):
        for j in range(i + 1, data.n_nodes):
            ni = data.node_names[i]
            nj = data.node_names[j]
            if ni in node_positions and nj in node_positions:
                xi, yi = node_positions[ni]
                xj, yj = node_positions[nj]
                ax.plot([xi, xj], [yi, yj], color="#cccccc", linewidth=0.8,
                        alpha=0.5, zorder=1)

    # Draw nodes
    for i in range(data.n_nodes):
        name = data.node_names[i]
        if name not in node_positions:
            continue
        x, y = node_positions[name]
        size = max(30, pops[i] / data.node_K[i] * 500)
        color = cmap(norm(dis_frac[i]))
        ax.scatter(x, y, s=size, c=[color], edgecolors="black", linewidth=1.5,
                   zorder=3)
        ax.annotate(f"{name}\nN={pops[i]}", (x, y),
                    textcoords="offset points", xytext=(12, 5),
                    fontsize=8, fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8))

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.5, label="Disease Mortality Fraction")

    ax.set_xlabel("Longitude", fontsize=10)
    ax.set_ylabel("Latitude", fontsize=10)
    ax.set_title(f"Network State — Year {year}", fontsize=14, fontweight="bold")
    ax.set_aspect(1.3)

    fig.tight_layout()
    return fig


def plot_genetic_heatmap(data: SimulationData, node_id: int) -> plt.Figure:
    """Heatmap of allele frequency change across all 52 loci.

    X = locus (0–51, locus 51 = EF1A), color = Δq (post − pre epidemic).

    Args:
        data: Loaded simulation data.
        node_id: Node index.

    Returns:
        matplotlib Figure.
    """
    fig, ax = plt.subplots(figsize=(14, 3))
    name = data.node_names[node_id]

    pre = data.pre_epidemic_allele_freq[node_id]   # (52,)
    post = data.post_epidemic_allele_freq[node_id]  # (52,)
    delta_q = post - pre

    # Reshape for imshow: 1 row × 52 columns
    im = ax.imshow(delta_q.reshape(1, -1), aspect="auto",
                   cmap="RdBu_r", vmin=-0.15, vmax=0.15,
                   interpolation="nearest")

    ax.set_xlabel("Locus Index", fontsize=10)
    ax.set_yticks([])
    ax.set_title(f"Allele Frequency Shift (Δq) — {name}", fontsize=12, fontweight="bold")

    # Mark EF1A (locus 51)
    ax.axvline(51, color="gold", linewidth=2, linestyle="-", alpha=0.8)
    ax.text(51, 0.6, "EF1A", fontsize=8, color="gold", ha="center",
            fontweight="bold", transform=ax.get_xaxis_transform())

    plt.colorbar(im, ax=ax, shrink=0.8, label="Δq (post − pre epidemic)")

    fig.tight_layout()
    return fig


def plot_genetic_heatmap_all_nodes(data: SimulationData) -> plt.Figure:
    """Heatmap of Δq across all nodes and loci."""
    fig, axes = plt.subplots(data.n_nodes, 1, figsize=(14, 2.5 * data.n_nodes),
                             sharex=True)
    if data.n_nodes == 1:
        axes = [axes]

    for i, ax in enumerate(axes):
        pre = data.pre_epidemic_allele_freq[i]
        post = data.post_epidemic_allele_freq[i]
        delta_q = post - pre

        im = ax.imshow(delta_q.reshape(1, -1), aspect="auto",
                       cmap="RdBu_r", vmin=-0.15, vmax=0.15,
                       interpolation="nearest")
        ax.set_yticks([])
        ax.set_ylabel(data.node_names[i], fontsize=8, rotation=0,
                      ha="right", va="center")
        ax.axvline(51, color="gold", linewidth=2, alpha=0.8)

    axes[-1].set_xlabel("Locus Index", fontsize=10)
    fig.suptitle("Allele Frequency Shifts (Δq) — All Nodes",
                 fontsize=14, fontweight="bold")

    # Shared colorbar
    fig.subplots_adjust(right=0.88)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    fig.colorbar(im, cax=cbar_ax, label="Δq")

    return fig


def plot_resistance_distribution(data: SimulationData, node_id: int,
                                 year: int) -> plt.Figure:
    """Mean resistance over time with bar showing value at specified year.

    Note: We only have yearly mean resistance (not per-individual distributions)
    from the recorded data. Shows the trajectory with a marker at the target year.
    """
    fig, ax = plt.subplots(figsize=(10, 4))
    years = data.years
    name = data.node_names[node_id]

    r_mean = data.yearly_mean_resistance[node_id]

    ax.plot(years, r_mean, color=NODE_COLORS[node_id], linewidth=2.5,
            marker="o", markersize=4, label=f"Mean r̄")
    ax.axvline(year, color="purple", linewidth=2, linestyle=":", alpha=0.7)
    ax.plot(year, r_mean[year], "s", color="purple", markersize=12, zorder=5,
            label=f"Year {year}: r̄ = {r_mean[year]:.4f}")

    _add_disease_vline(ax, data.disease_year, label=False)
    _style_axis(ax, "Year", "Mean Resistance (r̄)",
                f"Resistance Trajectory — {name}")
    ax.legend(loc="best", fontsize=9)
    ax.set_xlim(0, data.n_years - 1)

    fig.tight_layout()
    return fig


def plot_genetics_summary(data: SimulationData) -> plt.Figure:
    """Multi-panel genetics summary: resistance, EF1A, Ne/N, Va for all nodes."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 8), sharex=True)
    years = data.years

    # Panel 1: Mean resistance
    ax = axes[0, 0]
    for i in range(data.n_nodes):
        ax.plot(years, data.yearly_mean_resistance[i], color=NODE_COLORS[i],
                linewidth=2, label=data.node_names[i])
    ax.axvline(data.disease_year, color="#e41a1c", linewidth=1.5,
               linestyle="--", alpha=0.7)
    _style_axis(ax, "", "Mean r̄", "Mean Resistance")
    ax.legend(fontsize=7, loc="upper left")

    # Panel 2: EF1A frequency
    ax = axes[0, 1]
    for i in range(data.n_nodes):
        ax.plot(years, data.yearly_ef1a_freq[i], color=NODE_COLORS[i],
                linewidth=2, label=data.node_names[i])
    ax.axvline(data.disease_year, color="#e41a1c", linewidth=1.5,
               linestyle="--", alpha=0.7)
    _style_axis(ax, "", "q (EF1A)", "EF1A Allele Frequency")

    # Panel 3: Ne/N ratio
    ax = axes[1, 0]
    for i in range(data.n_nodes):
        ax.plot(years, data.yearly_ne_ratio[i], color=NODE_COLORS[i],
                linewidth=2, label=data.node_names[i])
    ax.axvline(data.disease_year, color="#e41a1c", linewidth=1.5,
               linestyle="--", alpha=0.7)
    _style_axis(ax, "Year", "Ne / N", "Effective Population Ratio")

    # Panel 4: Additive genetic variance
    ax = axes[1, 1]
    for i in range(data.n_nodes):
        ax.plot(years, data.yearly_va[i], color=NODE_COLORS[i],
                linewidth=2, label=data.node_names[i])
    ax.axvline(data.disease_year, color="#e41a1c", linewidth=1.5,
               linestyle="--", alpha=0.7)
    _style_axis(ax, "Year", "V_A", "Additive Genetic Variance")

    fig.suptitle("Genetics Summary — All Nodes", fontsize=14, fontweight="bold")
    fig.tight_layout()
    return fig
