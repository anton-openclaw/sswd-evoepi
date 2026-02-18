"""Spatial & metapopulation visualizations for SSWD-EvoEpi.

Every function:
  - Accepts model results (SpatialSimResult), MetapopulationNetwork, or config
  - Returns a matplotlib Figure
  - Has an optional ``save_path`` parameter (saves PNG when given)
  - Uses the shared dark theme from ``sswd_evoepi.viz.style``

matplotlib backend is forced to Agg (no display) on import.

Geographic plots: x = longitude, y = latitude (actual positions from NodeDefinition).

Authors: Anton 🔬 & Willem Weertman
"""

from __future__ import annotations

import matplotlib
matplotlib.use('Agg')

from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, TYPE_CHECKING

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches
import numpy as np

from sswd_evoepi.viz.style import (
    ACCENT_COLORS,
    DARK_BG,
    DARK_PANEL,
    GRID_COLOR,
    NODE_COLORS,
    TEXT_COLOR,
    apply_dark_theme,
    dark_figure,
    save_figure,
)

if TYPE_CHECKING:
    from sswd_evoepi.model import SpatialSimResult
    from sswd_evoepi.spatial import MetapopulationNetwork, NodeDefinition


# ═══════════════════════════════════════════════════════════════════════
# HELPER: extract geo coords from network
# ═══════════════════════════════════════════════════════════════════════

def _get_coords(network: 'MetapopulationNetwork') -> Tuple[np.ndarray, np.ndarray]:
    """Extract lon/lat arrays from a network's node definitions."""
    lons = np.array([n.definition.lon for n in network.nodes])
    lats = np.array([n.definition.lat for n in network.nodes])
    return lons, lats


def _get_names(network: 'MetapopulationNetwork') -> List[str]:
    """Extract node name list."""
    return [n.definition.name for n in network.nodes]


# ═══════════════════════════════════════════════════════════════════════
# 1. NETWORK MAP
# ═══════════════════════════════════════════════════════════════════════

def plot_network_map(
    network: 'MetapopulationNetwork',
    metric_by_node: Optional[np.ndarray] = None,
    title: str = 'Metapopulation Network',
    metric_label: str = 'Metric',
    cmap_name: str = 'plasma',
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Geographic scatter of nodes, coloured by a metric, with connectivity lines.

    Nodes are sized by carrying capacity, coloured by the supplied metric
    (e.g. mortality fraction, resistance, virulence).  Lines between nodes
    show larval connectivity strength (C matrix).

    Args:
        network: MetapopulationNetwork.
        metric_by_node: Array (n_nodes,) of values to colour nodes by.
            If None, colours by carrying capacity.
        title: Plot title.
        metric_label: Colourbar label.
        cmap_name: Matplotlib colourmap name.
        save_path: Optional path to save the figure.

    Returns:
        matplotlib Figure.
    """
    lons, lats = _get_coords(network)
    names = _get_names(network)
    N = network.n_nodes
    K = np.array([n.definition.carrying_capacity for n in network.nodes])

    if metric_by_node is None:
        metric_by_node = K.astype(float)
        metric_label = 'Carrying capacity (K)'

    fig, ax = dark_figure(figsize=(10, 10))

    # Draw connectivity lines (C matrix)
    C = network.C
    c_max = C.max()
    if c_max > 0:
        for i in range(N):
            for j in range(i + 1, N):
                strength = C[i, j] + C[j, i]
                if strength > 0:
                    alpha = min(1.0, 0.1 + 0.9 * strength / (2 * c_max))
                    width = 0.5 + 3.0 * strength / (2 * c_max)
                    ax.plot(
                        [lons[i], lons[j]], [lats[i], lats[j]],
                        color=GRID_COLOR, linewidth=width, alpha=alpha,
                        zorder=1,
                    )

    # Node scatter
    sizes = 100 + 400 * (K / K.max())
    cmap = plt.get_cmap(cmap_name)
    vmin, vmax = metric_by_node.min(), metric_by_node.max()
    if vmin == vmax:
        vmax = vmin + 1.0
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    scatter = ax.scatter(
        lons, lats, c=metric_by_node, cmap=cmap, norm=norm,
        s=sizes, edgecolors='white', linewidths=1.5, zorder=3,
    )

    # Node labels
    for i, name in enumerate(names):
        short = name.split(',')[0]  # "Sitka, AK" → "Sitka"
        ax.annotate(
            short, (lons[i], lats[i]),
            xytext=(8, 8), textcoords='offset points',
            color=TEXT_COLOR, fontsize=9, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.2', facecolor=DARK_PANEL,
                      edgecolor=GRID_COLOR, alpha=0.8),
        )

    cbar = fig.colorbar(scatter, ax=ax, pad=0.02, shrink=0.7)
    cbar.set_label(metric_label, color=TEXT_COLOR, fontsize=11)
    cbar.ax.yaxis.set_tick_params(color=TEXT_COLOR)
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color=TEXT_COLOR)

    ax.set_xlabel('Longitude', fontsize=12)
    ax.set_ylabel('Latitude', fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 2. CONNECTIVITY HEATMAP
# ═══════════════════════════════════════════════════════════════════════

def plot_connectivity_heatmap(
    network: 'MetapopulationNetwork',
    matrix: str = 'C',
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Heatmap of larval connectivity (C) or pathogen dispersal (D) matrix.

    Args:
        network: MetapopulationNetwork.
        matrix: 'C' for larval connectivity, 'D' for pathogen dispersal.
        save_path: Optional path to save the figure.

    Returns:
        matplotlib Figure.
    """
    mat = network.C if matrix == 'C' else network.D
    names = _get_names(network)
    short_names = [n.split(',')[0] for n in names]
    N = len(names)

    fig, ax = dark_figure(figsize=(8, 7))

    im = ax.imshow(mat, cmap='magma', interpolation='nearest')

    ax.set_xticks(range(N))
    ax.set_yticks(range(N))
    ax.set_xticklabels(short_names, rotation=45, ha='right', fontsize=10)
    ax.set_yticklabels(short_names, fontsize=10)

    # Annotate cells
    for i in range(N):
        for j in range(N):
            val = mat[i, j]
            text_color = 'white' if val < mat.max() * 0.5 else 'black'
            ax.text(j, i, f'{val:.3f}', ha='center', va='center',
                    fontsize=9, color=text_color)

    kind = 'Larval Connectivity (C)' if matrix == 'C' else 'Pathogen Dispersal (D)'
    ax.set_title(f'{kind} Matrix', fontsize=14, fontweight='bold')
    ax.set_xlabel('Destination node', fontsize=11)
    ax.set_ylabel('Source node', fontsize=11)

    cbar = fig.colorbar(im, ax=ax, pad=0.02, shrink=0.8)
    cbar.set_label('Probability', color=TEXT_COLOR, fontsize=11)
    cbar.ax.yaxis.set_tick_params(color=TEXT_COLOR)
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color=TEXT_COLOR)

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 3. NORTH-SOUTH GRADIENT
# ═══════════════════════════════════════════════════════════════════════

def plot_north_south_gradient(
    spatial_result: 'SpatialSimResult',
    metric: str = 'mortality',
    network: Optional['MetapopulationNetwork'] = None,
    lats: Optional[np.ndarray] = None,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Latitude (y) vs metric value (x) for each node.

    Shows the expected N-S gradient for mortality, resistance, or virulence.

    Args:
        spatial_result: SpatialSimResult.
        metric: 'mortality', 'resistance', or 'virulence'.
        network: MetapopulationNetwork (for lat/lon). If None, uses
            evenly spaced latitudes from node names.
        lats: Pre-computed latitudes (overrides network).
        save_path: Optional path to save the figure.

    Returns:
        matplotlib Figure.
    """
    sr = spatial_result
    N = sr.n_nodes
    names = sr.node_names if sr.node_names else [f'Node {i}' for i in range(N)]
    short_names = [n.split(',')[0] for n in names]

    # Get latitudes
    if lats is not None:
        lat_arr = lats
    elif network is not None:
        lat_arr = np.array([n.definition.lat for n in network.nodes])
    else:
        # Fallback: evenly spaced
        lat_arr = np.linspace(35, 57, N)

    # Compute metric per node (mean over all years post-epidemic)
    disease_yr = sr.disease_year if sr.disease_year is not None else 3
    post_years = slice(disease_yr, sr.n_years)

    if metric == 'mortality':
        total_dd = sr.yearly_disease_deaths[:, post_years].sum(axis=1).astype(float)
        total_pop = sr.yearly_pop[:, post_years].sum(axis=1).astype(float)
        values = np.where(total_pop > 0, total_dd / total_pop, 0)
        xlabel = 'Cumulative disease mortality fraction'
        color = ACCENT_COLORS[6]  # red
    elif metric == 'resistance':
        values = np.nanmean(sr.yearly_mean_resistance[:, post_years], axis=1)
        xlabel = 'Mean resistance (r̄)'
        color = ACCENT_COLORS[7]  # green
    elif metric == 'virulence':
        if sr.yearly_mean_virulence is not None:
            values = np.nanmean(sr.yearly_mean_virulence[:, post_years], axis=1)
        else:
            values = np.zeros(N)
        xlabel = 'Mean virulence (v̄)'
        color = ACCENT_COLORS[6]
    else:
        raise ValueError(f"Unknown metric: {metric}")

    fig, ax = dark_figure()

    ax.scatter(values, lat_arr, c=color, s=150, edgecolors='white',
               linewidths=1.5, zorder=3)

    for i in range(N):
        ax.annotate(
            short_names[i], (values[i], lat_arr[i]),
            xytext=(10, 0), textcoords='offset points',
            color=TEXT_COLOR, fontsize=10,
        )

    # Trend line
    if N >= 3:
        z = np.polyfit(values[values > 0], lat_arr[values > 0], 1) if np.sum(values > 0) >= 2 else None
        if z is not None:
            x_fit = np.linspace(values.min(), values.max(), 50)
            ax.plot(x_fit, np.polyval(z, x_fit), color=GRID_COLOR,
                    linestyle='--', linewidth=1, alpha=0.7)

    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel('Latitude (°N)', fontsize=12)
    ax.set_title(f'North–South Gradient: {metric.capitalize()}',
                 fontsize=14, fontweight='bold')

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 4. FJORD VS OPEN COAST
# ═══════════════════════════════════════════════════════════════════════

def plot_fjord_vs_open(
    spatial_result: 'SpatialSimResult',
    fjord_idx: int = 1,
    open_indices: Optional[List[int]] = None,
    disease_year: int = 3,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Paired comparison: fjord node vs open-coast nodes.

    4-panel: population, disease deaths, resistance, virulence.

    Args:
        spatial_result: SpatialSimResult.
        fjord_idx: Index of the fjord node (default 1 = Howe Sound).
        open_indices: Indices of open-coast nodes to compare.
        disease_year: Year disease was introduced.
        save_path: Optional path to save the figure.

    Returns:
        matplotlib Figure.
    """
    sr = spatial_result
    N = sr.n_nodes
    names = sr.node_names if sr.node_names else [f'Node {i}' for i in range(N)]

    if open_indices is None:
        open_indices = [i for i in range(N) if i != fjord_idx]

    years = np.arange(sr.n_years)

    fig, axes = dark_figure(nrows=2, ncols=2, figsize=(14, 10))

    fjord_name = names[fjord_idx].split(',')[0]
    fjord_color = ACCENT_COLORS[3]  # teal
    open_color = ACCENT_COLORS[0]   # crimson

    panels = [
        ('Population', sr.yearly_pop, False),
        ('Disease deaths/yr', sr.yearly_disease_deaths, False),
        ('Mean resistance', sr.yearly_mean_resistance, False),
        ('Mean virulence', sr.yearly_mean_virulence, True),
    ]

    for ax, (label, data, is_vir) in zip(axes.flat, panels):
        if data is None:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center',
                    color=TEXT_COLOR, fontsize=12, transform=ax.transAxes)
            ax.set_title(label, fontsize=12, fontweight='bold')
            continue

        # Fjord
        ax.plot(years, data[fjord_idx], color=fjord_color, linewidth=2.5,
                label=f'{fjord_name} (fjord)')

        # Open-coast mean ± range
        open_data = data[open_indices]
        mean_open = np.mean(open_data, axis=0)
        min_open = np.min(open_data, axis=0)
        max_open = np.max(open_data, axis=0)
        ax.plot(years, mean_open, color=open_color, linewidth=2.5,
                label=f'Open coast (mean, n={len(open_indices)})')
        ax.fill_between(years, min_open, max_open, color=open_color,
                        alpha=0.15)

        ax.axvline(disease_year, color=ACCENT_COLORS[4], linestyle=':',
                   linewidth=1.2, alpha=0.7)

        ax.set_title(label, fontsize=12, fontweight='bold')
        ax.set_xlabel('Year', fontsize=10)
        ax.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
                  labelcolor=TEXT_COLOR, fontsize=8, loc='best')

    fig.suptitle(f'Fjord vs Open Coast Comparison',
                 fontsize=15, fontweight='bold', color=TEXT_COLOR, y=1.01)

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 5. METAPOPULATION TIMESERIES
# ═══════════════════════════════════════════════════════════════════════

def plot_metapopulation_timeseries(
    spatial_result: 'SpatialSimResult',
    disease_year: int = 3,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """All nodes' population trajectories overlaid.

    Each node a different colour. Shows synchrony or asynchrony of
    population dynamics across the metapopulation.

    Args:
        spatial_result: SpatialSimResult with yearly_pop.
        disease_year: Year disease was introduced.
        save_path: Optional path to save the figure.

    Returns:
        matplotlib Figure.
    """
    sr = spatial_result
    N = sr.n_nodes
    names = sr.node_names if sr.node_names else [f'Node {i}' for i in range(N)]
    years = np.arange(sr.n_years)

    fig, ax = dark_figure(figsize=(12, 7))

    for i in range(N):
        color = NODE_COLORS[i % len(NODE_COLORS)]
        short = names[i].split(',')[0]
        ax.plot(years, sr.yearly_pop[i], color=color, linewidth=2,
                label=short, alpha=0.85)

    ax.axvline(disease_year, color='white', linestyle=':', linewidth=1.5,
               alpha=0.6, label='Disease intro')

    # Also show total (dashed)
    if sr.yearly_total_pop is not None:
        ax.plot(years, sr.yearly_total_pop, color='white', linewidth=2,
                linestyle='--', alpha=0.5, label='Total')

    ax.set_xlabel('Year', fontsize=12)
    ax.set_ylabel('Population', fontsize=12)
    ax.set_title('Metapopulation Trajectories', fontsize=14, fontweight='bold')
    ax.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
              labelcolor=TEXT_COLOR, fontsize=10, ncol=2)
    ax.set_xlim(0, sr.n_years - 1)
    ax.set_ylim(bottom=0)

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 6. LARVAL FLOW DIAGRAM
# ═══════════════════════════════════════════════════════════════════════

def plot_larval_flow_diagram(
    spatial_result: 'SpatialSimResult',
    network: 'MetapopulationNetwork',
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Arrow diagram showing net larval dispersal between nodes.

    Arrow width proportional to connectivity strength (C matrix).
    Geographic layout using actual lat/lon.

    Args:
        spatial_result: SpatialSimResult.
        network: MetapopulationNetwork with C matrix and node definitions.
        save_path: Optional path to save the figure.

    Returns:
        matplotlib Figure.
    """
    lons, lats = _get_coords(network)
    names = _get_names(network)
    short_names = [n.split(',')[0] for n in names]
    N = network.n_nodes
    C = network.C

    fig, ax = dark_figure(figsize=(10, 10))

    # Draw flow arrows (off-diagonal only)
    c_max = 0
    for i in range(N):
        for j in range(N):
            if i != j and C[i, j] > 0:
                c_max = max(c_max, C[i, j])

    if c_max > 0:
        for i in range(N):
            for j in range(N):
                if i != j and C[i, j] > 0:
                    frac = C[i, j] / c_max
                    alpha = 0.2 + 0.8 * frac
                    width = 0.5 + 4.0 * frac
                    # Slight offset to avoid overlapping opposing arrows
                    dx = lons[j] - lons[i]
                    dy = lats[j] - lats[i]
                    perp = np.array([-dy, dx])
                    norm = np.sqrt(perp[0]**2 + perp[1]**2)
                    if norm > 0:
                        perp = perp / norm * 0.15
                    else:
                        perp = np.array([0, 0])
                    ax.annotate(
                        '',
                        xy=(lons[j] + perp[0], lats[j] + perp[1]),
                        xytext=(lons[i] + perp[0], lats[i] + perp[1]),
                        arrowprops=dict(
                            arrowstyle='->', color=ACCENT_COLORS[3],
                            lw=width, alpha=alpha,
                            mutation_scale=12 + 10 * frac,
                            connectionstyle='arc3,rad=0.1',
                        ),
                    )

    # Draw self-recruitment as circles
    for i in range(N):
        if C[i, i] > 0:
            frac = C[i, i] / c_max if c_max > 0 else 0
            circle = plt.Circle(
                (lons[i], lats[i]), 0.3 + 0.5 * frac,
                fill=False, edgecolor=ACCENT_COLORS[4],
                linewidth=1 + 2 * frac, alpha=0.6,
            )
            ax.add_patch(circle)

    # Node scatter (fixed size for flow diagram)
    K = np.array([n.definition.carrying_capacity for n in network.nodes])
    sizes = 100 + 300 * (K / K.max())
    ax.scatter(lons, lats, s=sizes, color=ACCENT_COLORS[5],
               edgecolors='white', linewidths=1.5, zorder=5)

    for i, name in enumerate(short_names):
        ax.annotate(
            name, (lons[i], lats[i]),
            xytext=(10, 10), textcoords='offset points',
            color=TEXT_COLOR, fontsize=10, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.2', facecolor=DARK_PANEL,
                      edgecolor=GRID_COLOR, alpha=0.8),
        )

    ax.set_xlabel('Longitude', fontsize=12)
    ax.set_ylabel('Latitude', fontsize=12)
    ax.set_title('Larval Dispersal Flow Diagram', fontsize=14, fontweight='bold')

    # Legend for arrow meaning
    ax.annotate(
        'Arrow width ∝ C[i,j]\n(larval connectivity)',
        xy=(0.02, 0.02), xycoords='axes fraction',
        fontsize=9, color=GRID_COLOR,
    )

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 7. SPATIAL EPIDEMIC TIMELINE
# ═══════════════════════════════════════════════════════════════════════

def plot_spatial_epidemic_timeline(
    spatial_result: 'SpatialSimResult',
    disease_year: int = 3,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Gantt-chart style: each node as a horizontal bar showing disease state.

    Green = disease-free, Red = active epidemic, Blue = recovery.
    Classification per year: epidemic if disease_deaths > 5% of pop,
    recovery if pop is growing after epidemic.

    Args:
        spatial_result: SpatialSimResult.
        disease_year: Year disease was introduced.
        save_path: Optional path to save the figure.

    Returns:
        matplotlib Figure.
    """
    sr = spatial_result
    N = sr.n_nodes
    names = sr.node_names if sr.node_names else [f'Node {i}' for i in range(N)]
    short_names = [n.split(',')[0] for n in names]
    years = np.arange(sr.n_years)

    # Classify each node-year
    # 0 = disease-free (green), 1 = epidemic (red), 2 = recovery (blue)
    state_map = np.zeros((N, sr.n_years), dtype=int)

    for i in range(N):
        epidemic_started = False
        epidemic_ended = False
        for yr in range(sr.n_years):
            pop = sr.yearly_pop[i, yr]
            dd = sr.yearly_disease_deaths[i, yr] if yr < sr.n_years else 0
            # Epidemic: disease deaths > 2% of population
            if pop > 0 and dd > 0.02 * pop:
                state_map[i, yr] = 1
                epidemic_started = True
                epidemic_ended = False
            elif epidemic_started and not epidemic_ended:
                # Post-epidemic: check if recovering
                if pop > 0 and yr > 0 and sr.yearly_pop[i, yr] > sr.yearly_pop[i, yr - 1]:
                    state_map[i, yr] = 2  # recovery
                elif pop > 0 and dd > 0:
                    state_map[i, yr] = 1  # still epidemic (some deaths)
                else:
                    state_map[i, yr] = 2  # recovery
                    epidemic_ended = True
            elif epidemic_ended:
                state_map[i, yr] = 2  # post-recovery

    phase_colors = {0: '#2ecc71', 1: '#e74c3c', 2: '#3498db'}
    phase_labels = {0: 'Disease-free', 1: 'Epidemic', 2: 'Recovery'}

    fig, ax = dark_figure(figsize=(14, max(4, N * 1.2)))

    bar_height = 0.7
    for i in range(N):
        for yr in range(sr.n_years):
            ax.barh(
                i, 1, left=yr, height=bar_height,
                color=phase_colors[state_map[i, yr]],
                edgecolor='none', alpha=0.85,
            )

    ax.axvline(disease_year, color='white', linestyle=':', linewidth=1.5,
               alpha=0.7)

    ax.set_yticks(range(N))
    ax.set_yticklabels(short_names, fontsize=11)
    ax.set_xlabel('Year', fontsize=12)
    ax.set_title('Spatial Epidemic Timeline', fontsize=14, fontweight='bold')
    ax.set_xlim(0, sr.n_years)
    ax.set_ylim(-0.5, N - 0.5)
    ax.invert_yaxis()

    # Legend
    patches = [mpatches.Patch(color=phase_colors[k], label=phase_labels[k])
               for k in sorted(phase_colors.keys())]
    ax.legend(handles=patches, facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
              labelcolor=TEXT_COLOR, fontsize=10, loc='lower right')

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 8. NODE FATE MATRIX
# ═══════════════════════════════════════════════════════════════════════

def plot_node_fate_matrix(
    spatial_result: 'SpatialSimResult',
    disease_year: int = 3,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Multi-panel grid: one small plot per node showing pop + disease + resistance.

    Quick overview of all nodes at a glance.

    Args:
        spatial_result: SpatialSimResult.
        disease_year: Year disease was introduced.
        save_path: Optional path to save the figure.

    Returns:
        matplotlib Figure.
    """
    sr = spatial_result
    N = sr.n_nodes
    names = sr.node_names if sr.node_names else [f'Node {i}' for i in range(N)]
    years = np.arange(sr.n_years)

    # Layout: ceil(sqrt(N)) x ceil(N/cols)
    ncols = min(N, max(2, int(np.ceil(np.sqrt(N)))))
    nrows = int(np.ceil(N / ncols))

    fig, axes = dark_figure(nrows=nrows, ncols=ncols,
                            figsize=(5 * ncols, 4 * nrows))
    if isinstance(axes, np.ndarray):
        axes_flat = axes.flat
    else:
        axes_flat = [axes]

    for i, ax in enumerate(axes_flat):
        if i >= N:
            ax.set_visible(False)
            continue

        short = names[i].split(',')[0]

        # Population (left y-axis)
        pop = sr.yearly_pop[i]
        ax.plot(years, pop, color=ACCENT_COLORS[0], linewidth=1.8,
                label='Pop')
        ax.set_ylabel('Pop', fontsize=8, color=ACCENT_COLORS[0])
        ax.tick_params(axis='y', labelcolor=ACCENT_COLORS[0], labelsize=7)

        # Disease deaths (shaded area)
        dd = sr.yearly_disease_deaths[i]
        ax.fill_between(years, 0, dd * 5, color=ACCENT_COLORS[6],
                        alpha=0.2, label='Disease×5')

        # Resistance (right y-axis)
        ax2 = ax.twinx()
        apply_dark_theme(ax=ax2)
        res = sr.yearly_mean_resistance[i]
        ax2.plot(years, res, color=ACCENT_COLORS[7], linewidth=1.5,
                 linestyle='--', label='Resist.')
        ax2.set_ylabel('r̄', fontsize=8, color=ACCENT_COLORS[7])
        ax2.tick_params(axis='y', labelcolor=ACCENT_COLORS[7], labelsize=7)
        ax2.set_ylim(0, max(0.5, np.nanmax(res) * 1.3) if np.any(res > 0) else 0.5)

        # Virulence if available
        if sr.yearly_mean_virulence is not None:
            vir = sr.yearly_mean_virulence[i]
            valid_v = vir > 0
            if np.any(valid_v):
                ax2.plot(years[valid_v], vir[valid_v],
                         color=ACCENT_COLORS[6], linewidth=1.2,
                         linestyle=':', alpha=0.8, label='Vir.')

        ax.axvline(disease_year, color='white', linestyle=':', linewidth=0.8,
                   alpha=0.5)
        ax.set_title(short, fontsize=10, fontweight='bold', pad=3)
        ax.set_xlim(0, sr.n_years - 1)
        ax.set_ylim(bottom=0)

        if i >= N - ncols:
            ax.set_xlabel('Year', fontsize=8)

    fig.suptitle('Node Fate Matrix', fontsize=15, fontweight='bold',
                 color=TEXT_COLOR, y=1.01)

    if save_path:
        save_figure(fig, save_path)
    return fig
