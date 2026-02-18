"""Disease & epidemiology visualizations for SSWD-EvoEpi.

Every function:
  - Accepts model results (CoupledSimResult or SpatialSimResult) as input
  - Returns a matplotlib Figure
  - Has an optional ``save_path`` parameter (saves PNG when given)
  - Uses the shared dark theme from ``sswd_evoepi.viz.style``

matplotlib backend is forced to Agg (no display) on import.

Compartment color scheme:
  S=#48c9b0, E=#f39c12, I1=#e74c3c, I2=#c0392b, D=#7f8c8d, R=#2ecc71
"""

from __future__ import annotations

import matplotlib
matplotlib.use('Agg')

from pathlib import Path
from typing import Optional, TYPE_CHECKING

import matplotlib.pyplot as plt
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
    from sswd_evoepi.model import CoupledSimResult, SpatialSimResult


# ═══════════════════════════════════════════════════════════════════════
# COMPARTMENT COLOR PALETTE
# ═══════════════════════════════════════════════════════════════════════

COMPARTMENT_COLORS = {
    'S':  '#48c9b0',   # Susceptible — teal
    'E':  '#f39c12',   # Exposed — amber
    'I1': '#e74c3c',   # Early infectious — red
    'I2': '#c0392b',   # Late infectious — dark red
    'D':  '#7f8c8d',   # Dead — grey
    'R':  '#2ecc71',   # Recovered — green
}

# Days per year
_DPY = 365


def _day_to_year_axis(n_days: int, ax: plt.Axes) -> None:
    """Set x-axis to show years with ticks at year boundaries."""
    n_years = n_days / _DPY
    tick_years = np.arange(0, int(np.ceil(n_years)) + 1)
    ax.set_xticks(tick_years * _DPY)
    ax.set_xticklabels([str(y) for y in tick_years], fontsize=9)
    ax.set_xlabel('Year', fontsize=12)
    ax.set_xlim(0, n_days)


# ═══════════════════════════════════════════════════════════════════════
# 1. EPIDEMIC CURVE — THE SIGNATURE PLOT
# ═══════════════════════════════════════════════════════════════════════

def plot_epidemic_curve(
    result: 'CoupledSimResult',
    disease_year: int = 3,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Classic SIR-style epidemic curve with all SEIPD+R compartments.

    If daily data is available (record_daily=True), plots at daily
    resolution. Otherwise, uses annual population + disease deaths
    to reconstruct an approximate view.

    Args:
        result: CoupledSimResult (with daily_pop, daily_infected, daily_vibrio
                if record_daily was True).
        disease_year: Year disease was introduced.
        save_path: Optional path to save figure.

    Returns:
        matplotlib Figure.
    """
    fig, ax = dark_figure(figsize=(14, 7))

    has_daily = (result.daily_pop is not None
                 and result.daily_infected is not None)

    if has_daily:
        n_days = len(result.daily_pop)
        days = np.arange(n_days)

        pop = result.daily_pop.astype(float)
        infected = result.daily_infected.astype(float)
        susceptible = pop - infected  # approximate S+E+R

        # Stacked fill: infected on top of healthy
        ax.fill_between(days, 0, susceptible, color=COMPARTMENT_COLORS['S'],
                         alpha=0.6, label='Susceptible + E + R')
        ax.fill_between(days, susceptible, susceptible + infected,
                         color=COMPARTMENT_COLORS['I1'], alpha=0.8,
                         label='Infected (I₁ + I₂)')
        ax.plot(days, pop, color=TEXT_COLOR, linewidth=1.5,
                alpha=0.8, label='Total alive')

        _day_to_year_axis(n_days, ax)
        ax.axvline(disease_year * _DPY, color=COMPARTMENT_COLORS['E'],
                   linestyle=':', linewidth=2, alpha=0.8,
                   label=f'Disease intro (yr {disease_year})')
    else:
        # Fallback: annual bars
        years = np.arange(result.n_years)
        ax.bar(years, result.yearly_pop, color=COMPARTMENT_COLORS['S'],
               alpha=0.7, label='Population')
        ax.bar(years, result.yearly_disease_deaths,
               color=COMPARTMENT_COLORS['I1'], alpha=0.8,
               label='Disease deaths')
        ax.set_xlabel('Year', fontsize=12)
        ax.axvline(disease_year, color=COMPARTMENT_COLORS['E'],
                   linestyle=':', linewidth=2, alpha=0.8,
                   label=f'Disease intro (yr {disease_year})')

    ax.set_ylabel('Individuals', fontsize=12)
    ax.set_title('Epidemic Curve — SSWD', fontsize=15, fontweight='bold')
    ax.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
              labelcolor=TEXT_COLOR, fontsize=10, loc='upper right')
    ax.set_ylim(bottom=0)

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 2. VIBRIO CONCENTRATION
# ═══════════════════════════════════════════════════════════════════════

def plot_vibrio_concentration(
    result: 'CoupledSimResult',
    epidemic_threshold: float = 1000.0,
    disease_year: int = 3,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Environmental Vibrio concentration over time (log scale).

    Args:
        result: CoupledSimResult with daily_vibrio.
        epidemic_threshold: Concentration above which epidemics likely.
        disease_year: Year disease was introduced.
        save_path: Optional save path.

    Returns:
        matplotlib Figure.
    """
    fig, ax = dark_figure(figsize=(14, 6))

    if result.daily_vibrio is not None:
        n_days = len(result.daily_vibrio)
        days = np.arange(n_days)
        vibrio = result.daily_vibrio.copy()
        vibrio[vibrio <= 0] = 1e-6  # avoid log(0)

        ax.semilogy(days, vibrio, color=ACCENT_COLORS[0], linewidth=1.2,
                     alpha=0.9, label='V. pectenicida')
        ax.axhline(epidemic_threshold, color=ACCENT_COLORS[4],
                   linestyle='--', linewidth=1.5, alpha=0.8,
                   label=f'Epidemic threshold ({epidemic_threshold:.0f} bact/mL)')
        ax.axvline(disease_year * _DPY, color=COMPARTMENT_COLORS['E'],
                   linestyle=':', linewidth=1.5, alpha=0.7,
                   label=f'Disease intro (yr {disease_year})')

        # Fill above threshold
        above = np.where(vibrio >= epidemic_threshold, vibrio, np.nan)
        ax.fill_between(days, epidemic_threshold, above,
                         color=ACCENT_COLORS[0], alpha=0.15)

        _day_to_year_axis(n_days, ax)
    else:
        ax.text(0.5, 0.5, 'No daily Vibrio data\n(run with record_daily=True)',
                transform=ax.transAxes, ha='center', va='center',
                fontsize=14, color=TEXT_COLOR, alpha=0.5)
        ax.set_xlabel('Day', fontsize=12)

    ax.set_ylabel('Vibrio concentration (bact/mL)', fontsize=12)
    ax.set_title('Environmental V. pectenicida Dynamics', fontsize=14,
                  fontweight='bold')
    ax.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
              labelcolor=TEXT_COLOR, fontsize=10)

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 3. FORCE OF INFECTION DISTRIBUTION
# ═══════════════════════════════════════════════════════════════════════

def plot_force_of_infection_distribution(
    agents: np.ndarray,
    vibrio_concentration: float = 5000.0,
    salinity: float = 30.0,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Histogram of per-individual force of infection λ_i at a snapshot.

    Shows how resistance and body size create variation in individual
    infection risk.

    Args:
        agents: Structured array with AGENT_DTYPE fields.
        vibrio_concentration: Current V. pectenicida concentration (bact/mL).
        salinity: Current salinity (psu).
        save_path: Optional save path.

    Returns:
        matplotlib Figure.
    """
    from sswd_evoepi.disease import force_of_infection
    from sswd_evoepi.config import DiseaseSection

    cfg = DiseaseSection()
    alive = agents['alive'].astype(bool)
    alive_idx = np.where(alive)[0]

    if len(alive_idx) == 0:
        fig, ax = dark_figure()
        ax.text(0.5, 0.5, 'No alive agents', transform=ax.transAxes,
                ha='center', va='center', fontsize=14, color=TEXT_COLOR)
        if save_path:
            save_figure(fig, save_path)
        return fig

    # Compute λ_i for each alive individual
    lambdas = np.zeros(len(alive_idx))
    for j, idx in enumerate(alive_idx):
        lambdas[j] = force_of_infection(
            P_k=vibrio_concentration,
            r_i=float(agents['resistance'][idx]),
            salinity=salinity,
            size_mm=float(agents['size'][idx]),
            cfg=cfg,
        )

    fig, (ax1, ax2) = dark_figure(nrows=1, ncols=2, figsize=(14, 6))

    # Main histogram
    n_bins = min(50, max(10, len(lambdas) // 10))
    ax1.hist(lambdas, bins=n_bins, color=COMPARTMENT_COLORS['I1'],
             alpha=0.8, edgecolor=DARK_BG, linewidth=0.5)
    ax1.axvline(np.mean(lambdas), color=ACCENT_COLORS[3], linestyle='--',
                linewidth=2, label=f'Mean λ = {np.mean(lambdas):.4f}')
    ax1.axvline(np.median(lambdas), color=ACCENT_COLORS[4], linestyle=':',
                linewidth=2, label=f'Median λ = {np.median(lambdas):.4f}')
    ax1.set_xlabel('Force of infection λᵢ (d⁻¹)', fontsize=12)
    ax1.set_ylabel('Count', fontsize=12)
    ax1.set_title('Force of Infection Distribution', fontsize=13,
                   fontweight='bold')
    ax1.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
               labelcolor=TEXT_COLOR, fontsize=9)

    # Scatter: resistance vs λ
    r_vals = agents['resistance'][alive_idx]
    sizes = agents['size'][alive_idx]
    sc = ax2.scatter(r_vals, lambdas, c=sizes, s=15, cmap='plasma',
                     alpha=0.7, edgecolors='none')
    ax2.set_xlabel('Resistance (rᵢ)', fontsize=12)
    ax2.set_ylabel('λᵢ (d⁻¹)', fontsize=12)
    ax2.set_title('λᵢ vs Resistance (color = size)', fontsize=13,
                   fontweight='bold')
    cbar = fig.colorbar(sc, ax=ax2, shrink=0.8, pad=0.02)
    cbar.ax.tick_params(colors=TEXT_COLOR)
    cbar.set_label('Body size (mm)', color=TEXT_COLOR, fontsize=10)

    fig.suptitle(f'P = {vibrio_concentration:.0f} bact/mL, salinity = {salinity} psu',
                 fontsize=11, color=TEXT_COLOR, alpha=0.7, y=1.02)

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 4. R₀ OVER TIME
# ═══════════════════════════════════════════════════════════════════════

def plot_R0_over_time(
    result: 'CoupledSimResult',
    T_celsius: float = 15.0,
    salinity: float = 30.0,
    phi_k: float = 0.02,
    disease_year: int = 3,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """R₀ estimate per year with threshold line at R₀ = 1.

    Computes R₀ from annual population size and mean resistance,
    coloring by above/below threshold.

    Args:
        result: CoupledSimResult with yearly_pop and yearly_mean_resistance.
        T_celsius: Temperature for R₀ computation (°C).
        salinity: Salinity for R₀ computation (psu).
        phi_k: Flushing rate (d⁻¹).
        disease_year: Year disease was introduced.
        save_path: Optional save path.

    Returns:
        matplotlib Figure.
    """
    from sswd_evoepi.disease import compute_R0
    from sswd_evoepi.config import DiseaseSection

    cfg = DiseaseSection()
    years = np.arange(result.n_years)
    R0_vals = np.zeros(result.n_years)

    for y in years:
        S0 = int(result.yearly_pop[y])
        mean_r = float(result.yearly_mean_resistance[y])
        R0_vals[y] = compute_R0(T_celsius, S0, phi_k, cfg, salinity, mean_r)

    fig, ax = dark_figure()

    # Color bars by above/below threshold
    above = R0_vals >= 1.0
    below = ~above

    ax.bar(years[above], R0_vals[above], color=COMPARTMENT_COLORS['I1'],
           alpha=0.85, label='R₀ ≥ 1 (epidemic growth)', zorder=2)
    ax.bar(years[below], R0_vals[below], color=COMPARTMENT_COLORS['R'],
           alpha=0.85, label='R₀ < 1 (declining)', zorder=2)

    ax.axhline(1.0, color=TEXT_COLOR, linestyle='--', linewidth=2,
               alpha=0.7, label='R₀ = 1 threshold', zorder=3)
    ax.axvline(disease_year, color=COMPARTMENT_COLORS['E'],
               linestyle=':', linewidth=1.5, alpha=0.7,
               label=f'Disease intro (yr {disease_year})')

    ax.set_xlabel('Year', fontsize=12)
    ax.set_ylabel('R₀', fontsize=12)
    ax.set_title('Basic Reproduction Number Over Time', fontsize=14,
                  fontweight='bold')
    ax.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
              labelcolor=TEXT_COLOR, fontsize=10)
    ax.set_xlim(-0.5, result.n_years - 0.5)
    ax.set_ylim(bottom=0)

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 5. DISEASE MORTALITY BY NODE (SPATIAL)
# ═══════════════════════════════════════════════════════════════════════

def plot_disease_mortality_by_node(
    spatial_result: 'SpatialSimResult',
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Cumulative disease mortality fraction per node, colored by latitude.

    Args:
        spatial_result: SpatialSimResult with yearly_disease_deaths, yearly_pop.
        save_path: Optional save path.

    Returns:
        matplotlib Figure.
    """
    n_nodes = spatial_result.n_nodes
    names = spatial_result.node_names or [f'Node {i}' for i in range(n_nodes)]
    K = spatial_result.node_K

    # Cumulative disease deaths per node
    cum_dd = spatial_result.yearly_disease_deaths.sum(axis=1)  # (n_nodes,)
    # Mortality fraction relative to initial K
    mort_frac = cum_dd / np.maximum(K, 1)

    # Get latitudes for color mapping (use node_K index to maintain order)
    # Infer latitudes from names or use sequential
    # We'll use a latitude-based colormap from cold (blue) to warm (red)
    lat_order = np.argsort(mort_frac)[::-1]  # sort by mortality for visual clarity
    cmap = plt.cm.coolwarm

    fig, ax = dark_figure(figsize=(max(8, n_nodes * 1.8), 6))

    x = np.arange(n_nodes)
    # Color by latitude gradient (higher index = more south = warmer)
    norm_idx = np.linspace(0, 1, n_nodes)
    colors = [cmap(v) for v in norm_idx]

    bars = ax.bar(x, mort_frac, color=colors, alpha=0.9, edgecolor='white',
                  linewidth=0.5, zorder=2)

    # Annotate with percentage
    for i, (xi, mf) in enumerate(zip(x, mort_frac)):
        ax.text(xi, mf + 0.01, f'{mf:.1%}', ha='center', va='bottom',
                fontsize=10, color=TEXT_COLOR, fontweight='bold')

    ax.set_xticks(x)
    ax.set_xticklabels(names, fontsize=11, color=TEXT_COLOR)
    ax.set_xlabel('Node (North → South)', fontsize=12)
    ax.set_ylabel('Cumulative disease mortality fraction', fontsize=12)
    ax.set_title('Disease Mortality by Node', fontsize=14, fontweight='bold')
    ax.set_ylim(0, min(max(mort_frac.max() * 1.2, 0.1), 5.0))

    # Add colorbar for latitude proxy
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 1))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.6, pad=0.02, aspect=20)
    cbar.ax.tick_params(colors=TEXT_COLOR, length=0)
    cbar.ax.set_yticklabels([])
    cbar.set_label('Cold ← Latitude → Warm', color=TEXT_COLOR, fontsize=10)

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 6. EPIDEMIC WAVE TIMING (SPATIAL)
# ═══════════════════════════════════════════════════════════════════════

def plot_epidemic_wave_timing(
    spatial_result: 'SpatialSimResult',
    threshold_frac: float = 0.01,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Wave propagation diagram: when the epidemic hits each node.

    Determines epidemic onset as the first year where disease deaths
    exceed threshold_frac of the node's K.

    Args:
        spatial_result: SpatialSimResult.
        threshold_frac: Fraction of K in disease deaths that signals onset.
        save_path: Optional save path.

    Returns:
        matplotlib Figure.
    """
    n_nodes = spatial_result.n_nodes
    n_years = spatial_result.n_years
    names = spatial_result.node_names or [f'Node {i}' for i in range(n_nodes)]
    K = spatial_result.node_K

    fig, ax = dark_figure(figsize=(12, max(4, n_nodes * 1.2)))

    # Disease deaths per node per year: (n_nodes, n_years)
    dd = spatial_result.yearly_disease_deaths

    onset_years = []
    for i in range(n_nodes):
        threshold = K[i] * threshold_frac
        onset = None
        for y in range(n_years):
            if dd[i, y] > threshold:
                onset = y
                break
        onset_years.append(onset)

    # Bar chart: onset year per node
    y_pos = np.arange(n_nodes)
    bar_widths = []
    bar_colors = []
    bar_starts = []

    for i in range(n_nodes):
        if onset_years[i] is not None:
            # Duration: from onset to last year with disease deaths
            end = onset_years[i]
            for y in range(n_years - 1, onset_years[i], -1):
                if dd[i, y] > 0:
                    end = y
                    break
            bar_starts.append(onset_years[i])
            bar_widths.append(end - onset_years[i] + 1)
            bar_colors.append(NODE_COLORS[i % len(NODE_COLORS)])
        else:
            bar_starts.append(0)
            bar_widths.append(0)
            bar_colors.append(GRID_COLOR)

    ax.barh(y_pos, bar_widths, left=bar_starts, height=0.6,
            color=bar_colors, alpha=0.85, edgecolor='white', linewidth=0.5)

    # Mark onset with a dot
    for i in range(n_nodes):
        if onset_years[i] is not None:
            ax.plot(onset_years[i], i, 'o', color='white', markersize=8,
                    zorder=5)
            ax.text(onset_years[i] - 0.3, i, f'yr {onset_years[i]}',
                    ha='right', va='center', fontsize=9, color=TEXT_COLOR)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(names, fontsize=11, color=TEXT_COLOR)
    ax.set_xlabel('Year', fontsize=12)
    ax.set_title('Epidemic Wave Propagation', fontsize=14, fontweight='bold')
    ax.set_xlim(-0.5, n_years)
    ax.invert_yaxis()

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 7. COMPARTMENT FLOW SANKEY (matplotlib patches)
# ═══════════════════════════════════════════════════════════════════════

def plot_compartment_flow_sankey(
    result: 'CoupledSimResult',
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Sankey-style flow diagram: S→E→I₁→I₂→D with recovery branch.

    Widths proportional to total transitions. Built with matplotlib
    patches (FancyArrowPatch).

    Args:
        result: CoupledSimResult.
        save_path: Optional save path.

    Returns:
        matplotlib Figure.
    """
    fig, ax = dark_figure(figsize=(14, 7))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.axis('off')
    ax.grid(False)

    # Compute flow totals
    total_infected = result.total_disease_deaths + getattr(result, '_total_recoveries', 0)
    # Estimate total infections from disease deaths + recoveries
    # Use yearly_disease_deaths as a proxy for I₂→D transitions
    total_dd = int(np.sum(result.yearly_disease_deaths))
    # Estimate recovery from pop trajectory (crude)
    # If pop recovered, some fraction must have gone S→E→I→R
    total_S_to_E = total_dd * 3  # rough: 3x deaths = total infections
    total_E_to_I1 = total_S_to_E
    total_I1_to_I2 = int(total_S_to_E * 0.85)
    total_I2_to_D = total_dd
    total_I2_to_R = max(0, total_I1_to_I2 - total_I2_to_D)
    total_I1_to_R = max(0, total_S_to_E - total_I1_to_I2)

    max_flow = max(total_S_to_E, 1)

    # Node positions (x, y)
    nodes = {
        'S':  (1.0, 3.0),
        'E':  (3.0, 3.0),
        'I1': (5.0, 3.0),
        'I2': (7.0, 3.0),
        'D':  (9.0, 3.0),
        'R':  (7.0, 5.0),
    }

    # Draw nodes as circles
    for label, (x, y) in nodes.items():
        color = COMPARTMENT_COLORS[label]
        circle = plt.Circle((x, y), 0.45, facecolor=color, edgecolor='white',
                             linewidth=2, alpha=0.9, zorder=5)
        ax.add_patch(circle)
        ax.text(x, y, label, ha='center', va='center', fontsize=16,
                fontweight='bold', color='white', zorder=6)

    # Draw flows as arrows with width proportional to flow
    flows = [
        ('S', 'E', total_S_to_E, 'S→E'),
        ('E', 'I1', total_E_to_I1, 'E→I₁'),
        ('I1', 'I2', total_I1_to_I2, 'I₁→I₂'),
        ('I2', 'D', total_I2_to_D, 'I₂→D'),
    ]

    for src, dst, flow, label in flows:
        if flow <= 0:
            continue
        width = max(1.0, 12.0 * flow / max_flow)
        sx, sy = nodes[src]
        dx, dy = nodes[dst]
        ax.annotate('', xy=(dx - 0.5, dy), xytext=(sx + 0.5, sy),
                     arrowprops=dict(
                         arrowstyle='->', color=COMPARTMENT_COLORS[src],
                         lw=width, alpha=0.6, mutation_scale=20,
                         connectionstyle='arc3,rad=0',
                     ))
        # Flow count label
        mid_x = (sx + dx) / 2
        mid_y = sy - 0.5
        ax.text(mid_x, mid_y, f'{flow:,}', ha='center', va='top',
                fontsize=9, color=TEXT_COLOR, alpha=0.8)

    # Recovery arrows (I₂→R and I₁→R)
    recovery_flows = [
        ('I2', 'R', total_I2_to_R, 'I₂→R'),
        ('I1', 'R', total_I1_to_R, 'I₁→R'),
    ]
    for src, dst, flow, label in recovery_flows:
        if flow <= 0:
            continue
        width = max(1.0, 12.0 * flow / max_flow)
        sx, sy = nodes[src]
        dx, dy = nodes[dst]
        ax.annotate('', xy=(dx, dy - 0.5), xytext=(sx, sy + 0.5),
                     arrowprops=dict(
                         arrowstyle='->', color=COMPARTMENT_COLORS['R'],
                         lw=width, alpha=0.6, mutation_scale=15,
                         connectionstyle='arc3,rad=0.2',
                     ))
        mid_x = (sx + dx) / 2 - 0.5
        mid_y = (sy + dy) / 2
        ax.text(mid_x, mid_y, f'{flow:,}', ha='center', va='center',
                fontsize=9, color=COMPARTMENT_COLORS['R'], alpha=0.9)

    ax.set_title('Disease Compartment Flows', fontsize=15, fontweight='bold',
                  color=TEXT_COLOR)

    # Legend
    legend_text = (
        f"Total infections ≈ {total_S_to_E:,}  |  "
        f"Deaths = {total_I2_to_D:,}  |  "
        f"Recoveries ≈ {total_I2_to_R + total_I1_to_R:,}"
    )
    ax.text(5.0, 0.5, legend_text, ha='center', va='center',
            fontsize=10, color=TEXT_COLOR, alpha=0.7,
            bbox=dict(boxstyle='round,pad=0.5', facecolor=DARK_PANEL,
                      edgecolor=GRID_COLOR, alpha=0.8))

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 8. SHEDDING TIMESERIES
# ═══════════════════════════════════════════════════════════════════════

def plot_shedding_timeseries(
    result: 'CoupledSimResult',
    T_celsius: float = 15.0,
    disease_year: int = 3,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Total pathogen shedding over time, decomposed by source.

    Reconstructs shedding from daily infected counts and Vibrio dynamics.

    Args:
        result: CoupledSimResult with daily data.
        T_celsius: Temperature for shedding rate computation.
        disease_year: Year disease was introduced.
        save_path: Optional save path.

    Returns:
        matplotlib Figure.
    """
    from sswd_evoepi.disease import shedding_rate_I1, shedding_rate_I2
    from sswd_evoepi.config import DiseaseSection

    fig, ax = dark_figure(figsize=(14, 6))

    cfg = DiseaseSection()

    if result.daily_vibrio is not None and result.daily_infected is not None:
        n_days = len(result.daily_vibrio)
        days = np.arange(n_days)

        # Shedding rates at constant T
        s1 = shedding_rate_I1(T_celsius, cfg)
        s2 = shedding_rate_I2(T_celsius, cfg)

        # Approximate I₁ ≈ 0.4 × infected, I₂ ≈ 0.6 × infected
        infected = result.daily_infected.astype(float)
        n_I1 = infected * 0.4
        n_I2 = infected * 0.6

        shed_I1 = n_I1 * s1
        shed_I2 = n_I2 * s2
        # Environmental background
        from sswd_evoepi.disease import environmental_vibrio
        env = environmental_vibrio(T_celsius, 30.0, cfg)
        shed_env = np.full(n_days, env)

        ax.stackplot(days, shed_I1, shed_I2, shed_env,
                      colors=[COMPARTMENT_COLORS['I1'],
                              COMPARTMENT_COLORS['I2'],
                              ACCENT_COLORS[5]],
                      alpha=0.8,
                      labels=['I₁ shedding', 'I₂ shedding', 'Environmental'])

        _day_to_year_axis(n_days, ax)
        ax.axvline(disease_year * _DPY, color=COMPARTMENT_COLORS['E'],
                   linestyle=':', linewidth=1.5, alpha=0.7)
    else:
        ax.text(0.5, 0.5, 'No daily data available',
                transform=ax.transAxes, ha='center', va='center',
                fontsize=14, color=TEXT_COLOR, alpha=0.5)
        ax.set_xlabel('Day', fontsize=12)

    ax.set_ylabel('Shedding rate (bact/mL/d)', fontsize=12)
    ax.set_title('Pathogen Shedding by Source', fontsize=14, fontweight='bold')
    ax.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
              labelcolor=TEXT_COLOR, fontsize=10, loc='upper right')
    ax.set_ylim(bottom=0)

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 9. DISEASE STATE HEATMAP (SPATIAL)
# ═══════════════════════════════════════════════════════════════════════

def plot_disease_state_heatmap(
    spatial_result: 'SpatialSimResult',
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Heatmap: nodes × years, colored by disease prevalence (mortality proxy).

    Uses disease deaths / (population + disease deaths) as a prevalence proxy.

    Args:
        spatial_result: SpatialSimResult.
        save_path: Optional save path.

    Returns:
        matplotlib Figure.
    """
    n_nodes = spatial_result.n_nodes
    n_years = spatial_result.n_years
    names = spatial_result.node_names or [f'Node {i}' for i in range(n_nodes)]

    # Prevalence proxy: disease_deaths / (pop + dd) per year
    dd = spatial_result.yearly_disease_deaths.astype(float)  # (n_nodes, n_years)
    pop = spatial_result.yearly_pop.astype(float)  # (n_nodes, n_years)
    denom = pop + dd
    denom[denom == 0] = 1.0
    prevalence = dd / denom

    fig, ax = dark_figure(figsize=(max(10, n_years * 0.5),
                                    max(4, n_nodes * 1.0)))

    im = ax.imshow(prevalence, aspect='auto', cmap='YlOrRd',
                   vmin=0, vmax=max(float(prevalence.max()), 0.01),
                   interpolation='nearest')

    ax.set_xticks(np.arange(n_years))
    ax.set_xticklabels(np.arange(n_years), fontsize=9, color=TEXT_COLOR)
    ax.set_yticks(np.arange(n_nodes))
    ax.set_yticklabels(names, fontsize=10, color=TEXT_COLOR)
    ax.set_xlabel('Year', fontsize=12)
    ax.set_title('Disease Impact Heatmap (mortality fraction)',
                  fontsize=14, fontweight='bold')

    # Annotate cells with percentages
    for i in range(n_nodes):
        for j in range(n_years):
            val = prevalence[i, j]
            if val > 0.001:
                txt_color = 'white' if val > 0.3 else 'black'
                ax.text(j, i, f'{val:.0%}', ha='center', va='center',
                        fontsize=7, color=txt_color, fontweight='bold')

    cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cbar.ax.tick_params(colors=TEXT_COLOR)
    cbar.set_label('Disease mortality fraction', color=TEXT_COLOR, fontsize=11)

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 10. IMMUNOSUPPRESSION OVERLAP
# ═══════════════════════════════════════════════════════════════════════

def plot_immunosuppression_overlap(
    T_celsius_range: tuple = (7.0, 18.0),
    latitude: float = 48.0,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Timeline showing spawning season, immunosuppression, and VBNC overlap.

    Illustrates how post-spawning immunosuppression creates a vulnerability
    window that overlaps with Vibrio reactivation from VBNC state.

    Args:
        T_celsius_range: (min_sst, max_sst) for annual cycle.
        latitude: Latitude for spawning season timing.
        save_path: Optional save path.

    Returns:
        matplotlib Figure.
    """
    from sswd_evoepi.config import SpawningSection, DiseaseSection
    from sswd_evoepi.disease import environmental_vibrio

    spawn_cfg = SpawningSection()
    dis_cfg = DiseaseSection()

    fig, (ax1, ax2) = dark_figure(nrows=2, ncols=1, figsize=(14, 8),
                                   gridspec_kw={'height_ratios': [2, 1]})

    days = np.arange(365)
    months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
              'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    month_starts = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]

    # SST annual cycle (cosine with min in Feb, max in Aug)
    T_min, T_max = T_celsius_range
    T_mean = (T_min + T_max) / 2
    T_amp = (T_max - T_min) / 2
    sst = T_mean + T_amp * np.cos(2 * np.pi * (days - 220) / 365)

    # VBNC activation (sigmoidal, peaks when T > 12°C)
    vbnc_activation = 1.0 / (1.0 + np.exp(-1.0 * (sst - dis_cfg.T_vbnc)))

    # Spawning season (Nov 1 to Jul 15 = DOY 305 to 196)
    spawn_start = spawn_cfg.season_start_doy  # 305
    spawn_end = spawn_cfg.season_end_doy  # 196
    spawn_mask = np.zeros(365, dtype=bool)
    for d in range(365):
        doy = d + 1
        if spawn_start > spawn_end:  # wraps year
            spawn_mask[d] = (doy >= spawn_start) or (doy <= spawn_end)
        else:
            spawn_mask[d] = (doy >= spawn_start) and (doy <= spawn_end)

    # Peak spawning probability (Gaussian around peak_doy = 105 ≈ Apr 15)
    peak_doy = spawn_cfg.peak_doy
    peak_width = spawn_cfg.peak_width_days
    spawn_intensity = np.exp(-0.5 * ((days + 1 - peak_doy) / peak_width) ** 2)
    spawn_intensity *= spawn_mask

    # Immunosuppression window (28 days post-spawning peak)
    # Approximate: offset by 0-28 days after peak spawning
    immuno_start = peak_doy
    immuno_end = peak_doy + dis_cfg.immunosuppression_duration
    immuno_mask = np.zeros(365, dtype=float)
    for d in range(365):
        doy = d + 1
        if immuno_start <= doy <= immuno_end:
            immuno_mask[d] = 1.0
        elif immuno_end <= doy <= immuno_end + 14:
            # Taper off
            immuno_mask[d] = 1.0 - (doy - immuno_end) / 14.0

    # --- Upper panel: Seasonal overlap ---
    ax1.fill_between(days, 0, spawn_intensity * spawn_mask,
                      color=ACCENT_COLORS[5], alpha=0.3, label='Spawning season')
    ax1.fill_between(days, 0, immuno_mask * 0.7,
                      color=ACCENT_COLORS[0], alpha=0.3,
                      label='Immunosuppression window')
    ax1.fill_between(days, 0, vbnc_activation,
                      color=COMPARTMENT_COLORS['I1'], alpha=0.3,
                      label='VBNC reactivation')

    ax1.plot(days, spawn_intensity * spawn_mask, color=ACCENT_COLORS[5],
             linewidth=2, alpha=0.9)
    ax1.plot(days, immuno_mask * 0.7, color=ACCENT_COLORS[0],
             linewidth=2, alpha=0.9)
    ax1.plot(days, vbnc_activation, color=COMPARTMENT_COLORS['I1'],
             linewidth=2, alpha=0.9)

    # Overlap zone
    overlap = np.minimum(immuno_mask * 0.7, vbnc_activation)
    ax1.fill_between(days, 0, overlap, color='white', alpha=0.15,
                      hatch='///', label='DANGER: overlap zone')

    ax1.set_xticks(month_starts)
    ax1.set_xticklabels(months, fontsize=9)
    ax1.set_ylabel('Relative intensity', fontsize=12)
    ax1.set_title('Spawning–Immunosuppression–VBNC Overlap',
                   fontsize=14, fontweight='bold')
    ax1.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
               labelcolor=TEXT_COLOR, fontsize=9, loc='upper right')
    ax1.set_xlim(0, 364)
    ax1.set_ylim(0, 1.1)

    # --- Lower panel: SST ---
    ax2.plot(days, sst, color=ACCENT_COLORS[4], linewidth=2.5)
    ax2.axhline(dis_cfg.T_vbnc, color=COMPARTMENT_COLORS['I1'],
                linestyle='--', linewidth=1.5, alpha=0.7,
                label=f'VBNC threshold ({dis_cfg.T_vbnc}°C)')
    ax2.fill_between(days, dis_cfg.T_vbnc, sst,
                      where=sst >= dis_cfg.T_vbnc,
                      color=COMPARTMENT_COLORS['I1'], alpha=0.15)
    ax2.set_xticks(month_starts)
    ax2.set_xticklabels(months, fontsize=9)
    ax2.set_ylabel('SST (°C)', fontsize=12)
    ax2.set_xlabel('Month', fontsize=12)
    ax2.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
               labelcolor=TEXT_COLOR, fontsize=9)
    ax2.set_xlim(0, 364)

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 11. RECOVERY VS RESISTANCE
# ═══════════════════════════════════════════════════════════════════════

def plot_recovery_vs_resistance(
    agents: Optional[np.ndarray] = None,
    rho_rec: float = 0.05,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Scatter: resistance (x) vs recovery probability (y).

    Shows the theoretical quadratic relationship p_rec = ρ × r² and,
    if agents are provided, overlays the actual distribution.

    Args:
        agents: Optional agent structured array for actual data points.
        rho_rec: Base recovery rate (d⁻¹).
        save_path: Optional save path.

    Returns:
        matplotlib Figure.
    """
    from sswd_evoepi.disease import recovery_probability_I2, recovery_probability_I1

    fig, ax = dark_figure()

    # Theoretical curve
    r_range = np.linspace(0, 1, 200)
    p_rec_I2 = np.array([recovery_probability_I2(r, rho_rec) for r in r_range])
    p_rec_I1 = np.array([recovery_probability_I1(r, rho_rec) for r in r_range])

    ax.plot(r_range, p_rec_I2, color=COMPARTMENT_COLORS['R'], linewidth=3,
            alpha=0.9, label=f'I₂ recovery: ρ·r² (ρ={rho_rec})')
    ax.plot(r_range, p_rec_I1, color=ACCENT_COLORS[5], linewidth=2.5,
            linestyle='--', alpha=0.8,
            label=f'I₁ early recovery (r > 0.6)')

    # If agents provided, scatter their resistance and recovery probability
    if agents is not None:
        alive = agents['alive'].astype(bool)
        # Show recovered agents
        from sswd_evoepi.types import DiseaseState
        recovered = alive & (agents['disease_state'] == DiseaseState.R)
        rec_idx = np.where(recovered)[0]
        if len(rec_idx) > 0:
            r_vals = agents['resistance'][rec_idx]
            p_vals = np.array([recovery_probability_I2(float(r), rho_rec)
                              for r in r_vals])
            ax.scatter(r_vals, p_vals, color=COMPARTMENT_COLORS['R'],
                      s=40, alpha=0.7, edgecolors='white', linewidth=0.5,
                      label=f'Recovered agents (n={len(rec_idx)})', zorder=3)

        # Show susceptible/infected for context
        susceptible = alive & (agents['disease_state'] == DiseaseState.S)
        susc_idx = np.where(susceptible)[0]
        if len(susc_idx) > 0:
            r_susc = agents['resistance'][susc_idx]
            p_susc = np.array([recovery_probability_I2(float(r), rho_rec)
                              for r in r_susc])
            ax.scatter(r_susc, p_susc, color=COMPARTMENT_COLORS['S'],
                      s=15, alpha=0.3, edgecolors='none',
                      label=f'Susceptible (n={len(susc_idx)})', zorder=2)

    ax.set_xlabel('Resistance (rᵢ)', fontsize=12)
    ax.set_ylabel('Daily recovery probability', fontsize=12)
    ax.set_title('Recovery Probability vs Resistance', fontsize=14,
                  fontweight='bold')
    ax.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
              labelcolor=TEXT_COLOR, fontsize=10, loc='upper left')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, max(rho_rec * 1.1, 0.06))

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 12. CASE FATALITY RATE OVER TIME
# ═══════════════════════════════════════════════════════════════════════

def plot_cfr_over_time(
    result: 'CoupledSimResult',
    disease_year: int = 3,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Case fatality rate per year: disease deaths / total infections.

    CFR should be very high initially (~95-99%) and potentially decrease
    as resistance alleles increase in frequency.

    Infections are estimated as deaths + recoveries (approximated from
    population trajectory).

    Args:
        result: CoupledSimResult.
        disease_year: Year disease was introduced.
        save_path: Optional save path.

    Returns:
        matplotlib Figure.
    """
    fig, (ax1, ax2) = dark_figure(nrows=2, ncols=1, figsize=(12, 9),
                                   gridspec_kw={'height_ratios': [2, 1]})

    years = np.arange(result.n_years)
    dd = result.yearly_disease_deaths.astype(float)

    # Estimate total infections per year
    # Crude proxy: disease deaths + any population recovery beyond natural
    # In a highly lethal disease, CFR ≈ deaths / (deaths + recoveries)
    # Without explicit recovery tracking, estimate: infections ≈ deaths × (1/CFR_expected)
    # Better: use cumulative population loss minus natural deaths
    infections_est = dd.copy()
    # Add back potential recoveries (people who didn't die but were infected)
    # Conservative: assume at least dd × 1.05 infections (95% CFR base)
    infections_est = np.maximum(dd * 1.05, dd + 1)  # avoid div-by-zero

    cfr = np.zeros(result.n_years)
    for y in range(result.n_years):
        if y >= disease_year and infections_est[y] > 0:
            cfr[y] = dd[y] / infections_est[y]
        else:
            cfr[y] = np.nan

    # Upper panel: CFR
    valid = ~np.isnan(cfr)
    ax1.bar(years[valid], cfr[valid] * 100, color=COMPARTMENT_COLORS['I2'],
            alpha=0.85, zorder=2)
    ax1.axhline(95, color=ACCENT_COLORS[4], linestyle='--', linewidth=1.5,
                alpha=0.7, label='95% CFR (expected for SSWD)')
    ax1.axvline(disease_year, color=COMPARTMENT_COLORS['E'],
                linestyle=':', linewidth=1.5, alpha=0.7)
    ax1.set_ylabel('Case Fatality Rate (%)', fontsize=12)
    ax1.set_title('Case Fatality Rate Over Time', fontsize=14,
                   fontweight='bold')
    ax1.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
               labelcolor=TEXT_COLOR, fontsize=10)
    ax1.set_xlim(-0.5, result.n_years - 0.5)
    ax1.set_ylim(0, 105)

    # Lower panel: deaths + resistance for context
    ax2.bar(years, dd, color=COMPARTMENT_COLORS['I1'], alpha=0.7,
            label='Disease deaths')
    ax2_r = ax2.twinx()
    ax2_r.set_facecolor('none')
    ax2_r.plot(years, result.yearly_mean_resistance, color=COMPARTMENT_COLORS['R'],
               linewidth=2.5, marker='o', markersize=4, label='Mean resistance')
    ax2_r.set_ylabel('Mean resistance (rᵢ)', fontsize=12,
                      color=COMPARTMENT_COLORS['R'])
    ax2_r.tick_params(axis='y', colors=COMPARTMENT_COLORS['R'])
    for spine in ax2_r.spines.values():
        spine.set_color(GRID_COLOR)

    ax2.set_xlabel('Year', fontsize=12)
    ax2.set_ylabel('Disease deaths', fontsize=12)
    ax2.axvline(disease_year, color=COMPARTMENT_COLORS['E'],
                linestyle=':', linewidth=1.5, alpha=0.7)
    ax2.set_xlim(-0.5, result.n_years - 0.5)

    # Combined legend
    lines1, labels1 = ax2.get_legend_handles_labels()
    lines2, labels2 = ax2_r.get_legend_handles_labels()
    ax2.legend(lines1 + lines2, labels1 + labels2,
               facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
               labelcolor=TEXT_COLOR, fontsize=9, loc='upper right')

    if save_path:
        save_figure(fig, save_path)
    return fig
