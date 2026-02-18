"""Dashboard & composite view visualizations for SSWD-EvoEpi.

Multi-panel figures that combine insights from population, disease,
genetics, coevolution, spatial, and sensitivity analysis modules.

Every function:
  - Returns a matplotlib Figure
  - Has an optional ``save_path`` parameter (saves PNG when given)
  - Uses the shared dark theme from ``sswd_evoepi.viz.style``

matplotlib backend is forced to Agg (no display) on import.

Authors: Anton 🔬 & Willem Weertman
"""

from __future__ import annotations

import matplotlib
matplotlib.use('Agg')

from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, TYPE_CHECKING

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import numpy as np

from sswd_evoepi.viz.style import (
    ACCENT_COLORS,
    DARK_BG,
    DARK_PANEL,
    DEATH_COLORS,
    GRID_COLOR,
    NODE_COLORS,
    STAGE_COLORS,
    TEXT_COLOR,
    apply_dark_theme,
    dark_figure,
    save_figure,
)

if TYPE_CHECKING:
    from sswd_evoepi.model import CoupledSimResult, SpatialSimResult

# ═══════════════════════════════════════════════════════════════════════
# CONSTANTS
# ═══════════════════════════════════════════════════════════════════════

COMPARTMENT_COLORS = {
    'S': '#48c9b0', 'E': '#f39c12', 'I1': '#e74c3c',
    'I2': '#c0392b', 'D': '#7f8c8d', 'R': '#2ecc71',
}

VIRULENCE_LOW = '#2ecc71'
VIRULENCE_HIGH = '#e74c3c'

# Category colors for SA parameter grouping
SA_CATEGORY_COLORS = {
    'disease':            '#e74c3c',
    'population':         '#3498db',
    'genetics':           '#2ecc71',
    'spawning':           '#f39c12',
    'spatial':            '#9b59b6',
    'pathogen_evolution':  '#e67e22',
}

_DPY = 365


def _param_category(name: str) -> str:
    """Extract category from dotted param name (e.g. 'disease.a_exposure' → 'disease')."""
    return name.split('.')[0] if '.' in name else 'other'


def _legend(ax, **kwargs):
    """Standard dark legend. Skips if no labeled artists."""
    handles, labels = ax.get_legend_handles_labels()
    if not labels:
        return
    defaults = dict(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
                    labelcolor=TEXT_COLOR, fontsize=9)
    defaults.update(kwargs)
    ax.legend(**defaults)


# ═══════════════════════════════════════════════════════════════════════
# 1. SIMULATION DASHBOARD — THE MASTER OVERVIEW
# ═══════════════════════════════════════════════════════════════════════

def plot_simulation_dashboard(
    result: 'CoupledSimResult',
    title: str = 'SSWD-EvoEpi Simulation Dashboard',
    disease_year: int = 3,
    carrying_capacity: int = 500,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Master overview dashboard: 2×3 grid covering all model dimensions.

    Panels:
      [0,0] Population trajectory with K and disease-year markers
      [0,1] Epidemic curve (susceptible / infected stacked)
      [0,2] Mean resistance over time
      [1,0] Disease deaths per year (bar)
      [1,1] Virulence trajectory (PE) or Vibrio concentration
      [1,2] Co-evolution phase portrait (resistance vs virulence)

    Args:
        result: CoupledSimResult from a single-node run.
        title: Super-title for the figure.
        disease_year: Year disease was introduced.
        carrying_capacity: K for the node.
        save_path: Optional save path.

    Returns:
        matplotlib Figure (20×12 inches).
    """
    fig, axes = plt.subplots(2, 3, figsize=(20, 12))
    fig.patch.set_facecolor(DARK_BG)
    for ax in axes.flat:
        apply_dark_theme(ax=ax)

    years = np.arange(result.n_years)

    # ── [0,0] Population trajectory ──────────────────────────────────
    ax = axes[0, 0]
    ax.plot(years, result.yearly_pop, color=ACCENT_COLORS[3], linewidth=2.5,
            label='Population')
    ax.axhline(carrying_capacity, color=ACCENT_COLORS[4], linestyle='--',
               alpha=0.6, label=f'K={carrying_capacity}')
    ax.axvline(disease_year, color=COMPARTMENT_COLORS['E'], linestyle=':',
               alpha=0.7, label='Disease intro')
    ax.fill_between(years, result.yearly_pop, alpha=0.15, color=ACCENT_COLORS[3])
    ax.set_xlabel('Year')
    ax.set_ylabel('Population')
    ax.set_title('Population Trajectory', fontsize=12, fontweight='bold')
    _legend(ax)

    # ── [0,1] Epidemic curve ─────────────────────────────────────────
    ax = axes[0, 1]
    pop = result.yearly_pop
    infected_frac = np.zeros_like(years, dtype=float)
    if result.daily_infected is not None and len(result.daily_infected) > 0:
        # Compute yearly mean infected fraction
        daily_pop = result.daily_pop if result.daily_pop is not None else None
        for yr in years:
            start, end = yr * _DPY, (yr + 1) * _DPY
            if result.daily_infected is not None and end <= len(result.daily_infected):
                inf_slice = result.daily_infected[start:end]
                pop_slice = daily_pop[start:end] if daily_pop is not None else np.full(365, pop[yr])
                valid = pop_slice > 0
                if valid.any():
                    infected_frac[yr] = np.mean(inf_slice[valid] / pop_slice[valid])
    else:
        # Estimate from disease deaths
        if result.yearly_disease_deaths is not None:
            infected_frac = np.where(pop > 0,
                                     result.yearly_disease_deaths / np.maximum(pop, 1), 0)

    susceptible_frac = 1.0 - infected_frac
    ax.fill_between(years, 0, susceptible_frac, color=COMPARTMENT_COLORS['S'],
                    alpha=0.7, label='Susceptible')
    ax.fill_between(years, susceptible_frac, 1.0, color=COMPARTMENT_COLORS['I1'],
                    alpha=0.8, label='Infected')
    ax.axvline(disease_year, color=COMPARTMENT_COLORS['E'], linestyle=':',
               alpha=0.7)
    ax.set_xlabel('Year')
    ax.set_ylabel('Fraction')
    ax.set_ylim(0, 1)
    ax.set_title('Disease Prevalence', fontsize=12, fontweight='bold')
    _legend(ax)

    # ── [0,2] Mean resistance ────────────────────────────────────────
    ax = axes[0, 2]
    if result.yearly_mean_resistance is not None:
        ax.plot(years, result.yearly_mean_resistance, color=ACCENT_COLORS[7],
                linewidth=2.5, label='Mean r̄')
        ax.fill_between(years, result.yearly_mean_resistance, alpha=0.15,
                        color=ACCENT_COLORS[7])
        ax.axvline(disease_year, color=COMPARTMENT_COLORS['E'], linestyle=':',
                   alpha=0.7)
    ax.set_xlabel('Year')
    ax.set_ylabel('Mean Resistance (r̄)')
    ax.set_title('Host Resistance Evolution', fontsize=12, fontweight='bold')
    _legend(ax)

    # ── [1,0] Disease deaths per year ────────────────────────────────
    ax = axes[1, 0]
    if result.yearly_disease_deaths is not None:
        colors = [DEATH_COLORS['disease'] if yr >= disease_year else '#333'
                  for yr in years]
        ax.bar(years, result.yearly_disease_deaths, color=colors, alpha=0.85,
               edgecolor='white', linewidth=0.3)
    ax.set_xlabel('Year')
    ax.set_ylabel('Disease Deaths')
    ax.set_title('Annual Disease Mortality', fontsize=12, fontweight='bold')

    # ── [1,1] Virulence trajectory OR Vibrio ─────────────────────────
    ax = axes[1, 1]
    has_pe = (result.yearly_mean_virulence is not None and
              np.any(result.yearly_mean_virulence != 0))
    if has_pe:
        ax.plot(years, result.yearly_mean_virulence, color='#e74c3c',
                linewidth=2.5, label='Mean v̄')
        if result.yearly_virulence_new_infections is not None:
            ax.plot(years, result.yearly_virulence_new_infections,
                    color='#f39c12', linewidth=1.5, linestyle='--',
                    label='v̄ new infections')
        ax.set_ylabel('Virulence')
        ax.set_title('Pathogen Virulence Trajectory', fontsize=12,
                      fontweight='bold')
    else:
        # Fallback: Vibrio concentration from daily data
        if result.daily_vibrio is not None and len(result.daily_vibrio) > 0:
            days = np.arange(len(result.daily_vibrio))
            ax.plot(days / _DPY, result.daily_vibrio, color=ACCENT_COLORS[6],
                    linewidth=0.8, alpha=0.7)
            ax.axvline(disease_year, color=COMPARTMENT_COLORS['E'],
                       linestyle=':', alpha=0.7)
            ax.set_ylabel('Vibrio (cells/L)')
            ax.set_title('Environmental Vibrio', fontsize=12, fontweight='bold')
        else:
            ax.text(0.5, 0.5, 'No PE / Vibrio data', transform=ax.transAxes,
                    ha='center', va='center', color=TEXT_COLOR, fontsize=12)
            ax.set_title('Virulence / Vibrio', fontsize=12, fontweight='bold')
    ax.set_xlabel('Year')
    _legend(ax)

    # ── [1,2] Co-evolution phase portrait ────────────────────────────
    ax = axes[1, 2]
    if has_pe and result.yearly_mean_resistance is not None:
        r = result.yearly_mean_resistance
        v = result.yearly_mean_virulence
        # Colored by time
        cmap = plt.cm.viridis
        for i in range(len(years) - 1):
            ax.plot(r[i:i+2], v[i:i+2], color=cmap(i / max(len(years) - 1, 1)),
                    linewidth=2)
        ax.scatter(r[0], v[0], color='white', s=80, zorder=5, marker='o',
                   edgecolors='black', label='Start')
        ax.scatter(r[-1], v[-1], color='#e94560', s=80, zorder=5, marker='*',
                   edgecolors='black', label='End')
        ax.set_xlabel('Mean Host Resistance (r̄)')
        ax.set_ylabel('Mean Pathogen Virulence (v̄)')
        ax.set_title('Co-evolution Phase Portrait', fontsize=12,
                      fontweight='bold')
        # Add colorbar for time
        sm = plt.cm.ScalarMappable(cmap=cmap,
                                    norm=plt.Normalize(0, result.n_years))
        sm.set_array([])
        cb = fig.colorbar(sm, ax=ax, label='Year', shrink=0.8)
        cb.ax.yaxis.set_tick_params(color=TEXT_COLOR)
        cb.ax.yaxis.label.set_color(TEXT_COLOR)
        plt.setp(cb.ax.get_yticklabels(), color=TEXT_COLOR)
    else:
        ax.text(0.5, 0.5, 'Co-evolution requires PE', transform=ax.transAxes,
                ha='center', va='center', color=TEXT_COLOR, fontsize=12)
        ax.set_title('Co-evolution Phase Portrait', fontsize=12,
                      fontweight='bold')
    _legend(ax)

    fig.suptitle(title, fontsize=16, fontweight='bold', color=TEXT_COLOR, y=0.98)
    fig.subplots_adjust(hspace=0.35, wspace=0.3, top=0.92)
    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 2. SPATIAL DASHBOARD — MULTI-NODE OVERVIEW
# ═══════════════════════════════════════════════════════════════════════

def plot_spatial_dashboard(
    spatial_result: 'SpatialSimResult',
    title: str = 'Spatial Metapopulation Dashboard',
    network: Any = None,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Multi-node overview dashboard: 2×3 grid.

    Panels:
      [0,0] All nodes population trajectories
      [0,1] Disease mortality by node (bar)
      [0,2] Network map with population sizes
      [1,0] Mean resistance per node over time
      [1,1] N-S mortality gradient
      [1,2] Fjord vs open coast comparison

    Args:
        spatial_result: SpatialSimResult from run_spatial_simulation.
        title: Figure super-title.
        network: Optional MetapopulationNetwork (for geo coords).
        save_path: Optional save path.

    Returns:
        matplotlib Figure (20×12 inches).
    """
    fig, axes = plt.subplots(2, 3, figsize=(20, 12))
    fig.patch.set_facecolor(DARK_BG)
    for ax in axes.flat:
        apply_dark_theme(ax=ax)

    sr = spatial_result
    n_nodes = sr.n_nodes
    n_years = sr.n_years
    years = np.arange(n_years)
    names = sr.node_names or [f'Node {i}' for i in range(n_nodes)]

    # ── [0,0] All nodes population trajectories ─────────────────────
    ax = axes[0, 0]
    for i in range(n_nodes):
        ax.plot(years, sr.yearly_pop[i], color=NODE_COLORS[i % len(NODE_COLORS)],
                linewidth=2, label=names[i])
    ax.set_xlabel('Year')
    ax.set_ylabel('Population')
    ax.set_title('Population by Node', fontsize=12, fontweight='bold')
    _legend(ax, fontsize=8)

    # ── [0,1] Disease mortality by node (stacked bar) ────────────────
    ax = axes[0, 1]
    if sr.yearly_disease_deaths is not None:
        total_deaths = sr.yearly_disease_deaths.sum(axis=1)  # per node
        x = np.arange(n_nodes)
        ax.bar(x, total_deaths, color=[NODE_COLORS[i % len(NODE_COLORS)]
               for i in range(n_nodes)], alpha=0.85, edgecolor='white',
               linewidth=0.5)
        ax.set_xticks(x)
        ax.set_xticklabels(names, fontsize=9, rotation=30, ha='right',
                           color=TEXT_COLOR)
    ax.set_ylabel('Total Disease Deaths')
    ax.set_title('Cumulative Disease Mortality', fontsize=12, fontweight='bold')

    # ── [0,2] Network map ────────────────────────────────────────────
    ax = axes[0, 2]
    if network is not None:
        lons = np.array([n.definition.lon for n in network.nodes])
        lats = np.array([n.definition.lat for n in network.nodes])
    else:
        # Attempt evenly-spaced layout
        lons = np.linspace(-135, -122, n_nodes)
        lats = np.linspace(57, 37, n_nodes)

    final_pop = sr.yearly_pop[:, -1] if sr.yearly_pop is not None else np.ones(n_nodes) * 100
    sizes = np.clip(final_pop / np.maximum(final_pop.max(), 1) * 400, 50, 600)
    ax.scatter(lons, lats, s=sizes, c=[NODE_COLORS[i % len(NODE_COLORS)]
               for i in range(n_nodes)], alpha=0.85, edgecolors='white',
               linewidth=0.8, zorder=5)
    for i, name in enumerate(names):
        ax.annotate(name, (lons[i], lats[i]), textcoords='offset points',
                    xytext=(8, 4), fontsize=8, color=TEXT_COLOR)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title('Network Map (bubble = final pop)', fontsize=12,
                  fontweight='bold')

    # ── [1,0] Mean resistance per node over time ─────────────────────
    ax = axes[1, 0]
    if sr.yearly_mean_resistance is not None:
        for i in range(n_nodes):
            ax.plot(years, sr.yearly_mean_resistance[i],
                    color=NODE_COLORS[i % len(NODE_COLORS)],
                    linewidth=2, label=names[i])
    ax.set_xlabel('Year')
    ax.set_ylabel('Mean Resistance (r̄)')
    ax.set_title('Resistance Evolution by Node', fontsize=12, fontweight='bold')
    _legend(ax, fontsize=8)

    # ── [1,1] N-S mortality gradient ─────────────────────────────────
    ax = axes[1, 1]
    if sr.yearly_disease_deaths is not None and network is not None:
        lats_arr = np.array([n.definition.lat for n in network.nodes])
        total_dd = sr.yearly_disease_deaths.sum(axis=1)
        K = sr.node_K if sr.node_K is not None else sr.yearly_pop[:, 0]
        mort_frac = total_dd / np.maximum(K, 1)
        sort_idx = np.argsort(-lats_arr)  # north → south
        ax.barh(np.arange(n_nodes), mort_frac[sort_idx],
                color=[NODE_COLORS[sort_idx[i] % len(NODE_COLORS)]
                       for i in range(n_nodes)],
                alpha=0.85, edgecolor='white', linewidth=0.5)
        ax.set_yticks(np.arange(n_nodes))
        sorted_names = [names[i] for i in sort_idx]
        ax.set_yticklabels(sorted_names, fontsize=9, color=TEXT_COLOR)
        ax.set_xlabel('Mortality Fraction (deaths/K)')
    elif sr.yearly_disease_deaths is not None:
        total_dd = sr.yearly_disease_deaths.sum(axis=1)
        K = sr.node_K if sr.node_K is not None else sr.yearly_pop[:, 0]
        mort_frac = total_dd / np.maximum(K, 1)
        ax.barh(np.arange(n_nodes), mort_frac,
                color=[NODE_COLORS[i % len(NODE_COLORS)] for i in range(n_nodes)],
                alpha=0.85)
        ax.set_yticks(np.arange(n_nodes))
        ax.set_yticklabels(names, fontsize=9, color=TEXT_COLOR)
        ax.set_xlabel('Mortality Fraction')
    ax.set_title('North → South Mortality Gradient', fontsize=12,
                  fontweight='bold')

    # ── [1,2] Fjord vs open coast comparison ─────────────────────────
    ax = axes[1, 2]
    if network is not None:
        fjord_idx = [i for i, n in enumerate(network.nodes) if n.definition.is_fjord]
        open_idx = [i for i, n in enumerate(network.nodes) if not n.definition.is_fjord]
    else:
        fjord_idx, open_idx = [], list(range(n_nodes))

    if fjord_idx and open_idx and sr.yearly_pop is not None:
        fjord_mean = sr.yearly_pop[fjord_idx].mean(axis=0)
        open_mean = sr.yearly_pop[open_idx].mean(axis=0)
        # Normalize to initial
        fjord_norm = fjord_mean / np.maximum(fjord_mean[0], 1)
        open_norm = open_mean / np.maximum(open_mean[0], 1)
        ax.plot(years, fjord_norm, color='#48c9b0', linewidth=2.5,
                label=f'Fjord (n={len(fjord_idx)})')
        ax.plot(years, open_norm, color='#e94560', linewidth=2.5,
                label=f'Open coast (n={len(open_idx)})')
        ax.axhline(1.0, color=GRID_COLOR, linestyle='--', alpha=0.5)
        ax.set_ylabel('Population (fraction of initial)')
    else:
        ax.text(0.5, 0.5, 'No fjord/open distinction', transform=ax.transAxes,
                ha='center', va='center', color=TEXT_COLOR, fontsize=11)
    ax.set_xlabel('Year')
    ax.set_title('Fjord vs Open Coast', fontsize=12, fontweight='bold')
    _legend(ax)

    fig.suptitle(title, fontsize=16, fontweight='bold', color=TEXT_COLOR, y=0.98)
    fig.subplots_adjust(hspace=0.35, wspace=0.3, top=0.92)
    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 3. SCENARIO COMPARISON — SIDE-BY-SIDE N CONFIGS
# ═══════════════════════════════════════════════════════════════════════

def plot_scenario_comparison(
    results_dict: Dict[str, 'CoupledSimResult'],
    labels: Optional[List[str]] = None,
    title: str = 'Scenario Comparison',
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Side-by-side comparison of N scenarios, 4 rows × N columns.

    Rows: population, resistance, virulence, disease deaths.
    Shared y-axes within each row.

    Args:
        results_dict: Mapping label → CoupledSimResult.
        labels: Optional explicit labels (defaults to dict keys).
        title: Figure super-title.
        save_path: Optional save path.

    Returns:
        matplotlib Figure.
    """
    if labels is None:
        labels = list(results_dict.keys())
    results = [results_dict[k] for k in (labels if labels == list(results_dict.keys())
               else results_dict.keys())]
    # Ensure labels match results
    if labels != list(results_dict.keys()):
        results = list(results_dict.values())

    n_scenarios = len(results)
    row_names = ['Population', 'Mean Resistance', 'Virulence', 'Disease Deaths/yr']
    n_rows = len(row_names)

    fig, axes = plt.subplots(n_rows, n_scenarios,
                             figsize=(5 * n_scenarios, 3.5 * n_rows),
                             sharey='row')
    fig.patch.set_facecolor(DARK_BG)
    if n_scenarios == 1:
        axes = axes.reshape(-1, 1)
    for ax in axes.flat:
        apply_dark_theme(ax=ax)

    scenario_colors = [ACCENT_COLORS[i % len(ACCENT_COLORS)] for i in range(n_scenarios)]

    for col, (label, res) in enumerate(zip(labels, results)):
        years = np.arange(res.n_years)

        # Row 0: Population
        ax = axes[0, col]
        ax.plot(years, res.yearly_pop, color=scenario_colors[col], linewidth=2)
        ax.fill_between(years, res.yearly_pop, alpha=0.15, color=scenario_colors[col])
        if col == 0:
            ax.set_ylabel('Population')
        ax.set_title(label, fontsize=11, fontweight='bold')

        # Row 1: Resistance
        ax = axes[1, col]
        if res.yearly_mean_resistance is not None:
            ax.plot(years, res.yearly_mean_resistance,
                    color=scenario_colors[col], linewidth=2)
        if col == 0:
            ax.set_ylabel('Mean Resistance')

        # Row 2: Virulence
        ax = axes[2, col]
        has_pe = (res.yearly_mean_virulence is not None and
                  np.any(res.yearly_mean_virulence != 0))
        if has_pe:
            ax.plot(years, res.yearly_mean_virulence,
                    color=scenario_colors[col], linewidth=2)
        else:
            ax.text(0.5, 0.5, 'No PE', transform=ax.transAxes,
                    ha='center', va='center', color=TEXT_COLOR, fontsize=10)
        if col == 0:
            ax.set_ylabel('Virulence')

        # Row 3: Disease deaths
        ax = axes[3, col]
        if res.yearly_disease_deaths is not None:
            ax.bar(years, res.yearly_disease_deaths,
                   color=scenario_colors[col], alpha=0.8, edgecolor='white',
                   linewidth=0.3)
        if col == 0:
            ax.set_ylabel('Deaths/yr')
        ax.set_xlabel('Year')

    fig.suptitle(title, fontsize=16, fontweight='bold', color=TEXT_COLOR, y=0.98)
    fig.subplots_adjust(hspace=0.35, wspace=0.3, top=0.92)
    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 4. SENSITIVITY TORNADO — SOBOL TOTAL-ORDER
# ═══════════════════════════════════════════════════════════════════════

def plot_sensitivity_tornado(
    sobol_indices: Dict,
    metric: str,
    top_n: int = 15,
    title: Optional[str] = None,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Horizontal bar chart of Sobol S1 and ST indices for one metric.

    Bars color-coded by parameter category.

    Args:
        sobol_indices: Dict with 'param_names', 'indices' keys.
        metric: Metric name (e.g. 'pop_crash_pct').
        top_n: Show top N parameters by ST.
        title: Optional override title.
        save_path: Optional save path.

    Returns:
        matplotlib Figure.
    """
    params = sobol_indices['param_names']
    idx = sobol_indices['indices'][metric]
    ST = np.array(idx['ST'])
    S1 = np.array(idx['S1'])

    # Sort by ST descending
    order = np.argsort(-ST)[:top_n]

    fig, ax = dark_figure(figsize=(12, max(6, top_n * 0.45)))

    y = np.arange(len(order))
    bar_h = 0.35

    for i, pidx in enumerate(order):
        cat = _param_category(params[pidx])
        color = SA_CATEGORY_COLORS.get(cat, '#95a5a6')
        ax.barh(i + bar_h / 2, ST[pidx], bar_h, color=color, alpha=0.85,
                edgecolor='white', linewidth=0.3)
        ax.barh(i - bar_h / 2, S1[pidx], bar_h, color=color, alpha=0.5,
                edgecolor='white', linewidth=0.3)

    ax.set_yticks(y)
    ax.set_yticklabels([params[i].split('.')[-1] for i in order],
                       fontsize=10, color=TEXT_COLOR)
    ax.invert_yaxis()
    ax.set_xlabel('Sobol Index', fontsize=12)

    # Legend for S1/ST
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#888', alpha=0.85, label='ST (total)'),
        Patch(facecolor='#888', alpha=0.45, label='S1 (first-order)'),
    ]
    # Add category legend
    cats_present = set(_param_category(params[i]) for i in order)
    for cat in sorted(cats_present):
        legend_elements.append(
            Patch(facecolor=SA_CATEGORY_COLORS.get(cat, '#95a5a6'),
                  label=cat.replace('_', ' ').title())
        )
    ax.legend(handles=legend_elements, loc='lower right',
              facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
              labelcolor=TEXT_COLOR, fontsize=9)

    if title is None:
        title = f'Sobol Sensitivity — {metric}'
    ax.set_title(title, fontsize=14, fontweight='bold')

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 5. SENSITIVITY HEATMAP — PARAMS × METRICS
# ═══════════════════════════════════════════════════════════════════════

def plot_sensitivity_heatmap(
    sobol_indices: Dict,
    title: str = 'Sobol ST Heatmap: Parameters × Metrics',
    top_n: Optional[int] = None,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Heatmap of Sobol total-order indices: parameters (y) × metrics (x).

    Parameters sorted by mean ST. Colored by parameter category on y-axis.

    Args:
        sobol_indices: Dict with 'param_names', 'metric_names', 'indices'.
        title: Figure title.
        top_n: Optional limit to top N parameters.
        save_path: Optional save path.

    Returns:
        matplotlib Figure.
    """
    params = sobol_indices['param_names']
    metrics = sobol_indices['metric_names']
    n_params = len(params)
    n_metrics = len(metrics)

    # Build ST matrix
    ST_matrix = np.zeros((n_params, n_metrics))
    for j, m in enumerate(metrics):
        ST_matrix[:, j] = np.array(sobol_indices['indices'][m]['ST'])

    # Sort by mean ST
    mean_ST = ST_matrix.mean(axis=1)
    order = np.argsort(-mean_ST)
    if top_n is not None:
        order = order[:top_n]

    ST_sorted = ST_matrix[order]
    param_labels = [params[i].split('.')[-1] for i in order]
    metric_labels = [m.replace('_', '\n') for m in metrics]

    fig_h = max(6, len(order) * 0.5)
    fig, ax = dark_figure(figsize=(max(10, n_metrics * 0.9), fig_h))

    # Custom colormap: dark → bright
    cmap = mcolors.LinearSegmentedColormap.from_list(
        'sensitivity', ['#1a1a2e', '#2ecc71', '#f39c12', '#e74c3c'], N=256)
    im = ax.imshow(np.clip(ST_sorted, 0, 1), cmap=cmap, aspect='auto',
                   vmin=0, vmax=0.8)

    ax.set_xticks(np.arange(n_metrics))
    ax.set_xticklabels(metric_labels, fontsize=8, color=TEXT_COLOR,
                       rotation=45, ha='right')
    ax.set_yticks(np.arange(len(order)))
    ax.set_yticklabels(param_labels, fontsize=9, color=TEXT_COLOR)

    # Highlight top 5 cells
    top5_flat = np.argsort(-ST_sorted.ravel())[:5]
    for flat_idx in top5_flat:
        row, col = divmod(flat_idx, n_metrics)
        ax.add_patch(plt.Rectangle((col - 0.5, row - 0.5), 1, 1,
                                    fill=False, edgecolor='white',
                                    linewidth=2))

    # Category color bar on left
    for i, pidx in enumerate(order):
        cat = _param_category(params[pidx])
        color = SA_CATEGORY_COLORS.get(cat, '#95a5a6')
        ax.add_patch(plt.Rectangle((-1.5, i - 0.5), 0.8, 1,
                                    fill=True, color=color, alpha=0.8,
                                    clip_on=False))

    cb = fig.colorbar(im, ax=ax, label='Sobol ST', shrink=0.8)
    cb.ax.yaxis.set_tick_params(color=TEXT_COLOR)
    cb.ax.yaxis.label.set_color(TEXT_COLOR)
    plt.setp(cb.ax.get_yticklabels(), color=TEXT_COLOR)

    ax.set_title(title, fontsize=14, fontweight='bold')
    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 6. PARAMETER INTERACTION WEB
# ═══════════════════════════════════════════════════════════════════════

def plot_parameter_interaction_web(
    sobol_indices: Dict,
    metric: str,
    top_n: int = 12,
    title: Optional[str] = None,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Matrix plot showing interaction strength (ST − S1) between parameters.

    Diagonal shows ST. Off-diagonal shows estimated pairwise interaction
    contribution (proportional allocation of interaction budget).

    Args:
        sobol_indices: Dict with 'param_names', 'indices'.
        metric: Which metric to show.
        top_n: Limit to top N parameters by ST.
        title: Optional override title.
        save_path: Optional save path.

    Returns:
        matplotlib Figure.
    """
    params = sobol_indices['param_names']
    idx = sobol_indices['indices'][metric]
    ST = np.array(idx['ST'])
    S1 = np.array(idx['S1'])
    interaction = np.clip(ST - S1, 0, None)  # interaction budget per param

    # Top N by ST
    order = np.argsort(-ST)[:top_n]
    n = len(order)

    # Build approximate interaction matrix
    # Allocate interaction proportionally: I_ij ∝ interaction_i × interaction_j
    inter_vals = interaction[order]
    total_inter = inter_vals.sum()
    if total_inter > 0:
        inter_matrix = np.outer(inter_vals, inter_vals) / total_inter
    else:
        inter_matrix = np.zeros((n, n))

    # Set diagonal to ST
    np.fill_diagonal(inter_matrix, ST[order])

    param_labels = [params[i].split('.')[-1] for i in order]

    fig_size = max(8, n * 0.7)
    fig, ax = dark_figure(figsize=(fig_size, fig_size))

    cmap = mcolors.LinearSegmentedColormap.from_list(
        'interact', ['#1a1a2e', '#3498db', '#f39c12', '#e74c3c'], N=256)
    im = ax.imshow(inter_matrix, cmap=cmap, aspect='equal',
                   vmin=0, vmax=max(ST[order].max(), 0.5))

    ax.set_xticks(np.arange(n))
    ax.set_xticklabels(param_labels, fontsize=8, color=TEXT_COLOR,
                       rotation=45, ha='right')
    ax.set_yticks(np.arange(n))
    ax.set_yticklabels(param_labels, fontsize=8, color=TEXT_COLOR)

    # Annotate values
    for i in range(n):
        for j in range(n):
            val = inter_matrix[i, j]
            if val > 0.01:
                ax.text(j, i, f'{val:.2f}', ha='center', va='center',
                        fontsize=7, color='white' if val > 0.15 else TEXT_COLOR)

    cb = fig.colorbar(im, ax=ax, label='Index', shrink=0.8)
    cb.ax.yaxis.set_tick_params(color=TEXT_COLOR)
    cb.ax.yaxis.label.set_color(TEXT_COLOR)
    plt.setp(cb.ax.get_yticklabels(), color=TEXT_COLOR)

    if title is None:
        title = f'Parameter Interaction Web — {metric}'
    ax.set_title(title, fontsize=14, fontweight='bold')
    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 7. EVOLUTIONARY RESCUE ASSESSMENT
# ═══════════════════════════════════════════════════════════════════════

def plot_evolutionary_rescue_assessment(
    result: 'CoupledSimResult',
    carrying_capacity: int = 500,
    disease_year: int = 3,
    rescue_threshold: float = 0.5,
    title: str = 'Evolutionary Rescue Assessment',
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Combined panel assessing whether evolutionary rescue is occurring.

    4 subplots:
      [0] Population with rescue threshold line
      [1] Resistance increase rate vs required rate
      [2] Additive variance V_A depletion
      [3] Estimated generations to recovery

    Green/red coding for rescue likelihood.

    Args:
        result: CoupledSimResult.
        carrying_capacity: K for threshold calculation.
        disease_year: Year disease introduced.
        rescue_threshold: Fraction of K considered "rescued".
        title: Super-title.
        save_path: Optional save path.

    Returns:
        matplotlib Figure (16×10).
    """
    fig, axes = plt.subplots(2, 2, figsize=(16, 10))
    fig.patch.set_facecolor(DARK_BG)
    for ax in axes.flat:
        apply_dark_theme(ax=ax)

    years = np.arange(result.n_years)
    pop = result.yearly_pop
    rescue_pop = carrying_capacity * rescue_threshold

    # Determine rescue status
    post_disease = years >= disease_year
    final_pop = pop[-1] if len(pop) > 0 else 0
    rescued = final_pop >= rescue_pop

    status_color = '#2ecc71' if rescued else '#e74c3c'
    status_text = 'RESCUE LIKELY' if rescued else 'RESCUE UNLIKELY'

    # ── [0,0] Population with threshold ──────────────────────────────
    ax = axes[0, 0]
    ax.plot(years, pop, color=status_color, linewidth=2.5, label='Population')
    ax.axhline(rescue_pop, color='#f39c12', linestyle='--', linewidth=1.5,
               label=f'Rescue threshold ({rescue_threshold:.0%} K)')
    ax.axhline(carrying_capacity, color=GRID_COLOR, linestyle=':', alpha=0.5,
               label=f'K={carrying_capacity}')
    ax.axvline(disease_year, color=COMPARTMENT_COLORS['E'], linestyle=':',
               alpha=0.7)
    ax.fill_between(years, pop, alpha=0.1, color=status_color)
    ax.set_xlabel('Year')
    ax.set_ylabel('Population')
    ax.set_title('Population vs Rescue Threshold', fontsize=12, fontweight='bold')
    _legend(ax)

    # ── [0,1] Resistance increase rate ───────────────────────────────
    ax = axes[0, 1]
    if result.yearly_mean_resistance is not None:
        r = result.yearly_mean_resistance
        dr = np.diff(r)
        ax.bar(years[1:], dr, color=[('#2ecc71' if d > 0 else '#e74c3c')
               for d in dr], alpha=0.8, edgecolor='white', linewidth=0.3)
        ax.axhline(0, color=GRID_COLOR, linewidth=0.5)
        # Show cumulative shift
        total_shift = r[-1] - r[disease_year] if disease_year < len(r) else 0
        ax.text(0.95, 0.95, f'Δr̄ = {total_shift:+.4f}',
                transform=ax.transAxes, ha='right', va='top',
                fontsize=11, color=TEXT_COLOR,
                bbox=dict(boxstyle='round', facecolor=DARK_PANEL, alpha=0.8))
    ax.set_xlabel('Year')
    ax.set_ylabel('Δ Mean Resistance / yr')
    ax.set_title('Resistance Increase Rate', fontsize=12, fontweight='bold')

    # ── [1,0] V_A depletion ──────────────────────────────────────────
    ax = axes[1, 0]
    if result.yearly_va is not None:
        va = result.yearly_va
        va_norm = va / np.maximum(va[0], 1e-10)
        color = '#2ecc71' if va_norm[-1] > 0.3 else '#e74c3c'
        ax.plot(years, va_norm, color=color, linewidth=2.5)
        ax.fill_between(years, va_norm, alpha=0.15, color=color)
        ax.axhline(0.3, color='#f39c12', linestyle='--', alpha=0.6,
                   label='Critical V_A (30%)')
        ax.text(0.95, 0.05, f'V_A retention: {va_norm[-1]:.1%}',
                transform=ax.transAxes, ha='right', va='bottom',
                fontsize=11, color=TEXT_COLOR,
                bbox=dict(boxstyle='round', facecolor=DARK_PANEL, alpha=0.8))
    ax.set_xlabel('Year')
    ax.set_ylabel('V_A / V_A₀')
    ax.set_ylim(0, 1.1)
    ax.set_title('Additive Variance Depletion', fontsize=12, fontweight='bold')
    _legend(ax)

    # ── [1,1] Generations to recovery estimate ───────────────────────
    ax = axes[1, 1]
    if result.yearly_mean_resistance is not None:
        r = result.yearly_mean_resistance
        # Estimated generation time for Pycnopodia: ~10 years
        gen_time = 10
        post_years = np.arange(disease_year, result.n_years)
        if len(post_years) > 1 and disease_year < len(r):
            r_post = r[disease_year:]
            if len(r_post) > 1:
                # Δr per generation
                dr_per_gen = (r_post[-1] - r_post[0]) / max(len(r_post) / gen_time, 0.1)
                # Target: r̄ = 0.5 (arbitrary "resistant enough")
                target_r = 0.5
                current_r = r_post[-1]
                gap = target_r - current_r
                if dr_per_gen > 0:
                    gens_needed = gap / dr_per_gen
                    years_needed = gens_needed * gen_time
                else:
                    gens_needed = float('inf')
                    years_needed = float('inf')

                # Timeline bar
                stages = ['Current', 'Target']
                values = [current_r, target_r]
                colors_bar = [status_color, '#f39c12']
                ax.barh([0, 1], values, color=colors_bar, alpha=0.8,
                        edgecolor='white', linewidth=0.5)
                ax.set_yticks([0, 1])
                ax.set_yticklabels(stages, fontsize=11, color=TEXT_COLOR)
                ax.set_xlabel('Mean Resistance')

                est_text = (f'{gens_needed:.0f} gen ({years_needed:.0f} yr)'
                           if np.isfinite(gens_needed) else '∞ (no selection)')
                ax.text(0.95, 0.95, f'Estimate: {est_text}',
                        transform=ax.transAxes, ha='right', va='top',
                        fontsize=11, color=TEXT_COLOR,
                        bbox=dict(boxstyle='round', facecolor=DARK_PANEL,
                                  alpha=0.8))
    ax.set_title('Generations to Recovery', fontsize=12, fontweight='bold')

    # Overall assessment banner
    fig.suptitle(f'{title}  —  {status_text}', fontsize=16, fontweight='bold',
                 color=status_color, y=0.98)
    fig.subplots_adjust(hspace=0.35, wspace=0.3, top=0.92)
    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 8. CONSERVATION SCENARIO MATRIX
# ═══════════════════════════════════════════════════════════════════════

def plot_conservation_scenario_matrix(
    results_dict: Dict[str, 'CoupledSimResult'],
    title: str = 'Conservation Intervention Comparison',
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Grid comparing conservation scenarios.

    Rows: one per scenario (e.g. no intervention, captive breeding, gene flow).
    Columns: population, genetics (resistance), disease deaths.

    Args:
        results_dict: Mapping scenario_name → CoupledSimResult.
        title: Super-title.
        save_path: Optional save path.

    Returns:
        matplotlib Figure.
    """
    scenarios = list(results_dict.keys())
    n_rows = len(scenarios)
    col_names = ['Population', 'Mean Resistance', 'Disease Deaths']

    fig, axes = plt.subplots(n_rows, 3, figsize=(18, 4 * n_rows), sharey='col')
    fig.patch.set_facecolor(DARK_BG)
    if n_rows == 1:
        axes = axes.reshape(1, -1)
    for ax in axes.flat:
        apply_dark_theme(ax=ax)

    scenario_colors = ['#95a5a6', '#2ecc71', '#3498db', '#f39c12', '#e94560']

    for row, (name, res) in enumerate(results_dict.items()):
        years = np.arange(res.n_years)
        color = scenario_colors[row % len(scenario_colors)]

        # Col 0: Population
        ax = axes[row, 0]
        ax.plot(years, res.yearly_pop, color=color, linewidth=2.5)
        ax.fill_between(years, res.yearly_pop, alpha=0.15, color=color)
        ax.set_ylabel(name, fontsize=11, fontweight='bold', rotation=0,
                      labelpad=60, ha='right')
        if row == 0:
            ax.set_title(col_names[0], fontsize=12, fontweight='bold')
        if row == n_rows - 1:
            ax.set_xlabel('Year')

        # Col 1: Resistance
        ax = axes[row, 1]
        if res.yearly_mean_resistance is not None:
            ax.plot(years, res.yearly_mean_resistance, color=color, linewidth=2.5)
        if row == 0:
            ax.set_title(col_names[1], fontsize=12, fontweight='bold')
        if row == n_rows - 1:
            ax.set_xlabel('Year')

        # Col 2: Disease deaths
        ax = axes[row, 2]
        if res.yearly_disease_deaths is not None:
            ax.bar(years, res.yearly_disease_deaths, color=color, alpha=0.8,
                   edgecolor='white', linewidth=0.3)
        if row == 0:
            ax.set_title(col_names[2], fontsize=12, fontweight='bold')
        if row == n_rows - 1:
            ax.set_xlabel('Year')

    fig.suptitle(title, fontsize=16, fontweight='bold', color=TEXT_COLOR, y=0.98)
    fig.subplots_adjust(hspace=0.35, wspace=0.3, top=0.92)
    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 9. MODEL VALIDATION PANEL
# ═══════════════════════════════════════════════════════════════════════

def plot_model_validation_panel(
    result: 'CoupledSimResult',
    observed_data: Optional[Dict[str, Any]] = None,
    disease_year: int = 3,
    title: str = 'Model Validation Panel',
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Model vs observed data comparison (3 panels).

    Panels:
      [0] Mortality fraction (model vs Hamilton 2021 ~90-96%)
      [1] Recovery timeline (model vs first wild sightings)
      [2] Population decline pattern

    Accepts optional observed_data dict with keys:
      'hamilton_mortality': float (0.86-0.96)
      'sighting_year': int (year of first wild post-crash sighting)
      'field_mortality_by_lat': list of (lat, mortality) tuples

    Args:
        result: CoupledSimResult.
        observed_data: Optional observed data for comparison.
        disease_year: Year disease introduced.
        title: Super-title.
        save_path: Optional save path.

    Returns:
        matplotlib Figure.
    """
    if observed_data is None:
        observed_data = {}

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.patch.set_facecolor(DARK_BG)
    for ax in axes.flat:
        apply_dark_theme(ax=ax)

    years = np.arange(result.n_years)
    pop = result.yearly_pop

    # ── [0] Mortality fraction ───────────────────────────────────────
    ax = axes[0]
    if pop[0] > 0:
        min_pop = pop.min()
        model_mort = 1.0 - min_pop / pop[0]
    else:
        model_mort = 0.0

    # Hamilton 2021 range: 86-96%
    hamilton_low = observed_data.get('hamilton_mortality_low', 0.86)
    hamilton_high = observed_data.get('hamilton_mortality_high', 0.96)

    bar_colors = ['#e94560', '#f39c12']
    ax.bar([0], [model_mort], color=bar_colors[0], alpha=0.85, width=0.5,
           edgecolor='white', label='Model')
    ax.axhspan(hamilton_low, hamilton_high, color='#f39c12', alpha=0.2,
               label=f'Hamilton 2021 ({hamilton_low:.0%}–{hamilton_high:.0%})')
    ax.set_xticks([0])
    ax.set_xticklabels(['Peak crash'], color=TEXT_COLOR)
    ax.set_ylabel('Mortality Fraction')
    ax.set_ylim(0, 1.05)
    ax.set_title('Mortality vs Field Data', fontsize=12, fontweight='bold')
    # Status badge
    in_range = hamilton_low <= model_mort <= hamilton_high
    badge = '✓ Match' if in_range else '✗ Mismatch'
    badge_c = '#2ecc71' if in_range else '#e74c3c'
    ax.text(0.95, 0.95, badge, transform=ax.transAxes, ha='right', va='top',
            fontsize=13, color=badge_c, fontweight='bold',
            bbox=dict(boxstyle='round', facecolor=DARK_PANEL, alpha=0.8))
    _legend(ax)

    # ── [1] Recovery timeline ────────────────────────────────────────
    ax = axes[1]
    ax.plot(years, pop, color=ACCENT_COLORS[3], linewidth=2.5)
    ax.fill_between(years, pop, alpha=0.15, color=ACCENT_COLORS[3])
    # Mark nadir
    nadir_yr = np.argmin(pop)
    ax.scatter([nadir_yr], [pop[nadir_yr]], color='#e74c3c', s=100, zorder=5,
               marker='v', label=f'Nadir (yr {nadir_yr})')
    # Mark first sighting if provided
    sighting_yr = observed_data.get('sighting_year', None)
    if sighting_yr is not None and sighting_yr < len(years):
        ax.axvline(sighting_yr, color='#2ecc71', linestyle='--', linewidth=2,
                   label=f'Wild sighting (yr {sighting_yr})')
    ax.set_xlabel('Year')
    ax.set_ylabel('Population')
    ax.set_title('Recovery Timeline', fontsize=12, fontweight='bold')
    _legend(ax)

    # ── [2] Population decline pattern ───────────────────────────────
    ax = axes[2]
    if pop[0] > 0:
        pop_frac = pop / pop[0]
    else:
        pop_frac = np.ones_like(pop)
    ax.plot(years, pop_frac, color='#e94560', linewidth=2.5, label='Model')
    # Canonical SSWD timeline: crash within 2-3 years
    crash_window = [disease_year, disease_year + 3]
    ax.axvspan(crash_window[0], min(crash_window[1], years[-1]),
               color='#e74c3c', alpha=0.1, label='Expected crash window')
    ax.axhline(0.1, color='#e74c3c', linestyle=':', alpha=0.5,
               label='90% decline')
    ax.set_xlabel('Year')
    ax.set_ylabel('Population (fraction of initial)')
    ax.set_ylim(0, 1.1)
    ax.set_title('Decline Pattern', fontsize=12, fontweight='bold')
    _legend(ax)

    fig.suptitle(title, fontsize=16, fontweight='bold', color=TEXT_COLOR, y=0.98)
    fig.subplots_adjust(hspace=0.35, wspace=0.3, top=0.88)
    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 10. PARAMETER SPACE EXPLORATION — MORRIS VS SOBOL
# ═══════════════════════════════════════════════════════════════════════

def plot_parameter_space_exploration(
    morris_results: Dict,
    sobol_indices: Dict,
    metric: str = 'pop_crash_pct',
    title: Optional[str] = None,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Morris μ* vs Sobol ST scatter — each point is a parameter.

    Points near origin = unimportant.
    Points far from diagonal = strong nonlinear interactions.

    Args:
        morris_results: Dict from morris_screening.json.
        sobol_indices: Dict from sobol_indices.json.
        metric: Which metric to compare.
        title: Optional title override.
        save_path: Optional save path.

    Returns:
        matplotlib Figure.
    """
    morris_params = morris_results['param_names']
    sobol_params = sobol_indices['param_names']

    # Get Morris mu_star for the metric
    morris_data = morris_results['morris_results'][metric]
    mu_star = np.array(morris_data['mu_star'])

    # Get Sobol ST for same metric
    sobol_data = sobol_indices['indices'][metric]
    ST = np.array(sobol_data['ST'])

    # Match parameters (should be in same order for our data)
    assert len(morris_params) == len(sobol_params), \
        f"Param count mismatch: Morris={len(morris_params)}, Sobol={len(sobol_params)}"

    fig, ax = dark_figure(figsize=(10, 10))

    # Scatter with category colors
    for i, param in enumerate(morris_params):
        cat = _param_category(param)
        color = SA_CATEGORY_COLORS.get(cat, '#95a5a6')
        ax.scatter(mu_star[i], ST[i], color=color, s=100, alpha=0.85,
                   edgecolors='white', linewidth=0.5, zorder=5)
        # Label top parameters
        if mu_star[i] > np.percentile(mu_star, 75) or ST[i] > np.percentile(ST, 75):
            ax.annotate(param.split('.')[-1],
                        (mu_star[i], ST[i]),
                        textcoords='offset points', xytext=(6, 6),
                        fontsize=8, color=TEXT_COLOR, alpha=0.9)

    # Diagonal line (proportional relationship)
    max_mu = mu_star.max()
    max_st = ST.max()
    if max_mu > 0 and max_st > 0:
        scale = max_st / max_mu
        diag_x = np.linspace(0, max_mu * 1.1, 100)
        ax.plot(diag_x, diag_x * scale, color=GRID_COLOR, linestyle='--',
                alpha=0.4, label='Proportional')

    ax.set_xlabel('Morris μ* (elementary effect)', fontsize=12)
    ax.set_ylabel('Sobol ST (total-order)', fontsize=12)

    # Category legend
    handles = []
    cats_present = set(_param_category(p) for p in morris_params)
    for cat in sorted(cats_present):
        handles.append(mpatches.Patch(
            color=SA_CATEGORY_COLORS.get(cat, '#95a5a6'),
            label=cat.replace('_', ' ').title()))
    ax.legend(handles=handles, loc='upper left',
              facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
              labelcolor=TEXT_COLOR, fontsize=10)

    # Quadrant annotations
    ax.text(0.02, 0.98, 'Low μ*, High ST\n→ Strong interactions',
            transform=ax.transAxes, ha='left', va='top',
            fontsize=9, color='#f39c12', alpha=0.7,
            fontstyle='italic')
    ax.text(0.98, 0.02, 'High μ*, Low ST\n→ Linear effects',
            transform=ax.transAxes, ha='right', va='bottom',
            fontsize=9, color='#3498db', alpha=0.7,
            fontstyle='italic')

    if title is None:
        title = f'Parameter Space: Morris vs Sobol — {metric}'
    ax.set_title(title, fontsize=14, fontweight='bold')
    if save_path:
        save_figure(fig, save_path)
    return fig
