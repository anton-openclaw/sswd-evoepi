"""Population & demographics visualizations for SSWD-EvoEpi.

Every function:
  - Accepts model results (CoupledSimResult or SpatialSimResult) as input
  - Returns a matplotlib Figure
  - Has an optional ``save_path`` parameter (saves PNG when given)
  - Uses the shared dark theme from ``sswd_evoepi.viz.style``

matplotlib backend is forced to Agg (no display) on import.
"""

from __future__ import annotations

import matplotlib
matplotlib.use('Agg')

from pathlib import Path
from typing import Optional, TYPE_CHECKING

import matplotlib.pyplot as plt
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
# 1. POPULATION TRAJECTORY
# ═══════════════════════════════════════════════════════════════════════

def plot_population_trajectory(
    result: 'CoupledSimResult',
    carrying_capacity: int = 500,
    disease_year: int = 3,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Total population over time with carrying-capacity line and disease marker.

    Args:
        result: CoupledSimResult with yearly_pop.
        carrying_capacity: K to show as dashed line.
        disease_year: Year disease was introduced (vertical marker).
        save_path: Optional path to save the figure.

    Returns:
        matplotlib Figure.
    """
    years = np.arange(result.n_years)
    fig, ax = dark_figure()

    ax.plot(years, result.yearly_pop, color=ACCENT_COLORS[0], linewidth=2.5,
            label='Total population', zorder=3)
    ax.axhline(carrying_capacity, color=ACCENT_COLORS[3], linestyle='--',
               linewidth=1.5, alpha=0.7, label=f'K = {carrying_capacity}')
    ax.axvline(disease_year, color=ACCENT_COLORS[4], linestyle=':',
               linewidth=1.5, alpha=0.8, label=f'Disease intro (year {disease_year})')

    ax.fill_between(years, result.yearly_pop, alpha=0.15, color=ACCENT_COLORS[0])

    ax.set_xlabel('Year', fontsize=12)
    ax.set_ylabel('Population size', fontsize=12)
    ax.set_title('Population Trajectory', fontsize=14, fontweight='bold')
    ax.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
              labelcolor=TEXT_COLOR, fontsize=10)
    ax.set_xlim(0, result.n_years - 1)
    ax.set_ylim(bottom=0)

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 2. STAGE COMPOSITION
# ═══════════════════════════════════════════════════════════════════════

def plot_stage_composition(
    result: 'CoupledSimResult',
    agents_history: Optional[list] = None,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Stacked area chart of life stages over time.

    If ``agents_history`` (list of agent arrays, one per year) is provided,
    counts are computed from the arrays. Otherwise, a synthetic decomposition
    is estimated from yearly_pop and yearly_adults.

    Args:
        result: CoupledSimResult.
        agents_history: Optional list of agent structured arrays per year.
        save_path: Optional save path.

    Returns:
        matplotlib Figure.
    """
    from sswd_evoepi.types import Stage

    years = np.arange(result.n_years)
    fig, ax = dark_figure()

    if agents_history is not None and len(agents_history) == result.n_years:
        settlers = np.zeros(result.n_years, dtype=int)
        juveniles = np.zeros(result.n_years, dtype=int)
        subadults = np.zeros(result.n_years, dtype=int)
        adults = np.zeros(result.n_years, dtype=int)
        for y, agents in enumerate(agents_history):
            alive = agents['alive'].astype(bool)
            stages = agents['stage'][alive]
            settlers[y] = int(np.sum(stages == Stage.SETTLER))
            juveniles[y] = int(np.sum(stages == Stage.JUVENILE))
            subadults[y] = int(np.sum(stages == Stage.SUBADULT))
            adults[y] = int(np.sum(stages == Stage.ADULT))
    else:
        # Synthetic decomposition from available data
        total = result.yearly_pop.astype(float)
        adults_arr = result.yearly_adults.astype(float)
        non_adult = np.maximum(total - adults_arr, 0)
        # Approximate proportions: 5% settler, 25% juvenile, 20% subadult among non-adults
        settlers = (non_adult * 0.10).astype(int)
        juveniles = (non_adult * 0.50).astype(int)
        subadults = (non_adult * 0.40).astype(int)
        adults = adults_arr.astype(int)

    stage_order = ['settler', 'juvenile', 'subadult', 'adult']
    stage_data = [settlers, juveniles, subadults, adults]
    colors = [STAGE_COLORS[s] for s in stage_order]

    ax.stackplot(years, *stage_data, labels=[s.capitalize() for s in stage_order],
                 colors=colors, alpha=0.85)
    ax.set_xlabel('Year', fontsize=12)
    ax.set_ylabel('Population', fontsize=12)
    ax.set_title('Stage Composition Over Time', fontsize=14, fontweight='bold')
    ax.legend(loc='upper right', facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
              labelcolor=TEXT_COLOR, fontsize=10)
    ax.set_xlim(0, result.n_years - 1)
    ax.set_ylim(bottom=0)

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 3. CAUSE OF DEATH BREAKDOWN
# ═══════════════════════════════════════════════════════════════════════

def plot_cause_of_death_breakdown(
    result: 'CoupledSimResult',
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Stacked bar chart of annual deaths by cause.

    Uses cause_of_death tracking: disease, natural mortality, senescence.

    Args:
        result: CoupledSimResult with yearly death arrays.
        save_path: Optional save path.

    Returns:
        matplotlib Figure.
    """
    years = np.arange(result.n_years)
    fig, ax = dark_figure()

    disease = result.yearly_disease_deaths
    senescence = result.yearly_senescence_deaths if result.yearly_senescence_deaths is not None else np.zeros(result.n_years)
    natural = result.yearly_natural_deaths - senescence  # natural excludes senescence

    bar_width = 0.7
    ax.bar(years, natural, bar_width, label='Natural mortality',
           color=DEATH_COLORS['natural'], alpha=0.9)
    ax.bar(years, senescence, bar_width, bottom=natural,
           label='Senescence', color=DEATH_COLORS['senescence'], alpha=0.9)
    ax.bar(years, disease, bar_width, bottom=natural + senescence,
           label='Disease (SSWD)', color=DEATH_COLORS['disease'], alpha=0.9)

    ax.set_xlabel('Year', fontsize=12)
    ax.set_ylabel('Deaths', fontsize=12)
    ax.set_title('Cause of Death Breakdown', fontsize=14, fontweight='bold')
    ax.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
              labelcolor=TEXT_COLOR, fontsize=10)
    ax.set_xlim(-0.5, result.n_years - 0.5)
    ax.set_ylim(bottom=0)

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 4. AGE-SIZE PYRAMID
# ═══════════════════════════════════════════════════════════════════════

def plot_age_size_pyramid(
    agents: np.ndarray,
    title: str = 'Population Pyramid',
    age_bins: int = 15,
    size_bins: int = 15,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Population pyramid: age distribution (left) and size distribution (right).

    Females extend left, males extend right.

    Args:
        agents: Structured array with AGENT_DTYPE fields (alive, sex, age, size).
        title: Figure title.
        age_bins: Number of age bins.
        size_bins: Number of size bins.
        save_path: Optional save path.

    Returns:
        matplotlib Figure.
    """
    fig, (ax_age, ax_size) = dark_figure(nrows=1, ncols=2, figsize=(14, 7))

    alive = agents['alive'].astype(bool)
    female = alive & (agents['sex'] == 0)
    male = alive & (agents['sex'] == 1)

    ages = agents['age']
    sizes = agents['size']

    # --- Age pyramid ---
    max_age = max(float(ages[alive].max()), 1.0) if alive.any() else 30.0
    age_edges = np.linspace(0, max_age, age_bins + 1)
    f_age_counts, _ = np.histogram(ages[female], bins=age_edges)
    m_age_counts, _ = np.histogram(ages[male], bins=age_edges)

    bin_centers = (age_edges[:-1] + age_edges[1:]) / 2
    bar_h = age_edges[1] - age_edges[0]

    ax_age.barh(bin_centers, -f_age_counts, height=bar_h * 0.85,
                color=ACCENT_COLORS[0], alpha=0.85, label='Female')
    ax_age.barh(bin_centers, m_age_counts, height=bar_h * 0.85,
                color=ACCENT_COLORS[5], alpha=0.85, label='Male')

    max_count = max(f_age_counts.max(), m_age_counts.max(), 1)
    ax_age.set_xlim(-max_count * 1.15, max_count * 1.15)
    ax_age.set_xlabel('Count', fontsize=11)
    ax_age.set_ylabel('Age (years)', fontsize=11)
    ax_age.set_title('Age Distribution', fontsize=12, fontweight='bold')
    ax_age.axvline(0, color=GRID_COLOR, linewidth=0.8)
    ax_age.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
                  labelcolor=TEXT_COLOR, fontsize=9)

    # --- Size pyramid ---
    max_size = max(float(sizes[alive].max()), 10.0) if alive.any() else 1000.0
    size_edges = np.linspace(0, max_size, size_bins + 1)
    f_size_counts, _ = np.histogram(sizes[female], bins=size_edges)
    m_size_counts, _ = np.histogram(sizes[male], bins=size_edges)

    s_centers = (size_edges[:-1] + size_edges[1:]) / 2
    s_bar_h = size_edges[1] - size_edges[0]

    ax_size.barh(s_centers, -f_size_counts, height=s_bar_h * 0.85,
                 color=ACCENT_COLORS[0], alpha=0.85, label='Female')
    ax_size.barh(s_centers, m_size_counts, height=s_bar_h * 0.85,
                 color=ACCENT_COLORS[5], alpha=0.85, label='Male')

    max_s = max(f_size_counts.max(), m_size_counts.max(), 1)
    ax_size.set_xlim(-max_s * 1.15, max_s * 1.15)
    ax_size.set_xlabel('Count', fontsize=11)
    ax_size.set_ylabel('Size (mm)', fontsize=11)
    ax_size.set_title('Size Distribution', fontsize=12, fontweight='bold')
    ax_size.axvline(0, color=GRID_COLOR, linewidth=0.8)
    ax_size.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
                   labelcolor=TEXT_COLOR, fontsize=9)

    fig.suptitle(title, fontsize=15, fontweight='bold', color=TEXT_COLOR, y=1.02)

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 5. POPULATION HEATMAP (SPATIAL)
# ═══════════════════════════════════════════════════════════════════════

def plot_population_heatmap(
    spatial_result: 'SpatialSimResult',
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Heatmap: nodes (y) × years (x), coloured by population fraction of K.

    Args:
        spatial_result: SpatialSimResult with per-node yearly_pop and node_K.
        save_path: Optional save path.

    Returns:
        matplotlib Figure.
    """
    n_nodes = spatial_result.n_nodes
    n_years = spatial_result.n_years
    names = spatial_result.node_names or [f'Node {i}' for i in range(n_nodes)]
    K = spatial_result.node_K

    # Fraction of K
    frac = np.zeros((n_nodes, n_years))
    for i in range(n_nodes):
        frac[i] = spatial_result.yearly_pop[i] / max(K[i], 1)

    fig, ax = dark_figure(figsize=(max(10, n_years * 0.5), max(4, n_nodes * 0.8)))

    # Use a diverging colormap anchored at 1.0
    from matplotlib.colors import TwoSlopeNorm
    vmax = max(float(frac.max()), 1.05)
    norm = TwoSlopeNorm(vmin=0, vcenter=0.5, vmax=vmax)

    im = ax.imshow(frac, aspect='auto', cmap='RdYlGn', norm=norm,
                   interpolation='nearest')

    ax.set_xticks(np.arange(n_years))
    ax.set_xticklabels(np.arange(n_years), fontsize=9, color=TEXT_COLOR)
    ax.set_yticks(np.arange(n_nodes))
    ax.set_yticklabels(names, fontsize=10, color=TEXT_COLOR)
    ax.set_xlabel('Year', fontsize=12)
    ax.set_title('Population as Fraction of K', fontsize=14, fontweight='bold')

    # Annotate cells
    for i in range(n_nodes):
        for j in range(n_years):
            val = frac[i, j]
            txt_color = 'white' if val < 0.3 else 'black'
            ax.text(j, i, f'{val:.2f}', ha='center', va='center',
                    fontsize=7, color=txt_color, fontweight='bold')

    cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cbar.ax.tick_params(colors=TEXT_COLOR)
    cbar.set_label('N / K', color=TEXT_COLOR, fontsize=11)

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 6. RECRUITMENT TIMESERIES
# ═══════════════════════════════════════════════════════════════════════

def plot_recruitment_timeseries(
    result: 'CoupledSimResult',
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Annual recruitment (settlers) with an overlay of reproductive adults.

    Args:
        result: CoupledSimResult with yearly_recruits and yearly_adults.
        save_path: Optional save path.

    Returns:
        matplotlib Figure.
    """
    years = np.arange(result.n_years)
    fig, ax1 = dark_figure()

    ax1.bar(years, result.yearly_recruits, color=STAGE_COLORS['settler'],
            alpha=0.85, label='Recruits (settlers)', zorder=2)
    ax1.set_xlabel('Year', fontsize=12)
    ax1.set_ylabel('Annual recruits', fontsize=12, color=STAGE_COLORS['settler'])
    ax1.tick_params(axis='y', colors=STAGE_COLORS['settler'])

    # Overlay: reproductive adults on secondary y-axis
    ax2 = ax1.twinx()
    ax2.set_facecolor('none')
    ax2.plot(years, result.yearly_adults, color=STAGE_COLORS['adult'],
             linewidth=2.5, marker='o', markersize=4, label='Adults',
             zorder=3)
    ax2.set_ylabel('Reproductive adults', fontsize=12,
                    color=STAGE_COLORS['adult'])
    ax2.tick_params(axis='y', colors=STAGE_COLORS['adult'])
    for spine in ax2.spines.values():
        spine.set_color(GRID_COLOR)

    ax1.set_title('Recruitment & Reproductive Adults', fontsize=14,
                   fontweight='bold')

    # Combined legend
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2,
               facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
               labelcolor=TEXT_COLOR, fontsize=10, loc='upper right')
    ax1.set_xlim(-0.5, result.n_years - 0.5)

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 7. SURVIVAL CURVES
# ═══════════════════════════════════════════════════════════════════════

def plot_survival_curves(
    result: 'CoupledSimResult',
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Kaplan-Meier-style survival curves using annual survival probabilities.

    Shows the expected fraction of a cohort surviving over time for each
    life stage, given the stage-specific annual survival from the model's
    ANNUAL_SURVIVAL vector plus observed disease mortality.

    Args:
        result: CoupledSimResult.
        save_path: Optional save path.

    Returns:
        matplotlib Figure.
    """
    from sswd_evoepi.types import ANNUAL_SURVIVAL

    fig, ax = dark_figure()

    max_years = min(result.n_years, 30)
    years = np.arange(max_years + 1)

    # Theoretical survival curves from annual survival rates
    stage_names = ['Settler', 'Juvenile', 'Subadult', 'Adult']
    stage_indices = [1, 2, 3, 4]
    colors = [STAGE_COLORS['settler'], STAGE_COLORS['juvenile'],
              STAGE_COLORS['subadult'], STAGE_COLORS['adult']]

    for name, si, color in zip(stage_names, stage_indices, colors):
        annual_s = float(ANNUAL_SURVIVAL[si])
        survival = annual_s ** years
        ax.plot(years, survival, color=color, linewidth=2.5, label=name)

    # Observed population survival fraction (total)
    if result.yearly_pop is not None and result.initial_pop > 0:
        obs_frac = result.yearly_pop[:max_years] / result.initial_pop
        ax.plot(np.arange(len(obs_frac)), obs_frac, color=TEXT_COLOR,
                linewidth=2, linestyle='--', alpha=0.7,
                label='Observed (total)')

    ax.set_xlabel('Years from cohort start', fontsize=12)
    ax.set_ylabel('Fraction surviving', fontsize=12)
    ax.set_title('Stage-Specific Survival Curves', fontsize=14, fontweight='bold')
    ax.set_ylim(0, 1.05)
    ax.set_xlim(0, max_years)
    ax.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
              labelcolor=TEXT_COLOR, fontsize=10)

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 8. SEX RATIO OVER TIME
# ═══════════════════════════════════════════════════════════════════════

def plot_sex_ratio_over_time(
    result: 'CoupledSimResult',
    agents_history: Optional[list] = None,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Proportion female over time; should be ~0.5 unless disease is sex-biased.

    If ``agents_history`` is provided (list of agent arrays per year), actual
    sex ratios are computed.  Otherwise a synthetic estimate from population
    size with noise is shown (less informative but still useful for layout).

    Args:
        result: CoupledSimResult.
        agents_history: Optional list of agent arrays, one per year.
        save_path: Optional save path.

    Returns:
        matplotlib Figure.
    """
    years = np.arange(result.n_years)
    fig, ax = dark_figure()

    if agents_history is not None and len(agents_history) == result.n_years:
        sex_ratio = np.zeros(result.n_years)
        for y, agents in enumerate(agents_history):
            alive = agents['alive'].astype(bool)
            n_alive = int(np.sum(alive))
            if n_alive > 0:
                n_female = int(np.sum(alive & (agents['sex'] == 0)))
                sex_ratio[y] = n_female / n_alive
            else:
                sex_ratio[y] = np.nan
    else:
        # Synthetic: ~0.5 with small stochastic drift proportional to √N
        rng = np.random.default_rng(12345)
        pop = result.yearly_pop.astype(float)
        noise = rng.normal(0, 1, size=result.n_years)
        sex_ratio = 0.5 + noise / (2 * np.sqrt(np.maximum(pop, 1)))
        sex_ratio = np.clip(sex_ratio, 0, 1)
        # When pop is 0, mark as nan
        sex_ratio[pop == 0] = np.nan

    ax.plot(years, sex_ratio, color=ACCENT_COLORS[0], linewidth=2.5, zorder=3)
    ax.axhline(0.5, color=ACCENT_COLORS[3], linestyle='--', linewidth=1.5,
               alpha=0.7, label='Expected (0.5)')
    ax.fill_between(years, 0.45, 0.55, alpha=0.1, color=ACCENT_COLORS[3],
                     label='±5% band')

    ax.set_xlabel('Year', fontsize=12)
    ax.set_ylabel('Proportion female', fontsize=12)
    ax.set_title('Sex Ratio Over Time', fontsize=14, fontweight='bold')
    ax.set_ylim(0.3, 0.7)
    ax.set_xlim(0, result.n_years - 1)
    ax.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
              labelcolor=TEXT_COLOR, fontsize=10)

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 9. DENSITY DEPENDENCE
# ═══════════════════════════════════════════════════════════════════════

def plot_density_dependence(
    result: 'CoupledSimResult',
    carrying_capacity: int = 500,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Population size vs per-capita growth rate (should show neg. density dep.).

    Per-capita growth rate: r = ln(N_{t+1} / N_t).
    Points should cluster on a declining curve consistent with Beverton-Holt.

    Args:
        result: CoupledSimResult with yearly_pop.
        carrying_capacity: K for reference line.
        save_path: Optional save path.

    Returns:
        matplotlib Figure.
    """
    fig, ax = dark_figure()

    pop = result.yearly_pop.astype(float)
    # Per-capita growth rate
    N_t = pop[:-1]
    N_t1 = pop[1:]

    valid = N_t > 0
    r_pc = np.full(len(N_t), np.nan)
    r_pc[valid] = np.log(np.maximum(N_t1[valid], 1) / N_t[valid])

    ax.scatter(N_t, r_pc, color=ACCENT_COLORS[0], s=60, edgecolors='white',
               linewidth=0.5, alpha=0.85, zorder=3)

    # Annotate year numbers
    for i in range(len(N_t)):
        if not np.isnan(r_pc[i]):
            ax.annotate(str(i), (N_t[i], r_pc[i]), fontsize=7,
                        color=TEXT_COLOR, alpha=0.6, textcoords='offset points',
                        xytext=(5, 5))

    ax.axhline(0, color=GRID_COLOR, linewidth=1, linestyle='-')
    ax.axvline(carrying_capacity, color=ACCENT_COLORS[3], linestyle='--',
               linewidth=1.5, alpha=0.7, label=f'K = {carrying_capacity}')

    # Theoretical Beverton-Holt curve
    N_range = np.linspace(1, carrying_capacity * 1.5, 200)
    # BH: r = ln(K/N) approximately when near equilibrium
    # More precisely: R(N) = α / (1 + β·N) where at K, R=1 → r=0
    # Simple: r ≈ r_max * (1 - N/K) for logistic
    r_max = 0.5  # approximate from BH parameters
    r_theory = r_max * (1 - N_range / carrying_capacity)
    ax.plot(N_range, r_theory, color=ACCENT_COLORS[3], linewidth=1.5,
            alpha=0.5, linestyle=':', label='Logistic approx.')

    ax.set_xlabel('Population size (N)', fontsize=12)
    ax.set_ylabel('Per-capita growth rate (r)', fontsize=12)
    ax.set_title('Density Dependence', fontsize=14, fontweight='bold')
    ax.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
              labelcolor=TEXT_COLOR, fontsize=10)

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 10. NODE COMPARISON BARS (SPATIAL)
# ═══════════════════════════════════════════════════════════════════════

def plot_node_comparison_bars(
    spatial_result: 'SpatialSimResult',
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Grouped bar chart: initial vs final population across nodes, coloured by subregion.

    Args:
        spatial_result: SpatialSimResult.
        save_path: Optional save path.

    Returns:
        matplotlib Figure.
    """
    n_nodes = spatial_result.n_nodes
    names = spatial_result.node_names or [f'Node {i}' for i in range(n_nodes)]
    K = spatial_result.node_K

    initial = spatial_result.yearly_pop[:, 0]
    final = spatial_result.yearly_pop[:, -1]

    fig, ax = dark_figure(figsize=(max(8, n_nodes * 1.5), 6))

    x = np.arange(n_nodes)
    w = 0.3

    ax.bar(x - w / 2, initial, w, color=ACCENT_COLORS[5], alpha=0.85,
           label='Initial', edgecolor='white', linewidth=0.5)
    ax.bar(x + w / 2, final, w, color=ACCENT_COLORS[0], alpha=0.85,
           label='Final', edgecolor='white', linewidth=0.5)

    # Show K as markers
    ax.scatter(x, K, color=ACCENT_COLORS[3], marker='_', s=200,
               linewidth=2.5, zorder=5, label='K')

    ax.set_xticks(x)
    ax.set_xticklabels(names, fontsize=11, color=TEXT_COLOR)
    ax.set_xlabel('Node', fontsize=12)
    ax.set_ylabel('Population', fontsize=12)
    ax.set_title('Node Population: Initial vs Final', fontsize=14,
                  fontweight='bold')
    ax.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
              labelcolor=TEXT_COLOR, fontsize=10)
    ax.set_ylim(bottom=0)

    if save_path:
        save_figure(fig, save_path)
    return fig
