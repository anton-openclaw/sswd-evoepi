"""Host genetics & evolution visualizations for SSWD-EvoEpi.

Every function:
  - Accepts model results (CoupledSimResult or SpatialSimResult) as input
  - Returns a matplotlib Figure
  - Has an optional ``save_path`` parameter (saves PNG when given)
  - Uses the shared dark theme from ``sswd_evoepi.viz.style``

matplotlib backend is forced to Agg (no display) on import.

Resistance color gradient:
  low=#e74c3c (red) → mid=#f39c12 (orange) → high=#2ecc71 (green)

Authors: Anton 🔬 & Willem Weertman
"""

from __future__ import annotations

import matplotlib
matplotlib.use('Agg')

from pathlib import Path
from typing import Dict, List, Optional, Sequence, TYPE_CHECKING

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
import numpy as np
from scipy import stats as sp_stats

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
# RESISTANCE COLOUR MAP
# ═══════════════════════════════════════════════════════════════════════

RESISTANCE_LOW = '#e74c3c'   # red
RESISTANCE_MID = '#f39c12'   # orange
RESISTANCE_HIGH = '#2ecc71'  # green

_RESISTANCE_CMAP = mcolors.LinearSegmentedColormap.from_list(
    'resistance', [RESISTANCE_LOW, RESISTANCE_MID, RESISTANCE_HIGH], N=256,
)


def _resistance_color(value: float, vmin: float = 0.0,
                      vmax: float = 0.4) -> tuple:
    """Map a resistance value to an RGBA tuple."""
    norm = max(0.0, min(1.0, (value - vmin) / (vmax - vmin)))
    return _RESISTANCE_CMAP(norm)


# ═══════════════════════════════════════════════════════════════════════
# 1. RESISTANCE TRAJECTORY
# ═══════════════════════════════════════════════════════════════════════

def plot_resistance_trajectory(
    result: 'CoupledSimResult',
    disease_year: int = 3,
    std_band: Optional[np.ndarray] = None,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Mean resistance over time with ±1 SD shaded band.

    THE key evolutionary plot — shows whether selection is shifting
    resistance upward following disease introduction.

    Args:
        result: CoupledSimResult with yearly_mean_resistance.
        disease_year: Year disease was introduced (vertical marker).
        std_band: Optional (n_years,) array of resistance SD.  If None,
            a rough estimate of SD is computed from V_A if available.
        save_path: Path to save figure.

    Returns:
        matplotlib Figure.
    """
    years = np.arange(result.n_years)
    mean_r = result.yearly_mean_resistance

    fig, ax = dark_figure()

    # Compute SD band from V_A if available and std_band not provided
    if std_band is None and result.yearly_va is not None:
        # V_P ≈ V_A + V_E; conservative estimate SD ≈ sqrt(2 * V_A)
        std_band = np.sqrt(np.maximum(result.yearly_va * 2.0, 0.0))

    # Colour the line by resistance level using a line collection
    for i in range(len(years) - 1):
        color = _resistance_color(mean_r[i])
        ax.plot(years[i:i + 2], mean_r[i:i + 2], color=color,
                linewidth=2.5, solid_capstyle='round')

    # SD band
    if std_band is not None:
        ax.fill_between(years, mean_r - std_band, mean_r + std_band,
                        color=RESISTANCE_MID, alpha=0.15,
                        label='± 1 SD')

    # Disease marker
    ax.axvline(disease_year, color='#e94560', linestyle=':', linewidth=1.5,
               alpha=0.8, label=f'Disease intro (yr {disease_year})')

    # Annotations
    if len(mean_r) > disease_year:
        ax.annotate(
            f'Pre: {mean_r[disease_year]:.3f}',
            xy=(disease_year, mean_r[disease_year]),
            xytext=(disease_year + 1.5, mean_r[disease_year] + 0.02),
            color=TEXT_COLOR, fontsize=9,
            arrowprops=dict(arrowstyle='->', color=TEXT_COLOR, lw=0.8),
        )
    if len(mean_r) > 0:
        ax.annotate(
            f'Final: {mean_r[-1]:.3f}',
            xy=(years[-1], mean_r[-1]),
            xytext=(years[-1] - 3, mean_r[-1] + 0.02),
            color=TEXT_COLOR, fontsize=9,
            arrowprops=dict(arrowstyle='->', color=TEXT_COLOR, lw=0.8),
        )

    ax.set_xlabel('Year', fontsize=12)
    ax.set_ylabel('Mean resistance ($\\bar{r}_i$)', fontsize=12)
    ax.set_title('Host Resistance Evolution', fontsize=14, fontweight='bold')
    ax.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
              labelcolor=TEXT_COLOR, fontsize=10, loc='upper left')
    ax.set_xlim(0, result.n_years - 1)

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 2. RESISTANCE DISTRIBUTION (RIDGE / OVERLAID HISTOGRAMS)
# ═══════════════════════════════════════════════════════════════════════

def plot_resistance_distribution(
    snapshots: Dict[int, np.ndarray],
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Overlaid KDE curves of resistance at multiple timepoints.

    Shows the distribution shifting rightward under selection.

    Args:
        snapshots: Mapping {year: 1-D resistance array}.
        save_path: Path to save figure.

    Returns:
        matplotlib Figure.
    """
    fig, ax = dark_figure()

    sorted_years = sorted(snapshots.keys())
    n = len(sorted_years)
    cmap = plt.cm.viridis

    for idx, yr in enumerate(sorted_years):
        r_vals = snapshots[yr]
        if len(r_vals) < 5:
            continue
        color = cmap(idx / max(n - 1, 1))
        kde_x = np.linspace(0, max(0.6, r_vals.max() * 1.3), 300)
        kde = sp_stats.gaussian_kde(r_vals, bw_method=0.05)
        ax.fill_between(kde_x, kde(kde_x), alpha=0.20, color=color)
        ax.plot(kde_x, kde(kde_x), color=color, linewidth=2,
                label=f'Year {yr} (n={len(r_vals)}, '
                      f'$\\mu$={np.mean(r_vals):.3f})')

    ax.set_xlabel('Resistance ($r_i$)', fontsize=12)
    ax.set_ylabel('Density', fontsize=12)
    ax.set_title('Resistance Distribution Shift Under Selection',
                 fontsize=14, fontweight='bold')
    ax.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
              labelcolor=TEXT_COLOR, fontsize=9, loc='upper right')
    ax.set_xlim(left=0)

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 3. ALLELE FREQUENCY SPAGHETTI PLOT
# ═══════════════════════════════════════════════════════════════════════

def plot_allele_frequency_spaghetti(
    yearly_allele_freq: np.ndarray,
    effect_sizes: np.ndarray,
    disease_year: int = 3,
    n_highlight: int = 5,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Per-locus allele frequency trajectories over time.

    Each of 51 loci as a thin line, coloured by effect size.
    Top ``n_highlight`` loci by effect size are highlighted.

    Args:
        yearly_allele_freq: (n_years, N_ADDITIVE) allele freq per year.
        effect_sizes: (N_ADDITIVE,) float64.
        disease_year: Year disease was introduced (vertical marker).
        n_highlight: Number of top loci to highlight.
        save_path: Path to save figure.

    Returns:
        matplotlib Figure.
    """
    n_years, n_loci = yearly_allele_freq.shape
    years = np.arange(n_years)
    fig, ax = dark_figure()

    # Colour by effect size
    e_min, e_max = effect_sizes.min(), effect_sizes.max()
    norm = mcolors.Normalize(vmin=e_min, vmax=e_max)

    # Sort by effect size for layering (small effects behind)
    order = np.argsort(effect_sizes)  # ascending
    top_idx = set(np.argsort(-effect_sizes)[:n_highlight])

    for rank, l_idx in enumerate(order):
        if l_idx in top_idx:
            continue  # draw highlighted loci on top
        color = _RESISTANCE_CMAP(norm(effect_sizes[l_idx]))
        ax.plot(years, yearly_allele_freq[:, l_idx],
                color=color, alpha=0.25, linewidth=0.7)

    # Highlight top loci
    top_sorted = sorted(top_idx, key=lambda i: -effect_sizes[i])
    for rank, l_idx in enumerate(top_sorted):
        color = _RESISTANCE_CMAP(norm(effect_sizes[l_idx]))
        ax.plot(years, yearly_allele_freq[:, l_idx],
                color=color, linewidth=2.0, alpha=0.9,
                label=f'Locus {l_idx} (e={effect_sizes[l_idx]:.3f})')

    ax.axvline(disease_year, color='#e94560', linestyle=':', linewidth=1.5,
               alpha=0.8)

    # Colourbar
    sm = plt.cm.ScalarMappable(cmap=_RESISTANCE_CMAP, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label('Effect size', color=TEXT_COLOR, fontsize=10)
    cbar.ax.yaxis.set_tick_params(color=TEXT_COLOR)
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color=TEXT_COLOR)

    ax.set_xlabel('Year', fontsize=12)
    ax.set_ylabel('Allele frequency ($q_l$)', fontsize=12)
    ax.set_title('Per-Locus Allele Frequency Trajectories',
                 fontsize=14, fontweight='bold')
    ax.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
              labelcolor=TEXT_COLOR, fontsize=8, loc='upper left')
    ax.set_xlim(0, n_years - 1)
    ax.set_ylim(0, 1)

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 4. ADDITIVE VARIANCE OVER TIME
# ═══════════════════════════════════════════════════════════════════════

def plot_additive_variance_over_time(
    result: 'CoupledSimResult',
    disease_year: int = 3,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """V_A (additive genetic variance) over time.

    Selection erodes V_A, but drift and mutation can maintain it.
    Critical for judging evolutionary rescue potential.

    Args:
        result: CoupledSimResult with yearly_va.
        disease_year: Year disease was introduced.
        save_path: Path to save figure.

    Returns:
        matplotlib Figure.
    """
    years = np.arange(result.n_years)
    va = result.yearly_va

    fig, ax = dark_figure()

    ax.plot(years, va, color='#48c9b0', linewidth=2.5)
    ax.fill_between(years, va, alpha=0.15, color='#48c9b0')

    ax.axvline(disease_year, color='#e94560', linestyle=':', linewidth=1.5,
               alpha=0.8, label=f'Disease intro (yr {disease_year})')

    # Annotate initial and final V_A
    ax.annotate(f'V_A₀ = {va[0]:.4f}', xy=(0, va[0]),
                xytext=(2, va[0] + 0.001), fontsize=9, color=TEXT_COLOR,
                arrowprops=dict(arrowstyle='->', color=TEXT_COLOR, lw=0.8))
    ax.annotate(f'V_A_final = {va[-1]:.4f}', xy=(years[-1], va[-1]),
                xytext=(years[-1] - 4, va[-1] + 0.001),
                fontsize=9, color=TEXT_COLOR,
                arrowprops=dict(arrowstyle='->', color=TEXT_COLOR, lw=0.8))

    ax.set_xlabel('Year', fontsize=12)
    ax.set_ylabel('Additive genetic variance ($V_A$)', fontsize=12)
    ax.set_title('Additive Genetic Variance Over Time', fontsize=14,
                 fontweight='bold')
    ax.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
              labelcolor=TEXT_COLOR, fontsize=10)
    ax.set_xlim(0, result.n_years - 1)
    ax.set_ylim(bottom=0)

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 5. EF1A DYNAMICS
# ═══════════════════════════════════════════════════════════════════════

def plot_ef1a_dynamics(
    result: 'CoupledSimResult',
    disease_year: int = 3,
    expected_eq: float = 0.24,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """EF1A overdominant locus: allele frequency + heterozygosity.

    Shows expected equilibrium frequency.  Separate from additive loci
    because EF1A is under balancing selection (heterozygote advantage).

    Args:
        result: CoupledSimResult with yearly_ef1a_freq.
        disease_year: Year disease was introduced.
        expected_eq: Expected equilibrium allele frequency.
        save_path: Path to save figure.

    Returns:
        matplotlib Figure.
    """
    years = np.arange(result.n_years)
    ef1a_q = result.yearly_ef1a_freq

    fig, (ax1, ax2) = dark_figure(nrows=2, ncols=1, figsize=(10, 8))

    # Top panel: allele frequency
    ax1.plot(years, ef1a_q, color='#f39c12', linewidth=2.5,
             label='EF1A $q$')
    ax1.axhline(expected_eq, color='#2ecc71', linestyle='--',
                linewidth=1.5, alpha=0.7,
                label=f'Expected eq. ($q$ = {expected_eq})')
    ax1.axvline(disease_year, color='#e94560', linestyle=':', linewidth=1.5,
                alpha=0.8)
    ax1.set_ylabel('Allele frequency ($q$)', fontsize=12)
    ax1.set_title('EF1A Overdominant Locus Dynamics', fontsize=14,
                  fontweight='bold')
    ax1.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
               labelcolor=TEXT_COLOR, fontsize=10)
    ax1.set_xlim(0, result.n_years - 1)

    # Bottom panel: expected heterozygosity 2pq
    het = 2 * ef1a_q * (1 - ef1a_q)
    het_eq = 2 * expected_eq * (1 - expected_eq)
    ax2.plot(years, het, color='#3498db', linewidth=2.5,
             label='$H_e$ = 2pq')
    ax2.axhline(het_eq, color='#2ecc71', linestyle='--', linewidth=1.5,
                alpha=0.7, label=f'Expected $H_e$ = {het_eq:.3f}')
    ax2.axvline(disease_year, color='#e94560', linestyle=':', linewidth=1.5,
                alpha=0.8)
    ax2.set_xlabel('Year', fontsize=12)
    ax2.set_ylabel('Expected heterozygosity ($H_e$)', fontsize=12)
    ax2.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
               labelcolor=TEXT_COLOR, fontsize=10)
    ax2.set_xlim(0, result.n_years - 1)

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 6. SELECTION DIFFERENTIAL
# ═══════════════════════════════════════════════════════════════════════

def plot_selection_differential(
    yearly_mean_resistance: np.ndarray,
    yearly_pop: np.ndarray,
    disease_year: int = 3,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Selection differential: Δr per year (difference in mean r).

    Positive = directional selection for resistance.  The differential
    is computed as the year-over-year change in mean resistance, which
    captures the combined effect of differential survival and
    differential reproduction.

    Args:
        yearly_mean_resistance: (n_years,) mean resistance per year.
        yearly_pop: (n_years,) population per year.
        disease_year: Year disease was introduced.
        save_path: Path to save figure.

    Returns:
        matplotlib Figure.
    """
    n_years = len(yearly_mean_resistance)
    years = np.arange(1, n_years)
    delta_r = np.diff(yearly_mean_resistance)

    fig, (ax1, ax2) = dark_figure(nrows=2, ncols=1, figsize=(10, 8))

    # Top: bar chart of selection differential
    colors = [RESISTANCE_HIGH if d > 0 else RESISTANCE_LOW for d in delta_r]
    ax1.bar(years, delta_r, color=colors, alpha=0.8, edgecolor='none',
            width=0.8)
    ax1.axhline(0, color=TEXT_COLOR, linewidth=0.8, alpha=0.5)
    ax1.axvline(disease_year, color='#e94560', linestyle=':', linewidth=1.5,
                alpha=0.8, label=f'Disease intro')

    ax1.set_ylabel('$\\Delta \\bar{r}$ (selection differential)', fontsize=12)
    ax1.set_title('Annual Selection Differential', fontsize=14,
                  fontweight='bold')
    ax1.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
               labelcolor=TEXT_COLOR, fontsize=10)

    # Bottom: cumulative selection response
    cumulative = np.cumsum(delta_r)
    ax2.plot(years, cumulative, color='#48c9b0', linewidth=2.5)
    ax2.fill_between(years, cumulative, alpha=0.15, color='#48c9b0')
    ax2.axhline(0, color=TEXT_COLOR, linewidth=0.8, alpha=0.5)
    ax2.axvline(disease_year, color='#e94560', linestyle=':', linewidth=1.5,
                alpha=0.8)
    ax2.set_xlabel('Year', fontsize=12)
    ax2.set_ylabel('Cumulative $\\Delta \\bar{r}$', fontsize=12)
    ax2.set_title('Cumulative Selection Response', fontsize=14,
                  fontweight='bold')

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 7. HERITABILITY OVER TIME
# ═══════════════════════════════════════════════════════════════════════

def plot_heritability_over_time(
    yearly_va: np.ndarray,
    yearly_vp: np.ndarray,
    disease_year: int = 3,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Narrow-sense h² = V_A / V_P over time.

    Should start moderate, potentially decline as genetic variance
    is consumed by selection.

    Args:
        yearly_va: (n_years,) additive genetic variance per year.
        yearly_vp: (n_years,) phenotypic variance per year.
        disease_year: Year disease was introduced.
        save_path: Path to save figure.

    Returns:
        matplotlib Figure.
    """
    n_years = len(yearly_va)
    years = np.arange(n_years)

    # h² = V_A / V_P (guard against division by zero)
    h2 = np.where(yearly_vp > 1e-12, yearly_va / yearly_vp, 0.0)

    fig, (ax1, ax2) = dark_figure(nrows=2, ncols=1, figsize=(10, 8))

    # Top: h² over time
    ax1.plot(years, h2, color='#f39c12', linewidth=2.5, label='$h^2$')
    ax1.axhline(h2[0], color='#2ecc71', linestyle='--', linewidth=1,
                alpha=0.6, label=f'Initial $h^2$ = {h2[0]:.3f}')
    ax1.axvline(disease_year, color='#e94560', linestyle=':', linewidth=1.5,
                alpha=0.8)
    ax1.set_ylabel('$h^2 = V_A / V_P$', fontsize=12)
    ax1.set_title('Narrow-Sense Heritability Over Time',
                  fontsize=14, fontweight='bold')
    ax1.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
               labelcolor=TEXT_COLOR, fontsize=10)
    ax1.set_xlim(0, n_years - 1)
    ax1.set_ylim(0, 1)

    # Bottom: V_A and V_P decomposition
    ax2.plot(years, yearly_va, color='#48c9b0', linewidth=2,
             label='$V_A$ (additive)')
    ax2.plot(years, yearly_vp, color='#e94560', linewidth=2,
             label='$V_P$ (phenotypic)')
    ve = yearly_vp - yearly_va
    ax2.fill_between(years, yearly_va, yearly_vp, alpha=0.15,
                     color='#e94560', label='$V_E$ (environmental)')
    ax2.fill_between(years, 0, yearly_va, alpha=0.15,
                     color='#48c9b0')
    ax2.axvline(disease_year, color='#e94560', linestyle=':', linewidth=1.5,
                alpha=0.8)
    ax2.set_xlabel('Year', fontsize=12)
    ax2.set_ylabel('Variance', fontsize=12)
    ax2.set_title('Variance Decomposition', fontsize=14, fontweight='bold')
    ax2.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
               labelcolor=TEXT_COLOR, fontsize=9)
    ax2.set_xlim(0, n_years - 1)
    ax2.set_ylim(bottom=0)

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 8. GENOTYPE → PHENOTYPE MAP
# ═══════════════════════════════════════════════════════════════════════

def plot_genotype_phenotype_map(
    genotypes: np.ndarray,
    agents: np.ndarray,
    effect_sizes: np.ndarray,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Scatter: genotype score (sum of allele counts × effects) vs r_i.

    Shows the linear mapping from genotype to resistance phenotype.

    Args:
        genotypes: (max_agents, N_LOCI, 2) int8.
        agents: Agent structured array.
        effect_sizes: (N_ADDITIVE,) float64 effect sizes.
        save_path: Path to save figure.

    Returns:
        matplotlib Figure.
    """
    from sswd_evoepi.types import N_ADDITIVE, IDX_EF1A

    alive_mask = agents['alive']
    alive_idx = np.where(alive_mask)[0]
    if len(alive_idx) == 0:
        fig, ax = dark_figure()
        ax.text(0.5, 0.5, 'No alive agents', transform=ax.transAxes,
                ha='center', va='center', color=TEXT_COLOR, fontsize=14)
        if save_path:
            save_figure(fig, save_path)
        return fig

    # Additive genotype score
    alive_geno = genotypes[alive_idx]
    allele_sums = alive_geno[:, :N_ADDITIVE, :].sum(axis=2)  # 0, 1, or 2
    additive_score = (allele_sums.astype(np.float64) * 0.5) @ effect_sizes

    # EF1A bonus
    ef1a_sum = alive_geno[:, IDX_EF1A, :].sum(axis=1)
    has_od_bonus = ef1a_sum == 1  # heterozygotes

    r_i = agents['resistance'][alive_idx]

    fig, ax = dark_figure()

    # Plot non-heterozygotes and heterozygotes with different markers
    non_het = ~has_od_bonus
    ax.scatter(additive_score[non_het], r_i[non_het],
               c=[RESISTANCE_LOW], alpha=0.4, s=12, edgecolors='none',
               label='EF1A hom (no OD bonus)')
    ax.scatter(additive_score[has_od_bonus], r_i[has_od_bonus],
               c=[RESISTANCE_HIGH], alpha=0.5, s=15, edgecolors='none',
               marker='D', label='EF1A het (+OD bonus)')

    # Perfect-mapping reference line
    x_range = np.linspace(additive_score.min(), additive_score.max(), 100)
    ax.plot(x_range, x_range, color=TEXT_COLOR, linestyle='--',
            linewidth=1, alpha=0.5, label='1:1 additive')
    ax.plot(x_range, x_range + 0.160, color='#f39c12', linestyle='--',
            linewidth=1, alpha=0.5, label='1:1 + W_OD')

    ax.set_xlabel('Additive genotype score', fontsize=12)
    ax.set_ylabel('Resistance ($r_i$)', fontsize=12)
    ax.set_title('Genotype → Phenotype Map', fontsize=14, fontweight='bold')
    ax.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
              labelcolor=TEXT_COLOR, fontsize=9)

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 9. LOCUS EFFECT SIZE DISTRIBUTION
# ═══════════════════════════════════════════════════════════════════════

def plot_locus_effect_size_distribution(
    effect_sizes: np.ndarray,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Bar chart of effect sizes for all 51 additive loci, ranked.

    Shows the exponential decay characteristic of polygenic architecture.

    Args:
        effect_sizes: (N_ADDITIVE,) float64 (assumed sorted descending).
        save_path: Path to save figure.

    Returns:
        matplotlib Figure.
    """
    n = len(effect_sizes)
    fig, ax = dark_figure()

    # Sort descending for display
    sorted_e = np.sort(effect_sizes)[::-1]
    locus_idx = np.arange(n)

    colors = [_resistance_color(e, vmin=sorted_e.min(), vmax=sorted_e.max())
              for e in sorted_e]

    ax.bar(locus_idx, sorted_e, color=colors, edgecolor='none',
           width=0.8, alpha=0.9)

    # Overlay fitted exponential
    from scipy.optimize import curve_fit
    def _exp_func(x, a, b):
        return a * np.exp(-b * x)
    try:
        popt, _ = curve_fit(_exp_func, locus_idx, sorted_e,
                            p0=[sorted_e[0], 0.05], maxfev=5000)
        fit_x = np.linspace(0, n - 1, 200)
        ax.plot(fit_x, _exp_func(fit_x, *popt), color=TEXT_COLOR,
                linewidth=1.5, linestyle='--', alpha=0.7,
                label=f'Exp fit: {popt[0]:.3f}·e^(-{popt[1]:.3f}x)')
    except Exception:
        pass

    ax.set_xlabel('Locus rank (by effect size)', fontsize=12)
    ax.set_ylabel('Effect size ($\\tilde{e}_l$)', fontsize=12)
    ax.set_title('Additive Locus Effect Size Distribution',
                 fontsize=14, fontweight='bold')
    ax.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
              labelcolor=TEXT_COLOR, fontsize=10)
    ax.set_xlim(-0.5, n - 0.5)

    # Summary annotation
    ax.text(0.95, 0.85,
            f'Σ = {sorted_e.sum():.3f}\n'
            f'Top 5 = {sorted_e[:5].sum():.3f} '
            f'({sorted_e[:5].sum() / sorted_e.sum() * 100:.0f}%)\n'
            f'n = {n}',
            transform=ax.transAxes, fontsize=9, color=TEXT_COLOR,
            ha='right', va='top',
            bbox=dict(boxstyle='round,pad=0.3', facecolor=DARK_PANEL,
                      edgecolor=GRID_COLOR, alpha=0.8))

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 10. RESISTANCE BY NODE — VIOLIN PLOT
# ═══════════════════════════════════════════════════════════════════════

def plot_resistance_by_node_violin(
    node_resistance: Dict[str, np.ndarray],
    year: int = 0,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Violin plot of resistance distributions at each node at a given year.

    Shows node-to-node variation in evolutionary response.

    Args:
        node_resistance: Mapping {node_name: 1-D resistance array}.
        year: Year label for the title.
        save_path: Path to save figure.

    Returns:
        matplotlib Figure.
    """
    fig, ax = dark_figure(figsize=(12, 6))

    names = list(node_resistance.keys())
    data = [node_resistance[n] for n in names]
    positions = np.arange(len(names))

    parts = ax.violinplot(data, positions=positions, showmeans=True,
                          showextrema=True, showmedians=True)

    # Colour each violin by its mean resistance
    for idx, body in enumerate(parts['bodies']):
        mean_r = np.mean(data[idx])
        body.set_facecolor(_resistance_color(mean_r))
        body.set_edgecolor(TEXT_COLOR)
        body.set_alpha(0.7)

    for key in ('cmeans', 'cmedians', 'cmins', 'cmaxes', 'cbars'):
        if key in parts:
            parts[key].set_color(TEXT_COLOR)
            parts[key].set_linewidth(1.0)

    # Mean annotations
    for idx, name in enumerate(names):
        mean_r = np.mean(data[idx])
        ax.text(idx, mean_r + 0.01, f'{mean_r:.3f}',
                ha='center', fontsize=8, color=TEXT_COLOR)

    ax.set_xticks(positions)
    ax.set_xticklabels(names, rotation=30, ha='right', fontsize=10)
    ax.set_ylabel('Resistance ($r_i$)', fontsize=12)
    ax.set_title(f'Resistance Distribution by Node — Year {year}',
                 fontsize=14, fontweight='bold')

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 11. GENETIC DRIFT NULL COMPARISON
# ═══════════════════════════════════════════════════════════════════════

def plot_genetic_drift_null(
    disease_results: List['CoupledSimResult'],
    null_results: List['CoupledSimResult'],
    disease_year: int = 3,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Multi-seed drift null vs. disease comparison.

    No-disease runs show neutral drift envelope; disease runs show
    directional selection response.

    Args:
        disease_results: List of CoupledSimResult with disease.
        null_results: List of CoupledSimResult without disease (drift only).
        disease_year: Year disease was introduced.
        save_path: Path to save figure.

    Returns:
        matplotlib Figure.
    """
    fig, ax = dark_figure()

    n_years = max(
        max((r.n_years for r in disease_results), default=0),
        max((r.n_years for r in null_results), default=0),
    )
    years = np.arange(n_years)

    # Null (drift) runs — gray
    null_trajs = []
    for r in null_results:
        yrs = np.arange(r.n_years)
        ax.plot(yrs, r.yearly_mean_resistance,
                color='#7f8c8d', alpha=0.3, linewidth=0.8)
        null_trajs.append(r.yearly_mean_resistance)

    # Null envelope
    if null_trajs:
        min_len = min(len(t) for t in null_trajs)
        null_arr = np.array([t[:min_len] for t in null_trajs])
        null_mean = null_arr.mean(axis=0)
        null_std = null_arr.std(axis=0)
        yrs_null = np.arange(min_len)
        ax.fill_between(yrs_null, null_mean - 2 * null_std,
                        null_mean + 2 * null_std,
                        color='#7f8c8d', alpha=0.12, label='Drift ±2σ envelope')
        ax.plot(yrs_null, null_mean, color='#7f8c8d',
                linewidth=2, linestyle='--', label='Drift mean')

    # Disease runs — coloured
    disease_trajs = []
    for r in disease_results:
        yrs = np.arange(r.n_years)
        ax.plot(yrs, r.yearly_mean_resistance,
                color='#e94560', alpha=0.35, linewidth=0.8)
        disease_trajs.append(r.yearly_mean_resistance)

    if disease_trajs:
        min_len = min(len(t) for t in disease_trajs)
        dis_arr = np.array([t[:min_len] for t in disease_trajs])
        dis_mean = dis_arr.mean(axis=0)
        yrs_dis = np.arange(min_len)
        ax.plot(yrs_dis, dis_mean, color='#e94560', linewidth=2.5,
                label='Disease mean')

    ax.axvline(disease_year, color='#f39c12', linestyle=':', linewidth=1.5,
               alpha=0.8, label=f'Disease intro')

    ax.set_xlabel('Year', fontsize=12)
    ax.set_ylabel('Mean resistance ($\\bar{r}_i$)', fontsize=12)
    ax.set_title('Selection vs. Neutral Drift',
                 fontsize=14, fontweight='bold')
    ax.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
              labelcolor=TEXT_COLOR, fontsize=9)
    ax.set_xlim(0, n_years - 1)

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 12. BETA INITIALIZATION VISUALIZATION
# ═══════════════════════════════════════════════════════════════════════

def plot_beta_init_visualization(
    beta_a: float = 2.0,
    beta_b: float = 8.0,
    target_mean_r: float = 0.15,
    n_agents: int = 500,
    seed: int = 42,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Three-panel figure showing Beta(a,b) initialization.

    Panel 1: Beta(a,b) PDF
    Panel 2: Resulting per-locus q values (sorted by effect size)
    Panel 3: Population resistance distribution at t=0

    Args:
        beta_a: Beta shape parameter a.
        beta_b: Beta shape parameter b.
        target_mean_r: Target population mean resistance.
        n_agents: Population size for demonstration.
        seed: RNG seed.
        save_path: Path to save figure.

    Returns:
        matplotlib Figure.
    """
    from sswd_evoepi.genetics import (
        initialize_effect_sizes,
        initialize_genotypes,
        compute_resistance_batch,
        W_OD,
    )
    from sswd_evoepi.types import N_ADDITIVE, allocate_agents

    rng = np.random.default_rng(seed)
    effects = initialize_effect_sizes(rng, n_additive=N_ADDITIVE)

    fig, (ax1, ax2, ax3) = dark_figure(nrows=1, ncols=3, figsize=(18, 5))

    # Panel 1: Beta PDF
    x = np.linspace(0, 1, 500)
    pdf = sp_stats.beta.pdf(x, beta_a, beta_b)
    ax1.fill_between(x, pdf, alpha=0.3, color='#48c9b0')
    ax1.plot(x, pdf, color='#48c9b0', linewidth=2.5)
    ax1.set_xlabel('Allele frequency ($q$)', fontsize=11)
    ax1.set_ylabel('Density', fontsize=11)
    ax1.set_title(f'Beta({beta_a}, {beta_b}) Prior', fontsize=13,
                  fontweight='bold')
    # Stats
    mean_q = beta_a / (beta_a + beta_b)
    ax1.axvline(mean_q, color='#f39c12', linestyle='--', linewidth=1.5,
                label=f'E[q] = {mean_q:.3f}')
    ax1.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
               labelcolor=TEXT_COLOR, fontsize=9)

    # Panel 2: Per-locus q values
    # Draw raw qs and scale
    raw_q = rng.beta(beta_a, beta_b, size=N_ADDITIVE)
    ef1a_het_contrib = 2 * 0.24 * (1 - 0.24) * W_OD
    target_additive = max(0.0, target_mean_r - ef1a_het_contrib)
    current_additive = np.dot(effects, raw_q)
    if current_additive > 0:
        scale = target_additive / current_additive
        q_vals = np.clip(raw_q * scale, 0.001, 0.999)
    else:
        q_vals = np.full(N_ADDITIVE, 0.01)

    locus_idx = np.arange(N_ADDITIVE)
    colors = [_resistance_color(e, vmin=effects.min(), vmax=effects.max())
              for e in effects]
    ax2.bar(locus_idx, q_vals, color=colors, edgecolor='none',
            width=0.8, alpha=0.9)
    ax2.axhline(np.mean(q_vals), color='#f39c12', linestyle='--',
                linewidth=1.5, label=f'Mean q = {np.mean(q_vals):.3f}')
    ax2.set_xlabel('Locus (ranked by effect)', fontsize=11)
    ax2.set_ylabel('Allele frequency ($q_l$)', fontsize=11)
    ax2.set_title('Scaled Per-Locus Frequencies', fontsize=13,
                  fontweight='bold')
    ax2.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
               labelcolor=TEXT_COLOR, fontsize=9)

    # Panel 3: Resulting population r_i distribution
    geno = initialize_genotypes(
        n_agents=n_agents, effects=effects, rng=rng,
        target_mean_r=target_mean_r, q_additive=q_vals,
    )
    # Compute resistance
    agents = allocate_agents(n_agents)
    agents['alive'][:n_agents] = True
    alive_mask = agents['alive']
    r_vals = compute_resistance_batch(geno, effects, alive_mask[:n_agents])
    r_vals = r_vals[:n_agents]

    ax3.hist(r_vals, bins=40, color='#48c9b0', alpha=0.7,
             edgecolor=DARK_BG, linewidth=0.5)
    ax3.axvline(np.mean(r_vals), color='#e94560', linestyle='--',
                linewidth=2, label=f'Mean $r_i$ = {np.mean(r_vals):.3f}')
    ax3.axvline(target_mean_r, color='#f39c12', linestyle=':',
                linewidth=1.5, label=f'Target = {target_mean_r}')
    ax3.set_xlabel('Resistance ($r_i$)', fontsize=11)
    ax3.set_ylabel('Count', fontsize=11)
    ax3.set_title('Population $r_i$ at $t=0$', fontsize=13,
                  fontweight='bold')
    ax3.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
               labelcolor=TEXT_COLOR, fontsize=9)

    fig.suptitle(
        f'Genetic Initialization: Beta({beta_a},{beta_b}) → '
        f'target $\\bar{{r}}$ = {target_mean_r}',
        fontsize=15, fontweight='bold', color=TEXT_COLOR, y=1.02,
    )

    if save_path:
        save_figure(fig, save_path)
    return fig
