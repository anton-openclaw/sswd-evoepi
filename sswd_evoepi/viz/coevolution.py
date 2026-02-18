"""Pathogen co-evolution visualizations for SSWD-EvoEpi.

Every function:
  - Accepts model results (CoupledSimResult or SpatialSimResult) or config as input
  - Returns a matplotlib Figure
  - Has an optional ``save_path`` parameter (saves PNG when given)
  - Uses the shared dark theme from ``sswd_evoepi.viz.style``

matplotlib backend is forced to Agg (no display) on import.

Virulence color gradient:
  low=#2ecc71 (green/mild) → high=#e74c3c (red/virulent)

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
    from sswd_evoepi.config import (
        DiseaseSection,
        PathogenEvolutionSection,
        SimulationConfig,
    )


# ═══════════════════════════════════════════════════════════════════════
# VIRULENCE COLOUR MAP
# ═══════════════════════════════════════════════════════════════════════

VIRULENCE_LOW = '#2ecc71'    # green (mild)
VIRULENCE_HIGH = '#e74c3c'   # red (virulent)

_VIRULENCE_CMAP = mcolors.LinearSegmentedColormap.from_list(
    'virulence', [VIRULENCE_LOW, '#f39c12', VIRULENCE_HIGH], N=256,
)


def _virulence_color(value: float, vmin: float = 0.0,
                     vmax: float = 1.0) -> tuple:
    """Map a virulence value to an RGBA tuple."""
    norm = max(0.0, min(1.0, (value - vmin) / (vmax - vmin)))
    return _VIRULENCE_CMAP(norm)


# ═══════════════════════════════════════════════════════════════════════
# HELPER: compute virulence SD from yearly_virulence_new_infections
# ═══════════════════════════════════════════════════════════════════════

def _estimate_virulence_sd(result: 'CoupledSimResult') -> np.ndarray:
    """Estimate virulence SD from spread between new-infection and death virulences.

    When only the mean is available (single-node results only store means),
    we approximate SD from the difference between virulence of new infections
    vs virulence of deaths (higher-v strains die faster → selection gradient).
    Falls back to a fixed 0.05 if data is missing.
    """
    n = result.n_years
    sd = np.full(n, 0.05)
    if (result.yearly_virulence_new_infections is not None
            and result.yearly_virulence_of_deaths is not None):
        diff = np.abs(
            result.yearly_virulence_new_infections
            - result.yearly_virulence_of_deaths
        )
        # Scale: the spread between new infections and deaths ≈ 2 SD
        sd = np.maximum(diff / 2.0, 0.01)
    return sd


# ═══════════════════════════════════════════════════════════════════════
# 1. VIRULENCE TRAJECTORY
# ═══════════════════════════════════════════════════════════════════════

def plot_virulence_trajectory(
    result: 'CoupledSimResult',
    disease_year: int = 3,
    v_init: Optional[float] = None,
    v_ess: Optional[float] = None,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Mean pathogen virulence over time with ±SD band.

    THE key co-evolution plot.  Shows whether virulence evolves up or down
    following epidemic introduction.

    Args:
        result: CoupledSimResult with yearly_mean_virulence.
        disease_year: Year disease was introduced (vertical marker).
        v_init: Initial virulence (horizontal reference line).
        v_ess: Predicted ESS virulence (horizontal reference line).
        save_path: Optional path to save the figure.

    Returns:
        matplotlib Figure.
    """
    years = np.arange(result.n_years)
    mean_v = result.yearly_mean_virulence
    sd_v = _estimate_virulence_sd(result)

    fig, ax = dark_figure()

    # ±SD band
    valid = mean_v > 0
    ax.fill_between(
        years[valid],
        (mean_v - sd_v)[valid],
        (mean_v + sd_v)[valid],
        alpha=0.25, color=ACCENT_COLORS[6],  # red
        label='±1 SD',
    )

    # Mean line — color by virulence
    for i in range(len(years) - 1):
        if valid[i] and valid[i + 1]:
            ax.plot(
                years[i:i + 2], mean_v[i:i + 2],
                color=_virulence_color(mean_v[i]),
                linewidth=2.5, solid_capstyle='round',
            )

    # Disease introduction
    ax.axvline(disease_year, color=ACCENT_COLORS[4], linestyle=':',
               linewidth=1.5, alpha=0.8, label=f'Disease intro (year {disease_year})')

    # Reference lines
    if v_init is not None:
        ax.axhline(v_init, color=ACCENT_COLORS[3], linestyle='--',
                   linewidth=1.2, alpha=0.7, label=f'v_init = {v_init:.2f}')
    if v_ess is not None:
        ax.axhline(v_ess, color=ACCENT_COLORS[7], linestyle='--',
                   linewidth=1.2, alpha=0.7, label=f'Predicted ESS = {v_ess:.2f}')

    ax.set_xlabel('Year', fontsize=12)
    ax.set_ylabel('Mean virulence (v)', fontsize=12)
    ax.set_title('Pathogen Virulence Trajectory', fontsize=14, fontweight='bold')
    ax.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
              labelcolor=TEXT_COLOR, fontsize=10)
    ax.set_xlim(0, result.n_years - 1)
    ax.set_ylim(0, 1)

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 2. COEVOLUTION PHASE PORTRAIT
# ═══════════════════════════════════════════════════════════════════════

def plot_coevolution_phase_portrait(
    result: 'CoupledSimResult',
    disease_year: int = 3,
    cmap_name: str = 'plasma',
    save_path: Optional[str] = None,
) -> plt.Figure:
    """2D trajectory: mean host resistance vs mean pathogen virulence.

    Time encoded as a colour gradient along the path.  Arrows indicate
    direction of co-evolutionary dynamics (arms race).

    Args:
        result: CoupledSimResult with yearly_mean_resistance and
            yearly_mean_virulence.
        disease_year: Year disease was introduced (start marker).
        cmap_name: Matplotlib colourmap for time axis.
        save_path: Optional path to save the figure.

    Returns:
        matplotlib Figure.
    """
    res = result.yearly_mean_resistance
    vir = result.yearly_mean_virulence

    # Only plot years where both metrics are valid
    valid = (res > 0) & (vir > 0)
    idx = np.where(valid)[0]
    if len(idx) < 2:
        # Not enough data; return empty figure
        fig, ax = dark_figure()
        ax.text(0.5, 0.5, 'Insufficient co-evolution data',
                ha='center', va='center', color=TEXT_COLOR,
                fontsize=14, transform=ax.transAxes)
        if save_path:
            save_figure(fig, save_path)
        return fig

    r_vals = res[idx]
    v_vals = vir[idx]
    t_norm = np.linspace(0, 1, len(idx))

    cmap = plt.get_cmap(cmap_name)

    fig, ax = dark_figure()

    # Draw coloured line segments
    for i in range(len(idx) - 1):
        ax.plot(
            r_vals[i:i + 2], v_vals[i:i + 2],
            color=cmap(t_norm[i]),
            linewidth=2.5, solid_capstyle='round',
        )

    # Arrows every ~15% of trajectory to show direction
    n_arrows = max(3, len(idx) // 7)
    arrow_indices = np.linspace(0, len(idx) - 2, n_arrows, dtype=int)
    for ai in arrow_indices:
        dx = r_vals[ai + 1] - r_vals[ai]
        dy = v_vals[ai + 1] - v_vals[ai]
        ax.annotate(
            '', xy=(r_vals[ai + 1], v_vals[ai + 1]),
            xytext=(r_vals[ai], v_vals[ai]),
            arrowprops=dict(
                arrowstyle='->', color=cmap(t_norm[ai]),
                lw=1.8, mutation_scale=15,
            ),
        )

    # Start and end markers
    ax.scatter(r_vals[0], v_vals[0], color=cmap(0.0), s=120,
               edgecolors='white', linewidths=1.5, zorder=5,
               label=f'Year {idx[0]}')
    ax.scatter(r_vals[-1], v_vals[-1], color=cmap(1.0), s=120,
               edgecolors='white', linewidths=1.5, zorder=5, marker='*',
               label=f'Year {idx[-1]}')

    # Colourbar for time
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(idx[0], idx[-1]))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label('Year', color=TEXT_COLOR, fontsize=11)
    cbar.ax.yaxis.set_tick_params(color=TEXT_COLOR)
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color=TEXT_COLOR)

    ax.set_xlabel('Mean host resistance (r̄)', fontsize=12)
    ax.set_ylabel('Mean pathogen virulence (v̄)', fontsize=12)
    ax.set_title('Host–Pathogen Co-Evolution Phase Portrait',
                 fontsize=14, fontweight='bold')
    ax.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
              labelcolor=TEXT_COLOR, fontsize=10, loc='upper left')

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 3. VIRULENCE DISTRIBUTION OVER TIME
# ═══════════════════════════════════════════════════════════════════════

def plot_virulence_distribution_over_time(
    result: 'CoupledSimResult',
    disease_year: int = 3,
    n_timepoints: int = 5,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Violin plots of virulence distribution at multiple timepoints.

    Since individual-level virulence distributions aren't stored in
    CoupledSimResult (only means), we reconstruct approximate distributions
    using the mean and estimated SD, showing the spread of pathogen strains.

    For SpatialSimResult, uses per-node virulence as samples.

    Args:
        result: CoupledSimResult or SpatialSimResult with virulence data.
        disease_year: Year disease was introduced.
        n_timepoints: Number of timepoints to show.
        save_path: Optional path to save the figure.

    Returns:
        matplotlib Figure.
    """
    mean_v = result.yearly_mean_virulence
    sd_v = _estimate_virulence_sd(result)

    # Pick timepoints: from disease_year onward
    valid_years = np.where(mean_v > 0)[0]
    if len(valid_years) < 2:
        fig, ax = dark_figure()
        ax.text(0.5, 0.5, 'Insufficient virulence data',
                ha='center', va='center', color=TEXT_COLOR,
                fontsize=14, transform=ax.transAxes)
        if save_path:
            save_figure(fig, save_path)
        return fig

    # Sample evenly across valid years
    tp_idx = np.linspace(0, len(valid_years) - 1, min(n_timepoints, len(valid_years)),
                         dtype=int)
    timepoints = valid_years[tp_idx]

    fig, ax = dark_figure()

    # Generate synthetic distributions (normal, clipped to [0,1])
    violin_data = []
    positions = []
    for i, yr in enumerate(timepoints):
        samples = np.random.default_rng(yr).normal(
            mean_v[yr], max(sd_v[yr], 0.01), 200,
        )
        samples = np.clip(samples, 0, 1)
        violin_data.append(samples)
        positions.append(i)

    parts = ax.violinplot(
        violin_data, positions=positions, showmeans=True,
        showextrema=False, widths=0.7,
    )

    # Colour violins by mean virulence
    for i, body in enumerate(parts['bodies']):
        colour = _virulence_color(mean_v[timepoints[i]])
        body.set_facecolor(colour)
        body.set_edgecolor('white')
        body.set_alpha(0.7)
    parts['cmeans'].set_color('white')
    parts['cmeans'].set_linewidth(1.5)

    ax.set_xticks(positions)
    ax.set_xticklabels([f'Year {yr}' for yr in timepoints], fontsize=10)
    ax.set_ylabel('Pathogen virulence (v)', fontsize=12)
    ax.set_title('Virulence Distribution Over Time',
                 fontsize=14, fontweight='bold')
    ax.set_ylim(0, 1)

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 4. TRADEOFF CURVE
# ═══════════════════════════════════════════════════════════════════════

def plot_tradeoff_curve(
    pe_cfg: 'PathogenEvolutionSection',
    disease_cfg: 'DiseaseSection',
    T_celsius: float = 15.0,
    v_current: Optional[float] = None,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Static plot of the virulence-transmission tradeoff curve.

    x-axis: virulence (v)
    y1: shedding rate (σ₂, late infectious — primary transmission)
    y2: I₂→D death rate (µ_I2D — cost of virulence)
    y3: total lifetime output (TLO = σ₂ / µ_I2D — proxy for fitness)

    Shows how the mechanistic tradeoff constrains pathogen evolution.

    Args:
        pe_cfg: PathogenEvolutionSection config.
        disease_cfg: DiseaseSection config.
        T_celsius: Temperature for Arrhenius rates.
        v_current: Optional current mean virulence to mark.
        save_path: Optional path to save the figure.

    Returns:
        matplotlib Figure.
    """
    from sswd_evoepi.disease import (
        sigma_2_strain, mu_I2D_strain, sigma_1_strain, arrhenius,
    )

    v_range = np.linspace(pe_cfg.v_min + 0.01, pe_cfg.v_max - 0.01, 200)

    shed = np.array([sigma_2_strain(v, T_celsius, disease_cfg, pe_cfg)
                     for v in v_range])
    death = np.array([mu_I2D_strain(v, T_celsius, disease_cfg, pe_cfg)
                      for v in v_range])
    tlo = shed / death  # total lifetime output

    fig, axes = dark_figure(nrows=3, ncols=1, figsize=(10, 12))

    # Panel 1: Shedding rate
    ax1 = axes[0]
    ax1.plot(v_range, shed, color=ACCENT_COLORS[5], linewidth=2.5,
             label='σ₂(v) — late shedding rate')
    # Also show early shedding for comparison
    shed_early = np.array([sigma_1_strain(v, T_celsius, disease_cfg, pe_cfg)
                           for v in v_range])
    ax1.plot(v_range, shed_early, color=ACCENT_COLORS[3], linewidth=2,
             linestyle='--', alpha=0.7, label='σ₁(v) — early shedding rate')
    ax1.set_ylabel('Shedding rate (bact/mL/d)', fontsize=11)
    ax1.set_title('Virulence–Transmission Tradeoff', fontsize=14,
                  fontweight='bold')
    ax1.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
               labelcolor=TEXT_COLOR, fontsize=9)

    # Panel 2: Death rate
    ax2 = axes[1]
    ax2.plot(v_range, death, color=ACCENT_COLORS[6], linewidth=2.5,
             label='µ_I2D(v) — disease death rate')
    ax2.set_ylabel('I₂→D death rate (d⁻¹)', fontsize=11)
    ax2.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
               labelcolor=TEXT_COLOR, fontsize=9)

    # Panel 3: TLO (fitness proxy)
    ax3 = axes[2]
    ax3.plot(v_range, tlo, color=ACCENT_COLORS[7], linewidth=2.5,
             label='TLO = σ₂/µ_I2D')
    ax3.set_ylabel('Total lifetime output (bact/mL)', fontsize=11)
    ax3.set_xlabel('Virulence (v)', fontsize=12)
    ax3.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
               labelcolor=TEXT_COLOR, fontsize=9)

    # Mark TLO peak
    tlo_peak_idx = np.argmax(tlo)
    v_peak = v_range[tlo_peak_idx]
    ax3.axvline(v_peak, color='white', linestyle=':', linewidth=1,
                alpha=0.6)
    ax3.annotate(
        f'TLO peak\nv = {v_peak:.2f}',
        xy=(v_peak, tlo[tlo_peak_idx]),
        xytext=(v_peak + 0.1, tlo[tlo_peak_idx] * 0.85),
        arrowprops=dict(arrowstyle='->', color='white', lw=1.2),
        color='white', fontsize=10,
    )

    # Mark current virulence on all panels if provided
    if v_current is not None:
        for ax in axes:
            ax.axvline(v_current, color=ACCENT_COLORS[4], linestyle='--',
                       linewidth=1.5, alpha=0.8)
        axes[0].annotate(
            f'v̄ = {v_current:.2f}', xy=(v_current, shed[0]),
            xytext=(v_current + 0.08, shed.max() * 0.7),
            arrowprops=dict(arrowstyle='->', color=ACCENT_COLORS[4], lw=1.2),
            color=ACCENT_COLORS[4], fontsize=10,
        )

    # Mark v_anchor
    for ax in axes:
        ax.axvline(pe_cfg.v_anchor, color=GRID_COLOR, linestyle='-.',
                   linewidth=1, alpha=0.5)

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 5. R₀ BY VIRULENCE
# ═══════════════════════════════════════════════════════════════════════

def plot_R0_by_virulence(
    disease_cfg: 'DiseaseSection',
    pe_cfg: 'PathogenEvolutionSection',
    densities: Optional[List[int]] = None,
    T_celsius: float = 15.0,
    phi_k: float = 0.5,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """R₀ as a function of virulence for different host densities.

    Shows the intermediate virulence optimum and how it shifts with
    host density — a key prediction of virulence evolution theory.

    Args:
        disease_cfg: DiseaseSection config.
        pe_cfg: PathogenEvolutionSection config.
        densities: List of host population sizes to plot.
        T_celsius: Temperature for Arrhenius.
        phi_k: Flushing rate.
        save_path: Optional path to save the figure.

    Returns:
        matplotlib Figure.
    """
    from sswd_evoepi.disease import compute_R0

    if densities is None:
        densities = [100, 300, 500, 800, 1000]

    v_range = np.linspace(pe_cfg.v_min + 0.01, pe_cfg.v_max - 0.01, 150)

    fig, ax = dark_figure()

    colors = [VIRULENCE_LOW, ACCENT_COLORS[3], ACCENT_COLORS[4],
              ACCENT_COLORS[5], VIRULENCE_HIGH]

    for i, S0 in enumerate(densities):
        R0_vals = np.array([
            compute_R0(T_celsius, S0, phi_k, disease_cfg,
                       v=v, pe_cfg=pe_cfg)
            for v in v_range
        ])
        color = colors[i % len(colors)]
        ax.plot(v_range, R0_vals, color=color, linewidth=2,
                label=f'S₀ = {S0}')

        # Mark peak
        peak_idx = np.argmax(R0_vals)
        ax.scatter(v_range[peak_idx], R0_vals[peak_idx],
                   color=color, s=60, zorder=5, edgecolors='white',
                   linewidths=1)

    ax.axhline(1.0, color='white', linestyle=':', linewidth=1,
               alpha=0.5, label='R₀ = 1')

    ax.set_xlabel('Virulence (v)', fontsize=12)
    ax.set_ylabel('R₀', fontsize=12)
    ax.set_title('R₀ vs Virulence by Host Density',
                 fontsize=14, fontweight='bold')
    ax.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
              labelcolor=TEXT_COLOR, fontsize=10)
    ax.set_xlim(v_range[0], v_range[-1])
    ax.set_ylim(bottom=0)

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 6. VIRULENCE VS HOST DENSITY
# ═══════════════════════════════════════════════════════════════════════

def plot_virulence_vs_host_density(
    result: 'CoupledSimResult',
    disease_year: int = 3,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Scatter of mean virulence vs total population over time.

    Shows whether virulence tracks density changes — an expected outcome
    if density-dependent virulence evolution is operating.

    Args:
        result: CoupledSimResult with yearly_mean_virulence and yearly_pop.
        disease_year: Year disease was introduced.
        save_path: Optional path to save the figure.

    Returns:
        matplotlib Figure.
    """
    vir = result.yearly_mean_virulence
    pop = result.yearly_pop

    valid = vir > 0
    if np.sum(valid) < 2:
        fig, ax = dark_figure()
        ax.text(0.5, 0.5, 'Insufficient data',
                ha='center', va='center', color=TEXT_COLOR,
                fontsize=14, transform=ax.transAxes)
        if save_path:
            save_figure(fig, save_path)
        return fig

    years = np.arange(result.n_years)
    t_norm = np.zeros(result.n_years)
    valid_idx = np.where(valid)[0]
    if len(valid_idx) > 1:
        t_norm[valid_idx] = np.linspace(0, 1, len(valid_idx))

    fig, ax = dark_figure()

    cmap = plt.get_cmap('plasma')
    scatter = ax.scatter(
        pop[valid], vir[valid],
        c=t_norm[valid], cmap=cmap, s=80,
        edgecolors='white', linewidths=0.8, zorder=3,
    )

    # Connect with thin line
    ax.plot(pop[valid], vir[valid], color=GRID_COLOR, linewidth=0.8,
            alpha=0.5, zorder=2)

    # Mark start and end
    ax.scatter(pop[valid_idx[0]], vir[valid_idx[0]], color='white',
               s=150, marker='o', edgecolors=ACCENT_COLORS[0],
               linewidths=2, zorder=6, label=f'Year {valid_idx[0]}')
    ax.scatter(pop[valid_idx[-1]], vir[valid_idx[-1]], color='white',
               s=200, marker='*', edgecolors=ACCENT_COLORS[0],
               linewidths=2, zorder=6, label=f'Year {valid_idx[-1]}')

    cbar = fig.colorbar(scatter, ax=ax, pad=0.02)
    cbar.set_label('Time →', color=TEXT_COLOR, fontsize=11)
    cbar.ax.yaxis.set_tick_params(color=TEXT_COLOR)
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color=TEXT_COLOR)

    ax.set_xlabel('Population size', fontsize=12)
    ax.set_ylabel('Mean virulence (v)', fontsize=12)
    ax.set_title('Virulence vs Host Density Over Time',
                 fontsize=14, fontweight='bold')
    ax.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
              labelcolor=TEXT_COLOR, fontsize=10)

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 7. STRAIN COMPETITION
# ═══════════════════════════════════════════════════════════════════════

def plot_strain_competition(
    result: 'CoupledSimResult',
    n_bins: int = 3,
    disease_year: int = 3,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Stacked area plot of virulence bins (low/medium/high) over time.

    Shows which strain types dominate and how competitive dynamics
    shift between mild and virulent strains.

    Since CoupledSimResult only stores mean virulence, we approximate
    the bin fractions from a synthetic distribution using mean ± SD.

    Args:
        result: CoupledSimResult with virulence data.
        n_bins: Number of virulence bins (default 3: low/med/high).
        disease_year: Year disease was introduced.
        save_path: Optional path to save the figure.

    Returns:
        matplotlib Figure.
    """
    mean_v = result.yearly_mean_virulence
    sd_v = _estimate_virulence_sd(result)
    years = np.arange(result.n_years)

    # Bin boundaries
    bin_edges = np.linspace(0, 1, n_bins + 1)
    bin_labels = ['Low (mild)', 'Medium', 'High (virulent)'] if n_bins == 3 else [
        f'Bin {i+1}' for i in range(n_bins)
    ]
    bin_colors = [VIRULENCE_LOW, ACCENT_COLORS[4], VIRULENCE_HIGH] if n_bins == 3 else [
        _virulence_color((bin_edges[i] + bin_edges[i + 1]) / 2)
        for i in range(n_bins)
    ]

    # Compute bin fractions per year from synthetic distributions
    fractions = np.zeros((n_bins, result.n_years))
    rng = np.random.default_rng(42)

    for yr in range(result.n_years):
        if mean_v[yr] <= 0:
            continue
        samples = rng.normal(mean_v[yr], max(sd_v[yr], 0.01), 500)
        samples = np.clip(samples, 0, 1)
        for b in range(n_bins):
            fractions[b, yr] = np.mean(
                (samples >= bin_edges[b]) & (samples < bin_edges[b + 1])
            )
        # Normalize
        total = fractions[:, yr].sum()
        if total > 0:
            fractions[:, yr] /= total

    fig, ax = dark_figure()

    ax.stackplot(
        years, *fractions,
        labels=bin_labels, colors=bin_colors, alpha=0.85,
    )

    ax.axvline(disease_year, color='white', linestyle=':', linewidth=1.5,
               alpha=0.7, label=f'Disease intro')

    ax.set_xlabel('Year', fontsize=12)
    ax.set_ylabel('Fraction of strains', fontsize=12)
    ax.set_title('Strain Competition: Virulence Bin Composition',
                 fontsize=14, fontweight='bold')
    ax.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
              labelcolor=TEXT_COLOR, fontsize=10, loc='upper right')
    ax.set_xlim(0, result.n_years - 1)
    ax.set_ylim(0, 1)

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 8. COEVOLUTION MULTI-SEED
# ═══════════════════════════════════════════════════════════════════════

def plot_coevolution_multi_seed(
    results_dict: Dict[int, 'CoupledSimResult'],
    disease_year: int = 3,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Multiple seeds overlaid: resistance + virulence trajectories.

    Shows convergence (or lack thereof) to evolutionary equilibrium.
    Thin alpha=0.2 lines per seed, thick mean line.

    Args:
        results_dict: Dict mapping seed → CoupledSimResult.
        disease_year: Year disease was introduced.
        save_path: Optional path to save the figure.

    Returns:
        matplotlib Figure.
    """
    fig, axes = dark_figure(nrows=2, ncols=1, figsize=(12, 10))
    ax_r, ax_v = axes

    seeds = sorted(results_dict.keys())
    n_years = max(r.n_years for r in results_dict.values())
    years = np.arange(n_years)

    # Collect all trajectories
    all_resistance = []
    all_virulence = []

    for seed in seeds:
        r = results_dict[seed]
        ny = r.n_years
        yrs = np.arange(ny)

        # Resistance
        res = r.yearly_mean_resistance
        if res is not None:
            ax_r.plot(yrs, res, color=ACCENT_COLORS[7], alpha=0.2,
                      linewidth=1)
            padded = np.full(n_years, np.nan)
            padded[:ny] = res
            all_resistance.append(padded)

        # Virulence
        vir = r.yearly_mean_virulence
        if vir is not None:
            valid = vir > 0
            vir_plot = np.where(valid, vir, np.nan)
            ax_v.plot(yrs, vir_plot, color=ACCENT_COLORS[6], alpha=0.2,
                      linewidth=1)
            padded = np.full(n_years, np.nan)
            padded[:ny] = vir_plot
            all_virulence.append(padded)

    # Compute and plot means
    if all_resistance:
        mean_r = np.nanmean(np.array(all_resistance), axis=0)
        ax_r.plot(years, mean_r, color=ACCENT_COLORS[7], linewidth=3,
                  label='Mean resistance', zorder=5)

    if all_virulence:
        with np.errstate(all='ignore'):
            mean_v = np.nanmean(np.array(all_virulence), axis=0)
        ax_v.plot(years, mean_v, color=ACCENT_COLORS[6], linewidth=3,
                  label='Mean virulence', zorder=5)

    # Disease introduction line
    for ax in axes:
        ax.axvline(disease_year, color=ACCENT_COLORS[4], linestyle=':',
                   linewidth=1.5, alpha=0.8)

    ax_r.set_ylabel('Mean resistance (r̄)', fontsize=12)
    ax_r.set_title(f'Co-Evolution Across {len(seeds)} Seeds',
                   fontsize=14, fontweight='bold')
    ax_r.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
                labelcolor=TEXT_COLOR, fontsize=10)
    ax_r.set_xlim(0, n_years - 1)

    ax_v.set_xlabel('Year', fontsize=12)
    ax_v.set_ylabel('Mean virulence (v̄)', fontsize=12)
    ax_v.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
                labelcolor=TEXT_COLOR, fontsize=10)
    ax_v.set_xlim(0, n_years - 1)
    ax_v.set_ylim(0, 1)

    # Label seed count
    ax_r.text(0.98, 0.02, f'n = {len(seeds)} seeds',
              transform=ax_r.transAxes, ha='right', va='bottom',
              color=GRID_COLOR, fontsize=10)

    if save_path:
        save_figure(fig, save_path)
    return fig
