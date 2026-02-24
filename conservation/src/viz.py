"""Visualization functions for conservation genetics analyses.

Plotting functions for:
- Trait distribution histograms with ratio framing
- Exceedance curves
- Screening effort plots (expected max vs sample size)
- Breeding trajectory plots (mean/max/VA over generations)
- Gain–diversity Pareto frontier
- Complementarity heatmaps

All functions take data as input (no model dependencies) and return
matplotlib Figure objects. Caller is responsible for fig.savefig()
or plt.show().
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.figure import Figure
from typing import Optional, Sequence


# ═══════════════════════════════════════════════════════════════════════
# STYLE DEFAULTS
# ═══════════════════════════════════════════════════════════════════════

_COLORS = {
    'resistance': '#2196F3',
    'tolerance': '#FF9800',
    'recovery': '#4CAF50',
    'highlight': '#E91E63',
    'neutral': '#9E9E9E',
    'pareto': '#673AB7',
}

_LABEL_MAP = {
    'r': 'Resistance',
    't': 'Tolerance',
    'c': 'Recovery',
    'resistance': 'Resistance',
    'tolerance': 'Tolerance',
    'recovery': 'Recovery',
}


def _apply_style(ax, xlabel='', ylabel='', title=''):
    """Apply consistent styling to an axes."""
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    if title:
        ax.set_title(title, fontsize=13, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(labelsize=10)


# ═══════════════════════════════════════════════════════════════════════
# TRAIT DISTRIBUTION HISTOGRAMS
# ═══════════════════════════════════════════════════════════════════════


def plot_trait_distribution(
    scores: np.ndarray,
    trait_name: str = 'resistance',
    threshold: Optional[float] = None,
    ratio_label: Optional[str] = None,
    bins: int = 50,
    figsize: tuple = (8, 5),
) -> Figure:
    """Histogram of a trait distribution with optional threshold line.

    Uses "ratio framing" — annotates the fraction above/below threshold
    to communicate exceedance probability visually.

    Args:
        scores: (N,) trait scores.
        trait_name: 'resistance', 'tolerance', or 'recovery'.
        threshold: If given, draw a vertical line and annotate exceedance.
        ratio_label: Custom label for the exceedance ratio (e.g., "1 in 50").
        bins: Number of histogram bins.
        figsize: Figure size.

    Returns:
        matplotlib Figure.
    """
    fig, ax = plt.subplots(figsize=figsize)

    color = _COLORS.get(trait_name, _COLORS['neutral'])
    label = _LABEL_MAP.get(trait_name, trait_name.title())

    ax.hist(scores, bins=bins, color=color, alpha=0.7, edgecolor='white',
            linewidth=0.5, density=True, label=f'{label} (n={len(scores)})')

    if threshold is not None:
        ax.axvline(threshold, color=_COLORS['highlight'], linestyle='--',
                   linewidth=2, label=f'Threshold = {threshold:.3f}')

        n_above = np.sum(scores >= threshold)
        frac = n_above / len(scores) if len(scores) > 0 else 0

        if ratio_label is None:
            if frac > 0:
                ratio_label = f'1 in {int(round(1/frac))}' if frac < 1 else 'all'
            else:
                ratio_label = '< 1 in {}'.format(len(scores))

        ax.annotate(
            f'{n_above}/{len(scores)} above ({ratio_label})',
            xy=(threshold, ax.get_ylim()[1] * 0.9),
            fontsize=10, color=_COLORS['highlight'],
            fontweight='bold',
            ha='left', va='top',
        )

    ax.legend(fontsize=10, frameon=False)
    _apply_style(ax, xlabel=f'{label} Score', ylabel='Density',
                 title=f'{label} Distribution')

    fig.tight_layout()
    return fig


def plot_trait_distributions_overlay(
    scores_dict: dict,
    thresholds: Optional[dict] = None,
    bins: int = 50,
    figsize: tuple = (10, 5),
) -> Figure:
    """Overlaid histograms for multiple traits.

    Args:
        scores_dict: {'resistance': array, 'tolerance': array, ...}.
        thresholds: Optional {'resistance': 0.3, ...} threshold lines.
        bins: Number of histogram bins.
        figsize: Figure size.

    Returns:
        matplotlib Figure.
    """
    fig, ax = plt.subplots(figsize=figsize)

    for trait_name, scores in scores_dict.items():
        color = _COLORS.get(trait_name, _COLORS['neutral'])
        label = _LABEL_MAP.get(trait_name, trait_name.title())
        ax.hist(scores, bins=bins, color=color, alpha=0.4, edgecolor='white',
                linewidth=0.5, density=True, label=label)

    if thresholds:
        for trait_name, thresh in thresholds.items():
            color = _COLORS.get(trait_name, _COLORS['neutral'])
            ax.axvline(thresh, color=color, linestyle='--', linewidth=1.5)

    ax.legend(fontsize=10, frameon=False)
    _apply_style(ax, xlabel='Trait Score', ylabel='Density',
                 title='Trait Distributions')
    fig.tight_layout()
    return fig


# ═══════════════════════════════════════════════════════════════════════
# EXCEEDANCE CURVES
# ═══════════════════════════════════════════════════════════════════════


def plot_exceedance(
    thresholds: np.ndarray,
    exceedance_probs: np.ndarray,
    trait_name: str = 'resistance',
    target_prob: Optional[float] = None,
    figsize: tuple = (8, 5),
) -> Figure:
    """Plot exceedance curve P(τ ≥ τ*) vs threshold.

    Args:
        thresholds: (K,) threshold values.
        exceedance_probs: (K,) exceedance probabilities.
        trait_name: Trait being plotted.
        target_prob: If given, draw horizontal line and annotate.
        figsize: Figure size.

    Returns:
        matplotlib Figure.
    """
    fig, ax = plt.subplots(figsize=figsize)

    color = _COLORS.get(trait_name, _COLORS['neutral'])
    label = _LABEL_MAP.get(trait_name, trait_name.title())

    ax.plot(thresholds, exceedance_probs, color=color, linewidth=2)
    ax.fill_between(thresholds, exceedance_probs, alpha=0.15, color=color)

    if target_prob is not None:
        ax.axhline(target_prob, color=_COLORS['highlight'], linestyle=':',
                   linewidth=1.5, label=f'p = {target_prob}')
        # Find threshold at target probability
        idx = np.argmin(np.abs(exceedance_probs - target_prob))
        ax.axvline(thresholds[idx], color=_COLORS['highlight'], linestyle=':',
                   linewidth=1, alpha=0.5)
        ax.annotate(f'τ* ≈ {thresholds[idx]:.3f}',
                    xy=(thresholds[idx], target_prob),
                    fontsize=10, color=_COLORS['highlight'])
        ax.legend(fontsize=10, frameon=False)

    ax.set_yscale('log')
    ax.set_ylim(bottom=max(1e-4, exceedance_probs[exceedance_probs > 0].min() * 0.5))

    _apply_style(ax, xlabel=f'{label} Threshold (τ*)',
                 ylabel='P(τ ≥ τ*)',
                 title=f'{label} Exceedance Curve')
    fig.tight_layout()
    return fig


# ═══════════════════════════════════════════════════════════════════════
# SCREENING EFFORT PLOTS
# ═══════════════════════════════════════════════════════════════════════


def plot_screening_effort(
    sample_sizes: np.ndarray,
    expected_maxes: np.ndarray,
    trait_name: str = 'resistance',
    target_value: Optional[float] = None,
    figsize: tuple = (8, 5),
) -> Figure:
    """Expected maximum vs sample size (diminishing returns).

    Shows logarithmic saturation — relates to Eq. 4.6.

    Args:
        sample_sizes: (K,) sample sizes evaluated.
        expected_maxes: (K,) expected maximum at each sample size.
        trait_name: Trait being plotted.
        target_value: Target trait value line.
        figsize: Figure size.

    Returns:
        matplotlib Figure.
    """
    fig, ax = plt.subplots(figsize=figsize)

    color = _COLORS.get(trait_name, _COLORS['neutral'])
    label = _LABEL_MAP.get(trait_name, trait_name.title())

    ax.plot(sample_sizes, expected_maxes, 'o-', color=color, linewidth=2,
            markersize=4, label=f'E[max {label}]')

    if target_value is not None:
        ax.axhline(target_value, color=_COLORS['highlight'], linestyle='--',
                   linewidth=1.5, label=f'Target = {target_value:.3f}')
        ax.legend(fontsize=10, frameon=False)

    ax.set_xscale('log')
    _apply_style(ax, xlabel='Sample Size (n)',
                 ylabel=f'Expected Best {label}',
                 title='Screening Effort: Diminishing Returns')
    fig.tight_layout()
    return fig


def plot_multisite_allocation(
    site_names: Sequence[str],
    allocations: np.ndarray,
    site_means: np.ndarray,
    figsize: tuple = (10, 5),
) -> Figure:
    """Bar chart of screening allocation across sites.

    Args:
        site_names: (K,) site labels.
        allocations: (K,) samples per site.
        site_means: (K,) site mean traits (for coloring).
        figsize: Figure size.

    Returns:
        matplotlib Figure.
    """
    fig, ax = plt.subplots(figsize=figsize)

    # Color by mean trait value
    norm = plt.Normalize(site_means.min(), site_means.max())
    cmap = plt.cm.Blues
    colors = [cmap(norm(m)) for m in site_means]

    bars = ax.bar(site_names, allocations, color=colors, edgecolor='white',
                  linewidth=0.5)

    # Annotate with mean trait values
    for bar, mean_val in zip(bars, site_means):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                f'μ={mean_val:.3f}', ha='center', va='bottom', fontsize=8)

    _apply_style(ax, xlabel='Site', ylabel='Samples Allocated',
                 title='Optimal Screening Allocation')

    ax.tick_params(axis='x', rotation=45)
    fig.tight_layout()
    return fig


# ═══════════════════════════════════════════════════════════════════════
# BREEDING TRAJECTORY PLOTS
# ═══════════════════════════════════════════════════════════════════════


def plot_breeding_trajectory(
    generations: np.ndarray,
    mean_r: np.ndarray,
    max_r: np.ndarray,
    va_r: np.ndarray,
    he: np.ndarray,
    f_mean: np.ndarray,
    loci_fixed: Optional[np.ndarray] = None,
    n_loci_total: int = 17,
    figsize: tuple = (14, 10),
) -> Figure:
    """Multi-panel breeding trajectory plot.

    Four panels: (1) Mean & max resistance, (2) V_A,
    (3) Heterozygosity & F, (4) Loci fixed (optional).

    Args:
        generations: (G,) generation numbers.
        mean_r: (G,) mean resistance per generation.
        max_r: (G,) max resistance per generation.
        va_r: (G,) additive variance per generation.
        he: (G,) expected heterozygosity.
        f_mean: (G,) mean inbreeding coefficient.
        loci_fixed: (G,) number of fixed loci.
        n_loci_total: Total loci for trait (for % fixed axis).
        figsize: Figure size.

    Returns:
        matplotlib Figure.
    """
    n_panels = 4 if loci_fixed is not None else 3
    fig, axes = plt.subplots(n_panels, 1, figsize=figsize, sharex=True)

    # Panel 1: Trait values
    ax = axes[0]
    ax.plot(generations, mean_r, 'o-', color=_COLORS['resistance'],
            linewidth=2, markersize=5, label='Mean')
    ax.plot(generations, max_r, 's--', color=_COLORS['highlight'],
            linewidth=1.5, markersize=4, label='Max')
    ax.fill_between(generations, mean_r, max_r,
                    color=_COLORS['resistance'], alpha=0.1)
    ax.legend(fontsize=10, frameon=False)
    _apply_style(ax, ylabel='Resistance Score', title='Breeding Program Trajectory')

    # Panel 2: Additive variance
    ax = axes[1]
    ax.plot(generations, va_r, 'o-', color=_COLORS['recovery'],
            linewidth=2, markersize=5)
    ax.fill_between(generations, 0, va_r, color=_COLORS['recovery'], alpha=0.15)
    _apply_style(ax, ylabel='Additive Variance (V_A)')

    # Panel 3: Diversity
    ax = axes[2]
    ax.plot(generations, he, 'o-', color=_COLORS['tolerance'],
            linewidth=2, markersize=5, label='H_e')
    ax2 = ax.twinx()
    ax2.plot(generations, f_mean, 's--', color=_COLORS['highlight'],
             linewidth=1.5, markersize=4, label='Mean F')
    ax.legend(loc='upper left', fontsize=10, frameon=False)
    ax2.legend(loc='upper right', fontsize=10, frameon=False)
    _apply_style(ax, ylabel='Expected Heterozygosity')
    ax2.set_ylabel('Inbreeding (F)', fontsize=11, color=_COLORS['highlight'])
    ax2.tick_params(labelsize=10)

    # Panel 4: Loci fixed
    if loci_fixed is not None:
        ax = axes[3]
        pct_fixed = loci_fixed / n_loci_total * 100
        ax.bar(generations, pct_fixed, color=_COLORS['neutral'], alpha=0.7,
               edgecolor='white', width=0.8)
        _apply_style(ax, xlabel='Generation', ylabel='Loci Fixed (%)')
        ax.set_ylim(0, 100)

    if loci_fixed is None:
        axes[-1].set_xlabel('Generation', fontsize=11)

    fig.tight_layout()
    return fig


def plot_breeding_comparison(
    results: dict,
    metric: str = 'mean_r',
    figsize: tuple = (10, 6),
) -> Figure:
    """Compare multiple breeding schemes on one plot.

    Args:
        results: {scheme_name: {'generations': array, metric: array}}.
        metric: Which metric to plot ('mean_r', 'max_r', 'va_r', 'he', 'f_mean').
        figsize: Figure size.

    Returns:
        matplotlib Figure.
    """
    fig, ax = plt.subplots(figsize=figsize)

    colors = list(plt.cm.Set2.colors)
    for i, (name, data) in enumerate(results.items()):
        color = colors[i % len(colors)]
        ax.plot(data['generations'], data[metric], 'o-', color=color,
                linewidth=2, markersize=5, label=name)

    ax.legend(fontsize=10, frameon=False)
    _apply_style(ax, xlabel='Generation', ylabel=metric.replace('_', ' ').title(),
                 title=f'Breeding Scheme Comparison: {metric}')
    fig.tight_layout()
    return fig


# ═══════════════════════════════════════════════════════════════════════
# GAIN–DIVERSITY PARETO FRONTIER
# ═══════════════════════════════════════════════════════════════════════


def plot_pareto_frontier(
    gain: np.ndarray,
    diversity: np.ndarray,
    labels: Optional[Sequence[str]] = None,
    pareto_idx: Optional[np.ndarray] = None,
    xlabel: str = 'Genetic Gain (Δτ)',
    ylabel: str = 'Diversity Retained (H_e)',
    figsize: tuple = (8, 6),
) -> Figure:
    """Plot gain vs diversity with Pareto frontier highlighted.

    Args:
        gain: (K,) genetic gain values.
        diversity: (K,) diversity metric values.
        labels: (K,) labels for each point.
        pareto_idx: Indices of Pareto-optimal points. If None, computed.
        xlabel: X-axis label.
        ylabel: Y-axis label.
        figsize: Figure size.

    Returns:
        matplotlib Figure.
    """
    fig, ax = plt.subplots(figsize=figsize)

    if pareto_idx is None:
        pareto_idx = _compute_pareto(gain, diversity)

    # Non-Pareto points
    non_pareto = np.ones(len(gain), dtype=bool)
    non_pareto[pareto_idx] = False

    ax.scatter(gain[non_pareto], diversity[non_pareto],
               color=_COLORS['neutral'], s=60, alpha=0.5,
               edgecolors='white', zorder=2, label='Dominated')

    # Pareto points
    sorted_idx = pareto_idx[np.argsort(gain[pareto_idx])]
    ax.scatter(gain[pareto_idx], diversity[pareto_idx],
               color=_COLORS['pareto'], s=100, edgecolors='white',
               zorder=3, label='Pareto-optimal')
    ax.plot(gain[sorted_idx], diversity[sorted_idx],
            '--', color=_COLORS['pareto'], linewidth=1.5, alpha=0.7, zorder=2)

    if labels is not None:
        for i in pareto_idx:
            ax.annotate(labels[i], (gain[i], diversity[i]),
                        fontsize=9, ha='left', va='bottom',
                        xytext=(5, 5), textcoords='offset points')

    ax.legend(fontsize=10, frameon=False)
    _apply_style(ax, xlabel=xlabel, ylabel=ylabel,
                 title='Gain–Diversity Trade-off')
    fig.tight_layout()
    return fig


def _compute_pareto(gain: np.ndarray, diversity: np.ndarray) -> np.ndarray:
    """Find Pareto-optimal indices (maximize both gain and diversity)."""
    n = len(gain)
    is_pareto = np.ones(n, dtype=bool)

    for i in range(n):
        if not is_pareto[i]:
            continue
        for j in range(n):
            if i == j or not is_pareto[j]:
                continue
            # j dominates i if j is better in both
            if gain[j] >= gain[i] and diversity[j] >= diversity[i]:
                if gain[j] > gain[i] or diversity[j] > diversity[i]:
                    is_pareto[i] = False
                    break

    return np.where(is_pareto)[0]


# ═══════════════════════════════════════════════════════════════════════
# COMPLEMENTARITY HEATMAPS
# ═══════════════════════════════════════════════════════════════════════


def plot_complementarity_heatmap(
    matrix: np.ndarray,
    labels: Optional[Sequence[str]] = None,
    metric_name: str = 'Locus Union',
    cmap: str = 'YlOrRd',
    figsize: tuple = (10, 8),
) -> Figure:
    """Heatmap of pairwise complementarity/union/overlap.

    Args:
        matrix: (N, N) int or float pairwise matrix.
        labels: (N,) individual labels.
        metric_name: Name for the colorbar.
        cmap: Colormap name.
        figsize: Figure size.

    Returns:
        matplotlib Figure.
    """
    fig, ax = plt.subplots(figsize=figsize)

    n = len(matrix)
    if labels is None:
        labels = [str(i) for i in range(n)]

    im = ax.imshow(matrix, cmap=cmap, aspect='equal')
    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label(metric_name, fontsize=11)

    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
    ax.set_yticklabels(labels, fontsize=8)

    # Annotate cells if matrix is small enough
    if n <= 20:
        for i in range(n):
            for j in range(n):
                val = matrix[i, j]
                text_color = 'white' if val > (matrix.max() + matrix.min()) / 2 else 'black'
                ax.text(j, i, str(int(val)), ha='center', va='center',
                        fontsize=7, color=text_color)

    _apply_style(ax, title=f'Pairwise {metric_name}')
    fig.tight_layout()
    return fig


def plot_kinship_heatmap(
    kinship_matrix: np.ndarray,
    labels: Optional[Sequence[str]] = None,
    figsize: tuple = (10, 8),
) -> Figure:
    """Heatmap of pairwise kinship coefficients.

    Args:
        kinship_matrix: (N, N) float kinship/relationship matrix.
        labels: (N,) individual labels.
        figsize: Figure size.

    Returns:
        matplotlib Figure.
    """
    fig, ax = plt.subplots(figsize=figsize)

    n = len(kinship_matrix)
    if labels is None:
        labels = [str(i) for i in range(n)]

    # Diverging colormap centered at 0
    vmax = max(abs(kinship_matrix.min()), abs(kinship_matrix.max()))
    im = ax.imshow(kinship_matrix, cmap='RdBu_r', aspect='equal',
                   vmin=-vmax, vmax=vmax)
    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label('Kinship Coefficient', fontsize=11)

    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
    ax.set_yticklabels(labels, fontsize=8)

    _apply_style(ax, title='Genomic Kinship Matrix')
    fig.tight_layout()
    return fig
