"""Salinity & freshwater influence visualizations for SSWD-EvoEpi.

Comprehensive plots for the seasonal salinity mechanism:
  - Spatial maps of salinity depression and disease suppression
  - Temporal profiles showing seasonal cycles by region/latitude
  - Parameter sensitivity (fw_strength, fw_depth_exp)
  - Validation against DFO lighthouse observations
  - Asymmetry diagnostics (AK vs CA gradient)
  - Integration with disease module (sal_mod impact)

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
from typing import Dict, List, Optional, Sequence, Tuple, TYPE_CHECKING

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches
import matplotlib.cm as cm
import numpy as np

from sswd_evoepi.viz.style import (
    ACCENT_COLORS,
    DARK_BG,
    DARK_PANEL,
    GRID_COLOR,
    LIGHT_BG,
    LIGHT_GRID,
    LIGHT_PANEL,
    LIGHT_TEXT,
    NODE_COLORS,
    TEXT_COLOR,
    apply_dark_theme,
    apply_light_theme,
    dark_figure,
    pub_figure,
    save_figure,
    themed_figure,
    _theme_colors,
)
from sswd_evoepi.salinity import (
    ocean_baseline,
    freshwater_melt_pulse,
    latitude_melt_factor,
    compute_salinity_array,
    DAYS_PER_YEAR,
    _PEAK_DAY,
)

if TYPE_CHECKING:
    from sswd_evoepi.spatial import NodeDefinition


# ═══════════════════════════════════════════════════════════════════════
# COLOUR CONSTANTS
# ═══════════════════════════════════════════════════════════════════════

FRESH_WATER = '#3498db'     # blue — freshwater
SALT_WATER = '#e74c3c'      # red — high salinity / disease
SUPPRESSION = '#2ecc71'     # green — disease suppression
NEUTRAL = '#f39c12'         # amber — neutral / baseline
FJORD_COLOR = '#48c9b0'     # teal — fjord sites
OPEN_COLOR = '#e94560'      # crimson — open coast

# Region colour map (consistent with calibration plots)
REGION_COLORS = {
    'AK-AL': '#1f77b4', 'AK-WG': '#2ca02c', 'AK-OC': '#9467bd',
    'AK-EG': '#8c564b', 'AK-PWS': '#e377c2', 'AK-FN': '#17becf',
    'AK-FS': '#bcbd22', 'BC-N': '#ff7f0e', 'BC-C': '#d62728',
    'SS-N': '#7f7f7f', 'SS-S': '#aec7e8', 'JDF': '#ffbb78',
    'WA-O': '#98df8a', 'OR': '#ff9896', 'CA-N': '#c5b0d5',
    'CA-C': '#c49c94', 'CA-S': '#f7b6d2', 'BJ': '#dbdb8d',
}

# Month labels
MONTH_LABELS = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
MONTH_STARTS = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]


# ═══════════════════════════════════════════════════════════════════════
# HELPER: salinity modifier (mirrors disease.py salinity_modifier)
# ═══════════════════════════════════════════════════════════════════════

def _sal_mod(salinity: float, s_min: float = 10.0, s_full: float = 28.0,
             eta: float = 2.0) -> float:
    """Compute Vibrio salinity modifier in [0, 1].

    Mirrors disease.py::salinity_modifier().
    """
    if salinity >= s_full:
        return 1.0
    if salinity <= s_min:
        return 0.0
    x = (salinity - s_min) / (s_full - s_min)
    return x ** eta




def _tc(theme='light'):
    """Get text color for theme."""
    return LIGHT_TEXT if theme == 'light' else TEXT_COLOR

def _pc(theme='light'):
    """Get panel color for theme."""
    return LIGHT_PANEL if theme == 'light' else DARK_PANEL

def _gc(theme='light'):
    """Get grid color for theme."""
    return LIGHT_GRID if theme == 'light' else GRID_COLOR

def _sal_mod_array(salinity_arr: np.ndarray, s_min: float = 10.0,
                   s_full: float = 28.0, eta: float = 2.0) -> np.ndarray:
    """Vectorized salinity modifier."""
    x = np.clip((salinity_arr - s_min) / (s_full - s_min), 0.0, 1.0)
    return x ** eta


# ═══════════════════════════════════════════════════════════════════════
# 1. SALINITY SEASONAL HEATMAP (nodes × day-of-year)
# ═══════════════════════════════════════════════════════════════════════

def plot_salinity_heatmap(
    nodes: List['NodeDefinition'],
    fw_strength: float = 15.0,
    fw_depth_exp: float = 1.0,
    sort_by: str = 'latitude',
    title: Optional[str] = None,
    theme: str = 'light',
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Heatmap of daily salinity for all nodes across a year.

    Rows = nodes (sorted by latitude or fjord_depth_norm).
    Columns = day of year (0-364).
    Color = salinity (psu).

    Args:
        nodes: List of NodeDefinition objects.
        fw_strength: Freshwater depression strength (psu).
        fw_depth_exp: Exponent on fjord_depth_norm.
        sort_by: 'latitude' or 'fjord_depth' — row ordering.
        title: Plot title (auto-generated if None).
        save_path: Optional path to save figure.

    Returns:
        matplotlib Figure.
    """
    sal = compute_salinity_array(nodes, fw_strength, fw_depth_exp=fw_depth_exp)
    N = len(nodes)

    # Sort order
    if sort_by == 'fjord_depth':
        order = np.argsort([nd.fjord_depth_norm for nd in nodes])[::-1]
        ylabel = 'Node (sorted by fjord depth)'
    else:
        order = np.argsort([nd.lat for nd in nodes])[::-1]  # north at top
        ylabel = 'Node (sorted by latitude, N→S)'

    sal_sorted = sal[order, :]
    names_sorted = [f"{nodes[i].name} ({nodes[i].lat:.1f}°N)" for i in order]

    fig, ax = themed_figure(theme=theme, figsize=(14, max(6, N * 0.12)))

    im = ax.imshow(
        sal_sorted, aspect='auto', cmap='YlGnBu_r',
        interpolation='nearest',
        vmin=max(5.0, sal.min()), vmax=sal.max(),
    )

    # Month ticks on x-axis
    ax.set_xticks(MONTH_STARTS)
    ax.set_xticklabels(MONTH_LABELS, fontsize=9)

    # Node labels on y-axis (subsample if >30 nodes)
    if N <= 40:
        ax.set_yticks(range(N))
        ax.set_yticklabels(names_sorted, fontsize=7)
    else:
        step = max(1, N // 25)
        ticks = list(range(0, N, step))
        ax.set_yticks(ticks)
        ax.set_yticklabels([names_sorted[t] for t in ticks], fontsize=7)

    ax.set_xlabel('Month', fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)

    if title is None:
        title = f'Seasonal Salinity (fw_strength={fw_strength}, exp={fw_depth_exp})'
    ax.set_title(title, fontsize=14, fontweight='bold')

    cbar = fig.colorbar(im, ax=ax, pad=0.02, shrink=0.8)
    cbar.set_label('Salinity (psu)', color=_tc(theme), fontsize=11)
    cbar.ax.yaxis.set_tick_params(color=_tc(theme))
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color=_tc(theme))

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 2. SALINITY DEPRESSION HEATMAP (shows Δ from baseline)
# ═══════════════════════════════════════════════════════════════════════

def plot_depression_heatmap(
    nodes: List['NodeDefinition'],
    fw_strength: float = 15.0,
    fw_depth_exp: float = 1.0,
    theme: str = 'light',
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Heatmap of salinity DEPRESSION (baseline - actual) across nodes × months.

    Highlights where and when freshwater influence is strongest.
    Blue = strong depression, white = no depression.

    Args:
        nodes: List of NodeDefinition objects.
        fw_strength: Freshwater depression strength.
        fw_depth_exp: Exponent on fjord_depth_norm.
        save_path: Optional path to save figure.

    Returns:
        matplotlib Figure.
    """
    sal = compute_salinity_array(nodes, fw_strength, fw_depth_exp=fw_depth_exp)
    sal_baseline = compute_salinity_array(nodes, fw_strength=0.0)
    depression = sal_baseline - sal  # positive = more depression
    N = len(nodes)

    # Monthly averages
    monthly_depression = np.zeros((N, 12), dtype=np.float32)
    for m in range(12):
        start = MONTH_STARTS[m]
        end = MONTH_STARTS[m + 1] if m < 11 else 365
        monthly_depression[:, m] = depression[:, start:end].mean(axis=1)

    # Sort by max depression (most depressed at top)
    max_dep = monthly_depression.max(axis=1)
    order = np.argsort(max_dep)[::-1]
    monthly_sorted = monthly_depression[order, :]
    names_sorted = [f"{nodes[i].name}" for i in order]

    fig, ax = themed_figure(theme=theme, figsize=(12, max(6, N * 0.12)))

    im = ax.imshow(
        monthly_sorted, aspect='auto', cmap='Blues',
        interpolation='nearest', vmin=0, vmax=max(monthly_depression.max(), 1.0),
    )

    ax.set_xticks(range(12))
    ax.set_xticklabels(MONTH_LABELS, fontsize=10)

    if N <= 40:
        ax.set_yticks(range(N))
        ax.set_yticklabels(names_sorted, fontsize=7)
    else:
        step = max(1, N // 25)
        ticks = list(range(0, N, step))
        ax.set_yticks(ticks)
        ax.set_yticklabels([names_sorted[t] for t in ticks], fontsize=7)

    ax.set_xlabel('Month', fontsize=12)
    ax.set_ylabel('Node (sorted by max depression)', fontsize=12)
    ax.set_title(f'Salinity Depression Δ (fw_strength={fw_strength})',
                 fontsize=14, fontweight='bold')

    cbar = fig.colorbar(im, ax=ax, pad=0.02, shrink=0.8)
    cbar.set_label('Depression (psu)', color=_tc(theme), fontsize=11)
    cbar.ax.yaxis.set_tick_params(color=_tc(theme))
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color=_tc(theme))

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 3. DISEASE SUPPRESSION MAP (geographic scatter)
# ═══════════════════════════════════════════════════════════════════════

def plot_suppression_map(
    nodes: List['NodeDefinition'],
    fw_strength: float = 15.0,
    fw_depth_exp: float = 1.0,
    month: int = 5,  # June (0-indexed)
    s_min: float = 10.0,
    s_full: float = 28.0,
    title: Optional[str] = None,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Geographic map of disease suppression from salinity.

    Each node plotted at lat/lon, coloured by salinity modifier (sal_mod)
    for a given month. Green = strong suppression, red = no suppression.

    Args:
        nodes: List of NodeDefinition objects.
        fw_strength: Freshwater depression strength.
        fw_depth_exp: Exponent on fjord_depth_norm.
        month: Month index (0=Jan, 5=Jun, etc.).
        s_min, s_full: Salinity modifier parameters.
        title: Plot title.
        save_path: Optional path to save figure.

    Returns:
        matplotlib Figure.
    """
    sal = compute_salinity_array(nodes, fw_strength, fw_depth_exp=fw_depth_exp)
    N = len(nodes)

    # Monthly mean salinity for the target month
    start = MONTH_STARTS[month]
    end = MONTH_STARTS[month + 1] if month < 11 else 365
    monthly_sal = sal[:, start:end].mean(axis=1)

    # Compute salinity modifier
    sal_mod = _sal_mod_array(monthly_sal, s_min=s_min, s_full=s_full)
    suppression = 1.0 - sal_mod  # 0 = no suppression, 1 = full suppression

    lons = np.array([nd.lon for nd in nodes])
    lats = np.array([nd.lat for nd in nodes])

    fig, ax = themed_figure(theme=theme, figsize=(10, 12))

    # Custom diverging colourmap: red (no suppression) → green (full)
    cmap = mcolors.LinearSegmentedColormap.from_list(
        'suppression', [SALT_WATER, NEUTRAL, SUPPRESSION], N=256
    )

    scatter = ax.scatter(
        lons, lats, c=suppression, cmap=cmap,
        vmin=0, vmax=max(suppression.max(), 0.1),
        s=80, edgecolors='white', linewidths=0.5, zorder=3,
    )

    ax.set_xlabel('Longitude', fontsize=12)
    ax.set_ylabel('Latitude (°N)', fontsize=12)

    if title is None:
        title = (f'Disease Suppression from Salinity — {MONTH_LABELS[month]}\n'
                 f'(fw_strength={fw_strength}, exp={fw_depth_exp})')
    ax.set_title(title, fontsize=14, fontweight='bold')

    cbar = fig.colorbar(scatter, ax=ax, pad=0.02, shrink=0.7)
    cbar.set_label('Disease suppression (1 − sal_mod)', color=_tc(theme), fontsize=11)
    cbar.ax.yaxis.set_tick_params(color=_tc(theme))
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color=_tc(theme))

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 4. SEASONAL SALINITY PROFILES BY REGION
# ═══════════════════════════════════════════════════════════════════════

def plot_regional_salinity_profiles(
    nodes: List['NodeDefinition'],
    fw_strength: float = 15.0,
    fw_depth_exp: float = 1.0,
    regions: Optional[List[str]] = None,
    theme: str = 'light',
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Seasonal salinity curves averaged per region.

    One line per region showing mean daily salinity through the year.
    Highlights the latitude-dependent melt pulse.

    Args:
        nodes: List of NodeDefinition objects.
        fw_strength: Freshwater depression strength.
        fw_depth_exp: Exponent on fjord_depth_norm.
        regions: Subset of regions to plot (None = auto-select key regions).
        save_path: Optional path to save figure.

    Returns:
        matplotlib Figure.
    """
    sal = compute_salinity_array(nodes, fw_strength, fw_depth_exp=fw_depth_exp)
    N = len(nodes)

    # Group nodes by region
    region_nodes: Dict[str, List[int]] = {}
    for i, nd in enumerate(nodes):
        region_nodes.setdefault(nd.subregion, []).append(i)

    if regions is None:
        # Show key gradient regions
        regions = ['AK-PWS', 'AK-FN', 'BC-N', 'SS-S', 'JDF', 'OR', 'CA-N']
        regions = [r for r in regions if r in region_nodes]

    days = np.arange(DAYS_PER_YEAR)
    fig, ax = themed_figure(theme=theme, figsize=(12, 7))

    for region in regions:
        idxs = region_nodes.get(region, [])
        if not idxs:
            continue
        mean_sal = sal[idxs, :].mean(axis=0)
        color = REGION_COLORS.get(region, ACCENT_COLORS[len(regions) % len(ACCENT_COLORS)])
        ax.plot(days, mean_sal, color=color, linewidth=2.2, label=region, alpha=0.9)

        # Shade ±1 SD if more than 3 nodes
        if len(idxs) > 3:
            std_sal = sal[idxs, :].std(axis=0)
            ax.fill_between(days, mean_sal - std_sal, mean_sal + std_sal,
                            color=color, alpha=0.1)

    # Mark June peak
    ax.axvline(_PEAK_DAY, color='white', linestyle=':', linewidth=1, alpha=0.5)
    ax.text(_PEAK_DAY + 3, ax.get_ylim()[0] + 0.5, 'June 15\n(peak melt)',
            color=_tc(theme), fontsize=9, alpha=0.7)

    ax.set_xticks(MONTH_STARTS)
    ax.set_xticklabels(MONTH_LABELS, fontsize=10)
    ax.set_xlabel('Month', fontsize=12)
    ax.set_ylabel('Salinity (psu)', fontsize=12)
    ax.set_title(f'Regional Seasonal Salinity (fw_strength={fw_strength})',
                 fontsize=14, fontweight='bold')
    ax.legend(facecolor=_pc(theme), edgecolor=_gc(theme),
              labelcolor=_tc(theme), fontsize=10, loc='lower left')

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 5. LATITUDE VS PEAK DEPRESSION (asymmetry diagnostic)
# ═══════════════════════════════════════════════════════════════════════

def plot_latitude_asymmetry(
    nodes: List['NodeDefinition'],
    fw_strength: float = 15.0,
    fw_depth_exp: float = 1.0,
    theme: str = 'light',
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Scatter: latitude vs peak salinity depression at each node.

    Key diagnostic for the AK-CA asymmetry. Should show depression
    concentrated at high-latitude, high-fjord-depth nodes.

    Two panels:
      Left: latitude vs peak depression, coloured by fjord_depth_norm
      Right: latitude vs peak disease suppression (1 - sal_mod)

    Args:
        nodes: List of NodeDefinition objects.
        fw_strength: Freshwater depression strength.
        fw_depth_exp: Exponent on fjord_depth_norm.
        save_path: Optional path to save figure.

    Returns:
        matplotlib Figure.
    """
    sal = compute_salinity_array(nodes, fw_strength, fw_depth_exp=fw_depth_exp)
    sal_baseline = compute_salinity_array(nodes, fw_strength=0.0)
    depression = (sal_baseline - sal).max(axis=1)  # peak depression per node

    # Peak-month salinity modifier
    peak_sal = sal[:, _PEAK_DAY]
    sal_mod = _sal_mod_array(peak_sal)
    suppression = 1.0 - sal_mod

    lats = np.array([nd.lat for nd in nodes])
    fd = np.array([nd.fjord_depth_norm for nd in nodes])

    fig, (ax1, ax2) = themed_figure(theme=theme, nrows=1, ncols=2, figsize=(16, 7))

    # Left: depression
    sc1 = ax1.scatter(lats, depression, c=fd, cmap='viridis',
                      s=50, edgecolors='white', linewidths=0.3, alpha=0.8)
    ax1.set_xlabel('Latitude (°N)', fontsize=12)
    ax1.set_ylabel('Peak salinity depression (psu)', fontsize=12)
    ax1.set_title('Freshwater Depression by Latitude', fontsize=13, fontweight='bold')
    cb1 = fig.colorbar(sc1, ax=ax1, pad=0.02, shrink=0.8)
    cb1.set_label('fjord_depth_norm', color=_tc(theme), fontsize=10)
    cb1.ax.yaxis.set_tick_params(color=_tc(theme))
    plt.setp(cb1.ax.yaxis.get_ticklabels(), color=_tc(theme))

    # Trend line
    mask = depression > 0.1
    if mask.sum() > 5:
        z = np.polyfit(lats[mask], depression[mask], 1)
        x_fit = np.linspace(35, 62, 50)
        ax1.plot(x_fit, np.polyval(z, x_fit), color=FRESH_WATER,
                 linestyle='--', linewidth=1.5, alpha=0.7)

    # Right: disease suppression
    sc2 = ax2.scatter(lats, suppression * 100, c=fd, cmap='viridis',
                      s=50, edgecolors='white', linewidths=0.3, alpha=0.8)
    ax2.set_xlabel('Latitude (°N)', fontsize=12)
    ax2.set_ylabel('Peak disease suppression (%)', fontsize=12)
    ax2.set_title('Disease Suppression by Latitude', fontsize=13, fontweight='bold')
    cb2 = fig.colorbar(sc2, ax=ax2, pad=0.02, shrink=0.8)
    cb2.set_label('fjord_depth_norm', color=_tc(theme), fontsize=10)
    cb2.ax.yaxis.set_tick_params(color=_tc(theme))
    plt.setp(cb2.ax.yaxis.get_ticklabels(), color=_tc(theme))

    # Annotate asymmetry ratio
    ak_mask = lats > 55
    ca_mask = lats < 40
    if ak_mask.any() and ca_mask.any():
        ak_mean = suppression[ak_mask].mean() * 100
        ca_mean = suppression[ca_mask].mean() * 100
        ratio = ak_mean / ca_mean if ca_mean > 0.01 else float('inf')
        ax2.text(0.02, 0.95, f'AK mean: {ak_mean:.1f}%\nCA mean: {ca_mean:.1f}%\n'
                 f'Asymmetry: {ratio:.0f}×',
                 transform=ax2.transAxes, fontsize=10, color=SUPPRESSION,
                 verticalalignment='top',
                 bbox=dict(boxstyle='round', facecolor=_pc(theme),
                           edgecolor=_gc(theme), alpha=0.9))

    fig.suptitle(f'Latitude-Dependent Salinity Asymmetry (fw_strength={fw_strength})',
                 fontsize=15, fontweight='bold', color=_tc(theme), y=1.02)

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 6. FW_STRENGTH SENSITIVITY (multi-panel dose-response)
# ═══════════════════════════════════════════════════════════════════════

def plot_fw_strength_sensitivity(
    nodes: List['NodeDefinition'],
    fw_values: Optional[List[float]] = None,
    fw_depth_exp: float = 1.0,
    representative_nodes: Optional[Dict[str, int]] = None,
    theme: str = 'light',
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Show how fw_strength affects salinity and disease suppression.

    Three panels:
      Top: June salinity vs fw_strength for representative nodes
      Middle: June disease suppression vs fw_strength
      Bottom: Regional mean suppression gradient

    Args:
        nodes: List of NodeDefinition objects.
        fw_values: List of fw_strength values to sweep. Default [0, 5, ..., 30].
        fw_depth_exp: Exponent on fjord_depth_norm.
        representative_nodes: Dict mapping label → node index for
            individual node traces. Auto-selected if None.
        save_path: Optional path to save figure.

    Returns:
        matplotlib Figure.
    """
    if fw_values is None:
        fw_values = [0, 5, 8, 10, 12, 15, 18, 20, 25, 30]

    N = len(nodes)

    # Auto-select representative nodes (high/mid/low fjord depth × latitude)
    if representative_nodes is None:
        representative_nodes = {}
        # Find AK fjord, BC fjord, WA coast, CA coast
        for i, nd in enumerate(nodes):
            if nd.lat > 56 and nd.fjord_depth_norm > 0.7 and 'AK' not in representative_nodes:
                representative_nodes['AK fjord'] = i
            elif 49 < nd.lat < 52 and nd.fjord_depth_norm > 0.4 and 'BC fjord' not in representative_nodes:
                representative_nodes['BC fjord'] = i
            elif 47 < nd.lat < 49 and nd.fjord_depth_norm < 0.1 and 'WA coast' not in representative_nodes:
                representative_nodes['WA coast'] = i
            elif nd.lat < 38 and 'CA coast' not in representative_nodes:
                representative_nodes['CA coast'] = i
            if len(representative_nodes) >= 4:
                break

    fig, (ax1, ax2, ax3) = themed_figure(theme=theme, nrows=3, ncols=1, figsize=(12, 14))

    rep_colors = [FJORD_COLOR, FRESH_WATER, NEUTRAL, SALT_WATER, '#9b59b6']

    # Sweep fw_strength values
    for ci, (label, idx) in enumerate(representative_nodes.items()):
        june_sal = []
        june_sup = []
        for fw in fw_values:
            sal = compute_salinity_array([nodes[idx]], fw, fw_depth_exp=fw_depth_exp)
            s_june = float(sal[0, _PEAK_DAY])
            june_sal.append(s_june)
            june_sup.append((1.0 - _sal_mod(s_june)) * 100)

        color = rep_colors[ci % len(rep_colors)]
        nd = nodes[idx]
        lbl = f'{label} ({nd.lat:.1f}°N, fd={nd.fjord_depth_norm:.2f})'

        ax1.plot(fw_values, june_sal, color=color, linewidth=2.2,
                 marker='o', markersize=5, label=lbl)
        ax2.plot(fw_values, june_sup, color=color, linewidth=2.2,
                 marker='s', markersize=5, label=lbl)

    ax1.set_ylabel('June 15 salinity (psu)', fontsize=12)
    ax1.set_title('Salinity Dose-Response', fontsize=13, fontweight='bold')
    ax1.legend(facecolor=_pc(theme), edgecolor=_gc(theme),
               labelcolor=_tc(theme), fontsize=9, loc='best')
    ax1.axhline(28, color=_gc(theme), linestyle=':', alpha=0.5)
    ax1.text(fw_values[-1], 28.5, 's_full=28', color=_gc(theme), fontsize=8)

    ax2.set_ylabel('June disease suppression (%)', fontsize=12)
    ax2.set_title('Disease Suppression Dose-Response', fontsize=13, fontweight='bold')
    ax2.legend(facecolor=_pc(theme), edgecolor=_gc(theme),
               labelcolor=_tc(theme), fontsize=9, loc='best')

    # Bottom: regional mean suppression by fw_strength
    region_nodes: Dict[str, List[int]] = {}
    for i, nd in enumerate(nodes):
        region_nodes.setdefault(nd.subregion, []).append(i)

    key_regions = ['AK-PWS', 'AK-FN', 'BC-N', 'OR', 'CA-N']
    key_regions = [r for r in key_regions if r in region_nodes]

    for region in key_regions:
        idxs = region_nodes[region]
        mean_sup = []
        for fw in fw_values:
            sal = compute_salinity_array(nodes, fw, fw_depth_exp=fw_depth_exp)
            june_sal_region = sal[idxs, _PEAK_DAY].mean()
            mean_sup.append((1.0 - _sal_mod(float(june_sal_region))) * 100)

        color = REGION_COLORS.get(region, GRID_COLOR)
        ax3.plot(fw_values, mean_sup, color=color, linewidth=2.2,
                 marker='D', markersize=4, label=region)

    ax3.set_xlabel('fw_strength (psu)', fontsize=12)
    ax3.set_ylabel('Regional mean June suppression (%)', fontsize=12)
    ax3.set_title('Regional Suppression Sensitivity', fontsize=13, fontweight='bold')
    ax3.legend(facecolor=_pc(theme), edgecolor=_gc(theme),
               labelcolor=_tc(theme), fontsize=9, loc='best')

    fig.suptitle('Freshwater Strength Sensitivity Analysis',
                 fontsize=15, fontweight='bold', color=_tc(theme), y=1.01)

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 7. FW_DEPTH_EXP COMPARISON (linear vs sqrt vs square)
# ═══════════════════════════════════════════════════════════════════════

def plot_depth_exp_comparison(
    nodes: List['NodeDefinition'],
    fw_strength: float = 15.0,
    exponents: Optional[List[float]] = None,
    theme: str = 'light',
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Compare different fw_depth_exp values on the spatial gradient.

    Two panels:
      Left: fjord_depth_norm vs effective depth (fd^exp) for each exponent
      Right: latitude vs June suppression for each exponent

    Args:
        nodes: List of NodeDefinition objects.
        fw_strength: Freshwater depression strength.
        exponents: List of exponents to compare. Default [0.5, 1.0, 2.0].
        save_path: Optional path to save figure.

    Returns:
        matplotlib Figure.
    """
    if exponents is None:
        exponents = [0.5, 1.0, 2.0]

    lats = np.array([nd.lat for nd in nodes])
    fd = np.array([nd.fjord_depth_norm for nd in nodes])

    exp_colors = ['#3498db', '#f39c12', '#e74c3c', '#2ecc71', '#9b59b6']

    fig, (ax1, ax2) = themed_figure(theme=theme, nrows=1, ncols=2, figsize=(16, 7))

    # Left: transform curve
    x = np.linspace(0, 1, 100)
    for ci, exp in enumerate(exponents):
        color = exp_colors[ci % len(exp_colors)]
        ax1.plot(x, x ** exp, color=color, linewidth=2.5,
                 label=f'exp={exp}')

    ax1.set_xlabel('fjord_depth_norm (raw)', fontsize=12)
    ax1.set_ylabel('Effective depth (fd^exp)', fontsize=12)
    ax1.set_title('Depth Exponent Transform', fontsize=13, fontweight='bold')
    ax1.legend(facecolor=_pc(theme), edgecolor=_gc(theme),
               labelcolor=_tc(theme), fontsize=11)
    ax1.plot([0, 1], [0, 1], color=_gc(theme), linestyle=':', alpha=0.5)

    # Right: spatial impact
    for ci, exp in enumerate(exponents):
        sal = compute_salinity_array(nodes, fw_strength, fw_depth_exp=exp)
        june_sal = sal[:, _PEAK_DAY]
        sup = (1.0 - _sal_mod_array(june_sal)) * 100
        color = exp_colors[ci % len(exp_colors)]
        ax2.scatter(lats, sup, color=color, s=20, alpha=0.5,
                    label=f'exp={exp}')

    ax2.set_xlabel('Latitude (°N)', fontsize=12)
    ax2.set_ylabel('June peak suppression (%)', fontsize=12)
    ax2.set_title(f'Spatial Impact (fw_strength={fw_strength})',
                  fontsize=13, fontweight='bold')
    ax2.legend(facecolor=_pc(theme), edgecolor=_gc(theme),
               labelcolor=_tc(theme), fontsize=10, loc='best')

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 8. SAL_MOD TRANSFER FUNCTION
# ═══════════════════════════════════════════════════════════════════════

def plot_sal_mod_transfer(
    s_min: float = 10.0,
    s_full: float = 28.0,
    eta: float = 2.0,
    annotate_sites: Optional[Dict[str, float]] = None,
    theme: str = 'light',
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Plot the salinity → disease modifier transfer function.

    Shows how salinity maps to Vibrio activity (sal_mod). Optionally
    annotates specific site salinities.

    Args:
        s_min: Minimum salinity for Vibrio.
        s_full: Full-activity salinity.
        eta: Exponent for the power law.
        annotate_sites: Dict of {label: salinity} to mark on curve.
        save_path: Optional path to save figure.

    Returns:
        matplotlib Figure.
    """
    fig, ax = themed_figure(theme=theme, figsize=(10, 6))

    sal_range = np.linspace(0, 35, 500)
    mod = _sal_mod_array(sal_range, s_min=s_min, s_full=s_full, eta=eta)

    ax.plot(sal_range, mod, color=SALT_WATER, linewidth=3)
    ax.fill_between(sal_range, mod, 0, color=SALT_WATER, alpha=0.1)
    ax.fill_between(sal_range, mod, 1, color=SUPPRESSION, alpha=0.1)

    ax.axvline(s_min, color=FRESH_WATER, linestyle='--', linewidth=1, alpha=0.7)
    ax.axvline(s_full, color=SALT_WATER, linestyle='--', linewidth=1, alpha=0.7)
    ax.text(s_min - 1.5, 0.9, f's_min={s_min}', color=FRESH_WATER, fontsize=9,
            rotation=90)
    ax.text(s_full + 0.5, 0.9, f's_full={s_full}', color=SALT_WATER, fontsize=9,
            rotation=90)

    # Annotate sites
    if annotate_sites:
        for label, sal_val in annotate_sites.items():
            y_val = _sal_mod(sal_val, s_min, s_full, eta)
            ax.plot(sal_val, y_val, 'o', color='white', markersize=8, zorder=5)
            ax.annotate(
                f'{label}\n{sal_val:.1f} psu → {y_val:.2f}',
                (sal_val, y_val), xytext=(15, 15),
                textcoords='offset points', fontsize=9, color=_tc(theme),
                arrowprops=dict(arrowstyle='->', color=_gc(theme), lw=1),
                bbox=dict(boxstyle='round', facecolor=_pc(theme),
                          edgecolor=_gc(theme), alpha=0.9),
            )

    ax.set_xlabel('Salinity (psu)', fontsize=12)
    ax.set_ylabel('sal_mod (Vibrio activity)', fontsize=12)
    ax.set_title(f'Salinity → Disease Modifier Transfer Function (η={eta})',
                 fontsize=14, fontweight='bold')
    ax.set_xlim(0, 35)
    ax.set_ylim(-0.02, 1.05)

    # Dual annotation
    ax.text(5, 0.5, 'Disease\nSuppressed', color=SUPPRESSION,
            fontsize=14, fontweight='bold', alpha=0.5, ha='center')
    ax.text(31, 0.5, 'Full\nDisease', color=SALT_WATER,
            fontsize=14, fontweight='bold', alpha=0.5, ha='center')

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 9. MELT PULSE & LATITUDE FACTOR (mechanism components)
# ═══════════════════════════════════════════════════════════════════════

def plot_mechanism_components(
    theme: str = 'light',
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Three-panel plot showing the three multiplicative components.

    Top: freshwater_melt_pulse(day) — seasonal cosine
    Middle: latitude_melt_factor(lat) — latitude ramp
    Bottom: combined effect (fw_strength × fd_norm × f_melt × pulse)
        for a few example nodes

    Args:
        save_path: Optional path to save figure.

    Returns:
        matplotlib Figure.
    """
    fig, (ax1, ax2, ax3) = themed_figure(theme=theme, nrows=3, ncols=1, figsize=(12, 12))

    days = np.arange(DAYS_PER_YEAR)
    pulse = np.array([freshwater_melt_pulse(d) for d in days])

    # Top: melt pulse
    ax1.plot(days, pulse, color=FRESH_WATER, linewidth=2.5)
    ax1.fill_between(days, pulse, 0, color=FRESH_WATER, alpha=0.15)
    ax1.axvline(_PEAK_DAY, color='white', linestyle=':', linewidth=1, alpha=0.5)
    ax1.set_xticks(MONTH_STARTS)
    ax1.set_xticklabels(MONTH_LABELS)
    ax1.set_ylabel('Melt pulse', fontsize=12)
    ax1.set_title('Seasonal Meltwater Pulse', fontsize=13, fontweight='bold')
    ax1.set_ylim(-0.05, 1.1)
    ax1.text(_PEAK_DAY + 5, 0.95, 'June 15 peak', color=_tc(theme), fontsize=10)

    # Middle: latitude factor
    lats = np.linspace(30, 65, 200)
    f_melt = np.array([latitude_melt_factor(lat) for lat in lats])
    ax2.plot(lats, f_melt, color=FJORD_COLOR, linewidth=2.5)
    ax2.fill_between(lats, f_melt, 0, color=FJORD_COLOR, alpha=0.15)
    ax2.axhline(0, color=_gc(theme), linewidth=0.5)
    ax2.axhline(1, color=_gc(theme), linewidth=0.5, linestyle=':')
    ax2.axvline(35, color=SALT_WATER, linestyle='--', linewidth=1, alpha=0.5)
    ax2.axvline(60, color=SUPPRESSION, linestyle='--', linewidth=1, alpha=0.5)
    ax2.text(35.5, 0.85, 'CA: zero', color=SALT_WATER, fontsize=10)
    ax2.text(57, 0.85, 'AK: full', color=SUPPRESSION, fontsize=10)
    ax2.set_xlabel('Latitude (°N)', fontsize=12)
    ax2.set_ylabel('Latitude melt factor', fontsize=12)
    ax2.set_title('Latitude-Dependent Melt Amplitude', fontsize=13, fontweight='bold')

    # Bottom: combined depression for example nodes
    examples = [
        ('AK fjord (58°N, fd=0.8)', 58.0, 0.8, SUPPRESSION),
        ('BC fjord (50°N, fd=0.5)', 50.0, 0.5, FJORD_COLOR),
        ('WA coast (48°N, fd=0.0)', 48.0, 0.0, NEUTRAL),
        ('CA coast (36°N, fd=0.3)', 36.0, 0.3, SALT_WATER),
    ]

    fw_strength = 15.0
    for label, lat, fd_norm, color in examples:
        base = ocean_baseline(lat)
        f_m = latitude_melt_factor(lat)
        depression = fw_strength * fd_norm * f_m * pulse
        sal = np.maximum(5.0, base - depression)
        ax3.plot(days, sal, color=color, linewidth=2.2, label=label)

    ax3.axhline(28, color=_gc(theme), linestyle=':', alpha=0.5)
    ax3.text(5, 28.3, 's_full', color=_gc(theme), fontsize=8)
    ax3.axhline(10, color=_gc(theme), linestyle=':', alpha=0.5)
    ax3.text(5, 10.3, 's_min', color=_gc(theme), fontsize=8)

    ax3.set_xticks(MONTH_STARTS)
    ax3.set_xticklabels(MONTH_LABELS)
    ax3.set_xlabel('Month', fontsize=12)
    ax3.set_ylabel('Salinity (psu)', fontsize=12)
    ax3.set_title(f'Combined Salinity Profiles (fw_strength={fw_strength})',
                  fontsize=13, fontweight='bold')
    ax3.legend(facecolor=_pc(theme), edgecolor=_gc(theme),
               labelcolor=_tc(theme), fontsize=10, loc='lower left')

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 10. VALIDATION: MODEL VS DFO LIGHTHOUSE OBSERVATIONS
# ═══════════════════════════════════════════════════════════════════════

def plot_dfo_validation(
    dfo_csv_path: str,
    nodes: Optional[List['NodeDefinition']] = None,
    fw_strength: float = 15.0,
    fw_depth_exp: float = 1.0,
    theme: str = 'light',
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Compare model salinity predictions against DFO lighthouse observations.

    Loads DFO monthly climatology CSV and overlays model predictions for
    matching latitudes.

    Two panels:
      Top: monthly salinity for observed stations + model predictions
      Bottom: residual (model - observed) by month

    Args:
        dfo_csv_path: Path to dfo_monthly_climatology.csv.
        nodes: Optional list of NodeDefinition for model predictions.
            If None, creates synthetic nodes at DFO station latitudes.
        fw_strength: Model fw_strength.
        fw_depth_exp: Model fw_depth_exp.
        save_path: Optional path to save figure.

    Returns:
        matplotlib Figure.
    """
    import csv

    # Load DFO data
    stations = {}  # name → {lat, months: [12 values]}
    with open(dfo_csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row.get('station', row.get('Station', row.get('name', '?')))
            lat = float(row.get('latitude', row.get('Latitude', 49.0)))
            months_data = []
            for m in range(1, 13):
                key_candidates = [f'month_{m}', f'Month_{m}', f'm{m}',
                                  MONTH_LABELS[m-1], str(m)]
                val = None
                for k in key_candidates:
                    if k in row and row[k]:
                        try:
                            val = float(row[k])
                        except ValueError:
                            pass
                        break
                months_data.append(val)

            if any(v is not None for v in months_data):
                stations[name] = {'lat': lat, 'months': months_data}

    if not stations:
        fig, ax = themed_figure(theme=theme, )
        ax.text(0.5, 0.5, f'No data loaded from {dfo_csv_path}',
                ha='center', va='center', color=_tc(theme), fontsize=14,
                transform=ax.transAxes)
        if save_path:
            save_figure(fig, save_path)
        return fig

    fig, (ax1, ax2) = themed_figure(theme=theme, nrows=2, ncols=1, figsize=(14, 10))
    months_x = np.arange(12)

    from sswd_evoepi.spatial import NodeDefinition

    for si, (name, data) in enumerate(sorted(stations.items(),
                                              key=lambda x: -x[1]['lat'])):
        color = NODE_COLORS[si % len(NODE_COLORS)]
        lat = data['lat']
        obs = data['months']

        # Plot observed
        valid = [(m, v) for m, v in enumerate(obs) if v is not None]
        if valid:
            mx, mv = zip(*valid)
            ax1.plot(mx, mv, color=color, linewidth=2, marker='o',
                     markersize=5, label=f'{name} ({lat:.1f}°N) obs')

        # Model prediction at same latitude (assume open coast — fd=0)
        # DFO stations are on headlands, so fd_norm ≈ 0
        syn_node = NodeDefinition(
            node_id=0, name=name, lat=lat, lon=-125.0,
            subregion='DFO', habitat_area=1000, carrying_capacity=100,
            fjord_depth_norm=0.0,
        )
        sal = compute_salinity_array([syn_node], fw_strength=fw_strength,
                                     fw_depth_exp=fw_depth_exp)
        # Monthly mean from model
        model_monthly = []
        for m in range(12):
            start = MONTH_STARTS[m]
            end = MONTH_STARTS[m + 1] if m < 11 else 365
            model_monthly.append(float(sal[0, start:end].mean()))

        ax1.plot(months_x, model_monthly, color=color, linewidth=1.5,
                 linestyle='--', marker='s', markersize=3, alpha=0.7)

        # Residuals
        residuals = []
        for m in range(12):
            if obs[m] is not None:
                residuals.append(model_monthly[m] - obs[m])
            else:
                residuals.append(None)
        valid_res = [(m, r) for m, r in enumerate(residuals) if r is not None]
        if valid_res:
            rx, rv = zip(*valid_res)
            ax2.plot(rx, rv, color=color, linewidth=1.5, marker='o',
                     markersize=4, label=name)

    ax1.set_xticks(range(12))
    ax1.set_xticklabels(MONTH_LABELS)
    ax1.set_ylabel('Salinity (psu)', fontsize=12)
    ax1.set_title('Model vs DFO Lighthouse Observations', fontsize=13,
                  fontweight='bold')
    ax1.legend(facecolor=_pc(theme), edgecolor=_gc(theme),
               labelcolor=_tc(theme), fontsize=8, ncol=2, loc='lower left')

    ax2.axhline(0, color='white', linewidth=0.8, alpha=0.5)
    ax2.set_xticks(range(12))
    ax2.set_xticklabels(MONTH_LABELS)
    ax2.set_xlabel('Month', fontsize=12)
    ax2.set_ylabel('Residual (model − obs, psu)', fontsize=12)
    ax2.set_title('Residuals', fontsize=13, fontweight='bold')
    ax2.legend(facecolor=_pc(theme), edgecolor=_gc(theme),
               labelcolor=_tc(theme), fontsize=8, ncol=2, loc='best')

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 11. REGIONAL SUPPRESSION BAR CHART
# ═══════════════════════════════════════════════════════════════════════

def plot_regional_suppression_bars(
    nodes: List['NodeDefinition'],
    fw_strength: float = 15.0,
    fw_depth_exp: float = 1.0,
    month: int = 5,  # June
    theme: str = 'light',
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Bar chart of mean disease suppression per region.

    Sorted by latitude (north to south). Shows the N-S gradient
    in a single clear figure.

    Args:
        nodes: List of NodeDefinition objects.
        fw_strength: Freshwater depression strength.
        fw_depth_exp: Exponent on fjord_depth_norm.
        month: Month index (0=Jan, 5=Jun).
        save_path: Optional path to save figure.

    Returns:
        matplotlib Figure.
    """
    sal = compute_salinity_array(nodes, fw_strength, fw_depth_exp=fw_depth_exp)

    # Monthly mean salinity
    start = MONTH_STARTS[month]
    end = MONTH_STARTS[month + 1] if month < 11 else 365
    monthly_sal = sal[:, start:end].mean(axis=1)
    sal_mod = _sal_mod_array(monthly_sal)
    suppression = (1.0 - sal_mod) * 100

    # Group by region
    region_data: Dict[str, Dict] = {}
    for i, nd in enumerate(nodes):
        r = nd.subregion
        if r not in region_data:
            region_data[r] = {'lats': [], 'sup': []}
        region_data[r]['lats'].append(nd.lat)
        region_data[r]['sup'].append(suppression[i])

    # Sort regions by mean latitude (north to south)
    region_stats = []
    for r, d in region_data.items():
        region_stats.append({
            'region': r,
            'mean_lat': np.mean(d['lats']),
            'mean_sup': np.mean(d['sup']),
            'std_sup': np.std(d['sup']),
            'n_nodes': len(d['sup']),
        })
    region_stats.sort(key=lambda x: -x['mean_lat'])

    fig, ax = themed_figure(theme=theme, figsize=(14, 7))

    x = np.arange(len(region_stats))
    bars = ax.bar(
        x,
        [s['mean_sup'] for s in region_stats],
        yerr=[s['std_sup'] for s in region_stats],
        capsize=3,
        color=[REGION_COLORS.get(s['region'], ACCENT_COLORS[0])
               for s in region_stats],
        edgecolor='white', linewidth=0.5,
    )

    ax.set_xticks(x)
    ax.set_xticklabels(
        [f"{s['region']}\n({s['mean_lat']:.0f}°N, n={s['n_nodes']})"
         for s in region_stats],
        fontsize=9,
    )
    ax.set_ylabel(f'{MONTH_LABELS[month]} disease suppression (%)', fontsize=12)
    ax.set_title(f'Regional Disease Suppression from Salinity\n'
                 f'(fw_strength={fw_strength}, exp={fw_depth_exp})',
                 fontsize=14, fontweight='bold')
    ax.set_ylim(bottom=0)

    # Annotate bars
    for i, s in enumerate(region_stats):
        if s['mean_sup'] > 1:
            ax.text(i, s['mean_sup'] + s['std_sup'] + 1,
                    f"{s['mean_sup']:.1f}%", ha='center', fontsize=8,
                    color=_tc(theme))

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 12. FJORD DEPTH DISTRIBUTION BY REGION
# ═══════════════════════════════════════════════════════════════════════

def plot_fjord_depth_by_region(
    nodes: List['NodeDefinition'],
    theme: str = 'light',
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Box plot of fjord_depth_norm by region.

    Shows why salinity creates a latitude gradient: AK has many deep
    fjord sites (high fd_norm), CA has few.

    Args:
        nodes: List of NodeDefinition objects.
        save_path: Optional path to save figure.

    Returns:
        matplotlib Figure.
    """
    region_data: Dict[str, List[float]] = {}
    region_lats: Dict[str, float] = {}
    for nd in nodes:
        r = nd.subregion
        region_data.setdefault(r, []).append(nd.fjord_depth_norm)
        region_lats.setdefault(r, [])
        region_lats[r].append(nd.lat)

    # Sort by mean latitude
    regions_sorted = sorted(region_data.keys(),
                            key=lambda r: -np.mean(region_lats[r]))

    fig, ax = themed_figure(theme=theme, figsize=(14, 7))

    bp = ax.boxplot(
        [region_data[r] for r in regions_sorted],
        labels=[f'{r}\n({np.mean(region_lats[r]):.0f}°N)' for r in regions_sorted],
        patch_artist=True,
        medianprops=dict(color='white', linewidth=2),
        whiskerprops=dict(color=_tc(theme)),
        capprops=dict(color=_tc(theme)),
        flierprops=dict(markerfacecolor=_gc(theme), markersize=3),
    )

    for i, (patch, region) in enumerate(zip(bp['boxes'], regions_sorted)):
        color = REGION_COLORS.get(region, ACCENT_COLORS[i % len(ACCENT_COLORS)])
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
        patch.set_edgecolor('white')

    ax.set_ylabel('fjord_depth_norm', fontsize=12)
    ax.set_title('Fjord Depth Distribution by Region',
                 fontsize=14, fontweight='bold')
    ax.tick_params(axis='x', labelsize=8)

    # Annotate median
    for i, r in enumerate(regions_sorted):
        med = np.median(region_data[r])
        n = len(region_data[r])
        ax.text(i + 1, ax.get_ylim()[1] * 0.95,
                f'n={n}', ha='center', fontsize=7, color=_tc(theme))

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 13. MULTI-MONTH SUPPRESSION PANELS (spatial snapshots)
# ═══════════════════════════════════════════════════════════════════════

def plot_suppression_monthly_panels(
    nodes: List['NodeDefinition'],
    fw_strength: float = 15.0,
    fw_depth_exp: float = 1.0,
    months: Optional[List[int]] = None,
    theme: str = 'light',
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Multi-panel geographic maps of disease suppression for key months.

    Shows the spatial pattern of suppression evolving through the year.

    Args:
        nodes: List of NodeDefinition objects.
        fw_strength: Freshwater depression strength.
        fw_depth_exp: Exponent on fjord_depth_norm.
        months: List of month indices. Default [1, 3, 5, 7, 9, 11]
            (Feb, Apr, Jun, Aug, Oct, Dec).
        save_path: Optional path to save figure.

    Returns:
        matplotlib Figure.
    """
    if months is None:
        months = [1, 3, 5, 7, 9, 11]  # Feb, Apr, Jun, Aug, Oct, Dec

    sal = compute_salinity_array(nodes, fw_strength, fw_depth_exp=fw_depth_exp)
    lons = np.array([nd.lon for nd in nodes])
    lats = np.array([nd.lat for nd in nodes])

    n_panels = len(months)
    ncols = min(3, n_panels)
    nrows = int(np.ceil(n_panels / ncols))

    fig, axes = themed_figure(theme=theme, nrows=nrows, ncols=ncols,
                            figsize=(6 * ncols, 8 * nrows))
    if isinstance(axes, np.ndarray):
        axes_flat = axes.flat
    else:
        axes_flat = [axes]

    cmap = mcolors.LinearSegmentedColormap.from_list(
        'suppression', [SALT_WATER, NEUTRAL, SUPPRESSION], N=256
    )

    for idx, (ax, month) in enumerate(zip(axes_flat, months)):
        start = MONTH_STARTS[month]
        end = MONTH_STARTS[month + 1] if month < 11 else 365
        monthly_sal = sal[:, start:end].mean(axis=1)
        sup = (1.0 - _sal_mod_array(monthly_sal)) * 100

        vmax = max(sup.max(), 5.0)
        sc = ax.scatter(lons, lats, c=sup, cmap=cmap,
                        vmin=0, vmax=vmax,
                        s=30, edgecolors='white', linewidths=0.2, alpha=0.85)
        ax.set_title(MONTH_LABELS[month], fontsize=13, fontweight='bold')
        ax.set_xlabel('Lon', fontsize=9)
        ax.set_ylabel('Lat', fontsize=9)

        if idx == n_panels - 1 or idx == ncols - 1:
            cb = fig.colorbar(sc, ax=ax, pad=0.02, shrink=0.7)
            cb.set_label('Suppression %', color=_tc(theme), fontsize=9)
            cb.ax.yaxis.set_tick_params(color=_tc(theme))
            plt.setp(cb.ax.yaxis.get_ticklabels(), color=_tc(theme))

    # Hide unused axes
    for idx in range(n_panels, nrows * ncols):
        if idx < len(list(axes_flat)):
            list(axes_flat)[idx].set_visible(False)

    fig.suptitle(f'Seasonal Disease Suppression (fw={fw_strength}, exp={fw_depth_exp})',
                 fontsize=15, fontweight='bold', color=_tc(theme), y=1.02)

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 14. CALIBRATION IMPACT: SUPPRESSION × RECOVERY OVERLAY
# ═══════════════════════════════════════════════════════════════════════

def plot_suppression_vs_recovery(
    nodes: List['NodeDefinition'],
    region_recovery: Dict[str, float],
    fw_strength: float = 15.0,
    fw_depth_exp: float = 1.0,
    recovery_targets: Optional[Dict[str, float]] = None,
    theme: str = 'light',
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Scatter: regional mean suppression vs observed/simulated recovery.

    Key diagnostic: if salinity is the missing mechanism, high-suppression
    regions should have higher recovery.

    Args:
        nodes: List of NodeDefinition objects.
        region_recovery: Dict {region: recovery_fraction} from calibration.
        fw_strength: Freshwater depression strength.
        fw_depth_exp: Exponent on fjord_depth_norm.
        recovery_targets: Optional target recovery fractions for comparison.
        save_path: Optional path to save figure.

    Returns:
        matplotlib Figure.
    """
    sal = compute_salinity_array(nodes, fw_strength, fw_depth_exp=fw_depth_exp)
    june_sal = sal[:, _PEAK_DAY]

    # Regional means
    region_nodes: Dict[str, List[int]] = {}
    for i, nd in enumerate(nodes):
        region_nodes.setdefault(nd.subregion, []).append(i)

    fig, ax = themed_figure(theme=theme, figsize=(10, 8))

    regions_both = set(region_nodes.keys()) & set(region_recovery.keys())

    for region in sorted(regions_both):
        idxs = region_nodes[region]
        mean_sal = float(june_sal[idxs].mean())
        sup = (1.0 - _sal_mod(mean_sal)) * 100
        rec = region_recovery[region] * 100

        color = REGION_COLORS.get(region, ACCENT_COLORS[0])
        ax.scatter(sup, rec, color=color, s=200, edgecolors='white',
                   linewidths=1.5, zorder=3)
        ax.annotate(
            region, (sup, rec), xytext=(8, 8), textcoords='offset points',
            fontsize=10, color=_tc(theme), fontweight='bold',
            bbox=dict(boxstyle='round', facecolor=_pc(theme),
                      edgecolor=_gc(theme), alpha=0.8),
        )

        # Target if available
        if recovery_targets and region in recovery_targets:
            target_rec = recovery_targets[region] * 100
            ax.scatter(sup, target_rec, color=color, s=80,
                       marker='x', linewidths=2, zorder=4, alpha=0.6)

    # Trend line
    sups = []
    recs = []
    for region in regions_both:
        idxs = region_nodes[region]
        mean_sal = float(june_sal[idxs].mean())
        sups.append((1.0 - _sal_mod(mean_sal)) * 100)
        recs.append(region_recovery[region] * 100)
    sups = np.array(sups)
    recs = np.array(recs)

    if len(sups) > 3 and np.std(sups) > 1e-6 and np.std(recs) > 1e-6:
        try:
            z = np.polyfit(sups, recs, 1)
            x_fit = np.linspace(0, max(sups) * 1.1, 50)
            ax.plot(x_fit, np.polyval(z, x_fit), color=_gc(theme),
                    linestyle='--', linewidth=1.5, alpha=0.7)
            r_corr = np.corrcoef(sups, recs)[0, 1]
            ax.text(0.95, 0.05, f'r = {r_corr:.3f}', transform=ax.transAxes,
                    fontsize=12, color=_tc(theme), ha='right',
                    bbox=dict(boxstyle='round', facecolor=_pc(theme),
                              edgecolor=_gc(theme)))
        except (np.linalg.LinAlgError, ValueError):
            pass  # Polyfit can fail with edge-case data

    if recovery_targets:
        ax.plot([], [], 'x', color=_gc(theme), markersize=8, label='Target')
        ax.legend(facecolor=_pc(theme), edgecolor=_gc(theme),
                  labelcolor=_tc(theme), fontsize=10)

    ax.set_xlabel('June disease suppression from salinity (%)', fontsize=12)
    ax.set_ylabel('Population recovery (%)', fontsize=12)
    ax.set_title('Salinity Suppression vs Recovery',
                 fontsize=14, fontweight='bold')

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# REGION BOUNDING BOXES for zoom panels
# ═══════════════════════════════════════════════════════════════════════

REGION_BOUNDS = {
    'Alaska': {'lat': (54.5, 62.0), 'lon': (-155, -130), 'regions': ['AK-AL', 'AK-WG', 'AK-OC', 'AK-EG', 'AK-PWS', 'AK-FN', 'AK-FS']},
    'British Columbia': {'lat': (48.0, 55.0), 'lon': (-134, -122), 'regions': ['BC-N', 'BC-C']},
    'Salish Sea & WA': {'lat': (46.5, 49.5), 'lon': (-126, -121.5), 'regions': ['SS-N', 'SS-S', 'JDF', 'WA-O']},
    'Oregon & California': {'lat': (32.0, 47.0), 'lon': (-126, -116), 'regions': ['OR', 'CA-N', 'CA-C', 'CA-S', 'BJ']},
}


# ═══════════════════════════════════════════════════════════════════════
# 15. REGIONAL ZOOM: SUPPRESSION MAP (4-panel geographic zooms)
# ═══════════════════════════════════════════════════════════════════════

def plot_suppression_regional_zoom(
    nodes: List['NodeDefinition'],
    fw_strength: float = 15.0,
    fw_depth_exp: float = 1.0,
    month: int = 5,  # June
    s_min: float = 10.0,
    s_full: float = 28.0,
    theme: str = 'light',
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Four-panel regional zoom of disease suppression.

    Each panel zooms into a geographic region (Alaska, BC, Salish Sea,
    OR/CA) showing individual site suppression at high resolution.

    Args:
        nodes: List of NodeDefinition objects.
        fw_strength: Freshwater depression strength.
        fw_depth_exp: Exponent on fjord_depth_norm.
        month: Month index (0=Jan, 5=Jun).
        s_min, s_full: Salinity modifier parameters.
        theme: 'dark' or 'light'.
        save_path: Optional path to save figure.

    Returns:
        matplotlib Figure.
    """
    sal = compute_salinity_array(nodes, fw_strength, fw_depth_exp=fw_depth_exp)
    N = len(nodes)

    start = MONTH_STARTS[month]
    end = MONTH_STARTS[month + 1] if month < 11 else 365
    monthly_sal = sal[:, start:end].mean(axis=1)
    sal_mod = _sal_mod_array(monthly_sal, s_min=s_min, s_full=s_full)
    suppression = (1.0 - sal_mod) * 100

    lons = np.array([nd.lon for nd in nodes])
    lats = np.array([nd.lat for nd in nodes])
    names = [nd.name for nd in nodes]
    regions_arr = [nd.subregion for nd in nodes]

    bg, panel, text, grid = _theme_colors(theme)

    cmap = mcolors.LinearSegmentedColormap.from_list(
        'suppression', [SALT_WATER, NEUTRAL, SUPPRESSION], N=256
    )

    region_list = list(REGION_BOUNDS.items())
    fig, axes = themed_figure(theme=theme, nrows=2, ncols=2, figsize=(16, 14))

    for ax, (region_name, bounds) in zip(axes.flat, region_list):
        lat_min, lat_max = bounds['lat']
        lon_min, lon_max = bounds['lon']
        valid_regions = set(bounds['regions'])

        # Filter nodes in this region
        mask = np.array([
            (lat_min <= nd.lat <= lat_max and lon_min <= nd.lon <= lon_max)
            for nd in nodes
        ])
        idx = np.where(mask)[0]

        if len(idx) == 0:
            ax.text(0.5, 0.5, 'No nodes', ha='center', va='center',
                    color=text, transform=ax.transAxes)
            ax.set_title(region_name, fontsize=13, fontweight='bold')
            continue

        vmax = max(suppression[idx].max(), 5.0)
        sc = ax.scatter(
            lons[idx], lats[idx], c=suppression[idx], cmap=cmap,
            vmin=0, vmax=vmax,
            s=60, edgecolors='white' if theme == 'dark' else '#333333',
            linewidths=0.4, alpha=0.85, zorder=3,
        )

        # Label top-5 most suppressed sites in this region
        region_sup = suppression[idx]
        top5 = np.argsort(region_sup)[-5:]
        for ti in top5:
            gi = idx[ti]  # global index
            if suppression[gi] > 1.0:
                short = nodes[gi].name.split(',')[0][:15]
                ax.annotate(
                    f'{short}\n{suppression[gi]:.0f}%',
                    (lons[gi], lats[gi]),
                    xytext=(6, 6), textcoords='offset points',
                    fontsize=7, color=text,
                    bbox=dict(boxstyle='round,pad=0.15', facecolor=panel,
                              edgecolor=grid, alpha=0.85),
                )

        ax.set_xlim(lon_min, lon_max)
        ax.set_ylim(lat_min, lat_max)
        ax.set_xlabel('Longitude', fontsize=10)
        ax.set_ylabel('Latitude (°N)', fontsize=10)
        ax.set_title(f'{region_name}  ({len(idx)} sites)', fontsize=13,
                     fontweight='bold')

        cb = fig.colorbar(sc, ax=ax, pad=0.02, shrink=0.8)
        cb.set_label('Suppression %', color=text, fontsize=9)
        cb.ax.yaxis.set_tick_params(color=text)
        plt.setp(cb.ax.yaxis.get_ticklabels(), color=text)

    fig.suptitle(
        f'Regional Disease Suppression — {MONTH_LABELS[month]}\n'
        f'(fw_strength={fw_strength}, exp={fw_depth_exp})',
        fontsize=15, fontweight='bold', color=text, y=1.02,
    )

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 16. REGIONAL ZOOM: SALINITY PROFILES (4 panels, per-site lines)
# ═══════════════════════════════════════════════════════════════════════

def plot_salinity_regional_profiles(
    nodes: List['NodeDefinition'],
    fw_strength: float = 15.0,
    fw_depth_exp: float = 1.0,
    theme: str = 'light',
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Four-panel seasonal salinity profiles, one panel per region.

    Each panel shows individual site salinity curves (thin lines) and
    the regional mean (thick line), revealing within-region variability.

    Args:
        nodes: List of NodeDefinition objects.
        fw_strength: Freshwater depression strength.
        fw_depth_exp: Exponent on fjord_depth_norm.
        theme: 'dark' or 'light'.
        save_path: Optional path to save figure.

    Returns:
        matplotlib Figure.
    """
    sal = compute_salinity_array(nodes, fw_strength, fw_depth_exp=fw_depth_exp)
    days = np.arange(DAYS_PER_YEAR)
    bg, panel, text, grid = _theme_colors(theme)

    region_list = list(REGION_BOUNDS.items())
    fig, axes = themed_figure(theme=theme, nrows=2, ncols=2, figsize=(16, 12))

    for ax, (region_name, bounds) in zip(axes.flat, region_list):
        lat_min, lat_max = bounds['lat']
        lon_min, lon_max = bounds['lon']

        mask = np.array([
            (lat_min <= nd.lat <= lat_max and lon_min <= nd.lon <= lon_max)
            for nd in nodes
        ])
        idx = np.where(mask)[0]

        if len(idx) == 0:
            ax.set_title(region_name, fontsize=13, fontweight='bold')
            continue

        # Individual site lines (thin, low alpha)
        for i in idx:
            fd = nodes[i].fjord_depth_norm
            # Colour by fjord depth: blue=deep fjord, grey=open coast
            alpha = 0.15 + 0.3 * fd
            color = FRESH_WATER if fd > 0.3 else (grid if theme == 'dark' else '#aaaaaa')
            ax.plot(days, sal[i, :], color=color, linewidth=0.5, alpha=alpha)

        # Regional mean (thick)
        mean_sal = sal[idx, :].mean(axis=0)
        ax.plot(days, mean_sal, color=ACCENT_COLORS[0], linewidth=2.5,
                label=f'Mean (n={len(idx)})')

        # Highlight most-depressed site
        min_june = sal[idx, _PEAK_DAY]
        most_depressed = idx[np.argmin(min_june)]
        md_name = nodes[most_depressed].name.split(',')[0][:20]
        ax.plot(days, sal[most_depressed, :], color=FRESH_WATER, linewidth=2,
                linestyle='--', label=f'{md_name} (fd={nodes[most_depressed].fjord_depth_norm:.2f})')

        # Reference lines
        ax.axhline(28, color=grid, linestyle=':', alpha=0.5, linewidth=0.8)
        ax.axhline(10, color=grid, linestyle=':', alpha=0.5, linewidth=0.8)
        ax.axvline(_PEAK_DAY, color=grid, linestyle=':', alpha=0.3)

        ax.set_xticks(MONTH_STARTS)
        ax.set_xticklabels(MONTH_LABELS, fontsize=8)
        ax.set_ylabel('Salinity (psu)', fontsize=10)
        ax.set_title(region_name, fontsize=13, fontweight='bold')
        ax.legend(fontsize=8, loc='lower left',
                  facecolor=panel, edgecolor=grid, labelcolor=text)
        ax.set_ylim(bottom=max(0, sal[idx, :].min() - 2))

    fig.suptitle(
        f'Regional Salinity Profiles (fw_strength={fw_strength})',
        fontsize=15, fontweight='bold', color=text, y=1.02,
    )

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# 17. REGIONAL ZOOM: DEPRESSION SCATTER (lat vs depression, per region)
# ═══════════════════════════════════════════════════════════════════════

def plot_depression_regional_scatter(
    nodes: List['NodeDefinition'],
    fw_strength: float = 15.0,
    fw_depth_exp: float = 1.0,
    theme: str = 'light',
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Four-panel scatter: fjord_depth_norm vs peak depression, per region.

    Shows within-region structure — which individual sites get the most
    freshwater influence.

    Args:
        nodes: List of NodeDefinition objects.
        fw_strength: Freshwater depression strength.
        fw_depth_exp: Exponent on fjord_depth_norm.
        theme: 'dark' or 'light'.
        save_path: Optional path to save figure.

    Returns:
        matplotlib Figure.
    """
    sal = compute_salinity_array(nodes, fw_strength, fw_depth_exp=fw_depth_exp)
    sal_base = compute_salinity_array(nodes, fw_strength=0.0)
    peak_dep = (sal_base - sal).max(axis=1)  # max depression per node

    fd = np.array([nd.fjord_depth_norm for nd in nodes])
    bg, panel, text, grid = _theme_colors(theme)

    region_list = list(REGION_BOUNDS.items())
    fig, axes = themed_figure(theme=theme, nrows=2, ncols=2, figsize=(16, 12))

    for ax, (region_name, bounds) in zip(axes.flat, region_list):
        lat_min, lat_max = bounds['lat']
        lon_min, lon_max = bounds['lon']

        mask = np.array([
            (lat_min <= nd.lat <= lat_max and lon_min <= nd.lon <= lon_max)
            for nd in nodes
        ])
        idx = np.where(mask)[0]

        if len(idx) == 0:
            ax.set_title(region_name, fontsize=13, fontweight='bold')
            continue

        # Colour by subregion
        for i in idx:
            r = nodes[i].subregion
            c = REGION_COLORS.get(r, ACCENT_COLORS[0])
            ax.scatter(fd[i], peak_dep[i], color=c, s=40,
                       edgecolors='white' if theme == 'dark' else '#333',
                       linewidths=0.3, alpha=0.7, zorder=3)

        # Label top-3 most depressed
        region_dep = peak_dep[idx]
        top3 = idx[np.argsort(region_dep)[-3:]]
        for gi in top3:
            if peak_dep[gi] > 0.5:
                short = nodes[gi].name.split(',')[0][:18]
                ax.annotate(
                    f'{short}\n{peak_dep[gi]:.1f} psu',
                    (fd[gi], peak_dep[gi]),
                    xytext=(8, 5), textcoords='offset points',
                    fontsize=7, color=text,
                    bbox=dict(boxstyle='round,pad=0.15', facecolor=panel,
                              edgecolor=grid, alpha=0.85),
                )

        # Legend by subregion
        seen = set()
        for i in idx:
            r = nodes[i].subregion
            if r not in seen:
                seen.add(r)
                c = REGION_COLORS.get(r, ACCENT_COLORS[0])
                ax.scatter([], [], color=c, s=30, label=r)
        ax.legend(fontsize=7, loc='upper left', facecolor=panel,
                  edgecolor=grid, labelcolor=text, ncol=2)

        ax.set_xlabel('fjord_depth_norm', fontsize=10)
        ax.set_ylabel('Peak salinity depression (psu)', fontsize=10)
        ax.set_title(f'{region_name}  ({len(idx)} sites)', fontsize=13,
                     fontweight='bold')

    fig.suptitle(
        f'Fjord Depth vs Salinity Depression by Region (fw={fw_strength})',
        fontsize=15, fontweight='bold', color=text, y=1.02,
    )

    if save_path:
        save_figure(fig, save_path)
    return fig
