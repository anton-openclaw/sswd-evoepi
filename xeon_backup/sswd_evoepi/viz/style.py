"""Theme styling for SSWD-EvoEpi visualizations.

Provides consistent colors, palettes, and theme-application helpers.
Supports both dark (dashboard) and light (publication) themes.

Usage:
    # Dark theme (default — dashboards, interactive)
    fig, ax = dark_figure()

    # Light/publication theme (papers, reports)
    fig, ax = pub_figure()

    # Any function can accept theme='light' to switch
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

# ═══════════════════════════════════════════════════════════════════════
# DARK COLOR PALETTE
# ═══════════════════════════════════════════════════════════════════════

DARK_BG = '#1a1a2e'
DARK_PANEL = '#16213e'
TEXT_COLOR = '#e0e0e0'
GRID_COLOR = '#2a2a4a'

# ═══════════════════════════════════════════════════════════════════════
# LIGHT / PUBLICATION COLOR PALETTE
# ═══════════════════════════════════════════════════════════════════════

LIGHT_BG = '#ffffff'
LIGHT_PANEL = '#ffffff'
LIGHT_TEXT = '#222222'
LIGHT_GRID = '#cccccc'
LIGHT_SPINE = '#888888'

ACCENT_COLORS = [
    '#e94560',  # crimson
    '#0f3460',  # navy
    '#533483',  # purple
    '#48c9b0',  # teal
    '#f39c12',  # amber
    '#3498db',  # sky blue
    '#e74c3c',  # red
    '#2ecc71',  # green
]

STAGE_COLORS = {
    'settler':  '#f39c12',
    'juvenile': '#48c9b0',
    'subadult': '#3498db',
    'adult':    '#e94560',
}

DEATH_COLORS = {
    'disease':    '#e74c3c',
    'natural':    '#95a5a6',
    'senescence': '#8e44ad',
}

NODE_COLORS = [
    '#e94560', '#48c9b0', '#f39c12', '#3498db', '#2ecc71',
    '#533483', '#e74c3c', '#f1c40f', '#1abc9c', '#9b59b6',
]


# ═══════════════════════════════════════════════════════════════════════
# THEME HELPERS
# ═══════════════════════════════════════════════════════════════════════

def apply_dark_theme(fig=None, ax=None):
    """Apply dark theme to a matplotlib Figure and/or Axes."""
    if fig is not None:
        fig.patch.set_facecolor(DARK_BG)
    if ax is not None:
        ax.set_facecolor(DARK_PANEL)
        ax.tick_params(colors=TEXT_COLOR)
        ax.xaxis.label.set_color(TEXT_COLOR)
        ax.yaxis.label.set_color(TEXT_COLOR)
        ax.title.set_color(TEXT_COLOR)
        for spine in ax.spines.values():
            spine.set_color(GRID_COLOR)
        ax.grid(True, color=GRID_COLOR, alpha=0.3, linewidth=0.5)


def apply_light_theme(fig=None, ax=None):
    """Apply light/publication theme to a matplotlib Figure and/or Axes."""
    if fig is not None:
        fig.patch.set_facecolor(LIGHT_BG)
    if ax is not None:
        ax.set_facecolor(LIGHT_PANEL)
        ax.tick_params(colors=LIGHT_TEXT, labelsize=10)
        ax.xaxis.label.set_color(LIGHT_TEXT)
        ax.yaxis.label.set_color(LIGHT_TEXT)
        ax.title.set_color(LIGHT_TEXT)
        for spine in ax.spines.values():
            spine.set_color(LIGHT_SPINE)
        ax.grid(True, color=LIGHT_GRID, alpha=0.4, linewidth=0.5)


def _apply_theme(fig=None, ax=None, theme='dark'):
    """Apply the specified theme."""
    if theme == 'light':
        apply_light_theme(fig=fig, ax=ax)
    else:
        apply_dark_theme(fig=fig, ax=ax)


def _theme_colors(theme='dark'):
    """Return (bg, panel, text, grid) colours for the given theme."""
    if theme == 'light':
        return LIGHT_BG, LIGHT_PANEL, LIGHT_TEXT, LIGHT_GRID
    return DARK_BG, DARK_PANEL, TEXT_COLOR, GRID_COLOR


def dark_figure(nrows=1, ncols=1, figsize=None, **kwargs):
    """Create a Figure + Axes with the dark theme applied.

    Returns (fig, ax) where ax may be a single Axes or an ndarray.
    """
    if figsize is None:
        figsize = (10, 6) if (nrows == 1 and ncols == 1) else (14, 5 * nrows)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, **kwargs)
    apply_dark_theme(fig=fig)
    if isinstance(axes, np.ndarray):
        for a in axes.flat:
            apply_dark_theme(ax=a)
    else:
        apply_dark_theme(ax=axes)
    return fig, axes


def pub_figure(nrows=1, ncols=1, figsize=None, **kwargs):
    """Create a Figure + Axes with the light/publication theme.

    Returns (fig, ax) where ax may be a single Axes or an ndarray.
    """
    if figsize is None:
        figsize = (10, 6) if (nrows == 1 and ncols == 1) else (14, 5 * nrows)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, **kwargs)
    apply_light_theme(fig=fig)
    if isinstance(axes, np.ndarray):
        for a in axes.flat:
            apply_light_theme(ax=a)
    else:
        apply_light_theme(ax=axes)
    return fig, axes


def themed_figure(nrows=1, ncols=1, figsize=None, theme='light', **kwargs):
    """Create a Figure + Axes with the specified theme.

    Default is 'light' (publication-ready white backgrounds).

    Args:
        theme: 'light' (default) or 'dark'.

    Returns (fig, ax) where ax may be a single Axes or an ndarray.
    """
    if theme == 'light':
        return pub_figure(nrows=nrows, ncols=ncols, figsize=figsize, **kwargs)
    return dark_figure(nrows=nrows, ncols=ncols, figsize=figsize, **kwargs)


def save_figure(fig, save_path, dpi=150):
    """Save a figure with tight layout."""
    fig.tight_layout()
    fig.savefig(save_path, dpi=dpi, facecolor=fig.get_facecolor(),
                edgecolor='none', bbox_inches='tight')
    plt.close(fig)
