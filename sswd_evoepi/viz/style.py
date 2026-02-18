"""Dark theme styling for SSWD-EvoEpi visualizations.

Provides consistent colors, palettes, and a theme-application helper
so every plot has the same look.
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

# ═══════════════════════════════════════════════════════════════════════
# COLOR PALETTE
# ═══════════════════════════════════════════════════════════════════════

DARK_BG = '#1a1a2e'
DARK_PANEL = '#16213e'
TEXT_COLOR = '#e0e0e0'
GRID_COLOR = '#2a2a4a'

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


def dark_figure(nrows=1, ncols=1, figsize=None, **kwargs):
    """Create a Figure + Axes with the dark theme already applied.

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


def save_figure(fig, save_path, dpi=150):
    """Save a figure with tight layout and dark background."""
    fig.tight_layout()
    fig.savefig(save_path, dpi=dpi, facecolor=fig.get_facecolor(),
                edgecolor='none', bbox_inches='tight')
    plt.close(fig)
