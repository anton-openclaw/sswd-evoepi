#!/usr/bin/env python3
"""Wildfire visualization — watch SSWD spread through individual sea stars.

Design principles (Kornhauser et al. 2009, JASSS; Paul Tol 2021):
  - Paul Tol "bright" colorblind-safe qualitative palette
  - Dark navy background (marine context, better contrast)
  - Color = disease state (categorical, pre-attentive)
  - Size = body diameter (continuous, ordered)
  - Opacity = resistance score (continuous, ordered)
  - Shape: × for dead agents, ● for alive
  - Persistent legend, scale bar, time indicator
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Dict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np

from sswd_evoepi.snapshots import SnapshotRecorder, IndividualSnapshot
from sswd_evoepi.types import DiseaseState


# ═══════════════════════════════════════════════════════════════════════
# PAUL TOL "BRIGHT" COLORBLIND-SAFE PALETTE
# ═══════════════════════════════════════════════════════════════════════

DS_COLORS = {
    DiseaseState.S:  '#4477AA',   # steel blue — healthy baseline
    DiseaseState.E:  '#CCBB44',   # amber — early warning
    DiseaseState.I1: '#EE7733',   # orange — escalating
    DiseaseState.I2: '#EE6677',   # rose — severe
    DiseaseState.D:  '#888888',   # gray — gone
    DiseaseState.R:  '#66CCEE',   # cyan — recovered
}

DS_LABELS = {
    DiseaseState.S:  'Susceptible',
    DiseaseState.E:  'Exposed',
    DiseaseState.I1: 'Pre-symptomatic',
    DiseaseState.I2: 'Wasting',
    DiseaseState.D:  'Dead',
    DiseaseState.R:  'Recovered',
}

# Background
BG_COLOR = '#0d1b2a'        # dark navy
PANEL_BG = '#1b2838'        # slightly lighter panel
TEXT_COLOR = '#c8d6e5'       # soft white
GRID_COLOR = '#2d3f52'      # subtle grid
BORDER_COLOR = '#3d5a80'    # panel border

# Sizing
MIN_DOT_SIZE = 15            # pts² — settlers still visible
MAX_DOT_SIZE = 120           # pts² — large adults prominent
MAX_BODY_SIZE = 1000.0       # mm — L_inf reference for scaling


# ═══════════════════════════════════════════════════════════════════════
# FRAME RENDERING
# ═══════════════════════════════════════════════════════════════════════

def _render_node(
    ax: plt.Axes,
    snap: Optional[IndividualSnapshot],
    node_name: str,
    sim_day: int,
    habitat_side: float,
):
    """Render one node panel for one timestep."""
    ax.set_facecolor(PANEL_BG)

    if snap is None or snap.n_alive == 0:
        ax.text(0.5, 0.5, "Extinct", ha='center', va='center',
                transform=ax.transAxes, fontsize=16, color='#555555',
                fontweight='bold', style='italic')
        _style_axes(ax, habitat_side, node_name, sim_day, 0, {})
        return

    ds = snap.disease_state
    n_total = snap.n_alive

    # Separate alive (circle) from dead (×)
    alive_mask = ds != DiseaseState.D
    dead_mask = ds == DiseaseState.D

    # --- Alive agents (circles) ---
    if np.any(alive_mask):
        x_a = snap.x[alive_mask]
        y_a = snap.y[alive_mask]
        ds_a = ds[alive_mask]
        sz_a = snap.size[alive_mask]
        r_a = snap.resistance[alive_mask]

        colors_a = np.array([DS_COLORS.get(int(d), '#888888') for d in ds_a])
        sizes_a = MIN_DOT_SIZE + (MAX_DOT_SIZE - MIN_DOT_SIZE) * np.clip(sz_a / MAX_BODY_SIZE, 0, 1)
        # Opacity: base 0.5, scales to 1.0 with resistance
        alphas_a = 0.5 + 0.5 * np.clip(r_a, 0, 1)

        # Draw each alpha group (matplotlib scatter doesn't support per-point alpha easily)
        # Bin into 4 alpha levels for efficiency
        for alpha_lo, alpha_hi, alpha_val in [(0.0, 0.375, 0.35), (0.375, 0.625, 0.55),
                                               (0.625, 0.875, 0.75), (0.875, 1.01, 0.95)]:
            mask = (alphas_a >= alpha_lo) & (alphas_a < alpha_hi)
            if np.any(mask):
                ax.scatter(x_a[mask], y_a[mask], c=colors_a[mask], s=sizes_a[mask],
                           alpha=alpha_val, edgecolors='none', zorder=2, marker='o')

    # --- Dead agents (× markers, small, gray) ---
    if np.any(dead_mask):
        x_d = snap.x[dead_mask]
        y_d = snap.y[dead_mask]
        ax.scatter(x_d, y_d, c=DS_COLORS[DiseaseState.D], s=20,
                   alpha=0.3, marker='x', linewidths=0.8, zorder=1)

    # Compartment counts for title
    counts = {}
    for state in [DiseaseState.S, DiseaseState.E, DiseaseState.I1,
                  DiseaseState.I2, DiseaseState.R]:
        counts[state] = int(np.sum(ds == state))

    _style_axes(ax, habitat_side, node_name, sim_day, n_total, counts)


def _style_axes(ax, habitat_side, node_name, sim_day, n_total, counts):
    """Apply consistent styling to a node panel."""
    ax.set_xlim(-2, habitat_side + 2)
    ax.set_ylim(-2, habitat_side + 2)
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])

    # Border
    for spine in ax.spines.values():
        spine.set_color(BORDER_COLOR)
        spine.set_linewidth(1.5)

    # Title
    year = sim_day // 365
    doy = sim_day % 365
    ax.set_title(f"{node_name}", fontsize=11, color=TEXT_COLOR,
                 fontweight='bold', pad=8)

    # Stats bar at bottom
    if counts:
        s = counts.get(DiseaseState.S, 0)
        e = counts.get(DiseaseState.E, 0)
        i1 = counts.get(DiseaseState.I1, 0)
        i2 = counts.get(DiseaseState.I2, 0)
        r = counts.get(DiseaseState.R, 0)
        stats = f"N={n_total}  S={s}  E={e}  I₁={i1}  I₂={i2}  R={r}"
    else:
        stats = f"N={n_total}"
    ax.text(0.5, -0.04, stats, ha='center', va='top',
            transform=ax.transAxes, fontsize=8, color='#8899aa',
            fontfamily='monospace')

    # Scale bar (bottom-left, 50m)
    bar_len = min(50.0, habitat_side * 0.25)
    bar_y = habitat_side * 0.03
    bar_x = habitat_side * 0.03
    ax.plot([bar_x, bar_x + bar_len], [bar_y, bar_y],
            color='#556677', linewidth=2, solid_capstyle='butt', zorder=5)
    ax.text(bar_x + bar_len / 2, bar_y + habitat_side * 0.02,
            f"{bar_len:.0f}m", ha='center', va='bottom',
            fontsize=6, color='#556677')


def _make_legend(fig):
    """Add a persistent legend to the figure."""
    legend_handles = []
    for state in [DiseaseState.S, DiseaseState.E, DiseaseState.I1,
                  DiseaseState.I2, DiseaseState.R]:
        h = mlines.Line2D([], [], color=DS_COLORS[state], marker='o',
                          markersize=8, linestyle='None', label=DS_LABELS[state],
                          markeredgecolor='none', alpha=0.85)
        legend_handles.append(h)
    # Dead
    h = mlines.Line2D([], [], color=DS_COLORS[DiseaseState.D], marker='x',
                      markersize=7, linestyle='None', label='Dead',
                      markeredgewidth=1.5, alpha=0.5)
    legend_handles.append(h)

    leg = fig.legend(handles=legend_handles, loc='lower center',
                     ncol=6, fontsize=8, frameon=True, fancybox=True,
                     edgecolor=BORDER_COLOR, facecolor=BG_COLOR,
                     labelcolor=TEXT_COLOR, handletextpad=0.3,
                     columnspacing=1.0)
    return leg


def _make_time_label(fig, sim_day):
    """Add time indicator to figure."""
    year = sim_day // 365
    doy = sim_day % 365
    month_names = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                   'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    month = month_names[min(doy // 30, 11)]
    return fig.suptitle(f"Year {year}  ·  {month} (day {doy})",
                        fontsize=13, color=TEXT_COLOR, fontweight='bold',
                        y=0.98)


# ═══════════════════════════════════════════════════════════════════════
# ANIMATION (GIF / MP4)
# ═══════════════════════════════════════════════════════════════════════

def render_wildfire(
    recorder: SnapshotRecorder,
    output_path: str = "wildfire.gif",
    node_id: Optional[int] = None,
    node_names: Optional[List[str]] = None,
    fps: int = 12,
    dpi: int = 100,
    max_frames: int = 500,
    habitat_sides: Optional[Dict[int, float]] = None,
):
    """Render the wildfire animation with dark marine theme.

    Args:
        recorder: SnapshotRecorder with captured data.
        output_path: Output file path (.gif or .mp4).
        node_id: Animate only this node (None = all nodes).
        node_names: Node name list indexed by node_id.
        fps: Frames per second.
        dpi: Image resolution.
        max_frames: Max frames (subsamples if exceeded).
        habitat_sides: Dict of node_id -> habitat side (m).
    """
    import matplotlib.animation as animation

    days = recorder.get_days()
    nodes = recorder.get_nodes() if node_id is None else [node_id]

    if not days or not nodes:
        print("No snapshot data to animate.")
        return

    if len(days) > max_frames:
        step = max(1, len(days) // max_frames)
        days = days[::step]

    if habitat_sides is None:
        habitat_sides = {}

    # Layout
    ncols = min(len(nodes), 5)
    nrows = max(1, (len(nodes) + ncols - 1) // ncols)
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(4.5 * ncols, 4.8 * nrows + 0.8),
                             facecolor=BG_COLOR)
    if isinstance(axes, np.ndarray):
        axes = axes.reshape(nrows, ncols) if axes.ndim == 2 else axes.reshape(1, -1)
    else:
        axes = np.array([[axes]])

    # Hide unused
    for idx in range(len(nodes), nrows * ncols):
        r, c = divmod(idx, ncols)
        axes[r, c].set_visible(False)

    _make_legend(fig)
    title_text = _make_time_label(fig, days[0])

    def update(frame_idx):
        sim_day = days[frame_idx]
        title_text.set_text(
            f"Year {sim_day // 365}  ·  "
            f"{['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'][min(sim_day % 365 // 30, 11)]} "
            f"(day {sim_day % 365})"
        )
        for idx, nid in enumerate(nodes):
            r, c = divmod(idx, ncols)
            ax = axes[r, c]
            ax.clear()
            name = node_names[nid] if node_names and nid < len(node_names) else f"Node {nid}"
            hab = habitat_sides.get(nid, 100.0)
            snap = recorder.get_snapshot(sim_day, nid)
            _render_node(ax, snap, name, sim_day, hab)
        return []

    update(0)
    fig.subplots_adjust(left=0.02, right=0.98, top=0.93, bottom=0.08,
                        wspace=0.08, hspace=0.15)

    anim = animation.FuncAnimation(fig, update, frames=len(days),
                                   interval=1000 // fps, blit=False)

    output_path = str(output_path)
    if output_path.endswith('.mp4'):
        anim.save(output_path, writer='ffmpeg', fps=fps, dpi=dpi,
                  savefig_kwargs={'facecolor': BG_COLOR})
    else:
        anim.save(output_path, writer='pillow', fps=fps, dpi=dpi,
                  savefig_kwargs={'facecolor': BG_COLOR})

    plt.close(fig)
    print(f"Wildfire animation saved: {output_path} ({len(days)} frames, {len(nodes)} nodes)")


# ═══════════════════════════════════════════════════════════════════════
# STATIC SNAPSHOTS
# ═══════════════════════════════════════════════════════════════════════

def render_wildfire_static(
    recorder: SnapshotRecorder,
    output_dir: str,
    node_names: Optional[List[str]] = None,
    key_days: Optional[List[int]] = None,
    habitat_sides: Optional[Dict[int, float]] = None,
    dpi: int = 150,
):
    """Render static snapshots at key moments."""
    days = recorder.get_days()
    nodes = recorder.get_nodes()

    if not days or not nodes:
        return

    if key_days is None:
        indices = np.linspace(0, len(days) - 1, 8, dtype=int)
        key_days = [days[i] for i in indices]

    if habitat_sides is None:
        habitat_sides = {}

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    for sim_day in key_days:
        ncols = min(len(nodes), 5)
        nrows = max(1, (len(nodes) + ncols - 1) // ncols)
        fig, axes = plt.subplots(nrows, ncols,
                                 figsize=(4.5 * ncols, 4.8 * nrows + 0.8),
                                 facecolor=BG_COLOR)
        if isinstance(axes, np.ndarray):
            axes = axes.reshape(nrows, ncols) if axes.ndim == 2 else axes.reshape(1, -1)
        else:
            axes = np.array([[axes]])

        for idx, nid in enumerate(nodes):
            r, c = divmod(idx, ncols)
            ax = axes[r, c]
            name = node_names[nid] if node_names and nid < len(node_names) else f"Node {nid}"
            hab = habitat_sides.get(nid, 100.0)
            snap = recorder.get_snapshot(sim_day, nid)
            _render_node(ax, snap, name, sim_day, hab)

        for idx in range(len(nodes), nrows * ncols):
            r, c = divmod(idx, ncols)
            axes[r, c].set_visible(False)

        _make_legend(fig)
        _make_time_label(fig, sim_day)
        fig.subplots_adjust(left=0.02, right=0.98, top=0.93, bottom=0.08,
                            wspace=0.08, hspace=0.15)

        fig.savefig(out / f"wildfire_day{sim_day:05d}.png", dpi=dpi,
                    facecolor=BG_COLOR, bbox_inches='tight')
        plt.close(fig)

    print(f"Static snapshots saved to {output_dir}/ ({len(key_days)} frames)")
