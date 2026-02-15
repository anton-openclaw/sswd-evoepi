#!/usr/bin/env python3
"""Wildfire visualization — watch the epidemic spread through individuals.

Each frame shows all agents at a node as dots:
  - Color = disease state (green=S, yellow=E, orange=I1, red=I2, black=D, blue=R)
  - Size = body size (larger individuals = bigger dots)
  - Position = x, y within habitat
  
Can show a single node or a grid of all nodes.

Usage:
    from sswd_evoepi.snapshots import SnapshotRecorder
    from viz.wildfire import render_wildfire
    
    render_wildfire(recorder, output_path="wildfire.gif", node_id=2, fps=10)
    render_wildfire(recorder, output_path="wildfire_all.gif")  # all nodes
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
import numpy as np

from sswd_evoepi.snapshots import SnapshotRecorder
from sswd_evoepi.types import DiseaseState


# Disease state → color mapping
DS_COLORS = {
    DiseaseState.S: '#2ecc71',   # green — susceptible
    DiseaseState.E: '#f1c40f',   # yellow — exposed
    DiseaseState.I1: '#e67e22',  # orange — early infectious
    DiseaseState.I2: '#e74c3c',  # red — late infectious / wasting
    DiseaseState.D: '#2c3e50',   # dark gray — dead (shouldn't appear if alive-only)
    DiseaseState.R: '#3498db',   # blue — recovered
}

DS_LABELS = {
    DiseaseState.S: 'Susceptible',
    DiseaseState.E: 'Exposed',
    DiseaseState.I1: 'Early Infectious',
    DiseaseState.I2: 'Late Infectious',
    DiseaseState.D: 'Dead',
    DiseaseState.R: 'Recovered',
}


def _frame_single_node(
    ax: plt.Axes,
    recorder: SnapshotRecorder,
    sim_day: int,
    node_id: int,
    node_name: str = "",
    max_size: float = 1000.0,
    habitat_side: float = 100.0,
):
    """Render one frame for one node."""
    ax.clear()
    snap = recorder.get_snapshot(sim_day, node_id)
    
    if snap is None or snap.n_alive == 0:
        ax.text(0.5, 0.5, "No data", ha='center', va='center',
                transform=ax.transAxes, fontsize=14, color='gray')
        ax.set_xlim(0, habitat_side)
        ax.set_ylim(0, habitat_side)
        ax.set_aspect('equal')
        ax.set_title(f"{node_name} — Day {sim_day} (N=0)")
        return
    
    # Map disease states to colors
    colors = np.array([DS_COLORS.get(int(ds), '#999999') for ds in snap.disease_state])
    
    # Scale dot size by body size (clamped)
    sizes = np.clip(snap.size, 1.0, max_size)
    dot_sizes = 5 + 50 * (sizes / max_size)  # 5-55 point range
    
    ax.scatter(snap.x, snap.y, c=colors, s=dot_sizes, alpha=0.7, edgecolors='none')
    
    # Compute stats
    n_total = snap.n_alive
    ds = snap.disease_state
    n_s = int(np.sum(ds == DiseaseState.S))
    n_e = int(np.sum(ds == DiseaseState.E))
    n_i = int(np.sum((ds == DiseaseState.I1) | (ds == DiseaseState.I2)))
    n_r = int(np.sum(ds == DiseaseState.R))
    
    year = sim_day // 365
    day_of_year = sim_day % 365
    
    title = f"{node_name} — Year {year}, Day {day_of_year}"
    subtitle = f"N={n_total}  S={n_s}  E={n_e}  I={n_i}  R={n_r}"
    ax.set_title(f"{title}\n{subtitle}", fontsize=10)
    
    ax.set_xlim(0, habitat_side)
    ax.set_ylim(0, habitat_side)
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])


def render_wildfire(
    recorder: SnapshotRecorder,
    output_path: str = "wildfire.gif",
    node_id: Optional[int] = None,
    node_names: Optional[List[str]] = None,
    fps: int = 10,
    dpi: int = 100,
    max_frames: int = 500,
    habitat_sides: Optional[dict] = None,
):
    """Render the wildfire animation.
    
    Args:
        recorder: SnapshotRecorder with captured data.
        output_path: Output file path (.gif or .mp4).
        node_id: If specified, animate only this node. Otherwise all nodes in a grid.
        node_names: List of node names indexed by node_id.
        fps: Frames per second.
        dpi: Image resolution.
        max_frames: Maximum number of frames (subsample if exceeded).
        habitat_sides: Dict of node_id -> habitat side length for axis limits.
    """
    import matplotlib.animation as animation
    
    days = recorder.get_days()
    nodes = recorder.get_nodes() if node_id is None else [node_id]
    
    if not days or not nodes:
        print("No snapshot data to animate.")
        return
    
    # Subsample if too many frames
    if len(days) > max_frames:
        step = len(days) // max_frames
        days = days[::step]
    
    # Figure layout
    if len(nodes) == 1:
        fig, axes = plt.subplots(1, 1, figsize=(8, 8))
        axes = np.array([[axes]])
    elif len(nodes) <= 3:
        fig, axes = plt.subplots(1, len(nodes), figsize=(6 * len(nodes), 6))
        axes = axes.reshape(1, -1)
    elif len(nodes) <= 6:
        ncols = 3
        nrows = (len(nodes) + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 6 * nrows))
        axes = axes.reshape(nrows, ncols)
    else:
        ncols = 4
        nrows = (len(nodes) + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 5 * nrows))
        axes = axes.reshape(nrows, ncols)
    
    # Default habitat sides
    if habitat_sides is None:
        habitat_sides = {}
    
    # Legend
    legend_elements = [
        Patch(facecolor=DS_COLORS[DiseaseState.S], label='Susceptible'),
        Patch(facecolor=DS_COLORS[DiseaseState.E], label='Exposed'),
        Patch(facecolor=DS_COLORS[DiseaseState.I1], label='Early Infectious'),
        Patch(facecolor=DS_COLORS[DiseaseState.I2], label='Late Infectious'),
        Patch(facecolor=DS_COLORS[DiseaseState.R], label='Recovered'),
    ]
    
    def update(frame_idx):
        sim_day = days[frame_idx]
        for idx, nid in enumerate(nodes):
            row = idx // axes.shape[1]
            col = idx % axes.shape[1]
            ax = axes[row, col]
            name = node_names[nid] if node_names and nid < len(node_names) else f"Node {nid}"
            hab = habitat_sides.get(nid, 100.0)
            _frame_single_node(ax, recorder, sim_day, nid, name, habitat_side=hab)
        
        # Hide unused subplots
        for idx in range(len(nodes), axes.shape[0] * axes.shape[1]):
            row = idx // axes.shape[1]
            col = idx % axes.shape[1]
            axes[row, col].set_visible(False)
        
        return []
    
    # First frame
    update(0)
    fig.legend(handles=legend_elements, loc='lower center', ncol=5,
               fontsize=9, frameon=True, fancybox=True)
    fig.tight_layout(rect=[0, 0.05, 1, 1])
    
    anim = animation.FuncAnimation(
        fig, update, frames=len(days),
        interval=1000 // fps, blit=False,
    )
    
    output_path = str(output_path)
    if output_path.endswith('.gif'):
        anim.save(output_path, writer='pillow', fps=fps, dpi=dpi)
    elif output_path.endswith('.mp4'):
        anim.save(output_path, writer='ffmpeg', fps=fps, dpi=dpi)
    else:
        anim.save(output_path, fps=fps, dpi=dpi)
    
    plt.close(fig)
    print(f"Wildfire animation saved: {output_path} ({len(days)} frames, {len(nodes)} nodes)")


def render_wildfire_static(
    recorder: SnapshotRecorder,
    output_dir: str,
    node_names: Optional[List[str]] = None,
    key_days: Optional[List[int]] = None,
    habitat_sides: Optional[dict] = None,
    dpi: int = 150,
):
    """Render static snapshots at key moments (for quick inspection).
    
    Args:
        recorder: SnapshotRecorder with data.
        output_dir: Directory for output PNGs.
        node_names: Node name list.
        key_days: Specific days to render. If None, picks 8 evenly spaced.
        habitat_sides: Dict of node_id -> habitat side.
        dpi: Resolution.
    """
    days = recorder.get_days()
    nodes = recorder.get_nodes()
    
    if not days or not nodes:
        return
    
    if key_days is None:
        # Pick 8 evenly spaced frames
        indices = np.linspace(0, len(days) - 1, 8, dtype=int)
        key_days = [days[i] for i in indices]
    
    if habitat_sides is None:
        habitat_sides = {}
    
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    
    for sim_day in key_days:
        ncols = min(len(nodes), 5)
        nrows = (len(nodes) + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 5 * nrows))
        if len(nodes) == 1:
            axes = np.array([[axes]])
        elif nrows == 1:
            axes = axes.reshape(1, -1)
        
        for idx, nid in enumerate(nodes):
            row = idx // ncols
            col = idx % ncols
            ax = axes[row, col]
            name = node_names[nid] if node_names and nid < len(node_names) else f"Node {nid}"
            hab = habitat_sides.get(nid, 100.0)
            _frame_single_node(ax, recorder, sim_day, nid, name, habitat_side=hab)
        
        # Hide unused
        for idx in range(len(nodes), nrows * ncols):
            row = idx // ncols
            col = idx % ncols
            axes[row, col].set_visible(False)
        
        year = sim_day // 365
        doy = sim_day % 365
        fig.suptitle(f"Year {year}, Day {doy} (sim_day={sim_day})", fontsize=14)
        fig.tight_layout()
        fig.savefig(out / f"wildfire_day{sim_day:05d}.png", dpi=dpi, bbox_inches='tight')
        plt.close(fig)
    
    print(f"Static wildfire snapshots saved to {output_dir}/ ({len(key_days)} frames)")
