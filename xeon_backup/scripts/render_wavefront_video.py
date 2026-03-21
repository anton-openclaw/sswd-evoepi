#!/usr/bin/env python3
"""Render wavefront visualization video from monthly snapshot data.

Creates a coastline map showing disease spread as the wavefront moves
north. Each site is a circle:
  - Size: population relative to carrying capacity
  - Color: green (no sick) → red (all sick)

Usage:
    python3 scripts/render_wavefront_video.py \
        results/calibration/W29/monthly_seed42.npz \
        --output results/calibration/W29/wavefront.gif \
        --fps 4
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap

# Try cartopy for proper map projection
try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    HAS_CARTOPY = True
except ImportError:
    HAS_CARTOPY = False
    print("Warning: cartopy not available, using simple lat/lon plot")

# Try geopandas for coastline shapefile
try:
    import geopandas as gpd
    HAS_GEOPANDAS = True
except ImportError:
    HAS_GEOPANDAS = False


def load_snapshots(path: str) -> dict:
    """Load monthly snapshot data from .npz file."""
    data = np.load(path, allow_pickle=True)
    return {
        'sim_days': data['sim_days'],
        'populations': data['populations'],
        'infected': data['infected'],
        'site_lats': data['site_lats'],
        'site_lons': data['site_lons'],
        'site_names': data['site_names'],
        'K': int(data['K']),
        'sst_start_year': int(data['sst_start_year']),
    }


def sim_day_to_date(sim_day: int, start_year: int) -> str:
    """Convert simulation day to approximate calendar date string."""
    year = start_year + sim_day // 365
    month = (sim_day % 365) // 30 + 1
    if month > 12:
        month = 12
    return f"{year}-{month:02d}"


def make_sick_colormap():
    """Green (healthy) → Yellow → Red (all sick)."""
    return LinearSegmentedColormap.from_list(
        'sick_ratio',
        [(0.2, 0.8, 0.2),   # green at 0
         (1.0, 1.0, 0.2),   # yellow at 0.5
         (0.9, 0.1, 0.1)],  # red at 1.0
        N=256,
    )


def render_frames(snap: dict, output_dir: Path, fps: int = 4,
                  min_circle_size: float = 5.0, max_circle_size: float = 120.0):
    """Render all monthly frames and save as individual PNGs + GIF."""
    from PIL import Image

    output_dir.mkdir(parents=True, exist_ok=True)

    lats = snap['site_lats']
    lons = snap['site_lons']
    K = snap['K']
    start_year = snap['sst_start_year']
    populations = snap['populations']
    infected = snap['infected']
    sim_days = snap['sim_days']
    n_frames = len(sim_days)

    cmap = make_sick_colormap()
    norm = mcolors.Normalize(vmin=0.0, vmax=1.0)

    # Map extent: NE Pacific coast (Baja to Aleutians)
    lon_min, lon_max = -180, -112
    lat_min, lat_max = 25, 62

    frame_paths = []

    for fi in range(n_frames):
        day = sim_days[fi]
        pop = populations[fi]
        inf = infected[fi]
        date_str = sim_day_to_date(day, start_year)

        # Compute ratios
        pop_frac = pop.astype(float) / max(K, 1)  # population / K
        sick_ratio = np.where(pop > 0, inf.astype(float) / pop, 0.0)

        # Circle sizes: proportional to pop_frac, with min size for alive sites
        sizes = np.clip(pop_frac, 0, 1) * max_circle_size
        # Sites with pop>0 but tiny get minimum visible size
        alive_mask = pop > 0
        sizes[alive_mask] = np.maximum(sizes[alive_mask], min_circle_size)
        # Dead sites get very small gray dot
        sizes[~alive_mask] = 2.0

        # Colors
        colors = np.zeros((len(lats), 4))
        for i in range(len(lats)):
            if pop[i] > 0:
                colors[i] = cmap(norm(sick_ratio[i]))
            else:
                colors[i] = (0.3, 0.3, 0.3, 0.4)  # gray for extinct

        # Create figure
        if HAS_CARTOPY:
            fig = plt.figure(figsize=(8, 14))
            # Albers Equal Area — good for north-south extent
            proj = ccrs.AlbersEqualArea(
                central_longitude=-135,
                central_latitude=45,
                standard_parallels=(30, 60),
            )
            ax = fig.add_subplot(1, 1, 1, projection=proj)
            ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

            # Minimal background
            ax.add_feature(cfeature.LAND, facecolor='#f0efe8', edgecolor='none')
            ax.add_feature(cfeature.OCEAN, facecolor='#d4e6f1', edgecolor='none')
            ax.add_feature(cfeature.COASTLINE, linewidth=0.5, color='#555555')
            ax.add_feature(cfeature.BORDERS, linewidth=0.3, color='#888888', linestyle='--')

            # Plot sites
            ax.scatter(
                lons, lats,
                s=sizes,
                c=colors,
                transform=ccrs.PlateCarree(),
                zorder=5,
                edgecolors='#333333',
                linewidths=0.2,
                alpha=0.85,
            )
        else:
            fig, ax = plt.subplots(figsize=(8, 14))
            ax.set_xlim(lon_min, lon_max)
            ax.set_ylim(lat_min, lat_max)
            ax.set_facecolor('#d4e6f1')
            ax.scatter(lons, lats, s=sizes, c=colors, zorder=5,
                       edgecolors='#333333', linewidths=0.2, alpha=0.85)

        # Title with date
        ax.set_title(f'SSWD Wavefront — {date_str}', fontsize=14, fontweight='bold', pad=10)

        # Colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, shrink=0.3, pad=0.02, aspect=20)
        cbar.set_label('Infection ratio (sick/total)', fontsize=9)
        cbar.set_ticks([0, 0.25, 0.5, 0.75, 1.0])
        cbar.set_ticklabels(['0% sick', '25%', '50%', '75%', '100% sick'])

        # Size legend
        for pf, label in [(1.0, '100% K'), (0.5, '50% K'), (0.1, '10% K')]:
            ax.scatter([], [], s=pf * max_circle_size, c='gray', alpha=0.5,
                       edgecolors='#333', linewidths=0.3, label=label)
        ax.legend(loc='lower left', title='Pop. Size', fontsize=7, title_fontsize=8,
                  framealpha=0.8, borderpad=0.5)

        # Save frame
        frame_path = output_dir / f'frame_{fi:04d}.png'
        fig.savefig(frame_path, dpi=100, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        frame_paths.append(frame_path)

        if fi % 12 == 0:
            print(f"  Frame {fi+1}/{n_frames} ({date_str})")

    # Compile video — prefer MP4 (ffmpeg) over GIF
    import shutil
    has_ffmpeg = shutil.which('ffmpeg') is not None

    if has_ffmpeg:
        mp4_path = output_dir.parent / (output_dir.stem + '_wavefront.mp4')
        print(f"Compiling {n_frames} frames into MP4 via ffmpeg...")
        import subprocess
        cmd = [
            'ffmpeg', '-y',
            '-framerate', str(fps),
            '-i', str(output_dir / 'frame_%04d.png'),
            '-c:v', 'libx264',
            '-pix_fmt', 'yuv420p',
            '-crf', '23',
            '-preset', 'medium',
            str(mp4_path),
        ]
        subprocess.run(cmd, capture_output=True)
        print(f"  MP4 saved: {mp4_path} ({mp4_path.stat().st_size / 1024 / 1024:.1f} MB)")
        out_path = str(mp4_path)
    else:
        print(f"Compiling {n_frames} frames into GIF (install ffmpeg for MP4)...")
        images = [Image.open(p) for p in frame_paths]
        duration_ms = int(1000 / fps)
        gif_path = output_dir.parent / (output_dir.stem + '_wavefront.gif')
        images[0].save(
            gif_path,
            save_all=True,
            append_images=images[1:],
            duration=duration_ms,
            loop=0,
            optimize=True,
        )
        print(f"  GIF saved: {gif_path} ({gif_path.stat().st_size / 1024 / 1024:.1f} MB)")
        out_path = str(gif_path)

    # Clean up frames
    for p in frame_paths:
        p.unlink()
    if output_dir.exists() and not any(output_dir.iterdir()):
        output_dir.rmdir()

    return out_path


def main():
    parser = argparse.ArgumentParser(description='Render wavefront visualization')
    parser.add_argument('snapshot', help='Path to monthly_seed*.npz file')
    parser.add_argument('--output-dir', default=None,
                        help='Output directory for frames (default: alongside snapshot)')
    parser.add_argument('--fps', type=int, default=4, help='Frames per second in GIF')
    parser.add_argument('--max-size', type=float, default=120,
                        help='Max circle size for full K population')
    args = parser.parse_args()

    snap_path = Path(args.snapshot)
    if not snap_path.exists():
        print(f"Error: {snap_path} not found")
        sys.exit(1)

    print(f"Loading {snap_path}...")
    snap = load_snapshots(str(snap_path))
    print(f"  {len(snap['sim_days'])} frames, {len(snap['site_lats'])} sites, K={snap['K']}")

    output_dir = Path(args.output_dir) if args.output_dir else snap_path.parent / 'frames'
    gif_path = render_frames(snap, output_dir, fps=args.fps, max_circle_size=args.max_size)
    print(f"\nDone! GIF at: {gif_path}")


if __name__ == '__main__':
    main()
