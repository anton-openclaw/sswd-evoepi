#!/usr/bin/env python3
"""Generate disease dynamics video from monthly calibration NPZ data.

Coastline map with circles at each site:
- Size: population / K (bigger = more stars)
- Color: green (healthy) → red (sick), based on infected fraction
- Animated monthly over 13 years
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.animation import FuncAnimation, PillowWriter
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from pathlib import Path
import sys

def load_data(npz_path):
    data = np.load(npz_path, allow_pickle=True)
    return {
        'populations': data['populations'],     # (T, N)
        'infected': data['infected'],           # (T, N)
        'site_lats': data['site_lats'],         # (N,)
        'site_lons': data['site_lons'],         # (N,)
        'site_names': data['site_names'],       # (N,)
        'K': int(data['K']),
        'sst_start_year': int(data['sst_start_year']),
    }

def make_video(npz_path, output_path, fps=6, dpi=120):
    print(f"Loading {npz_path}...")
    d = load_data(npz_path)
    
    pops = d['populations']       # (T, N)
    infected = d['infected']      # (T, N)
    lats = d['site_lats']
    lons = d['site_lons']
    K = d['K']
    start_year = d['sst_start_year']
    T, N = pops.shape
    
    print(f"  {T} timesteps, {N} sites, K={K}, start_year={start_year}")
    
    # Compute metrics per frame
    pop_frac = pops / K                                          # population / K
    inf_frac = np.divide(infected, pops, out=np.zeros_like(infected, dtype=float),
                         where=pops > 0)                         # infected / pop
    
    # Size scaling: pop_frac * max_size
    max_marker = 120  # max marker size in points²
    min_marker = 5    # minimum visible size
    
    # Color: green (healthy) → yellow → red (sick)
    cmap = mcolors.LinearSegmentedColormap.from_list('health', 
        ['#2ecc71', '#f1c40f', '#e74c3c'], N=256)
    
    # Map extent (Pacific coast: Alaska to Baja)
    lon_min, lon_max = lons.min() - 3, min(lons.max() + 3, -110)
    lat_min, lat_max = max(lats.min() - 2, 22), min(lats.max() + 2, 62)
    
    # Use Lambert Conformal — great for N America west coast
    proj = ccrs.LambertConformal(
        central_longitude=-135, central_latitude=45,
        standard_parallels=(30, 60)
    )
    
    # Create figure — compact, map fills the frame
    fig = plt.figure(figsize=(7, 10))
    ax = fig.add_axes([0.01, 0.05, 0.85, 0.90], projection=proj)
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    
    # Background — use 50m resolution for better coastline detail
    ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m',
                   facecolor='#f0ebe3', edgecolor='none'))
    ax.add_feature(cfeature.NaturalEarthFeature('physical', 'ocean', '50m',
                   facecolor='#dbe9f4', edgecolor='none'))
    ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '50m',
                   facecolor='none', edgecolor='#888888', linewidth=0.5))
    ax.add_feature(cfeature.BORDERS, linewidth=0.3, color='#cccccc')
    
    # Initial scatter — use RGBA for per-point alpha (extinct sites = transparent)
    sizes = pop_frac[0] * max_marker + min_marker
    colors_raw = inf_frac[0]
    alive_mask = pops[0] > 0
    
    # Convert scalar colors to RGBA via colormap, then set alpha per-point
    norm = plt.Normalize(vmin=0, vmax=0.8)
    rgba = cmap(norm(colors_raw))
    rgba[:, 3] = np.where(alive_mask, 0.85, 0.0)  # transparent if extinct
    
    scatter = ax.scatter(lons, lats, s=sizes, c=rgba,
                        edgecolors='black', linewidths=0.4,
                        transform=ccrs.PlateCarree(), zorder=5)
    
    # Need a ScalarMappable for the colorbar since scatter uses direct RGBA
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    
    # Colorbar — manually positioned tight to right edge (uses ScalarMappable)
    cax = fig.add_axes([0.88, 0.20, 0.025, 0.30])
    cbar = plt.colorbar(sm, cax=cax, label='Infected fraction')
    
    # Title and stats
    title = ax.set_title('', fontsize=14, fontweight='bold', pad=10)
    
    # Stats text box
    stats_text = ax.text(0.02, 0.02, '', transform=ax.transAxes, fontsize=9,
                         verticalalignment='bottom', fontfamily='monospace',
                         bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Size legend
    for frac, label in [(1.0, '100%'), (0.5, '50%'), (0.1, '10%')]:
        ax.scatter([], [], s=frac * max_marker + min_marker, c='gray', alpha=0.5,
                   edgecolors='black', linewidths=0.3, label=f'Pop/K = {label}')
    ax.legend(loc='upper left', fontsize=9, title='Population', title_fontsize=10,
              framealpha=0.9)
    
    def update(frame):
        # Sizes — extinct sites get minimum size (won't show anyway due to alpha=0)
        s = pop_frac[frame] * max_marker + min_marker
        scatter.set_sizes(s)
        
        # Colors — RGBA with per-point alpha for extinction transparency
        c = inf_frac[frame]
        alive = pops[frame] > 0
        frame_rgba = cmap(norm(c))
        frame_rgba[:, 3] = np.where(alive, 0.85, 0.0)  # extinct = fully transparent
        scatter.set_facecolors(frame_rgba)
        
        # Edge colors — also transparent for extinct sites
        edge_rgba = np.zeros((N, 4))
        edge_rgba[:, :3] = 0  # black edges
        edge_rgba[:, 3] = np.where(alive, 1.0, 0.0)
        scatter.set_edgecolors(edge_rgba)
        
        # Time label
        month = frame % 12
        year = frame // 12
        calendar_year = start_year + year
        month_names = ['Jan','Feb','Mar','Apr','May','Jun',
                       'Jul','Aug','Sep','Oct','Nov','Dec']
        title.set_text(f'SSWD Disease Dynamics — {month_names[month]} {calendar_year}')
        
        # Stats
        total_pop = pops[frame].sum()
        total_inf = infected[frame].sum()
        total_K = N * K
        overall_recovery = total_pop / total_K * 100
        overall_inf = total_inf / max(total_pop, 1) * 100
        extinct = (pops[frame] == 0).sum()
        
        stats_text.set_text(
            f'Pop: {total_pop:,.0f} / {total_K:,.0f} ({overall_recovery:.1f}%)\n'
            f'Infected: {overall_inf:.1f}%\n'
            f'Extinct sites: {extinct}/{N}'
        )
        
        return scatter, title, stats_text
    
    print(f"Rendering {T} frames at {fps} fps...")
    anim = FuncAnimation(fig, update, frames=T, interval=1000//fps, blit=False)
    
    # Save as GIF, then auto-crop whitespace from all frames
    output = Path(output_path)
    tmp_output = output.with_stem(output.stem + '_raw')
    
    writer = PillowWriter(fps=fps)
    anim.save(str(tmp_output), writer=writer, dpi=dpi)
    plt.close()
    
    # Auto-crop: find bounding box of non-white content across frames
    from PIL import Image as PILImage, ImageChops
    print("  Auto-cropping whitespace...")
    raw = PILImage.open(str(tmp_output))
    
    # Find crop bbox from first frame
    raw.seek(0)
    first = raw.convert('RGB')
    bg = PILImage.new('RGB', first.size, (255, 255, 255))
    diff = ImageChops.difference(first, bg)
    bbox = diff.getbbox()
    
    # Check a few more frames to get the union bbox
    for f in [raw.n_frames // 4, raw.n_frames // 2, raw.n_frames * 3 // 4]:
        raw.seek(f)
        frame = raw.convert('RGB')
        diff = ImageChops.difference(frame, bg)
        fb = diff.getbbox()
        if fb:
            bbox = (min(bbox[0], fb[0]), min(bbox[1], fb[1]),
                    max(bbox[2], fb[2]), max(bbox[3], fb[3]))
    
    # Add small padding
    pad = 10
    w, h = first.size
    bbox = (max(0, bbox[0] - pad), max(0, bbox[1] - pad),
            min(w, bbox[2] + pad), min(h, bbox[3] + pad))
    
    # Crop all frames
    frames = []
    durations = []
    for i in range(raw.n_frames):
        raw.seek(i)
        cropped = raw.convert('RGBA').crop(bbox)
        frames.append(cropped)
        durations.append(raw.info.get('duration', 1000 // fps))
    
    frames[0].save(str(output), save_all=True, append_images=frames[1:],
                   duration=durations, loop=0, optimize=True)
    
    # Clean up raw
    tmp_output.unlink()
    
    print(f"  Saved: {output} ({output.stat().st_size / 1e6:.1f} MB), "
          f"cropped to {frames[0].size[0]}x{frames[0].size[1]}")
    return str(output)


if __name__ == '__main__':
    npz = sys.argv[1] if len(sys.argv) > 1 else 'results/calibration/W142/monthly_seed42.npz'
    out = sys.argv[2] if len(sys.argv) > 2 else 'results/calibration/W142/disease_dynamics.gif'
    
    import os
    os.chdir(os.path.join(os.path.dirname(__file__), '..'))
    
    make_video(npz, out)
