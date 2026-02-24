#!/usr/bin/env python3
"""One-shot overwater distance matrix computation for SSWD-EvoEpi.

Computes shortest overwater distances between all 489 sites using:
  1. Coastline rasterization → binary land/sea grid
  2. Vectorized sparse graph construction (fixes the bottleneck)
  3. Single-source Dijkstra, one site at a time, saving incrementally

Run once, save forever. Resumes from where it left off if interrupted.

Usage:
    python compute_overwater_distances.py [--resolution 2000] [--resume]
"""

import json
import time
import sys
from pathlib import Path

import numpy as np
import geopandas as gpd
from shapely.geometry import Point
from rasterio import features, Affine
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.csgraph import dijkstra

# ── Config ────────────────────────────────────────────────────────────

LAND_POLYGONS_PATH = Path("data/shorelines/ne_10m_land/ne_10m_land.shp")
SITES_PATH = Path("data/nodes/all_sites.json")
OUTPUT_DIR = Path("results/overwater")
CHECKPOINT_PATH = OUTPUT_DIR / "checkpoint.npz"
FINAL_PATH = OUTPUT_DIR / "distance_matrix.npz"

# Alaska Albers — good equal-area projection for the full NE Pacific
CRS = "EPSG:3338"
DEFAULT_RESOLUTION = 2000  # 2 km — good enough, keeps grid manageable
BUFFER_KM = 100  # buffer around sites

# Clip box for NE Pacific (lat/lon before projection)
CLIP_BOX_LON = (-190, -110)  # Extended past -180 for dateline-crossing Aleutians
CLIP_BOX_LAT = (25, 65)

COASTAL_PENALTY = 1.2
SQRT2 = np.sqrt(2)


# ── Load sites ────────────────────────────────────────────────────────

def load_sites():
    """Load all sites, return (N,2) array of (lat, lon) and names."""
    with open(SITES_PATH) as f:
        sites = json.load(f)
    coords = []
    for s in sites:
        lat = s["latitude"]
        lon = s["longitude"]
        # Normalize positive longitudes (western Aleutians crossing dateline)
        if lon > 0:
            lon = lon - 360.0
        coords.append([lat, lon])
    coords = np.array(coords)
    names = [s["name"] for s in sites]
    regions = [s.get("region", "?") for s in sites]
    return coords, names, regions


# ── Rasterize ─────────────────────────────────────────────────────────

def rasterize_land(land_gdf, bounds, resolution_m):
    """Rasterize land polygons into a binary sea/land grid.
    
    Returns (sea_raster, affine_transform) where sea=True, land=False.
    """
    minx, miny, maxx, maxy = bounds
    width = int(np.ceil((maxx - minx) / resolution_m))
    height = int(np.ceil((maxy - miny) / resolution_m))
    
    transform = Affine.translation(minx, maxy) * Affine.scale(resolution_m, -resolution_m)
    
    land_shapes = [(geom, 1) for geom in land_gdf.geometry if geom is not None]
    land_raster = features.rasterize(
        land_shapes,
        out_shape=(height, width),
        transform=transform,
        fill=0,
        dtype=np.uint8,
    )
    sea_raster = land_raster == 0  # True = sea
    
    print(f"  Grid: {width}×{height} = {width*height:,} cells")
    print(f"  Land: {(land_raster==1).sum():,} ({100*(land_raster==1).mean():.1f}%)")
    print(f"  Sea:  {sea_raster.sum():,} ({100*sea_raster.mean():.1f}%)")
    
    return sea_raster, transform


# ── Build sparse graph (VECTORIZED — fixes the bottleneck) ───────────

def build_sea_graph(sea_raster, resolution_m):
    """Build 8-connected weighted sparse graph of sea cells.
    
    Fully vectorized — no Python loops over cells.
    Returns (csr_matrix, cell_to_idx mapping).
    """
    H, W = sea_raster.shape
    
    # Assign sequential indices to sea cells
    # cell_idx[r, c] = index into sea-cell array, or -1 for land
    cell_idx = np.full((H, W), -1, dtype=np.int32)
    sea_rows, sea_cols = np.where(sea_raster)
    n_sea = len(sea_rows)
    cell_idx[sea_rows, sea_cols] = np.arange(n_sea, dtype=np.int32)
    
    print(f"  Building graph for {n_sea:,} sea cells...")
    
    # Check which sea cells are coastal (adjacent to land) — vectorized
    # Pad with land (False) on edges
    padded = np.pad(sea_raster, 1, constant_values=False)
    # A cell is coastal if it's sea AND any of its 8 neighbors is land
    is_coastal_grid = sea_raster & ~(
        padded[0:H, 0:W] & padded[0:H, 1:W+1] & padded[0:H, 2:W+2] &
        padded[1:H+1, 0:W] &                     padded[1:H+1, 2:W+2] &
        padded[2:H+2, 0:W] & padded[2:H+2, 1:W+1] & padded[2:H+2, 2:W+2]
    )
    is_coastal = is_coastal_grid[sea_rows, sea_cols]  # per sea-cell boolean
    
    # 8 neighbor offsets: (dr, dc, base_cost)
    neighbors = [
        (-1, -1, SQRT2), (-1, 0, 1.0), (-1, 1, SQRT2),
        ( 0, -1, 1.0),                  ( 0, 1, 1.0),
        ( 1, -1, SQRT2), ( 1, 0, 1.0),  ( 1, 1, SQRT2),
    ]
    
    row_list, col_list, cost_list = [], [], []
    
    for dr, dc, base_cost in neighbors:
        # Shifted neighbor coordinates
        nr = sea_rows + dr
        nc = sea_cols + dc
        
        # In-bounds mask
        valid = (nr >= 0) & (nr < H) & (nc >= 0) & (nc < W)
        
        # Neighbor must also be sea
        src_idx = np.arange(n_sea)[valid]
        nbr_idx = cell_idx[nr[valid], nc[valid]]
        sea_mask = nbr_idx >= 0
        
        src_idx = src_idx[sea_mask]
        nbr_idx = nbr_idx[sea_mask]
        
        # Cost: base × coastal penalty if either endpoint is coastal
        cost = np.full(len(src_idx), base_cost * resolution_m / 1000.0)  # km
        coastal_edge = is_coastal[src_idx] | is_coastal[nbr_idx]
        cost[coastal_edge] *= COASTAL_PENALTY
        
        # Map to sea-cell indices
        src_sea = cell_idx[sea_rows[src_idx], sea_cols[src_idx]]
        
        row_list.append(src_sea)
        col_list.append(nbr_idx)
        cost_list.append(cost)
    
    rows = np.concatenate(row_list)
    cols = np.concatenate(col_list)
    costs = np.concatenate(cost_list)
    
    graph = csr_matrix((costs, (rows, cols)), shape=(n_sea, n_sea))
    
    print(f"  Graph: {n_sea:,} nodes, {len(costs):,} edges")
    
    return graph, cell_idx


# ── Map sites to sea cells ────────────────────────────────────────────

def map_sites_to_cells(coords, cell_idx, transform, crs, sea_raster):
    """Map (lat, lon) site coordinates to sea-cell graph indices.
    
    Sites on land get snapped to nearest sea cell.
    Returns array of graph indices (length N), -1 for unmappable.
    """
    H, W = sea_raster.shape
    
    # Project lat/lon → CRS
    gdf = gpd.GeoDataFrame(
        geometry=[Point(lon, lat) for lat, lon in coords],
        crs="EPSG:4326",
    ).to_crs(crs)
    
    inv = ~transform
    result = np.full(len(coords), -1, dtype=np.int32)
    
    for i, geom in enumerate(gdf.geometry):
        col_f, row_f = inv * (geom.x, geom.y)
        c, r = int(round(col_f)), int(round(row_f))
        
        if 0 <= r < H and 0 <= c < W and cell_idx[r, c] >= 0:
            result[i] = cell_idx[r, c]
        else:
            # Snap to nearest sea cell (expanding square search)
            found = False
            for radius in range(1, 50):
                for dr in range(-radius, radius + 1):
                    for dc in range(-radius, radius + 1):
                        if abs(dr) != radius and abs(dc) != radius:
                            continue
                        nr, nc = r + dr, c + dc
                        if 0 <= nr < H and 0 <= nc < W and cell_idx[nr, nc] >= 0:
                            result[i] = cell_idx[nr, nc]
                            found = True
                            break
                    if found:
                        break
                if found:
                    break
            if not found:
                print(f"  WARNING: Could not map site {i} ({coords[i]}) to any sea cell")
    
    n_mapped = (result >= 0).sum()
    print(f"  Mapped {n_mapped}/{len(coords)} sites to sea cells")
    return result


# ── Main computation ──────────────────────────────────────────────────

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Compute overwater distance matrix")
    parser.add_argument("--resolution", type=int, default=DEFAULT_RESOLUTION,
                        help=f"Grid resolution in meters (default: {DEFAULT_RESOLUTION})")
    parser.add_argument("--resume", action="store_true",
                        help="Resume from checkpoint")
    args = parser.parse_args()
    
    resolution_m = args.resolution
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    t0 = time.time()
    
    # ── 1. Load sites ──
    print("1. Loading sites...")
    coords, names, regions = load_sites()
    N = len(coords)
    print(f"  {N} sites loaded")
    
    # ── 2. Load and clip land polygons ──
    print("2. Loading land polygons...")
    from shapely.geometry import box
    land_gdf = gpd.read_file(LAND_POLYGONS_PATH)
    clip = box(CLIP_BOX_LON[0], CLIP_BOX_LAT[0], CLIP_BOX_LON[1], CLIP_BOX_LAT[1])
    land_clipped = gpd.clip(land_gdf, clip)
    print(f"  {len(land_clipped)} land features (clipped to NE Pacific)")
    land_proj = land_clipped.to_crs(CRS)
    
    # ── 3. Compute bounds (land + sites + buffer) ──
    print("3. Computing bounds...")
    site_gdf = gpd.GeoDataFrame(
        geometry=[Point(lon, lat) for lat, lon in coords],
        crs="EPSG:4326",
    ).to_crs(CRS)
    
    all_geoms = list(land_proj.geometry) + list(site_gdf.geometry)
    all_bounds = gpd.GeoSeries(all_geoms, crs=CRS).total_bounds
    buf = BUFFER_KM * 1000
    bounds = (all_bounds[0] - buf, all_bounds[1] - buf,
              all_bounds[2] + buf, all_bounds[3] + buf)
    print(f"  Bounds (m): {bounds[0]:.0f}, {bounds[1]:.0f} → {bounds[2]:.0f}, {bounds[3]:.0f}")
    
    # ── 4. Rasterize ──
    print("4. Rasterizing land polygons...")
    sea_raster, affine = rasterize_land(land_proj, bounds, resolution_m)
    
    # ── 5. Build graph ──
    print("5. Building sea graph...")
    t_graph = time.time()
    graph, cell_idx = build_sea_graph(sea_raster, resolution_m)
    print(f"  Graph built in {time.time() - t_graph:.1f}s")
    
    # ── 6. Map sites ──
    print("6. Mapping sites to graph nodes...")
    site_nodes = map_sites_to_cells(coords, cell_idx, affine, CRS, sea_raster)
    
    unmapped = np.where(site_nodes < 0)[0]
    if len(unmapped) > 0:
        print(f"  WARNING: {len(unmapped)} sites could not be mapped:")
        for idx in unmapped:
            print(f"    {idx}: {names[idx]} ({coords[idx][0]:.2f}, {coords[idx][1]:.2f})")
    
    # ── 7. Compute distances, one source at a time ──
    print(f"7. Computing distances ({N} sources)...")
    
    # Load checkpoint or initialize
    if args.resume and CHECKPOINT_PATH.exists():
        ckpt = np.load(CHECKPOINT_PATH)
        dist_matrix = ckpt["distances"]
        completed = set(ckpt["completed"].tolist())
        print(f"  Resumed from checkpoint: {len(completed)}/{N} sources done")
    else:
        dist_matrix = np.full((N, N), np.inf, dtype=np.float32)
        np.fill_diagonal(dist_matrix, 0.0)
        completed = set()
    
    mapped_mask = site_nodes >= 0
    mapped_indices = np.where(mapped_mask)[0]
    mapped_nodes = site_nodes[mapped_mask]
    
    n_todo = sum(1 for i in mapped_indices if i not in completed)
    print(f"  {n_todo} sources remaining...")
    
    done_count = len(completed)
    t_start = time.time()
    
    for count, src_i in enumerate(mapped_indices):
        if src_i in completed:
            continue
        
        src_node = site_nodes[src_i]
        
        # Single-source Dijkstra → distances to ALL sea cells
        dists = dijkstra(graph, directed=False, indices=src_node,
                         return_predecessors=False)
        
        # Extract distances to other mapped sites
        for dst_i, dst_node in zip(mapped_indices, mapped_nodes):
            d = dists[dst_node]
            if np.isfinite(d):
                dist_matrix[src_i, dst_i] = d
                dist_matrix[dst_i, src_i] = d  # symmetric
        
        completed.add(int(src_i))
        done_count += 1
        
        # Progress report
        elapsed = time.time() - t_start
        rate = (count + 1) / elapsed if elapsed > 0 else 0
        eta = (n_todo - count - 1) / rate if rate > 0 else 0
        
        if (done_count % 10 == 0) or (done_count == N):
            print(f"  [{done_count}/{N}] {names[src_i][:30]:30s} "
                  f"({elapsed:.0f}s elapsed, ~{eta:.0f}s remaining)")
        
        # Checkpoint every 50 sources
        if done_count % 50 == 0:
            np.savez_compressed(CHECKPOINT_PATH,
                                distances=dist_matrix,
                                completed=np.array(list(completed)))
            print(f"  ✓ Checkpoint saved ({done_count} sources)")
    
    # ── 8. Save final result ──
    print("8. Saving final distance matrix...")
    
    finite = dist_matrix[(dist_matrix > 0) & np.isfinite(dist_matrix)]
    
    np.savez_compressed(
        FINAL_PATH,
        distances=dist_matrix,
        coordinates=coords,
        names=np.array(names),
        regions=np.array(regions),
        resolution_m=resolution_m,
        crs=CRS,
    )
    
    # Human-readable summary
    with open(OUTPUT_DIR / "summary.txt", "w") as f:
        f.write(f"Overwater Distance Matrix — SSWD-EvoEpi\n")
        f.write(f"{'='*50}\n\n")
        f.write(f"Sites: {N}\n")
        f.write(f"Resolution: {resolution_m}m\n")
        f.write(f"CRS: {CRS}\n")
        f.write(f"Mapped: {mapped_mask.sum()}/{N}\n")
        
        n_pairs = N * (N - 1) // 2
        upper = dist_matrix[np.triu_indices(N, k=1)]
        n_connected = np.sum(np.isfinite(upper))
        f.write(f"Connected pairs: {n_connected:,}/{n_pairs:,} ({100*n_connected/n_pairs:.1f}%)\n\n")
        
        if len(finite) > 0:
            f.write(f"Distance range: {finite.min():.1f} – {finite.max():.1f} km\n")
            f.write(f"Mean: {finite.mean():.1f} km\n")
            f.write(f"Median: {np.median(finite):.1f} km\n")
    
    total = time.time() - t0
    
    print(f"\n{'='*50}")
    print(f"DONE in {total:.0f}s ({total/60:.1f} min)")
    print(f"Output: {FINAL_PATH}")
    if len(finite) > 0:
        print(f"Distance range: {finite.min():.1f} – {finite.max():.1f} km")
        print(f"Connected pairs: {np.sum(np.isfinite(upper)):,}/{n_pairs:,}")
    print(f"Unmapped sites: {len(unmapped)}")


if __name__ == "__main__":
    main()
