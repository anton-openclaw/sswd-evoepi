#!/usr/bin/env python3
"""
Compute per-site "enclosedness" metrics for the SSWD-EvoEpi model's 896-node coastal network.

Two complementary metrics:
1. Tortuosity Ratio: over-water distance vs. straight-line distance to nearest neighbors
2. Horizon Exposure: fraction of rays that don't hit land within 50km (ray-casting)

Combined into a single enclosedness score for mapping to flushing rates.
"""
import os
import sys
import json
import numpy as np
import pandas as pd
from pathlib import Path
import time
from typing import Dict, List, Tuple, Optional
import warnings

# Suppress some warnings for cleaner output
warnings.filterwarnings('ignore', category=UserWarning)

try:
    import geopandas as gpd
    import shapely
    from shapely.geometry import Point, LineString
    from shapely.ops import unary_union
    HAS_GEOPANDAS = True
except ImportError:
    print("Warning: geopandas/shapely not available. Ray-casting will be skipped.")
    HAS_GEOPANDAS = False

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    from cartopy.feature import NaturalEarthFeature
    HAS_CARTOPY = True
except ImportError:
    print("Warning: cartopy not available. Using alternative land polygon source.")
    HAS_CARTOPY = False


def haversine_distance(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """
    Calculate the great-circle distance between two points on Earth in km.
    """
    from math import radians, cos, sin, asin, sqrt
    
    # Convert to radians
    lat1, lon1, lat2, lon2 = map(radians, [lat1, lon1, lat2, lon2])
    
    # Haversine formula
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a))
    
    # Radius of Earth in km
    r = 6371
    return c * r


def compute_tortuosity_ratio(sites: List[Dict], distance_matrix: np.ndarray, k: int = 20) -> np.ndarray:
    """
    Compute tortuosity ratio for each site using k nearest neighbors.
    
    tortuosity_i = mean(overwater_dist[i,j] / haversine_dist[i,j]) for j in k_nearest
    """
    n_sites = len(sites)
    tortuosity_ratios = np.full(n_sites, np.nan)
    
    print(f"Computing tortuosity ratios for {n_sites} sites...")
    
    for i in range(n_sites):
        site_i = sites[i]
        # Handle both 'lat'/'lon' and 'latitude'/'longitude' fields
        lat_i = site_i.get('lat', site_i.get('latitude'))
        lon_i = site_i.get('lon', site_i.get('longitude'))
        
        if lat_i is None or lon_i is None:
            print(f"Warning: Site {i} missing coordinates: {site_i}")
            continue
        
        # Compute haversine distances to all other sites
        haversine_dists = []
        indices = []
        
        for j in range(n_sites):
            if i == j:
                continue  # Skip self
            site_j = sites[j]
            lat_j = site_j.get('lat', site_j.get('latitude'))
            lon_j = site_j.get('lon', site_j.get('longitude'))
            if lat_j is None or lon_j is None:
                continue
            haversine_dist = haversine_distance(lat_i, lon_i, lat_j, lon_j)
            haversine_dists.append(haversine_dist)
            indices.append(j)
        
        haversine_dists = np.array(haversine_dists)
        indices = np.array(indices)
        
        # Find k nearest neighbors by haversine distance
        k_actual = min(k, len(haversine_dists))
        nearest_idx = np.argsort(haversine_dists)[:k_actual]
        
        # Compute tortuosity ratios for k nearest neighbors
        ratios = []
        for idx in nearest_idx:
            j = indices[idx]
            overwater_dist = distance_matrix[i, j]
            haversine_dist = haversine_dists[idx]
            
            # Handle disconnected pairs (inf/nan in distance matrix)
            if np.isfinite(overwater_dist) and overwater_dist > 0 and haversine_dist > 0:
                ratio = overwater_dist / haversine_dist
                if np.isfinite(ratio) and ratio >= 1.0:  # Sanity check: overwater >= haversine
                    ratios.append(ratio)
        
        # Compute mean ratio
        if len(ratios) > 0:
            tortuosity_ratios[i] = np.mean(ratios)
        
        if (i + 1) % 100 == 0:
            print(f"  Processed {i + 1}/{n_sites} sites")
    
    # Handle any remaining NaNs with the median value
    valid_ratios = tortuosity_ratios[np.isfinite(tortuosity_ratios)]
    if len(valid_ratios) > 0:
        median_ratio = np.median(valid_ratios)
        tortuosity_ratios[np.isnan(tortuosity_ratios)] = median_ratio
    else:
        tortuosity_ratios[:] = 1.5  # Fallback value
    
    print(f"Tortuosity ratio stats: min={np.min(tortuosity_ratios):.2f}, "
          f"mean={np.mean(tortuosity_ratios):.2f}, max={np.max(tortuosity_ratios):.2f}")
    
    return tortuosity_ratios


def load_land_polygons(shapefile_path: str) -> Optional[gpd.GeoDataFrame]:
    """
    Load land polygons from Natural Earth shapefile.
    """
    if not HAS_GEOPANDAS:
        return None
    
    if os.path.exists(shapefile_path):
        print(f"Loading land polygons from {shapefile_path}", flush=True)
        try:
            land_gdf = gpd.read_file(shapefile_path)
            print(f"Loaded {len(land_gdf)} land polygons", flush=True)
            # Ensure we're in WGS84
            if land_gdf.crs != 'EPSG:4326':
                print("Converting to WGS84...", flush=True)
                land_gdf = land_gdf.to_crs('EPSG:4326')
            
            # Create a unified land polygon for faster intersection checking
            print("Creating unified land polygon (this may take a while)...", flush=True)
            unified_land = unary_union(land_gdf.geometry)
            print("Unified land polygon created successfully", flush=True)
            
            return unified_land
        except Exception as e:
            print(f"Error loading shapefile: {e}")
            return None
    else:
        print(f"Shapefile not found at {shapefile_path}")
        return None


def compute_horizon_exposure(sites: List[Dict], land_polygon, max_distance_km: float = 50.0) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute horizon exposure using ray-casting.
    
    For each site, cast 72 rays at 5° intervals outward to max_distance_km.
    horizon_exposure = fraction of rays that DON'T hit land within max_distance_km
    enclosedness = 1 - horizon_exposure
    """
    if not HAS_GEOPANDAS or land_polygon is None:
        print("Skipping ray-casting due to missing dependencies or land data")
        n_sites = len(sites)
        return np.full(n_sites, 0.5), np.full(n_sites, 0.5)  # Default values
    
    n_sites = len(sites)
    horizon_exposures = np.zeros(n_sites)
    enclosedness_rays = np.zeros(n_sites)
    
    # Ray directions: 72 rays at 5° intervals
    n_rays = 72
    ray_angles = np.linspace(0, 360, n_rays, endpoint=False)
    
    print(f"Computing horizon exposure for {n_sites} sites using {n_rays} rays...")
    
    for i, site in enumerate(sites):
        lat = site.get('lat', site.get('latitude'))
        lon = site.get('lon', site.get('longitude'))
        
        if lat is None or lon is None:
            print(f"Warning: Site {i} missing coordinates, setting default enclosedness")
            horizon_exposures[i] = 0.5
            enclosedness_rays[i] = 0.5
            continue
        
        # Count rays that don't hit land
        rays_clear = 0
        
        for angle_deg in ray_angles:
            # Convert angle to radians
            angle_rad = np.radians(angle_deg)
            
            # Approximate endpoint (simple lat/lon offset for 50km)
            # This is a rough approximation; more accurate would use proper geodesic
            km_per_deg_lat = 111.0
            km_per_deg_lon = 111.0 * np.cos(np.radians(lat))
            
            delta_lat = (max_distance_km * np.sin(angle_rad)) / km_per_deg_lat
            delta_lon = (max_distance_km * np.cos(angle_rad)) / km_per_deg_lon
            
            end_lat = lat + delta_lat
            end_lon = lon + delta_lon
            
            # Create ray line
            ray_line = LineString([(lon, lat), (end_lon, end_lat)])
            
            # Check if ray intersects land
            try:
                intersects_land = ray_line.intersects(land_polygon)
                if not intersects_land:
                    rays_clear += 1
            except Exception:
                # If there's an error, assume ray hits land (conservative)
                pass
        
        # Compute horizon exposure
        horizon_exposure = rays_clear / n_rays
        horizon_exposures[i] = horizon_exposure
        enclosedness_rays[i] = 1.0 - horizon_exposure
        
        if (i + 1) % 100 == 0:
            print(f"  Processed {i + 1}/{n_sites} sites")
    
    print(f"Horizon exposure stats: min={np.min(horizon_exposures):.3f}, "
          f"mean={np.mean(horizon_exposures):.3f}, max={np.max(horizon_exposures):.3f}")
    print(f"Ray-cast enclosedness stats: min={np.min(enclosedness_rays):.3f}, "
          f"mean={np.mean(enclosedness_rays):.3f}, max={np.max(enclosedness_rays):.3f}")
    
    return horizon_exposures, enclosedness_rays


def normalize_tortuosity(tortuosity_ratios: np.ndarray) -> np.ndarray:
    """
    Normalize tortuosity ratios to [0,1].
    tortuosity=1 → 0, tortuosity≥3 → 1, linear between.
    """
    normalized = np.clip((tortuosity_ratios - 1.0) / (3.0 - 1.0), 0.0, 1.0)
    return normalized


def combine_metrics(normalized_tortuosity: np.ndarray, enclosedness_rays: np.ndarray, 
                   alpha: float = 0.5) -> np.ndarray:
    """
    Combine tortuosity and ray-cast enclosedness metrics.
    combined = alpha * normalized_tortuosity + (1-alpha) * enclosedness_rays
    """
    combined = alpha * normalized_tortuosity + (1.0 - alpha) * enclosedness_rays
    return np.clip(combined, 0.0, 1.0)


def compute_fjord_depth(sites: List[Dict], distance_matrix: np.ndarray,
                        enclosedness: np.ndarray,
                        open_threshold: float = 0.3) -> np.ndarray:
    """Compute 'fjord depth' — over-water distance to nearest open-coast site.

    Sites classified as open coast (enclosedness < open_threshold) get depth 0.
    All other sites get the minimum over-water distance to any open-coast site.
    This distinguishes fjord mouths (close to open water) from fjord heads
    (far from open water along the waterway).

    Returns array of fjord depths in km, one per site.
    """
    n = len(sites)
    open_mask = enclosedness < open_threshold
    open_indices = np.where(open_mask)[0]
    print(f"  Open-coast sites (enclosedness < {open_threshold}): {len(open_indices)}/{n}")

    fjord_depth = np.zeros(n)
    if len(open_indices) == 0:
        print("  WARNING: No open-coast sites found! Using all zeros.")
        return fjord_depth

    for i in range(n):
        if open_mask[i]:
            fjord_depth[i] = 0.0
        else:
            # Min over-water distance to any open-coast site
            dists = distance_matrix[i, open_indices]
            valid = dists[np.isfinite(dists) & (dists > 0)]
            if len(valid) > 0:
                fjord_depth[i] = np.min(valid)
            else:
                fjord_depth[i] = 0.0  # disconnected — treat as open

    print(f"  Fjord depth stats: min={np.min(fjord_depth):.1f} km, "
          f"median={np.median(fjord_depth):.1f} km, "
          f"max={np.max(fjord_depth):.1f} km, "
          f"mean={np.mean(fjord_depth):.1f} km")
    return fjord_depth


def normalize_fjord_depth(fjord_depth: np.ndarray,
                          max_depth_km: float = 200.0) -> np.ndarray:
    """Normalize fjord depth to [0, 1].
    0 km → 0.0, ≥ max_depth_km → 1.0, linear between.
    """
    return np.clip(fjord_depth / max_depth_km, 0.0, 1.0)


def map_to_flushing_rate(combined_enclosedness: np.ndarray,
                         fjord_depth_norm: np.ndarray = None,
                         phi_open: float = 0.8, phi_fjord: float = 0.03) -> np.ndarray:
    """Map enclosedness + fjord depth to flushing rate.

    If fjord_depth_norm is provided, we use a two-factor model:
        effective_enclosedness = enclosedness * (0.5 + 0.5 * fjord_depth_norm)
    This means:
    - Open coast (low enclosedness): φ ≈ phi_open regardless of depth
    - Enclosed + at mouth (depth~0): φ intermediate (encl * 0.5)
    - Enclosed + deep in fjord (depth~1): φ ≈ phi_fjord (full enclosedness)
    """
    if fjord_depth_norm is not None:
        # Fjord depth amplifies enclosedness: mouth gets 50% effect, head gets 100%
        effective = combined_enclosedness * (0.5 + 0.5 * fjord_depth_norm)
        effective = np.clip(effective, 0.0, 1.0)
    else:
        effective = combined_enclosedness

    flushing_rates = phi_open * (1.0 - effective) + phi_fjord * effective
    return flushing_rates


def main():
    """
    Main function to compute enclosedness metrics and save results.
    """
    # Set up paths
    script_dir = Path(__file__).parent
    repo_root = script_dir.parent
    
    # Data paths
    sites_path = repo_root / "data" / "nodes" / "all_sites.json"
    distance_matrix_path = repo_root / "results" / "overwater" / "distance_matrix.npz"
    land_shapefile_path = repo_root / "data" / "shorelines" / "ne_10m_land" / "ne_10m_land.shp"
    
    # Output paths
    output_dir = repo_root / "data" / "nodes"
    output_json_path = output_dir / "site_enclosedness.json"
    output_csv_path = output_dir / "site_enclosedness.csv"
    
    # Load sites data
    print("Loading sites data...", flush=True)
    if not sites_path.exists():
        print(f"Error: Sites file not found at {sites_path}", flush=True)
        return 1
    
    with open(sites_path, 'r') as f:
        sites = json.load(f)
    
    print(f"Loaded {len(sites)} sites", flush=True)
    
    # Load distance matrix
    print("Loading distance matrix...", flush=True)
    if not distance_matrix_path.exists():
        print(f"Error: Distance matrix not found at {distance_matrix_path}", flush=True)
        return 1
    
    distance_data = np.load(distance_matrix_path)
    distance_matrix = distance_data['distances']
    print(f"Distance matrix shape: {distance_matrix.shape}", flush=True)
    
    if len(sites) != distance_matrix.shape[0]:
        print(f"Error: Site count mismatch. Sites: {len(sites)}, Matrix: {distance_matrix.shape[0]}", flush=True)
        return 1
    
    # Load land polygons
    print("Loading land polygons...", flush=True)
    land_polygon = load_land_polygons(str(land_shapefile_path))
    
    # Compute tortuosity ratios
    print("\n=== Computing Tortuosity Ratios ===")
    tortuosity_ratios = compute_tortuosity_ratio(sites, distance_matrix, k=20)
    
    # Compute horizon exposure
    print("\n=== Computing Horizon Exposure ===")
    horizon_exposures, enclosedness_rays = compute_horizon_exposure(sites, land_polygon)
    
    # Normalize and combine metrics
    print("\n=== Combining Metrics ===")
    normalized_tortuosity = normalize_tortuosity(tortuosity_ratios)
    combined_enclosedness = combine_metrics(normalized_tortuosity, enclosedness_rays)

    # Compute fjord depth — distance to open ocean along waterway
    print("\n=== Computing Fjord Depth ===")
    fjord_depth = compute_fjord_depth(sites, distance_matrix, combined_enclosedness,
                                      open_threshold=0.3)
    fjord_depth_norm = normalize_fjord_depth(fjord_depth, max_depth_km=200.0)
    print(f"Normalized fjord depth stats: min={np.min(fjord_depth_norm):.3f}, "
          f"mean={np.mean(fjord_depth_norm):.3f}, max={np.max(fjord_depth_norm):.3f}")

    flushing_rates = map_to_flushing_rate(combined_enclosedness, fjord_depth_norm)

    print(f"\nCombined enclosedness stats: min={np.min(combined_enclosedness):.3f}, "
          f"mean={np.mean(combined_enclosedness):.3f}, max={np.max(combined_enclosedness):.3f}")
    print(f"Flushing rate stats: min={np.min(flushing_rates):.3f}, "
          f"mean={np.mean(flushing_rates):.3f}, max={np.max(flushing_rates):.3f}")
    
    # Prepare output data
    print("\n=== Saving Results ===")
    output_data = []
    csv_rows = []
    
    for i, site in enumerate(sites):
        lat = site.get('lat', site.get('latitude'))
        lon = site.get('lon', site.get('longitude'))
        
        site_result = {
            "name": site["name"],
            "lat": lat,
            "lon": lon,
            "region": site["region"],
            "tortuosity_ratio": float(tortuosity_ratios[i]),
            "horizon_exposure": float(horizon_exposures[i]),
            "enclosedness_rays": float(enclosedness_rays[i]),
            "enclosedness_combined": float(combined_enclosedness[i]),
            "fjord_depth_km": float(fjord_depth[i]),
            "fjord_depth_norm": float(fjord_depth_norm[i]),
            "flushing_rate": float(flushing_rates[i])
        }
        output_data.append(site_result)
        
        # CSV row
        csv_row = [
            site["name"], lat, lon, site["region"],
            tortuosity_ratios[i], horizon_exposures[i], enclosedness_rays[i],
            combined_enclosedness[i], fjord_depth[i], fjord_depth_norm[i],
            flushing_rates[i]
        ]
        csv_rows.append(csv_row)
    
    # Save JSON
    with open(output_json_path, 'w') as f:
        json.dump(output_data, f, indent=2)
    print(f"Saved JSON results to {output_json_path}")
    
    # Save CSV
    csv_df = pd.DataFrame(csv_rows, columns=[
        'name', 'lat', 'lon', 'region', 'tortuosity_ratio', 'horizon_exposure',
        'enclosedness_rays', 'enclosedness_combined',
        'fjord_depth_km', 'fjord_depth_norm', 'flushing_rate'
    ])
    csv_df.to_csv(output_csv_path, index=False)
    print(f"Saved CSV results to {output_csv_path}")
    
    print("\n=== Summary by Region ===")
    # Group by region and show summary stats
    region_stats = {}
    for i, site in enumerate(sites):
        region = site["region"]
        if region not in region_stats:
            region_stats[region] = []
        region_stats[region].append(combined_enclosedness[i])
    
    for region in sorted(region_stats.keys()):
        values = np.array(region_stats[region])
        print(f"{region:>8}: mean={np.mean(values):.3f} ± {np.std(values):.3f} "
              f"(n={len(values):2d}) [min={np.min(values):.3f}, max={np.max(values):.3f}]")
    
    print("\nEnclosedness computation completed successfully!")
    return 0


if __name__ == "__main__":
    sys.exit(main())