"""Overwater distance calculation module for SSWD-EvoEpi.

Computes shortest overwater distances between coastal sites using:
1. Coastline rasterization to land/sea binary grid
2. 8-connected graph of sea cells
3. Dijkstra pathfinding algorithm
4. Proper handling of CRS projection and antimeridian

Key classes:
  - OverwaterDistanceCalculator: Main computation engine
  - RasterCache: Cached grid storage and retrieval

Key functions:  
  - compute_overwater_distance: Single pair distance
  - compute_distance_matrix: Full N×N matrix via parallel Dijkstra

References:
  - OVERWATER_DISTANCE_DESIGN.md (complete specification)
  - spatial.py compute_distance_matrix() (API compatibility)

Build target: Step 2 of 3 (overwater-2-implement).
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
import time

import numpy as np
import geopandas as gpd
from shapely.geometry import Point, Polygon
from shapely.ops import transform
import rasterio
from rasterio import features, Affine
from rasterio.crs import CRS
from scipy import sparse
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra
from skimage import measure
import pyproj
import warnings

# Suppress rasterio warnings during rasterization
warnings.filterwarnings("ignore", category=UserWarning, module="rasterio")


# ═══════════════════════════════════════════════════════════════════════
# CONSTANTS
# ═══════════════════════════════════════════════════════════════════════

# CRS definitions per design doc
ALASKA_CRS = "EPSG:3338"  # Alaska Albers Equal Area
PACIFIC_COAST_CRS = "+proj=lcc +lat_0=42.5 +lon_0=-125 +lat_1=35 +lat_2=50 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# Default parameters
DEFAULT_RESOLUTION = 500.0  # meters
DEFAULT_MAX_RANGE = 2000.0  # km - maximum distance to compute
ANTIMERIDIAN_SPLIT = 170.0  # degrees west - split longitude for Aleutians
COASTAL_COST_PENALTY = 1.2  # cost multiplier for cells adjacent to land
DEEP_WATER_DISCOUNT = 0.9   # cost multiplier for cells far from shore
DEEP_WATER_THRESHOLD = 20000.0  # meters from shore to be "deep water"

# Progress reporting thresholds
PROGRESS_INTERVAL = 10  # report progress every N nodes for batch computation


# ═══════════════════════════════════════════════════════════════════════
# RASTER CACHE
# ═══════════════════════════════════════════════════════════════════════

class RasterCache:
    """Manages caching of rasterized grids to .npy files."""
    
    def __init__(self, cache_dir: Union[str, Path] = "cache/overwater"):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
    
    def _get_cache_key(self, coastline_path: str, bounds: Tuple[float, ...], 
                      resolution: float, crs: str) -> str:
        """Generate cache key from parameters."""
        coastline_file = Path(coastline_path).name
        bounds_str = "_".join(f"{x:.6f}" for x in bounds)
        crs_str = crs.replace("+", "").replace(" ", "_").replace("=", "")[:20]
        return f"{coastline_file}_{bounds_str}_{resolution:.1f}m_{crs_str}"
    
    def get_cached_raster(self, coastline_path: str, bounds: Tuple[float, ...],
                         resolution: float, crs: str) -> Optional[Tuple[np.ndarray, Affine]]:
        """Load cached raster if available."""
        cache_key = self._get_cache_key(coastline_path, bounds, resolution, crs)
        raster_file = self.cache_dir / f"{cache_key}_raster.npy"
        metadata_file = self.cache_dir / f"{cache_key}_metadata.json"
        
        if raster_file.exists() and metadata_file.exists():
            try:
                raster = np.load(raster_file)
                with open(metadata_file) as f:
                    metadata = json.load(f)
                
                # Reconstruct affine transform
                transform_data = metadata['transform']
                transform = Affine(transform_data[0], transform_data[1], transform_data[2],
                                 transform_data[3], transform_data[4], transform_data[5])
                
                return raster, transform
            except Exception:
                # Cache corrupted, ignore
                pass
        return None
    
    def save_raster(self, raster: np.ndarray, transform: Affine,
                   coastline_path: str, bounds: Tuple[float, ...], 
                   resolution: float, crs: str) -> None:
        """Save raster and metadata to cache."""
        cache_key = self._get_cache_key(coastline_path, bounds, resolution, crs)
        raster_file = self.cache_dir / f"{cache_key}_raster.npy"
        metadata_file = self.cache_dir / f"{cache_key}_metadata.json"
        
        try:
            np.save(raster_file, raster)
            metadata = {
                'transform': list(transform)[:6],
                'shape': raster.shape,
                'resolution': resolution,
                'crs': crs,
                'bounds': bounds,
                'created': time.time()
            }
            with open(metadata_file, 'w') as f:
                json.dump(metadata, f, indent=2)
        except Exception as e:
            print(f"Warning: Failed to cache raster: {e}")


# ═══════════════════════════════════════════════════════════════════════
# MAIN CALCULATOR CLASS
# ═══════════════════════════════════════════════════════════════════════

class OverwaterDistanceCalculator:
    """Compute shortest overwater distances between coastal sites.
    
    Uses raster-based approach: coastline → binary grid → graph → Dijkstra.
    
    Args:
        coastline_path: Path to coastline shapefile/GeoJSON/GPKG
        resolution_m: Grid resolution in meters
        crs: Target coordinate reference system
        adaptive_resolution: Enable adaptive resolution (not yet implemented)
        cache_dir: Directory for caching rasterized grids
    """
    
    def __init__(self, 
                 coastline_path: Union[str, Path],
                 resolution_m: float = DEFAULT_RESOLUTION,
                 crs: Optional[str] = None,
                 adaptive_resolution: bool = False,
                 cache_dir: Optional[Union[str, Path]] = None):
        
        self.coastline_path = Path(coastline_path)
        if not self.coastline_path.exists():
            raise FileNotFoundError(f"Coastline file not found: {coastline_path}")
        
        self.resolution_m = resolution_m
        self.adaptive_resolution = adaptive_resolution
        
        # Initialize cache
        if cache_dir is None:
            cache_dir = self.coastline_path.parent / "cache" / "overwater"
        self.cache = RasterCache(cache_dir)
        
        # Load coastline data
        print(f"Loading coastline from {self.coastline_path}")
        self.coastline_gdf = gpd.read_file(self.coastline_path)
        print(f"Loaded {len(self.coastline_gdf)} coastline features")
        
        # Set CRS - auto-detect best projection
        if crs is None:
            crs = self._select_crs()
        self.crs = crs
        
        # Reproject coastline to working CRS
        self.coastline_projected = self.coastline_gdf.to_crs(self.crs)
        
        # Cached computation products
        self._raster = None
        self._raster_transform = None
        self._sea_graph = None
        self._raster_bounds = None
    
    def _select_crs(self) -> str:
        """Select appropriate CRS based on coastline extent."""
        bounds = self.coastline_gdf.total_bounds
        min_lon, min_lat, max_lon, max_lat = bounds
        
        # Alaska region (north of 55°N)
        if max_lat > 55.0:
            print(f"Using Alaska Albers (EPSG:3338) for northern extent ({max_lat:.1f}°N)")
            return ALASKA_CRS
        else:
            print(f"Using Pacific Coast Lambert for extent ({min_lat:.1f}°N to {max_lat:.1f}°N)")
            return PACIFIC_COAST_CRS
    
    def _get_processing_bounds(self, site_coords: Optional[np.ndarray] = None,
                             buffer_km: float = 50.0) -> Tuple[float, float, float, float]:
        """Get processing bounds, optionally expanded to include sites.
        
        Args:
            site_coords: (N, 2) array of (lat, lon) coordinates
            buffer_km: Buffer distance in km
            
        Returns:
            (minx, miny, maxx, maxy) in projected CRS
        """
        bounds = self.coastline_projected.total_bounds
        minx, miny, maxx, maxy = bounds
        
        if site_coords is not None:
            # Convert site coords to projected CRS and expand bounds
            site_gdf = gpd.GeoDataFrame(
                geometry=[Point(lon, lat) for lat, lon in site_coords],
                crs="EPSG:4326"
            )
            site_proj = site_gdf.to_crs(self.crs)
            site_bounds = site_proj.total_bounds
            
            minx = min(minx, site_bounds[0])
            miny = min(miny, site_bounds[1]) 
            maxx = max(maxx, site_bounds[2])
            maxy = max(maxy, site_bounds[3])
        
        # Add buffer
        buffer_m = buffer_km * 1000.0
        minx -= buffer_m
        miny -= buffer_m
        maxx += buffer_m  
        maxy += buffer_m
        
        return minx, miny, maxx, maxy
    
    def rasterize_coastline(self, bounds: Tuple[float, float, float, float]) -> Tuple[np.ndarray, Affine]:
        """Convert coastline to binary land/sea raster.
        
        Args:
            bounds: (minx, miny, maxx, maxy) in projected CRS
            
        Returns:
            (raster, transform) where raster is uint8 with 0=land, 1=sea
        """
        # Check cache first
        cached = self.cache.get_cached_raster(
            str(self.coastline_path), bounds, self.resolution_m, self.crs
        )
        if cached is not None:
            print("Using cached raster")
            return cached
        
        print(f"Rasterizing coastline at {self.resolution_m}m resolution...")
        
        minx, miny, maxx, maxy = bounds
        
        # Calculate grid dimensions
        width = int(np.ceil((maxx - minx) / self.resolution_m))
        height = int(np.ceil((maxy - miny) / self.resolution_m))
        
        print(f"Grid dimensions: {width} × {height} = {width*height:,} cells ({width*height*4/1e6:.1f} MB)")
        
        # Create affine transform
        transform = Affine.translation(minx, maxy) * Affine.scale(self.resolution_m, -self.resolution_m)
        
        # Rasterize land polygons (land = 1, sea = 0)  
        land_shapes = [(geom, 1) for geom in self.coastline_projected.geometry if geom is not None]
        
        if not land_shapes:
            raise ValueError("No valid coastline geometries found")
        
        land_raster = features.rasterize(
            land_shapes,
            out_shape=(height, width),
            transform=transform,
            fill=0,
            dtype=np.uint8
        )
        
        # Invert: sea = 1, land = 0
        sea_raster = 1 - land_raster
        
        print(f"Rasterization complete: {np.sum(sea_raster):,} sea cells ({100*np.sum(sea_raster)/(width*height):.1f}%)")
        
        # Cache result
        self.cache.save_raster(sea_raster, transform, str(self.coastline_path), bounds, self.resolution_m, self.crs)
        
        return sea_raster, transform
    
    def build_sea_graph(self, raster: np.ndarray, add_costs: bool = True) -> csr_matrix:
        """Build 8-connected adjacency graph of sea cells.
        
        Args:
            raster: Binary raster with 1=sea, 0=land
            add_costs: Add coastal/deep-water cost modifiers
            
        Returns:
            Sparse CSR adjacency matrix
        """
        print("Building sea cell graph...")
        
        height, width = raster.shape
        n_cells = height * width
        
        # Find sea cells
        sea_mask = (raster == 1)
        sea_indices = np.where(sea_mask)
        n_sea_cells = len(sea_indices[0])
        
        print(f"Sea cells: {n_sea_cells:,} / {n_cells:,} ({100*n_sea_cells/n_cells:.1f}%)")
        
        if n_sea_cells == 0:
            raise ValueError("No sea cells found - check coastline data")
        
        # Create mapping from (row, col) to sea-cell index
        cell_to_index = np.full((height, width), -1, dtype=np.int32)
        for i, (row, col) in enumerate(zip(sea_indices[0], sea_indices[1])):
            cell_to_index[row, col] = i
        
        # 8-connected neighbors (includes diagonals)
        neighbors = [
            (-1, -1), (-1, 0), (-1, 1),
            (0, -1),           (0, 1),
            (1, -1),  (1, 0),  (1, 1)
        ]
        
        # Distance costs for neighbors
        neighbor_costs = [
            np.sqrt(2), 1.0, np.sqrt(2),  # diagonal, straight, diagonal
            1.0,             1.0,         # straight, straight  
            np.sqrt(2), 1.0, np.sqrt(2)  # diagonal, straight, diagonal
        ]
        
        # Build adjacency lists
        rows, cols, data = [], [], []
        
        for i, (row, col) in enumerate(zip(sea_indices[0], sea_indices[1])):
            for (dr, dc), base_cost in zip(neighbors, neighbor_costs):
                nr, nc = row + dr, col + dc
                
                # Check bounds
                if 0 <= nr < height and 0 <= nc < width:
                    neighbor_idx = cell_to_index[nr, nc]
                    if neighbor_idx >= 0:  # neighbor is also sea
                        cost = base_cost
                        
                        if add_costs:
                            # Coastal penalty: check if current or neighbor cell is adjacent to land
                            is_coastal = self._is_coastal_cell(row, col, raster) or self._is_coastal_cell(nr, nc, raster)
                            if is_coastal:
                                cost *= COASTAL_COST_PENALTY
                        
                        rows.append(i)
                        cols.append(neighbor_idx)
                        data.append(cost)
        
        # Create sparse matrix
        graph = csr_matrix((data, (rows, cols)), shape=(n_sea_cells, n_sea_cells))
        
        print(f"Graph built: {n_sea_cells:,} nodes, {len(data):,} edges")
        
        return graph
    
    def _is_coastal_cell(self, row: int, col: int, raster: np.ndarray) -> bool:
        """Check if a sea cell is adjacent to land."""
        height, width = raster.shape
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue
                nr, nc = row + dr, col + dc
                if 0 <= nr < height and 0 <= nc < width:
                    if raster[nr, nc] == 0:  # land
                        return True
        return False
    
    def map_sites_to_cells(self, site_coords: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Map lat/lon coordinates to raster cell indices.
        
        Args:
            site_coords: (N, 2) array of (lat, lon) coordinates
            
        Returns:
            Tuple of (valid_indices, cell_indices) where:
            - valid_indices: boolean mask of sites that map to sea cells
            - cell_indices: sea-cell indices for valid sites
        """
        if self._raster is None or self._raster_transform is None:
            raise ValueError("Must call rasterize_coastline() first")
        
        print(f"Mapping {len(site_coords)} sites to grid cells...")
        
        # Convert coordinates to projected CRS
        site_gdf = gpd.GeoDataFrame(
            geometry=[Point(lon, lat) for lat, lon in site_coords],
            crs="EPSG:4326"
        )
        site_proj = site_gdf.to_crs(self.crs)
        
        # Get projected coordinates
        proj_coords = np.array([[geom.x, geom.y] for geom in site_proj.geometry])
        
        # Convert to raster indices
        inv_transform = ~self._raster_transform
        raster_coords = np.array([inv_transform * (x, y) for x, y in proj_coords])
        
        # Round to integer indices
        col_indices = np.round(raster_coords[:, 0]).astype(int)
        row_indices = np.round(raster_coords[:, 1]).astype(int)
        
        # Check bounds and sea cells
        height, width = self._raster.shape
        valid_mask = (
            (0 <= row_indices) & (row_indices < height) &
            (0 <= col_indices) & (col_indices < width)
        )
        
        sea_mask = np.zeros(len(site_coords), dtype=bool)
        sea_mask[valid_mask] = self._raster[row_indices[valid_mask], col_indices[valid_mask]] == 1
        
        # Create mapping from grid coordinates to sea-cell indices
        sea_indices = np.where(self._raster == 1)
        cell_to_index = np.full((height, width), -1, dtype=np.int32)
        for i, (row, col) in enumerate(zip(sea_indices[0], sea_indices[1])):
            cell_to_index[row, col] = i
        
        # Get sea-cell indices for valid sites
        cell_indices = np.full(len(site_coords), -1, dtype=np.int32)
        for i, (row, col) in enumerate(zip(row_indices, col_indices)):
            if sea_mask[i]:
                cell_indices[i] = cell_to_index[row, col]
        
        n_valid = np.sum(sea_mask)
        n_land = np.sum(valid_mask) - n_valid
        n_oob = len(site_coords) - np.sum(valid_mask)
        
        print(f"Site mapping: {n_valid} → sea cells, {n_land} → land (will snap), {n_oob} out of bounds")
        
        # For sites on land, find nearest sea cell
        if n_land > 0:
            print("Snapping land-based sites to nearest sea cells...")
            for i in range(len(site_coords)):
                if valid_mask[i] and not sea_mask[i]:
                    # Find nearest sea cell within reasonable distance
                    row, col = row_indices[i], col_indices[i]
                    nearest_sea_idx = self._find_nearest_sea_cell(row, col, cell_to_index)
                    if nearest_sea_idx >= 0:
                        cell_indices[i] = nearest_sea_idx
                        sea_mask[i] = True
        
        n_final = np.sum(sea_mask)
        print(f"Final mapping: {n_final} sites successfully mapped to sea cells")
        
        return sea_mask, cell_indices[sea_mask]
    
    def _find_nearest_sea_cell(self, row: int, col: int, cell_to_index: np.ndarray,
                             max_search_radius: int = 20) -> int:
        """Find nearest sea cell using expanding search."""
        height, width = cell_to_index.shape
        
        for radius in range(1, max_search_radius + 1):
            for dr in range(-radius, radius + 1):
                for dc in range(-radius, radius + 1):
                    if abs(dr) == radius or abs(dc) == radius:  # only check perimeter
                        nr, nc = row + dr, col + dc
                        if 0 <= nr < height and 0 <= nc < width:
                            cell_idx = cell_to_index[nr, nc]
                            if cell_idx >= 0:
                                return cell_idx
        return -1
    
    def compute_all_pairs_distances(self, site_coords: np.ndarray,
                                   progress_callback: Optional[callable] = None) -> np.ndarray:
        """Compute shortest overwater distances between all site pairs.
        
        Args:
            site_coords: (N, 2) array of (lat, lon) coordinates
            progress_callback: Optional function called with (current, total)
            
        Returns:
            (N, N) distance matrix in kilometers
        """
        print(f"Computing overwater distances for {len(site_coords)} sites...")
        
        # Set up raster and graph if not already done
        bounds = self._get_processing_bounds(site_coords)
        if self._raster is None:
            self._raster, self._raster_transform = self.rasterize_coastline(bounds)
            self._raster_bounds = bounds
        
        if self._sea_graph is None:
            self._sea_graph = self.build_sea_graph(self._raster)
        
        # Map sites to cells
        site_mask, cell_indices = self.map_sites_to_cells(site_coords)
        
        if len(cell_indices) == 0:
            raise ValueError("No sites could be mapped to sea cells")
        
        n_valid_sites = len(cell_indices)
        
        print(f"Running Dijkstra for {n_valid_sites} source nodes...")
        start_time = time.time()
        
        # Compute all-pairs shortest paths using Dijkstra
        # This returns distances in terms of cells
        try:
            distance_matrix_cells = dijkstra(
                self._sea_graph, 
                directed=False,
                indices=cell_indices,
                return_predecessors=False
            )
        except Exception as e:
            raise RuntimeError(f"Dijkstra computation failed: {e}")
        
        # Convert cell distances to metric distances
        distance_matrix_km = distance_matrix_cells * (self.resolution_m / 1000.0)
        
        # Create full N×N matrix for all input sites
        n_sites = len(site_coords)
        full_matrix = np.full((n_sites, n_sites), np.inf, dtype=np.float64)
        
        # Fill valid entries
        valid_indices = np.where(site_mask)[0]
        for i, vi in enumerate(valid_indices):
            for j, vj in enumerate(valid_indices):
                full_matrix[vi, vj] = distance_matrix_km[i, j]
        
        # Set diagonal to zero
        np.fill_diagonal(full_matrix, 0.0)
        
        elapsed = time.time() - start_time
        print(f"Distance computation complete in {elapsed:.1f}s")
        print(f"Distance range: {np.nanmin(distance_matrix_km[distance_matrix_km > 0]):.1f} - {np.nanmax(distance_matrix_km):.1f} km")
        
        return full_matrix
    
    def compute_single_distance(self, coord1: Tuple[float, float], 
                              coord2: Tuple[float, float]) -> float:
        """Compute overwater distance between two coordinates.
        
        Args:
            coord1: (lat, lon) of first point
            coord2: (lat, lon) of second point
            
        Returns:
            Distance in kilometers (np.inf if no path exists)
        """
        site_coords = np.array([coord1, coord2])
        matrix = self.compute_all_pairs_distances(site_coords)
        return matrix[0, 1]
    
    @classmethod
    def from_sites_file(cls, sites_file: Union[str, Path], 
                       coastline_path: Union[str, Path],
                       **kwargs) -> 'OverwaterDistanceCalculator':
        """Factory method to create calculator from site definitions file.
        
        Args:
            sites_file: CSV/YAML file containing site definitions
            coastline_path: Path to coastline data
            **kwargs: Additional arguments passed to constructor
            
        Returns:
            Configured OverwaterDistanceCalculator
        """
        sites_path = Path(sites_file)
        
        if sites_path.suffix.lower() == '.csv':
            import pandas as pd
            df = pd.read_csv(sites_path)
            if 'lat' not in df.columns or 'lon' not in df.columns:
                raise ValueError("Sites CSV must contain 'lat' and 'lon' columns")
        elif sites_path.suffix.lower() in ['.yaml', '.yml']:
            # Assume SSWD-EvoEpi format from spatial.py
            from .spatial import load_node_definitions_yaml
            node_defs = load_node_definitions_yaml(sites_path)
            # Extract coordinates
            import pandas as pd
            df = pd.DataFrame([{
                'lat': nd.lat,
                'lon': nd.lon,
                'name': nd.name
            } for nd in node_defs])
        else:
            raise ValueError("Sites file must be CSV or YAML")
        
        return cls(coastline_path, **kwargs)


# ═══════════════════════════════════════════════════════════════════════
# UTILITY FUNCTIONS  
# ═══════════════════════════════════════════════════════════════════════

def validate_distance_matrix(distances: np.ndarray, 
                            site_coords: np.ndarray,
                            tolerance: float = 0.1) -> Dict[str, float]:
    """Validate distance matrix against Haversine baselines.
    
    Args:
        distances: (N, N) overwater distance matrix (km)
        site_coords: (N, 2) array of (lat, lon) coordinates  
        tolerance: Maximum fractional difference to report as error
        
    Returns:
        Dictionary of validation statistics
    """
    from .spatial import haversine_km
    
    n_sites = len(site_coords)
    n_pairs = n_sites * (n_sites - 1) // 2
    
    haversine_diffs = []
    overwater_valid = 0
    
    for i in range(n_sites):
        for j in range(i + 1, n_sites):
            lat1, lon1 = site_coords[i]
            lat2, lon2 = site_coords[j]
            
            haversine_dist = haversine_km(lat1, lon1, lat2, lon2)
            overwater_dist = distances[i, j]
            
            if np.isfinite(overwater_dist):
                overwater_valid += 1
                if haversine_dist > 0:
                    ratio = overwater_dist / haversine_dist
                    haversine_diffs.append(ratio)
    
    haversine_diffs = np.array(haversine_diffs)
    
    return {
        'n_pairs': n_pairs,
        'overwater_valid': overwater_valid,
        'coverage': overwater_valid / n_pairs,
        'mean_ratio': np.mean(haversine_diffs) if len(haversine_diffs) > 0 else np.nan,
        'median_ratio': np.median(haversine_diffs) if len(haversine_diffs) > 0 else np.nan,
        'min_ratio': np.min(haversine_diffs) if len(haversine_diffs) > 0 else np.nan,
        'max_ratio': np.max(haversine_diffs) if len(haversine_diffs) > 0 else np.nan,
        'std_ratio': np.std(haversine_diffs) if len(haversine_diffs) > 0 else np.nan
    }


def export_distance_matrix(distances: np.ndarray,
                          site_coords: np.ndarray, 
                          site_ids: Optional[np.ndarray] = None,
                          output_path: Union[str, Path] = "overwater_distances.npz",
                          metadata: Optional[Dict] = None) -> None:
    """Export distance matrix in compressed NumPy format.
    
    Args:
        distances: (N, N) distance matrix (km)
        site_coords: (N, 2) array of (lat, lon) coordinates
        site_ids: Optional (N,) array of site identifiers  
        output_path: Output file path (.npz)
        metadata: Optional metadata dictionary
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    if site_ids is None:
        site_ids = np.arange(len(site_coords))
    
    if metadata is None:
        metadata = {}
    
    metadata.update({
        'created': time.time(),
        'n_sites': len(site_coords),
        'format_version': '1.0'
    })
    
    np.savez_compressed(
        output_path,
        distances=distances.astype(np.float32),
        site_ids=site_ids,
        coordinates=site_coords.astype(np.float64),
        metadata=np.array([metadata], dtype=object)[0]  # ensure object array
    )
    
    print(f"Distance matrix exported to {output_path}")


# ═══════════════════════════════════════════════════════════════════════
# INTEGRATION WITH SSWD-EVOEPI
# ═══════════════════════════════════════════════════════════════════════

def compute_distance_matrix_overwater(
    lats: np.ndarray,
    lons: np.ndarray, 
    coastline_path: Optional[Union[str, Path]] = None,
    resolution_m: float = DEFAULT_RESOLUTION,
    use_cache: bool = True,
    **kwargs
) -> np.ndarray:
    """Drop-in replacement for spatial.compute_distance_matrix() using overwater distances.
    
    Args:
        lats: (N,) node latitudes
        lons: (N,) node longitudes  
        coastline_path: Path to coastline data (required)
        resolution_m: Raster resolution in meters
        use_cache: Enable raster caching
        **kwargs: Additional arguments (ignored for compatibility)
        
    Returns:
        (N, N) waterway distance matrix (km)
    """
    if coastline_path is None:
        raise ValueError("coastline_path is required for overwater distance computation")
    
    site_coords = np.column_stack([lats, lons])
    
    calc = OverwaterDistanceCalculator(
        coastline_path=coastline_path,
        resolution_m=resolution_m
    )
    
    return calc.compute_all_pairs_distances(site_coords)


# ═══════════════════════════════════════════════════════════════════════
# COMMAND LINE INTERFACE
# ═══════════════════════════════════════════════════════════════════════

def main():
    """Command line interface for overwater distance computation."""
    import argparse
    import json
    
    parser = argparse.ArgumentParser(
        description="Compute overwater distances between coastal sites",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Compute distances for 5-node test network
  python -m sswd_evoepi.overwater_distance \\
    --coastline data/shorelines/ne_pacific_coastline.geojson \\
    --sites-json '[{"lat": 57.06, "lon": -135.34}, {"lat": 49.52, "lon": -123.25}]' \\
    --output results/

  # Use sites from file  
  python -m sswd_evoepi.overwater_distance \\
    --coastline data/shorelines/ne_pacific_coastline.geojson \\
    --sites-file data/nodes/5node_network.yaml \\
    --resolution 250 --output results/
        """
    )
    
    parser.add_argument("--coastline", required=True,
                       help="Path to coastline shapefile/GeoJSON/GPKG")
    parser.add_argument("--sites-json", 
                       help="JSON array of site coordinates: [{'lat': ..., 'lon': ...}, ...]")
    parser.add_argument("--sites-file",
                       help="CSV/YAML file with site definitions")
    parser.add_argument("--output", default="results/",
                       help="Output directory")
    parser.add_argument("--resolution", type=float, default=DEFAULT_RESOLUTION,
                       help="Grid resolution in meters")
    parser.add_argument("--crs", 
                       help="Coordinate reference system (auto-detect if not specified)")
    parser.add_argument("--no-cache", action="store_true",
                       help="Disable raster caching")
    parser.add_argument("--validate", action="store_true",
                       help="Validate results against Haversine distances")
    
    args = parser.parse_args()
    
    # Load site coordinates
    if args.sites_json:
        site_data = json.loads(args.sites_json)
        site_coords = np.array([[s['lat'], s['lon']] for s in site_data])
        site_names = [s.get('name', f'Site_{i}') for i, s in enumerate(site_data)]
    elif args.sites_file:
        sites_path = Path(args.sites_file)
        if sites_path.suffix.lower() == '.csv':
            import pandas as pd
            df = pd.read_csv(sites_path)
            if 'lat' not in df.columns or 'lon' not in df.columns:
                raise ValueError("Sites CSV must contain 'lat' and 'lon' columns")
            site_coords = df[['lat', 'lon']].values
            site_names = df.get('name', [f'Site_{i}' for i in range(len(df))]).tolist()
        elif sites_path.suffix.lower() in ['.yaml', '.yml']:
            # Import here to avoid circular imports
            from .spatial import load_node_definitions_yaml
            node_defs = load_node_definitions_yaml(sites_path)
            site_coords = np.array([[nd.lat, nd.lon] for nd in node_defs])
            site_names = [nd.name for nd in node_defs]
        else:
            raise ValueError("Sites file must be CSV or YAML")
    else:
        raise ValueError("Must specify either --sites-json or --sites-file")
    
    print(f"Loaded {len(site_coords)} sites")
    for i, (name, (lat, lon)) in enumerate(zip(site_names, site_coords)):
        print(f"  {i}: {name} ({lat:.3f}°N, {lon:.3f}°E)")
    
    # Set up calculator
    calc = OverwaterDistanceCalculator(
        coastline_path=args.coastline,
        resolution_m=args.resolution,
        crs=args.crs
    )
    
    # Compute distances  
    print("\nComputing overwater distances...")
    distances = calc.compute_all_pairs_distances(site_coords)
    
    # Export results
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    matrix_file = output_dir / "overwater_distances.npz"
    export_distance_matrix(
        distances, site_coords, 
        site_ids=np.arange(len(site_coords)),
        output_path=matrix_file,
        metadata={
            'coastline_file': str(args.coastline),
            'resolution_m': args.resolution,
            'crs': calc.crs,
            'site_names': site_names
        }
    )
    
    # Save human-readable summary
    summary_file = output_dir / "distance_summary.txt" 
    with open(summary_file, 'w') as f:
        f.write("Overwater Distance Matrix Summary\n")
        f.write("=" * 35 + "\n\n")
        f.write(f"Sites: {len(site_coords)}\n")
        f.write(f"Coastline: {args.coastline}\n")
        f.write(f"Resolution: {args.resolution}m\n")
        f.write(f"CRS: {calc.crs}\n\n")
        
        # Distance statistics
        finite_distances = distances[np.isfinite(distances) & (distances > 0)]
        if len(finite_distances) > 0:
            f.write(f"Distance range: {np.min(finite_distances):.1f} - {np.max(finite_distances):.1f} km\n")
            f.write(f"Mean distance: {np.mean(finite_distances):.1f} km\n")
            f.write(f"Median distance: {np.median(finite_distances):.1f} km\n")
        
        n_pairs = len(site_coords) * (len(site_coords) - 1) // 2
        n_connected = np.sum(np.isfinite(distances[np.triu_indices(len(site_coords), k=1)]))
        f.write(f"\nConnected pairs: {n_connected} / {n_pairs} ({100*n_connected/n_pairs:.1f}%)\n")
        
        f.write("\nSite Pairs:\n")
        for i in range(len(site_coords)):
            for j in range(i + 1, len(site_coords)):
                dist = distances[i, j]
                if np.isfinite(dist):
                    f.write(f"  {site_names[i]} ↔ {site_names[j]}: {dist:.1f} km\n")
                else:
                    f.write(f"  {site_names[i]} ↔ {site_names[j]}: no path\n")
    
    print(f"Results saved to {output_dir}")
    
    # Validation
    if args.validate:
        print("\nValidating against Haversine distances...")
        stats = validate_distance_matrix(distances, site_coords)
        
        print(f"  Coverage: {stats['coverage']:.1%}")
        if not np.isnan(stats['mean_ratio']):
            print(f"  Overwater/Haversine ratio: {stats['mean_ratio']:.2f} ± {stats['std_ratio']:.2f}")
            print(f"  Range: {stats['min_ratio']:.2f} - {stats['max_ratio']:.2f}")


if __name__ == "__main__":
    main()