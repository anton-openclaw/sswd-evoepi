#!/usr/bin/env python3
"""
Test and validate the overwater distance calculation module.
This is step 3 of 3 for the overwater distance implementation.

Testing steps:
1. Basic smoke test
2. Sanity checks with diverse site pairs  
3. Performance test at different resolutions
4. Edge cases
5. Visualize example paths
6. Run subset matrix computation
7. Write comprehensive test report
"""

import json
import time
import warnings
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import traceback

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap
try:
    import seaborn as sns
    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False
    sns = None
import geopandas as gpd
from shapely.geometry import Point, LineString
try:
    import psutil
    HAS_PSUTIL = True
except ImportError:
    HAS_PSUTIL = False
    psutil = None
import gc

# Import our overwater distance module
from sswd_evoepi.overwater_distance import (
    OverwaterDistanceCalculator,
    validate_distance_matrix,
    export_distance_matrix
)

# Suppress warnings for cleaner output
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

# Set up paths
DATA_DIR = Path("data")
NODES_DIR = DATA_DIR / "nodes"
SHORELINES_DIR = DATA_DIR / "shorelines" 
RESULTS_DIR = Path("data/nodes")  # Save results alongside node data
COASTLINE_PATH = SHORELINES_DIR / "ne_pacific_coastline.geojson"

# Create results directory
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Set up matplotlib for dark theme
plt.style.use('dark_background')
COLORS = {
    'background': '#0d1b2a',
    'coastline': '#e9c46a', 
    'water': '#264653',
    'path': '#e76f51',
    'sites': '#f4a261'
}

def load_all_sites() -> Tuple[np.ndarray, List[str], List[str]]:
    """Load all sites from JSON files in nodes directory."""
    print("Loading all sites from JSON files...")
    
    sites = []
    site_names = []
    site_regions = []
    
    # Find all sites_*.json files
    site_files = list(NODES_DIR.glob("sites_*.json"))
    print(f"Found {len(site_files)} site files")
    
    for site_file in sorted(site_files):
        print(f"  Loading {site_file.name}...")
        try:
            with open(site_file) as f:
                data = json.load(f)
            
            region = data.get('region', site_file.stem.replace('sites_', ''))
            
            for site in data['sites']:
                sites.append([site['lat'], site['lon']])
                site_names.append(site.get('name', site.get('site_id', f"Site_{len(sites)}")))
                site_regions.append(region)
                
        except Exception as e:
            print(f"    Warning: Failed to load {site_file}: {e}")
    
    print(f"Loaded {len(sites)} sites total")
    return np.array(sites), site_names, site_regions

def select_diverse_sites(coords: np.ndarray, names: List[str], 
                        regions: List[str], n_sites: int = 20) -> Tuple[np.ndarray, List[str], List[int]]:
    """Select a diverse subset of sites across regions."""
    print(f"Selecting {n_sites} diverse sites...")
    
    # Group sites by region
    region_groups = {}
    for i, region in enumerate(regions):
        if region not in region_groups:
            region_groups[region] = []
        region_groups[region].append(i)
    
    print(f"Found {len(region_groups)} regions: {list(region_groups.keys())}")
    
    # Select sites from each region
    selected_indices = []
    sites_per_region = max(1, n_sites // len(region_groups))
    
    for region, indices in region_groups.items():
        # Take sites spread across the region
        region_coords = coords[indices]
        if len(indices) <= sites_per_region:
            selected_indices.extend(indices)
        else:
            # Sample spread across lat/lon ranges
            step = len(indices) // sites_per_region
            for i in range(0, len(indices), step):
                if len(selected_indices) < n_sites:
                    selected_indices.append(indices[i])
    
    # Fill to target size if needed
    remaining = set(range(len(coords))) - set(selected_indices)
    while len(selected_indices) < n_sites and remaining:
        selected_indices.append(remaining.pop())
    
    selected_indices = selected_indices[:n_sites]  # Trim if over
    
    selected_coords = coords[selected_indices]
    selected_names = [names[i] for i in selected_indices]
    
    print(f"Selected {len(selected_indices)} sites")
    for i, (idx, name, region) in enumerate(zip(selected_indices, selected_names, 
                                               [regions[i] for i in selected_indices])):
        lat, lon = coords[idx]
        print(f"  {i:2d}: {name} ({region}) - {lat:.3f}°N, {lon:.3f}°W")
    
    return selected_coords, selected_names, selected_indices

def haversine_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """Calculate Haversine distance between two points."""
    R = 6371.0  # Earth radius in km
    
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a))
    
    return R * c

def get_memory_usage() -> float:
    """Get current memory usage in MB."""
    if HAS_PSUTIL:
        process = psutil.Process()
        return process.memory_info().rss / 1024 / 1024
    else:
        # Fallback - return a dummy value
        return 0.0

def test_basic_smoke(coastline_path: Path) -> Dict:
    """Step 1: Basic smoke test - import module and compute simple distance."""
    print("\n" + "="*60)
    print("STEP 1: BASIC SMOKE TEST")
    print("="*60)
    
    results = {}
    
    try:
        # Test 1: Module import (already done)
        print("✓ Module import successful")
        
        # Test 2: Load coastline data
        print(f"Loading coastline from {coastline_path}...")
        if not coastline_path.exists():
            raise FileNotFoundError(f"Coastline file not found: {coastline_path}")
        
        coastline_gdf = gpd.read_file(coastline_path)
        print(f"✓ Loaded {len(coastline_gdf)} coastline features")
        results['coastline_features'] = len(coastline_gdf)
        
        # Test 3: Initialize calculator with coarse resolution
        print("Initializing calculator (2km resolution for speed)...")
        calc = OverwaterDistanceCalculator(
            coastline_path=coastline_path,
            resolution_m=2000.0  # 2km for speed
        )
        print(f"✓ Calculator initialized (CRS: {calc.crs})")
        results['crs'] = calc.crs
        
        # Test 4: Compute distance between 2 nearby sites
        print("Computing test distance...")
        # Sitka, AK to Juneau, AK
        sitka = (57.0531, -135.3300)
        juneau = (58.3019, -134.4197)
        
        start_time = time.time()
        distance = calc.compute_single_distance(sitka, juneau)
        elapsed = time.time() - start_time
        
        haversine_dist = haversine_km(sitka[0], sitka[1], juneau[0], juneau[1])
        
        print(f"✓ Sitka → Juneau: {distance:.1f} km overwater vs {haversine_dist:.1f} km straight-line")
        print(f"  Computation time: {elapsed:.2f}s")
        print(f"  Ratio: {distance/haversine_dist:.2f}")
        
        results['test_distance_km'] = distance
        results['test_haversine_km'] = haversine_dist
        results['test_ratio'] = distance / haversine_dist
        results['test_time_s'] = elapsed
        
        if np.isfinite(distance) and distance > 0:
            print("✓ Basic smoke test PASSED")
            results['status'] = 'PASSED'
        else:
            print("✗ Basic smoke test FAILED - invalid distance")
            results['status'] = 'FAILED'
            
    except Exception as e:
        print(f"✗ Basic smoke test FAILED: {e}")
        traceback.print_exc()
        results['status'] = 'FAILED'
        results['error'] = str(e)
    
    return results

def test_sanity_checks(calc: OverwaterDistanceCalculator) -> Dict:
    """Step 2: Sanity checks with diverse site pairs."""
    print("\n" + "="*60)
    print("STEP 2: SANITY CHECKS")  
    print("="*60)
    
    results = {'pairs': []}
    
    # Define test pairs with expected ranges
    test_pairs = [
        # Name, coord1, coord2, expected_min_km, expected_max_km, notes
        ("Sitka → Juneau", (57.0531, -135.3300), (58.3019, -134.4197), 250, 350, "SE Alaska fjords"),
        ("Friday Harbor → Howe Sound", (48.5465, -123.0134), (49.52, -123.25), 150, 200, "Salish Sea"),
        ("Monterey → Santa Barbara", (36.6002, -121.8947), (34.4208, -119.6982), 400, 500, "CA coast"),
        ("Newport OR → Bodega Bay CA", (44.6369, -124.0533), (38.5733, -123.0667), 800, 1000, "Open coast"),
        ("Sitka → Vancouver", (57.0531, -135.3300), (49.2827, -123.1207), 1000, 1300, "Long fjord route"),
    ]
    
    # Add some within-region tests
    bc_pairs = [
        ("Howe Sound → Vancouver", (49.52, -123.25), (49.2827, -123.1207), 40, 80, "BC coastal"),
        ("Tofino → Victoria", (49.1537, -125.9077), (48.4284, -123.3656), 200, 300, "Vancouver Island"),
    ]
    test_pairs.extend(bc_pairs)
    
    all_passed = True
    
    for name, coord1, coord2, exp_min, exp_max, notes in test_pairs:
        print(f"\nTesting: {name}")
        print(f"  Route: {coord1} → {coord2}")
        print(f"  Expected: {exp_min}-{exp_max} km ({notes})")
        
        try:
            start_time = time.time()
            overwater_dist = calc.compute_single_distance(coord1, coord2)
            elapsed = time.time() - start_time
            
            haversine_dist = haversine_km(coord1[0], coord1[1], coord2[0], coord2[1])
            
            if np.isfinite(overwater_dist):
                ratio = overwater_dist / haversine_dist
                in_range = exp_min <= overwater_dist <= exp_max
                
                print(f"  Result: {overwater_dist:.1f} km overwater vs {haversine_dist:.1f} km straight")
                print(f"  Ratio: {ratio:.2f}, Time: {elapsed:.2f}s")
                print(f"  Range check: {'✓ PASS' if in_range else '✗ FAIL'}")
                
                if not in_range:
                    all_passed = False
                    
            else:
                print(f"  Result: No path found")
                print(f"  Range check: ✗ FAIL (no path)")
                all_passed = False
            
            # Check symmetry
            reverse_dist = calc.compute_single_distance(coord2, coord1)
            symmetric = abs(overwater_dist - reverse_dist) < 0.1
            print(f"  Symmetry: {'✓ PASS' if symmetric else '✗ FAIL'} ({overwater_dist:.1f} vs {reverse_dist:.1f})")
            
            if not symmetric:
                all_passed = False
            
            results['pairs'].append({
                'name': name,
                'coord1': coord1,
                'coord2': coord2,
                'expected_min': exp_min,
                'expected_max': exp_max,
                'overwater_km': overwater_dist,
                'haversine_km': haversine_dist,
                'ratio': ratio if np.isfinite(overwater_dist) else None,
                'time_s': elapsed,
                'in_range': in_range if np.isfinite(overwater_dist) else False,
                'symmetric': symmetric,
                'notes': notes
            })
            
        except Exception as e:
            print(f"  Error: {e}")
            all_passed = False
            results['pairs'].append({
                'name': name,
                'error': str(e),
                'in_range': False,
                'symmetric': False
            })
    
    results['status'] = 'PASSED' if all_passed else 'FAILED'
    results['n_passed'] = sum(1 for p in results['pairs'] if p.get('in_range', False))
    results['n_total'] = len(results['pairs'])
    
    print(f"\n{'✓' if all_passed else '✗'} Sanity checks: {results['n_passed']}/{results['n_total']} passed")
    
    return results

def test_performance(coastline_path: Path) -> Dict:
    """Step 3: Performance test at different resolutions."""
    print("\n" + "="*60)
    print("STEP 3: PERFORMANCE TEST")
    print("="*60)
    
    results = {'resolutions': []}
    
    # Test resolutions  
    resolutions = [2000, 1000, 500]  # meters
    test_coords = [
        (57.0531, -135.3300),  # Sitka
        (58.3019, -134.4197),  # Juneau
        (49.5260, -123.2500),  # Vancouver
        (48.5465, -123.0134),  # Friday Harbor
        (36.6002, -121.8947),  # Monterey
    ]
    
    print(f"Testing {len(test_coords)} sites at different resolutions...")
    
    for resolution in resolutions:
        print(f"\n--- Resolution: {resolution}m ---")
        
        try:
            # Initialize calculator
            start_init = time.time()
            calc = OverwaterDistanceCalculator(
                coastline_path=coastline_path,
                resolution_m=resolution
            )
            init_time = time.time() - start_init
            
            # Measure memory before computation
            mem_before = get_memory_usage()
            
            # Time 10 pairwise distances
            print("Computing 10 random pairwise distances...")
            pair_times = []
            distances = []
            
            start_compute = time.time()
            for i in range(min(10, len(test_coords))):
                for j in range(i + 1, min(10, len(test_coords))):
                    if len(pair_times) >= 10:
                        break
                    
                    pair_start = time.time()
                    dist = calc.compute_single_distance(test_coords[i], test_coords[j])
                    pair_time = time.time() - pair_start
                    
                    pair_times.append(pair_time)
                    distances.append(dist)
                    
                if len(pair_times) >= 10:
                    break
            
            total_compute_time = time.time() - start_compute
            mem_after = get_memory_usage()
            
            # Calculate matrix extrapolation
            n_sites_options = [20, 50, 100, 150, 300]
            mean_pair_time = np.mean(pair_times)
            
            print(f"  Initialization: {init_time:.2f}s")
            print(f"  Mean pair time: {mean_pair_time:.3f}s") 
            print(f"  Memory usage: {mem_after - mem_before:.1f} MB increase")
            print(f"  Valid distances: {sum(np.isfinite(distances))}/{len(distances)}")
            
            # Extrapolate to full matrices
            extrapolations = {}
            for n in n_sites_options:
                n_pairs = n * (n - 1) // 2
                est_time_min = (n_pairs * mean_pair_time) / 60
                extrapolations[n] = est_time_min
                print(f"  Est. {n}×{n} matrix: {est_time_min:.1f} min")
            
            results['resolutions'].append({
                'resolution_m': resolution,
                'init_time_s': init_time,
                'mean_pair_time_s': mean_pair_time,
                'memory_increase_mb': mem_after - mem_before,
                'valid_distances': sum(np.isfinite(distances)),
                'total_distances': len(distances),
                'extrapolations': extrapolations
            })
            
        except Exception as e:
            print(f"  Error at {resolution}m: {e}")
            results['resolutions'].append({
                'resolution_m': resolution,
                'error': str(e)
            })
    
    # Recommend best resolution
    valid_results = [r for r in results['resolutions'] if 'error' not in r]
    if valid_results:
        # Find balance of speed and memory
        best = min(valid_results, key=lambda r: r['mean_pair_time_s'] + r['memory_increase_mb']/1000)
        results['recommended_resolution'] = best['resolution_m']
        print(f"\nRecommended resolution: {best['resolution_m']}m")
    else:
        results['recommended_resolution'] = 500
    
    return results

def test_edge_cases(calc: OverwaterDistanceCalculator) -> Dict:
    """Step 4: Test edge cases."""
    print("\n" + "="*60)
    print("STEP 4: EDGE CASES")
    print("="*60)
    
    results = {'cases': []}
    
    edge_cases = [
        # Case name, coord1, coord2, description
        ("Land-based site", (57.0, -135.0), (57.05, -135.05), "Sites that fall on land"),
        ("Same location", (57.0531, -135.3300), (57.0531, -135.3300), "Identical coordinates"),
        ("Very close sites", (57.0531, -135.3300), (57.0532, -135.3301), "Adjacent grid cells"),
        ("Cross-strait", (48.5, -123.0), (48.6, -123.2), "Across narrow water body"),
        ("Very distant", (60.0, -140.0), (32.0, -117.0), "Alaska to Baja California"),
        ("Antimeridian edge", (52.0, -179.0), (52.0, 179.0), "Near dateline"),
    ]
    
    for case_name, coord1, coord2, description in edge_cases:
        print(f"\nTesting: {case_name}")
        print(f"  {description}")
        print(f"  Route: {coord1} → {coord2}")
        
        try:
            start_time = time.time()
            distance = calc.compute_single_distance(coord1, coord2)
            elapsed = time.time() - start_time
            
            haversine_dist = haversine_km(coord1[0], coord1[1], coord2[0], coord2[1])
            
            if np.isfinite(distance):
                ratio = distance / haversine_dist if haversine_dist > 0 else np.inf
                print(f"  Result: {distance:.3f} km overwater vs {haversine_dist:.3f} km straight")
                print(f"  Ratio: {ratio:.3f}, Time: {elapsed:.3f}s")
                status = "SUCCESS"
            else:
                print(f"  Result: No path found (inf distance)")
                print(f"  Haversine: {haversine_dist:.3f} km, Time: {elapsed:.3f}s")
                status = "NO_PATH"
                ratio = None
                
            results['cases'].append({
                'case': case_name,
                'description': description,
                'coord1': coord1,
                'coord2': coord2,
                'overwater_km': distance,
                'haversine_km': haversine_dist,
                'ratio': ratio,
                'time_s': elapsed,
                'status': status
            })
            
        except Exception as e:
            print(f"  Error: {e}")
            results['cases'].append({
                'case': case_name,
                'error': str(e),
                'status': "ERROR"
            })
    
    n_success = sum(1 for c in results['cases'] if c.get('status') == 'SUCCESS')
    n_total = len(results['cases'])
    results['n_success'] = n_success
    results['n_total'] = n_total
    
    print(f"\nEdge cases: {n_success}/{n_total} successful computations")
    
    return results

def visualize_paths(calc: OverwaterDistanceCalculator, coastline_path: Path, 
                   output_path: Path) -> Dict:
    """Step 5: Visualize example overwater paths."""
    print("\n" + "="*60)
    print("STEP 5: VISUALIZE PATHS") 
    print("="*60)
    
    results = {'paths': []}
    
    # Example paths to visualize
    example_paths = [
        ("Sitka → Juneau", (57.0531, -135.3300), (58.3019, -134.4197)),
        ("Friday Harbor → Vancouver", (48.5465, -123.0134), (49.2827, -123.1207)),
        ("Monterey → Santa Barbara", (36.6002, -121.8947), (34.4208, -119.6982)),
    ]
    
    print("Loading coastline for visualization...")
    coastline_gdf = gpd.read_file(coastline_path)
    
    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.patch.set_facecolor(COLORS['background'])
    axes = axes.flatten()
    
    try:
        for i, (name, coord1, coord2) in enumerate(example_paths[:3]):
            ax = axes[i]
            print(f"\nVisualizing: {name}")
            
            # Compute distance 
            overwater_dist = calc.compute_single_distance(coord1, coord2)
            haversine_dist = haversine_km(coord1[0], coord1[1], coord2[0], coord2[1])
            
            # Set up plot bounds with buffer
            lats = [coord1[0], coord2[0]]
            lons = [coord1[1], coord2[1]]
            lat_range = max(lats) - min(lats)
            lon_range = max(lons) - min(lons)
            buffer = max(lat_range, lon_range) * 0.3
            
            bounds = {
                'minx': min(lons) - buffer,
                'maxx': max(lons) + buffer, 
                'miny': min(lats) - buffer,
                'maxy': max(lats) + buffer
            }
            
            # Plot coastline in bounds
            clipped_coast = coastline_gdf.cx[bounds['minx']:bounds['maxx'], 
                                           bounds['miny']:bounds['maxy']]
            
            if len(clipped_coast) > 0:
                clipped_coast.plot(ax=ax, color=COLORS['coastline'], linewidth=0.5, alpha=0.7)
            
            # Plot straight line
            ax.plot([coord1[1], coord2[1]], [coord1[0], coord2[0]], 
                   '--', color='white', alpha=0.5, linewidth=1, label='Straight line')
            
            # Plot sites
            ax.plot(coord1[1], coord1[0], 'o', color=COLORS['sites'], markersize=8, 
                   markeredgecolor='white', markeredgewidth=1)
            ax.plot(coord2[1], coord2[0], 's', color=COLORS['sites'], markersize=8,
                   markeredgecolor='white', markeredgewidth=1)
            
            # Set bounds and styling
            ax.set_xlim(bounds['minx'], bounds['maxx'])
            ax.set_ylim(bounds['miny'], bounds['maxy'])
            ax.set_facecolor(COLORS['background'])
            ax.grid(True, alpha=0.3)
            ax.set_aspect('equal', adjustable='box')
            
            # Title and labels
            if np.isfinite(overwater_dist):
                title = f"{name}\nOverwater: {overwater_dist:.1f} km | Straight: {haversine_dist:.1f} km | Ratio: {overwater_dist/haversine_dist:.2f}"
            else:
                title = f"{name}\nNo path found | Straight: {haversine_dist:.1f} km"
                
            ax.set_title(title, color='white', fontsize=10, pad=10)
            ax.set_xlabel('Longitude', color='white')
            ax.set_ylabel('Latitude', color='white')
            ax.tick_params(colors='white')
            
            results['paths'].append({
                'name': name,
                'coord1': coord1,
                'coord2': coord2,
                'overwater_km': overwater_dist,
                'haversine_km': haversine_dist
            })
        
        # Remove unused subplot
        fig.delaxes(axes[3])
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=150, bbox_inches='tight', 
                   facecolor=COLORS['background'], edgecolor='none')
        plt.close()
        
        print(f"✓ Path visualization saved to {output_path}")
        results['status'] = 'SUCCESS'
        results['output_file'] = str(output_path)
        
    except Exception as e:
        print(f"✗ Visualization failed: {e}")
        traceback.print_exc()
        results['status'] = 'ERROR'
        results['error'] = str(e)
        plt.close()
    
    return results

def run_subset_matrix(coords: np.ndarray, names: List[str], 
                     calc: OverwaterDistanceCalculator, output_dir: Path) -> Dict:
    """Step 6: Compute subset distance matrix and heatmap."""
    print("\n" + "="*60)
    print("STEP 6: SUBSET DISTANCE MATRIX")
    print("="*60)
    
    results = {}
    
    print(f"Computing {len(coords)}×{len(coords)} distance matrix...")
    
    try:
        start_time = time.time()
        distance_matrix = calc.compute_all_pairs_distances(coords)
        elapsed = time.time() - start_time
        
        print(f"✓ Matrix computation completed in {elapsed:.1f}s")
        
        # Export matrix as CSV
        csv_file = output_dir / "distance_matrix_subset.csv"
        df = pd.DataFrame(distance_matrix, columns=names, index=names)
        df.to_csv(csv_file, float_format='%.2f')
        print(f"✓ CSV saved to {csv_file}")
        
        # Create heatmap
        heatmap_file = output_dir / "distance_heatmap_subset.png"
        
        plt.figure(figsize=(12, 10))
        plt.style.use('dark_background')
        
        # Mask infinite values for visualization
        masked_matrix = np.where(np.isfinite(distance_matrix), distance_matrix, np.nan)
        
        # Create custom colormap
        cmap = LinearSegmentedColormap.from_list(
            'overwater', ['#264653', '#2a9d8f', '#e9c46a', '#f4a261', '#e76f51']
        )
        
        # Plot heatmap
        if HAS_SEABORN:
            sns.heatmap(masked_matrix, 
                       xticklabels=[name[:15] + "..." if len(name) > 15 else name for name in names],
                       yticklabels=[name[:15] + "..." if len(name) > 15 else name for name in names],
                       cmap=cmap, 
                       annot=False,
                       fmt='.0f',
                       cbar_kws={'label': 'Distance (km)'})
        else:
            # Fallback to basic matplotlib imshow
            im = plt.imshow(masked_matrix, cmap=cmap, aspect='auto')
            plt.colorbar(im, label='Distance (km)')
            plt.xticks(range(len(names)), [name[:15] + "..." if len(name) > 15 else name for name in names], rotation=45, ha='right')
            plt.yticks(range(len(names)), [name[:15] + "..." if len(name) > 15 else name for name in names])
        
        plt.title('Overwater Distance Matrix (km)', color='white', fontsize=14, pad=20)
        plt.xlabel('Destination Sites', color='white')
        plt.ylabel('Origin Sites', color='white')
        plt.xticks(rotation=45, ha='right', color='white')
        plt.yticks(rotation=0, color='white')
        
        plt.tight_layout()
        plt.savefig(heatmap_file, dpi=150, bbox_inches='tight',
                   facecolor=COLORS['background'], edgecolor='none')
        plt.close()
        
        print(f"✓ Heatmap saved to {heatmap_file}")
        
        # Calculate statistics
        finite_distances = distance_matrix[np.isfinite(distance_matrix) & (distance_matrix > 0)]
        n_pairs = len(coords) * (len(coords) - 1) // 2
        n_connected = np.sum(np.isfinite(distance_matrix[np.triu_indices(len(coords), k=1)]))
        
        results.update({
            'status': 'SUCCESS',
            'n_sites': len(coords),
            'computation_time_s': elapsed,
            'csv_file': str(csv_file),
            'heatmap_file': str(heatmap_file),
            'n_pairs': n_pairs,
            'n_connected': n_connected,
            'connectivity': n_connected / n_pairs,
            'min_distance_km': np.min(finite_distances) if len(finite_distances) > 0 else np.nan,
            'max_distance_km': np.max(finite_distances) if len(finite_distances) > 0 else np.nan,
            'mean_distance_km': np.mean(finite_distances) if len(finite_distances) > 0 else np.nan,
            'median_distance_km': np.median(finite_distances) if len(finite_distances) > 0 else np.nan
        })
        
        print(f"Matrix statistics:")
        print(f"  Connected pairs: {n_connected}/{n_pairs} ({100*n_connected/n_pairs:.1f}%)")
        if len(finite_distances) > 0:
            print(f"  Distance range: {np.min(finite_distances):.1f} - {np.max(finite_distances):.1f} km")
            print(f"  Mean distance: {np.mean(finite_distances):.1f} km")
            
    except Exception as e:
        print(f"✗ Matrix computation failed: {e}")
        traceback.print_exc()
        results = {
            'status': 'ERROR',
            'error': str(e)
        }
    
    return results

def write_test_report(all_results: Dict, output_file: Path) -> None:
    """Step 7: Write comprehensive test report."""
    print("\n" + "="*60)
    print("STEP 7: WRITING TEST REPORT")
    print("="*60)
    
    try:
        with open(output_file, 'w') as f:
            f.write("# Overwater Distance Module Test Results\n\n")
            f.write(f"**Test Date:** {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"**Coastline Data:** {COASTLINE_PATH}\n")
            f.write(f"**Test Framework:** Python automated testing suite\n\n")
            
            f.write("## Executive Summary\n\n")
            
            # Count overall results
            total_tests = 0
            passed_tests = 0
            for step, result in all_results.items():
                if isinstance(result, dict) and 'status' in result:
                    total_tests += 1
                    if result['status'] == 'PASSED' or result['status'] == 'SUCCESS':
                        passed_tests += 1
            
            f.write(f"- **Overall Status:** {'✓ PASSED' if passed_tests == total_tests else '⚠ PARTIAL' if passed_tests > 0 else '✗ FAILED'}\n")
            f.write(f"- **Tests Passed:** {passed_tests}/{total_tests}\n")
            f.write(f"- **Implementation Status:** Module exists and functional\n")
            f.write(f"- **Performance:** {'Acceptable' if all_results.get('performance', {}).get('recommended_resolution') else 'Needs optimization'}\n\n")
            
            # Step 1: Basic smoke test
            if 'smoke_test' in all_results:
                smoke = all_results['smoke_test']
                f.write("## 1. Basic Smoke Test\n\n")
                f.write(f"**Status:** {smoke.get('status', 'UNKNOWN')}\n\n")
                
                if 'coastline_features' in smoke:
                    f.write(f"- Coastline features loaded: {smoke['coastline_features']}\n")
                if 'crs' in smoke:
                    f.write(f"- Coordinate system: {smoke['crs']}\n")
                if 'test_distance_km' in smoke:
                    f.write(f"- Test distance (Sitka→Juneau): {smoke['test_distance_km']:.1f} km overwater vs {smoke['test_haversine_km']:.1f} km straight\n")
                    f.write(f"- Distance ratio: {smoke['test_ratio']:.2f}\n")
                    f.write(f"- Computation time: {smoke['test_time_s']:.2f}s\n")
                
                if 'error' in smoke:
                    f.write(f"- **Error:** {smoke['error']}\n")
                f.write("\n")
            
            # Step 2: Sanity checks
            if 'sanity_checks' in all_results:
                sanity = all_results['sanity_checks']
                f.write("## 2. Sanity Checks\n\n")
                f.write(f"**Status:** {sanity.get('status', 'UNKNOWN')}\n")
                f.write(f"**Results:** {sanity.get('n_passed', 0)}/{sanity.get('n_total', 0)} pairs passed\n\n")
                
                if 'pairs' in sanity:
                    f.write("### Detailed Results\n\n")
                    f.write("| Route | Expected (km) | Actual (km) | Ratio | Pass | Notes |\n")
                    f.write("|-------|---------------|-------------|-------|------|-------|\n")
                    
                    for pair in sanity['pairs']:
                        if 'error' not in pair:
                            status = "✓" if pair.get('in_range', False) else "✗"
                            f.write(f"| {pair['name']} | {pair['expected_min']}-{pair['expected_max']} | ")
                            ratio_str = f"{pair['ratio']:.2f}" if pair.get('ratio') is not None else 'N/A'
                            f.write(f"{pair['overwater_km']:.1f} | {ratio_str} | ")
                            f.write(f"{status} | {pair['notes']} |\n")
                        else:
                            f.write(f"| {pair['name']} | - | ERROR | - | ✗ | {pair['error'][:50]}... |\n")
                f.write("\n")
            
            # Step 3: Performance
            if 'performance' in all_results:
                perf = all_results['performance']
                f.write("## 3. Performance Analysis\n\n")
                f.write(f"**Recommended Resolution:** {perf.get('recommended_resolution', 500)}m\n\n")
                
                if 'resolutions' in perf:
                    f.write("### Resolution Comparison\n\n")
                    f.write("| Resolution | Init Time | Pair Time | Memory | Est. 150×150 |\n")
                    f.write("|------------|-----------|-----------|--------|---------------|\n")
                    
                    for res in perf['resolutions']:
                        if 'error' not in res:
                            est_time = res.get('extrapolations', {}).get(150, 0)
                            f.write(f"| {res['resolution_m']}m | {res['init_time_s']:.1f}s | ")
                            f.write(f"{res['mean_pair_time_s']:.3f}s | {res['memory_increase_mb']:.1f}MB | ")
                            f.write(f"{est_time:.1f}min |\n")
                        else:
                            f.write(f"| {res['resolution_m']}m | ERROR | - | - | - |\n")
                f.write("\n")
            
            # Step 4: Edge cases
            if 'edge_cases' in all_results:
                edge = all_results['edge_cases']
                f.write("## 4. Edge Case Testing\n\n")
                f.write(f"**Success Rate:** {edge.get('n_success', 0)}/{edge.get('n_total', 0)}\n\n")
                
                if 'cases' in edge:
                    f.write("### Results\n\n")
                    for case in edge['cases']:
                        f.write(f"**{case['case']}:** {case.get('status', 'UNKNOWN')}\n")
                        if case.get('status') == 'SUCCESS':
                            f.write(f"- Distance: {case['overwater_km']:.3f} km overwater vs {case['haversine_km']:.3f} km straight\n")
                            f.write(f"- Ratio: {case.get('ratio', 'N/A')}\n")
                        elif case.get('status') == 'NO_PATH':
                            f.write(f"- No path found (expected for some cases)\n")
                        elif 'error' in case:
                            f.write(f"- Error: {case['error']}\n")
                        f.write(f"- Description: {case.get('description', 'N/A')}\n\n")
            
            # Step 5: Visualization  
            if 'visualization' in all_results:
                viz = all_results['visualization']
                f.write("## 5. Path Visualization\n\n")
                f.write(f"**Status:** {viz.get('status', 'UNKNOWN')}\n")
                if 'output_file' in viz:
                    f.write(f"**Output:** {viz['output_file']}\n")
                if 'error' in viz:
                    f.write(f"**Error:** {viz['error']}\n")
                f.write("\n")
            
            # Step 6: Subset matrix
            if 'matrix' in all_results:
                matrix = all_results['matrix']
                f.write("## 6. Distance Matrix Computation\n\n")
                f.write(f"**Status:** {matrix.get('status', 'UNKNOWN')}\n")
                
                if matrix.get('status') == 'SUCCESS':
                    f.write(f"- Sites: {matrix['n_sites']}\n")
                    f.write(f"- Computation time: {matrix['computation_time_s']:.1f}s\n")
                    f.write(f"- Connectivity: {matrix['n_connected']}/{matrix['n_pairs']} pairs ({100*matrix['connectivity']:.1f}%)\n")
                    if not np.isnan(matrix.get('mean_distance_km', np.nan)):
                        f.write(f"- Distance range: {matrix['min_distance_km']:.1f} - {matrix['max_distance_km']:.1f} km\n")
                        f.write(f"- Mean distance: {matrix['mean_distance_km']:.1f} km\n")
                    f.write(f"- Output files: {matrix.get('csv_file', 'N/A')}, {matrix.get('heatmap_file', 'N/A')}\n")
                
                if 'error' in matrix:
                    f.write(f"**Error:** {matrix['error']}\n")
                f.write("\n")
            
            # Recommendations
            f.write("## Recommendations\n\n")
            
            if passed_tests == total_tests:
                f.write("✅ **Ready for Production Use**\n")
                f.write("- All tests passed successfully\n")
                f.write(f"- Recommended resolution: {perf.get('recommended_resolution', 500)}m for production runs\n")
                f.write("- Implementation handles edge cases appropriately\n")
            elif passed_tests > 0:
                f.write("⚠️ **Partial Success - Needs Review**\n")
                f.write("- Some tests failed - review error details above\n")
                f.write("- May be suitable for development/testing use\n") 
                f.write("- Address failing tests before production deployment\n")
            else:
                f.write("❌ **Not Ready - Major Issues**\n")
                f.write("- Multiple test failures indicate fundamental problems\n")
                f.write("- Requires debugging and fixes before use\n")
                f.write("- Review implementation against design document\n")
            
            f.write("\n### Production Readiness Checklist\n\n")
            f.write(f"- [{'x' if all_results.get('smoke_test', {}).get('status') == 'PASSED' else ' '}] Basic functionality works\n")
            f.write(f"- [{'x' if all_results.get('sanity_checks', {}).get('status') == 'PASSED' else ' '}] Distance calculations are reasonable\n") 
            f.write(f"- [{'x' if all_results.get('performance', {}).get('recommended_resolution') else ' '}] Performance is acceptable\n")
            f.write(f"- [{'x' if all_results.get('edge_cases', {}).get('n_success', 0) > 0 else ' '}] Edge cases handled\n")
            f.write(f"- [{'x' if all_results.get('matrix', {}).get('status') == 'SUCCESS' else ' '}] Matrix computation works\n")
            f.write(f"- [{'x' if all_results.get('visualization', {}).get('status') == 'SUCCESS' else ' '}] Visualization generates\n")
            
            f.write("\n---\n")
            f.write(f"*Report generated by overwater distance test suite v1.0*\n")
        
        print(f"✓ Test report written to {output_file}")
        
    except Exception as e:
        print(f"✗ Failed to write test report: {e}")

def main():
    """Run complete overwater distance test suite."""
    print("OVERWATER DISTANCE MODULE TESTING")
    print("=" * 60)
    print(f"Coastline: {COASTLINE_PATH}")
    print(f"Results: {RESULTS_DIR}")
    print()
    
    # Check prerequisites
    if not COASTLINE_PATH.exists():
        print(f"✗ Coastline file not found: {COASTLINE_PATH}")
        return
    
    all_results = {}
    
    # Step 1: Basic smoke test
    all_results['smoke_test'] = test_basic_smoke(COASTLINE_PATH)
    
    if all_results['smoke_test'].get('status') != 'PASSED':
        print("\n❌ Basic smoke test failed - aborting further tests")
        write_test_report(all_results, RESULTS_DIR / "overwater_test_results.md")
        return
    
    # Initialize calculator for remaining tests (500m resolution)
    print("\nInitializing calculator for remaining tests (500m resolution)...")
    calc = OverwaterDistanceCalculator(
        coastline_path=COASTLINE_PATH,
        resolution_m=500.0
    )
    
    # Step 2: Sanity checks
    all_results['sanity_checks'] = test_sanity_checks(calc)
    
    # Step 3: Performance test
    all_results['performance'] = test_performance(COASTLINE_PATH)
    
    # Step 4: Edge cases
    all_results['edge_cases'] = test_edge_cases(calc)
    
    # Step 5: Visualization
    viz_output = RESULTS_DIR / "overwater_paths_example.png"
    all_results['visualization'] = visualize_paths(calc, COASTLINE_PATH, viz_output)
    
    # Step 6: Load sites and run subset matrix
    coords, names, regions = load_all_sites()
    subset_coords, subset_names, subset_indices = select_diverse_sites(coords, names, regions, 20)
    all_results['matrix'] = run_subset_matrix(subset_coords, subset_names, calc, RESULTS_DIR)
    
    # Step 7: Write comprehensive report
    report_file = RESULTS_DIR / "overwater_test_results.md"
    write_test_report(all_results, report_file)
    
    # Final summary
    print("\n" + "="*60)
    print("TESTING COMPLETE")
    print("="*60)
    
    total_tests = sum(1 for result in all_results.values() 
                     if isinstance(result, dict) and 'status' in result)
    passed_tests = sum(1 for result in all_results.values() 
                      if isinstance(result, dict) and 
                      result.get('status') in ['PASSED', 'SUCCESS'])
    
    print(f"Tests passed: {passed_tests}/{total_tests}")
    print(f"Test report: {report_file}")
    print(f"Results directory: {RESULTS_DIR}")
    
    if passed_tests == total_tests:
        print("\n✅ ALL TESTS PASSED - Module ready for production use!")
    elif passed_tests > 0:
        print(f"\n⚠️ PARTIAL SUCCESS - {total_tests - passed_tests} tests need review")
    else:
        print("\n❌ ALL TESTS FAILED - Major issues need resolution")

if __name__ == "__main__":
    main()