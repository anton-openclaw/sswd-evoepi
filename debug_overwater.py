#!/usr/bin/env python3
"""
Debug the overwater distance implementation step by step.
"""

import time
import numpy as np
from pathlib import Path
from sswd_evoepi.overwater_distance import OverwaterDistanceCalculator

def debug_step_by_step():
    coastline_path = Path('data/shorelines/ne_pacific_coastline.geojson')
    
    print("1. Initialize calculator...")
    start = time.time()
    calc = OverwaterDistanceCalculator(coastline_path, resolution_m=5000.0)  # Very coarse for speed
    print(f"   Initialization took {time.time() - start:.2f}s")
    print(f"   CRS: {calc.crs}")
    
    # Test coordinates (close together for simple test)
    sitka = (57.0531, -135.3300)
    juneau = (58.3019, -134.4197)
    
    print(f"\n2. Test rasterization for region around sites...")
    site_coords = np.array([sitka, juneau])
    
    start = time.time()
    bounds = calc._get_processing_bounds(site_coords)
    print(f"   Processing bounds: {bounds}")
    print(f"   Bounds calculation took {time.time() - start:.2f}s")
    
    start = time.time()
    raster, transform = calc.rasterize_coastline(bounds)
    print(f"   Rasterization took {time.time() - start:.2f}s")
    print(f"   Raster shape: {raster.shape}")
    print(f"   Sea cells: {np.sum(raster):,} / {raster.size:,} ({100*np.sum(raster)/raster.size:.1f}%)")
    
    # Cache the raster
    calc._raster = raster
    calc._raster_transform = transform
    calc._raster_bounds = bounds
    
    print(f"\n3. Test sea graph construction...")
    start = time.time()
    sea_graph = calc.build_sea_graph(raster)
    print(f"   Graph construction took {time.time() - start:.2f}s")
    print(f"   Graph shape: {sea_graph.shape}")
    print(f"   Graph edges: {sea_graph.nnz:,}")
    
    # Cache the graph
    calc._sea_graph = sea_graph
    
    print(f"\n4. Test site mapping...")
    start = time.time()
    site_mask, cell_indices = calc.map_sites_to_cells(site_coords)
    print(f"   Site mapping took {time.time() - start:.2f}s")
    print(f"   Valid sites: {np.sum(site_mask)}/{len(site_coords)}")
    print(f"   Cell indices: {cell_indices}")
    
    if len(cell_indices) >= 2:
        print(f"\n5. Test Dijkstra computation...")
        from scipy.sparse.csgraph import dijkstra
        
        start = time.time()
        # Just compute distances from first site
        distances = dijkstra(sea_graph, directed=False, indices=[cell_indices[0]])
        print(f"   Dijkstra took {time.time() - start:.2f}s")
        
        # Distance between the two sites
        dist_cells = distances[0, cell_indices[1]]
        dist_km = dist_cells * (calc.resolution_m / 1000.0)
        print(f"   Distance: {dist_km:.1f} km ({dist_cells:.1f} cells)")
        
        # Compare to Haversine
        from sswd_evoepi.overwater_distance import haversine_km
        haver_dist = haversine_km(sitka[0], sitka[1], juneau[0], juneau[1])
        print(f"   Haversine: {haver_dist:.1f} km")
        print(f"   Ratio: {dist_km / haver_dist:.2f}")
        
    else:
        print("   Cannot test Dijkstra - insufficient valid sites")
    
    print("\nDebug complete!")

if __name__ == "__main__":
    debug_step_by_step()