#!/usr/bin/env python3
"""
Debug exactly where the code is hanging.
"""

import sys
import time
import signal
import numpy as np
from pathlib import Path

def timeout_handler(signum, frame):
    print(f"\n⏰ TIMEOUT after {time.time() - start_time:.1f}s")
    sys.exit(1)

# Set up timeout
signal.signal(signal.SIGALRM, timeout_handler)
signal.alarm(30)  # 30 second timeout

start_time = time.time()

print("DEBUGGING OVERWATER DISTANCE HANG")
print("=" * 40)

try:
    print("1. Importing modules...", end=' ', flush=True)
    from sswd_evoepi.overwater_distance import OverwaterDistanceCalculator
    print(f"OK ({time.time() - start_time:.2f}s)")
    
    coastline_path = Path('data/shorelines/ne_pacific_coastline.geojson')
    
    print("2. Creating calculator...", end=' ', flush=True)
    step_start = time.time()
    calc = OverwaterDistanceCalculator(coastline_path, resolution_m=5000.0)
    print(f"OK ({time.time() - step_start:.2f}s)")
    
    print("3. Getting processing bounds...", end=' ', flush=True)
    step_start = time.time()
    site_coords = np.array([(57.0531, -135.3300), (58.3019, -134.4197)])
    bounds = calc._get_processing_bounds(site_coords)
    print(f"OK ({time.time() - step_start:.2f}s)")
    print(f"   Bounds: {bounds}")
    
    # Calculate expected grid size
    minx, miny, maxx, maxy = bounds
    width = int(np.ceil((maxx - minx) / calc.resolution_m))
    height = int(np.ceil((maxy - miny) / calc.resolution_m))
    total_cells = width * height
    print(f"   Expected grid: {width} × {height} = {total_cells:,} cells ({total_cells*4/1e6:.1f} MB)")
    
    if total_cells > 10_000_000:  # 10M cells
        print("   ⚠️  Grid is very large - this could cause hanging!")
    
    print("4. Rasterizing coastline...", end=' ', flush=True)
    step_start = time.time()
    raster, transform = calc.rasterize_coastline(bounds)
    print(f"OK ({time.time() - step_start:.2f}s)")
    print(f"   Actual grid shape: {raster.shape}")
    
    print("5. Building sea graph...", end=' ', flush=True)
    step_start = time.time()
    sea_graph = calc.build_sea_graph(raster)
    print(f"OK ({time.time() - step_start:.2f}s)")
    print(f"   Graph: {sea_graph.shape[0]:,} nodes, {sea_graph.nnz:,} edges")
    
    print("6. Mapping sites to cells...", end=' ', flush=True)
    step_start = time.time()
    calc._raster = raster
    calc._raster_transform = transform
    calc._raster_bounds = bounds
    site_mask, cell_indices = calc.map_sites_to_cells(site_coords)
    print(f"OK ({time.time() - step_start:.2f}s)")
    print(f"   Mapped sites: {np.sum(site_mask)}/{len(site_coords)}")
    
    if len(cell_indices) >= 2:
        print("7. Running Dijkstra...", end=' ', flush=True)
        step_start = time.time()
        from scipy.sparse.csgraph import dijkstra
        distances = dijkstra(sea_graph, directed=False, indices=[cell_indices[0]])
        print(f"OK ({time.time() - step_start:.2f}s)")
        
        dist_cells = distances[0, cell_indices[1]]
        dist_km = dist_cells * (calc.resolution_m / 1000.0)
        print(f"   Distance: {dist_km:.1f} km")
        
        print("✅ ALL STEPS COMPLETED - no hanging detected!")
    else:
        print("❌ Cannot test Dijkstra - sites not mapped to sea cells")
        
except Exception as e:
    print(f"❌ ERROR at {time.time() - start_time:.1f}s: {e}")
    import traceback
    traceback.print_exc()
    
signal.alarm(0)  # Cancel timeout