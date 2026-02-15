#!/usr/bin/env python3
"""
Very simple test of overwater distance with just the minimal components.
"""

import time
import numpy as np
import geopandas as gpd
from pathlib import Path
from shapely.geometry import Point

def simple_test():
    print("Simple overwater distance test")
    print("=" * 40)
    
    # Load a small subset of coastline around Sitka/Juneau
    coastline_path = Path('data/shorelines/ne_pacific_coastline.geojson')
    print(f"Loading coastline: {coastline_path}")
    
    start = time.time()
    full_coastline = gpd.read_file(coastline_path)
    print(f"Loaded {len(full_coastline)} features in {time.time() - start:.2f}s")
    
    # Test coordinates
    sitka = (-135.3300, 57.0531)  # lon, lat for bounds
    juneau = (-134.4197, 58.3019)
    
    # Create a bounding box around our test sites
    buffer = 1.0  # degrees
    minx, miny = min(sitka[0], juneau[0]) - buffer, min(sitka[1], juneau[1]) - buffer  
    maxx, maxy = max(sitka[0], juneau[0]) + buffer, max(sitka[1], juneau[1]) + buffer
    
    print(f"Clipping to bounds: {minx:.3f}, {miny:.3f}, {maxx:.3f}, {maxy:.3f}")
    
    start = time.time()
    clipped_coastline = full_coastline.cx[minx:maxx, miny:maxy]
    print(f"Clipped to {len(clipped_coastline)} features in {time.time() - start:.2f}s")
    
    if len(clipped_coastline) == 0:
        print("ERROR: No coastline features in test region!")
        return
    
    # Save clipped coastline for debugging
    test_coastline_path = Path('data/shorelines/test_region.geojson')
    clipped_coastline.to_file(test_coastline_path)
    print(f"Saved test region to {test_coastline_path}")
    
    # Now test with the overwater calculator on the small region
    try:
        print("\nTesting OverwaterDistanceCalculator with small region...")
        from sswd_evoepi.overwater_distance import OverwaterDistanceCalculator
        
        start = time.time()
        calc = OverwaterDistanceCalculator(
            coastline_path=test_coastline_path,
            resolution_m=1000.0  # 1km resolution
        )
        print(f"Calculator initialized in {time.time() - start:.2f}s")
        
        # Test distance calculation
        print(f"\nComputing Sitka → Juneau distance...")
        start = time.time()
        distance = calc.compute_single_distance((57.0531, -135.3300), (58.3019, -134.4197))
        elapsed = time.time() - start
        
        print(f"Distance: {distance:.1f} km")
        print(f"Computation time: {elapsed:.2f}s")
        
        # Compare to Haversine
        from math import radians, sin, cos, asin, sqrt
        def haversine_simple(lat1, lon1, lat2, lon2):
            R = 6371
            lat1, lon1, lat2, lon2 = map(radians, [lat1, lon1, lat2, lon2])
            dlat, dlon = lat2 - lat1, lon2 - lon1
            a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
            return 2 * R * asin(sqrt(a))
        
        haversine_dist = haversine_simple(57.0531, -135.3300, 58.3019, -134.4197)
        print(f"Haversine distance: {haversine_dist:.1f} km")
        print(f"Ratio: {distance / haversine_dist:.2f}")
        
        if np.isfinite(distance) and distance > 0:
            print("✅ SUCCESS: Basic overwater calculation works!")
            return True
        else:
            print("❌ FAILED: Invalid distance result")
            return False
        
    except Exception as e:
        print(f"❌ ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = simple_test()
    if success:
        print("\n🎉 Basic test passed - module is functional!")
    else:
        print("\n💥 Basic test failed - needs debugging")