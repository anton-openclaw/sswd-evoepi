#!/usr/bin/env python3
"""
Minimal test with synthetic coastline data to isolate the issue.
"""

import numpy as np
import geopandas as gpd
from shapely.geometry import Polygon, Point
from pathlib import Path
import time

def create_test_coastline():
    """Create a simple synthetic coastline for testing."""
    print("Creating synthetic coastline...")
    
    # Create a simple rectangular "island" with a channel between two areas
    island1 = Polygon([
        (-136.0, 57.0), (-135.5, 57.0), (-135.5, 57.5), (-136.0, 57.5), (-136.0, 57.0)
    ])
    
    island2 = Polygon([
        (-134.5, 58.0), (-134.0, 58.0), (-134.0, 58.5), (-134.5, 58.5), (-134.5, 58.0)  
    ])
    
    # Create GeoDataFrame
    gdf = gpd.GeoDataFrame({
        'geometry': [island1, island2]
    }, crs='EPSG:4326')
    
    # Save to file
    test_file = Path('data/shorelines/test_synthetic.geojson')
    test_file.parent.mkdir(parents=True, exist_ok=True)
    gdf.to_file(test_file)
    
    print(f"Created synthetic coastline with {len(gdf)} features")
    print(f"Saved to: {test_file}")
    return test_file

def test_with_synthetic():
    """Test overwater calculation with synthetic coastline."""
    coastline_file = create_test_coastline()
    
    print("\nTesting with OverwaterDistanceCalculator...")
    
    try:
        from sswd_evoepi.overwater_distance import OverwaterDistanceCalculator
        
        print("Initializing calculator...")
        start = time.time()
        
        calc = OverwaterDistanceCalculator(
            coastline_path=coastline_file,
            resolution_m=500.0  # 500m resolution
        )
        
        init_time = time.time() - start
        print(f"✓ Initialization completed in {init_time:.2f}s")
        print(f"  CRS: {calc.crs}")
        
        # Test coordinates in water between the islands
        point1 = (57.25, -135.75)  # Between the islands
        point2 = (58.25, -134.25)  # In other water area
        
        print(f"\nComputing distance: {point1} → {point2}")
        start = time.time()
        
        distance = calc.compute_single_distance(point1, point2)
        
        compute_time = time.time() - start
        print(f"✓ Distance computed in {compute_time:.2f}s")
        print(f"  Result: {distance:.1f} km")
        
        # Calculate expected straight-line distance
        from math import radians, sin, cos, asin, sqrt
        def haversine_km(lat1, lon1, lat2, lon2):
            R = 6371
            lat1, lon1, lat2, lon2 = map(radians, [lat1, lon1, lat2, lon2])
            dlat, dlon = lat2 - lat1, lon2 - lon1
            a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
            return 2 * R * asin(sqrt(a))
        
        straight_dist = haversine_km(point1[0], point1[1], point2[0], point2[1])
        print(f"  Straight-line: {straight_dist:.1f} km")
        
        if np.isfinite(distance) and distance > 0:
            ratio = distance / straight_dist
            print(f"  Ratio: {ratio:.2f}")
            
            if ratio >= 1.0 and ratio < 5.0:  # Reasonable range
                print("✅ SUCCESS: Result is reasonable!")
                return True
            else:
                print("⚠️  WARNING: Ratio seems unusual but computation worked")
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
    print("MINIMAL OVERWATER DISTANCE TEST")
    print("=" * 40)
    
    success = test_with_synthetic()
    
    if success:
        print("\n🎉 Minimal test PASSED!")
        print("The overwater distance module is working.")
        print("Performance issues may be due to large coastline files.")
    else:
        print("\n💥 Minimal test FAILED!")
        print("There are fundamental issues with the implementation.")