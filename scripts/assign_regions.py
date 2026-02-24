#!/usr/bin/env python3
"""
Assign regions to all sites using Natural Earth state/province shapefiles
+ Salish Sea polygon + geographic subdivision rules.

Strategy:
1. Load Natural Earth 10m admin-1 (states/provinces) boundaries
2. For each site, find the nearest state/province (handles marine sites)
3. Map state/province to our 13 model regions using geographic sub-rules
4. Salish Sea polygon overrides WA/BC assignments for inland waterway sites

The 13 model regions:
  AK-AL  Aleutian Islands (< ~55°N or < ~-164°W in Alaska)
  AK-WG  Western Gulf of Alaska (Kodiak / Cook Inlet, -159° to -150°)
  AK-EG  Eastern Gulf of Alaska (PWS / Yakutat, -150° to -138°)
  AK-SE  Southeast Alaska panhandle (south of Cross Sound)
  BC-N   Northern BC outer coast (> ~51.5°N, outside Salish Sea)
  BC-C   Central BC outer coast (≤ ~51.5°N, outside Salish Sea)
  SS     Salish Sea (Puget Sound + JdF + San Juans + Georgia Strait)
  WA-O   Washington outer coast (Pacific side)
  OR     Oregon coast
  CA-N   Northern California (> ~38.5°N)
  CA-C   Central California (~34.5° to ~38.5°N)
  CA-S   Southern California (< ~34.5°N)
  BJ     Baja California (Mexico)
"""

import json
import sys
from pathlib import Path
from collections import Counter

import geopandas as gpd
import numpy as np
from shapely.geometry import Point, Polygon
from shapely.ops import nearest_points

BASE = Path(__file__).parent.parent
SITES_FILE = BASE / "data" / "nodes" / "all_sites.json"
SHAPEFILE = BASE / "data" / "shapefiles" / "ne_10m_admin_1" / "ne_10m_admin_1_states_provinces.shp"


def get_lat_lon(site):
    """Extract lat/lon from site dict (handles both schemas)."""
    lat = site.get('lat', site.get('latitude', None))
    lon = site.get('lon', site.get('longitude', None))
    return lat, lon


def define_salish_sea_polygon():
    """
    The Salish Sea boundary polygon (Puget Sound + JdF + Georgia Strait
    + Desolation Sound). Takes priority over WA/BC state assignments.
    Uses (lon, lat) order for shapely compatibility.
    """
    # Define in (lon, lat) for shapely
    vertices = [
        # JdF entrance south side, going clockwise
        (-124.30, 48.15),
        (-123.70, 48.05),
        (-123.40, 48.00),
        (-123.10, 47.95),
        (-122.90, 47.95),
        # Hood Canal
        (-123.00, 47.80),
        (-123.30, 47.35),
        (-123.15, 47.00),
        # South Puget Sound
        (-122.95, 46.85),
        (-122.35, 46.85),
        # East side going north
        (-122.30, 47.10),
        (-122.15, 47.30),
        (-122.05, 47.55),
        (-122.05, 47.70),
        (-122.05, 47.90),
        (-122.00, 48.10),
        (-122.10, 48.30),
        (-122.15, 48.50),
        (-122.30, 48.75),
        (-122.40, 49.05),
        # Lower mainland
        (-122.40, 49.10),
        (-122.50, 49.25),
        (-122.70, 49.40),
        (-123.00, 49.50),
        (-123.10, 49.60),
        (-123.30, 49.75),
        # Up Georgia Strait mainland
        (-123.80, 49.85),
        (-123.80, 49.95),
        (-123.70, 50.15),
        (-124.30, 50.30),
        (-124.50, 50.60),
        # Desolation Sound / Discovery Islands
        (-125.00, 50.80),
        (-125.70, 50.80),
        # Down VI side
        (-125.40, 50.20),
        (-125.30, 50.00),
        (-125.00, 49.80),
        (-124.90, 49.55),
        (-124.40, 49.30),
        (-124.15, 49.10),
        (-123.60, 48.85),
        (-123.60, 48.75),
        (-123.50, 48.55),
        (-123.90, 48.50),
        (-124.20, 48.45),
        (-124.50, 48.55),
        # Back to JdF entrance north
        (-124.50, 48.45),
        (-124.30, 48.15),  # close
    ]
    return Polygon(vertices)


def load_admin_boundaries():
    """Load and filter Natural Earth admin-1 boundaries for our region."""
    gdf = gpd.read_file(SHAPEFILE)
    
    states = {
        'Alaska': gdf[gdf['name'] == 'Alaska'].iloc[0].geometry,
        'Washington': gdf[gdf['name'] == 'Washington'].iloc[0].geometry,
        'Oregon': gdf[gdf['name'] == 'Oregon'].iloc[0].geometry,
        'California': gdf[gdf['name'] == 'California'].iloc[0].geometry,
        'British Columbia': gdf[(gdf['name'] == 'British Columbia')].iloc[0].geometry,
        'Baja California': gdf[gdf['name'] == 'Baja California'].iloc[0].geometry,
        'Baja California Sur': gdf[gdf['name'] == 'Baja California Sur'].iloc[0].geometry,
    }
    
    return states


def find_nearest_state(lon, lat, states):
    """
    Find the nearest state/province for a marine or land point.
    
    Returns (state_name, distance_degrees).
    For points crossing the 180° meridian (western Aleutians),
    we normalize longitude.
    """
    point = Point(lon, lat)
    
    # First check direct containment (fast path for land points)
    for name, geom in states.items():
        if geom.contains(point):
            return name, 0.0
    
    # Marine point — find nearest land
    min_dist = float('inf')
    nearest_state = None
    
    for name, geom in states.items():
        dist = point.distance(geom)
        if dist < min_dist:
            min_dist = dist
            nearest_state = name
    
    return nearest_state, min_dist


def alaska_subregion(lat, lon):
    """
    Subdivide Alaska into AK-AL, AK-WG, AK-EG, AK-SE.
    
    Geographic rules:
    - AK-AL: Aleutian Islands chain (lat < 55° AND lon < -164°, 
              OR positive longitude / crossing dateline)
    - AK-SE: Southeast panhandle (lat < 58.5° AND lon > -138°,
              plus Yakutat area)
    - AK-EG: Eastern Gulf (lon > -150° AND not SE, OR Yakutat corridor)
    - AK-WG: Western Gulf (everything else: Kodiak, Cook Inlet, Peninsula)
    """
    # Normalize positive longitudes (western Aleutians)
    norm_lon = lon if lon < 0 else lon - 360.0
    
    # Aleutian Islands: far west or crossing dateline
    if norm_lon < -164.0 and lat < 56.0:
        return 'AK-AL'
    if lon > 0:  # Crossing dateline
        return 'AK-AL'
    
    # Southeast Alaska panhandle  
    # The panhandle is roughly: south of Cross Sound (~58.5°N) 
    # and east of ~-138°W. Glacier Bay / Icy Strait transitions.
    if lat < 58.5 and lon > -138.0:
        return 'AK-SE'
    # Yakutat is at ~59.5°N, -139.7°W — traditionally SE Alaska
    if lat < 60.0 and lon > -140.5 and lon < -138.0:
        return 'AK-SE'
    
    # Eastern Gulf: PWS, Copper River, Kenai east side
    # Roughly -150° to -138°
    if lon > -150.0:
        return 'AK-EG'
    
    # Western Gulf: Kodiak, Cook Inlet, Alaska Peninsula
    return 'AK-WG'


def bc_subregion(lat, lon, salish_sea):
    """
    Subdivide British Columbia into BC-N, BC-C, SS.
    """
    point = Point(lon, lat)
    
    if salish_sea.contains(point):
        return 'SS'
    
    if lat > 51.5:
        return 'BC-N'
    
    return 'BC-C'


def wa_subregion(lat, lon, salish_sea):
    """
    Subdivide Washington into WA-O or SS.
    """
    point = Point(lon, lat)
    
    if salish_sea.contains(point):
        return 'SS'
    
    return 'WA-O'


def ca_subregion(lat, lon):
    """
    Subdivide California into CA-N, CA-C, CA-S.
    """
    if lat > 38.5:
        return 'CA-N'
    if lat > 34.5:
        return 'CA-C'
    return 'CA-S'


def assign_region(lat, lon, states, salish_sea):
    """
    Full region assignment pipeline:
    1. Find nearest state/province
    2. Apply sub-regional rules
    3. Salish Sea override
    """
    # For dateline-crossing sites (positive lon), try BOTH the original
    # positive lon (NE Alaska polygon extends to 179.78°) AND normalized negative.
    # Use whichever gives a closer match.
    state_name, dist = find_nearest_state(lon, lat, states)
    
    if lon > 0:
        # Also try negative normalization
        alt_state, alt_dist = find_nearest_state(lon - 360.0, lat, states)
        if alt_dist < dist:
            state_name, dist = alt_state, alt_dist
    
    if state_name is None:
        return 'UNKNOWN'
    
    # Distance sanity check (> ~5° away from any land → UNKNOWN)
    if dist > 5.0:
        return 'UNKNOWN'
    
    # Map state → region with sub-rules
    if state_name == 'Alaska':
        return alaska_subregion(lat, lon)
    elif state_name == 'British Columbia':
        return bc_subregion(lat, lon, salish_sea)
    elif state_name == 'Washington':
        return wa_subregion(lat, lon, salish_sea)
    elif state_name == 'Oregon':
        return 'OR'
    elif state_name == 'California':
        return ca_subregion(lat, lon)
    elif state_name in ('Baja California', 'Baja California Sur'):
        return 'BJ'
    else:
        return 'UNKNOWN'


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Reassign regions using shapefiles')
    parser.add_argument('--dry-run', action='store_true',
                        help='Show changes without writing')
    parser.add_argument('--apply', action='store_true',
                        help='Write changes to all_sites.json')
    parser.add_argument('--show-changes', action='store_true',
                        help='Show individual site reassignments')
    parser.add_argument('--show-unknown', action='store_true',
                        help='Show sites that could not be assigned')
    args = parser.parse_args()
    
    if not args.dry_run and not args.apply:
        args.dry_run = True
        print("(No --apply flag, running in dry-run mode)\n")
    
    # Load sites
    sites = json.load(open(SITES_FILE))
    print(f"Loaded {len(sites)} sites from {SITES_FILE.name}")
    
    # Load boundaries
    print("Loading Natural Earth admin-1 boundaries...")
    states = load_admin_boundaries()
    print(f"  Loaded {len(states)} state/province boundaries")
    
    salish_sea = define_salish_sea_polygon()
    print(f"  Salish Sea polygon defined")
    
    # Current distribution
    old_counts = Counter(s.get('region', 'NONE') for s in sites)
    
    # Reassign all sites
    print(f"\nAssigning regions to {len(sites)} sites...")
    changes = []
    unknowns = []
    new_assignments = []
    
    for i, site in enumerate(sites):
        lat, lon = get_lat_lon(site)
        if lat is None or lon is None:
            unknowns.append((i, site.get('name', '?'), 'NO_COORDS'))
            new_assignments.append(site.get('region', 'UNKNOWN'))
            continue
        
        new_region = assign_region(lat, lon, states, salish_sea)
        old_region = site.get('region', 'NONE')
        new_assignments.append(new_region)
        
        if new_region == 'UNKNOWN':
            unknowns.append((i, site.get('name', '?'), f'{lat:.3f}, {lon:.3f}'))
        
        if new_region != old_region:
            changes.append((i, site.get('name', '?'), old_region, new_region, lat, lon))
    
    # New distribution
    new_counts = Counter(new_assignments)
    
    # Report
    print("\n=== Region counts: before → after ===")
    all_region_names = sorted(set(list(old_counts.keys()) + list(new_counts.keys())))
    for r in all_region_names:
        old = old_counts.get(r, 0)
        new = new_counts.get(r, 0)
        delta = new - old
        marker = " ←" if delta != 0 else ""
        print(f"  {r:6s}: {old:4d} → {new:4d}  ({delta:+d}){marker}")
    
    print(f"\n{len(changes)} sites would change region")
    print(f"{len(unknowns)} sites unassigned (UNKNOWN)")
    
    if args.show_changes and changes:
        print("\n=== Changes ===")
        for i, name, old_r, new_r, lat, lon in changes:
            print(f"  [{i}] {name[:50]:50s} {old_r:6s} → {new_r:6s}  ({lat:.3f}, {lon:.3f})")
    
    if args.show_unknown and unknowns:
        print("\n=== Unknown ===")
        for i, name, coords in unknowns:
            print(f"  [{i}] {name[:50]:50s} {coords}")
    
    if args.apply and changes:
        for i, name, old_r, new_r, lat, lon in changes:
            sites[i]['region'] = new_r
        
        with open(SITES_FILE, 'w') as f:
            json.dump(sites, f, indent=2)
        print(f"\n✓ Updated {len(changes)} sites in {SITES_FILE.name}")
    elif not args.apply:
        print("\nRun with --apply to write changes.")


if __name__ == '__main__':
    main()
