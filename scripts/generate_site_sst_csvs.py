#!/usr/bin/env python3
"""Generate per-site monthly SST CSVs from grid_cell_monthly.json + site_to_grid.json.

The model expects files like data/sst/site_sst/AK-AL-004_monthly.csv for each node.
Our SST data is stored as grid cell data. This script creates the per-site CSVs
by mapping each site to its OISST grid cell.
"""

import json
import csv
import os
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent
SST_DIR = PROJECT_ROOT / 'data' / 'sst' / 'site_sst'

def main():
    # Load mappings
    with open(SST_DIR / 'site_to_grid.json') as f:
        site_to_grid = json.load(f)  # index → grid_cell_key
    
    with open(SST_DIR / 'grid_cell_monthly.json') as f:
        grid_data = json.load(f)  # grid_cell_key → {monthly: {year: {month: sst}}, ...}
    
    with open(PROJECT_ROOT / 'data' / 'nodes' / 'all_sites.json') as f:
        sites = json.load(f)
    
    print(f"Sites: {len(sites)}")
    print(f"Grid cells: {len(grid_data)}")
    print(f"Mappings: {len(site_to_grid)}")
    
    # Build a spatial lookup: for each site, find the nearest OISST grid cell (0.25°)
    # OISST cells are at 0.25° resolution, snap lat/lon to nearest 0.25
    def snap_to_grid(lat, lon):
        """Snap to nearest 0.25° OISST grid cell."""
        glat = round(lat * 4) / 4
        glon = round(lon * 4) / 4
        return f"{glat:.4f},{glon:.4f}"
    
    generated = 0
    skipped = 0
    
    for i, site in enumerate(sites):
        name = site['name']
        
        # Find grid cell by snapping coordinates
        grid_key = snap_to_grid(site['latitude'], site['longitude'])
        
        cell = grid_data.get(grid_key)
        if cell is None:
            print(f"  SKIP {name}: grid cell {grid_key} not in data")
            skipped += 1
            continue
        
        monthly = cell.get('monthly', [])
        if not monthly:
            print(f"  SKIP {name}: no monthly data for cell {grid_key}")
            skipped += 1
            continue
        
        # Write CSV — monthly is a list of {year, month, sst_mean, ...}
        csv_path = SST_DIR / f"{name}_monthly.csv"
        with open(csv_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['year', 'month', 'sst'])
            
            # Sort by year then month
            sorted_records = sorted(monthly, key=lambda r: (r['year'], r['month']))
            for rec in sorted_records:
                sst = rec.get('sst_mean')
                if sst is not None:
                    writer.writerow([rec['year'], rec['month'], f"{sst:.3f}"])
        
        generated += 1
    
    print(f"\nGenerated: {generated} CSVs")
    print(f"Skipped: {skipped}")
    print(f"Output dir: {SST_DIR}")


if __name__ == '__main__':
    main()
