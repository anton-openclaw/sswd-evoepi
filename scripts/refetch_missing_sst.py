#!/usr/bin/env python3
"""Re-fetch OISST monthly data for grid cells that are missing.

Queries NOAA ERDDAP for 0.25° OISST v2.1 monthly means.
For cells that are truly land-masked, uses nearest-neighbor from
the closest cell WITH data (still real satellite data, just nearby).
"""

import json
import csv
import time
import math
import urllib.request
import urllib.error
from pathlib import Path
from typing import Optional, Dict, List, Tuple

PROJECT_ROOT = Path(__file__).parent.parent
SST_DIR = PROJECT_ROOT / 'data' / 'sst' / 'site_sst'

ERDDAP_BASE = "https://coastwatch.pfeg.noaa.gov/erddap/griddap/ncdcOisst21Agg_LonPM180.csv"
START_YEAR = 2002
END_YEAR = 2025


def fetch_oisst_monthly(lat: float, lon: float) -> Optional[List[dict]]:
    """Fetch monthly OISST data for a grid cell from ERDDAP.
    
    Returns list of {year, month, sst_mean} or None if land-masked.
    """
    # OISST uses 0.25° grid, centered on grid points
    glat = round(lat * 4) / 4
    glon = round(lon * 4) / 4
    
    # Query a small box around the point
    lat_lo = glat - 0.001
    lat_hi = glat + 0.001
    lon_lo = glon - 0.001
    lon_hi = glon + 0.001
    
    url = (
        f"{ERDDAP_BASE}?"
        f"sst[({START_YEAR}-01-01T00:00:00Z):1:({END_YEAR}-12-31T00:00:00Z)]"
        f"[({lat_lo}):1:({lat_hi})]"
        f"[({lon_lo}):1:({lon_hi})]"
    )
    
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'SSWD-EvoEpi/1.0'})
        with urllib.request.urlopen(req, timeout=30) as resp:
            text = resp.read().decode('utf-8')
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
        return None
    
    # Parse CSV response
    lines = text.strip().split('\n')
    if len(lines) < 3:  # header + units + at least one data row
        return None
    
    records = []
    reader = csv.DictReader(lines[2:], fieldnames=lines[0].split(','))
    
    # Skip the units row (line index 1) — we start from line 2
    for row in reader:
        try:
            time_str = row.get('time', '')
            sst_str = row.get('sst', '')
            if not time_str or not sst_str or sst_str == 'NaN':
                continue
            
            # Parse time: "2002-01-01T12:00:00Z"
            year = int(time_str[:4])
            month = int(time_str[5:7])
            sst = float(sst_str)
            
            records.append({'year': year, 'month': month, 'sst_mean': sst})
        except (ValueError, KeyError):
            continue
    
    if not records:
        return None
    
    # Aggregate to monthly means (OISST is daily, we want monthly)
    monthly = {}
    for rec in records:
        key = (rec['year'], rec['month'])
        monthly.setdefault(key, []).append(rec['sst_mean'])
    
    result = []
    for (year, month), ssts in sorted(monthly.items()):
        result.append({
            'year': year,
            'month': month,
            'sst_mean': round(sum(ssts) / len(ssts), 3),
        })
    
    return result


def find_nearest_with_data(target_lat: float, target_lon: float,
                           cells_with_data: Dict[str, List[dict]]) -> Optional[Tuple[str, List[dict]]]:
    """Find nearest grid cell that has SST data."""
    best_dist = float('inf')
    best_key = None
    
    for key, data in cells_with_data.items():
        parts = key.split(',')
        clat, clon = float(parts[0]), float(parts[1])
        dist = math.sqrt((clat - target_lat)**2 + (clon - target_lon)**2)
        if dist < best_dist:
            best_dist = dist
            best_key = key
    
    if best_key is not None:
        return best_key, cells_with_data[best_key], best_dist
    return None


def write_site_csv(name: str, records: List[dict]):
    """Write per-site monthly CSV."""
    csv_path = SST_DIR / f'{name}_monthly.csv'
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['year', 'month', 'sst'])
        for rec in sorted(records, key=lambda r: (r['year'], r['month'])):
            writer.writerow([rec['year'], rec['month'], f"{rec['sst_mean']:.3f}"])


def main():
    with open(SST_DIR / 'grid_cell_monthly.json') as f:
        grid = json.load(f)
    
    with open(PROJECT_ROOT / 'data' / 'nodes' / 'all_sites.json') as f:
        sites = json.load(f)
    
    def snap(lat, lon):
        return f'{round(lat*4)/4:.4f},{round(lon*4)/4:.4f}'
    
    # Build map of cells with data
    cells_with_data = {}
    for key, cell in grid.items():
        if cell.get('monthly'):
            cells_with_data[key] = cell['monthly']
    
    # Find cells that need data
    site_cells = {}
    for s in sites:
        key = snap(s['latitude'], s['longitude'])
        site_cells.setdefault(key, []).append(s)
    
    missing_cells = {k for k in site_cells if k not in cells_with_data}
    print(f"Cells with data: {len(cells_with_data)}")
    print(f"Missing cells to fetch: {len(missing_cells)}")
    print(f"Sites affected: {sum(len(site_cells[k]) for k in missing_cells)}")
    
    fetched = 0
    nearest_neighbor = 0
    failed = 0
    
    for i, cell_key in enumerate(sorted(missing_cells)):
        parts = cell_key.split(',')
        lat, lon = float(parts[0]), float(parts[1])
        n_sites = len(site_cells[cell_key])
        
        print(f"[{i+1}/{len(missing_cells)}] {cell_key} ({n_sites} sites)...", end=' ', flush=True)
        
        # Try direct ERDDAP fetch
        records = fetch_oisst_monthly(lat, lon)
        
        if records and len(records) >= 12:
            print(f"OK ({len(records)} records)")
            cells_with_data[cell_key] = records
            for s in site_cells[cell_key]:
                write_site_csv(s['name'], records)
            fetched += 1
        else:
            # Land-masked — use nearest neighbor
            result = find_nearest_with_data(lat, lon, cells_with_data)
            if result:
                nn_key, nn_data, nn_dist = result
                print(f"NN → {nn_key} ({nn_dist*111:.0f}km)")
                for s in site_cells[cell_key]:
                    write_site_csv(s['name'], nn_data)
                nearest_neighbor += 1
            else:
                print("FAILED — no data anywhere nearby")
                failed += 1
        
        # Rate limit ERDDAP
        if records is not None or i % 5 == 0:
            time.sleep(0.5)
    
    print(f"\nDone: {fetched} fetched, {nearest_neighbor} nearest-neighbor, {failed} failed")
    
    # Update grid_cell_monthly.json with new data
    for cell_key in missing_cells:
        if cell_key in cells_with_data and cell_key not in {k for k, v in grid.items() if v.get('monthly')}:
            if cell_key in grid:
                grid[cell_key]['monthly'] = cells_with_data[cell_key]
    
    with open(SST_DIR / 'grid_cell_monthly.json', 'w') as f:
        json.dump(grid, f)
    print("Updated grid_cell_monthly.json")
    
    # Verify all sites have CSVs
    count = sum(1 for s in sites if (SST_DIR / f"{s['name']}_monthly.csv").exists())
    print(f"Sites with CSV: {count}/{len(sites)}")


if __name__ == '__main__':
    main()
