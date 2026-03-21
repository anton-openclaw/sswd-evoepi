#!/usr/bin/env python3
"""
Fetch NOAA OISST v2.1 data for all sites in the network.

Strategy:
- OISST is on a 0.25° grid → many sites share grid cells
- Group sites by their nearest OISST grid cell
- Fetch SST for each unique grid cell from ERDDAP (once)
- Save per-site SST by mapping site → grid cell

Outputs:
  data/sst/site_sst/grid_cell_monthly.json     — all unique grid cells with monthly SST
  data/sst/site_sst/site_to_grid.json           — mapping: site_index → grid_cell_key
  data/sst/site_sst/all_sites_climatology.csv   — 365-day climatology for all sites
  data/sst/site_sst/fetch_summary.json          — metadata

Usage:
    python3 scripts/fetch_sst_all_sites.py [--max-concurrent 5]
"""

import json
import csv
import io
import sys
import time
import urllib.request
import urllib.error
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

import numpy as np

BASE = Path(__file__).parent.parent
SITES_FILE = BASE / "data" / "nodes" / "all_sites.json"
OUTPUT_DIR = BASE / "data" / "sst" / "site_sst"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ERDDAP configuration
ERDDAP_BASE = "https://coastwatch.pfeg.noaa.gov/erddap"
DATASET_ID = "ncdcOisst21Agg_LonPM180"
START_DATE = "2002-01-01"
END_DATE = "2025-12-31"

# OISST grid resolution
GRID_RES = 0.25


def snap_to_grid(lat, lon, res=GRID_RES):
    """Snap lat/lon to nearest OISST grid center."""
    grid_lat = round(round(lat / res) * res, 4)
    grid_lon = round(round(lon / res) * res, 4)
    return grid_lat, grid_lon


def grid_key(lat, lon):
    """Create a string key for a grid cell."""
    return f"{lat:.4f},{lon:.4f}"


def load_and_group_sites():
    """Load all sites and group by OISST grid cell."""
    with open(SITES_FILE) as f:
        sites = json.load(f)
    
    site_to_grid = {}  # site_index → grid_key
    grid_cells = {}    # grid_key → (grid_lat, grid_lon, [site_indices])
    
    for i, s in enumerate(sites):
        lat = s.get('latitude', s.get('lat'))
        lon = s.get('longitude', s.get('lon'))
        
        # Normalize positive longitudes
        if lon > 0:
            lon = lon - 360.0
        
        glat, glon = snap_to_grid(lat, lon)
        key = grid_key(glat, glon)
        
        site_to_grid[str(i)] = key
        
        if key not in grid_cells:
            grid_cells[key] = {
                'lat': glat,
                'lon': glon,
                'sites': []
            }
        grid_cells[key]['sites'].append(i)
    
    return sites, site_to_grid, grid_cells


def fetch_sst_by_year(lat, lon, year, retries=3):
    """
    Fetch daily SST for a grid cell for ONE year from ERDDAP.
    Returns list of records or None on failure.
    """
    start = f"{year}-01-01T12:00:00Z"
    end = f"{year}-12-31T12:00:00Z"
    
    url = (
        f"{ERDDAP_BASE}/griddap/{DATASET_ID}.csv?"
        f"sst[({start}):1:({end})]"
        f"[(0.0):1:(0.0)]"
        f"[({lat}):1:({lat})]"
        f"[({lon}):1:({lon})]"
    )
    
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url)
            req.add_header('User-Agent', 'SSWD-EvoEpi-Research/1.0')
            with urllib.request.urlopen(req, timeout=60) as response:
                text = response.read().decode('utf-8')
            
            reader = csv.reader(io.StringIO(text))
            header = next(reader)
            next(reader)  # units row
            
            sst_col = time_col = None
            for j, h in enumerate(header):
                if 'sst' in h.lower():
                    sst_col = j
                if 'time' in h.lower():
                    time_col = j
            
            if sst_col is None:
                return None
            
            records = []
            for row in reader:
                try:
                    time_str = row[time_col]
                    sst_val = float(row[sst_col])
                    dt = datetime.fromisoformat(time_str.replace('Z', '+00:00'))
                    records.append({
                        'year': dt.year,
                        'month': dt.month,
                        'day': dt.day,
                        'sst': sst_val
                    })
                except (ValueError, IndexError):
                    continue
            
            return records
            
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
            if attempt < retries - 1:
                time.sleep(2 ** (attempt + 1))
            else:
                return None
    
    return None


def fetch_monthly_sst(lat, lon, retries=3):
    """
    Fetch daily SST for a grid cell, year by year, then aggregate to monthly.
    Returns list of daily records or None if all years fail.
    """
    all_records = []
    years = list(range(int(START_DATE[:4]), int(END_DATE[:4]) + 1))
    
    for year in years:
        records = fetch_sst_by_year(lat, lon, year, retries=retries)
        if records:
            all_records.extend(records)
    
    if not all_records:
        return None
    
    return all_records


def compute_monthly_means(daily_records):
    """Aggregate daily records into monthly means."""
    monthly = defaultdict(list)
    for rec in daily_records:
        key = (rec['year'], rec['month'])
        if not np.isnan(rec['sst']):
            monthly[key].append(rec['sst'])
    
    result = []
    for (year, month), vals in sorted(monthly.items()):
        result.append({
            'year': year,
            'month': month,
            'sst_mean': round(np.mean(vals), 3),
            'sst_std': round(np.std(vals), 3),
            'n_days': len(vals)
        })
    
    return result


def compute_climatology(monthly_data):
    """Compute 12-month climatology from monthly means."""
    by_month = defaultdict(list)
    for rec in monthly_data:
        by_month[rec['month']].append(rec['sst_mean'])
    
    clim = {}
    for month in range(1, 13):
        vals = by_month.get(month, [])
        if vals:
            clim[month] = {
                'mean': round(np.mean(vals), 3),
                'std': round(np.std(vals), 3),
                'n_years': len(vals)
            }
        else:
            clim[month] = {'mean': None, 'std': None, 'n_years': 0}
    
    return clim


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--max-concurrent', type=int, default=2,
                        help='Max concurrent ERDDAP requests (ERDDAP rate-limits at >2)')
    parser.add_argument('--resume', action='store_true',
                        help='Skip already-fetched grid cells')
    args = parser.parse_args()
    
    print("Loading and grouping sites...")
    sites, site_to_grid, grid_cells = load_and_group_sites()
    
    print(f"  {len(sites)} sites → {len(grid_cells)} unique OISST grid cells")
    print(f"  Grid cells per site: min={min(len(g['sites']) for g in grid_cells.values())}, "
          f"max={max(len(g['sites']) for g in grid_cells.values())}, "
          f"mean={np.mean([len(g['sites']) for g in grid_cells.values()]):.1f}")
    
    # Save site-to-grid mapping
    mapping_file = OUTPUT_DIR / "site_to_grid.json"
    with open(mapping_file, 'w') as f:
        json.dump(site_to_grid, f, indent=2)
    print(f"  Saved site→grid mapping to {mapping_file.name}")
    
    # Load existing results if resuming
    results_file = OUTPUT_DIR / "grid_cell_monthly.json"
    if args.resume and results_file.exists():
        with open(results_file) as f:
            grid_results = json.load(f)
        print(f"  Resuming: {len(grid_results)} grid cells already fetched")
    else:
        grid_results = {}
    
    # Fetch SST for each unique grid cell
    keys_to_fetch = [k for k in grid_cells if k not in grid_results]
    print(f"\nFetching SST for {len(keys_to_fetch)} grid cells "
          f"({len(grid_cells) - len(keys_to_fetch)} already done)...")
    
    completed = 0
    failed = 0
    nan_cells = 0
    
    def fetch_one(key):
        cell = grid_cells[key]
        daily = fetch_monthly_sst(cell['lat'], cell['lon'])
        if daily is None:
            return key, None
        monthly = compute_monthly_means(daily)
        clim = compute_climatology(monthly)
        return key, {
            'lat': cell['lat'],
            'lon': cell['lon'],
            'n_sites': len(cell['sites']),
            'monthly': monthly,
            'climatology': clim,
        }
    
    start = time.time()
    save_interval = 20  # save progress every N completions
    
    with ThreadPoolExecutor(max_workers=args.max_concurrent) as pool:
        futures = {pool.submit(fetch_one, k): k for k in keys_to_fetch}
        
        for future in as_completed(futures):
            key, result = future.result()
            completed += 1
            
            if result is None:
                failed += 1
                # Try with slightly adjusted coordinates
                cell = grid_cells[key]
                print(f"  FAILED: {key} ({cell['lat']}, {cell['lon']}) "
                      f"- {len(cell['sites'])} sites affected")
            else:
                grid_results[key] = result
                # Check for NaN months
                if any(m['sst_mean'] is None or np.isnan(m['sst_mean']) 
                       for m in result.get('monthly', [])):
                    nan_cells += 1
            
            if completed % 10 == 0:
                elapsed = time.time() - start
                rate = completed / elapsed if elapsed > 0 else 0
                remaining = len(keys_to_fetch) - completed
                eta = remaining / rate if rate > 0 else 0
                print(f"  {completed}/{len(keys_to_fetch)} "
                      f"({failed} failed, {nan_cells} NaN) "
                      f"- {rate:.1f}/s, ETA {eta/60:.0f}m")
            
            # Save periodically
            if completed % save_interval == 0:
                with open(results_file, 'w') as f:
                    json.dump(grid_results, f)
    
    # Final save
    with open(results_file, 'w') as f:
        json.dump(grid_results, f)
    
    elapsed = time.time() - start
    print(f"\nDone: {completed} cells in {elapsed/60:.1f} min, "
          f"{failed} failed, {nan_cells} with NaN months")
    
    # Generate climatology CSV for all sites
    print("\nGenerating all-sites climatology CSV...")
    clim_file = OUTPUT_DIR / "all_sites_climatology.csv"
    with open(clim_file, 'w', newline='') as f:
        writer = csv.writer(f)
        header = ['site_index', 'site_name', 'region', 'grid_lat', 'grid_lon']
        for m in range(1, 13):
            header.extend([f'month_{m:02d}_mean', f'month_{m:02d}_std'])
        writer.writerow(header)
        
        for i, site in enumerate(sites):
            gkey = site_to_grid.get(str(i))
            gdata = grid_results.get(gkey, {})
            clim = gdata.get('climatology', {})
            
            row = [i, site.get('name', '?'), site.get('region', '?'),
                   gdata.get('lat', ''), gdata.get('lon', '')]
            for m in range(1, 13):
                mc = clim.get(str(m), clim.get(m, {}))
                row.append(mc.get('mean', ''))
                row.append(mc.get('std', ''))
            writer.writerow(row)
    
    print(f"  Saved {clim_file.name}")
    
    # Summary
    summary = {
        'timestamp': datetime.utcnow().isoformat(),
        'n_sites': len(sites),
        'n_grid_cells': len(grid_cells),
        'n_fetched': len(grid_results),
        'n_failed': failed,
        'n_nan_cells': nan_cells,
        'elapsed_min': round(elapsed / 60, 1),
        'source': 'NOAA OISST v2.1 via ERDDAP',
        'dataset': DATASET_ID,
        'period': f'{START_DATE} to {END_DATE}',
    }
    with open(OUTPUT_DIR / "fetch_summary.json", 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"\n✓ Complete: {len(grid_results)}/{len(grid_cells)} grid cells, "
          f"covering {sum(len(g['sites']) for g in grid_results.values() if isinstance(g, dict) and 'sites' in g)}"
          f"/{len(sites)} sites")


if __name__ == '__main__':
    main()
