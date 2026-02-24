#!/usr/bin/env python3
"""
Bulk-fetch OISST data for all sites using xarray + OPeNDAP.

Instead of making hundreds of ERDDAP requests, this:
1. Opens the OISST dataset remotely via OPeNDAP
2. Extracts data for all unique grid cells in one efficient operation
3. Computes monthly means and climatology locally

Much faster than individual ERDDAP CSV requests.

Usage:
    python3 scripts/fetch_sst_bulk.py
"""

import json
import csv
import sys
import time
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import numpy as np
import xarray as xr

BASE = Path(__file__).parent.parent
SITES_FILE = BASE / "data" / "nodes" / "all_sites.json"
OUTPUT_DIR = BASE / "data" / "sst" / "site_sst"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# OISST OPeNDAP URL
OPENDAP_URL = "https://coastwatch.pfeg.noaa.gov/erddap/griddap/ncdcOisst21Agg_LonPM180"

GRID_RES = 0.25
START_YEAR = 2002
END_YEAR = 2025


def snap_to_grid(val, res=GRID_RES):
    return round(round(val / res) * res, 4)


def load_sites():
    with open(SITES_FILE) as f:
        sites = json.load(f)
    return sites


def get_unique_grid_cells(sites):
    """Get unique OISST grid cells and site mapping."""
    site_to_grid = {}
    grid_cells = {}  # (glat, glon) → list of site indices

    for i, s in enumerate(sites):
        lat = s.get('latitude', s.get('lat'))
        lon = s.get('longitude', s.get('lon'))
        if lon > 0:
            lon -= 360.0
        
        glat = snap_to_grid(lat)
        glon = snap_to_grid(lon)
        
        key = (glat, glon)
        site_to_grid[i] = key
        
        if key not in grid_cells:
            grid_cells[key] = []
        grid_cells[key].append(i)
    
    return site_to_grid, grid_cells


def fetch_sst_for_cells(grid_cells, start_year=START_YEAR, end_year=END_YEAR):
    """
    Fetch SST for all grid cells using xarray OPeNDAP.
    
    Strategy: batch cells by latitude bands to avoid loading the entire dataset.
    """
    lats = sorted(set(glat for glat, glon in grid_cells.keys()))
    lons = sorted(set(glon for glat, glon in grid_cells.keys()))
    
    print(f"  Lat range: {min(lats):.2f} to {max(lats):.2f} ({len(lats)} unique)")
    print(f"  Lon range: {min(lons):.2f} to {max(lons):.2f} ({len(lons)} unique)")
    
    # Open dataset
    print(f"\n  Opening OPeNDAP dataset...", flush=True)
    ds = xr.open_dataset(OPENDAP_URL)
    print(f"  Dataset opened. Variables: {list(ds.data_vars)}", flush=True)
    
    # Select our time range
    time_sel = slice(f"{start_year}-01-01", f"{end_year}-12-31")
    
    results = {}
    
    # Process in latitude bands (10° chunks) to manage memory
    lat_bands = []
    band_start = min(lats)
    while band_start <= max(lats):
        band_end = min(band_start + 10.0, max(lats) + 0.1)
        band_lats = [l for l in lats if band_start <= l < band_end]
        band_lons_for_band = sorted(set(
            glon for glat, glon in grid_cells.keys() 
            if glat in band_lats
        ))
        if band_lats:
            lat_bands.append((band_start, band_end, band_lats, band_lons_for_band))
        band_start = band_end
    
    print(f"\n  Processing {len(lat_bands)} latitude bands...", flush=True)
    
    for band_idx, (lat_start, lat_end, band_lats, band_lons) in enumerate(lat_bands):
        t0 = time.time()
        print(f"\n  Band {band_idx+1}/{len(lat_bands)}: "
              f"lat {lat_start:.1f}–{lat_end:.1f}, "
              f"{len(band_lats)} lats × {len(band_lons)} lons", flush=True)
        
        # Select this band's data
        try:
            # Use .sel with method='nearest' to find the actual grid cells
            sst_band = ds['sst'].sel(
                time=time_sel,
                zlev=0.0,
                latitude=band_lats,
                longitude=band_lons,
                method='nearest'
            ).load()
            
            t1 = time.time()
            print(f"    Loaded in {t1-t0:.1f}s, shape: {sst_band.shape}", flush=True)
            
            # Extract monthly means for each cell in this band
            for glat, glon in grid_cells.keys():
                if glat not in band_lats:
                    continue
                
                try:
                    cell_data = sst_band.sel(
                        latitude=glat, 
                        longitude=glon,
                        method='nearest'
                    )
                    
                    # Compute monthly climatology
                    monthly = cell_data.groupby('time.month').mean(skipna=True)
                    monthly_std = cell_data.groupby('time.month').std(skipna=True)
                    
                    clim = {}
                    for month in range(1, 13):
                        try:
                            mean_val = float(monthly.sel(month=month).values)
                            std_val = float(monthly_std.sel(month=month).values)
                            clim[month] = {
                                'mean': round(mean_val, 3) if not np.isnan(mean_val) else None,
                                'std': round(std_val, 3) if not np.isnan(std_val) else None,
                            }
                        except:
                            clim[month] = {'mean': None, 'std': None}
                    
                    # Also get yearly monthly means for time series
                    yearly_monthly = cell_data.groupby('time.year').apply(
                        lambda x: x.groupby('time.month').mean(skipna=True)
                    )
                    
                    monthly_ts = []
                    for year in range(start_year, end_year + 1):
                        for month in range(1, 13):
                            try:
                                val = float(yearly_monthly.sel(year=year, month=month).values)
                                if not np.isnan(val):
                                    monthly_ts.append({
                                        'year': year, 'month': month,
                                        'sst_mean': round(val, 3)
                                    })
                            except:
                                pass
                    
                    results[(glat, glon)] = {
                        'lat': glat, 'lon': glon,
                        'n_sites': len(grid_cells[(glat, glon)]),
                        'climatology': clim,
                        'monthly_ts': monthly_ts,
                        'n_valid_months': len(monthly_ts),
                    }
                    
                except Exception as e:
                    print(f"    Warning: failed for ({glat}, {glon}): {e}", flush=True)
            
            cells_done = sum(1 for k in results if k[0] in band_lats)
            cells_total = sum(1 for glat, _ in grid_cells if glat in band_lats)
            print(f"    Extracted {cells_done}/{cells_total} cells from this band", flush=True)
            
        except Exception as e:
            print(f"    ERROR loading band: {e}", flush=True)
    
    ds.close()
    return results


def main():
    print("Loading sites...", flush=True)
    sites = load_sites()
    site_to_grid, grid_cells = get_unique_grid_cells(sites)
    
    print(f"  {len(sites)} sites → {len(grid_cells)} unique grid cells", flush=True)
    
    # Save mapping
    mapping = {str(k): f"{v[0]},{v[1]}" for k, v in site_to_grid.items()}
    with open(OUTPUT_DIR / "site_to_grid.json", 'w') as f:
        json.dump(mapping, f, indent=2)
    
    # Fetch
    print("\nFetching SST data via OPeNDAP...", flush=True)
    t0 = time.time()
    results = fetch_sst_for_cells(grid_cells)
    elapsed = time.time() - t0
    
    print(f"\n  Fetched {len(results)}/{len(grid_cells)} cells in {elapsed/60:.1f} min")
    
    # Save results
    serializable = {}
    for (glat, glon), data in results.items():
        key = f"{glat},{glon}"
        serializable[key] = data
    
    with open(OUTPUT_DIR / "grid_cell_data.json", 'w') as f:
        json.dump(serializable, f)
    print(f"  Saved grid cell data")
    
    # Generate climatology CSV
    clim_file = OUTPUT_DIR / "all_sites_climatology.csv"
    with open(clim_file, 'w', newline='') as f:
        writer = csv.writer(f)
        header = ['site_index', 'site_name', 'region', 'grid_lat', 'grid_lon']
        for m in range(1, 13):
            header.extend([f'month_{m:02d}_mean', f'month_{m:02d}_std'])
        writer.writerow(header)
        
        for i, site in enumerate(sites):
            gkey = site_to_grid[i]
            gdata = results.get(gkey, {})
            clim = gdata.get('climatology', {})
            
            row = [i, site.get('name', '?'), site.get('region', '?'),
                   gkey[0], gkey[1]]
            for m in range(1, 13):
                mc = clim.get(m, {})
                row.append(mc.get('mean', ''))
                row.append(mc.get('std', ''))
            writer.writerow(row)
    
    # Check coverage
    nan_sites = sum(1 for i in range(len(sites)) 
                    if site_to_grid[i] not in results or
                    all(results[site_to_grid[i]]['climatology'].get(m, {}).get('mean') is None 
                        for m in range(1, 13)))
    
    print(f"\n✓ Complete:")
    print(f"  {len(results)} grid cells with data")
    print(f"  {len(sites) - nan_sites}/{len(sites)} sites with valid SST")
    print(f"  {nan_sites} sites with no SST (land-masked or data gap)")
    print(f"  Output: {clim_file}")
    
    # Summary
    summary = {
        'timestamp': datetime.utcnow().isoformat(),
        'n_sites': len(sites),
        'n_grid_cells': len(grid_cells),
        'n_with_data': len(results),
        'n_nan_sites': nan_sites,
        'elapsed_min': round(elapsed / 60, 1),
        'source': 'NOAA OISST v2.1 via OPeNDAP',
        'period': f'{START_YEAR}-{END_YEAR}',
    }
    with open(OUTPUT_DIR / "fetch_summary.json", 'w') as f:
        json.dump(summary, f, indent=2)


if __name__ == '__main__':
    main()
