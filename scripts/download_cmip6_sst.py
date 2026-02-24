#!/usr/bin/env python3
"""
Download CMIP6 SST projections for 11 stepping-stone network nodes.

Uses the Pangeo CMIP6 catalog on Google Cloud Storage (publicly accessible, no API key).
Extracts monthly 'tos' (sea surface temperature) for SSP2-4.5 and SSP5-8.5
from multiple models. Uses xarray .sel(method='nearest') to download only
the grid cells we need (fast, minimal data transfer).

IMPORTANT: For curvilinear grids (gn), we can't use .sel(method='nearest')
directly. We prefer regridded (gr) grids, or crop to a spatial subset first.
"""

import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings('ignore')

import intake
import xarray as xr

# ─── Config ──────────────────────────────────────────────────────────
NODES = {
    'Sitka':         (57.05, -135.33),
    'Ketchikan':     (55.34, -131.64),
    'Haida_Gwaii':   (53.25, -132.07),
    'Bella_Bella':   (52.16, -128.14),
    'Howe_Sound':    (49.38, -123.23),
    'SJI':           (48.53, -123.01),
    'Westport':      (46.89, -124.10),
    'Newport':       (44.63, -124.05),
    'Crescent_City': (41.74, -124.18),
    'Fort_Bragg':    (39.43, -123.80),
    'Monterey':      (36.60, -121.90),
}

TARGET_MODELS = [
    'GFDL-ESM4', 'IPSL-CM6A-LR', 'MRI-ESM2-0', 'CanESM5', 'MIROC6',
    'CESM2', 'UKESM1-0-LL', 'NorESM2-LM', 'EC-Earth3',
]

SCENARIOS = ['ssp245', 'ssp585']
OUTPUT_DIR = Path(__file__).parent.parent / 'data' / 'sst' / 'cmip6'

# NE Pacific bounding box (generous) for spatial subsetting
LAT_MIN, LAT_MAX = 34.0, 60.0
LON_MIN_180, LON_MAX_180 = -140.0, -118.0  # in -180/180
LON_MIN_360, LON_MAX_360 = 220.0, 242.0    # same in 0/360


def extract_regular_grid(ds, variable='tos'):
    """
    Extract SST time series from a regular (gr) grid dataset using .sel().
    Returns (DataFrame, metadata_dict).
    """
    tos = ds[variable]
    
    # Squeeze singleton dims
    for dim in list(tos.dims):
        if dim != 'time' and tos.sizes[dim] == 1:
            tos = tos.isel({dim: 0})
    
    # Determine lon convention
    lon_vals = ds['lon'].values
    use_360 = np.max(lon_vals) > 180
    
    times = []
    node_data = {n: [] for n in NODES}
    metadata = {}
    
    for node_name, (tlat, tlon) in NODES.items():
        tlon_use = (tlon % 360) if use_360 else tlon
        
        # Select nearest grid cell — this downloads only that cell's time series!
        point = tos.sel(lat=tlat, lon=tlon_use, method='nearest')
        values = point.values.astype(float)
        
        # Get actual grid coordinates
        alat = float(point.lat) if 'lat' in point.coords else float(ds['lat'].sel(lat=tlat, method='nearest'))
        alon = float(point.lon) if 'lon' in point.coords else float(ds['lon'].sel(lon=tlon_use, method='nearest'))
        alon_display = alon - 360 if alon > 180 else alon
        
        # Convert Kelvin to Celsius
        if np.nanmean(values) > 100:
            values = values - 273.15
        
        node_data[node_name] = values
        metadata[node_name] = {
            'target_lat': tlat, 'target_lon': tlon,
            'grid_lat': round(alat, 2), 'grid_lon': round(alon_display, 2),
        }
        print(f"      {node_name}: → ({alat:.2f}, {alon_display:.2f}), range: {np.nanmin(values):.1f}–{np.nanmax(values):.1f}°C")
    
    # Get time axis
    time_vals = pd.to_datetime([str(t) for t in tos.time.values])
    
    result = {'time': time_vals}
    result.update(node_data)
    return pd.DataFrame(result), metadata


def extract_curvilinear_grid(ds, variable='tos'):
    """
    Extract SST time series from a curvilinear (gn) grid by:
    1. Loading lat/lon arrays to find nearest ocean cells
    2. Subsetting spatially to NE Pacific first
    3. Extracting point time series
    """
    tos = ds[variable]
    
    # Squeeze singleton dims
    for dim in list(tos.dims):
        if dim != 'time' and tos.sizes[dim] == 1:
            tos = tos.isel({dim: 0})
    
    spatial_dims = [d for d in tos.dims if d != 'time']
    print(f"      Spatial dims: {spatial_dims}")
    
    # Find lat/lon arrays
    lat_name = lon_name = None
    for c in ds.coords:
        cl = c.lower()
        if 'lat' in cl and ds[c].ndim >= 1:
            lat_name = c
        if 'lon' in cl and ds[c].ndim >= 1:
            lon_name = c
    
    if not lat_name or not lon_name:
        raise ValueError(f"Can't find lat/lon in coords: {list(ds.coords)}")
    
    lat_arr = ds[lat_name].values
    lon_arr = ds[lon_name].values
    print(f"      lat ({lat_name}): shape {lat_arr.shape}, lon ({lon_name}): shape {lon_arr.shape}")
    
    # Normalize longitude
    lon_arr_180 = np.where(lon_arr > 180, lon_arr - 360, lon_arr)
    
    # Find NE Pacific mask to narrow our search
    in_region = (
        (lat_arr >= LAT_MIN) & (lat_arr <= LAT_MAX) &
        (lon_arr_180 >= LON_MIN_180) & (lon_arr_180 <= LON_MAX_180)
    )
    
    if not np.any(in_region):
        print(f"      ⚠️ No grid points in NE Pacific region!")
        return None, {}
    
    # Get a single time step for land masking (only NE Pacific region)
    # Find bounding box in index space
    if lat_arr.ndim == 2:
        y_idx, x_idx = np.where(in_region)
    else:
        # 1D lat but curvilinear context (unusual)
        y_idx = np.where(np.any(in_region, axis=1))[0] if in_region.ndim == 2 else np.where(in_region)[0]
        x_idx = np.arange(lon_arr.shape[-1])
    
    y_slice = slice(max(0, y_idx.min() - 2), min(lat_arr.shape[0], y_idx.max() + 3))
    x_slice = slice(max(0, x_idx.min() - 2), min(lat_arr.shape[-1], x_idx.max() + 3))
    
    # Load just one timestep of the subregion
    sub_tos_sample = tos.isel(time=0, **{spatial_dims[0]: y_slice, spatial_dims[1]: x_slice}).values
    sub_lat = lat_arr[y_slice, x_slice] if lat_arr.ndim == 2 else lat_arr[y_slice]
    sub_lon_180 = lon_arr_180[y_slice, x_slice] if lon_arr_180.ndim == 2 else lon_arr_180[x_slice]
    
    print(f"      Subregion: {sub_tos_sample.shape}, lat: {np.nanmin(sub_lat):.1f}–{np.nanmax(sub_lat):.1f}, lon: {np.nanmin(sub_lon_180):.1f}–{np.nanmax(sub_lon_180):.1f}")
    
    metadata = {}
    node_data = {}
    
    for node_name, (tlat, tlon) in NODES.items():
        # Distance in the subregion
        if sub_lat.ndim == 2:
            dist = np.sqrt((sub_lat - tlat)**2 + (sub_lon_180 - tlon)**2)
        else:
            lon_2d, lat_2d = np.meshgrid(sub_lon_180, sub_lat)
            dist = np.sqrt((lat_2d - tlat)**2 + (lon_2d - tlon)**2)
        
        # Mask land
        dist = np.where(np.isnan(sub_tos_sample), np.inf, dist)
        
        if np.all(np.isinf(dist)):
            print(f"      {node_name}: ⚠️ All land in subregion")
            node_data[node_name] = np.full(tos.sizes['time'], np.nan)
            continue
        
        local_idx = np.unravel_index(np.argmin(dist), dist.shape)
        
        # Convert back to global indices
        global_y = y_slice.start + local_idx[0]
        global_x = x_slice.start + local_idx[1]
        
        actual_lat = float(lat_arr[global_y, global_x]) if lat_arr.ndim == 2 else float(lat_arr[global_y])
        actual_lon = float(lon_arr_180[global_y, global_x]) if lon_arr_180.ndim == 2 else float(lon_arr_180[global_x])
        
        # Extract full time series for this single cell
        ts = tos.isel(**{spatial_dims[0]: global_y, spatial_dims[1]: global_x}).values.astype(float)
        
        if np.nanmean(ts) > 100:
            ts = ts - 273.15
        
        node_data[node_name] = ts
        metadata[node_name] = {
            'target_lat': tlat, 'target_lon': tlon,
            'grid_lat': round(actual_lat, 2), 'grid_lon': round(actual_lon, 2),
        }
        print(f"      {node_name}: → ({actual_lat:.2f}, {actual_lon:.2f}), range: {np.nanmin(ts):.1f}–{np.nanmax(ts):.1f}°C")
    
    time_vals = pd.to_datetime([str(t) for t in tos.time.values])
    result = {'time': time_vals}
    result.update(node_data)
    return pd.DataFrame(result), metadata


def process_model(col, model, scenario):
    """Process a single model/scenario combination."""
    print(f"\n--- {model} / {scenario} ---")
    
    # Prefer 'gr' (regridded) grid, fall back to 'gn'
    for grid in ['gr', 'gr1', 'gn']:
        q = col.search(
            variable_id='tos', table_id='Omon',
            experiment_id=scenario, source_id=model,
            member_id='r1i1p1f1', grid_label=grid,
        )
        if len(q.df) > 0:
            break
    
    if len(q.df) == 0:
        # Try any member
        q = col.search(
            variable_id='tos', table_id='Omon',
            experiment_id=scenario, source_id=model,
        )
        if len(q.df) == 0:
            print(f"  No data found")
            return None, {}
    
    grid_used = q.df['grid_label'].iloc[0]
    is_regular = grid_used in ('gr', 'gr1')
    print(f"  Grid: {grid_used} ({'regular' if is_regular else 'curvilinear'})")
    
    try:
        t0 = time.time()
        ddict = q.to_dataset_dict(storage_options={'token': 'anon'})
        key = list(ddict.keys())[0]
        ds = ddict[key]
        print(f"  Loaded in {time.time()-t0:.1f}s. Dims: {dict(ds.sizes)}")
        
        if is_regular:
            df, meta = extract_regular_grid(ds)
        else:
            df, meta = extract_curvilinear_grid(ds)
        
        if df is not None:
            # Filter to 2015-2100
            df = df[(df['time'] >= '2015-01-01') & (df['time'] <= '2100-12-31')].copy()
            print(f"  ✓ {len(df)} months ({df['time'].min().strftime('%Y-%m')} to {df['time'].max().strftime('%Y-%m')})")
        
        return df, meta
    
    except Exception as e:
        print(f"  ERROR: {e}")
        import traceback
        traceback.print_exc()
        return None, {}


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    print("=" * 70)
    print("CMIP6 SST Projection Download")
    print("Pangeo CMIP6 catalog → Google Cloud Storage (no API key)")
    print("=" * 70)
    
    # Open catalog
    print("\nLoading catalog...")
    catalog_url = "https://storage.googleapis.com/cmip6/pangeo-cmip6.json"
    col = intake.open_esm_datastore(catalog_url)
    print(f"Catalog: {len(col.df)} assets\n")
    
    # Check availability
    print("Model availability:")
    available = {}
    for model in TARGET_MODELS:
        ok = []
        for scen in SCENARIOS:
            q = col.search(variable_id='tos', table_id='Omon', experiment_id=scen, source_id=model)
            if len(q.df) > 0:
                ok.append(scen)
        if set(SCENARIOS).issubset(set(ok)):
            available[model] = ok
            print(f"  ✓ {model}")
        else:
            print(f"  ✗ {model} (has: {ok})")
    
    print(f"\n{len(available)} models with both scenarios")
    
    if not available:
        print("ERROR: No models found!")
        return 1
    
    # Process each model/scenario
    all_data = {s: {} for s in SCENARIOS}
    all_meta = {s: {} for s in SCENARIOS}
    
    for scenario in SCENARIOS:
        print(f"\n{'='*70}")
        print(f"SCENARIO: {scenario}")
        print(f"{'='*70}")
        
        for model in available:
            df, meta = process_model(col, model, scenario)
            if df is not None and len(df) > 0:
                all_data[scenario][model] = df
                all_meta[scenario][model] = meta
                
                fname = f"{model}_{scenario}_tos_monthly.csv"
                df.to_csv(OUTPUT_DIR / fname, index=False, float_format='%.4f')
                print(f"  Saved: {fname}")
            
            time.sleep(0.5)
    
    # Ensemble means
    print(f"\n{'='*70}")
    print("ENSEMBLE MEANS")
    print(f"{'='*70}")
    
    ensemble_info = {}
    
    for scenario in SCENARIOS:
        models_ok = list(all_data[scenario].keys())
        if not models_ok:
            continue
        
        print(f"\n{scenario}: {len(models_ok)} models")
        
        target_times = pd.date_range('2026-01-01', '2100-12-01', freq='MS')
        node_stacks = {n: [] for n in NODES}
        models_used = []
        
        for model in models_ok:
            df = all_data[scenario][model].set_index('time')
            usable = True
            temp_stacks = {}
            
            for node in NODES:
                if node not in df.columns:
                    usable = False
                    break
                ts = df[node].reindex(target_times, method='nearest', tolerance=pd.Timedelta('20D'))
                if ts.notna().mean() < 0.5:
                    usable = False
                    break
                temp_stacks[node] = ts.values
            
            if usable:
                models_used.append(model)
                for node in NODES:
                    node_stacks[node].append(temp_stacks[node])
        
        if not models_used:
            print(f"  No models with sufficient coverage")
            continue
        
        print(f"  Using {len(models_used)} models: {models_used}")
        
        # Mean
        ens_df = pd.DataFrame({'time': target_times})
        ens_std_df = pd.DataFrame({'time': target_times})
        
        for node in NODES:
            stacked = np.array(node_stacks[node])
            ens_df[node] = np.nanmean(stacked, axis=0)
            ens_std_df[node] = np.nanstd(stacked, axis=0)
        
        ens_df.to_csv(OUTPUT_DIR / f'ensemble_mean_{scenario}_tos_monthly.csv', index=False, float_format='%.4f')
        ens_std_df.to_csv(OUTPUT_DIR / f'ensemble_std_{scenario}_tos_monthly.csv', index=False, float_format='%.4f')
        
        for node in NODES:
            v = ens_df[node].dropna()
            if len(v) > 0:
                print(f"    {node}: {v.min():.1f}–{v.max():.1f}°C (μ={v.mean():.1f})")
        
        ensemble_info[scenario] = {
            'n_models': len(models_used),
            'models': models_used,
            'time_range': '2026-01 to 2100-12',
            'n_months': len(target_times),
        }
    
    # Metadata
    meta_out = {
        'source': 'Pangeo CMIP6 catalog (Google Cloud Storage)',
        'catalog_url': catalog_url,
        'variable': 'tos (sea surface temperature)',
        'frequency': 'monthly (Omon table)',
        'units': 'degrees Celsius',
        'scenarios': SCENARIOS,
        'models_queried': TARGET_MODELS,
        'models_with_data': {s: list(all_data[s].keys()) for s in SCENARIOS},
        'ensemble_info': ensemble_info,
        'grid_point_metadata': {
            s: {m: all_meta[s][m] for m in all_meta[s]} for s in SCENARIOS
        },
        'node_coordinates': {n: {'lat': c[0], 'lon': c[1]} for n, c in NODES.items()},
        'download_date': pd.Timestamp.now().isoformat(),
        'notes': [
            'SSP scenarios run from 2015-2100; ensemble means cover 2026-2100',
            'Each model: nearest ocean grid point to target coordinates',
            'Ensemble std files = inter-model spread (uncertainty)',
            'All temperatures in Celsius (Kelvin-origin data converted)',
            'Preferred regridded (gr) grids; used curvilinear (gn) where gr unavailable',
        ],
    }
    
    with open(OUTPUT_DIR / 'download_metadata.json', 'w') as f:
        json.dump(meta_out, f, indent=2, default=str)
    
    print(f"\n{'='*70}")
    print(f"COMPLETE — {OUTPUT_DIR}")
    n_files = sum(len(v) for v in all_data.values())
    print(f"  {n_files} model CSVs + {len(ensemble_info)} ensemble means")
    print(f"{'='*70}")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
