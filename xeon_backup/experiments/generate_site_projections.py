#!/usr/bin/env python3
"""Generate SST projection files for all model sites using delta method.

For each non-reference site:
  1. Find the nearest reference site (from 11 CMIP6 projection sites)
  2. Compute the site's observed monthly climatology (2021-2025 mean)
  3. Compute the reference site's observed monthly climatology (same period)
  4. Apply delta method:
     site_proj(yr, mo) = site_clim(mo) + [ref_proj(yr, mo) - ref_clim(mo)]
  5. Write as {site_name}_{scenario}_monthly.csv

This preserves each site's base temperature while inheriting regional
warming trends from the nearest CMIP6 reference site.
"""

import csv
import json
import math
import os
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SITES_PATH = PROJECT_ROOT / 'data' / 'nodes' / 'all_sites.json'
OBS_DIR = PROJECT_ROOT / 'data' / 'sst' / 'site_sst'
PROJ_DIR = PROJECT_ROOT / 'data' / 'sst' / 'projections'
CLIM_YEARS = range(2021, 2026)  # 5-year baseline for climatology

# Reference sites with approximate coordinates
REFERENCE_SITES = {
    "Ketchikan": (55.34, -131.64),
    "Sitka": (57.05, -135.33),
    "Haida_Gwaii": (53.25, -132.07),
    "Bella_Bella": (52.16, -128.14),
    "Howe_Sound": (49.37, -123.23),
    "SJI": (48.53, -123.01),
    "Westport": (46.89, -124.10),
    "Newport": (44.63, -124.05),
    "Crescent_City": (41.75, -124.20),
    "Fort_Bragg": (39.44, -123.80),
    "Monterey": (36.62, -121.90),
}

SCENARIOS = ['ssp245', 'ssp585']


def haversine(lat1, lon1, lat2, lon2):
    """Haversine distance in km."""
    R = 6371
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = (math.sin(dlat / 2) ** 2 +
         math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) *
         math.sin(dlon / 2) ** 2)
    return R * 2 * math.asin(math.sqrt(a))


def normalize_name(name):
    """Normalize node name for file lookup (spaces → underscores only)."""
    return name.strip().replace(' ', '_')


def load_monthly_obs(site_name, data_dir):
    """Load monthly SST observations for a site. Returns {(yr,mo): sst}."""
    norm = normalize_name(site_name)
    path = os.path.join(data_dir, f"{norm}_monthly.csv")
    if not os.path.isfile(path):
        return None
    data = {}
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            data[(int(row['year']), int(row['month']))] = float(row['sst'])
    return data


def compute_monthly_climatology(monthly_data, years=CLIM_YEARS):
    """Compute mean SST for each month over the given years."""
    clim = {}
    for mo in range(1, 13):
        vals = [monthly_data.get((yr, mo)) for yr in years
                if (yr, mo) in monthly_data]
        if vals:
            clim[mo] = sum(vals) / len(vals)
        else:
            # Fall back to all available data
            all_vals = [v for (y, m), v in monthly_data.items() if m == mo]
            clim[mo] = sum(all_vals) / len(all_vals) if all_vals else 10.0
    return clim


def load_projection(ref_name, scenario, proj_dir):
    """Load reference site projection data. Returns {(yr,mo): sst}."""
    # Try exact then placeholder
    for suffix in [f"_{scenario}_monthly.csv",
                   f"_{scenario}_placeholder_monthly.csv"]:
        path = os.path.join(proj_dir, ref_name + suffix)
        if os.path.isfile(path):
            data = {}
            with open(path) as f:
                reader = csv.DictReader(f)
                for row in reader:
                    data[(int(row['year']), int(row['month']))] = float(row['sst'])
            return data
    raise FileNotFoundError(f"No projection for {ref_name}/{scenario}")


def nearest_reference(lat, lon):
    """Find nearest reference site. Returns (name, distance_km)."""
    best = None
    best_dist = float('inf')
    for name, (rlat, rlon) in REFERENCE_SITES.items():
        d = haversine(lat, lon, rlat, rlon)
        if d < best_dist:
            best_dist = d
            best = name
    return best, best_dist


def main():
    # Load all model sites
    with open(SITES_PATH) as f:
        sites = json.load(f)
    print(f"Loaded {len(sites)} sites")

    # Pre-load reference site observations and projections
    ref_obs = {}
    ref_clim = {}
    ref_proj = {s: {} for s in SCENARIOS}

    for ref_name in REFERENCE_SITES:
        obs = load_monthly_obs(ref_name, str(OBS_DIR))
        if obs is None:
            # Reference sites might use different naming in obs
            print(f"WARNING: No obs for reference site {ref_name}, "
                  f"using projection baseline")
            ref_obs[ref_name] = None
        else:
            ref_obs[ref_name] = obs
            ref_clim[ref_name] = compute_monthly_climatology(obs)

        for scenario in SCENARIOS:
            ref_proj[scenario][ref_name] = load_projection(
                ref_name, scenario, str(PROJ_DIR))

    # For reference sites without obs, compute clim from first years of projection
    for ref_name in REFERENCE_SITES:
        if ref_name not in ref_clim:
            # Use mean of first 5 projection years as baseline
            proj = ref_proj['ssp245'][ref_name]
            clim = {}
            for mo in range(1, 13):
                vals = [proj.get((yr, mo), 10.0) for yr in range(2026, 2031)]
                clim[mo] = sum(vals) / len(vals)
            ref_clim[ref_name] = clim
            print(f"  Using projection baseline for {ref_name}")

    # Generate projection files for each site
    generated = 0
    skipped = 0
    errors = 0

    for site in sites:
        site_name = site['name']
        site_lat = site['latitude']
        site_lon = site['longitude']

        # Load site observations
        site_obs = load_monthly_obs(site_name, str(OBS_DIR))
        if site_obs is None:
            print(f"WARNING: No obs for site {site_name}, skipping")
            errors += 1
            continue

        site_clim = compute_monthly_climatology(site_obs)

        # Find nearest reference site
        ref_name, ref_dist = nearest_reference(site_lat, site_lon)

        for scenario in SCENARIOS:
            norm_name = normalize_name(site_name)
            out_path = PROJ_DIR / f"{norm_name}_{scenario}_monthly.csv"

            # Skip if already a real projection file (reference site)
            if out_path.exists() and site_name in [
                n.replace('_', '-') for n in REFERENCE_SITES]:
                skipped += 1
                continue

            # Delta method: site_proj = site_clim + (ref_proj - ref_clim)
            proj_data = ref_proj[scenario][ref_name]
            rc = ref_clim[ref_name]

            rows = []
            for (yr, mo), ref_sst in sorted(proj_data.items()):
                delta = ref_sst - rc.get(mo, ref_sst)
                site_sst = site_clim.get(mo, 10.0) + delta
                rows.append((yr, mo, round(site_sst, 4)))

            with open(out_path, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['year', 'month', 'sst'])
                for yr, mo, sst in rows:
                    writer.writerow([yr, mo, sst])

            generated += 1

    print(f"\nDone: {generated} files generated, {skipped} skipped (existing), "
          f"{errors} errors")
    print(f"Output directory: {PROJ_DIR}")


if __name__ == '__main__':
    main()
