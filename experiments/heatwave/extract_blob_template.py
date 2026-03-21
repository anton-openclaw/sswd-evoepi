#!/usr/bin/env python3
"""Extract the 2013-2016 Pacific Blob anomaly template from OISST site data.

Computes 2002-2012 monthly climatology per site, then extracts 2014/2015/2016
monthly anomalies. Saves as NPZ for use by generate_mhw_sst.py.
"""

import json
import os
import sys
import csv
import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SST_DIR = os.path.join(REPO, "data", "sst", "site_sst")
SITES_JSON = os.path.join(REPO, "data", "nodes", "all_sites.json")
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "templates")

CLIM_YEARS = range(2002, 2013)  # 2002-2012 inclusive
ANOMALY_YEARS = [2014, 2015, 2016]


def load_site_sst(site_name):
    """Load monthly SST for a site. Returns dict: (year, month) -> sst."""
    path = os.path.join(SST_DIR, f"{site_name}_monthly.csv")
    data = {}
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            yr, mo, sst = int(row["year"]), int(row["month"]), float(row["sst"])
            data[(yr, mo)] = sst
    return data


def main():
    # Load site list
    with open(SITES_JSON) as f:
        sites = json.load(f)

    site_names = [s["name"] for s in sites]
    n_sites = len(site_names)
    print(f"Processing {n_sites} sites...")

    climatology = np.full((n_sites, 12), np.nan)
    anomaly_2014 = np.full((n_sites, 12), np.nan)
    anomaly_2015 = np.full((n_sites, 12), np.nan)
    anomaly_2016 = np.full((n_sites, 12), np.nan)

    for i, name in enumerate(site_names):
        if i % 100 == 0:
            print(f"  Site {i}/{n_sites}: {name}")
        sst = load_site_sst(name)

        # Compute monthly climatology (2002-2012)
        for mo in range(1, 13):
            vals = [sst[(yr, mo)] for yr in CLIM_YEARS if (yr, mo) in sst]
            if vals:
                climatology[i, mo - 1] = np.mean(vals)

        # Extract anomalies for 2014, 2015, 2016
        for yr, arr in [(2014, anomaly_2014), (2015, anomaly_2015), (2016, anomaly_2016)]:
            for mo in range(1, 13):
                if (yr, mo) in sst and not np.isnan(climatology[i, mo - 1]):
                    arr[i, mo - 1] = sst[(yr, mo)] - climatology[i, mo - 1]

    # Verify no NaNs
    for label, arr in [("climatology", climatology), ("anomaly_2014", anomaly_2014),
                       ("anomaly_2015", anomaly_2015), ("anomaly_2016", anomaly_2016)]:
        n_nan = np.isnan(arr).sum()
        if n_nan > 0:
            print(f"  WARNING: {label} has {n_nan} NaN values")
        print(f"  {label}: shape={arr.shape}, range=[{np.nanmin(arr):.3f}, {np.nanmax(arr):.3f}], mean={np.nanmean(arr):.3f}")

    # Save
    os.makedirs(OUT_DIR, exist_ok=True)
    out_path = os.path.join(OUT_DIR, "blob_anomalies.npz")
    np.savez(out_path,
             site_names=np.array(site_names),
             anomaly_2014=anomaly_2014,
             anomaly_2015=anomaly_2015,
             anomaly_2016=anomaly_2016,
             climatology=climatology)
    print(f"\nSaved to {out_path}")
    print(f"  File size: {os.path.getsize(out_path) / 1024:.1f} KB")


if __name__ == "__main__":
    main()
