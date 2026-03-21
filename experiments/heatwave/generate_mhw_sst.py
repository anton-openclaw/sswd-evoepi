#!/usr/bin/env python3
"""Generate synthetic SST scenario CSV files for marine heatwave experiments.

Reads the Blob anomaly template and generates 12 scenario directories,
each with 896 site CSV files covering 2002-2050.

For years < 2025: use original OISST data unchanged.
For years >= 2025: use 2002-2012 climatology + heatwave anomalies at event years.
"""

import json
import os
import sys
import csv
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SST_DIR = os.path.join(REPO, "data", "sst", "site_sst")
SITES_JSON = os.path.join(REPO, "data", "nodes", "all_sites.json")
TEMPLATE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "templates", "blob_anomalies.npz")
OUT_BASE = os.path.join(REPO, "data", "sst", "scenarios")

SIM_START_YEAR = 2012
SIM_END_YEAR = 2050
FUTURE_START = 2025
MAX_SST = 30.0

# Scenario definitions
# Each scenario: dict with 'events' list of (start_year, intensity) or special type
SCENARIOS = {
    "MHW_CTRL": {"type": "control"},
    "MHW_S1_LOW": {"type": "pulse", "events": [(2030, 0.5)]},
    "MHW_S1_MED": {"type": "pulse", "events": [(2030, 1.0)]},
    "MHW_S1_HIGH": {"type": "pulse", "events": [(2030, 1.5)]},
    "MHW_R10_MED": {"type": "pulse", "events": [(2030, 1.0), (2040, 1.0), (2050, 1.0)]},
    "MHW_R7_MED": {"type": "pulse", "events": [(2028, 1.0), (2035, 1.0), (2042, 1.0), (2049, 1.0)]},
    "MHW_R5_MED": {"type": "pulse", "events": [(2028, 1.0), (2033, 1.0), (2038, 1.0), (2043, 1.0), (2048, 1.0)]},
    "MHW_R5_HIGH": {"type": "pulse", "events": [(2028, 1.5), (2033, 1.5), (2038, 1.5), (2043, 1.5), (2048, 1.5)]},
    "MHW_R3_MED": {"type": "pulse", "events": [(2028, 1.0), (2031, 1.0), (2034, 1.0), (2037, 1.0),
                                                  (2040, 1.0), (2043, 1.0), (2046, 1.0), (2049, 1.0)]},
    "MHW_RAMP": {"type": "pulse", "events": [(2028, 0.5), (2033, 0.75), (2038, 1.0), (2043, 1.25), (2048, 1.5)]},
    "MHW_CHRONIC": {"type": "chronic", "offset": 1.0},
    "MHW_COMPOUND": {"type": "compound", "events": [(2030, 1.0, 36), (2033, 1.0, 24)]},
}


def compute_heatwave_anomaly(month_offset, intensity, anom_2014, anom_2015, anom_2016, duration_months=24):
    """Compute the heatwave anomaly for a given month offset from event start.

    Standard 24-month lifecycle:
      Months 1-6: onset (2014 anomalies × 0.5 × intensity)
      Months 7-18: peak (2015 anomalies × 1.0 × intensity)
      Months 19-24: dissipation (2016 anomalies × 0.5 × intensity)

    For compound 36-month full event:
      Months 1-6: onset
      Months 7-30: peak (extended)
      Months 31-36: dissipation

    Args:
        month_offset: 0-based month offset from event start
        intensity: scaling factor
        anom_2014/2015/2016: (12,) arrays of monthly anomalies for this site
        duration_months: 24 (standard) or 36 (extended)

    Returns:
        anomaly value (float), or 0 if outside event window
    """
    if month_offset < 0 or month_offset >= duration_months:
        return 0.0

    # Calendar month (0-11) for this offset
    # The event starts in January of the event year
    cal_month = month_offset % 12

    if duration_months == 36:
        # Extended event: 6 onset, 24 peak, 6 dissipation
        if month_offset < 6:
            return anom_2014[cal_month] * 0.5 * intensity
        elif month_offset < 30:
            return anom_2015[cal_month] * 1.0 * intensity
        else:
            return anom_2016[cal_month] * 0.5 * intensity
    else:
        # Standard 24-month
        if month_offset < 6:
            return anom_2014[cal_month] * 0.5 * intensity
        elif month_offset < 18:
            return anom_2015[cal_month] * 1.0 * intensity
        else:
            return anom_2016[cal_month] * 0.5 * intensity


def generate_site_scenario(site_idx, site_name, scenario_id, scenario_def,
                           climatology, anom_2014, anom_2015, anom_2016,
                           original_data):
    """Generate SST time series for one site under one scenario.

    Returns list of (year, month, sst) tuples.
    """
    clim = climatology[site_idx]
    a14 = anom_2014[site_idx]
    a15 = anom_2015[site_idx]
    a16 = anom_2016[site_idx]

    rows = []

    for year in range(SIM_START_YEAR, SIM_END_YEAR + 1):
        for month in range(1, 13):
            if year < FUTURE_START:
                # Use original data
                sst = original_data.get((year, month), clim[month - 1])
            else:
                # Start with climatology
                sst = clim[month - 1]

                stype = scenario_def["type"]

                if stype == "chronic":
                    sst += scenario_def["offset"]

                elif stype == "pulse":
                    for event_year, intensity in scenario_def["events"]:
                        # Month offset from event start (event starts Jan of event_year)
                        mo_offset = (year - event_year) * 12 + (month - 1)
                        sst += compute_heatwave_anomaly(mo_offset, intensity, a14, a15, a16)

                elif stype == "compound":
                    for event_year, intensity, duration in scenario_def["events"]:
                        mo_offset = (year - event_year) * 12 + (month - 1)
                        sst += compute_heatwave_anomaly(mo_offset, intensity, a14, a15, a16, duration)

                # control: just climatology, no addition

            sst = min(sst, MAX_SST)
            rows.append((year, month, round(sst, 3)))

    return rows


def write_site_csv(out_dir, site_name, rows):
    """Write a site CSV file."""
    path = os.path.join(out_dir, f"{site_name}_monthly.csv")
    with open(path, "w", newline="") as f:
        f.write("year,month,sst\n")
        for yr, mo, sst in rows:
            f.write(f"{yr},{mo},{sst}\n")


def load_original_data(site_name):
    """Load original SST data for a site."""
    path = os.path.join(SST_DIR, f"{site_name}_monthly.csv")
    data = {}
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            data[(int(row["year"]), int(row["month"]))] = float(row["sst"])
    return data


def process_scenario(scenario_id, scenario_def, site_names, climatology,
                     anom_2014, anom_2015, anom_2016, all_original_data):
    """Generate all site CSVs for one scenario."""
    out_dir = os.path.join(OUT_BASE, scenario_id)
    os.makedirs(out_dir, exist_ok=True)

    for i, name in enumerate(site_names):
        rows = generate_site_scenario(
            i, name, scenario_id, scenario_def,
            climatology, anom_2014, anom_2015, anom_2016,
            all_original_data[name]
        )
        write_site_csv(out_dir, name, rows)

    return scenario_id


def main():
    print("Loading blob template...")
    tmpl = np.load(TEMPLATE_PATH, allow_pickle=True)
    site_names = list(tmpl["site_names"])
    climatology = tmpl["climatology"]
    anom_2014 = tmpl["anomaly_2014"]
    anom_2015 = tmpl["anomaly_2015"]
    anom_2016 = tmpl["anomaly_2016"]
    print(f"  {len(site_names)} sites, climatology shape: {climatology.shape}")

    print("Loading original SST data...")
    all_original_data = {}
    for i, name in enumerate(site_names):
        if i % 200 == 0:
            print(f"  Loading site {i}/{len(site_names)}")
        all_original_data[name] = load_original_data(name)

    print(f"\nGenerating {len(SCENARIOS)} scenarios...")
    os.makedirs(OUT_BASE, exist_ok=True)

    for sid, sdef in SCENARIOS.items():
        print(f"  Generating {sid}...")
        process_scenario(sid, sdef, site_names, climatology, anom_2014, anom_2015, anom_2016, all_original_data)
        # Quick validation
        out_dir = os.path.join(OUT_BASE, sid)
        n_files = len([f for f in os.listdir(out_dir) if f.endswith("_monthly.csv")])
        print(f"    -> {n_files} site files written")

    print("\nDone! Scenario directories:")
    for sid in SCENARIOS:
        d = os.path.join(OUT_BASE, sid)
        size_mb = sum(os.path.getsize(os.path.join(d, f)) for f in os.listdir(d)) / 1e6
        print(f"  {sid}: {size_mb:.1f} MB")


if __name__ == "__main__":
    main()
