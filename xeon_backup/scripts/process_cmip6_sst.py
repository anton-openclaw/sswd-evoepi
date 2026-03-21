#!/usr/bin/env python3
"""
Process CMIP6 SST data into model-ready per-node monthly CSVs.

Takes raw CMIP6 CSV files (from download_cmip6_sst.py or manual download)
and produces bias-corrected, per-node, per-scenario monthly CSV files
compatible with the SSWD-EvoEpi model.

Workflow:
  1. Load OISST observational monthly data (2002-2025) as ground truth
  2. Load CMIP6 model outputs (already extracted to node locations)
  3. Compute per-node, per-month bias correction (CMIP6 historical - OISST)
     over the overlap period (2015-2025)
  4. Multi-model ensemble mean (if multiple models available)
  5. Apply bias correction to projections
  6. Output per-node, per-scenario CSVs: {Node}_{scenario}_monthly.csv

Usage:
  python scripts/process_cmip6_sst.py [--cmip6-dir data/sst/cmip6]
                                      [--oisst-dir data/sst]
                                      [--output-dir data/sst/projections]
                                      [--scenarios ssp245,ssp585]

Output format (matches OISST *_monthly.csv schema):
  year,month,sst
  2026,1,7.23
  2026,2,6.01
  ...
"""

import argparse
import csv
import os
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import numpy as np


# ─── Node definitions ────────────────────────────────────────────────
NODES = [
    'Sitka', 'Ketchikan', 'Haida_Gwaii', 'Bella_Bella', 'Howe_Sound',
    'SJI', 'Westport', 'Newport', 'Crescent_City', 'Fort_Bragg', 'Monterey',
]

# Overlap period for bias correction
OVERLAP_START = 2015
OVERLAP_END = 2025  # inclusive


def load_oisst_monthly(node_name: str, data_dir: str) -> dict:
    """Load OISST monthly data for a node.

    Returns dict: (year, month) → SST (°C).
    """
    csv_path = os.path.join(data_dir, f"{node_name}_monthly.csv")
    if not os.path.isfile(csv_path):
        raise FileNotFoundError(f"OISST file not found: {csv_path}")

    monthly = {}
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            yr = int(row['year'])
            mo = int(row['month'])
            monthly[(yr, mo)] = float(row['sst'])
    return monthly


def load_cmip6_csv(csv_path: str) -> dict:
    """Load a CMIP6 CSV file (model_scenario_tos_monthly.csv).

    Expected format: time,Sitka,Ketchikan,...,Monterey
    Time format: 'YYYY-MM-DD HH:MM:SS'

    Returns dict: node_name → {(year, month) → SST (°C)}.
    """
    data = defaultdict(dict)
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            dt = datetime.strptime(row['time'][:10], '%Y-%m-%d')
            yr, mo = dt.year, dt.month
            for node in NODES:
                if node in row:
                    data[node][(yr, mo)] = float(row[node])
    return dict(data)


def discover_cmip6_files(cmip6_dir: str) -> dict:
    """Discover CMIP6 CSV files organized by (model, scenario).

    File naming convention: {Model}_{scenario}_tos_monthly.csv

    Returns dict: scenario → [list of (model_name, file_path)].
    """
    files_by_scenario = defaultdict(list)
    for fname in sorted(os.listdir(cmip6_dir)):
        if not fname.endswith('_tos_monthly.csv'):
            continue
        parts = fname.replace('_tos_monthly.csv', '').rsplit('_', 1)
        if len(parts) != 2:
            print(f"  Skipping unrecognized file: {fname}", file=sys.stderr)
            continue
        model_name, scenario = parts
        files_by_scenario[scenario].append(
            (model_name, os.path.join(cmip6_dir, fname))
        )
    return dict(files_by_scenario)


def compute_monthly_bias(oisst: dict, cmip6: dict,
                         start_year: int = OVERLAP_START,
                         end_year: int = OVERLAP_END) -> np.ndarray:
    """Compute per-month bias: mean(CMIP6 - OISST) over overlap period.

    Returns array of shape (12,) with bias for months 1-12.
    """
    monthly_diffs = defaultdict(list)
    for yr in range(start_year, end_year + 1):
        for mo in range(1, 13):
            if (yr, mo) in oisst and (yr, mo) in cmip6:
                monthly_diffs[mo].append(cmip6[(yr, mo)] - oisst[(yr, mo)])

    bias = np.zeros(12)
    for mo in range(1, 13):
        if monthly_diffs[mo]:
            bias[mo - 1] = np.mean(monthly_diffs[mo])
        else:
            print(f"  Warning: no overlap data for month {mo}", file=sys.stderr)
    return bias


def process_scenario(scenario: str, model_files: list, oisst_dir: str,
                     output_dir: str, projection_start: int = 2026,
                     projection_end: int = 2100) -> None:
    """Process all models for one scenario, producing bias-corrected ensemble CSVs.

    For each node:
      1. Load each model's data
      2. Compute per-model, per-month bias correction from overlap period
      3. Apply bias correction to projection years
      4. Average across models (ensemble mean)
      5. Write output CSV
    """
    print(f"\nProcessing scenario: {scenario}")
    print(f"  Models: {[m for m, _ in model_files]}")

    os.makedirs(output_dir, exist_ok=True)

    for node in NODES:
        # Load OISST for this node
        oisst = load_oisst_monthly(node, oisst_dir)

        # Collect bias-corrected projections from each model
        corrected_by_model = []

        for model_name, fpath in model_files:
            cmip6_data = load_cmip6_csv(fpath)
            if node not in cmip6_data:
                print(f"  Warning: {node} not found in {model_name}, skipping",
                      file=sys.stderr)
                continue

            node_cmip6 = cmip6_data[node]

            # Compute per-month bias for this model at this node
            bias = compute_monthly_bias(oisst, node_cmip6)

            # Apply bias correction to projection years
            corrected = {}
            for yr in range(projection_start, projection_end + 1):
                for mo in range(1, 13):
                    if (yr, mo) in node_cmip6:
                        corrected[(yr, mo)] = node_cmip6[(yr, mo)] - bias[mo - 1]

            corrected_by_model.append(corrected)

        if not corrected_by_model:
            print(f"  Warning: no models available for {node}/{scenario}, "
                  f"skipping", file=sys.stderr)
            continue

        # Ensemble mean across models
        ensemble = {}
        all_keys = set()
        for c in corrected_by_model:
            all_keys.update(c.keys())

        for key in sorted(all_keys):
            vals = [c[key] for c in corrected_by_model if key in c]
            ensemble[key] = np.mean(vals)

        # Write output CSV
        out_path = os.path.join(output_dir, f"{node}_{scenario}_monthly.csv")
        with open(out_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['year', 'month', 'sst'])
            for (yr, mo) in sorted(ensemble.keys()):
                writer.writerow([yr, mo, f"{ensemble[(yr, mo)]:.4f}"])

        n_years = len(set(yr for yr, _ in ensemble.keys()))
        print(f"  {node}: {len(ensemble)} months ({n_years} years), "
              f"{len(corrected_by_model)} models → {out_path}")


def generate_trend_extrapolation(oisst_dir: str, output_dir: str,
                                 scenario_name: str,
                                 trend_multiplier: float = 1.0,
                                 projection_start: int = 2026,
                                 projection_end: int = 2100) -> None:
    """Generate placeholder projections using linear trend extrapolation from OISST.

    CLEARLY LABELED AS PLACEHOLDER — not real CMIP6 projections.

    For each node:
      1. Fit linear trend to OISST monthly data (per month-of-year)
      2. Extrapolate forward with optional trend multiplier
      3. Write output CSV with same schema

    The trend_multiplier scales the observed trend:
      1.0 = linear extrapolation of observed warming
      ~1.5 = approximate SSP2-4.5 acceleration
      ~2.5 = approximate SSP5-8.5 acceleration

    Args:
        oisst_dir: Directory with OISST *_monthly.csv files.
        output_dir: Output directory.
        scenario_name: Scenario label (e.g., 'ssp585_placeholder').
        trend_multiplier: Scaling factor for observed trend.
        projection_start: First projection year.
        projection_end: Last projection year.
    """
    print(f"\nGenerating trend extrapolation: {scenario_name} "
          f"(trend_multiplier={trend_multiplier})")

    os.makedirs(output_dir, exist_ok=True)

    for node in NODES:
        oisst = load_oisst_monthly(node, oisst_dir)
        if not oisst:
            print(f"  Warning: no OISST data for {node}, skipping",
                  file=sys.stderr)
            continue

        # Fit per-month linear trends
        projections = {}
        for mo in range(1, 13):
            years = []
            ssts = []
            for (yr, m), sst in sorted(oisst.items()):
                if m == mo:
                    years.append(yr)
                    ssts.append(sst)

            if len(years) < 5:
                continue

            years_arr = np.array(years, dtype=np.float64)
            ssts_arr = np.array(ssts, dtype=np.float64)

            # Linear regression
            coeffs = np.polyfit(years_arr, ssts_arr, 1)
            slope, intercept = coeffs

            # Extrapolate with scaled trend
            for yr in range(projection_start, projection_end + 1):
                # Base = climatological mean + observed trend to 2025
                base_2025 = slope * 2025 + intercept
                # Additional warming beyond 2025 scaled by multiplier
                delta = slope * trend_multiplier * (yr - 2025)
                projections[(yr, mo)] = base_2025 + delta

        # Write output
        out_path = os.path.join(output_dir, f"{node}_{scenario_name}_monthly.csv")
        with open(out_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['year', 'month', 'sst'])
            for (yr, mo) in sorted(projections.keys()):
                writer.writerow([yr, mo, f"{projections[(yr, mo)]:.4f}"])

        n_years = len(set(yr for yr, _ in projections.keys()))
        print(f"  {node}: {len(projections)} months ({n_years} years) → {out_path}")


def write_bias_report(cmip6_dir: str, oisst_dir: str, output_dir: str) -> None:
    """Write a bias correction report showing per-node, per-model offsets."""
    report_path = os.path.join(output_dir, "bias_correction_report.txt")
    files_by_scenario = discover_cmip6_files(cmip6_dir)

    with open(report_path, 'w') as f:
        f.write("CMIP6 Bias Correction Report\n")
        f.write(f"Overlap period: {OVERLAP_START}-{OVERLAP_END}\n")
        f.write("=" * 70 + "\n\n")

        for scenario, model_files in sorted(files_by_scenario.items()):
            f.write(f"Scenario: {scenario}\n")
            f.write("-" * 40 + "\n")

            for model_name, fpath in model_files:
                f.write(f"\n  Model: {model_name}\n")
                cmip6_data = load_cmip6_csv(fpath)

                for node in NODES:
                    oisst = load_oisst_monthly(node, oisst_dir)
                    if node not in cmip6_data:
                        f.write(f"    {node}: NO DATA\n")
                        continue

                    bias = compute_monthly_bias(oisst, cmip6_data[node])
                    mean_bias = np.mean(bias)
                    f.write(f"    {node}: mean bias = {mean_bias:+.2f}°C "
                            f"(range {np.min(bias):+.2f} to {np.max(bias):+.2f})\n")

            f.write("\n")

    print(f"Bias correction report: {report_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Process CMIP6 SST data into model-ready format"
    )
    parser.add_argument('--cmip6-dir', default='data/sst/cmip6',
                        help='Directory with raw CMIP6 CSV files')
    parser.add_argument('--oisst-dir', default='data/sst',
                        help='Directory with OISST *_monthly.csv files')
    parser.add_argument('--output-dir', default='data/sst/projections',
                        help='Output directory for processed projections')
    parser.add_argument('--scenarios', default='ssp245,ssp585',
                        help='Comma-separated list of scenarios to process')
    parser.add_argument('--projection-start', type=int, default=2026,
                        help='First projection year')
    parser.add_argument('--projection-end', type=int, default=2100,
                        help='Last projection year')
    args = parser.parse_args()

    scenarios = [s.strip() for s in args.scenarios.split(',')]

    # Discover available CMIP6 files
    print("Discovering CMIP6 files...")
    files_by_scenario = discover_cmip6_files(args.cmip6_dir)
    print(f"Found scenarios: {list(files_by_scenario.keys())}")
    for sc, files in files_by_scenario.items():
        print(f"  {sc}: {[m for m, _ in files]}")

    # Write bias report
    write_bias_report(args.cmip6_dir, args.oisst_dir, args.output_dir)

    # Process each requested scenario
    for scenario in scenarios:
        if scenario in files_by_scenario:
            # Real CMIP6 data available
            process_scenario(
                scenario=scenario,
                model_files=files_by_scenario[scenario],
                oisst_dir=args.oisst_dir,
                output_dir=args.output_dir,
                projection_start=args.projection_start,
                projection_end=args.projection_end,
            )
        else:
            # No real data — generate trend extrapolation placeholder
            print(f"\n  No CMIP6 data for {scenario}, generating placeholder "
                  f"(trend extrapolation)")

            # SSP2-4.5 ≈ 1.5× observed trend, SSP5-8.5 ≈ 2.5× observed
            multipliers = {
                'ssp126': 0.7,
                'ssp245': 1.5,
                'ssp370': 2.0,
                'ssp585': 2.5,
            }
            mult = multipliers.get(scenario, 1.5)

            generate_trend_extrapolation(
                oisst_dir=args.oisst_dir,
                output_dir=args.output_dir,
                scenario_name=f"{scenario}_placeholder",
                trend_multiplier=mult,
                projection_start=args.projection_start,
                projection_end=args.projection_end,
            )

    print("\nDone.")


if __name__ == '__main__':
    main()
