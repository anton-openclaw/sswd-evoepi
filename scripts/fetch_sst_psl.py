#!/usr/bin/env python3
"""
Fetch satellite SST data from NOAA PSL OPeNDAP (OISST v2.1 monthly means).

Coastwatch ERDDAP is currently unresponsive, so we use PSL's THREDDS OPeNDAP
for the monthly mean SST dataset and interpolate to daily climatology.

Source: NOAA OI SST V2 High Resolution Dataset
  - 0.25° daily, monthly means
  - https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.highres.html
  - OPeNDAP: https://psl.noaa.gov/thredds/dodsC/Datasets/noaa.oisst.v2.highres/sst.mon.mean.nc

Period: Sep 1981 – present (we use 2002-2025 for climatology)
"""

import csv
import io
import json
import os
import sys
import time
import urllib.request
import urllib.error
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import numpy as np

DATA_DIR = Path(__file__).resolve().parent.parent / "data" / "sst"

# Node coordinates
NODES = [
    ("Sitka",         57.06, -135.34),
    ("Ketchikan",     55.34, -131.64),
    ("Haida_Gwaii",   53.25, -132.07),
    ("Bella_Bella",   52.16, -128.15),
    ("Howe_Sound",    49.52, -123.25),
    ("SJI",           48.53, -123.02),
    ("Westport",      46.89, -124.10),
    ("Newport",       44.63, -124.05),
    ("Crescent_City", 41.76, -124.20),
    ("Fort_Bragg",    39.45, -123.80),
    ("Monterey",      36.62, -121.90),
]

# OPeNDAP base
OPENDAP_BASE = "https://psl.noaa.gov/thredds/dodsC/Datasets/noaa.oisst.v2.highres/sst.mon.mean.nc"

# Grid resolution = 0.25°
# lat: -89.875 to 89.875, step 0.25 → index = (lat + 89.875) / 0.25
# lon: 0.125 to 359.875, step 0.25 → index = (lon % 360) / 0.25
# time: monthly since 1981-09-01. Index 0 = Sep 1981.
# We want 2002-01 to 2025-12 → that's month indices for those dates.

# Sep 1981 = index 0.  Jan 1982 = index 4.
# Jan 2002 = (2002-1981)*12 + (1-9) = 21*12 - 8 = 252 - 8 = 244
# Actually: Sep 1981 = 0, Oct 1981 = 1, ..., Dec 1981 = 3, Jan 1982 = 4
# Jan Y = (Y - 1981)*12 + (1-9) = (Y-1981)*12 - 8
# Jan 2002 = 21*12 - 8 = 252 - 8 = 244
# Dec 2025 = (2025-1981)*12 + (12-9) = 44*12 + 3 = 528 + 3 = 531
# But dataset may not extend to Dec 2025 yet; we'll try and handle truncation.

START_YEAR = 2002
END_YEAR = 2025

def time_index(year, month):
    """Convert year/month to OPeNDAP time index (0-based, Sep 1981 = 0)."""
    return (year - 1981) * 12 + (month - 9)

def lat_index(lat):
    """Convert latitude to grid index."""
    return round((lat + 89.875) / 0.25)

def lon_index(lon):
    """Convert longitude to grid index (handles negative lon)."""
    lon360 = lon % 360
    return round((lon360 - 0.125) / 0.25)


def fetch_opendap_ascii(url, timeout=60, retries=3):
    """Fetch OPeNDAP ASCII data with retries."""
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "sswd-evoepi/1.0"})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read().decode("utf-8")
        except Exception as e:
            print(f"    Attempt {attempt+1}/{retries}: {e}", flush=True)
            if attempt < retries - 1:
                time.sleep(5 * (attempt + 1))
    raise ConnectionError(f"Failed after {retries} attempts")


def parse_opendap_sst(text):
    """Parse OPeNDAP ASCII response for SST array.
    
    Format is:
    Dataset { ... } name;
    ---
    sst.sst[T][1][1]
    [t], value
    ...
    sst.time[T]
    t0, t1, ...
    sst.lat[1]
    latval
    sst.lon[1]
    lonval
    """
    lines = text.strip().split("\n")
    
    # Find the data section after "-----"
    data_start = None
    for i, line in enumerate(lines):
        if line.startswith("---"):
            data_start = i + 1
            break
    
    if data_start is None:
        return [], []
    
    # Parse SST values  
    sst_values = []
    time_values = []
    
    parsing_sst = False
    parsing_time = False
    
    for i in range(data_start, len(lines)):
        line = lines[i].strip()
        if not line:
            continue
            
        if line.startswith("sst.sst["):
            parsing_sst = True
            parsing_time = False
            continue
        elif line.startswith("sst.time["):
            parsing_sst = False
            parsing_time = True
            continue
        elif line.startswith("sst.lat[") or line.startswith("sst.lon["):
            parsing_sst = False
            parsing_time = False
            continue
        
        if parsing_sst:
            # Format: [idx][0][0], value  or  [idx], value
            parts = line.split(",")
            if len(parts) >= 2:
                val_str = parts[-1].strip()
                try:
                    val = float(val_str)
                    # Fill value is -9.96921E36
                    if val < -1000:
                        val = float("nan")
                    sst_values.append(val)
                except ValueError:
                    sst_values.append(float("nan"))
        elif parsing_time:
            # Comma-separated time values (days since 1800-01-01)
            for part in line.split(","):
                part = part.strip()
                if part:
                    try:
                        time_values.append(float(part))
                    except ValueError:
                        pass
    
    return time_values, sst_values


def days_since_1800_to_date(days):
    """Convert 'days since 1800-01-01' to (year, month)."""
    from datetime import datetime, timedelta
    base = datetime(1800, 1, 1)
    dt = base + timedelta(days=days)
    return dt.year, dt.month


def fetch_node_monthly(name, lat, lon):
    """Fetch monthly SST for one node from PSL OPeNDAP."""
    t0 = time.time()
    
    li = lat_index(lat)
    loi = lon_index(lon)
    t_start = time_index(START_YEAR, 1)
    t_end = time_index(END_YEAR, 12)
    
    print(f"  [{name}] lat_idx={li}, lon_idx={loi}, t=[{t_start}:{t_end}]", flush=True)
    
    # Try original coordinates first
    url = f"{OPENDAP_BASE}.ascii?sst[{t_start}:{t_end}][{li}:{li}][{loi}:{loi}]"
    
    try:
        text = fetch_opendap_ascii(url, timeout=120)
        times, values = parse_opendap_sst(text)
    except Exception as e:
        print(f"  [{name}] Failed: {e}", flush=True)
        return None
    
    # Check for fill values (coastal point hitting land)
    nan_count = sum(1 for v in values if np.isnan(v))
    nan_rate = nan_count / max(len(values), 1)
    
    if nan_rate > 0.5:
        print(f"  [{name}] {nan_rate:.0%} NaN — trying offshore offsets...", flush=True)
        
        # Try offsets (move seaward — generally westward for west coast)
        offsets = [
            (0, -1),   # 0.25° west
            (0, -2),   # 0.50° west
            (-1, 0),   # 0.25° south
            (1, 0),    # 0.25° north
            (-1, -1),  # SW
            (0, 1),    # 0.25° east (for fjords that may be inside)
        ]
        
        best_values = values
        best_times = times
        best_nan = nan_rate
        best_offset = (0, 0)
        
        for dlat, dlon in offsets:
            try:
                url2 = f"{OPENDAP_BASE}.ascii?sst[{t_start}:{t_end}][{li+dlat}:{li+dlat}][{loi+dlon}:{loi+dlon}]"
                text2 = fetch_opendap_ascii(url2, timeout=120)
                t2, v2 = parse_opendap_sst(text2)
                nr = sum(1 for v in v2 if np.isnan(v)) / max(len(v2), 1)
                print(f"  [{name}]   offset ({dlat},{dlon}): {nr:.0%} NaN ({len(v2)} months)", flush=True)
                if nr < best_nan:
                    best_values = v2
                    best_times = t2
                    best_nan = nr
                    best_offset = (dlat, dlon)
                    if nr < 0.05:
                        break
            except Exception as e:
                print(f"  [{name}]   offset ({dlat},{dlon}): error {e}", flush=True)
        
        values = best_values
        times = best_times
        nan_rate = best_nan
        if best_offset != (0, 0):
            print(f"  [{name}] Using offset {best_offset} ({best_nan:.0%} NaN)", flush=True)
    
    # Convert times to year-month
    months = []
    for t in times:
        y, m = days_since_1800_to_date(t)
        months.append((y, m))
    
    elapsed = time.time() - t0
    valid = [v for v in values if not np.isnan(v)]
    
    print(f"  [{name}] {len(values)} months, {nan_rate:.1%} NaN, "
          f"{min(valid):.1f}–{max(valid):.1f}°C (mean {np.mean(valid):.1f}), "
          f"{elapsed:.1f}s", flush=True)
    
    return {
        "name": name,
        "lat": lat, "lon": lon,
        "months": months,
        "values": values,
        "nan_rate": nan_rate,
        "elapsed": elapsed,
    }


def monthly_to_daily_climatology(months, values):
    """Convert monthly time series to 365-day climatology.
    
    1. Compute mean for each calendar month across years
    2. Interpolate monthly means (assigned to mid-month) to daily via cubic-like interp
    """
    # Group by calendar month
    month_means = defaultdict(list)
    for (y, m), v in zip(months, values):
        if not np.isnan(v):
            month_means[m].append(v)
    
    # Monthly climatology (12 values)
    monthly_clim = {}
    for m in range(1, 13):
        vals = month_means.get(m, [])
        if vals:
            monthly_clim[m] = {
                "mean": float(np.mean(vals)),
                "std": float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0,
                "n": len(vals),
            }
        else:
            monthly_clim[m] = {"mean": float("nan"), "std": 0.0, "n": 0}
    
    # Mid-month day-of-year for interpolation anchor points
    # Jan=16, Feb=45, ..., Dec=350
    mid_month_doy = [16, 45, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349]
    
    # Get monthly means in order, handling any missing with linear interpolation
    monthly_means = [monthly_clim[m]["mean"] for m in range(1, 13)]
    monthly_stds = [monthly_clim[m]["std"] for m in range(1, 13)]
    monthly_ns = [monthly_clim[m]["n"] for m in range(1, 13)]
    
    # Interpolate to daily using periodic (wrap-around) linear interpolation
    daily_clim = {}
    
    # Extended arrays for wrap-around
    ext_doys = [mid_month_doy[-1] - 365] + mid_month_doy + [mid_month_doy[0] + 365]
    ext_means = [monthly_means[-1]] + monthly_means + [monthly_means[0]]
    ext_stds = [monthly_stds[-1]] + monthly_stds + [monthly_stds[0]]
    
    for doy in range(1, 366):
        # Find bracketing months
        for i in range(len(ext_doys) - 1):
            if ext_doys[i] <= doy <= ext_doys[i + 1]:
                frac = (doy - ext_doys[i]) / (ext_doys[i + 1] - ext_doys[i])
                mean_val = ext_means[i] + frac * (ext_means[i + 1] - ext_means[i])
                std_val = ext_stds[i] + frac * (ext_stds[i + 1] - ext_stds[i])
                # n = average across years for this month
                break
        else:
            mean_val = monthly_means[0]
            std_val = monthly_stds[0]
        
        daily_clim[doy] = {
            "mean": float(mean_val),
            "std": float(std_val),
            "n": int(np.mean(monthly_ns)),  # approximate
        }
    
    return daily_clim, monthly_clim


def save_monthly_csv(name, months, values):
    """Save raw monthly data."""
    path = DATA_DIR / f"{name}_monthly.csv"
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["year", "month", "sst"])
        for (y, m), v in zip(months, values):
            w.writerow([y, m, f"{v:.4f}" if not np.isnan(v) else "NaN"])
    return path


def save_climatology_csv(name, clim):
    """Save 365-day climatology."""
    path = DATA_DIR / f"{name}_climatology.csv"
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["day_of_year", "sst_mean", "sst_std", "n_years"])
        for doy in range(1, 366):
            c = clim[doy]
            w.writerow([
                doy,
                f"{c['mean']:.4f}" if not np.isnan(c['mean']) else "NaN",
                f"{c['std']:.4f}" if not np.isnan(c['std']) else "NaN",
                c['n'],
            ])
    return path


def save_combined_climatology(all_clim):
    """Save combined climatology for all nodes."""
    path = DATA_DIR / "all_nodes_climatology.csv"
    names = list(all_clim.keys())
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        header = ["day_of_year"]
        for n in names:
            header.extend([f"{n}_mean", f"{n}_std"])
        w.writerow(header)
        for doy in range(1, 366):
            row = [doy]
            for n in names:
                c = all_clim[n].get(doy, {"mean": float("nan"), "std": float("nan")})
                row.append(f"{c['mean']:.4f}" if not np.isnan(c['mean']) else "NaN")
                row.append(f"{c['std']:.4f}" if not np.isnan(c['std']) else "NaN")
            w.writerow(row)
    return path


def plot_climatologies(all_clim, all_monthly, model_params):
    """Plot all climatologies + sinusoidal comparison."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  matplotlib not available, skipping plot", flush=True)
        return None
    
    fig, axes = plt.subplots(3, 4, figsize=(20, 14))
    axes = axes.flatten()
    
    cmap = plt.cm.tab10
    doys = list(range(1, 366))
    
    for i, (name, clim) in enumerate(all_clim.items()):
        ax = axes[i]
        
        # Satellite climatology
        means = [clim[d]["mean"] for d in doys]
        ax.plot(doys, means, color="steelblue", linewidth=2, label="OISST satellite")
        
        # Monthly points
        if name in all_monthly:
            mc = all_monthly[name]
            mid_doys = [16, 45, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349]
            m_means = [mc[m]["mean"] for m in range(1, 13)]
            m_stds = [mc[m]["std"] for m in range(1, 13)]
            ax.errorbar(mid_doys, m_means, yerr=m_stds, fmt='o', color="navy",
                       markersize=4, capsize=3, label="Monthly ±1σ")
        
        # Sinusoidal model
        if name in model_params:
            mp = model_params[name]
            sin_sst = [mp["mean"] + mp["amp"] * np.cos(2 * np.pi * (d - 1 - 45) / 365)
                       for d in doys]
            ax.plot(doys, sin_sst, color="orangered", linewidth=1.5, linestyle="--",
                   label=f"Sinusoidal (μ={mp['mean']}, A={mp['amp']})")
        
        ax.set_title(name.replace("_", " "), fontsize=11, fontweight="bold")
        ax.set_xlim(1, 365)
        ax.grid(True, alpha=0.3)
        
        month_starts = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]
        month_labels = ["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"]
        ax.set_xticks(month_starts)
        ax.set_xticklabels(month_labels, fontsize=8)
        ax.set_ylabel("SST (°C)", fontsize=9)
        
        if i == 0:
            ax.legend(fontsize=7, loc="upper left")
    
    # Hide last subplot (12th, we only have 11)
    axes[11].set_visible(False)
    
    fig.suptitle("NOAA OISST v2.1 Climatology (2002–2025) vs Model Sinusoidal\n"
                 "11 Stepping-Stone Network Nodes", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    
    path = DATA_DIR / "validation_plots.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  Plot saved: {path}", flush=True)
    
    # Also make overlay plot
    fig2, ax2 = plt.subplots(figsize=(14, 8))
    colors = plt.cm.coolwarm(np.linspace(0, 1, len(all_clim)))
    
    for idx, (name, clim) in enumerate(all_clim.items()):
        means = [clim[d]["mean"] for d in doys]
        ax2.plot(doys, means, color=colors[idx], linewidth=1.5,
                label=name.replace("_", " "))
    
    ax2.set_xlabel("Day of Year", fontsize=12)
    ax2.set_ylabel("SST (°C)", fontsize=12)
    ax2.set_title("Satellite SST Climatology (2002–2025) — All 11 Nodes", fontsize=14)
    ax2.legend(loc="upper left", fontsize=9, ncol=2)
    ax2.set_xlim(1, 365)
    ax2.grid(True, alpha=0.3)
    ax2.set_xticks(month_starts)
    ax2.set_xticklabels(["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                         "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
    
    path2 = DATA_DIR / "climatology_curves.png"
    fig2.tight_layout()
    fig2.savefig(path2, dpi=150)
    plt.close(fig2)
    print(f"  Overlay plot saved: {path2}", flush=True)
    
    return path


def main():
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    t_start = time.time()
    
    # Model sinusoidal parameters for comparison
    model_params = {
        "Sitka":         {"mean": 8.0,  "amp": 3.5},
        "Ketchikan":     {"mean": 8.5,  "amp": 3.5},
        "Haida_Gwaii":   {"mean": 9.0,  "amp": 3.0},
        "Bella_Bella":   {"mean": 9.5,  "amp": 3.5},
        "Howe_Sound":    {"mean": 10.0, "amp": 4.0},
        "SJI":           {"mean": 10.5, "amp": 4.0},
        "Westport":      {"mean": 11.0, "amp": 3.5},
        "Newport":       {"mean": 11.5, "amp": 3.0},
        "Crescent_City": {"mean": 12.0, "amp": 2.5},
        "Fort_Bragg":    {"mean": 12.5, "amp": 2.5},
        "Monterey":      {"mean": 13.0, "amp": 2.5},
    }
    
    print("=" * 60, flush=True)
    print("  SST Data Extraction — NOAA PSL OPeNDAP", flush=True)
    print(f"  Period: {START_YEAR}–{END_YEAR}", flush=True)
    print(f"  Source: OISST v2.1 monthly means", flush=True)
    print(f"  Output: {DATA_DIR}", flush=True)
    print("=" * 60, flush=True)
    
    all_clim = {}
    all_monthly = {}
    summaries = []
    
    for name, lat, lon in NODES:
        print(f"\n--- {name} ({lat}, {lon}) ---", flush=True)
        result = fetch_node_monthly(name, lat, lon)
        
        if result is None:
            summaries.append({
                "name": name, "source": "FAILED",
                "total_months": 0, "nan_rate": 1.0,
                "sst_min": None, "sst_max": None, "sst_mean": None,
            })
            continue
        
        # Save monthly
        save_monthly_csv(name, result["months"], result["values"])
        
        # Compute climatology
        daily_clim, monthly_clim = monthly_to_daily_climatology(
            result["months"], result["values"]
        )
        save_climatology_csv(name, daily_clim)
        all_clim[name] = daily_clim
        all_monthly[name] = monthly_clim
        
        valid = [v for v in result["values"] if not np.isnan(v)]
        years_covered = len(set(y for (y, m) in result["months"] 
                               if not np.isnan(result["values"][result["months"].index((y,m))])))
        
        summaries.append({
            "name": name,
            "source": "OISST_v2.1_PSL",
            "total_months": len(result["values"]),
            "nan_rate": result["nan_rate"],
            "sst_min": float(min(valid)) if valid else None,
            "sst_max": float(max(valid)) if valid else None,
            "sst_mean": float(np.mean(valid)) if valid else None,
            "years_covered": years_covered,
            "elapsed_s": result["elapsed"],
        })
        
        # Brief pause between nodes
        time.sleep(1)
    
    # Combined climatology
    if all_clim:
        save_combined_climatology(all_clim)
        print(f"\n  Combined climatology saved", flush=True)
    
    # Plots
    if all_clim:
        plot_climatologies(all_clim, all_monthly, model_params)
    
    # Summary JSON
    summary_path = DATA_DIR / "fetch_summary.json"
    with open(summary_path, "w") as f:
        json.dump({
            "fetch_time": datetime.now().isoformat(),
            "period": f"{START_YEAR}-{END_YEAR}",
            "source": "NOAA PSL OPeNDAP (OISST v2.1 monthly means)",
            "note": "Coastwatch ERDDAP was unresponsive; used PSL THREDDS instead. "
                    "Monthly means interpolated to daily climatology.",
            "total_elapsed_s": time.time() - t_start,
            "nodes": summaries,
        }, f, indent=2, default=str)
    
    # Print summary
    elapsed = time.time() - t_start
    print(f"\n{'='*80}", flush=True)
    print(f"  TOTAL TIME: {elapsed:.0f}s ({elapsed/60:.1f} min)", flush=True)
    print(f"  {'Node':<15} {'Months':>7} {'NaN%':>6} {'Min':>6} {'Max':>6} {'Mean':>6} {'Years':>6}", flush=True)
    print("-" * 80, flush=True)
    for s in summaries:
        if s.get("sst_min") is not None:
            print(f"  {s['name']:<15} {s['total_months']:>7} {s['nan_rate']:>5.1%} "
                  f"{s['sst_min']:>6.1f} {s['sst_max']:>6.1f} {s['sst_mean']:>6.1f} "
                  f"{s.get('years_covered', '?'):>6}", flush=True)
        else:
            print(f"  {s['name']:<15} {'FAILED':>7}", flush=True)
    print("=" * 80, flush=True)
    
    failures = [s for s in summaries if s.get("total_months", 0) == 0]
    if failures:
        print(f"\n  ⚠️  {len(failures)} node(s) failed", flush=True)
    else:
        print(f"\n  ✅ All {len(summaries)} nodes fetched successfully!", flush=True)
    
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())
