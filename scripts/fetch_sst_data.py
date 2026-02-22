#!/usr/bin/env python3
"""
Fetch satellite SST data for all 11 stepping-stone network nodes.

Primary: NOAA OISST v2.1 (0.25°, daily, 1981-present) via ERDDAP
Fallback: MUR SST (1km, daily, 2002-present) for coastal/fjord nodes with NaN issues

Outputs:
  data/sst/{node}_daily.csv         — raw daily SST time series
  data/sst/{node}_climatology.csv   — 365-day climatology (mean, std, n)
  data/sst/all_nodes_climatology.csv — combined climatology
  data/sst/climatology_curves.png   — comparison plot (if matplotlib available)
  data/sst/fetch_summary.json       — metadata about the fetch
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
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Node definitions
# ---------------------------------------------------------------------------
NODES = [
    ("Sitka",        57.06, -135.34),
    ("Ketchikan",    55.34, -131.64),
    ("Haida_Gwaii",  53.25, -132.07),
    ("Bella_Bella",  52.16, -128.15),
    ("Howe_Sound",   49.52, -123.25),
    ("SJI",          48.53, -123.02),
    ("Westport",     46.89, -124.10),
    ("Newport",      44.63, -124.05),
    ("Crescent_City",41.76, -124.20),
    ("Fort_Bragg",   39.45, -123.80),
    ("Monterey",     36.62, -121.90),
]

# ERDDAP base URLs
OISST_BASE = "https://coastwatch.pfeg.noaa.gov/erddap/griddap/ncdcOisst21Agg_LonPM180.csv"
MUR_BASE   = "https://coastwatch.pfeg.noaa.gov/erddap/griddap/jplMURSST41.csv"

DATA_DIR = Path(__file__).resolve().parent.parent / "data" / "sst"

START_YEAR = 2002
END_YEAR   = 2025
CHUNK_YEARS = 8  # years per request — large chunks are efficient

MAX_RETRIES = 3
RETRY_DELAY = 10
CONCURRENT_WORKERS = 2  # parallel node fetches (ERDDAP 408s with >2)


def fetch_url(url: str, retries: int = MAX_RETRIES) -> str:
    """Fetch URL with retry logic."""
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "sswd-evoepi/1.0"})
            with urllib.request.urlopen(req, timeout=300) as resp:
                return resp.read().decode("utf-8")
        except urllib.error.HTTPError as e:
            if e.code in (413, 404):
                raise
            print(f"    HTTP {e.code} attempt {attempt+1}/{retries}")
        except (urllib.error.URLError, TimeoutError, OSError) as e:
            print(f"    Connection error attempt {attempt+1}/{retries}: {e}")
        if attempt < retries - 1:
            time.sleep(RETRY_DELAY * (attempt + 1))
    raise ConnectionError(f"Failed after {retries} attempts: {url[:120]}...")


def parse_erddap_csv(text: str, sst_var: str = "sst") -> list[tuple[str, float]]:
    """Parse ERDDAP CSV → list of (date_str, sst_value)."""
    lines = text.strip().split("\n")
    if len(lines) < 3:
        return []
    reader = csv.reader(io.StringIO("\n".join(lines)))
    header = next(reader)
    _units = next(reader)

    time_idx = sst_idx = None
    for i, col in enumerate(header):
        c = col.strip().lower()
        if c == "time":
            time_idx = i
        if c == sst_var.lower() or c == "analysed_sst":
            sst_idx = i

    if time_idx is None or sst_idx is None:
        return []

    records = []
    for row in reader:
        if len(row) <= max(time_idx, sst_idx):
            continue
        date_str = row[time_idx].strip()[:10]
        sst_raw = row[sst_idx].strip()
        if sst_raw == "" or sst_raw.lower() == "nan":
            sst_val = float("nan")
        else:
            try:
                sst_val = float(sst_raw)
            except ValueError:
                sst_val = float("nan")
        records.append((date_str, sst_val))
    return records


def make_year_chunks(start: int, end: int, chunk_size: int) -> list[tuple[int, int]]:
    """Split year range into chunks."""
    chunks = []
    y = start
    while y <= end:
        y_end = min(y + chunk_size - 1, end)
        chunks.append((y, y_end))
        y = y_end + 1
    return chunks


def fetch_oisst_chunk(lat: float, lon: float, y_start: int, y_end: int) -> list[tuple[str, float]]:
    """Fetch a multi-year chunk of OISST data."""
    start_d = f"{y_start}-01-01"
    end_d = f"{y_end}-12-31"
    url = (
        f"{OISST_BASE}?sst[({start_d}):1:({end_d})]"
        f"[(0.0):1:(0.0)]"
        f"[({lat}):1:({lat})]"
        f"[({lon}):1:({lon})]"
    )
    try:
        text = fetch_url(url)
        return parse_erddap_csv(text, "sst")
    except urllib.error.HTTPError as e:
        if e.code == 404 and y_end == END_YEAR:
            # Try up to today
            today = datetime.now().strftime("%Y-%m-%d")
            url2 = (
                f"{OISST_BASE}?sst[({start_d}):1:({today})]"
                f"[(0.0):1:(0.0)]"
                f"[({lat}):1:({lat})]"
                f"[({lon}):1:({lon})]"
            )
            text = fetch_url(url2)
            return parse_erddap_csv(text, "sst")
        raise


def fetch_mur_chunk(lat: float, lon: float, y_start: int, y_end: int) -> list[tuple[str, float]]:
    """Fetch a multi-year chunk of MUR SST data."""
    start_d = f"{y_start}-01-01"
    end_d = f"{y_end}-12-31"
    url = (
        f"{MUR_BASE}?analysed_sst[({start_d}):1:({end_d})]"
        f"[({lat}):1:({lat})]"
        f"[({lon}):1:({lon})]"
    )
    try:
        text = fetch_url(url)
        return parse_erddap_csv(text, "analysed_sst")
    except urllib.error.HTTPError as e:
        if e.code == 404 and y_end == END_YEAR:
            today = datetime.now().strftime("%Y-%m-%d")
            url2 = (
                f"{MUR_BASE}?analysed_sst[({start_d}):1:({today})]"
                f"[({lat}):1:({lat})]"
                f"[({lon}):1:({lon})]"
            )
            text = fetch_url(url2)
            return parse_erddap_csv(text, "analysed_sst")
        raise


def interpolate_short_gaps(values: list[float], max_gap: int = 3) -> list[float]:
    """Linearly interpolate NaN gaps of max_gap days or fewer."""
    arr = np.array(values, dtype=float)
    n = len(arr)
    if n == 0:
        return values
    is_nan = np.isnan(arr)
    result = arr.copy()
    i = 0
    while i < n:
        if is_nan[i]:
            j = i
            while j < n and is_nan[j]:
                j += 1
            gap_len = j - i
            if gap_len <= max_gap:
                left = arr[i - 1] if i > 0 and not np.isnan(arr[i - 1]) else None
                right = arr[j] if j < n and not np.isnan(arr[j]) else None
                if left is not None and right is not None:
                    for k in range(i, j):
                        frac = (k - i + 1) / (gap_len + 1)
                        result[k] = left + frac * (right - left)
                elif left is not None:
                    result[i:j] = left
                elif right is not None:
                    result[i:j] = right
            i = j
        else:
            i += 1
    return result.tolist()


def compute_climatology(dates: list[str], values: list[float]) -> dict:
    """Compute 365-day climatology."""
    day_values = defaultdict(list)
    for date_str, val in zip(dates, values):
        if np.isnan(val):
            continue
        dt = datetime.strptime(date_str, "%Y-%m-%d")
        doy = min(dt.timetuple().tm_yday, 365)
        day_values[doy].append(val)

    clim = {}
    for doy in range(1, 366):
        vals = day_values.get(doy, [])
        if vals:
            clim[doy] = {
                "mean": float(np.mean(vals)),
                "std": float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0,
                "n": len(vals),
            }
        else:
            clim[doy] = {"mean": float("nan"), "std": float("nan"), "n": 0}
    return clim


# Coordinate offsets for coastal nodes hitting land in OISST
OFFSET_ATTEMPTS = [
    (0, 0),
    (0, -0.125),
    (0, -0.25),
    (-0.125, 0),
    (0.125, 0),
    (-0.125, -0.125),
    (0, 0.125),
]


def fetch_node_all(name: str, lat: float, lon: float) -> dict:
    """Fetch complete SST time series for one node."""
    chunks = make_year_chunks(START_YEAR, END_YEAR, CHUNK_YEARS)
    all_records = []
    source = "OISST"
    used_lat, used_lon = lat, lon
    t0 = time.time()

    print(f"  [{name}] Fetching OISST ({lat}, {lon}) — {len(chunks)} chunks...")

    for ci, (y_start, y_end) in enumerate(chunks):
        if ci > 0:
            time.sleep(2)  # breathing room between chunks
        label = f"{y_start}-{y_end}"
        try:
            recs = fetch_oisst_chunk(lat, lon, y_start, y_end)
            nan_c = sum(1 for _, v in recs if np.isnan(v))
            print(f"  [{name}] OISST {label}: {len(recs)} days, {nan_c} NaN")
            all_records.extend(recs)
        except Exception as e:
            print(f"  [{name}] OISST {label}: ERROR {e}")

    # Evaluate NaN rate
    total = len(all_records)
    nan_total = sum(1 for _, v in all_records if np.isnan(v))
    nan_rate = nan_total / total if total > 0 else 1.0

    # If high NaN, try MUR with coordinate offsets
    if nan_rate > 0.5 or total == 0:
        print(f"  [{name}] OISST NaN rate {nan_rate:.1%} — trying MUR fallback...")
        for lat_off, lon_off in OFFSET_ATTEMPTS:
            mur_lat = lat + lat_off
            mur_lon = lon + lon_off
            off_label = f"({lat_off:+.3f},{lon_off:+.3f})" if (lat_off, lon_off) != (0, 0) else "(orig)"

            # Quick test with one chunk
            try:
                test = fetch_mur_chunk(mur_lat, mur_lon, 2010, 2010)
                test_nan = sum(1 for _, v in test if np.isnan(v))
                if not test or test_nan / max(len(test), 1) > 0.5:
                    print(f"  [{name}] MUR {off_label} test: {len(test)} days, {test_nan} NaN — skip")
                    continue
            except Exception as e:
                print(f"  [{name}] MUR {off_label} test error: {e}")
                continue

            # Fetch all chunks from MUR
            mur_records = []
            ok = True
            for y_start, y_end in chunks:
                try:
                    recs = fetch_mur_chunk(mur_lat, mur_lon, y_start, y_end)
                    mur_records.extend(recs)
                    print(f"  [{name}] MUR {off_label} {y_start}-{y_end}: {len(recs)} days")
                except Exception as e:
                    print(f"  [{name}] MUR {off_label} {y_start}-{y_end}: ERROR {e}")
                    ok = False
                    break

            if ok and mur_records:
                mur_nan = sum(1 for _, v in mur_records if np.isnan(v))
                mur_nan_rate = mur_nan / len(mur_records)
                if mur_nan_rate < nan_rate:
                    print(f"  [{name}] MUR better: {mur_nan_rate:.1%} vs OISST {nan_rate:.1%}")
                    all_records = mur_records
                    source = f"MUR{' '+off_label if (lat_off,lon_off)!=(0,0) else ''}"
                    used_lat, used_lon = mur_lat, mur_lon
                    nan_rate = mur_nan_rate
                    break

        # If MUR didn't help, try OISST with offsets
        if nan_rate > 0.5:
            print(f"  [{name}] Trying OISST with coordinate offsets...")
            for lat_off, lon_off in OFFSET_ATTEMPTS[1:]:
                off_lat, off_lon = lat + lat_off, lon + lon_off
                try:
                    test = fetch_oisst_chunk(off_lat, off_lon, 2010, 2010)
                    test_nan = sum(1 for _, v in test if np.isnan(v))
                    if test and test_nan / len(test) < 0.1:
                        # Fetch all
                        off_records = []
                        for y_start, y_end in chunks:
                            try:
                                recs = fetch_oisst_chunk(off_lat, off_lon, y_start, y_end)
                                off_records.extend(recs)
                            except Exception:
                                pass
                        off_nan = sum(1 for _, v in off_records if np.isnan(v))
                        off_nan_rate = off_nan / max(len(off_records), 1)
                        if off_nan_rate < nan_rate:
                            off_label = f"({lat_off:+.3f},{lon_off:+.3f})"
                            print(f"  [{name}] OISST offset {off_label}: {off_nan_rate:.1%} NaN — using it")
                            all_records = off_records
                            source = f"OISST {off_label}"
                            used_lat, used_lon = off_lat, off_lon
                            nan_rate = off_nan_rate
                            break
                except Exception:
                    pass

    # Separate dates and values
    dates = [r[0] for r in all_records]
    values = [r[1] for r in all_records]

    # Interpolate short gaps
    if values:
        values = interpolate_short_gaps(values, max_gap=3)

    # Stats
    final_nan = sum(1 for v in values if np.isnan(v))
    final_nan_rate = final_nan / max(len(values), 1)
    valid = [v for v in values if not np.isnan(v)]
    vmin = min(valid) if valid else float("nan")
    vmax = max(valid) if valid else float("nan")
    vmean = float(np.mean(valid)) if valid else float("nan")

    elapsed = time.time() - t0
    print(f"  [{name}] DONE in {elapsed:.0f}s — {source}, {len(values)} days, "
          f"{final_nan_rate:.1%} NaN, {vmin:.1f}–{vmax:.1f}°C (mean {vmean:.1f})")

    return {
        "name": name, "lat": lat, "lon": lon,
        "used_lat": used_lat, "used_lon": used_lon,
        "source": source, "dates": dates, "values": values,
        "total_days": len(values), "nan_count": final_nan,
        "nan_rate": final_nan_rate,
        "sst_min": vmin, "sst_max": vmax, "sst_mean": vmean,
        "suspicious": (vmin < -3 or vmax > 35) if valid else True,
        "elapsed_s": elapsed,
    }


def save_daily_csv(name: str, dates: list, values: list):
    path = DATA_DIR / f"{name}_daily.csv"
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["date", "sst"])
        for d, v in zip(dates, values):
            w.writerow([d, f"{v:.4f}" if not np.isnan(v) else "NaN"])
    return path


def save_climatology_csv(name: str, clim: dict):
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


def save_combined_climatology(all_clim: dict):
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


def plot_climatologies(all_clim: dict):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  matplotlib not available, skipping plot")
        return None

    fig, ax = plt.subplots(figsize=(14, 8))
    cmap = plt.cm.coolwarm
    n = len(all_clim)

    for i, (name, clim) in enumerate(all_clim.items()):
        doys = list(range(1, 366))
        means = [clim[d]["mean"] for d in doys]
        color = cmap(i / max(n - 1, 1))
        ax.plot(doys, means, label=name.replace("_", " "), color=color, linewidth=1.5)

    ax.set_xlabel("Day of Year", fontsize=12)
    ax.set_ylabel("SST (°C)", fontsize=12)
    ax.set_title("Satellite SST Climatology (2002–2025) — 11 Stepping-Stone Nodes", fontsize=14)
    ax.legend(loc="upper left", fontsize=9, ncol=2)
    ax.set_xlim(1, 365)
    ax.grid(True, alpha=0.3)

    month_starts = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]
    month_names = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                   "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    ax.set_xticks(month_starts)
    ax.set_xticklabels(month_names)

    path = DATA_DIR / "climatology_curves.png"
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  Plot saved: {path}")
    return path


def main():
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    t_start = time.time()

    print("=" * 60)
    print("  SST Data Extraction — 11 Stepping-Stone Nodes")
    print(f"  Period: {START_YEAR}–{END_YEAR} ({CHUNK_YEARS}-yr chunks)")
    print(f"  Workers: {CONCURRENT_WORKERS} concurrent")
    print(f"  Primary: NOAA OISST v2.1 | Fallback: MUR SST")
    print(f"  Output: {DATA_DIR}")
    print("=" * 60)

    # Fetch all nodes concurrently
    results = {}
    with ThreadPoolExecutor(max_workers=CONCURRENT_WORKERS) as pool:
        futures = {
            pool.submit(fetch_node_all, name, lat, lon): name
            for name, lat, lon in NODES
        }
        for future in as_completed(futures):
            name = futures[future]
            try:
                results[name] = future.result()
            except Exception as e:
                print(f"  [{name}] FATAL ERROR: {e}")
                results[name] = {
                    "name": name, "dates": [], "values": [],
                    "total_days": 0, "nan_count": 0, "nan_rate": 1.0,
                    "sst_min": float("nan"), "sst_max": float("nan"),
                    "sst_mean": float("nan"), "suspicious": True,
                    "source": "FAILED", "lat": 0, "lon": 0,
                    "used_lat": 0, "used_lon": 0, "elapsed_s": 0,
                }

    # Process results in node order
    all_clim = {}
    summaries = []

    for name, lat, lon in NODES:
        r = results[name]
        summaries.append({k: v for k, v in r.items() if k not in ("dates", "values")})

        if not r["dates"]:
            print(f"\n  ❌ NO DATA for {name}")
            continue

        # Save daily
        dp = save_daily_csv(name, r["dates"], r["values"])
        print(f"  Saved: {dp}")

        # Climatology
        clim = compute_climatology(r["dates"], r["values"])
        cp = save_climatology_csv(name, clim)
        print(f"  Saved: {cp}")
        all_clim[name] = clim

    # Combined climatology
    if all_clim:
        p = save_combined_climatology(all_clim)
        print(f"\n  Combined: {p}")

    # Plot
    if all_clim:
        plot_climatologies(all_clim)

    # Summary JSON
    summary_path = DATA_DIR / "fetch_summary.json"
    with open(summary_path, "w") as f:
        json.dump({
            "fetch_time": datetime.now().isoformat(),
            "period": f"{START_YEAR}-{END_YEAR}",
            "total_elapsed_s": time.time() - t_start,
            "nodes": summaries,
        }, f, indent=2, default=str)

    # Summary table
    elapsed = time.time() - t_start
    print(f"\n{'='*80}")
    print(f"  TOTAL TIME: {elapsed:.0f}s ({elapsed/60:.1f} min)")
    print(f"  {'Node':<15} {'Source':<18} {'Days':>6} {'NaN%':>6} {'Min':>6} {'Max':>6} {'Mean':>6}")
    print("-" * 80)
    for s in summaries:
        src = s['source'][:17]
        print(f"  {s['name']:<15} {src:<18} {s['total_days']:>6} {s['nan_rate']:>5.1%} "
              f"{s['sst_min']:>6.1f} {s['sst_max']:>6.1f} {s['sst_mean']:>6.1f}")
    print("=" * 80)

    failures = [s for s in summaries if s["total_days"] == 0 or s["nan_rate"] > 0.5]
    if failures:
        print(f"\n  ⚠️  {len(failures)} node(s) with issues:")
        for fn in failures:
            print(f"    - {fn['name']}: {fn['total_days']} days, {fn['nan_rate']:.1%} NaN")
    else:
        print(f"\n  ✅ All {len(summaries)} nodes fetched successfully!")

    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())
