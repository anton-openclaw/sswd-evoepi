#!/usr/bin/env python3
"""
Generate figure: Pathogen Thermal Tolerance Evolution
Shows T_vbnc_local adaptation trajectories by region alongside SST context.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from pathlib import Path
from datetime import datetime, timedelta
from scipy.interpolate import interp1d

# ── Paths ──────────────────────────────────────────────────────────────────
BASE = Path(__file__).resolve().parent.parent.parent.parent  # sswd-evoepi root
SST_DIR = BASE / "data" / "sst"
SNAP_PATH = BASE / "results" / "calibration" / "W142" / "monthly_seed123.npz"
OUT_PATH = Path(__file__).resolve().parent / "fig_thermal_adaptation.png"

# ── Model parameters ──────────────────────────────────────────────────────
T_VBNC_INITIAL = 12.0   # Starting threshold °C
T_VBNC_MIN = 9.0        # Hard floor °C
ADAPT_RATE = 0.001       # °C/day per unit selection pressure
P_ADAPT_HALF = 500.0     # Half-saturation for pool-driven adaptation
REVERT_RATE = 0.0        # Reversion speed (disabled)
P_ENV_MAX = 2000.0       # Max environmental pathogen pool
K = 5000                 # Carrying capacity

# ── Region definitions ────────────────────────────────────────────────────
# Region label → (SST source file, site prefix patterns)
REGIONS = {
    "Alaska (Ketchikan)":       ("Ketchikan",    ["AK-"]),
    "BC North (Bella Bella)":   ("Bella_Bella",  ["BC-"]),
    "Salish Sea (SJI)":         ("SJI",          ["SS-", "JDF-"]),
    "Oregon (Newport)":         ("Newport",      ["OR-", "WA-"]),
    "N. California (Crescent City)": ("Crescent_City", ["CA-N"]),
    "S. California (Monterey)": ("Monterey",     ["CA-S", "CA-C"]),
}

# Colors: cool (north) → warm (south)
COLORS = {
    "Alaska (Ketchikan)":       "#2166ac",   # deep blue
    "BC North (Bella Bella)":   "#4393c3",   # medium blue
    "Salish Sea (SJI)":         "#92c5de",   # light blue
    "Oregon (Newport)":         "#f4a582",   # light salmon
    "N. California (Crescent City)": "#d6604d",  # salmon red
    "S. California (Monterey)": "#b2182b",   # deep red
}


def load_sst(source_name):
    """Load monthly SST for a source, return DataFrame with year, month, sst."""
    path = SST_DIR / f"{source_name}_monthly.csv"
    df = pd.read_csv(path)
    return df[["year", "month", "sst"]].dropna()


def monthly_sst_to_daily(df_sst, start_year=2012, end_year=2025):
    """Interpolate monthly SST to daily resolution via linear interp at month midpoints."""
    # Build month-midpoint dates and SST values
    dates = []
    temps = []
    for _, row in df_sst.iterrows():
        yr, mo, sst = int(row["year"]), int(row["month"]), row["sst"]
        if yr < start_year or yr > end_year:
            continue
        # Midpoint of each month (day 15)
        dates.append(datetime(yr, mo, 15))
        temps.append(sst)

    if not dates:
        return None, None

    # Convert to ordinal for interpolation
    ordinals = np.array([d.toordinal() for d in dates])
    temps = np.array(temps)

    # Daily date range
    d0 = datetime(start_year, 1, 1)
    d1 = datetime(end_year, 12, 31)
    daily_ordinals = np.arange(d0.toordinal(), d1.toordinal() + 1)

    # Interpolate (extrapolate at edges)
    f = interp1d(ordinals, temps, kind="linear", fill_value="extrapolate")
    daily_sst = f(daily_ordinals)

    daily_dates = [datetime.fromordinal(int(o)) for o in daily_ordinals]
    return daily_dates, daily_sst


def get_region_site_indices(site_names, prefixes):
    """Get indices of sites matching any prefix."""
    indices = []
    for i, name in enumerate(site_names):
        if any(name.startswith(p) for p in prefixes):
            indices.append(i)
    return np.array(indices)


def compute_monthly_pathogen_pool(infected, populations, site_indices, K):
    """
    Compute regional P_env_pool estimate from monthly snapshot data.
    Returns array of shape (n_months,) with pool values scaled to 0-P_ENV_MAX.
    """
    if len(site_indices) == 0:
        return np.zeros(infected.shape[0])

    # Sum infected across all sites in region
    region_infected = infected[:, site_indices].sum(axis=1).astype(float)
    region_pop = populations[:, site_indices].sum(axis=1).astype(float)

    # Fraction infected → scale to pathogen pool
    frac_infected = np.where(region_pop > 0, region_infected / region_pop, 0.0)
    # Scale: at 100% infection, pool = P_ENV_MAX
    pool = frac_infected * P_ENV_MAX
    return pool


def interpolate_monthly_to_daily(sim_days, monthly_values, total_days):
    """Interpolate monthly snapshot values to daily resolution."""
    f = interp1d(sim_days, monthly_values, kind="linear",
                 bounds_error=False, fill_value=(monthly_values[0], monthly_values[-1]))
    return f(np.arange(total_days))


def run_adaptation(daily_sst, daily_pool, total_days):
    """Run the thermal adaptation algorithm forward, return daily T_vbnc_local."""
    T_local = np.full(total_days, T_VBNC_INITIAL)

    for day in range(1, total_days):
        T_local[day] = T_local[day - 1]
        T_celsius = daily_sst[day] if day < len(daily_sst) else daily_sst[-1]
        P_env = daily_pool[day] if day < len(daily_pool) else 0.0

        if P_env > 0 and T_celsius < T_local[day]:
            pool_factor = min(P_env / max(P_ADAPT_HALF, 1e-9), 1.0)
            temp_gap = T_local[day] - T_celsius
            delta_T = ADAPT_RATE * pool_factor * temp_gap
            T_local[day] = max(T_local[day] - delta_T, T_VBNC_MIN)
        elif P_env == 0 and T_local[day] < T_VBNC_INITIAL:
            T_local[day] = min(T_local[day] + REVERT_RATE, T_VBNC_INITIAL)

    return T_local


def compute_mean_annual_sst(daily_dates, daily_sst):
    """Compute mean annual SST from daily data."""
    df = pd.DataFrame({"date": daily_dates, "sst": daily_sst})
    df["year"] = df["date"].apply(lambda d: d.year)
    return df.groupby("year")["sst"].mean()


def main():
    # ── Load snapshot data ────────────────────────────────────────────────
    snap = np.load(SNAP_PATH, allow_pickle=True)
    sim_days = snap["sim_days"]
    infected = snap["infected"]
    populations = snap["populations"]
    site_names = snap["site_names"]
    sst_start_year = int(snap["sst_start_year"])

    # Simulation time span
    total_sim_days = int(sim_days[-1]) + 1
    sim_start = datetime(sst_start_year, 1, 1)

    # ── Process each region ───────────────────────────────────────────────
    results = {}

    for region_label, (sst_source, prefixes) in REGIONS.items():
        # Load & interpolate SST
        df_sst = load_sst(sst_source)
        daily_dates, daily_sst = monthly_sst_to_daily(df_sst, sst_start_year, 2025)
        if daily_dates is None:
            print(f"WARNING: No SST data for {region_label}, skipping")
            continue

        # Get site indices for this region
        site_idx = get_region_site_indices(site_names, prefixes)
        print(f"{region_label}: {len(site_idx)} sites, SST source={sst_source}")

        # Compute monthly pathogen pool
        monthly_pool = compute_monthly_pathogen_pool(infected, populations, site_idx, K)

        # Interpolate pool to daily
        daily_pool = interpolate_monthly_to_daily(sim_days, monthly_pool, total_sim_days)

        # Trim SST to simulation length
        sst_for_sim = np.array(daily_sst[:total_sim_days])
        if len(sst_for_sim) < total_sim_days:
            sst_for_sim = np.pad(sst_for_sim, (0, total_sim_days - len(sst_for_sim)),
                                 mode="edge")

        # Run adaptation
        T_local = run_adaptation(sst_for_sim, daily_pool, total_sim_days)

        # Build date array for simulation
        sim_dates = [sim_start + timedelta(days=int(d)) for d in range(total_sim_days)]

        # Mean annual SST
        mean_annual = compute_mean_annual_sst(daily_dates, daily_sst)

        results[region_label] = {
            "sim_dates": sim_dates,
            "T_local": T_local,
            "daily_dates": daily_dates,
            "daily_sst": daily_sst,
            "mean_annual_sst": mean_annual,
        }

    # ── Plot ──────────────────────────────────────────────────────────────
    plt.style.use("seaborn-v0_8-whitegrid")
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), dpi=150,
                                    sharex=True, height_ratios=[1.2, 1])

    # ── Top panel: T_vbnc_local trajectories ──────────────────────────────
    for region_label in REGIONS:
        if region_label not in results:
            continue
        r = results[region_label]
        ax1.plot(r["sim_dates"], r["T_local"],
                 color=COLORS[region_label], linewidth=1.8, label=region_label,
                 alpha=0.9)

    # Reference lines
    ax1.axhline(y=T_VBNC_INITIAL, color="gray", linestyle="--", linewidth=1.0, alpha=0.7)
    ax1.text(datetime(2012, 2, 1), T_VBNC_INITIAL + 0.15, "Initial threshold (12°C)",
             fontsize=9, color="gray", style="italic")

    ax1.axhline(y=T_VBNC_MIN, color="gray", linestyle="--", linewidth=1.0, alpha=0.7)
    ax1.text(datetime(2012, 2, 1), T_VBNC_MIN - 0.35, "Biophysical floor (9°C)",
             fontsize=9, color="gray", style="italic")

    # Disease introduction line
    ax1.axvline(x=datetime(2013, 6, 1), color="#555555", linestyle=":", linewidth=1.2, alpha=0.6)
    ax1.text(datetime(2013, 7, 1), T_VBNC_INITIAL + 0.5,
             "Disease\nintroduction", fontsize=8, color="#555555",
             ha="left", va="bottom")

    ax1.set_ylabel("Thermal Activation Threshold (°C)", fontsize=12, fontweight="bold")
    ax1.set_ylim(8.0, 13.0)
    ax1.legend(loc="upper right", fontsize=9, framealpha=0.9, ncol=2)
    ax1.set_title("Pathogen Cold-Adaptation: Thermal Threshold Evolution",
                   fontsize=14, fontweight="bold", pad=12)
    ax1.tick_params(axis="both", labelsize=10)

    # ── Bottom panel: Mean annual SST ─────────────────────────────────────
    for region_label in REGIONS:
        if region_label not in results:
            continue
        r = results[region_label]
        annual = r["mean_annual_sst"]
        # Filter to simulation years
        annual = annual[(annual.index >= sst_start_year) & (annual.index <= 2024)]
        years_dt = [datetime(int(yr), 7, 1) for yr in annual.index]
        ax2.plot(years_dt, annual.values,
                 color=COLORS[region_label], linewidth=2.0, marker="o",
                 markersize=5, label=region_label, alpha=0.9)

    # Reference line at T_vbnc_min
    ax2.axhline(y=T_VBNC_MIN, color="gray", linestyle="--", linewidth=1.0, alpha=0.7)
    ax2.text(datetime(2012, 2, 1), T_VBNC_MIN - 0.3, "Biophysical floor (9°C)",
             fontsize=9, color="gray", style="italic")

    ax2.set_ylabel("Mean Annual SST (°C)", fontsize=12, fontweight="bold")
    ax2.set_xlabel("Year", fontsize=12, fontweight="bold")
    ax2.legend(loc="upper left", fontsize=9, framealpha=0.9, ncol=2)
    ax2.tick_params(axis="both", labelsize=10)

    # X-axis formatting
    ax2.xaxis.set_major_locator(mdates.YearLocator())
    ax2.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax2.set_xlim(datetime(2012, 1, 1), datetime(2025, 1, 1))
    plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45, ha="right")

    # ── Layout & save ─────────────────────────────────────────────────────
    fig.tight_layout(h_pad=1.5)

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PATH, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"\nSaved: {OUT_PATH}")


if __name__ == "__main__":
    main()
