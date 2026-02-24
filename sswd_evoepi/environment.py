"""Environmental forcing module.

Provides SST, salinity, and flushing-rate time series per spatial node.

SST forcing uses a sinusoidal annual cycle with:
  - Per-node mean SST and amplitude (latitude-dependent)
  - Linear warming trend (configurable per node)
  - Phase shift (peak in late summer, ~Aug 15 = day 227)

Salinity is constant per node (fjord nodes lower; open coast ~32 psu).

Flushing rate is constant or seasonally modulated for fjords.

References:
  - spatial-connectivity-spec.md §2, §6
  - CODE_ERRATA CE-6: σ_D reduced
  - Hood Canal φ ≈ 0.007 (ERRATA E3)

Build target: Phase 9 (Spatial Network).
"""

from __future__ import annotations

import os
import re
from typing import Optional

import numpy as np


# ═══════════════════════════════════════════════════════════════════════
# SST FORCING — sinusoidal annual cycle + linear trend
# ═══════════════════════════════════════════════════════════════════════

# Peak SST day-of-year (mid-August ≈ day 227)
_SST_PEAK_DOY = 227


def sinusoidal_sst(day_of_year: int, mean_sst: float,
                   amplitude: float, peak_doy: int = 227) -> float:
    """Compute SST from sinusoidal annual cycle.

    T(d) = T_mean + A × cos(2π × (d − d_peak) / 365)

    Args:
        day_of_year: 0-indexed day of year [0, 364].
        mean_sst: Annual mean SST (°C).
        amplitude: Half-range of annual cycle (°C).
        peak_doy: Day of year for SST maximum (default 227 ≈ Aug 15).

    Returns:
        SST (°C).
    """
    phase = 2.0 * np.pi * (day_of_year - peak_doy) / 365.0
    return mean_sst + amplitude * np.cos(phase)



def sst_with_trend(day_of_year: int, year: int, mean_sst: float,
                   amplitude: float, trend_per_year: float,
                   reference_year: int = 2000,
                   peak_doy: int = 227) -> float:
    """SST with sinusoidal cycle + linear warming trend.

    T(d, y) = T_mean + trend × (y − y_ref) + A × cos(2π × (d − d_peak) / 365)

    Args:
        day_of_year: 0-indexed [0, 364].
        year: Calendar year (e.g. 2013).
        mean_sst: Baseline annual mean SST (°C) at reference_year.
        amplitude: Annual cycle half-range (°C).
        trend_per_year: °C per year warming (e.g. 0.02).
        reference_year: Year at which mean_sst applies (default 2000).
        peak_doy: Day of peak SST.

    Returns:
        SST (°C).
    """
    base = mean_sst + trend_per_year * (year - reference_year)
    phase = 2.0 * np.pi * (day_of_year - peak_doy) / 365.0
    return base + amplitude * np.cos(phase)


def make_sst_timeseries(n_years: int, start_year: int,
                        mean_sst: float, amplitude: float,
                        trend_per_year: float = 0.02,
                        reference_year: int = 2000) -> np.ndarray:
    """Generate a daily SST time series for one node.

    Args:
        n_years: Number of years.
        start_year: First calendar year.
        mean_sst: Baseline mean SST (°C).
        amplitude: Annual cycle amplitude (°C).
        trend_per_year: Linear warming rate (°C/yr).
        reference_year: Year at which mean_sst applies.

    Returns:
        1-D array of shape (n_years * 365,) with daily SST values.
    """
    total_days = n_years * 365
    
    # Vectorized SST generation - create all day and year arrays at once
    day_array = np.tile(np.arange(365), n_years)  # [0,1,...,364,0,1,...,364,...]
    year_array = np.repeat(np.arange(start_year, start_year + n_years), 365)  # [2020,2020,...,2021,2021,...]
    
    # Vectorized SST calculation for all days at once
    sst = sst_with_trend(
        day_array, year_array, mean_sst, amplitude, trend_per_year, reference_year
    )
    return sst


# ═══════════════════════════════════════════════════════════════════════
# SST FORCING — satellite climatology
# ═══════════════════════════════════════════════════════════════════════


def _normalize_node_name(name: str) -> str:
    """Normalize a node name for climatology filename matching.

    Converts spaces to underscores, strips leading/trailing whitespace.
    Preserves original case (climatology files use Title_Case).

    Examples:
        "Howe Sound" → "Howe_Sound"
        "SJI"        → "SJI"
        " Fort Bragg " → "Fort_Bragg"
    """
    return re.sub(r'\s+', '_', name.strip())


def load_sst_climatology(node_name: str,
                         data_dir: str = 'data/sst') -> np.ndarray:
    """Load satellite-derived daily SST climatology for a node.

    Expects a CSV file at ``{data_dir}/{normalized_name}_climatology.csv``
    with columns ``day_of_year`` (1-indexed) and ``sst_mean``.

    Args:
        node_name: Node name (spaces and case are normalized).
        data_dir: Directory containing climatology CSV files.

    Returns:
        1-D array of shape (365,) with daily mean SST (°C).

    Raises:
        FileNotFoundError: If the climatology CSV does not exist.
        ValueError: If the file doesn't contain exactly 365 rows.
    """
    norm_name = _normalize_node_name(node_name)
    csv_path = os.path.join(data_dir, f"{norm_name}_climatology.csv")

    if not os.path.isfile(csv_path):
        raise FileNotFoundError(
            f"SST climatology file not found: {csv_path} "
            f"(node '{node_name}', normalized '{norm_name}')"
        )

    # Read CSV — expects columns: day_of_year, sst_mean, [sst_std, n_years]
    import csv
    sst_values = []
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            sst_values.append(float(row['sst_mean']))

    if len(sst_values) != 365:
        raise ValueError(
            f"SST climatology for '{node_name}' has {len(sst_values)} rows, "
            f"expected 365"
        )

    return np.array(sst_values, dtype=np.float64)


def load_sst_monthly(node_name: str,
                     data_dir: str = 'data/sst') -> dict:
    """Load monthly SST data for a node from *_monthly.csv.

    Expects a CSV file at ``{data_dir}/{normalized_name}_monthly.csv``
    with columns ``year``, ``month``, ``sst``.

    Args:
        node_name: Node name (spaces and case are normalized).
        data_dir: Directory containing monthly CSV files.

    Returns:
        Dict mapping (year, month) → SST value (°C).

    Raises:
        FileNotFoundError: If the monthly CSV does not exist.
    """
    norm_name = _normalize_node_name(node_name)
    csv_path = os.path.join(data_dir, f"{norm_name}_monthly.csv")

    if not os.path.isfile(csv_path):
        raise FileNotFoundError(
            f"SST monthly file not found: {csv_path} "
            f"(node '{node_name}', normalized '{norm_name}')"
        )

    import csv
    monthly = {}
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            yr = int(row['year'])
            mo = int(row['month'])
            monthly[(yr, mo)] = float(row['sst'])

    if len(monthly) == 0:
        raise ValueError(
            f"SST monthly file for '{node_name}' contains no data rows"
        )

    return monthly


def interpolate_monthly_to_daily(monthly_values_12: np.ndarray,
                                  year: int) -> np.ndarray:
    """Interpolate 12 monthly SST values to daily resolution.

    Places each monthly value at the midpoint of its month, then
    uses linear interpolation (with wrap-around for continuity)
    to produce a daily array.

    Args:
        monthly_values_12: Array of shape (12,) with monthly SST (°C),
            indexed 0=January .. 11=December.
        year: Calendar year (used to determine if leap year → 366 days).

    Returns:
        1-D array of shape (365,) or (366,) with daily SST values.
    """
    import calendar
    is_leap = calendar.isleap(year)
    n_days = 366 if is_leap else 365

    # Compute midpoint day-of-year (0-indexed) for each month
    month_midpoints = np.empty(12, dtype=np.float64)
    cumulative = 0
    for m in range(12):
        days_in_month = calendar.monthrange(year, m + 1)[1]
        month_midpoints[m] = cumulative + (days_in_month - 1) / 2.0
        cumulative += days_in_month

    # Wrap around for smooth interpolation at year boundaries:
    # prepend December of "previous year" and append January of "next year"
    extended_midpoints = np.empty(14, dtype=np.float64)
    extended_values = np.empty(14, dtype=np.float64)

    # December before: midpoint shifted back by n_days
    extended_midpoints[0] = month_midpoints[11] - n_days
    extended_values[0] = monthly_values_12[11]

    # The 12 months
    extended_midpoints[1:13] = month_midpoints
    extended_values[1:13] = monthly_values_12

    # January after: midpoint shifted forward by n_days
    extended_midpoints[13] = month_midpoints[0] + n_days
    extended_values[13] = monthly_values_12[0]

    # Linear interpolation to daily resolution
    days = np.arange(n_days, dtype=np.float64)
    daily_sst = np.interp(days, extended_midpoints, extended_values)

    return daily_sst


def generate_yearly_sst_series(
    node_name: str,
    start_year: int,
    n_years: int,
    data_dir: str = 'data/sst',
    climatology_dir: Optional[str] = None,
    trend_per_year: float = 0.0,
    reference_year: int = 2015,
) -> np.ndarray:
    """Generate a full daily SST array for multi-year simulation from monthly data.

    For each year in the range [start_year, start_year + n_years), loads
    the 12 monthly SST values and interpolates to daily resolution.

    For years beyond the available monthly data (past 2025), falls back
    to the satellite climatology with a logged warning. If no climatology
    is available either, the last available year of monthly data is
    repeated.

    Note: the model uses a fixed 365-day year internally. For leap years,
    the interpolated array is truncated to 365 days.

    Args:
        node_name: Node name for file lookup.
        start_year: First calendar year of the series.
        n_years: Number of years.
        data_dir: Directory containing ``*_monthly.csv`` files.
        climatology_dir: Directory for climatology fallback. If None,
            defaults to ``data_dir``.
        trend_per_year: Linear warming rate (°C/yr) applied on top of
            data (default 0.0 — raw observations).
        reference_year: Year at which trend contribution is zero.

    Returns:
        1-D array of shape (n_years * 365,) with daily SST values.
    """
    import logging
    logger = logging.getLogger(__name__)

    if climatology_dir is None:
        climatology_dir = data_dir

    # Load monthly data
    monthly_data = load_sst_monthly(node_name, data_dir)

    # Determine available year range from monthly data
    available_years = sorted(set(yr for yr, _ in monthly_data.keys()))
    max_data_year = max(available_years) if available_years else start_year

    # Try loading climatology for fallback
    climatology = None
    try:
        climatology = load_sst_climatology(node_name, climatology_dir)
    except (FileNotFoundError, ValueError):
        pass

    # Find last complete year of monthly data (all 12 months present)
    last_complete_monthly = None
    for yr in reversed(available_years):
        if all((yr, m) in monthly_data for m in range(1, 13)):
            last_complete_monthly = yr
            break

    sst_chunks = []

    for yr_offset in range(n_years):
        cal_year = start_year + yr_offset

        # Check if all 12 months available for this year
        has_all_months = all(
            (cal_year, m) in monthly_data for m in range(1, 13)
        )

        if has_all_months:
            # Use actual monthly data
            monthly_vals = np.array(
                [monthly_data[(cal_year, m)] for m in range(1, 13)],
                dtype=np.float64,
            )
            daily = interpolate_monthly_to_daily(monthly_vals, cal_year)
        elif climatology is not None:
            # Fall back to climatology
            if cal_year > max_data_year:
                logger.warning(
                    "SST monthly data for '%s' year %d not available; "
                    "falling back to climatology.",
                    node_name, cal_year,
                )
            daily = climatology.copy()  # already 365 days
        elif last_complete_monthly is not None:
            # No climatology — repeat last complete year
            logger.warning(
                "SST monthly data for '%s' year %d not available and no "
                "climatology found; repeating year %d.",
                node_name, cal_year, last_complete_monthly,
            )
            monthly_vals = np.array(
                [monthly_data[(last_complete_monthly, m)] for m in range(1, 13)],
                dtype=np.float64,
            )
            daily = interpolate_monthly_to_daily(
                monthly_vals, last_complete_monthly
            )
        else:
            raise ValueError(
                f"No SST data available for '{node_name}' year {cal_year}: "
                f"no monthly data, no climatology, no fallback."
            )

        # Truncate to 365 days (model uses fixed 365-day years)
        daily_365 = daily[:365]

        # Apply linear trend if requested
        if trend_per_year != 0.0:
            daily_365 = daily_365 + trend_per_year * (cal_year - reference_year)

        sst_chunks.append(daily_365)

    return np.concatenate(sst_chunks)


def generate_satellite_sst_series(
    n_years: int,
    start_year: int,
    node_name: str,
    trend_per_year: float = 0.0,
    reference_year: int = 2015,
    data_dir: str = 'data/sst',
) -> np.ndarray:
    """Generate a daily SST time series from satellite climatology.

    Loads the 365-day climatology for the node and tiles it across
    ``n_years``, optionally applying a linear warming trend.

    Same output interface as ``make_sst_timeseries()``: returns a 1-D
    array of shape ``(n_years * 365,)``.

    Args:
        n_years: Number of years.
        start_year: First calendar year.
        node_name: Node name (for climatology file lookup).
        trend_per_year: Linear warming rate (°C/yr) on top of climatology.
        reference_year: Year at which climatology applies unchanged.
        data_dir: Directory containing climatology CSV files.

    Returns:
        1-D array of shape (n_years * 365,) with daily SST values.
    """
    clim = load_sst_climatology(node_name, data_dir)  # shape (365,)

    # Tile the 365-day climatology across n_years
    sst = np.tile(clim, n_years)  # shape (n_years * 365,)

    # Apply linear warming trend if non-zero
    if trend_per_year != 0.0:
        year_array = np.repeat(
            np.arange(start_year, start_year + n_years), 365
        )
        sst = sst + trend_per_year * (year_array - reference_year)

    return sst


# ═══════════════════════════════════════════════════════════════════════
# SALINITY
# ═══════════════════════════════════════════════════════════════════════


def salinity_modifier(salinity: float, s_min: float = 10.0,
                      s_full: float = 28.0) -> float:
    """Vibrio viability modifier from salinity (S_sal).

    S_sal = clip( ((sal − s_min) / (s_full − s_min))², 0, 1 )

    Returns 0 at low salinity (< s_min), 1 at full marine (≥ s_full).
    """
    if salinity <= s_min:
        return 0.0
    if salinity >= s_full:
        return 1.0
    frac = (salinity - s_min) / (s_full - s_min)
    return frac * frac


# ═══════════════════════════════════════════════════════════════════════
# FLUSHING RATE
# ═══════════════════════════════════════════════════════════════════════


def seasonal_flushing(phi_annual: float, month: int,
                      is_fjord: bool) -> float:
    """Seasonally modulated flushing rate.

    Fjords: ±30% seasonal amplitude (peak Jun–Aug, trough Dec–Feb).
    Open coast: ±20% amplitude.

    Args:
        phi_annual: Mean annual flushing rate (d⁻¹).
        month: 0-indexed month [0, 11].
        is_fjord: True if fjord node.

    Returns:
        Flushing rate for this month (d⁻¹).
    """
    # month 5 = June peak
    if is_fjord:
        amplitude = 0.3
    else:
        amplitude = 0.2
    seasonal = 1.0 + amplitude * np.cos(2.0 * np.pi * (month - 5) / 12.0)
    return phi_annual * seasonal


# ═══════════════════════════════════════════════════════════════════════
# VBNC SIGMOID
# ═══════════════════════════════════════════════════════════════════════


def vbnc_fraction(T_celsius: float, T_vbnc: float = 12.0,
                  k_vbnc: float = 1.0) -> float:
    """Fraction of Vibrio that is culturable (non-VBNC).

    f_vbnc(T) = 1 / (1 + exp(−k × (T − T_vbnc)))

    Below T_vbnc, most Vibrio enters VBNC state (non-pathogenic).
    """
    x = -k_vbnc * (T_celsius - T_vbnc)
    if x > 50.0:
        return 0.0
    if x < -50.0:
        return 1.0
    return 1.0 / (1.0 + np.exp(x))


# ═══════════════════════════════════════════════════════════════════════
# COMBINED ENVIRONMENTAL STATE (convenience)
# ═══════════════════════════════════════════════════════════════════════

def compute_environmental_state(sst: float, salinity: float,
                                phi: float, is_fjord: bool,
                                month: int) -> dict:
    """Compute derived environmental quantities for one node.

    Returns dict with keys:
      sst, salinity, flushing_rate, S_sal, f_vbnc, vibrio_viability
    """
    phi_now = seasonal_flushing(phi, month, is_fjord)
    s_sal = salinity_modifier(salinity)
    f_v = vbnc_fraction(sst)
    return {
        'sst': sst,
        'salinity': salinity,
        'flushing_rate': phi_now,
        'S_sal': s_sal,
        'f_vbnc': f_v,
        'vibrio_viability': s_sal * f_v,
    }
