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
