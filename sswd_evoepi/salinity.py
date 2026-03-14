"""Seasonal salinity model for freshwater influence.

Two-layer model:
  Layer 1: WOA23 monthly climatology (0.25° grid) — captures real seasonal
           cycles at each site via nearest-ocean-cell interpolation.
  Layer 2: Fjord freshwater depression — fw_strength × fd_norm^exp adds
           extra freshwater that WOA23 can't resolve at coarse resolution.

Model:
    S(site, day) = max(s_floor,
        S_woa23(site, month(day)) - fw_strength × fd_norm^exp × pulse(day))

Where:
    S_woa23    — WOA23 monthly surface salinity at nearest ocean grid cell
    fw_strength — calibration parameter (psu), default 0.0 (OFF)
    fd_norm    — fjord_depth_norm from site_enclosedness.csv, [0, 1]
    fw_depth_exp — exponent on fd_norm, default 1.0
    pulse(day) — cosine seasonal pulse, peaks at day 166 (~June 15)
    s_floor    = 5.0 psu minimum

Fallback: if WOA23 data is not available, falls back to the parametric
baseline S_ocean(lat) = 31.32 + 0.054 × (lat - 50.0).

When fw_strength=0 (default), salinity equals the WOA23 baseline
(or parametric fallback) — backward compatible.
"""

from __future__ import annotations

import math
import os
from pathlib import Path
from typing import List, Optional, TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from sswd_evoepi.spatial import NodeDefinition

DAYS_PER_YEAR = 365
_PEAK_DAY = 166  # ~June 15

# Month boundaries (day-of-year, 0-indexed)
_MONTH_STARTS = np.array([0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334])
_MONTH_ENDS = np.array([31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365])

# WOA23 data paths (searched in order)
_WOA23_SEARCH_PATHS = [
    Path(__file__).resolve().parents[1] / 'data' / 'salinity' / 'woa23_surface_nepac.npz',
    Path(__file__).resolve().parents[2] / 'salinity_validation' / 'woa23_surface_nepac.npz',
    Path(os.environ.get('WOA23_SALINITY_PATH', '/dev/null')),
]

# ─── WOA23 data cache ────────────────────────────────────────────────
_woa23_cache: Optional[dict] = None
_woa23_load_attempted: bool = False


def _load_woa23() -> Optional[dict]:
    """Lazy-load WOA23 surface salinity data. Returns dict or None."""
    global _woa23_cache, _woa23_load_attempted
    if _woa23_load_attempted:
        return _woa23_cache
    _woa23_load_attempted = True

    for path in _WOA23_SEARCH_PATHS:
        if path.exists() and path.stat().st_size > 0:
            try:
                d = np.load(str(path))
                _woa23_cache = {
                    'lat': d['lat'],
                    'lon': d['lon'],
                    'salinity': d['salinity_filled'],  # (12, lat, lon) — land cells NN-filled
                }
                return _woa23_cache
            except Exception:
                continue
    return None


def _extract_woa23_monthly(lat: float, lon: float, woa: dict) -> np.ndarray:
    """Extract 12-month salinity cycle for a single site from WOA23.

    Uses nearest-grid-cell lookup (0.25° grid → max ~18km error).

    Args:
        lat: Site latitude (°N).
        lon: Site longitude (°E, west negative).
        woa: WOA23 data dict with 'lat', 'lon', 'salinity' keys.

    Returns:
        (12,) array of monthly salinity values (psu).
    """
    lat_i = int(np.argmin(np.abs(woa['lat'] - lat)))
    lon_i = int(np.argmin(np.abs(woa['lon'] - lon)))
    return woa['salinity'][:, lat_i, lon_i].astype(np.float64)


# ─── Parametric fallback (preserved for backward compat & testing) ───

def ocean_baseline(lat: float) -> float:
    """Open-ocean baseline salinity from DFO lighthouse regression.

    Used as fallback when WOA23 data is unavailable.

    Args:
        lat: Latitude in decimal degrees N.

    Returns:
        Baseline salinity (psu).
    """
    return 31.32 + 0.054 * (lat - 50.0)


def freshwater_melt_pulse(day_of_year: int) -> float:
    """Cosine meltwater pulse peaking at day 166 (~June 15).

    Returns 0 during winter (when cosine is negative), and peaks at 1.0
    on day 166.

    Args:
        day_of_year: Day of year (0-indexed, 0 = Jan 1).

    Returns:
        Pulse strength in [0, 1].
    """
    return max(0.0, math.cos(2.0 * math.pi * (day_of_year - _PEAK_DAY) / DAYS_PER_YEAR))


def latitude_melt_factor(lat: float) -> float:
    """Latitude-dependent melt amplitude.

    0 at 35°N (no glacial melt), 1 at 60°N (full glacial influence).

    Args:
        lat: Latitude in decimal degrees N.

    Returns:
        Melt factor in [0, 1].
    """
    return max(0.0, min(1.0, (lat - 35.0) / 25.0))


# ─── Daily interpolation from monthly WOA23 ─────────────────────────

def _monthly_to_daily(monthly: np.ndarray) -> np.ndarray:
    """Interpolate 12 monthly values to 365 daily values.

    Uses linear interpolation between month midpoints to avoid
    step discontinuities at month boundaries.

    Args:
        monthly: (12,) array of monthly values.

    Returns:
        (365,) array of daily values.
    """
    # Month midpoints (day of year, 0-indexed)
    midpoints = (_MONTH_STARTS + _MONTH_ENDS) / 2.0  # [15.5, 45.0, ...]
    days = np.arange(DAYS_PER_YEAR, dtype=np.float64)

    # Wrap for interpolation across Dec→Jan boundary
    mid_ext = np.concatenate([midpoints[-1:] - DAYS_PER_YEAR, midpoints,
                              midpoints[:1] + DAYS_PER_YEAR])
    val_ext = np.concatenate([monthly[-1:], monthly, monthly[:1]])

    return np.interp(days, mid_ext, val_ext).astype(np.float32)


# ─── Main entry point ────────────────────────────────────────────────

def compute_salinity_array(
    nodes: List["NodeDefinition"],
    fw_strength: float,
    fw_depth_exp: float = 1.0,
    s_floor: float = 5.0,
    woa23_path: Optional[str] = None,
) -> np.ndarray:
    """Pre-compute daily salinity for all nodes over one year.

    Two-layer model:
      1. WOA23 monthly climatology → daily interpolation (or parametric fallback)
      2. Fjord freshwater depression: fw_strength × fd_norm^exp × pulse(day)

    When fw_strength=0, every node gets its WOA23 baseline for all 365 days.
    When WOA23 is unavailable, falls back to parametric ocean_baseline(lat).

    Args:
        nodes: List of NodeDefinition objects (need .lat, .lon, .fjord_depth_norm).
        fw_strength: Freshwater strength parameter (psu). 0 = no fjord depression.
        fw_depth_exp: Exponent applied to fjord_depth_norm before use.
            1.0 = linear (default), 0.5 = sqrt (boosts moderate-depth sites).
        s_floor: Minimum physically reasonable salinity (psu).
        woa23_path: Optional explicit path to WOA23 npz file (overrides search).

    Returns:
        (N, 365) float32 array of daily salinity values.
    """
    N = len(nodes)
    salinity = np.empty((N, DAYS_PER_YEAR), dtype=np.float32)

    # Try loading WOA23
    woa = None
    if woa23_path is not None:
        # Explicit path given — use it or fall back to parametric (no global cache)
        try:
            d = np.load(woa23_path)
            woa = {'lat': d['lat'], 'lon': d['lon'], 'salinity': d['salinity_filled']}
        except Exception:
            pass  # woa stays None → parametric fallback
    else:
        # Default: try global WOA23 cache
        woa = _load_woa23()

    using_woa23 = woa is not None

    # Pre-compute daily pulse (for fjord depression layer)
    pulse = np.array(
        [freshwater_melt_pulse(d) for d in range(DAYS_PER_YEAR)],
        dtype=np.float64,
    )

    for i, nd in enumerate(nodes):
        # Layer 1: Baseline seasonal cycle
        if using_woa23:
            monthly = _extract_woa23_monthly(nd.lat, nd.lon, woa)
            baseline = _monthly_to_daily(monthly)
        else:
            # Parametric fallback: constant across all days
            baseline = np.full(DAYS_PER_YEAR, ocean_baseline(nd.lat), dtype=np.float32)

        if fw_strength == 0.0:
            # No fjord depression — just the baseline
            salinity[i, :] = baseline
        else:
            # Layer 2: Fjord freshwater depression
            # Depression scales with: fjord depth × latitude melt factor × seasonal pulse
            fd_norm = getattr(nd, 'fjord_depth_norm', 0.0)
            fd_effective = fd_norm ** fw_depth_exp if fd_norm > 0.0 else 0.0
            f_melt = latitude_melt_factor(nd.lat)
            depression = fw_strength * fd_effective * f_melt * pulse
            daily = np.maximum(s_floor, baseline - depression)
            salinity[i, :] = daily

    return salinity
