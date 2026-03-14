"""Seasonal salinity model for freshwater influence.

Computes daily salinity at each node as a function of latitude,
fjord enclosedness, and a cosine meltwater pulse peaking mid-June.

Model:
    S(site, day) = max(s_floor,
        S_ocean(lat) - fw_strength × fd_norm × f_melt(lat) × melt_pulse(day))

Where:
    S_ocean(lat) = 31.32 + 0.054 × (lat - 50.0)
    fw_strength  — single calibration parameter (psu), default 0.0 (OFF)
    fd_norm      — fjord_depth_norm from site_enclosedness.csv, [0, 1]
    f_melt(lat)  = clip((lat - 35) / 25, 0, 1)
    melt_pulse(day) = max(0, cos(2π × (day - 166) / 365))
    s_floor      = 5.0 psu minimum

When fw_strength=0 (default), salinity is exactly the ocean baseline
for each node — backward compatible.
"""

from __future__ import annotations

import math
from typing import List, TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from sswd_evoepi.spatial import NodeDefinition

DAYS_PER_YEAR = 365
_PEAK_DAY = 166  # ~June 15


def ocean_baseline(lat: float) -> float:
    """Open-ocean baseline salinity from DFO lighthouse regression.

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


def compute_salinity_array(
    nodes: List["NodeDefinition"],
    fw_strength: float,
    fw_depth_exp: float = 1.0,
    s_floor: float = 5.0,
) -> np.ndarray:
    """Pre-compute daily salinity for all nodes over one year.

    When fw_strength=0, every node gets its ocean baseline for all 365 days
    (exact backward compatibility: no depression, no floor clipping unless
    the baseline itself is below s_floor which cannot happen for realistic lats).

    Args:
        nodes: List of NodeDefinition objects (need .lat, .fjord_depth_norm).
        fw_strength: Freshwater strength parameter (psu). 0 = mechanism OFF.
        fw_depth_exp: Exponent applied to fjord_depth_norm before use.
            1.0 = linear (default), 0.5 = sqrt (boosts moderate-depth sites).
        s_floor: Minimum physically reasonable salinity (psu).

    Returns:
        (N, 365) float32 array of daily salinity values.
    """
    N = len(nodes)
    salinity = np.empty((N, DAYS_PER_YEAR), dtype=np.float32)

    # Pre-compute daily melt pulse (shared across all nodes)
    pulse = np.array(
        [freshwater_melt_pulse(d) for d in range(DAYS_PER_YEAR)],
        dtype=np.float64,
    )

    for i, nd in enumerate(nodes):
        base = ocean_baseline(nd.lat)
        if fw_strength == 0.0:
            # Exact backward compatibility: constant baseline, no floor clipping
            salinity[i, :] = base
        else:
            fd_norm = getattr(nd, 'fjord_depth_norm', 0.0)
            fd_effective = fd_norm ** fw_depth_exp if fd_norm > 0.0 else 0.0
            f_melt = latitude_melt_factor(nd.lat)
            depression = fw_strength * fd_effective * f_melt * pulse
            daily = np.maximum(s_floor, base - depression)
            salinity[i, :] = daily

    return salinity
