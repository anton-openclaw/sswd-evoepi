"""Unified result loading for SSWD-EvoEpi calibration runs.

Single source of truth for:
- Loading individual seed results (result_seed*.json)
- Loading multi-seed runs (mean across seeds)
- Loading sweeps (multiple configs)
- Region definitions and ordering
- Monthly snapshot NPZ loading
"""
from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

__all__ = [
    "REGION_ORDER",
    "COASTLINE_ORDER_SN",
    "SCORED_REGIONS",
    "site_to_region",
    "sites_to_regions",
    "SeedResult",
    "RunResult",
    "load_seed",
    "load_run",
    "load_sweep",
    "load_monthly_npz",
]

# Canonical north-to-south ordering of all 18 regions
# Coastline order (N→S): follows the coast from the Aleutians southward to Baja.
# Within Alaska, the order follows the coastline from the Aleutians east through
# the Western Gulf, then south through SE Alaska.
REGION_ORDER: List[str] = [
    "AK-AL", "AK-WG", "AK-OC", "AK-EG", "AK-PWS", "AK-FN", "AK-FS",
    "BC-N", "BC-C", "SS-N", "SS-S", "JDF", "WA-O", "OR",
    "CA-N", "CA-C", "CA-S", "BJ",
]

# Same ordering, south → north (for plots with S at bottom, N at top)
COASTLINE_ORDER_SN: List[str] = list(reversed(REGION_ORDER))

# The 8 regions with calibration targets
SCORED_REGIONS: List[str] = [
    "AK-PWS", "AK-FN", "AK-FS", "BC-N", "SS-S", "JDF", "OR", "CA-N",
]


def site_to_region(site_name: str) -> str:
    """Map a site name to its region.

    Site names follow two formats:
    - 2-part: 'BJ-001', 'JDF-023', 'OR-006' → region is first part
    - 3-part: 'AK-PWS-001', 'CA-N-042' → region is first two parts

    The rule: region = everything before the final numeric segment.

    >>> site_to_region('AK-PWS-001')
    'AK-PWS'
    >>> site_to_region('JDF-023')
    'JDF'
    >>> site_to_region('BJ-007')
    'BJ'
    >>> site_to_region('OR-006')
    'OR'
    """
    parts = str(site_name).split('-')
    # Strip the trailing numeric site ID
    if len(parts) >= 2 and parts[-1].isdigit():
        return '-'.join(parts[:-1])
    # Fallback: strip any trailing all-digit segments
    while len(parts) > 1 and parts[-1].isdigit():
        parts = parts[:-1]
    return '-'.join(parts)


def sites_to_regions(site_names) -> List[str]:
    """Map an array of site names to region names.

    Args:
        site_names: array-like of site name strings (e.g., from NPZ site_names)

    Returns:
        List of region names, one per site.
    """
    return [site_to_region(n) for n in site_names]


def _extract_seed_number(path: Path) -> int:
    """Extract seed number from filename like result_seed42.json."""
    m = re.search(r"result_seed(\d+)", path.stem)
    if m:
        return int(m.group(1))
    return 0


@dataclass
class SeedResult:
    """Typed wrapper around a single seed's result JSON."""

    seed: int
    wall_time: float
    rmsle: float
    region_recovery: Dict[str, float]
    region_details: Dict[str, dict]
    scoring: dict
    overall: dict
    arrival_timing: dict
    raw: dict = field(repr=False)

    @property
    def rmse_log(self) -> float:
        """Backward-compatible alias for rmsle."""
        return self.rmsle

    @classmethod
    def from_json(cls, data: dict) -> SeedResult:
        """Create a SeedResult from a parsed JSON dict."""
        scoring = data.get("scoring", {})
        return cls(
            seed=data.get("seed", 0),
            wall_time=data.get("wall_time_seconds", 0.0),
            rmsle=scoring.get("rmsle", scoring.get("rmse_log", float("inf"))),
            region_recovery=data.get("region_recovery", {}),
            region_details=data.get("region_details", {}),
            scoring=scoring,
            overall=data.get("overall", {}),
            arrival_timing=data.get("arrival_timing", {}),
            raw=data,
        )

    @classmethod
    def from_file(cls, path: Path) -> SeedResult:
        """Load a SeedResult from a JSON file path."""
        path = Path(path)
        with open(path) as f:
            data = json.load(f)
        # If seed not in JSON, try to extract from filename
        if "seed" not in data:
            data["seed"] = _extract_seed_number(path)
        return cls.from_json(data)

    def recovery(self, region: str) -> float:
        """Get recovery fraction for a region."""
        return self.region_recovery.get(region, 0.0)

    def yearly_pop(self, region: str) -> List[int]:
        """Get yearly total population for a region."""
        details = self.region_details.get(region, {})
        return details.get("yearly_totals", [])

    def yearly_trait(self, region: str, trait: str) -> List[float]:
        """Get yearly mean trait values for a region.

        Parameters
        ----------
        region : str
            Region name (e.g. 'AK-PWS').
        trait : str
            One of 'resistance', 'tolerance', 'recovery'.
        """
        if trait not in ("resistance", "tolerance", "recovery"):
            raise ValueError(
                f"trait must be 'resistance', 'tolerance', or 'recovery', got {trait!r}"
            )
        details = self.region_details.get(region, {})
        return details.get(f"yearly_mean_{trait}", [])


@dataclass
class RunResult:
    """Aggregated results across multiple seeds for one configuration."""

    config_name: str
    seeds: List[SeedResult]
    description: str = ""
    group: str = "default"

    @classmethod
    def from_dir(cls, run_dir: Path, config_name: str = "") -> RunResult:
        """Load all result_seed*.json files from a directory.

        Parameters
        ----------
        run_dir : Path
            Directory containing result_seed*.json files.
        config_name : str
            Name for this configuration. Defaults to directory name.
        """
        run_dir = Path(run_dir)
        if not config_name:
            config_name = run_dir.name

        files = sorted(run_dir.glob("result_seed*.json"), key=_extract_seed_number)
        seeds = [SeedResult.from_file(f) for f in files]
        return cls(config_name=config_name, seeds=seeds)

    def mean_rmsle(self) -> float:
        """Mean RMSLE across seeds."""
        if not self.seeds:
            return float("inf")
        return float(np.mean([s.rmsle for s in self.seeds]))

    def mean_rmse(self) -> float:
        """Backward-compatible alias for mean_rmsle."""
        return self.mean_rmsle()

    def mean_recovery(self, region: str) -> float:
        """Mean recovery fraction for a region across seeds."""
        if not self.seeds:
            return 0.0
        values = [s.recovery(region) for s in self.seeds]
        return float(np.mean(values))

    def recovery_std(self, region: str) -> float:
        """Standard deviation of recovery fraction across seeds."""
        if len(self.seeds) < 2:
            return 0.0
        values = [s.recovery(region) for s in self.seeds]
        return float(np.std(values, ddof=1))

    def n_seeds(self) -> int:
        """Number of seeds loaded."""
        return len(self.seeds)


# ---------------------------------------------------------------------------
# Convenience functions
# ---------------------------------------------------------------------------

def load_seed(path: Path) -> SeedResult:
    """Load a single seed result from a JSON file."""
    return SeedResult.from_file(path)


def load_run(run_dir: Path, config_name: str = "") -> RunResult:
    """Load all seed results from a run directory."""
    return RunResult.from_dir(run_dir, config_name=config_name)


def _natural_sort_key(name: str):
    """Sort key that handles numeric parts naturally (W1 < W10 < W100)."""
    parts = re.split(r"(\d+)", name)
    return [int(p) if p.isdigit() else p.lower() for p in parts]


def load_sweep(
    sweep_dir: Path,
    configs: Optional[List[str]] = None,
) -> Dict[str, RunResult]:
    """Load multiple configurations from a sweep directory.

    Parameters
    ----------
    sweep_dir : Path
        Parent directory containing subdirectories (e.g. W1/, W2/, ...).
    configs : list of str, optional
        Specific subdirectory names to load. If None, loads all W* dirs.

    Returns
    -------
    dict
        OrderedDict-like mapping config_name -> RunResult, sorted by mean RMSE.
    """
    sweep_dir = Path(sweep_dir)

    if configs is not None:
        dirs = [(sweep_dir / c) for c in configs if (sweep_dir / c).is_dir()]
    else:
        dirs = sorted(
            [d for d in sweep_dir.iterdir() if d.is_dir() and d.name.startswith("W")],
            key=lambda d: _natural_sort_key(d.name),
        )

    runs: Dict[str, RunResult] = {}
    for d in dirs:
        run = RunResult.from_dir(d)
        if run.seeds:  # skip empty directories
            runs[run.config_name] = run

    # Sort by mean RMSE (best first)
    sorted_runs = dict(sorted(runs.items(), key=lambda kv: kv[1].mean_rmse()))
    return sorted_runs


def load_monthly_npz(npz_path: Path) -> dict:
    """Load a monthly snapshot NPZ file.

    Returns a dict with:
    - region_name -> population array for each region found
    - '_sim_days' -> simulation days array (if present)
    - '_years' -> year labels (if present)
    """
    npz_path = Path(npz_path)
    data = np.load(npz_path, allow_pickle=True)

    result: dict = {}
    for key in data.files:
        if key == "sim_days":
            result["_sim_days"] = data[key]
        elif key == "years":
            result["_years"] = data[key]
        else:
            # Map array keys to regions using REGION_ORDER if they are
            # integer-indexed, otherwise use the key as-is
            result[key] = data[key]

    return result
