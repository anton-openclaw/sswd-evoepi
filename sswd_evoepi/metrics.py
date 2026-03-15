"""Unified scoring and comparison for SSWD-EvoEpi calibration.

Single source of truth for:
- Recovery targets (the AUTHORITATIVE values)
- Arrival timing targets
- RMSE computation
- Run comparison and formatting
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Dict, List, Optional

from .results import SCORED_REGIONS, RunResult, SeedResult, load_run

__all__ = [
    "RECOVERY_TARGETS",
    "ARRIVAL_TARGETS",
    "rmse_log",
    "score_regions",
    "score_run",
    "score_seed",
    "compare_runs",
    "format_comparison_table",
    "quick_score",
]

# ---------------------------------------------------------------------------
# Authoritative calibration targets
# ---------------------------------------------------------------------------

RECOVERY_TARGETS: Dict[str, float] = {
    "AK-PWS": 0.30,
    "AK-FN": 0.30,
    "AK-FS": 0.20,
    "BC-N": 0.20,
    "SS-S": 0.05,
    "JDF": 0.02,
    "OR": 0.0025,
    "CA-N": 0.001,
}

ARRIVAL_TARGETS: Dict[str, int] = {
    "CA-S": 0,
    "CA-C": 6,
    "CA-N": 6,
    "OR": 15,
    "WA-O": 15,
    "JDF": 26,
    "SS-S": 26,
    "SS-N": 26,
    "BC-C": 26,
    "BC-N": 26,
    "AK-FS": 26,
    "AK-FN": 26,
    "AK-PWS": 33,
    "AK-EG": 33,
    "AK-OC": 33,
    "AK-WG": 42,
    "AK-AL": 42,
}

# Small epsilon to avoid log(0)
_EPS = 1e-6


# ---------------------------------------------------------------------------
# Pure scoring functions (no dependency on dataclasses)
# ---------------------------------------------------------------------------

def rmse_log(
    actual: Dict[str, float],
    targets: Optional[Dict[str, float]] = None,
) -> float:
    """Compute log-space RMSE between actual and target recovery fractions.

    Parameters
    ----------
    actual : dict
        Region name -> actual recovery fraction.
    targets : dict, optional
        Region name -> target recovery fraction. Defaults to RECOVERY_TARGETS.

    Returns
    -------
    float
        Root mean squared error in log space.
    """
    if targets is None:
        targets = RECOVERY_TARGETS

    sq_errors: list[float] = []
    for region, target in targets.items():
        act = actual.get(region, 0.0)
        log_err = math.log(max(act, _EPS)) - math.log(max(target, _EPS))
        sq_errors.append(log_err * log_err)

    if not sq_errors:
        return float("inf")

    return math.sqrt(sum(sq_errors) / len(sq_errors))


def score_regions(
    actual: Dict[str, float],
    targets: Optional[Dict[str, float]] = None,
) -> Dict[str, dict]:
    """Per-region scoring with grade strings.

    Parameters
    ----------
    actual : dict
        Region name -> actual recovery fraction.
    targets : dict, optional
        Region name -> target recovery fraction. Defaults to RECOVERY_TARGETS.

    Returns
    -------
    dict
        Region name -> {target, actual, log_error, within_2x, within_5x, grade_str}.
    """
    if targets is None:
        targets = RECOVERY_TARGETS

    results: Dict[str, dict] = {}
    for region, target in targets.items():
        act = actual.get(region, 0.0)
        log_act = math.log(max(act, _EPS))
        log_tgt = math.log(max(target, _EPS))
        log_error = log_act - log_tgt

        ratio = max(act, _EPS) / max(target, _EPS)
        fold = max(ratio, 1.0 / ratio) if ratio > 0 else float("inf")

        w2x = fold <= 2.0
        w5x = fold <= 5.0

        if w2x:
            grade = "✓ 2×"
        elif w5x:
            grade = "~ 5×"
        else:
            grade = f"✗ {fold:.0f}×"

        results[region] = {
            "target": target,
            "actual": act,
            "log_error": log_error,
            "within_2x": w2x,
            "within_5x": w5x,
            "grade_str": grade,
        }

    return results


# ---------------------------------------------------------------------------
# Dataclass-aware scoring
# ---------------------------------------------------------------------------

def _build_score_dict(
    actual: Dict[str, float],
    targets: Optional[Dict[str, float]] = None,
) -> dict:
    """Internal helper: build a full score dictionary from recovery values."""
    if targets is None:
        targets = RECOVERY_TARGETS

    per_region = score_regions(actual, targets)
    rmse = rmse_log(actual, targets)
    n_2x = sum(1 for v in per_region.values() if v["within_2x"])
    n_5x = sum(1 for v in per_region.values() if v["within_5x"])

    return {
        "rmse_log": rmse,
        "per_region": per_region,
        "within_2x": n_2x,
        "within_5x": n_5x,
        "n_targets": len(targets),
    }


def score_seed(seed: SeedResult) -> dict:
    """Score a single SeedResult against recovery targets.

    Returns
    -------
    dict
        {rmse_log, per_region, within_2x, within_5x, n_targets}
    """
    return _build_score_dict(seed.region_recovery)


def score_run(run: RunResult) -> dict:
    """Score a RunResult (mean recovery across seeds) against targets.

    Returns
    -------
    dict
        {rmse_log, per_region, within_2x, within_5x, n_targets}
    """
    if not run.seeds:
        return _build_score_dict({})

    # Compute mean recovery across seeds for all regions
    all_regions: set[str] = set()
    for s in run.seeds:
        all_regions.update(s.region_recovery.keys())

    mean_recovery: Dict[str, float] = {}
    for region in all_regions:
        mean_recovery[region] = run.mean_recovery(region)

    return _build_score_dict(mean_recovery)


# ---------------------------------------------------------------------------
# Multi-run comparison
# ---------------------------------------------------------------------------

def compare_runs(
    runs: Dict[str, RunResult],
    sort_by: str = "rmse",
) -> List[dict]:
    """Compare multiple runs and return sorted results.

    Parameters
    ----------
    runs : dict
        Config name -> RunResult.
    sort_by : str
        Sort metric: 'rmse' (default) or 'within_2x'.

    Returns
    -------
    list of dict
        Each dict has: config, rmse, recovery, within_2x, within_5x, description.
    """
    comparisons: list[dict] = []

    for config, run in runs.items():
        scored = score_run(run)
        recovery = {
            region: run.mean_recovery(region)
            for region in RECOVERY_TARGETS
        }

        comparisons.append({
            "config": config,
            "rmse": scored["rmse_log"],
            "recovery": recovery,
            "within_2x": scored["within_2x"],
            "within_5x": scored["within_5x"],
            "description": run.description,
        })

    if sort_by == "within_2x":
        comparisons.sort(key=lambda x: (-x["within_2x"], x["rmse"]))
    else:
        comparisons.sort(key=lambda x: x["rmse"])

    return comparisons


def format_comparison_table(runs: Dict[str, RunResult]) -> str:
    """Format a comparison table for multiple runs.

    Returns a fixed-width text table suitable for reports or agent output.
    """
    comparisons = compare_runs(runs)
    if not comparisons:
        return "(no runs to compare)"

    regions = SCORED_REGIONS
    # Header
    header_parts = [f"{'Config':<8}", f"{'RMSE':>6}"]
    for r in regions:
        header_parts.append(f"{r:>7}")
    header_parts.append(f"{'Grade':>8}")
    header = "  ".join(header_parts)

    lines = [header]
    lines.append("-" * len(header))

    for comp in comparisons:
        parts = [f"{comp['config']:<8}", f"{comp['rmse']:>6.3f}"]
        for r in regions:
            val = comp["recovery"].get(r, 0.0) * 100
            parts.append(f"{val:>6.1f}%")
        parts.append(f"{comp['within_2x']}/{comp.get('within_5x', len(regions))} 2×")
        lines.append("  ".join(parts))

    return "\n".join(lines)


def quick_score(run_dir: Path) -> str:
    """One-liner score string for a run directory.

    Example: "W201: RMSE=0.663 | AK-PWS 3.1% | 3/8 within 2×"
    """
    run_dir = Path(run_dir)
    run = load_run(run_dir)
    scored = score_run(run)

    pws_recovery = run.mean_recovery("AK-PWS") * 100

    return (
        f"{run.config_name}: RMSE={scored['rmse_log']:.3f}"
        f" | AK-PWS {pws_recovery:.1f}%"
        f" | {scored['within_2x']}/{scored['n_targets']} within 2×"
    )
