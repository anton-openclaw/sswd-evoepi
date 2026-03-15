#!/usr/bin/env python3
"""Format calibration results for the tuning agent.

Reads completed round results and produces a concise summary
suitable for passing to the tuning agent as context.

Usage:
    python3 format_results.py results/calibration/round_00/ [round_01/ round_02/ ...]
"""

import json
import sys
import numpy as np
from pathlib import Path

# Add project root to path (3 levels up: experiments/calibration/format_results.py)
PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from sswd_evoepi.metrics import RECOVERY_TARGETS, score_regions
from sswd_evoepi.results import load_run


def format_round(round_dir: Path) -> str:
    """Format a single round's results."""
    round_name = round_dir.name
    
    # Load params
    params_file = round_dir / 'params.json'
    if params_file.exists():
        with open(params_file) as f:
            params = json.load(f)
    else:
        params = {}
    
    # Load seed results via consolidated module
    run = load_run(round_dir)
    if not run.seeds:
        return f"## {round_name}: NO RESULTS YET\n"
    
    # Filter out early-stopped / empty results
    valid_seeds = [s for s in run.seeds if s.scoring.get('per_region')]
    
    # Compute means across seeds
    lines = []
    lines.append(f"## {round_name}")
    
    # Load config.json as fallback for params
    config_file = round_dir / 'config.json'
    if not params and config_file.exists():
        with open(config_file) as f:
            params = json.load(f)
    
    if params:
        lines.append(f"Parameters changed: {json.dumps(params)}")
    else:
        lines.append("Parameters: ALL DEFAULTS")
    
    lines.append(f"Seeds: {[s.seed for s in run.seeds]}")
    early_stopped = len(run.seeds) - len(valid_seeds)
    if early_stopped:
        lines.append(f"Early-stopped seeds: {early_stopped}/{len(run.seeds)}")
    
    finite_rmse = [s.rmse_log for s in run.seeds if s.rmse_log != float('inf')]
    mean_rmse_val = np.mean(finite_rmse) if finite_rmse else float('inf')
    mean_crash = np.mean([s.raw['overall']['pop_crash_pct'] for s in run.seeds])
    lines.append(f"Mean RMSE(log): {mean_rmse_val:.3f}" if mean_rmse_val != float('inf') else "Mean RMSE(log): INF")
    lines.append(f"Mean crash: {mean_crash:.1f}%")
    lines.append(f"Wall time per seed: {np.mean([s.wall_time for s in run.seeds]):.0f}s")
    
    if not valid_seeds:
        lines.append("\nNo valid regional results (all early-stopped)")
        return "\n".join(lines)
    
    # Compute mean recovery across valid seeds and score with consolidated module
    mean_recovery = {}
    std_recovery = {}
    for region in RECOVERY_TARGETS:
        actuals = [s.recovery(region) for s in valid_seeds]
        mean_recovery[region] = np.mean(actuals)
        std_recovery[region] = np.std(actuals) if len(actuals) > 1 else 0
    
    scored = score_regions(mean_recovery)
    
    # Per-region table
    lines.append("")
    lines.append(f"{'Region':<10s} {'Target':>8s} {'Actual':>10s} {'LogErr':>8s} {'Grade'}")
    lines.append("-" * 45)
    
    for region in RECOVERY_TARGETS:
        target = RECOVERY_TARGETS[region]
        mean_actual = mean_recovery[region]
        std_actual = std_recovery[region]
        info = scored[region]
        
        actual_str = f"{mean_actual*100:.2f}%"
        if std_actual > 0:
            actual_str += f"±{std_actual*100:.2f}"
        
        lines.append(f"{region:<10s} {target*100:>7.2f}% {actual_str:>10s} {info['log_error']:>+7.2f}  {info['grade_str']}")
    
    # Unconstrained regions
    all_regions = set()
    for s in run.seeds:
        all_regions.update(s.region_recovery.keys())
    unconstrained = sorted(all_regions - set(RECOVERY_TARGETS.keys()))
    
    if unconstrained:
        lines.append("")
        lines.append("Unconstrained regions:")
        for region in unconstrained:
            vals = [s.recovery(region) for s in run.seeds]
            lines.append(f"  {region:<10s}: {np.mean(vals)*100:.2f}%")
    
    return "\n".join(lines)


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 format_results.py round_dir1 [round_dir2 ...]")
        sys.exit(1)
    
    all_text = []
    for arg in sys.argv[1:]:
        round_dir = Path(arg)
        if round_dir.exists():
            all_text.append(format_round(round_dir))
    
    print("\n\n".join(all_text))


if __name__ == '__main__':
    main()
