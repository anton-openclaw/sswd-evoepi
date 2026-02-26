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

TARGETS = {
    'AK-PWS': 0.50, 'AK-FN': 0.50, 'AK-FS': 0.20, 'BC-N': 0.20,
    'SS-S': 0.05, 'JDF': 0.02, 'OR': 0.0025, 'CA-N': 0.001,
}

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
    
    # Load individual seed results
    seed_files = sorted(round_dir.glob('result_seed*.json'))
    if not seed_files:
        return f"## {round_name}: NO RESULTS YET\n"
    
    results = []
    for sf in seed_files:
        with open(sf) as f:
            results.append(json.load(f))
    
    # Compute means across seeds
    lines = []
    lines.append(f"## {round_name}")
    
    if params:
        lines.append(f"Parameters changed: {json.dumps(params)}")
    else:
        lines.append("Parameters: ALL DEFAULTS")
    
    lines.append(f"Seeds: {[r['seed'] for r in results]}")
    
    mean_rmse = np.mean([r['scoring']['rmse_log'] for r in results])
    mean_crash = np.mean([r['overall']['pop_crash_pct'] for r in results])
    lines.append(f"Mean RMSE(log10): {mean_rmse:.3f}")
    lines.append(f"Mean crash: {mean_crash:.1f}%")
    lines.append(f"Wall time per seed: {np.mean([r['wall_time_seconds'] for r in results]):.0f}s")
    
    # Per-region table
    lines.append("")
    lines.append(f"{'Region':<10s} {'Target':>8s} {'Actual':>10s} {'LogErr':>8s} {'Grade'}")
    lines.append("-" * 45)
    
    for region in TARGETS:
        target = TARGETS[region]
        actuals = [r['scoring']['per_region'][region]['actual'] for r in results]
        mean_actual = np.mean(actuals)
        std_actual = np.std(actuals) if len(actuals) > 1 else 0
        
        log_err = np.log10(max(mean_actual, 1e-6)) - np.log10(max(target, 1e-6))
        
        if abs(log_err) < 0.301:
            grade = "✓ 2×"
        elif abs(log_err) < 0.699:
            grade = "~ 5×"
        else:
            grade = f"✗ {10**abs(log_err):.0f}×"
        
        actual_str = f"{mean_actual*100:.2f}%"
        if std_actual > 0:
            actual_str += f"±{std_actual*100:.2f}"
        
        lines.append(f"{region:<10s} {target*100:>7.2f}% {actual_str:>10s} {log_err:>+7.2f}  {grade}")
    
    # Unconstrained regions
    all_regions = set()
    for r in results:
        all_regions.update(r['region_recovery'].keys())
    unconstrained = sorted(all_regions - set(TARGETS.keys()))
    
    if unconstrained:
        lines.append("")
        lines.append("Unconstrained regions:")
        for region in unconstrained:
            vals = [r['region_recovery'].get(region, 0) for r in results]
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
