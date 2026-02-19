#!/usr/bin/env python3
"""Phase 5: Daily population trajectory plot — should be smooth, not staircase."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from sswd_evoepi.config import default_config
from sswd_evoepi.model import run_coupled_simulation


def main():
    config = default_config()
    print("Running single-node 10yr sim (no disease, record_daily=True)...")

    result = run_coupled_simulation(
        n_individuals=5000,
        n_years=10,
        disease_year=None,  # No disease
        seed=42,
        config=config,
        record_daily=True,
    )

    daily_pop = result.daily_pop
    total_days = len(daily_pop)
    days = np.arange(total_days)

    # Plot
    fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True)

    # Panel 1: Full daily trajectory
    ax1 = axes[0]
    ax1.plot(days / 365.0, daily_pop, color='#2196F3', linewidth=0.5, alpha=0.8)
    ax1.set_ylabel('Population', fontsize=12)
    ax1.set_title('Daily Population Trajectory — Continuous Mortality + Growth\n'
                   '(Single node, K=5000, no disease, 10 years)',
                   fontsize=13, fontweight='bold')

    # Mark year boundaries
    for y in range(1, 11):
        ax1.axvline(y, color='gray', linestyle='--', alpha=0.3, linewidth=0.5)

    ax1.grid(True, alpha=0.3)

    # Panel 2: Zoom into years 2-4 to show smoothness clearly
    ax2 = axes[1]
    start_day = 2 * 365
    end_day = 4 * 365
    zoom_days = days[start_day:end_day]
    zoom_pop = daily_pop[start_day:end_day]
    ax2.plot(zoom_days / 365.0, zoom_pop, color='#4CAF50', linewidth=1.0)
    ax2.set_xlabel('Year', fontsize=12)
    ax2.set_ylabel('Population', fontsize=12)
    ax2.set_title('Zoom: Years 2–4 (smooth daily changes, no year-boundary drops)',
                   fontsize=11)

    # Mark year boundaries in zoom
    for y in range(2, 5):
        ax2.axvline(y, color='red', linestyle='--', alpha=0.5, linewidth=1,
                     label='Year boundary' if y == 2 else None)

    ax2.legend(loc='upper right')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    os.makedirs("results/continuous_mortality", exist_ok=True)
    outpath = "results/continuous_mortality/daily_pop_smooth.png"
    fig.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {outpath}")

    # Check for staircase pattern: compute day-to-day changes
    diffs = np.diff(daily_pop[:end_day])
    # In a staircase, you'd see long runs of zero change then big jumps
    zero_runs = np.sum(diffs == 0)
    total_diffs = len(diffs)
    print(f"\nSmoothness check (days 0-{end_day}):")
    print(f"  Day-to-day diffs: {total_diffs}")
    print(f"  Zero-change days: {zero_runs} ({zero_runs/total_diffs*100:.1f}%)")
    print(f"  Non-zero days: {total_diffs - zero_runs} ({(total_diffs-zero_runs)/total_diffs*100:.1f}%)")
    big_jumps = np.sum(np.abs(diffs) > 50)
    print(f"  Large jumps (|Δ|>50): {big_jumps}")

    # Check year boundaries specifically
    print(f"\nYear-boundary population changes:")
    for y in range(1, min(10, len(daily_pop) // 365)):
        boundary_day = y * 365
        if boundary_day < len(daily_pop) - 1:
            before = daily_pop[boundary_day - 1]
            after = daily_pop[boundary_day]
            change = after - before
            print(f"  Year {y}→{y+1}: {before} → {after} (Δ={change:+d})")


if __name__ == "__main__":
    main()
