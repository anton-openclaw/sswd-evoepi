#!/usr/bin/env python3
"""Run all W325 figure generation scripts."""
import subprocess
import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(SCRIPT_DIR, "..", "..", "reports", "w325", "figures")

scripts = [
    "fig_3_1_infection_heatmap.py",
    "fig_3_2_population_heatmap.py",
    "fig_3_3_recovery_vs_latitude.py",
    "fig_3_4_cascade_panels.py",
    "fig_3_5_recovery_map.py",
    "fig_4_1_resistance_evolution.py",
    "fig_4_2_three_trait_panel.py",
    "fig_4_3_va_trajectory.py",
    "fig_4_4_pathogen_evolution.py",
    "fig_5_1_recovery_vs_targets.py",
    "fig_5_2_scoring_breakdown.py",
    "fig_A_1_report_card.py",
]

failed = []
for script in scripts:
    path = os.path.join(SCRIPT_DIR, script)
    print(f"\n{'='*60}")
    print(f"Running {script}...")
    print(f"{'='*60}")
    result = subprocess.run([sys.executable, path], capture_output=True, text=True)
    if result.stdout:
        print(result.stdout)
    if result.stderr:
        print(result.stderr, file=sys.stderr)
    if result.returncode != 0:
        failed.append(script)
        print(f"  ✗ FAILED (exit code {result.returncode})")
    else:
        print(f"  ✓ OK")

print(f"\n{'='*60}")
if failed:
    print(f"FAILED: {len(failed)}/{len(scripts)}")
    for s in failed:
        print(f"  ✗ {s}")
else:
    print(f"ALL {len(scripts)} figures generated successfully!")
print(f"{'='*60}")
