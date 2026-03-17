#!/usr/bin/env python3
"""Disease removal counterfactual experiment.

Runs two simulations with identical parameters:
  1. CONTROL: disease persists throughout (normal run)
  2. EXPERIMENT: disease removed at disease_end_year (all pathogen zeroed, infected cured)

This tests Augie's hypothesis: if populations recover after disease removal,
it proves ongoing pathogen pressure (not Allee effects or demographic debt)
is what keeps populations suppressed post-SSWD.

Usage:
    python experiments/disease_removal_experiment.py \
        --config experiments/calibration/W285_config.json \
        --disease-end-year 6 \
        --output-dir results/disease_removal/

    # disease_end_year is in simulation years (0-indexed).
    # Year 3 = epidemic start, so year 6 = 3 years post-epidemic.
"""

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

# Project root
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from experiments.calibration_runner import (
    build_full_network, load_sites, apply_param_overrides,
    compute_regional_recovery, score_against_targets,
)
from sswd_evoepi.config import default_config
from sswd_evoepi.model import run_spatial_simulation
from sswd_evoepi.metrics import RECOVERY_TARGETS


def run_experiment(config_path: str, disease_end_year: int, output_dir: str,
                   K: int = 5000, n_years: int = 13, disease_year: int = 3,
                   seed: int = 137):
    """Run control + disease-removal pair."""

    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)

    # Load config
    with open(config_path) as f:
        config_data = json.load(f)

    if 'param_overrides' in config_data:
        param_overrides = config_data['param_overrides']
    else:
        param_overrides = config_data

    # Extract spatial params
    n_connectivity = float(param_overrides.pop('spatial.n_connectivity', 0.3))
    alpha_self_open = float(param_overrides.pop('spatial.alpha_self_open', 0.05))
    alpha_self_fjord = float(param_overrides.pop('spatial.alpha_self_fjord', 0.50))
    r_total = float(param_overrides.pop('spatial.r_total', 0.02))
    phi_open = float(param_overrides.pop('spatial.phi_open', 0.8))
    phi_fjord = float(param_overrides.pop('spatial.phi_fjord', 0.03))
    D_P = param_overrides.pop('spatial.D_P', 15.0)
    D_L = param_overrides.pop('spatial.D_L', 400.0)

    print(f"Building 896-node network...")
    sites, network = build_full_network(
        K=K, seed=seed, D_L=D_L, D_P=D_P,
        n_connectivity=n_connectivity,
        alpha_self_open=alpha_self_open,
        alpha_self_fjord=alpha_self_fjord,
        r_total=r_total,
        phi_open=phi_open,
        phi_fjord=phi_fjord,
    )
    print(f"  {len(sites)} nodes ready")

    results = {}
    for label, end_year in [('control', None), ('disease_off', disease_end_year)]:
        print(f"\n{'='*60}")
        print(f"Running: {label.upper()}" +
              (f" (disease removed at sim year {end_year})" if end_year else " (disease persists)"))
        print(f"{'='*60}")

        config = default_config()
        config.simulation.sst_source = 'monthly'
        config.simulation.sst_data_dir = str(PROJECT_ROOT / 'data' / 'sst' / 'site_sst')
        config.simulation.sst_start_year = 2012
        config.simulation.seed = seed

        # Apply param overrides (make a copy so we don't mutate)
        overrides_copy = dict(param_overrides)
        if end_year is not None:
            overrides_copy['disease.disease_end_year'] = end_year
        apply_param_overrides(config, overrides_copy)

        t0 = time.time()
        result = run_spatial_simulation(
            network=network,
            n_years=n_years,
            disease_year=disease_year,
            seed=seed,
            config=config,
        )
        elapsed = time.time() - t0
        print(f"  Completed in {elapsed:.0f}s")

        # Score
        region_recovery, region_details = compute_regional_recovery(result, sites)
        scoring = score_against_targets(region_recovery)

        results[label] = {
            'scoring': scoring,
            'region_recovery': {k: float(v) for k, v in region_recovery.items()},
            'region_details': region_details,
            'wall_time': elapsed,
        }

        # Print recovery
        print(f"\n  Regional recovery ({label}):")
        for reg in sorted(region_recovery.keys()):
            target = RECOVERY_TARGETS.get(reg, None)
            val = region_recovery[reg] * 100
            marker = ""
            if target is not None:
                marker = f" (target: {target*100:.1f}%)"
            print(f"    {reg:8s}: {val:6.1f}%{marker}")
        print(f"  RMSLE: {scoring['rmsle']:.3f}")

    # Compare
    print(f"\n{'='*60}")
    print("COMPARISON: Disease ON vs Disease OFF")
    print(f"{'='*60}")
    print(f"{'Region':>10s} {'Control':>10s} {'Disease OFF':>12s} {'Δ':>10s}")
    print("-" * 45)
    for reg in sorted(set(results['control']['region_recovery'].keys()) |
                      set(results['disease_off']['region_recovery'].keys())):
        ctrl = results['control']['region_recovery'].get(reg, 0) * 100
        expr = results['disease_off']['region_recovery'].get(reg, 0) * 100
        delta = expr - ctrl
        print(f"{reg:>10s} {ctrl:>9.1f}% {expr:>11.1f}% {delta:>+9.1f}%")

    ctrl_rmsle = results['control']['scoring']['rmsle']
    expr_rmsle = results['disease_off']['scoring']['rmsle']
    print(f"\n  RMSLE: {ctrl_rmsle:.3f} → {expr_rmsle:.3f}")

    # Save
    out_path = output / 'disease_removal_results.json'
    with open(out_path, 'w') as f:
        json.dump({
            'config_path': str(config_path),
            'disease_end_year': disease_end_year,
            'seed': seed,
            'K': K,
            'n_years': n_years,
            'control': {k: v for k, v in results['control'].items() if k != 'region_details'},
            'disease_off': {k: v for k, v in results['disease_off'].items() if k != 'region_details'},
        }, f, indent=2, default=str)
    print(f"\nResults saved to {out_path}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Disease removal counterfactual experiment')
    parser.add_argument('--config', required=True, help='Config JSON path')
    parser.add_argument('--disease-end-year', type=int, default=6,
                        help='Simulation year to remove disease (default: 6 = 3yr post-epidemic)')
    parser.add_argument('--output-dir', default='results/disease_removal/',
                        help='Output directory')
    parser.add_argument('--K', type=int, default=5000, help='Carrying capacity')
    parser.add_argument('--n-years', type=int, default=13, help='Simulation years')
    parser.add_argument('--disease-year', type=int, default=3, help='Disease intro year')
    parser.add_argument('--seed', type=int, default=137, help='RNG seed')
    args = parser.parse_args()

    run_experiment(
        config_path=args.config,
        disease_end_year=args.disease_end_year,
        output_dir=args.output_dir,
        K=args.K,
        n_years=args.n_years,
        disease_year=args.disease_year,
        seed=args.seed,
    )
