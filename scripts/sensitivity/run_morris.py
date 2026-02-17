#!/usr/bin/env python3
"""
Phase 1: Morris Screening
 
Run Morris method (Elementary Effects) to identify which of the 25 
parameters most influence model outputs. Cheap screening to reduce 
the parameter set before expensive Sobol analysis.

Usage: python3 scripts/sensitivity/run_morris.py [--trajectories 15] [--cores 8]
"""

import argparse
import json
import time
import sys
import os
import numpy as np
from multiprocessing import Pool

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from SALib.sample import morris as morris_sample
from SALib.analyze import morris as morris_analyze

from scripts.sensitivity.param_spec import get_param_names, get_salib_problem
from scripts.sensitivity.spatial_runner import run_single_spatial, METRIC_NAMES


def main():
    parser = argparse.ArgumentParser(description="Morris screening")
    parser.add_argument("--trajectories", type=int, default=20,
                        help="Number of Morris trajectories (default: 20)")
    parser.add_argument("--cores", type=int, default=8,
                        help="Number of parallel workers (default: 8)")
    parser.add_argument("--seed", type=int, default=12345,
                        help="Base RNG seed")
    parser.add_argument("--outdir", type=str, 
                        default="results/sensitivity",
                        help="Output directory")
    args = parser.parse_args()
    
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)
    
    param_names = get_param_names()
    problem = get_salib_problem(param_names)
    n_params = len(param_names)
    
    print(f"Morris screening: {n_params} parameters, {args.trajectories} trajectories")
    print(f"Total runs: {args.trajectories * (n_params + 1)}")
    
    # Generate Morris samples
    t0 = time.time()
    X = morris_sample.sample(problem, N=args.trajectories, seed=args.seed)
    n_runs = len(X)
    print(f"Generated {n_runs} sample points in {time.time()-t0:.1f}s")
    
    # Run all simulations
    t0 = time.time()
    tasks = [
        (i, X[i], param_names, args.seed * 1000)
        for i in range(n_runs)
    ]
    
    print(f"Running {n_runs} spatial simulations on {args.cores} cores (~{n_runs*6.5/args.cores/60:.0f} min)...")
    with Pool(processes=args.cores) as pool:
        results_list = pool.map(run_single_spatial, tasks)
    
    elapsed = time.time() - t0
    print(f"Completed {n_runs} runs in {elapsed:.1f}s ({elapsed/n_runs:.2f}s/run)")
    
    # Count errors
    errors = [r for r in results_list if r["error"] is not None]
    if errors:
        print(f"WARNING: {len(errors)}/{n_runs} runs failed:")
        for e in errors[:5]:
            print(f"  Run {e['run_index']}: {e['error']}")
    
    # Build output matrix Y: (n_runs, n_metrics)
    Y = np.zeros((n_runs, len(METRIC_NAMES)))
    for r in results_list:
        idx = r["run_index"]
        for j, m in enumerate(METRIC_NAMES):
            Y[idx, j] = r["metrics"][m]
    
    # Analyze each metric
    print("\n" + "="*70)
    print("MORRIS SCREENING RESULTS")
    print("="*70)
    
    morris_results = {}
    
    for j, metric in enumerate(METRIC_NAMES):
        y = Y[:, j]
        
        # Skip if all NaN or constant
        valid = ~np.isnan(y)
        if np.sum(valid) < n_runs * 0.5:
            print(f"\n{metric}: >50% NaN — skipping")
            continue
        
        # Replace NaN with median for analysis
        if np.any(~valid):
            y_clean = y.copy()
            y_clean[~valid] = np.nanmedian(y)
        else:
            y_clean = y
        
        try:
            Si = morris_analyze.analyze(
                problem, X, y_clean,
                num_resamples=100,
                conf_level=0.95,
                seed=args.seed,
            )
            
            morris_results[metric] = {
                "mu_star": Si["mu_star"].tolist(),
                "mu_star_conf": Si["mu_star_conf"].tolist(),
                "sigma": Si["sigma"].tolist(),
                "names": list(Si["names"]),
            }
            
            # Print top 10
            ranked = np.argsort(Si["mu_star"])[::-1]
            print(f"\n{'─'*60}")
            print(f"Metric: {metric}")
            print(f"{'─'*60}")
            print(f"  {'Rank':<5} {'Parameter':<40} {'μ*':>8} {'σ':>8}")
            for rank, idx in enumerate(ranked[:10], 1):
                print(f"  {rank:<5} {param_names[idx]:<40} {Si['mu_star'][idx]:>8.3f} {Si['sigma'][idx]:>8.3f}")
        
        except Exception as e:
            print(f"\n{metric}: Analysis failed — {e}")
    
    # Determine screened parameters
    # Criterion: parameter is "important" if it ranks in top 10 for ANY metric
    # by having μ* > 10% of max μ* for that metric
    important_params = set()
    
    for metric, data in morris_results.items():
        mu_star = np.array(data["mu_star"])
        max_mu = np.max(mu_star)
        if max_mu > 0:
            threshold = max_mu * 0.05  # 5% of max effect
            for i, ms in enumerate(mu_star):
                if ms >= threshold:
                    important_params.add(param_names[i])
    
    important_list = sorted(important_params)
    
    print(f"\n{'='*70}")
    print(f"SCREENING SUMMARY")
    print(f"{'='*70}")
    print(f"Important parameters ({len(important_list)}/{n_params}):")
    for p in important_list:
        print(f"  ✓ {p}")
    
    eliminated = set(param_names) - important_params
    if eliminated:
        print(f"\nEliminated parameters ({len(eliminated)}):")
        for p in sorted(eliminated):
            print(f"  ✗ {p}")
    
    # Save results
    output = {
        "method": "Morris Elementary Effects",
        "n_params": n_params,
        "n_trajectories": args.trajectories,
        "n_runs": n_runs,
        "n_errors": len(errors),
        "elapsed_s": float(elapsed),
        "param_names": list(param_names),
        "metric_names": list(METRIC_NAMES),
        "morris_results": morris_results,
        "screened_params": list(important_list),
        "eliminated_params": list(sorted(eliminated)),
        "screening_threshold": "mu_star >= 5% of max mu_star for any metric",
        "seed": args.seed,
    }
    
    outpath = os.path.join(outdir, "morris_screening.json")
    with open(outpath, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nResults saved: {outpath}")
    
    # Save screened param list for downstream Sobol
    screened_path = os.path.join(outdir, "screened_params.json")
    with open(screened_path, "w") as f:
        json.dump({"screened_params": important_list}, f, indent=2)
    print(f"Screened params saved: {screened_path}")
    
    return output


if __name__ == "__main__":
    main()
