#!/usr/bin/env python3
"""
Phase 2: Sobol Sensitivity Analysis

Run Sobol variance decomposition on the screened parameter subset.
Generates Saltelli samples, runs simulations in parallel with 
checkpointing, and computes first-order (S1) and total-order (ST) indices.

Usage: 
    python3 scripts/sensitivity/run_sobol.py [--N 512] [--cores 8] [--batch 0] [--n-batches 1]
    
For parallel batch execution:
    python3 scripts/sensitivity/run_sobol.py --batch 0 --n-batches 4  # runs batch 0 of 4
    python3 scripts/sensitivity/run_sobol.py --batch 1 --n-batches 4  # runs batch 1 of 4
"""

import argparse
import json
import time
import sys
import os
import numpy as np
from multiprocessing import Pool

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from SALib.sample import saltelli
from SALib.analyze import sobol as sobol_analyze

from scripts.sensitivity.param_spec import get_param_names, get_salib_problem, PARAM_SPEC
from scripts.sensitivity.spatial_runner import run_single_spatial, METRIC_NAMES


def load_screened_params(outdir):
    """Load screened parameter list from Morris phase."""
    path = os.path.join(outdir, "screened_params.json")
    if os.path.exists(path):
        with open(path) as f:
            data = json.load(f)
        return data["screened_params"]
    else:
        print("WARNING: No screened_params.json found — using all 25 parameters")
        return get_param_names()


def main():
    parser = argparse.ArgumentParser(description="Sobol sensitivity analysis")
    parser.add_argument("--N", type=int, default=512,
                        help="Base sample size (default: 512)")
    parser.add_argument("--cores", type=int, default=8,
                        help="Number of parallel workers (default: 8)")
    parser.add_argument("--seed", type=int, default=54321,
                        help="Base RNG seed")
    parser.add_argument("--batch", type=int, default=0,
                        help="Batch index (0-indexed)")
    parser.add_argument("--n-batches", type=int, default=1,
                        help="Total number of batches")
    parser.add_argument("--outdir", type=str,
                        default="results/sensitivity",
                        help="Output directory")
    parser.add_argument("--checkpoint-every", type=int, default=500,
                        help="Checkpoint interval (runs)")
    parser.add_argument("--params", type=str, default=None,
                        help="Comma-separated param names (override screened list)")
    args = parser.parse_args()
    
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)
    
    # Load parameter list
    if args.params:
        param_names = args.params.split(",")
    else:
        param_names = load_screened_params(outdir)
    
    n_params = len(param_names)
    problem = get_salib_problem(param_names)
    
    # Generate full Saltelli sample matrix
    # Total samples = N * (2k + 2)
    X = saltelli.sample(problem, N=args.N, seed=args.seed)
    n_total = len(X)
    
    print(f"Sobol analysis: {n_params} parameters, N={args.N}")
    print(f"Total Saltelli samples: {n_total}")
    print(f"Parameters: {param_names}")
    
    # Save the full sample matrix (needed for analysis)
    sample_path = os.path.join(outdir, "sobol_samples.npz")
    if not os.path.exists(sample_path) or args.batch == 0:
        np.savez_compressed(sample_path, X=X, param_names=param_names)
        print(f"Sample matrix saved: {sample_path} ({X.shape})")
    
    # Determine this batch's slice
    batch_size = n_total // args.n_batches
    start_idx = args.batch * batch_size
    end_idx = n_total if args.batch == args.n_batches - 1 else start_idx + batch_size
    batch_runs = end_idx - start_idx
    
    print(f"\nBatch {args.batch}/{args.n_batches}: runs {start_idx}-{end_idx-1} ({batch_runs} runs)")
    
    # Check for existing checkpoint
    ckpt_path = os.path.join(outdir, f"sobol_batch{args.batch}_checkpoint.npz")
    completed_from = start_idx
    Y_batch = np.full((batch_runs, len(METRIC_NAMES)), np.nan)
    runtimes = np.zeros(batch_runs)
    errors_list = []
    
    if os.path.exists(ckpt_path):
        ckpt = np.load(ckpt_path)
        Y_prev = ckpt["Y"]
        n_completed = int(ckpt["n_completed"])
        completed_from = start_idx + n_completed
        Y_batch[:n_completed] = Y_prev[:n_completed]
        runtimes[:n_completed] = ckpt.get("runtimes", np.zeros(n_completed))[:n_completed]
        print(f"Resumed from checkpoint: {n_completed} runs already done")
    
    # Run remaining
    remaining = end_idx - completed_from
    if remaining <= 0:
        print("All runs in this batch already complete!")
    else:
        t0 = time.time()
        
        # Build tasks for remaining runs
        tasks = []
        for i in range(completed_from, end_idx):
            local_i = i - start_idx
            tasks.append((i, X[i], param_names, args.seed * 1000))
        
        print(f"Running {remaining} simulations on {args.cores} cores...")
        
        # Process in chunks for checkpointing
        chunk_size = args.checkpoint_every
        tasks_done = 0
        
        with Pool(processes=args.cores) as pool:
            for chunk_start in range(0, len(tasks), chunk_size):
                chunk = tasks[chunk_start:chunk_start + chunk_size]
                results = pool.map(run_single_spatial, chunk)
                
                for r in results:
                    local_i = r["run_index"] - start_idx
                    for j, m in enumerate(METRIC_NAMES):
                        Y_batch[local_i, j] = r["metrics"][m]
                    runtimes[local_i] = r["runtime"]
                    if r["error"]:
                        errors_list.append({
                            "run_index": r["run_index"],
                            "error": r["error"],
                        })
                
                tasks_done += len(chunk)
                n_done = completed_from - start_idx + tasks_done
                
                # Checkpoint
                np.savez_compressed(ckpt_path,
                    Y=Y_batch,
                    n_completed=n_done,
                    runtimes=runtimes,
                )
                
                elapsed = time.time() - t0
                rate = tasks_done / elapsed if elapsed > 0 else 0
                eta = (remaining - tasks_done) / rate if rate > 0 else 0
                print(f"  Progress: {n_done}/{batch_runs} | "
                      f"{rate:.1f} runs/s | ETA: {eta/60:.1f} min | "
                      f"Errors: {len(errors_list)}")
        
        total_elapsed = time.time() - t0
        print(f"\nBatch {args.batch} complete: {remaining} runs in {total_elapsed:.1f}s "
              f"({total_elapsed/remaining:.2f}s/run)")
    
    # Save batch results
    batch_path = os.path.join(outdir, f"sobol_batch{args.batch}_results.npz")
    np.savez_compressed(batch_path,
        Y=Y_batch,
        start_idx=start_idx,
        end_idx=end_idx,
        runtimes=runtimes,
        metric_names=METRIC_NAMES,
        param_names=param_names,
    )
    print(f"Batch results saved: {batch_path}")
    
    if errors_list:
        err_path = os.path.join(outdir, f"sobol_batch{args.batch}_errors.json")
        with open(err_path, "w") as f:
            json.dump(errors_list, f, indent=2)
        print(f"Errors logged: {err_path} ({len(errors_list)} errors)")
    
    # If single batch (n_batches=1) or this is the analysis job, do full analysis
    if args.n_batches == 1:
        analyze_sobol(outdir, X, Y_batch, problem, param_names)
    
    return Y_batch


def analyze_sobol(outdir, X, Y, problem, param_names):
    """Run Sobol analysis on complete Y matrix."""
    print(f"\n{'='*70}")
    print("SOBOL SENSITIVITY ANALYSIS")
    print(f"{'='*70}")
    
    sobol_results = {}
    
    for j, metric in enumerate(METRIC_NAMES):
        y = Y[:, j]
        
        # Check for enough valid data
        valid = ~np.isnan(y)
        valid_pct = np.sum(valid) / len(y) * 100
        
        if valid_pct < 50:
            print(f"\n{metric}: Only {valid_pct:.0f}% valid — skipping")
            continue
        
        # Replace NaN with median
        if np.any(~valid):
            y_clean = y.copy()
            y_clean[~valid] = np.nanmedian(y)
        else:
            y_clean = y
        
        # Check variance
        if np.std(y_clean) < 1e-10:
            print(f"\n{metric}: No variance — skipping")
            continue
        
        try:
            Si = sobol_analyze.analyze(
                problem, y_clean,
                calc_second_order=False,  # saves compute; add later if needed
                num_resamples=1000,
                conf_level=0.95,
                seed=54321,
            )
            
            sobol_results[metric] = {
                "S1": Si["S1"].tolist(),
                "S1_conf": Si["S1_conf"].tolist(),
                "ST": Si["ST"].tolist(),
                "ST_conf": Si["ST_conf"].tolist(),
            }
            
            # Print ranked by ST
            ranked = np.argsort(Si["ST"])[::-1]
            print(f"\n{'─'*70}")
            print(f"Metric: {metric}")
            print(f"{'─'*70}")
            print(f"  {'Rank':<5} {'Parameter':<40} {'S1':>8} {'ST':>8} {'Inter':>8}")
            for rank, idx in enumerate(ranked[:10], 1):
                s1 = Si["S1"][idx]
                st = Si["ST"][idx]
                inter = st - s1  # interaction contribution
                print(f"  {rank:<5} {param_names[idx]:<40} {s1:>8.3f} {st:>8.3f} {inter:>8.3f}")
        
        except Exception as e:
            print(f"\n{metric}: Analysis failed — {e}")
            import traceback
            traceback.print_exc()
    
    # Save
    output = {
        "method": "Sobol (Saltelli)",
        "n_params": len(param_names),
        "n_samples": len(Y),
        "param_names": param_names,
        "metric_names": METRIC_NAMES,
        "sobol_results": sobol_results,
    }
    
    outpath = os.path.join(outdir, "sobol_results.json")
    with open(outpath, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nSobol results saved: {outpath}")
    
    return sobol_results


if __name__ == "__main__":
    main()
