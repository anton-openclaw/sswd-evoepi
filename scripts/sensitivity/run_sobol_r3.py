#!/usr/bin/env python3
"""SA Round 3 — Sobol N=512 execution script.

Runs 45,056 spatial simulations across 8 cores with checkpointing.
Resume-safe: checks for existing checkpoint and continues from there.

Usage:
    setsid python3 scripts/sensitivity/run_sobol_r3.py > results/sensitivity_r3/sobol_log.txt 2>&1 &
    disown
"""
import numpy as np
import json
import time
import os
import sys
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from scripts.sensitivity.param_spec import get_param_names
from scripts.sensitivity.spatial_runner import run_single_spatial, METRIC_NAMES

# ── Config ────────────────────────────────────────────────────────
N_CORES = 8
BASE_SEED = 31415
CHECKPOINT_INTERVAL = 500
RESULTS_DIR = "results/sensitivity_r3"

def main():
    names = get_param_names()
    
    # Load pre-generated samples
    data = np.load(f"{RESULTS_DIR}/sobol_samples.npz")
    X = data['X']
    total_runs = X.shape[0]
    
    print(f"SA Round 3 Sobol — {total_runs} runs, {len(names)} params, {len(METRIC_NAMES)} metrics")
    print(f"Cores: {N_CORES}, checkpoint every {CHECKPOINT_INTERVAL}")
    print(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    sys.stdout.flush()
    
    # Check for existing checkpoint
    checkpoint_path = f"{RESULTS_DIR}/sobol_checkpoint.npz"
    if os.path.exists(checkpoint_path):
        ckpt = np.load(checkpoint_path)
        Y = ckpt['Y']
        completed_mask = ckpt['completed_mask']
        n_done = int(np.sum(completed_mask))
        print(f"Resuming from checkpoint: {n_done}/{total_runs} done")
    else:
        Y = np.full((total_runs, len(METRIC_NAMES)), np.nan)
        completed_mask = np.zeros(total_runs, dtype=bool)
        n_done = 0
    
    # Build list of remaining runs
    remaining = [i for i in range(total_runs) if not completed_mask[i]]
    print(f"Runs remaining: {len(remaining)}")
    sys.stdout.flush()
    
    if len(remaining) == 0:
        print("All runs complete!")
        return
    
    # Build args for remaining runs
    args_list = [(i, X[i], names, BASE_SEED) for i in remaining]
    
    t_start = time.time()
    n_errors = 0
    batch_done = 0
    
    with Pool(N_CORES) as pool:
        for result in pool.imap_unordered(run_single_spatial, args_list):
            idx = result['run_index']
            
            # Store metrics
            for j, m in enumerate(METRIC_NAMES):
                Y[idx, j] = result['metrics'].get(m, np.nan)
            
            completed_mask[idx] = True
            if result.get('error'):
                n_errors += 1
            
            batch_done += 1
            total_done = n_done + batch_done
            
            # Progress + checkpoint
            if batch_done % CHECKPOINT_INTERVAL == 0:
                elapsed = time.time() - t_start
                rate = batch_done / elapsed
                eta_s = (len(remaining) - batch_done) / rate if rate > 0 else 0
                eta_h = eta_s / 3600
                
                # Save checkpoint
                np.savez(checkpoint_path, Y=Y, completed_mask=completed_mask)
                
                print(f"[{time.strftime('%H:%M:%S')}] {total_done}/{total_runs} "
                      f"({total_done/total_runs*100:.1f}%) | "
                      f"{rate:.2f} runs/s | "
                      f"ETA: {eta_h:.1f}h | "
                      f"errors: {n_errors}")
                sys.stdout.flush()
            
            # Also print every 2000 for coarser tracking
            elif batch_done % 2000 == 0:
                elapsed = time.time() - t_start
                rate = batch_done / elapsed
                eta_h = (len(remaining) - batch_done) / rate / 3600 if rate > 0 else 0
                print(f"[{time.strftime('%H:%M:%S')}] {total_done}/{total_runs} "
                      f"({total_done/total_runs*100:.1f}%) | "
                      f"{rate:.2f} runs/s | "
                      f"ETA: {eta_h:.1f}h | "
                      f"errors: {n_errors}")
                sys.stdout.flush()
    
    # Final save
    elapsed = time.time() - t_start
    np.savez(checkpoint_path, Y=Y, completed_mask=completed_mask)
    
    # Also save final Y matrix separately
    np.savez(f"{RESULTS_DIR}/sobol_Y_final.npz", Y=Y, X=X)
    
    total_done = n_done + batch_done
    print(f"\n{'='*60}")
    print(f"COMPLETE: {total_done}/{total_runs} runs")
    print(f"Errors: {n_errors} ({n_errors/total_runs*100:.1f}%)")
    print(f"Wall time: {elapsed/3600:.1f}h ({elapsed:.0f}s)")
    print(f"Mean rate: {batch_done/elapsed:.2f} runs/s")
    print(f"Finished: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Write completion marker
    with open(f"{RESULTS_DIR}/sobol_complete.json", 'w') as f:
        json.dump({
            'total_runs': total_runs,
            'completed': int(total_done),
            'errors': n_errors,
            'wall_time_hours': elapsed / 3600,
            'finished_at': time.strftime('%Y-%m-%d %H:%M:%S'),
        }, f, indent=2)


if __name__ == "__main__":
    main()
