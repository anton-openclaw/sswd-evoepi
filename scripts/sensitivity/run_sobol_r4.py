#!/usr/bin/env python3
"""SA Round 4 — Sobol N=512 execution script.

47 params, 23 metrics, 11-node stepping-stone network.
49,152 runs on 48 Xeon cores with checkpointing.

Usage:
    setsid nohup python3 -u scripts/sensitivity/run_sobol_r4.py \
        > results/sensitivity_r4/sobol_log.txt 2>&1 &
    disown
"""
import numpy as np
import json
import time
import os
import sys
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from SALib.sample import saltelli
from scripts.sensitivity.param_spec import get_param_names, get_salib_problem
from scripts.sensitivity.spatial_runner import run_single_spatial, METRIC_NAMES

# ── Config ────────────────────────────────────────────────────────
N_SOBOL = 512          # Saltelli N → N×(2D+2) = 512×96 = 49,152 runs
N_CORES = 48
BASE_SEED = 27182      # Different from Morris (31415) and R3 (31415)
CHECKPOINT_INTERVAL = 500
RESULTS_DIR = "results/sensitivity_r4"
PROGRESS_INTERVAL = 20  # Log every N completions


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    
    names = get_param_names()
    problem = get_salib_problem()
    n_params = len(names)
    n_metrics = len(METRIC_NAMES)
    
    # ── Generate or load samples ──────────────────────────────────
    samples_path = f"{RESULTS_DIR}/sobol_samples.npz"
    if os.path.exists(samples_path):
        data = np.load(samples_path)
        X = data['X']
        print(f"Loaded existing samples: {X.shape}")
    else:
        X = saltelli.sample(problem, N_SOBOL, calc_second_order=False)
        np.savez(samples_path, X=X)
        print(f"Generated Saltelli samples: {X.shape}")
    
    total_runs = X.shape[0]
    expected = N_SOBOL * (2 * n_params + 2)
    
    print(f"SA Round 4 Sobol — {total_runs} runs, {n_params} params, {n_metrics} metrics")
    print(f"Network: 11-node stepping-stone chain, K=5000/node (55K total)")
    print(f"Cores: {N_CORES}, checkpoint every {CHECKPOINT_INTERVAL}")
    print(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    if total_runs != expected:
        print(f"WARNING: expected {expected} runs, got {total_runs}")
    sys.stdout.flush()
    
    # ── Check for existing checkpoint ─────────────────────────────
    checkpoint_path = f"{RESULTS_DIR}/sobol_checkpoint.npz"
    if os.path.exists(checkpoint_path):
        ckpt = np.load(checkpoint_path)
        Y = ckpt['Y']
        completed_mask = ckpt['completed_mask']
        n_done = int(np.sum(completed_mask))
        print(f"Resuming from checkpoint: {n_done}/{total_runs} done")
    else:
        Y = np.full((total_runs, n_metrics), np.nan)
        completed_mask = np.zeros(total_runs, dtype=bool)
        n_done = 0
    
    # ── Build remaining work ──────────────────────────────────────
    remaining = [i for i in range(total_runs) if not completed_mask[i]]
    print(f"Runs remaining: {len(remaining)}")
    sys.stdout.flush()
    
    if len(remaining) == 0:
        print("All runs already complete!")
        _run_analysis(X, Y, names, problem)
        return
    
    args_list = [(i, X[i], names, BASE_SEED) for i in remaining]
    
    t_start = time.time()
    n_errors = 0
    batch_done = 0
    
    with Pool(N_CORES) as pool:
        for result in pool.imap_unordered(run_single_spatial, args_list):
            idx = result['run_index']
            
            for j, m in enumerate(METRIC_NAMES):
                Y[idx, j] = result['metrics'].get(m, np.nan)
            
            completed_mask[idx] = True
            if result.get('error'):
                n_errors += 1
            
            batch_done += 1
            total_done = n_done + batch_done
            
            # Progress logging
            if batch_done % PROGRESS_INTERVAL == 0:
                elapsed = time.time() - t_start
                rate = batch_done / elapsed if elapsed > 0 else 0
                eta_s = (len(remaining) - batch_done) / rate if rate > 0 else 0
                print(f"  {total_done}/{total_runs} done, "
                      f"{elapsed:.0f}s elapsed, {rate:.1f} runs/s, "
                      f"ETA {eta_s:.0f}s, errors={n_errors}")
                sys.stdout.flush()
            
            # Checkpoint
            if batch_done % CHECKPOINT_INTERVAL == 0:
                np.savez(checkpoint_path, Y=Y, completed_mask=completed_mask)
                print(f"  Checkpoint saved ({total_done}/{total_runs})")
                sys.stdout.flush()
    
    # ── Final save ────────────────────────────────────────────────
    elapsed = time.time() - t_start
    np.savez(checkpoint_path, Y=Y, completed_mask=completed_mask)
    np.savez(f"{RESULTS_DIR}/sobol_Y_final.npz", Y=Y, X=X)
    
    total_done = n_done + batch_done
    print(f"\nAll runs complete: {elapsed:.0f}s ({elapsed/3600:.1f}h), errors={n_errors}")
    sys.stdout.flush()
    
    # ── Run Sobol analysis ────────────────────────────────────────
    _run_analysis(X, Y, names, problem)
    
    # ── Write completion marker ───────────────────────────────────
    with open(f"{RESULTS_DIR}/sobol_complete.json", 'w') as f:
        json.dump({
            'total_runs': total_runs,
            'completed': int(total_done),
            'errors': n_errors,
            'wall_time_hours': elapsed / 3600,
            'finished_at': time.strftime('%Y-%m-%d %H:%M:%S'),
        }, f, indent=2)
    
    print("Done.")


def _run_analysis(X, Y, names, problem):
    """Run SALib Sobol analysis and save indices."""
    from SALib.analyze import sobol as sobol_analyze
    
    n_metrics = len(METRIC_NAMES)
    results = {}
    
    print(f"\n=== Sobol Analysis ===")
    
    for j, metric in enumerate(METRIC_NAMES):
        y_col = Y[:, j]
        
        # Skip if too many NaN
        valid_frac = np.sum(~np.isnan(y_col)) / len(y_col)
        if valid_frac < 0.8:
            print(f"  {metric}: skipped ({valid_frac:.0%} valid)")
            continue
        
        # Replace NaN with column mean for SALib
        y_clean = y_col.copy()
        y_clean[np.isnan(y_clean)] = np.nanmean(y_clean)
        
        # Check for zero variance
        if np.std(y_clean) < 1e-12:
            print(f"  {metric}: skipped (zero variance)")
            continue
        
        try:
            Si = sobol_analyze.analyze(
                problem, y_clean, calc_second_order=False,
                print_to_console=False,
            )
            results[metric] = {
                'S1': Si['S1'].tolist(),
                'ST': Si['ST'].tolist(),
                'S1_conf': Si['S1_conf'].tolist(),
                'ST_conf': Si['ST_conf'].tolist(),
            }
        except Exception as e:
            print(f"  {metric}: analysis failed — {e}")
    
    # Save indices
    out_path = f"{RESULTS_DIR}/sobol_indices.json"
    with open(out_path, 'w') as f:
        json.dump({
            'param_names': names,
            'metric_names': list(METRIC_NAMES),
            'indices': results,
        }, f, indent=2)
    print(f"\nSaved: {out_path}")
    
    # Print top-10 by mean ST
    all_ST = np.zeros(len(names))
    n_valid = 0
    for metric, si in results.items():
        all_ST += np.array(si['ST'])
        n_valid += 1
    
    if n_valid > 0:
        all_ST /= n_valid
        ranking = np.argsort(-all_ST)
        
        print(f"\n=== Top Parameters (mean ST across {n_valid} metrics) ===\n")
        print(f"{'Rank':>4} {'Parameter':<45} {'Mean ST':>8}")
        print("-" * 60)
        for r, idx in enumerate(ranking[:20]):
            print(f"{r+1:4d} {names[idx]:<45} {all_ST[idx]:.4f}")
    
    sys.stdout.flush()


if __name__ == "__main__":
    main()
