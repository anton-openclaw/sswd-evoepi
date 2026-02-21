#!/usr/bin/env python3
"""SA Round 4 — Sobol N=512 execution script on 11-node network.

Key changes from R3 Sobol:
- 11-node stepping-stone chain (was 3-node)
- 47 params (was 43), 23 metrics (was 21)
- 48 cores on Xeon (was 12)
- N=512 → 49,152 total runs (N × (2D+2) = 512 × 96)
- Three-trait genetic architecture (resistance/tolerance/recovery)
- Pathogen evolution parameters included

Usage:
    setsid nohup python3 -u scripts/sensitivity/run_sobol_r4.py \
        > results/sensitivity_r4/sobol_log.txt 2>&1 &
"""
import numpy as np
import json
import time
import os
import sys
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from scripts.sensitivity.param_spec import get_param_names, get_salib_problem
from scripts.sensitivity.spatial_runner import run_single_spatial, METRIC_NAMES
from SALib.sample import saltelli

# ── Config ────────────────────────────────────────────────────────
N_SOBOL = 512           # Saltelli base sample size
N_CORES = 48            # Xeon W-3365: 64 cores / 128 threads
BASE_SEED = 88888       # Distinct from R3 (31415) and R4 Morris (77777)
CHECKPOINT_INTERVAL = 500
PROGRESS_INTERVAL = 20
RESULTS_DIR = "results/sensitivity_r4"


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)

    names = get_param_names()
    problem = get_salib_problem(names)
    n_params = len(names)

    # Generate Saltelli samples
    # calc_second_order=True → N × (2D+2) samples
    # Set global seed (saltelli.sample doesn't accept seed kwarg in this SALib version)
    np.random.seed(BASE_SEED)
    X = saltelli.sample(problem, N_SOBOL, calc_second_order=True)
    total_runs = X.shape[0]
    expected = N_SOBOL * (2 * n_params + 2)

    print(f"SA Round 4 Sobol — {total_runs} runs (expected {expected})")
    print(f"  {n_params} params, {len(METRIC_NAMES)} metrics")
    print(f"  Network: 11-node stepping-stone chain, K=5000/node (55K total)")
    print(f"  Cores: {N_CORES}, checkpoint every {CHECKPOINT_INTERVAL}")
    print(f"  Progress log every {PROGRESS_INTERVAL} runs")
    print(f"  Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    sys.stdout.flush()

    # Save samples for reproducibility and later analysis
    samples_path = f"{RESULTS_DIR}/sobol_samples.npz"
    if not os.path.exists(samples_path):
        np.savez(samples_path, X=X, param_names=names)
        print(f"  Saved samples: {samples_path} ({X.shape})")

    # Check for existing checkpoint
    checkpoint_path = f"{RESULTS_DIR}/sobol_checkpoint.npz"
    if os.path.exists(checkpoint_path):
        ckpt = np.load(checkpoint_path)
        Y = ckpt['Y']
        completed_mask = ckpt['completed_mask']
        n_done = int(np.sum(completed_mask))
        print(f"  Resuming from checkpoint: {n_done}/{total_runs} done")
    else:
        Y = np.full((total_runs, len(METRIC_NAMES)), np.nan)
        completed_mask = np.zeros(total_runs, dtype=bool)
        n_done = 0

    # Build list of remaining runs
    remaining = [i for i in range(total_runs) if not completed_mask[i]]
    print(f"  Runs remaining: {len(remaining)}")
    sys.stdout.flush()

    if len(remaining) == 0:
        print("All runs complete!")
        _save_final(Y, X, names, total_runs, n_done, 0, 0.0)
        return

    # Build args for remaining runs
    args_list = [(i, X[i], names, BASE_SEED) for i in remaining]

    t_start = time.time()
    last_ckpt_time = t_start
    n_errors = 0
    batch_done = 0

    with Pool(N_CORES) as pool:
        for result in pool.imap_unordered(run_single_spatial, args_list):
            idx = result['run_index']

            # Store metrics
            if result.get('error'):
                n_errors += 1
            else:
                for j, m in enumerate(METRIC_NAMES):
                    Y[idx, j] = result['metrics'].get(m, np.nan)

            completed_mask[idx] = True
            batch_done += 1
            total_done = n_done + batch_done

            # Progress logging
            if batch_done % PROGRESS_INTERVAL == 0:
                elapsed = time.time() - t_start
                rate = batch_done / elapsed
                eta_s = (len(remaining) - batch_done) / rate if rate > 0 else 0
                eta_h = eta_s / 3600
                pct = total_done / total_runs * 100
                print(f"[{time.strftime('%H:%M:%S')}] {total_done}/{total_runs} "
                      f"({pct:.1f}%) | "
                      f"{rate:.2f} r/s | "
                      f"ETA: {eta_h:.1f}h | "
                      f"errors: {n_errors}")
                sys.stdout.flush()

            # Checkpoint
            if batch_done % CHECKPOINT_INTERVAL == 0:
                np.savez(checkpoint_path, Y=Y, completed_mask=completed_mask)
                ckpt_elapsed = time.time() - last_ckpt_time
                last_ckpt_time = time.time()
                ckpt_rate = CHECKPOINT_INTERVAL / ckpt_elapsed if ckpt_elapsed > 0 else 0
                print(f"  ── Checkpoint saved ({total_done}/{total_runs}) | "
                      f"last {CHECKPOINT_INTERVAL}: {ckpt_rate:.2f} r/s ──")
                sys.stdout.flush()

    # Final save
    elapsed = time.time() - t_start
    np.savez(checkpoint_path, Y=Y, completed_mask=completed_mask)
    _save_final(Y, X, names, total_runs, n_done + batch_done, n_errors, elapsed)


def _save_final(Y, X, names, total_runs, total_done, n_errors, elapsed):
    """Save final results and completion marker."""
    # Save Y matrix + samples for analysis
    np.savez(f"{RESULTS_DIR}/sobol_Y_final.npz", Y=Y, X=X)

    # Save completion marker
    completion = {
        'round': 4,
        'method': 'sobol',
        'N': N_SOBOL,
        'n_params': len(names),
        'n_metrics': len(METRIC_NAMES),
        'total_runs': total_runs,
        'completed': total_done,
        'errors': n_errors,
        'wall_time_hours': elapsed / 3600,
        'cores': N_CORES,
        'network': '11-node stepping-stone chain',
        'base_seed': BASE_SEED,
        'finished_at': time.strftime('%Y-%m-%d %H:%M:%S'),
        'param_names': names,
        'metric_names': METRIC_NAMES,
    }

    with open(f"{RESULTS_DIR}/sobol_complete.json", 'w') as f:
        json.dump(completion, f, indent=2)

    print(f"\n{'='*60}")
    print(f"COMPLETE: {total_done}/{total_runs} runs")
    print(f"Errors: {n_errors} ({n_errors/total_runs*100:.1f}%)")
    print(f"Wall time: {elapsed/3600:.1f}h ({elapsed:.0f}s)")
    if elapsed > 0:
        print(f"Mean rate: {total_done/elapsed:.2f} runs/s")
    print(f"Finished: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Results: {RESULTS_DIR}/sobol_complete.json")


if __name__ == "__main__":
    main()
