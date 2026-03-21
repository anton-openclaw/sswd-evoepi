#!/usr/bin/env python3
"""SA Round 3 — Morris screening (re-run with fixed metrics).

880 runs across 8 cores with checkpointing, then analysis + JSON export.
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
from SALib.sample import morris as morris_sample
from SALib.analyze import morris as morris_analyze

N_CORES = 8
BASE_SEED = 54321
N_TRAJECTORIES = 20
CHECKPOINT_INTERVAL = 200
RESULTS_DIR = "results/sensitivity_r3"


def main():
    names = get_param_names()
    problem = get_salib_problem(names)

    # Generate Morris samples
    np.random.seed(BASE_SEED)
    X = morris_sample.sample(problem, N=N_TRAJECTORIES, seed=BASE_SEED)
    total_runs = X.shape[0]

    print(f"SA Round 3 Morris — {total_runs} runs, {len(names)} params, {len(METRIC_NAMES)} metrics")
    print(f"Cores: {N_CORES}, checkpoint every {CHECKPOINT_INTERVAL}")
    print(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    sys.stdout.flush()

    # Check for checkpoint
    ckpt_path = f"{RESULTS_DIR}/morris_checkpoint.npz"
    if os.path.exists(ckpt_path):
        ckpt = np.load(ckpt_path)
        Y = ckpt['Y']
        completed_mask = ckpt['completed_mask']
        n_done = int(np.sum(completed_mask))
        print(f"Resuming from checkpoint: {n_done}/{total_runs}")
    else:
        Y = np.full((total_runs, len(METRIC_NAMES)), np.nan)
        completed_mask = np.zeros(total_runs, dtype=bool)
        n_done = 0

    remaining = [i for i in range(total_runs) if not completed_mask[i]]
    print(f"Runs remaining: {len(remaining)}")
    sys.stdout.flush()

    if len(remaining) == 0:
        print("All runs already complete, skipping to analysis.")
    else:
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

                if batch_done % CHECKPOINT_INTERVAL == 0:
                    elapsed = time.time() - t_start
                    rate = batch_done / elapsed
                    eta_h = (len(remaining) - batch_done) / rate / 3600 if rate > 0 else 0
                    np.savez(ckpt_path, Y=Y, X=X, completed_mask=completed_mask)
                    print(f"[{time.strftime('%H:%M:%S')}] {total_done}/{total_runs} "
                          f"({total_done/total_runs*100:.1f}%) | {rate:.2f} runs/s | "
                          f"ETA: {eta_h:.1f}h | errors: {n_errors}")
                    sys.stdout.flush()

        # Final checkpoint
        np.savez(ckpt_path, Y=Y, X=X, completed_mask=completed_mask)
        elapsed = time.time() - t_start
        total_done = n_done + batch_done
        print(f"\nRuns complete: {total_done}/{total_runs}, errors: {n_errors}, "
              f"wall time: {elapsed/60:.1f}min")

    # ── Analysis ──────────────────────────────────────────────────
    print("\nRunning Morris analysis...")
    valid_rows = ~np.isnan(Y).all(axis=1)
    n_valid = int(np.sum(valid_rows))
    print(f"Valid rows: {n_valid}/{total_runs}")

    morris_results = {}
    for j, metric in enumerate(METRIC_NAMES):
        y = Y[:, j]
        n_metric_valid = int(np.sum(~np.isnan(y)))
        if n_metric_valid < total_runs // 2:
            print(f"  {metric}: SKIPPED ({n_metric_valid} valid)")
            continue
        y_clean = y.copy()
        y_clean[np.isnan(y_clean)] = np.nanmedian(y_clean)

        # Check for zero variance
        if np.std(y_clean) < 1e-12:
            print(f"  {metric}: ZERO VARIANCE (all values = {y_clean[0]:.6f})")
            continue

        Si = morris_analyze.analyze(problem, X, y_clean, seed=BASE_SEED)
        morris_results[metric] = {
            'mu_star': [float(x) for x in Si['mu_star']],
            'sigma': [float(x) for x in Si['sigma']],
            'mu_star_conf': [float(x) for x in Si['mu_star_conf']],
            'names': list(Si['names']),
        }
        print(f"  {metric}: analyzed ({n_metric_valid} valid, std={np.std(y_clean):.4g})")

    with open(f'{RESULTS_DIR}/morris_results.json', 'w') as f:
        json.dump(morris_results, f, indent=2)
    print(f"\nSaved {len(morris_results)} metric analyses to morris_results.json")

    # Completion marker
    with open(f'{RESULTS_DIR}/morris_complete.json', 'w') as f:
        json.dump({
            'total_runs': total_runs,
            'valid': n_valid,
            'metrics_analyzed': len(morris_results),
            'finished_at': time.strftime('%Y-%m-%d %H:%M:%S'),
        }, f, indent=2)
    print(f"Done: {time.strftime('%Y-%m-%d %H:%M:%S')}")


if __name__ == "__main__":
    main()
