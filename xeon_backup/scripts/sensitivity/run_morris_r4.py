#!/usr/bin/env python3
"""SA Round 4 — Morris screening on 11-node stepping-stone network.

Key changes from R3:
- 11-node chain (was 3): Sitka→Ketchikan→Haida Gwaii→Bella Bella→
  Howe Sound→SJI→Westport→Newport→Crescent City→Fort Bragg→Monterey
- 47 params, 23 metrics (added total_recovery_events, recovery_rate)
- target_mean_c = 0.02 (was 0.08), rho_rec = 0.05
- Recovery recording fixed (yearly_recoveries in SpatialSimResult)
- Designed for Xeon (48 cores)
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

N_CORES = 48
BASE_SEED = 77777
N_TRAJECTORIES = 20
CHECKPOINT_INTERVAL = 200
RESULTS_DIR = "results/sensitivity_r4"


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)

    names = get_param_names()
    problem = get_salib_problem(names)

    # Generate Morris samples
    np.random.seed(BASE_SEED)
    X = morris_sample.sample(problem, N=N_TRAJECTORIES, seed=BASE_SEED)
    total_runs = X.shape[0]

    print(f"SA Round 4 Morris — {total_runs} runs, {len(names)} params, {len(METRIC_NAMES)} metrics")
    print(f"Network: 11-node stepping-stone chain, K=5000/node (55K total)")
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
        # Build work items
        work = [(idx, X[idx], names, BASE_SEED) for idx in remaining]

        # Run in parallel
        n_errors = 0
        t0 = time.time()
        batch_done = 0

        with Pool(N_CORES) as pool:
            for result in pool.imap_unordered(run_single_spatial, work):
                idx = result["run_index"]
                if result["error"] is not None:
                    n_errors += 1
                    print(f"  ERROR run {idx}: {result['error']}")
                else:
                    for mi, mname in enumerate(METRIC_NAMES):
                        Y[idx, mi] = result["metrics"].get(mname, np.nan)

                completed_mask[idx] = True
                batch_done += 1

                if batch_done % 20 == 0:
                    elapsed = time.time() - t0
                    rate = batch_done / elapsed
                    eta = (len(remaining) - batch_done) / rate if rate > 0 else 0
                    print(f"  {batch_done}/{len(remaining)} done, "
                          f"{elapsed:.0f}s elapsed, {rate:.1f} runs/s, "
                          f"ETA {eta:.0f}s, errors={n_errors}")
                    sys.stdout.flush()

                # Checkpoint
                if batch_done % CHECKPOINT_INTERVAL == 0:
                    np.savez(ckpt_path, Y=Y, completed_mask=completed_mask)
                    print(f"  Checkpoint saved ({int(np.sum(completed_mask))}/{total_runs})")

        # Final checkpoint
        np.savez(ckpt_path, Y=Y, completed_mask=completed_mask)
        elapsed_total = time.time() - t0
        print(f"\nAll runs complete: {elapsed_total:.0f}s ({elapsed_total/60:.1f} min), errors={n_errors}")

    # ── Analysis ─────────────────────────────────────────────────────
    print("\n=== Morris Analysis ===")

    morris_results = {}
    for mi, mname in enumerate(METRIC_NAMES):
        y = Y[:, mi]
        valid = ~np.isnan(y)
        n_valid = int(np.sum(valid))

        if n_valid < total_runs * 0.8:
            print(f"  SKIP {mname}: only {n_valid}/{total_runs} valid")
            continue

        # Replace NaN with column median for analysis
        y_clean = y.copy()
        y_clean[np.isnan(y_clean)] = np.nanmedian(y_clean)

        try:
            result = morris_analyze.analyze(
                problem, X, y_clean,
                num_resamples=1000,
                conf_level=0.95,
                seed=BASE_SEED,
            )
            morris_results[mname] = {
                "mu_star": result["mu_star"].tolist(),
                "mu_star_conf": result["mu_star_conf"].tolist(),
                "sigma": result["sigma"].tolist(),
                "mu": result["mu"].tolist(),
                "n_valid": n_valid,
                "y_std": float(np.std(y_clean)),
                "y_mean": float(np.mean(y_clean)),
            }
        except Exception as e:
            print(f"  ERROR analyzing {mname}: {e}")

    # Save results
    output = {
        "round": 4,
        "n_params": len(names),
        "n_metrics": len(METRIC_NAMES),
        "n_trajectories": N_TRAJECTORIES,
        "total_runs": total_runs,
        "n_cores": N_CORES,
        "network": "11-node stepping-stone chain",
        "param_names": names,
        "metric_names": METRIC_NAMES,
        "results": morris_results,
    }

    with open(f"{RESULTS_DIR}/morris_results.json", "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved: {RESULTS_DIR}/morris_results.json")

    # ── Summary table ────────────────────────────────────────────────
    print("\n=== Top Parameters (mean μ* across metrics, normalized) ===")
    param_scores = np.zeros(len(names))
    n_metrics_used = 0

    for mname, mdata in morris_results.items():
        y_std = mdata["y_std"]
        if y_std < 1e-10:
            continue
        mu_stars = np.array(mdata["mu_star"])
        param_scores += mu_stars / y_std
        n_metrics_used += 1

    if n_metrics_used > 0:
        param_scores /= n_metrics_used
        ranking = np.argsort(param_scores)[::-1]
        print(f"\n{'Rank':>4} {'Parameter':<45} {'Norm μ*':>8}")
        print("-" * 60)
        for rank, idx in enumerate(ranking):
            print(f"{rank+1:>4} {names[idx]:<45} {param_scores[idx]:>8.4f}")

    print(f"\nDone. Results at {RESULTS_DIR}/")


if __name__ == "__main__":
    main()
