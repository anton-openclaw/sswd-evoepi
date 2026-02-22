#!/usr/bin/env python3
"""Benchmark the vectorized disease loop on 11-node network.

Runs 3 simulations (seeds 42, 43, 44) with 20 years, disease at year 2,
movement enabled, 1 substep/day. Reports wall-clock time per run.
"""
import json
import time
import sys
import os
import traceback

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from sswd_evoepi.config import default_config, validate_config
from sswd_evoepi.model import run_spatial_simulation
from scripts.sensitivity.spatial_runner import build_network

SEEDS = [42, 43, 44]
N_YEARS = 20
DISEASE_YEAR = 2

def main():
    # Build config with movement enabled, 1 substep/day
    cfg = default_config()
    cfg.movement.enabled = True
    cfg.movement.substeps_per_day = 1
    validate_config(cfg)

    results = []
    for seed in SEEDS:
        print(f"\n{'='*60}")
        print(f"Running seed={seed}, {N_YEARS} years, disease_year={DISEASE_YEAR}")
        print(f"{'='*60}")
        
        network = build_network(rng_seed=seed)
        
        try:
            t0 = time.perf_counter()
            result = run_spatial_simulation(
                config=cfg,
                network=network,
                n_years=N_YEARS,
                disease_year=DISEASE_YEAR,
                seed=seed,
            )
            elapsed = time.perf_counter() - t0
            
            info = {
                "seed": seed,
                "elapsed_s": round(elapsed, 2),
                "initial_pop": int(result.initial_total_pop),
                "final_pop": int(result.final_total_pop),
                "status": "ok",
            }
            print(f"  Done in {elapsed:.2f}s  |  pop {info['initial_pop']} → {info['final_pop']}")
            
        except Exception as e:
            elapsed = time.perf_counter() - t0
            tb = traceback.format_exc()
            info = {
                "seed": seed,
                "elapsed_s": round(elapsed, 2),
                "status": "error",
                "error": str(e),
                "traceback": tb,
            }
            print(f"  ERROR after {elapsed:.2f}s: {e}")
            print(tb)
        
        results.append(info)

    # Summary
    ok = [r for r in results if r["status"] == "ok"]
    mean_time = sum(r["elapsed_s"] for r in ok) / len(ok) if ok else 0
    
    xeon_baseline = 293.5  # pre-vectorization Xeon baseline
    local_est_baseline = 200.0  # estimated pre-vectorization local baseline
    
    summary = {
        "benchmark": "vectorized_disease_11node",
        "n_nodes": 11,
        "k_per_node": 5000,
        "n_years": N_YEARS,
        "disease_year": DISEASE_YEAR,
        "movement_enabled": True,
        "runs": results,
        "mean_time_s": round(mean_time, 2),
        "xeon_pre_vectorization_baseline_s": xeon_baseline,
        "local_est_pre_vectorization_baseline_s": local_est_baseline,
        "speedup_vs_xeon_baseline": round(xeon_baseline / mean_time, 2) if mean_time > 0 else None,
        "speedup_vs_local_est_baseline": round(local_est_baseline / mean_time, 2) if mean_time > 0 else None,
    }
    
    outpath = os.path.join(os.path.dirname(__file__), '..', 'results', 'vectorization_benchmark.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"\nResults saved to {outpath}")
    
    print(f"\n{'='*60}")
    print(f"SUMMARY")
    print(f"{'='*60}")
    print(f"  Runs OK:        {len(ok)}/{len(results)}")
    print(f"  Mean time:      {mean_time:.2f}s")
    print(f"  Xeon baseline:  {xeon_baseline}s → {summary['speedup_vs_xeon_baseline']}× faster")
    print(f"  Local baseline: ~{local_est_baseline}s → {summary['speedup_vs_local_est_baseline']}× faster")
    for r in results:
        t = r['elapsed_s']
        s = r['status']
        print(f"  seed={r['seed']}: {t:.2f}s ({s})")


if __name__ == "__main__":
    main()
