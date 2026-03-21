#!/usr/bin/env python3
"""Compute Sobol indices from R4 raw output (sobol_Y_final.npz)."""
import json
import numpy as np
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from scripts.sensitivity.param_spec import get_param_names, get_salib_problem
from scripts.sensitivity.spatial_runner import METRIC_NAMES
from SALib.analyze import sobol as sobol_analyze

RESULTS_DIR = "results/sensitivity_r4"

def main():
    # Load raw data
    data = np.load(f"{RESULTS_DIR}/sobol_Y_final.npz")
    X = data["X"]
    Y = data["Y"]
    
    # Load metadata
    with open(f"{RESULTS_DIR}/sobol_complete.json") as f:
        meta = json.load(f)
    
    param_names = meta["param_names"]
    metric_names = meta["metric_names"]
    problem = get_salib_problem(param_names)
    
    print(f"X shape: {X.shape}, Y shape: {Y.shape}")
    print(f"Params: {len(param_names)}, Metrics: {len(metric_names)}")
    
    sobol_results = {}
    for j, metric in enumerate(metric_names):
        y = Y[:, j]
        n_valid = np.sum(~np.isnan(y))
        n_total = len(y)
        
        if n_valid < n_total * 0.5:
            print(f"  {metric}: skipped ({n_valid}/{n_total} valid)")
            continue
        
        y_clean = y.copy()
        nan_mask = np.isnan(y_clean)
        if np.any(nan_mask):
            y_clean[nan_mask] = np.nanmedian(y)
        
        if np.std(y_clean) < 1e-10:
            print(f"  {metric}: skipped (zero variance)")
            continue
        
        try:
            Si = sobol_analyze.analyze(
                problem, y_clean,
                calc_second_order=False,
                num_resamples=1000,
                conf_level=0.95,
                seed=54321
            )
            sobol_results[metric] = {
                "S1": Si["S1"].tolist(),
                "S1_conf": Si["S1_conf"].tolist(),
                "ST": Si["ST"].tolist(),
                "ST_conf": Si["ST_conf"].tolist(),
            }
            # Print top 3 for this metric
            top3 = np.argsort(Si["ST"])[::-1][:3]
            top_str = ", ".join(f"{param_names[i].split('.')[-1]}={Si['ST'][i]:.3f}" for i in top3)
            print(f"  {metric}: OK (top ST: {top_str})")
        except Exception as e:
            print(f"  {metric}: ERROR — {e}")
    
    # Save
    out = {
        "sobol_results": sobol_results,
        "param_names": param_names,
        "metric_names": metric_names,
        "n_runs": meta["total_runs"],
        "n_params": len(param_names),
    }
    outpath = f"{RESULTS_DIR}/sobol_results.json"
    with open(outpath, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved: {outpath}")
    
    # Print aggregate top-10
    print("\n" + "="*60)
    print("AGGREGATE TOP-10 (max ST across all metrics)")
    print("="*60)
    max_ST = np.zeros(len(param_names))
    max_S1 = np.zeros(len(param_names))
    best_metric = [""] * len(param_names)
    for metric, mdata in sobol_results.items():
        ST = np.array(mdata["ST"])
        S1 = np.array(mdata["S1"])
        for i in range(len(param_names)):
            if ST[i] > max_ST[i]:
                max_ST[i] = ST[i]
                best_metric[i] = metric
            if S1[i] > max_S1[i]:
                max_S1[i] = S1[i]
    
    ranked = np.argsort(max_ST)[::-1]
    print(f"\n{'Rank':>4} {'Parameter':<45} {'S1':>8} {'ST':>8} {'Best metric'}")
    print("-"*100)
    for rank, idx in enumerate(ranked[:10], 1):
        print(f"{rank:>4} {param_names[idx]:<45} {max_S1[idx]:>8.3f} {max_ST[idx]:>8.3f} {best_metric[idx]}")

if __name__ == "__main__":
    main()
