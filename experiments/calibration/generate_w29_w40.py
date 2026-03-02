#!/usr/bin/env python3
"""Generate W29-W40 calibration configs.

Design: activation_threshold [50, 100] × K_half [200K, 400K, 600K] × T_vbnc [12, 14]
= 2 × 3 × 2 = 12 runs

Fixed: k_vbnc=2.0, s0=0.002, P_env_max=2000, D_P=50, D_P_max_range=500,
       5 Channel Islands origins, wavefront enabled

Changes from W17-W28:
- activation_threshold: 500 → [50, 100]  (let wavefront reach Alaska)
- D_P_max_range: 175 → 500 km  (longer-range waterborne dispersal)
- K_half: explore [200K, 400K, 600K]  (allow partial recovery)
- T_vbnc: explore [12, 14]  (shift VBNC midpoint for more winter relief)
"""
import json
import os

ORIGIN_NODES = [322, 319, 632, 633, 634]  # Channel Islands

configs = []
w = 29
for act_thresh in [50, 100]:
    for k_half in [200000, 400000, 600000]:
        for t_vbnc in [12.0, 14.0]:
            cfg = {
                "disease.K_half": float(k_half),
                "disease.P_env_max": 2000.0,
                "disease.k_vbnc": 2.0,
                "disease.T_vbnc": t_vbnc,
                "disease.activation_threshold": float(act_thresh),
                "disease.wavefront_enabled": True,
                "disease.disease_origin_nodes": ORIGIN_NODES,
                "population.settler_survival": 0.002,
                "spatial.D_P": 50.0,
                "spatial.D_P_max_range": 500.0,
            }
            fname = f"W{w}_config.json"
            with open(fname, 'w') as f:
                json.dump(cfg, f, indent=2)
            print(f"W{w}: act={act_thresh}, K_half={k_half//1000}K, T_vbnc={t_vbnc}")
            configs.append((w, cfg))
            w += 1

print(f"\nGenerated {len(configs)} configs (W29-W{w-1})")
