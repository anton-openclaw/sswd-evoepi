#!/usr/bin/env python3
"""Generate W53-W64 calibration configs.

Design rationale:
- W42 base: CDT=1000, K_half=400K, T_vbnc=12, k_vbnc=2.0, s0=0.002, 5 Channel Islands origins
- NEW: wavefront_D_P and wavefront_D_P_max_range for long-range dispersal kernel
- Also test P_env_max increases to fix southern overshoot (from W45-W52 design)

Sweep:
  wavefront_D_P:           [150, 300]     (km; standard D_P is 50km)
  wavefront_D_P_max_range: [1500, 3000]   (km; standard is 500km)
  P_env_max:               [2000, 6000, 10000]  (fix southern overshoot)

= 2 × 2 × 3 = 12 configs (W53-W64)
"""
import json
import os

BASE = {
    "disease.K_half": 400000.0,
    "disease.P_env_max": 2000.0,
    "disease.k_vbnc": 2.0,
    "disease.T_vbnc": 12.0,
    "disease.activation_threshold": 50.0,
    "disease.wavefront_enabled": True,
    "disease.disease_origin_nodes": [322, 319, 632, 633, 634],
    "population.settler_survival": 0.002,
    "spatial.D_P": 50.0,
    "spatial.D_P_max_range": 500.0,
    "disease.cumulative_dose_threshold": 1000.0,
}

# Sweep values
wf_dp_values = [150, 300]
wf_max_range_values = [1500, 3000]
p_env_max_values = [2000, 6000, 10000]

configs = []
round_num = 53
for wf_dp in wf_dp_values:
    for wf_max in wf_max_range_values:
        for p_env in p_env_max_values:
            name = f"W{round_num:02d}"
            cfg = dict(BASE)
            cfg["disease.wavefront_D_P"] = float(wf_dp)
            cfg["disease.wavefront_D_P_max_range"] = float(wf_max)
            cfg["disease.P_env_max"] = float(p_env)
            configs.append((name, cfg))
            round_num += 1

# Write configs
out_dir = os.path.dirname(os.path.abspath(__file__))
for name, cfg in configs:
    path = os.path.join(out_dir, f"{name}_config.json")
    with open(path, "w") as f:
        json.dump(cfg, f, indent=2)
    print(f"  {name}: wf_D_P={cfg['disease.wavefront_D_P']:.0f}km, "
          f"wf_max={cfg['disease.wavefront_D_P_max_range']:.0f}km, "
          f"P_env_max={cfg['disease.P_env_max']:.0f}")

# Summary table
print(f"\n{len(configs)} configs written (W53-W{52+len(configs)})")
print("\nDesign matrix:")
print(f"  wavefront_D_P:           {wf_dp_values}")
print(f"  wavefront_D_P_max_range: {wf_max_range_values}")
print(f"  P_env_max:               {p_env_max_values}")
