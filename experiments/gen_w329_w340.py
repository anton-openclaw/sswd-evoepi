#!/usr/bin/env python3
"""Generate W329-W340 configs — tuning around W325 baseline."""
import json, os

OUT = "/home/starbot/.openclaw/workspace/sswd-evoepi/experiments/calibration"

# W325 baseline
BASE = {
    "disease.K_half": 800000.0,
    "disease.T_vbnc_min": 10.0,
    "disease.seed_vibrio": 2000.0,
    "disease.fw_strength": 25.0,
    "disease.alpha_env": 0.18,
    "disease.k_vbnc": 2.0,
    "disease.pathogen_adapt_rate": 0.001,
    "disease.v_max_warm": 0.7,
    "disease.cumulative_dose_threshold": 1000.0,
    "disease.wavefront_D_P": 300.0,
    "spatial.n_connectivity": 0.3,
    "spatial.r_total": 0.003,
    "spatial.alpha_self_open": 0.05,
    "spatial.alpha_self_fjord": 0.5,
    "disease.delta_env_T_dependent": True,
    "disease.delta_env_cold": 0.05,
    "disease.delta_env_warm": 0.01,
    "disease.recovery_P_env_gated": True,
    "disease.recovery_P_env_half": 200.0,
    "disease.recovery_P_env_hill": 2.0,
}

# W325 problems: BC-N=4.3% (target 20%), AK-PWS=30% (target 20%), JDF=8.2% (target 2%), SS-S=8.6% (target 5%)
# Strategy: recovery_P_env_half ↑ (help BC-N), delta_env_cold ↓ (reduce AK overshoot), K_half bracket, hill ↓

configs = {
    # Single-parameter variations from W325
    "W329": {"disease.recovery_P_env_half": 400.0},      # relax P_env gating
    "W330": {"disease.recovery_P_env_half": 800.0},      # relax more
    "W331": {"disease.recovery_P_env_half": 1600.0},     # relax a lot
    "W332": {"disease.recovery_P_env_hill": 1.0},        # gentler sigmoid
    "W333": {"disease.delta_env_cold": 0.04},            # less pathogen death in cold → reduce AK
    "W334": {"disease.delta_env_cold": 0.03},            # even less → reduce AK more
    "W335": {"disease.K_half": 700000.0},                # more disease pressure
    "W336": {"disease.K_half": 900000.0},                # less disease pressure → help BC-N
    
    # Combination configs targeting BC-N + AK simultaneously
    "W337": {"disease.recovery_P_env_half": 400.0, "disease.delta_env_cold": 0.04},
    "W338": {"disease.recovery_P_env_half": 800.0, "disease.delta_env_cold": 0.04},
    "W339": {"disease.recovery_P_env_half": 400.0, "disease.K_half": 900000.0},
    "W340": {"disease.recovery_P_env_half": 800.0, "disease.recovery_P_env_hill": 1.0},
}

for name, overrides in configs.items():
    cfg = dict(BASE)
    cfg.update(overrides)
    path = os.path.join(OUT, f"{name}_config.json")
    with open(path, 'w') as f:
        json.dump({"param_overrides": cfg}, f, indent=2)
    # Show what changed
    changes = {k.split('.')[-1]: v for k, v in overrides.items()}
    print(f"{name}: {changes}")

print(f"\nGenerated {len(configs)} configs in {OUT}")
