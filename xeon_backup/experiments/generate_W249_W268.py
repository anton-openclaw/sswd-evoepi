#!/usr/bin/env python3
"""Generate W249-W268: First clean sweep after audit bug fixes + tooling consolidation.

This is the first calibration sweep with:
- 4 bug fixes (immunosuppression, ghost cascade, K overshoot, growth noise)
- Two-layer salinity model (WOA23 + fjord freshwater)
- Config defaults = W201 baseline
- Consolidated metrics (log10 RMSE, single source of truth)

Sweep strategy:
- W249: Pure baseline (W201 params, fw_strength=0) → measure bug fix impact
- W250-W253: fw_strength sweep (15, 20, 25, 30)
- W254-W255: alpha_env exploration with fw=25
- W256-W257: K_half exploration with fw=25
- W258-W259: alpha_env × K_half cross with fw=25
- W260-W263: n_connectivity exploration with fw=25
- W264-W265: alpha_self_fjord with fw=25
- W266-W268: r_total fine-tuning with fw=25

All configs: 1 seed (137), K=5000, 13 years. Seed CV < 1% from prior sweeps.
"""

import json
from pathlib import Path

OUTPUT_DIR = Path(__file__).parent / "calibration"
OUTPUT_DIR.mkdir(exist_ok=True)

# W201 baseline — every param explicit (no reliance on runner defaults)
BASELINE = {
    "disease.K_half": 800000.0,
    "disease.P_env_dynamic": True,
    "disease.P_env_floor": 500.0,
    "disease.P_env_max": 2000.0,
    "disease.T_vbnc": 12.0,
    "disease.T_vbnc_initial": 12.0,
    "disease.T_vbnc_min": 9.0,
    "disease.activation_threshold": 50.0,
    "disease.alpha_env": 0.18,
    "disease.cumulative_dose_threshold": 1000.0,
    "disease.delta_env": 0.02,
    "disease.disease_origin_nodes": [322, 319, 632, 633, 634],
    "disease.k_vbnc": 2.0,
    "disease.pathogen_adapt_rate": 0.001,
    "disease.pathogen_adaptation": True,
    "disease.pathogen_revert_rate": 0.0,
    "disease.v_adapt_rate": 0.001,
    "disease.v_max_warm": 0.7,
    "disease.virulence_evolution": True,
    "disease.wavefront_D_P": 300.0,
    "disease.wavefront_D_P_max_range": 3000.0,
    "disease.wavefront_enabled": True,
    "disease.fw_strength": 0.0,
    "disease.fw_depth_exp": 1.0,
    "disease.fw_lat_min": 48.0,
    "disease.fw_lat_max": 60.0,
    "population.settler_survival": 1.0,
    "spatial.n_connectivity": 0.3,
    "spatial.alpha_self_open": 0.02,
    "spatial.alpha_self_fjord": 0.7,
    "spatial.r_total": 0.003,
}

COMMON = {
    "K": 5000,
    "K_cv": 0.0,
    "years": 13,
    "disease_year": 0,
    "seeds": [137],
    "sst_start_year": 2012,
}


def make_config(name: str, overrides: dict, description: str = ""):
    """Generate a config file with W201 baseline + overrides."""
    params = dict(BASELINE)
    params.update(overrides)
    config = {"param_overrides": params, **COMMON}
    
    path = OUTPUT_DIR / f"{name}_config.json"
    with open(path, "w") as f:
        json.dump(config, f, indent=2)
    print(f"  {name}: {description}")
    return path


configs = []

print("=== W249-W268 Clean Sweep ===\n")

# Phase 1: Bug-fix baseline
print("Phase 1: Baseline")
configs.append(make_config("W249", {},
    "W201 baseline (bug fixes only, no salinity)"))

# Phase 2: fw_strength sweep
print("\nPhase 2: Salinity (fw_strength)")
for i, fw in enumerate([15, 20, 25, 30]):
    configs.append(make_config(f"W{250+i}", {"disease.fw_strength": fw},
        f"fw_strength={fw}"))

# Phase 3: alpha_env with optimal salinity
print("\nPhase 3: alpha_env × fw=25")
configs.append(make_config("W254", {"disease.fw_strength": 25, "disease.alpha_env": 0.22},
    "fw=25, alpha_env=0.22"))
configs.append(make_config("W255", {"disease.fw_strength": 25, "disease.alpha_env": 0.25},
    "fw=25, alpha_env=0.25"))

# Phase 4: K_half with fw=25
print("\nPhase 4: K_half × fw=25")
configs.append(make_config("W256", {"disease.fw_strength": 25, "disease.K_half": 400000},
    "fw=25, K_half=400K"))
configs.append(make_config("W257", {"disease.fw_strength": 25, "disease.K_half": 1200000},
    "fw=25, K_half=1.2M"))

# Phase 5: alpha_env × K_half cross
print("\nPhase 5: alpha_env × K_half × fw=25")
configs.append(make_config("W258",
    {"disease.fw_strength": 25, "disease.alpha_env": 0.25, "disease.K_half": 400000},
    "fw=25, alpha=0.25, K=400K"))
configs.append(make_config("W259",
    {"disease.fw_strength": 25, "disease.alpha_env": 0.25, "disease.K_half": 1200000},
    "fw=25, alpha=0.25, K=1.2M"))

# Phase 6: n_connectivity
print("\nPhase 6: n_connectivity × fw=25")
configs.append(make_config("W260", {"disease.fw_strength": 25, "spatial.n_connectivity": 0.2},
    "fw=25, n_conn=0.2"))
configs.append(make_config("W261", {"disease.fw_strength": 25, "spatial.n_connectivity": 0.5},
    "fw=25, n_conn=0.5"))
configs.append(make_config("W262",
    {"disease.fw_strength": 25, "disease.alpha_env": 0.25, "spatial.n_connectivity": 0.2},
    "fw=25, alpha=0.25, n_conn=0.2"))
configs.append(make_config("W263",
    {"disease.fw_strength": 25, "disease.alpha_env": 0.25, "spatial.n_connectivity": 0.5},
    "fw=25, alpha=0.25, n_conn=0.5"))

# Phase 7: alpha_self_fjord
print("\nPhase 7: alpha_self_fjord × fw=25")
configs.append(make_config("W264",
    {"disease.fw_strength": 25, "spatial.alpha_self_fjord": 0.5},
    "fw=25, alpha_fjord=0.5"))
configs.append(make_config("W265",
    {"disease.fw_strength": 25, "spatial.alpha_self_fjord": 0.9},
    "fw=25, alpha_fjord=0.9"))

# Phase 8: r_total fine-tuning
print("\nPhase 8: r_total × fw=25")
configs.append(make_config("W266", {"disease.fw_strength": 25, "spatial.r_total": 0.002},
    "fw=25, r_total=0.002"))
configs.append(make_config("W267", {"disease.fw_strength": 25, "spatial.r_total": 0.004},
    "fw=25, r_total=0.004"))
configs.append(make_config("W268",
    {"disease.fw_strength": 25, "disease.alpha_env": 0.25, "spatial.r_total": 0.002},
    "fw=25, alpha=0.25, r_total=0.002"))

print(f"\n=== Generated {len(configs)} configs ===")
