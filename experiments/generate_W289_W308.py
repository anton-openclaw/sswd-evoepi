#!/usr/bin/env python3
"""Generate W289-W308 calibration configs.

Unified sweep synthesized from three council recommendations:
  - Calibration strategist (epidemiological levers)
  - Spatial modeler (wavefront/connectivity levers)
  - Disease ecologist (reservoir/dormancy levers)

Core diagnosis: The latitudinal recovery gradient is too flat because all
previously tested levers act uniformly across latitudes. This sweep targets
temperature-gated and latitude-dependent mechanisms.

Excluded: delta_env changes (deemed a hack by Willem).

Baseline: W285 (RMSLE=0.666, AK-PWS=19.5%)
Target: RMSLE < 0.60, AK-PWS >= 18%, OR <= 1.0%, CA-N <= 0.5%
"""

import json
from pathlib import Path
from copy import deepcopy

# W285 baseline — all configs start here
BASELINE = {
    "K": 5000,
    "K_cv": 0.0,
    "disease_year": 0,
    "seeds": [137],
    "sst_start_year": 2012,
    "years": 13,
    "param_overrides": {
        "disease.K_half": 1500000,
        "disease.P_env_dynamic": True,
        "disease.P_env_floor": 500.0,
        "disease.P_env_max": 2000.0,
        "disease.T_vbnc": 12.0,
        "disease.T_vbnc_initial": 12.0,
        "disease.T_vbnc_min": 10.0,
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
        "disease.fw_strength": 25,
        "disease.fw_depth_exp": 1.0,
        "disease.fw_lat_min": 48.0,
        "disease.fw_lat_max": 60.0,
        "disease.seed_vibrio": 2000.0,
        "population.settler_survival": 1.0,
        "spatial.n_connectivity": 0.3,
        "spatial.alpha_self_open": 0.02,
        "spatial.alpha_self_fjord": 0.7,
        "spatial.r_total": 0.003,
    },
}


def make_config(name: str, overrides: dict, description: str) -> dict:
    """Create a config by applying overrides to the W285 baseline."""
    cfg = deepcopy(BASELINE)
    for k, v in overrides.items():
        cfg["param_overrides"][k] = v
    cfg["_name"] = name
    cfg["_description"] = description
    return cfg


configs = []

# ═══════════════════════════════════════════════════════════════════════
# BLOCK A: Single-lever isolates (W289-W295)
# Measure the standalone gradient effect of each promising parameter.
# ═══════════════════════════════════════════════════════════════════════

configs.append(make_config("W289",
    {"disease.v_max_warm": 0.85},
    "v_max_warm=0.85 — increase warm-water virulence by 21%. "
    "Southern SST 12-16°C hits virulence ceiling; northern SST 6-10°C stays below ramp. "
    "Unanimous #1 pick across all three councils."))

configs.append(make_config("W290",
    {"disease.v_max_warm": 1.0},
    "v_max_warm=1.0 — aggressive upper bound. Tests whether near-total "
    "warm-water lethality can suppress southern recovery to target levels."))

configs.append(make_config("W291",
    {"disease.pathogen_adapt_rate": 0.0003},
    "pathogen_adapt_rate=0.0003 — 3x slower cold adaptation. "
    "Pathogen reaches ~12.5°C instead of 9°C floor over 13 years. "
    "Northern populations retain wider VBNC dormancy window."))

configs.append(make_config("W292",
    {"disease.k_vbnc": 4.0},
    "k_vbnc=4.0 — sharper VBNC transition (2x steeper). "
    "Dormancy becomes near-binary: fully active or fully dormant. "
    "Sites near T_vbnc threshold see amplified protective effect."))

configs.append(make_config("W293",
    {"disease.T_vbnc_min": 11.5},
    "T_vbnc_min=11.5 — higher dormancy floor. Even fully cold-adapted pathogen "
    "enters dormancy below 11.5°C. Alaska gets 7-8mo disease-free; OR gets 4mo; "
    "southern CA gets almost none."))

configs.append(make_config("W294",
    {"disease.wavefront_D_P": 150.0},
    "wavefront_D_P=150km — halve dispersal kernel. Roughly doubles wavefront "
    "transit time to Alaska, giving northern populations an extra reproductive "
    "season before crash."))

configs.append(make_config("W295",
    {"disease.P_env_max": 4000.0},
    "P_env_max=4000 — 2x environmental reservoir capacity. "
    "More Vibrio accumulates in the water column. Should increase disease "
    "pressure everywhere but especially in warm water where Vibrio thrives."))

# ═══════════════════════════════════════════════════════════════════════
# BLOCK B: Two-lever compounds (W296-W301)
# Test synergies between the most promising individual levers.
# ═══════════════════════════════════════════════════════════════════════

configs.append(make_config("W296",
    {"disease.v_max_warm": 0.85, "disease.pathogen_adapt_rate": 0.0005},
    "v_max_warm + slow adaptation — dual gradient attack. Higher warm virulence "
    "(suppress south) AND slower cold adaptation (boost north). "
    "Calibration strategist's #1 compound recommendation."))

configs.append(make_config("W297",
    {"disease.v_max_warm": 0.85, "disease.k_vbnc": 4.0},
    "v_max_warm + sharp dormancy — virulence ceiling suppresses south; "
    "knife-edge VBNC transition amplifies northern dormancy protection."))

configs.append(make_config("W298",
    {"disease.T_vbnc_min": 11.5, "disease.pathogen_adapt_rate": 0.0003},
    "High dormancy floor + slow adaptation — compound dormancy enhancement. "
    "Pathogen stays at ~12.5°C, can't reach 11.5°C floor. Northern populations "
    "retain maximum disease-free window for entire simulation."))

configs.append(make_config("W299",
    {"disease.T_vbnc_min": 12.0, "disease.k_vbnc": 4.0},
    "T_vbnc_min=12 + k_vbnc=4.0 — high floor + sharp transition. "
    "Disease ecologist Config 8. Dormancy below 12°C with sharp cutoff. "
    "OR (winter SST 8-9°C) and AK (4-6°C) gain substantial dormancy windows."))

configs.append(make_config("W300",
    {"disease.v_max_warm": 0.85, "disease.wavefront_D_P": 200.0},
    "v_max_warm + slower wavefront — cross-domain compound. Temperature-gated "
    "virulence plus delayed Alaska arrival. Tests additive gradient effects "
    "from independent mechanisms."))

configs.append(make_config("W301",
    {"disease.wavefront_D_P": 150.0, "disease.cumulative_dose_threshold": 2000.0},
    "Slow wavefront + high CDT — spatial compound. Delayed pathogen arrival "
    "AND higher activation threshold. Cold-water dose accumulation is slower, "
    "so CDT creates an implicit latitude-dependent lethality gradient."))

# ═══════════════════════════════════════════════════════════════════════
# BLOCK C: Three-lever compounds (W302-W305)
# Stack the most promising mechanisms for synergistic gradient effects.
# ═══════════════════════════════════════════════════════════════════════

configs.append(make_config("W302",
    {"disease.v_max_warm": 0.85, "disease.pathogen_adapt_rate": 0.0003,
     "disease.k_vbnc": 4.0},
    "v_max_warm + slow adapt + sharp dormancy — triple disease ecology stack. "
    "Suppresses south via virulence, boosts north via dormancy sharpness "
    "and slow pathogen evolution."))

configs.append(make_config("W303",
    {"disease.v_max_warm": 0.85, "disease.T_vbnc_min": 11.5,
     "disease.k_vbnc": 4.0},
    "v_max_warm + high floor + sharp dormancy — temperature-gated trifecta. "
    "All three levers are SST-mediated through different mechanisms."))

configs.append(make_config("W304",
    {"disease.v_max_warm": 0.85, "disease.pathogen_adapt_rate": 0.0003,
     "disease.wavefront_D_P": 200.0},
    "v_max_warm + slow adapt + slower wavefront — cross-domain triple. "
    "Combines best epidemiological and spatial levers."))

configs.append(make_config("W305",
    {"disease.v_max_warm": 0.85, "disease.T_vbnc_min": 11.5,
     "disease.fw_lat_min": 46.0, "disease.fw_lat_max": 63.0},
    "v_max_warm + high dormancy + broad freshwater — multi-mechanism gradient. "
    "Freshwater band now covers Columbia River (46°N) to upper Alaska (63°N)."))

# ═══════════════════════════════════════════════════════════════════════
# BLOCK D: Multi-lever best guesses (W306-W308)
# Informed combinations targeting RMSLE < 0.60.
# ═══════════════════════════════════════════════════════════════════════

configs.append(make_config("W306",
    {"disease.v_max_warm": 0.85, "disease.T_vbnc_min": 11.5,
     "disease.k_vbnc": 4.0, "disease.pathogen_adapt_rate": 0.0003},
    "'Goldilocks' — disease ecology best guess (adapted from disease ecologist "
    "Config 16, delta_env removed). Four temperature-gated levers at moderate "
    "intensity. Predicted: OR ~0.5-1.0%, CA-N ~0.3-0.8%, AK-PWS ~15-20%."))

configs.append(make_config("W307",
    {"disease.v_max_warm": 0.9, "disease.T_vbnc_min": 12.0,
     "disease.k_vbnc": 4.5, "disease.pathogen_adapt_rate": 0.0002},
    "'Aggressive' — stress test for gradient steepness. All disease ecology "
    "levers pushed to plausible extremes. Risk: may over-recover in north "
    "or suppress BC-N (summer SST 10-14°C near T_vbnc threshold)."))

configs.append(make_config("W308",
    {"disease.v_max_warm": 0.85, "disease.T_vbnc_min": 11.5,
     "disease.k_vbnc": 4.0, "disease.pathogen_adapt_rate": 0.0003,
     "disease.K_half": 1800000,
     "disease.wavefront_D_P": 200.0,
     "disease.fw_lat_min": 46.0, "disease.fw_lat_max": 63.0},
    "'Kitchen sink' — all promising levers combined. Disease ecology Goldilocks "
    "plus K_half bump, slower wavefront, broader freshwater. "
    "If this doesn't steepen the gradient, spatial/disease mechanisms alone "
    "are insufficient and we need structural model changes."))

# ═══════════════════════════════════════════════════════════════════════
# Write configs
# ═══════════════════════════════════════════════════════════════════════

out_dir = Path("experiments/calibration")
out_dir.mkdir(parents=True, exist_ok=True)

for cfg in configs:
    name = cfg.pop("_name")
    desc = cfg.pop("_description")
    path = out_dir / f"{name}_config.json"
    with open(path, "w") as f:
        json.dump(cfg, f, indent=2)
    print(f"  {name}: {desc[:80]}...")

print(f"\nGenerated {len(configs)} configs in {out_dir}/")

# Also write a summary for reference
summary_path = out_dir / "W289_W308_summary.md"
with open(summary_path, "w") as f:
    f.write("# W289-W308 Sweep Summary\n\n")
    f.write("Baseline: W285 (RMSLE=0.666, AK-PWS=19.5%)\n")
    f.write("Target: RMSLE < 0.60, AK-PWS >= 18%, OR <= 1.0%, CA-N <= 0.5%\n\n")
    f.write("| Config | Changes vs W285 | Block |\n")
    f.write("|--------|-----------------|-------|\n")
    for i, cfg in enumerate(configs):
        name = f"W{289+i}"
        # Rebuild overrides description
        changed = {}
        for k, v in cfg["param_overrides"].items():
            bv = BASELINE["param_overrides"].get(k)
            if bv != v:
                short_k = k.replace("disease.", "").replace("spatial.", "")
                changed[short_k] = v
        block = "A:Isolate" if i < 7 else "B:Pair" if i < 13 else "C:Triple" if i < 17 else "D:Best"
        changes_str = ", ".join(f"{k}={v}" for k, v in changed.items())
        f.write(f"| {name} | {changes_str} | {block} |\n")

print(f"Summary: {summary_path}")
