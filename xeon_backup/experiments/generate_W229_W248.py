import json, os

OUT_DIR = os.path.expanduser("~/projects/sswd-evoepi/experiments/configs/W229-W248")
os.makedirs(OUT_DIR, exist_ok=True)

# Shared baseline: r_total=0.003, s0=1.0, K=5000, 13 years, disease_year=0
BASE = {
    "spatial.r_total": 0.003,
    "population.settler_survival": 1.0,
    "simulation.years": 13,
    "simulation.disease_year": 0,
    "simulation.K": 5000,
    "disease.P_env_dynamic": True,
}

configs = []

def cfg(name, fw, exp=1.0, lat_min=48.0, lat_max=60.0, alpha_env=0.10,
         phi_fjord=None, n_conn=None, notes=""):
    c = dict(BASE)
    c["disease.fw_strength"] = fw
    c["disease.fw_depth_exp"] = exp
    c["disease.fw_lat_min"] = lat_min
    c["disease.fw_lat_max"] = lat_max
    c["disease.alpha_env"] = alpha_env
    if phi_fjord is not None:
        c["spatial.phi_fjord"] = phi_fjord
    if n_conn is not None:
        c["spatial.n_connectivity"] = n_conn
    configs.append((name, c, notes))

# Group A: fw_strength dose-response, f_melt(48-60)
cfg("W229", fw=8,  notes="Conservative fw")
cfg("W230", fw=12, notes="Moderate fw")
cfg("W231", fw=16, notes="Strong fw")
cfg("W232", fw=20, notes="Aggressive fw")
cfg("W233", fw=25, notes="Upper bound fw")

# Group B: Latitude bound sensitivity
cfg("W234", fw=12, lat_min=50, lat_max=60, notes="Steep glacial 50-60")
cfg("W235", fw=16, lat_min=50, lat_max=60, notes="Steep glacial 50-60")
cfg("W236", fw=20, lat_min=50, lat_max=60, notes="Steep glacial 50-60")
cfg("W237", fw=16, lat_min=45, lat_max=58, notes="Broader, lower ceiling")

# Group C: fw_depth_exp
cfg("W238", fw=16, exp=0.5, notes="Broader fjord protection (sqrt)")
cfg("W239", fw=16, exp=1.5, notes="Deep fjords only")
cfg("W240", fw=20, exp=0.5, notes="Strong + broad")

# Group D: Synergy with alpha_env=0.25
cfg("W241", fw=12, alpha_env=0.25, notes="Moderate fw + alpha_env")
cfg("W242", fw=16, alpha_env=0.25, notes="Strong fw + alpha_env")
cfg("W243", fw=20, alpha_env=0.25, notes="Aggressive fw + alpha_env")
cfg("W244", fw=16, exp=0.5, alpha_env=0.25, notes="Broad + alpha_env")
cfg("W245", fw=20, exp=0.5, alpha_env=0.25, notes="Broad aggressive + alpha_env")

# Group E: Full stack
cfg("W246", fw=16, exp=0.5, alpha_env=0.25, phi_fjord=0.15, notes="Kitchen sink moderate")
cfg("W247", fw=20, exp=0.5, alpha_env=0.25, phi_fjord=0.15, notes="Kitchen sink aggressive")
cfg("W248", fw=16, alpha_env=0.25, lat_min=50, lat_max=60, n_conn=0.5, notes="Steep + connected")

# Write configs
for name, params, notes in configs:
    fpath = os.path.join(OUT_DIR, f"{name}.json")
    data = {"param_overrides": params}
    if notes:
        data["_notes"] = notes
    with open(fpath, "w") as f:
        json.dump(data, f, indent=2)
    print(f"{name}: fw={params[disease.fw_strength]}, exp={params[disease.fw_depth_exp]}, "
          f"lat=[{params[disease.fw_lat_min]}-{params[disease.fw_lat_max]}], "
          f"aenv={params[disease.alpha_env]}"
          + (f", phi_fjord={params.get(spatial.phi_fjord,)}" if "spatial.phi_fjord" in params else "")
          + (f", n_conn={params.get(spatial.n_connectivity,)}" if "spatial.n_connectivity" in params else "")
          + f"  # {notes}")

print(f"\n{len(configs)} configs written to {OUT_DIR}")
