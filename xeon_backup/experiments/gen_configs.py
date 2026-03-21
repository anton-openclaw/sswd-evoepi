import json, os

OUT_DIR = os.path.expanduser("~/projects/sswd-evoepi/experiments/configs/W229-W248")
os.makedirs(OUT_DIR, exist_ok=True)

# W201 baseline (27 params from W201 log, the r_total=0.003 best run)
BASE = {
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
    "population.settler_survival": 1.0,
    "spatial.alpha_self_fjord": 0.7,
    "spatial.alpha_self_open": 0.02,
    "spatial.n_connectivity": 0.3,
    "spatial.r_total": 0.003,
}

configs = []

def cfg(name, fw, exp=1.0, lat_min=48.0, lat_max=60.0, alpha_env=None,
         phi_fjord=None, n_conn=None, notes=""):
    c = dict(BASE)
    c["disease.fw_strength"] = fw
    c["disease.fw_depth_exp"] = exp
    c["disease.fw_lat_min"] = lat_min
    c["disease.fw_lat_max"] = lat_max
    if alpha_env is not None:
        c["disease.alpha_env"] = alpha_env
    # else keeps W201 baseline alpha_env=0.18
    if phi_fjord is not None:
        c["spatial.phi_fjord"] = phi_fjord
    if n_conn is not None:
        c["spatial.n_connectivity"] = n_conn
    configs.append((name, c, notes))

# Group A: fw_strength dose-response, f_melt(48-60), W201 baseline
cfg("W229", fw=8,  notes="Conservative fw")
cfg("W230", fw=12, notes="Moderate fw")
cfg("W231", fw=16, notes="Strong fw")
cfg("W232", fw=20, notes="Aggressive fw")
cfg("W233", fw=25, notes="Upper bound fw")

# Group B: Latitude bound sensitivity
cfg("W234", fw=12, lat_min=50, lat_max=60, notes="Steep glacial 50-60")
cfg("W235", fw=16, lat_min=50, lat_max=60, notes="Steep glacial 50-60")
cfg("W236", fw=20, lat_min=50, lat_max=60, notes="Steep glacial 50-60")
cfg("W237", fw=16, lat_min=45, lat_max=58, notes="Broader lower ceiling")

# Group C: fw_depth_exp
cfg("W238", fw=16, exp=0.5, notes="Broader fjord protection sqrt")
cfg("W239", fw=16, exp=1.5, notes="Deep fjords only")
cfg("W240", fw=20, exp=0.5, notes="Strong plus broad")

# Group D: Synergy with alpha_env=0.25
cfg("W241", fw=12, alpha_env=0.25, notes="Moderate fw plus alpha_env")
cfg("W242", fw=16, alpha_env=0.25, notes="Strong fw plus alpha_env")
cfg("W243", fw=20, alpha_env=0.25, notes="Aggressive fw plus alpha_env")
cfg("W244", fw=16, exp=0.5, alpha_env=0.25, notes="Broad plus alpha_env")
cfg("W245", fw=20, exp=0.5, alpha_env=0.25, notes="Broad aggressive plus alpha_env")

# Group E: Full stack
cfg("W246", fw=16, exp=0.5, alpha_env=0.25, phi_fjord=0.15, notes="Kitchen sink moderate")
cfg("W247", fw=20, exp=0.5, alpha_env=0.25, phi_fjord=0.15, notes="Kitchen sink aggressive")
cfg("W248", fw=16, alpha_env=0.25, lat_min=50, lat_max=60, n_conn=0.5, notes="Steep plus connected")

for name, params, notes in configs:
    fpath = os.path.join(OUT_DIR, name + ".json")
    data = {"param_overrides": params}
    if notes:
        data["_notes"] = notes
    with open(fpath, "w") as f:
        json.dump(data, f, indent=2)
    fw = params["disease.fw_strength"]
    exp = params["disease.fw_depth_exp"]
    lmin = params["disease.fw_lat_min"]
    lmax = params["disease.fw_lat_max"]
    aenv = params["disease.alpha_env"]
    extra = ""
    if "spatial.phi_fjord" in params:
        extra += " phi_fjord=" + str(params["spatial.phi_fjord"])
    nc = params.get("spatial.n_connectivity", 0.3)
    if nc != 0.3:
        extra += " n_conn=" + str(nc)
    print("%s: fw=%s, exp=%s, lat=[%s-%s], aenv=%s%s  # %s" % (
        name, fw, exp, lmin, lmax, aenv, extra, notes))

print("\n%d configs written (with %d base params each)" % (len(configs), len(BASE)))
