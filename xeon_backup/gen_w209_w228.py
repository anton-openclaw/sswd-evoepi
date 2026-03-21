import json, copy

# W201 base config
with open("experiments/calibration/W201_config.json") as f:
    base = json.load(f)

base["seeds"] = [137]

configs = {}

# --- H1: Fjord refuge ---
c = copy.deepcopy(base); c["param_overrides"]["spatial.phi_fjord"] = 0.15; configs["W209"] = c
c = copy.deepcopy(base); c["param_overrides"]["spatial.phi_fjord"] = 0.25; configs["W210"] = c
c = copy.deepcopy(base); c["param_overrides"]["spatial.phi_fjord"] = 0.15; c["param_overrides"]["spatial.alpha_self_fjord"] = 0.85; configs["W211"] = c

# --- H2: Southern persistence ---
c = copy.deepcopy(base); c["param_overrides"]["spatial.phi_open"] = 0.50; configs["W212"] = c
c = copy.deepcopy(base); c["param_overrides"]["spatial.phi_open"] = 0.50; c["param_overrides"]["spatial.phi_fjord"] = 0.15; configs["W213"] = c

# --- H3: Temperature amplification ---
c = copy.deepcopy(base); c["param_overrides"]["disease.alpha_env"] = 0.25; configs["W214"] = c
c = copy.deepcopy(base); c["param_overrides"]["disease.alpha_env"] = 0.30; configs["W215"] = c
c = copy.deepcopy(base); c["param_overrides"]["disease.alpha_env"] = 0.35; configs["W216"] = c

# --- H4: Disease gradient + fjord refuge ---
c = copy.deepcopy(base); c["param_overrides"]["disease.alpha_env"] = 0.25; c["param_overrides"]["spatial.phi_fjord"] = 0.15; configs["W217"] = c
c = copy.deepcopy(base); c["param_overrides"]["disease.alpha_env"] = 0.25; c["param_overrides"]["spatial.phi_open"] = 0.50; configs["W218"] = c
c = copy.deepcopy(base); c["param_overrides"]["disease.alpha_env"] = 0.25; c["param_overrides"]["spatial.phi_open"] = 0.60; c["param_overrides"]["spatial.phi_fjord"] = 0.15; configs["W219"] = c

# --- H5: Recruitment isolation ---
c = copy.deepcopy(base); c["param_overrides"]["spatial.alpha_self_fjord"] = 0.85; configs["W220"] = c
c = copy.deepcopy(base); c["param_overrides"]["spatial.alpha_self_fjord"] = 0.85; c["param_overrides"]["spatial.alpha_self_open"] = 0.01; configs["W221"] = c
c = copy.deepcopy(base); c["param_overrides"]["spatial.alpha_self_fjord"] = 0.85; c["param_overrides"]["disease.alpha_env"] = 0.25; configs["W222"] = c

# --- H6: Full asymmetry ---
c = copy.deepcopy(base); c["param_overrides"]["disease.alpha_env"] = 0.25; c["param_overrides"]["spatial.phi_open"] = 0.60; c["param_overrides"]["spatial.phi_fjord"] = 0.15; c["param_overrides"]["spatial.alpha_self_fjord"] = 0.80; configs["W223"] = c
c = copy.deepcopy(base); c["param_overrides"]["disease.alpha_env"] = 0.30; c["param_overrides"]["spatial.phi_open"] = 0.50; c["param_overrides"]["spatial.phi_fjord"] = 0.20; c["param_overrides"]["spatial.alpha_self_fjord"] = 0.85; configs["W224"] = c
c = copy.deepcopy(base); c["param_overrides"]["disease.alpha_env"] = 0.30; c["param_overrides"]["spatial.phi_open"] = 0.50; c["param_overrides"]["spatial.phi_fjord"] = 0.20; c["param_overrides"]["spatial.alpha_self_fjord"] = 0.85; c["param_overrides"]["spatial.n_connectivity"] = 0.5; configs["W225"] = c

# --- H7: Enclosedness amplification ---
c = copy.deepcopy(base); c["param_overrides"]["spatial.n_connectivity"] = 0.5; configs["W226"] = c
c = copy.deepcopy(base); c["param_overrides"]["spatial.n_connectivity"] = 0.5; c["param_overrides"]["disease.alpha_env"] = 0.25; configs["W227"] = c
c = copy.deepcopy(base); c["param_overrides"]["spatial.n_connectivity"] = 0.7; c["param_overrides"]["spatial.alpha_self_fjord"] = 0.85; configs["W228"] = c

for name, cfg in sorted(configs.items()):
    path = f"experiments/calibration/{name}_config.json"
    with open(path, "w") as f:
        json.dump(cfg, f, indent=2)
    p = cfg["param_overrides"]
    aenv = p.get("disease.alpha_env", 0.18)
    po = p.get("spatial.phi_open", 0.80)
    pf = p.get("spatial.phi_fjord", 0.03)
    asf = p.get("spatial.alpha_self_fjord", 0.70)
    nc = p.get("spatial.n_connectivity", 0.3)
    print(f"{name}: aenv={aenv}, phi_o={po}, phi_f={pf}, a_self_f={asf}, n_conn={nc}")

print(f"\nGenerated {len(configs)} configs")
