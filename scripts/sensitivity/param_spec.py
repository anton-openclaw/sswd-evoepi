"""
Parameter specification for sensitivity analysis.

Defines the 25 uncertain parameters, their ranges, distributions,
and how they map to SimulationConfig overrides.
"""

import numpy as np
from collections import OrderedDict


# ─── Parameter definitions ────────────────────────────────────────────
# Each entry: (config_path, low, high, distribution, description)
# config_path uses dot notation: "disease.a_exposure" → config.disease.a_exposure
# distribution: "uniform" (default), "loguniform", "discrete"

PARAM_SPEC = OrderedDict([
    # Disease (10)
    ("disease.a_exposure", {
        "low": 0.30, "high": 1.50, "dist": "uniform",
        "desc": "Exposure rate (d⁻¹)", "confidence": "★☆☆",
    }),
    ("disease.K_half", {
        "low": 20000.0, "high": 200000.0, "dist": "loguniform",
        "desc": "Half-infective dose (bact/mL)", "confidence": "★☆☆",
    }),
    ("disease.sigma_1_eff", {
        "low": 1.0, "high": 25.0, "dist": "loguniform",
        "desc": "I1 shedding rate (field-effective)", "confidence": "★☆☆",
    }),
    ("disease.sigma_2_eff", {
        "low": 10.0, "high": 250.0, "dist": "loguniform",
        "desc": "I2 shedding rate (field-effective)", "confidence": "★☆☆",
    }),
    ("disease.sigma_D", {
        "low": 3.0, "high": 75.0, "dist": "loguniform",
        "desc": "Saprophytic burst (field-effective)", "confidence": "★☆☆",
    }),
    ("disease.rho_rec", {
        "low": 0.0, "high": 0.20, "dist": "uniform",
        "desc": "Recovery rate (d⁻¹)", "confidence": "★☆☆",
    }),
    ("disease.mu_EI1_ref", {
        "low": 0.20, "high": 1.00, "dist": "uniform",
        "desc": "E→I1 progression rate at T_ref", "confidence": "★★☆",
    }),
    ("disease.mu_I2D_ref", {
        "low": 0.08, "high": 0.35, "dist": "uniform",
        "desc": "I2→Death rate at T_ref", "confidence": "★★☆",
    }),
    ("disease.P_env_max", {
        "low": 50.0, "high": 5000.0, "dist": "loguniform",
        "desc": "Background Vibrio input (bact/mL/d)", "confidence": "★☆☆",
    }),
    ("disease.T_ref", {
        "low": 17.0, "high": 23.0, "dist": "uniform",
        "desc": "V. pectenicida T_opt (°C)", "confidence": "★★☆",
    }),
    # Population (7)
    ("population.F0", {
        "low": 1e6, "high": 1e8, "dist": "loguniform",
        "desc": "Reference fecundity (eggs)", "confidence": "★☆☆",
    }),
    ("population.gamma_fert", {
        "low": 1.0, "high": 10.0, "dist": "uniform",
        "desc": "Fertilization kinetics parameter", "confidence": "★☆☆",
    }),
    ("population.settler_survival", {
        "low": 0.005, "high": 0.10, "dist": "loguniform",
        "desc": "B-H settler survival s0", "confidence": "★☆☆",
    }),
    ("population.alpha_srs", {
        "low": 1.0, "high": 1.8, "dist": "uniform",
        "desc": "SRS Pareto shape", "confidence": "★★☆",
    }),
    ("population.senescence_age", {
        "low": 20.0, "high": 80.0, "dist": "uniform",
        "desc": "Senescence onset age (yr)", "confidence": "★☆☆",
    }),
    ("population.k_growth", {
        "low": 0.03, "high": 0.15, "dist": "uniform",
        "desc": "VB growth rate (yr⁻¹)", "confidence": "★☆☆",
    }),
    ("population.L_min_repro", {
        "low": 200.0, "high": 500.0, "dist": "uniform",
        "desc": "Min reproductive size (mm)", "confidence": "★☆☆",
    }),
    # Genetics (1)
    # NOTE: s_het and q_ef1a_init (EF1A) dropped — Wares 2016 finding is
    # Pisaster ochraceus, not Pycnopodia. Overdominant locus is a structural
    # assumption to test separately, not a continuous parameter to sweep.
    ("genetics.n_additive", {
        "low": 10, "high": 51, "dist": "discrete",
        "values": [10, 20, 30, 40, 51],
        "desc": "Number of additive resistance loci", "confidence": "★★☆",
    }),
    # Spawning (3)
    ("spawning.p_spontaneous_female", {
        "low": 0.005, "high": 0.025, "dist": "uniform",
        "desc": "Daily spontaneous female spawning prob", "confidence": "★★☆",
    }),
    ("spawning.induction_female_to_male", {
        "low": 0.40, "high": 0.95, "dist": "uniform",
        "desc": "Female→male cascade induction", "confidence": "★★☆",
    }),
    ("disease.susceptibility_multiplier", {
        "low": 1.0, "high": 4.0, "dist": "uniform",
        "desc": "Post-spawning immunosuppression multiplier", "confidence": "★☆☆",
    }),
    # Environmental (2)
    ("disease.T_vbnc", {
        "low": 8.0, "high": 15.0, "dist": "uniform",
        "desc": "VBNC midpoint temperature (°C)", "confidence": "★★☆",
    }),
    ("disease.s_min", {
        "low": 5.0, "high": 15.0, "dist": "uniform",
        "desc": "Salinity minimum for Vibrio (psu)", "confidence": "★★☆",
    }),
])


def get_param_names():
    """Return ordered list of parameter names."""
    return list(PARAM_SPEC.keys())


def get_salib_problem(param_names=None):
    """Build SALib problem dict for given parameters (or all)."""
    if param_names is None:
        param_names = get_param_names()
    
    names = []
    bounds = []
    dists = []
    
    for name in param_names:
        spec = PARAM_SPEC[name]
        names.append(name)
        
        if spec["dist"] == "loguniform":
            # SALib expects bounds; we handle log-transform ourselves
            bounds.append([np.log10(spec["low"]), np.log10(spec["high"])])
            dists.append("unif")  # uniform in log-space
        elif spec["dist"] == "discrete":
            # Map to [0, len(values)-1] uniform, then discretize
            bounds.append([0.0, len(spec["values"]) - 0.001])
            dists.append("unif")
        else:
            bounds.append([spec["low"], spec["high"]])
            dists.append("unif")
    
    return {
        "num_vars": len(names),
        "names": names,
        "bounds": bounds,
        "dists": dists,
    }


def sample_to_config_overrides(sample_row, param_names=None):
    """Convert a single sample row to config override dict.
    
    Handles log-uniform back-transform and discrete parameter mapping.
    Returns a nested dict suitable for deep_merge with config YAML.
    """
    if param_names is None:
        param_names = get_param_names()
    
    overrides = {}
    
    for i, name in enumerate(param_names):
        spec = PARAM_SPEC[name]
        val = sample_row[i]
        
        # Transform
        if spec["dist"] == "loguniform":
            val = 10.0 ** val
        elif spec["dist"] == "discrete":
            idx = int(val)
            idx = min(idx, len(spec["values"]) - 1)
            val = spec["values"][idx]
        
        # Parse dotted path → nested dict
        parts = name.split(".")
        section = parts[0]
        field = parts[1]
        
        if section not in overrides:
            overrides[section] = {}
        
        # Type coercion
        if isinstance(val, (np.integer, np.int64)):
            val = int(val)
        elif isinstance(val, (np.floating, np.float64)):
            val = float(val)
        
        overrides[section][field] = val
    
    # Enforce constraint: n_loci = n_additive + 1
    if "genetics" in overrides and "n_additive" in overrides["genetics"]:
        overrides["genetics"]["n_loci"] = overrides["genetics"]["n_additive"] + 1
    
    return overrides


def apply_overrides_to_config(config, overrides):
    """Apply override dict to a SimulationConfig object in-place.
    
    overrides: {"disease": {"a_exposure": 0.5}, "population": {"F0": 1e7}, ...}
    """
    for section_name, fields in overrides.items():
        section = getattr(config, section_name)
        for field_name, value in fields.items():
            setattr(section, field_name, value)
    return config
