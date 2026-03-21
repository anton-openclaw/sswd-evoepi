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
        "low": 0.155, "high": 0.465, "dist": "uniform",
        "desc": "E→I1 progression rate at T_ref (Prentice 2025)", "confidence": "★★★",
    }),
    ("disease.mu_I2D_ref", {
        "low": 0.375, "high": 1.125, "dist": "uniform",
        "desc": "I2→Death rate at T_ref (Prentice 2025)", "confidence": "★★★",
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
    # Genetics — three-trait architecture (Phase 4)
    # n_resistance + n_tolerance + n_recovery = 51 (constrained)
    # Parameterize as two free variables; n_recovery = 51 - n_resistance - n_tolerance
    ("genetics.n_resistance", {
        "low": 5, "high": 30, "dist": "discrete",
        "values": [5, 10, 17, 25, 30],
        "desc": "Number of resistance loci", "confidence": "★★☆",
    }),
    ("genetics.n_tolerance", {
        "low": 5, "high": 30, "dist": "discrete",
        "values": [5, 10, 17, 25, 30],
        "desc": "Number of tolerance loci", "confidence": "★★☆",
    }),
    ("genetics.target_mean_t", {
        "low": 0.02, "high": 0.30, "dist": "uniform",
        "desc": "Target mean tolerance at t=0", "confidence": "★☆☆",
    }),
    ("genetics.target_mean_c", {
        "low": 0.02, "high": 0.25, "dist": "uniform",
        "desc": "Target mean recovery at t=0", "confidence": "★☆☆",
    }),
    ("genetics.tau_max", {
        "low": 0.3, "high": 0.95, "dist": "uniform",
        "desc": "Max tolerance mortality reduction", "confidence": "★☆☆",
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
    # ── NEW PARAMS (SA Round 2) ──────────────────────────────────────
    # Disease — missing middle transition + immunosuppression duration
    ("disease.mu_I1I2_ref", {
        "low": 0.289, "high": 0.867, "dist": "uniform",
        "desc": "I1→I2 progression rate at T_ref (Prentice 2025)", "confidence": "★★☆",
    }),
    ("disease.immunosuppression_duration", {
        "low": 7, "high": 56, "dist": "uniform",
        "desc": "Post-spawning immunosuppression duration (days)", "confidence": "★★☆",
    }),
    # Spawning — reverse induction
    ("spawning.induction_male_to_female", {
        "low": 0.10, "high": 0.60, "dist": "uniform",
        "desc": "Male→female cascade induction", "confidence": "★★☆",
    }),
    # Spatial — dispersal kernel
    ("spatial.D_L", {
        "low": 100.0, "high": 1000.0, "dist": "loguniform",
        "desc": "Larval dispersal scale (km)", "confidence": "★☆☆",
    }),
    # Genetics initialization
    ("genetics.target_mean_r", {
        "low": 0.05, "high": 0.30, "dist": "uniform",
        "desc": "Target mean resistance at t=0", "confidence": "★☆☆",
    }),
    ("genetics.q_init_beta_a", {
        "low": 1.0, "high": 5.0, "dist": "uniform",
        "desc": "Beta shape a for per-locus allele freq", "confidence": "★☆☆",
    }),
    ("genetics.q_init_beta_b", {
        "low": 3.0, "high": 15.0, "dist": "uniform",
        "desc": "Beta shape b for per-locus allele freq", "confidence": "★☆☆",
    }),
    # Larval retention
    ("spatial.alpha_self_fjord", {
        "low": 0.10, "high": 0.50, "dist": "uniform",
        "desc": "Larval self-recruitment fraction (fjord)", "confidence": "★☆☆",
    }),
    ("spatial.alpha_self_open", {
        "low": 0.02, "high": 0.20, "dist": "uniform",
        "desc": "Larval self-recruitment fraction (open coast)", "confidence": "★☆☆",
    }),
    # Pathogen evolution (6)
    ("pathogen_evolution.alpha_kill", {
        "low": 1.0, "high": 4.0, "dist": "uniform",
        "desc": "Death rate scaling exponent", "confidence": "★☆☆",
    }),
    ("pathogen_evolution.alpha_shed", {
        "low": 0.5, "high": 3.0, "dist": "uniform",
        "desc": "Shedding rate scaling exponent", "confidence": "★☆☆",
    }),
    ("pathogen_evolution.alpha_prog", {
        "low": 0.5, "high": 2.0, "dist": "uniform",
        "desc": "I1→I2 progression scaling exponent", "confidence": "★☆☆",
    }),
    ("pathogen_evolution.gamma_early", {
        "low": 0.0, "high": 1.0, "dist": "uniform",
        "desc": "Early shedding attenuation factor", "confidence": "★☆☆",
    }),
    ("pathogen_evolution.sigma_v_mutation", {
        "low": 0.005, "high": 0.10, "dist": "loguniform",
        "desc": "Mutation step size (std dev)", "confidence": "★☆☆",
    }),
    ("pathogen_evolution.v_init", {
        "low": 0.2, "high": 0.8, "dist": "uniform",
        "desc": "Initial pathogen virulence", "confidence": "★☆☆",
    }),
    # ── NEW PARAMS (SA Round 3 — Phase 11) ───────────────────────────
    # Juvenile immunity
    ("disease.min_susceptible_age_days", {
        "low": 0, "high": 180, "dist": "uniform",
        "desc": "Days post-settlement before susceptible to infection", "confidence": "★☆☆",
    }),
    # ── NEW PARAMS (SA Round 3 — spawning overhaul) ──────────────────
    ("spawning.p_spontaneous_male", {
        "low": 0.005, "high": 0.025, "dist": "uniform",
        "desc": "Daily spontaneous male spawning prob", "confidence": "★★☆",
    }),
    ("spawning.peak_width_days", {
        "low": 30.0, "high": 90.0, "dist": "uniform",
        "desc": "Spawning season peak width (std dev days)", "confidence": "★★☆",
    }),
    ("spawning.readiness_induction_prob", {
        "low": 0.1, "high": 0.8, "dist": "uniform",
        "desc": "Social spawning readiness induction probability", "confidence": "★☆☆",
    }),
    ("spawning.female_max_bouts", {
        "low": 1, "high": 3, "dist": "discrete",
        "values": [1, 2, 3],
        "desc": "Maximum spawning bouts per female per season", "confidence": "★★☆",
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
    
    # Enforce constraint: n_resistance + n_tolerance + n_recovery = 51
    if "genetics" in overrides:
        g = overrides["genetics"]
        n_r = g.get("n_resistance", 17)
        n_t = g.get("n_tolerance", 17)
        # Derive n_recovery = 51 - n_r - n_t (clamped to ≥1)
        n_c = max(1, 51 - n_r - n_t)
        # If partition overflows, adjust tolerance down
        if n_r + n_t + n_c != 51:
            n_t = max(1, 51 - n_r - n_c)
        g["n_recovery"] = n_c
        g["n_resistance"] = n_r
        g["n_tolerance"] = n_t
        # Sync tau_max to disease section (disease reads it from there too)
        if "tau_max" in g:
            if "disease" not in overrides:
                overrides["disease"] = {}
            overrides["disease"]["tau_max"] = g["tau_max"]
    
    # Integer coercion for fields that must be int
    int_fields = {
        ("disease", "immunosuppression_duration"),
        ("disease", "min_susceptible_age_days"),
        ("spawning", "female_max_bouts"),
    }
    for (sec, fld) in int_fields:
        if sec in overrides and fld in overrides[sec]:
            overrides[sec][fld] = int(round(overrides[sec][fld]))
    
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
