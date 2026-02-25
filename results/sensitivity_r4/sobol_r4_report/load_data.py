"""Shared data loader for Sobol R4 report figures."""
import json
import numpy as np
import os

DATA_DIR = os.path.join(os.path.dirname(__file__), '..')
FIG_DIR = os.path.join(os.path.dirname(__file__), 'figures')

def load_sobol():
    with open(os.path.join(DATA_DIR, 'sobol_results.json')) as f:
        d = json.load(f)
    return d

def get_param_names(d):
    return d['param_names']

def get_metric_names(d):
    return d['metric_names']

def get_short_names(param_names):
    """Strip module prefix for display."""
    return [p.split('.')[-1] for p in param_names]

def get_param_groups(param_names):
    """Return group name for each parameter."""
    return [p.split('.')[0] for p in param_names]

# Consistent color scheme for parameter groups
GROUP_COLORS = {
    'disease': '#e74c3c',
    'genetics': '#3498db',
    'population': '#2ecc71',
    'spawning': '#f39c12',
    'spatial': '#9b59b6',
    'pathogen_evolution': '#e67e22',
}

GROUP_ORDER = ['disease', 'genetics', 'population', 'spawning', 'spatial', 'pathogen_evolution']

GROUP_LABELS = {
    'disease': 'Disease',
    'genetics': 'Genetics',
    'population': 'Population',
    'spawning': 'Spawning',
    'spatial': 'Spatial',
    'pathogen_evolution': 'Pathogen Evol.',
}

def setup_style():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    try:
        plt.style.use('seaborn-v0_8-whitegrid')
    except:
        plt.style.use('seaborn-whitegrid')
    plt.rcParams.update({
        'font.size': 10,
        'axes.titlesize': 12,
        'axes.labelsize': 11,
        'xtick.labelsize': 9,
        'ytick.labelsize': 9,
        'legend.fontsize': 9,
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.1,
    })
    return plt

# Parameter ranges (from the Sobol problem definition, R4)
PARAM_RANGES = {
    'disease.a_exposure': (0.5, 5.0),
    'disease.K_half': (50, 5000),
    'disease.sigma_1_eff': (1e4, 1e7),
    'disease.sigma_2_eff': (1e5, 1e8),
    'disease.sigma_D': (0.01, 0.5),
    'disease.rho_rec': (0.001, 0.05),
    'disease.mu_EI1_ref': (0.05, 0.5),
    'disease.mu_I2D_ref': (0.005, 0.1),
    'disease.P_env_max': (1e3, 1e7),
    'disease.T_ref': (10.0, 18.0),
    'population.F0': (1e4, 1e6),
    'population.gamma_fert': (0.5, 2.0),
    'population.settler_survival': (0.001, 0.05),
    'population.alpha_srs': (0.0001, 0.01),
    'population.senescence_age': (15, 35),
    'population.k_growth': (0.1, 0.5),
    'population.L_min_repro': (10, 30),
    'genetics.n_resistance': (2, 30),
    'genetics.n_tolerance': (2, 30),
    'genetics.target_mean_t': (0.1, 0.9),
    'genetics.target_mean_c': (0.1, 0.9),
    'genetics.tau_max': (0.1, 0.9),
    'spawning.p_spontaneous_female': (0.001, 0.1),
    'spawning.induction_female_to_male': (0.01, 0.5),
    'disease.susceptibility_multiplier': (0.5, 2.0),
    'disease.T_vbnc': (2.0, 15.0),
    'disease.s_min': (0.01, 0.5),
    'disease.mu_I1I2_ref': (0.01, 0.2),
    'disease.immunosuppression_duration': (30, 365),
    'spawning.induction_male_to_female': (0.01, 0.5),
    'spatial.D_L': (1.0, 100.0),
    'genetics.target_mean_r': (0.1, 0.9),
    'genetics.q_init_beta_a': (0.5, 5.0),
    'genetics.q_init_beta_b': (0.5, 5.0),
    'spatial.alpha_self_fjord': (0.5, 0.99),
    'spatial.alpha_self_open': (0.1, 0.8),
    'pathogen_evolution.alpha_kill': (0.0, 1.0),
    'pathogen_evolution.alpha_shed': (0.0, 1.0),
    'pathogen_evolution.alpha_prog': (0.0, 1.0),
    'pathogen_evolution.gamma_early': (0.0, 0.5),
    'pathogen_evolution.sigma_v_mutation': (0.001, 0.1),
    'pathogen_evolution.v_init': (0.1, 0.9),
    'disease.min_susceptible_age_days': (30, 365),
    'spawning.p_spontaneous_male': (0.001, 0.1),
    'spawning.peak_width_days': (10, 60),
    'spawning.readiness_induction_prob': (0.01, 0.5),
    'spawning.female_max_bouts': (1, 5),
}

def get_range_widths(param_names):
    """Return normalized range widths (range / midpoint) for each parameter."""
    widths = []
    for p in param_names:
        if p in PARAM_RANGES:
            lo, hi = PARAM_RANGES[p]
            mid = (lo + hi) / 2.0
            if mid > 0:
                widths.append((hi - lo) / mid)
            else:
                widths.append(hi - lo)
        else:
            widths.append(1.0)
    return np.array(widths)

# Morris R4 ranking (from the analysis MD)
MORRIS_RANKING = {
    'disease.rho_rec': 1,
    'population.k_growth': 2,
    'disease.K_half': 3,
    'disease.P_env_max': 4,
    'genetics.n_resistance': 5,
    'population.settler_survival': 6,
    'disease.sigma_2_eff': 7,
    'disease.mu_I2D_ref': 8,
    'spawning.peak_width_days': 9,
    'genetics.target_mean_c': 10,
    'disease.a_exposure': 11,
    'disease.T_vbnc': 12,
    'genetics.tau_max': 13,
    'pathogen_evolution.sigma_v_mutation': 14,
    'pathogen_evolution.alpha_kill': 15,
    'disease.sigma_1_eff': 16,
    'genetics.target_mean_r': 17,
    'disease.T_ref': 18,
    'disease.min_susceptible_age_days': 19,
    'disease.sigma_D': 20,
    'spawning.female_max_bouts': 21,
    'spawning.readiness_induction_prob': 22,
    'genetics.target_mean_t': 23,
    'genetics.n_tolerance': 24,
    'spatial.alpha_self_open': 25,
    'spatial.D_L': 26,
    'spawning.induction_male_to_female': 27,
    'disease.s_min': 28,
    'pathogen_evolution.v_init': 29,
    'spawning.p_spontaneous_male': 30,
    'disease.mu_I1I2_ref': 31,
    'genetics.q_init_beta_a': 32,
    'spatial.alpha_self_fjord': 33,
    'population.senescence_age': 34,
    'pathogen_evolution.gamma_early': 35,
    'population.alpha_srs': 36,
    'pathogen_evolution.alpha_prog': 37,
    'disease.mu_EI1_ref': 38,
    'population.L_min_repro': 39,
    'pathogen_evolution.alpha_shed': 40,
    'spawning.induction_female_to_male': 41,
    'disease.immunosuppression_duration': 42,
    'population.gamma_fert': 43,
    'disease.susceptibility_multiplier': 44,
    'spawning.p_spontaneous_female': 45,
    'genetics.q_init_beta_b': 46,
    'population.F0': 47,
}
