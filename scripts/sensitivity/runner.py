"""
Core sensitivity analysis runner.

Runs a single simulation with overridden parameters and extracts metrics.
Designed to be called from multiprocessing pool.
"""

import time
import traceback
import numpy as np

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from sswd_evoepi.config import default_config, validate_config
from sswd_evoepi.model import run_coupled_simulation
from scripts.sensitivity.param_spec import (
    sample_to_config_overrides,
    apply_overrides_to_config,
    get_param_names,
)


# ─── Simulation configuration ────────────────────────────────────────
SIM_N = 200          # agents
SIM_K = 200          # carrying capacity
SIM_YEARS = 20       # total years
SIM_DISEASE_YEAR = 3 # disease introduction year
SIM_T = 15.0         # temperature (°C)
SIM_SAL = 30.0       # salinity (psu)
SIM_AREA = 10000.0   # habitat area (m²)
SIM_INFECTED = 5     # initial infected


# ─── Metric extraction ────────────────────────────────────────────────

def extract_metrics(result):
    """Extract sensitivity analysis metrics from CoupledSimResult.
    
    Returns dict of scalar metric values.
    """
    metrics = {}
    
    # Population dynamics
    initial = result.initial_pop
    final = result.final_pop
    minimum = result.min_pop
    
    metrics["pop_crash_pct"] = (
        (initial - minimum) / initial * 100.0 if initial > 0 else 0.0
    )
    metrics["final_pop_frac"] = final / initial if initial > 0 else 0.0
    metrics["recovery"] = 1.0 if final > initial * 0.5 else 0.0
    metrics["extinction"] = 1.0 if final == 0 else 0.0
    metrics["peak_mortality"] = result.peak_mortality_fraction
    metrics["time_to_nadir"] = float(result.min_pop_year)
    metrics["total_disease_deaths"] = float(result.total_disease_deaths)
    
    # Resistance evolution
    mr = result.yearly_mean_resistance
    if mr is not None and len(mr) > SIM_DISEASE_YEAR:
        pre_resist = mr[SIM_DISEASE_YEAR]  # just before epidemic wave
        post_resist = mr[-1]
        metrics["resistance_shift"] = post_resist - pre_resist
    else:
        metrics["resistance_shift"] = 0.0
    
    # Additive variance retention
    va = result.yearly_va
    if va is not None and len(va) > SIM_DISEASE_YEAR and va[SIM_DISEASE_YEAR] > 0:
        metrics["va_retention"] = va[-1] / va[SIM_DISEASE_YEAR]
    else:
        metrics["va_retention"] = 0.0 if (va is None or len(va) == 0) else 1.0
    
    # EF1A dynamics
    ef = result.yearly_ef1a_freq
    if ef is not None and len(ef) > 1:
        metrics["ef1a_shift"] = abs(float(ef[-1]) - float(ef[0]))
    else:
        metrics["ef1a_shift"] = 0.0
    
    return metrics


METRIC_NAMES = [
    "pop_crash_pct",
    "final_pop_frac",
    "recovery",
    "extinction",
    "peak_mortality",
    "time_to_nadir",
    "total_disease_deaths",
    "resistance_shift",
    "va_retention",
    "ef1a_shift",
]


def run_single(args):
    """Run a single sensitivity analysis simulation.
    
    Args:
        args: tuple (run_index, sample_row, param_names, base_seed)
    
    Returns:
        dict with run_index, metrics, runtime, error
    """
    run_index, sample_row, param_names, base_seed = args
    
    t0 = time.time()
    
    try:
        # Build config with overrides
        config = default_config()
        overrides = sample_to_config_overrides(sample_row, param_names)
        apply_overrides_to_config(config, overrides)
        
        # Validate
        validate_config(config)
        
        # Unique seed per run
        seed = base_seed + run_index
        
        # Run simulation
        result = run_coupled_simulation(
            n_individuals=SIM_N,
            carrying_capacity=SIM_K,
            habitat_area=SIM_AREA,
            T_celsius=SIM_T,
            salinity=SIM_SAL,
            n_years=SIM_YEARS,
            disease_year=SIM_DISEASE_YEAR,
            initial_infected=SIM_INFECTED,
            seed=seed,
            config=config,
        )
        
        metrics = extract_metrics(result)
        elapsed = time.time() - t0
        
        return {
            "run_index": run_index,
            "metrics": metrics,
            "runtime": elapsed,
            "error": None,
        }
    
    except Exception as e:
        elapsed = time.time() - t0
        return {
            "run_index": run_index,
            "metrics": {m: np.nan for m in METRIC_NAMES},
            "runtime": elapsed,
            "error": f"{type(e).__name__}: {str(e)[:200]}",
        }
