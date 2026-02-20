"""
Spatial sensitivity analysis runner.

Runs a 5-node spatial simulation with overridden parameters and extracts metrics.
Designed to be called from multiprocessing pool.
"""

import time
import numpy as np
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from sswd_evoepi.config import default_config, validate_config
from sswd_evoepi.model import run_spatial_simulation
from sswd_evoepi.spatial import NodeDefinition, SpatialNode, MetapopulationNetwork
from scripts.sensitivity.param_spec import (
    sample_to_config_overrides,
    apply_overrides_to_config,
    get_param_names,
)


# ─── 5-node prototype network ────────────────────────────────────────
# Realistic K values from the prototype runs

def build_network(rng_seed=42):
    """Build a 3-node NE Pacific network with realistic K.
    
    3 nodes capturing:
    - Sitka (cold, high-latitude Alaska)
    - Howe Sound (fjord refugium, BC) 
    - Monterey (warm, southern range)
    
    K=5000 per site based on field densities of 0.1-0.2 ind/m²
    at healthy sites (Willem, pers. comm.), with ~33,000 m² habitat patches.
    """
    node_defs = [
        NodeDefinition(
            node_id=0, name="Sitka", lat=57.05, lon=-135.33,
            subregion="AK-SE", habitat_area=33333.0, carrying_capacity=5000,
            mean_sst=8.5, sst_amplitude=4.0, salinity=30.0,
            flushing_rate=0.3,
        ),
        NodeDefinition(
            node_id=1, name="Howe_Sound", lat=49.55, lon=-123.23,
            subregion="BC", habitat_area=33333.0, carrying_capacity=5000,
            is_fjord=True, sill_depth=50.0, flushing_rate=0.15,
            mean_sst=10.0, sst_amplitude=5.0, salinity=28.0,
        ),
        NodeDefinition(
            node_id=2, name="Monterey", lat=36.60, lon=-121.90,
            subregion="CA", habitat_area=33333.0, carrying_capacity=5000,
            mean_sst=13.0, sst_amplitude=3.0, salinity=33.5,
            flushing_rate=0.6,
        ),
    ]
    
    # Haversine distances (km)
    n = len(node_defs)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                dlat = np.radians(node_defs[i].lat - node_defs[j].lat)
                dlon = np.radians(node_defs[i].lon - node_defs[j].lon)
                a = (np.sin(dlat/2)**2 + 
                     np.cos(np.radians(node_defs[i].lat)) * 
                     np.cos(np.radians(node_defs[j].lat)) * 
                     np.sin(dlon/2)**2)
                D[i,j] = 6371 * 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
    
    # Wrap in SpatialNodes
    nodes = []
    for nd in node_defs:
        sn = SpatialNode(definition=nd)
        nodes.append(sn)
    
    rng = np.random.default_rng(rng_seed)
    network = MetapopulationNetwork(
        nodes=nodes,
        C=np.zeros((n, n)),
        D=np.zeros((n, n)),
        distances=D,
        rng=rng,
    )
    
    return network


# Keep old name for backward compatibility
build_5node_network = build_network


# ─── Simulation config ────────────────────────────────────────────────
SIM_YEARS = 20
SIM_DISEASE_YEAR = 3
SIM_INFECTED_PER_NODE = 3


# ─── Metric extraction ────────────────────────────────────────────────

METRIC_NAMES = [
    "pop_crash_pct",
    "final_pop_frac",
    "recovery",
    "extinction",
    "peak_mortality",
    "time_to_nadir",
    "total_disease_deaths",
    "resistance_shift_mean",
    "resistance_shift_max",
    "va_retention_mean",
    "ef1a_shift_mean",
    "n_extinct_nodes",
    "north_south_mortality_gradient",
    "fjord_protection_effect",
    # SA Round 2: pathogen evolution + cause-of-death metrics
    "mean_final_virulence",
    "virulence_shift",
    "disease_death_fraction",
    # SA Round 3: spawning + composite metrics
    "spawning_participation",
    "mean_recruitment_rate",
    "evolutionary_rescue_index",
]


def extract_spatial_metrics(result):
    """Extract sensitivity analysis metrics from SpatialSimResult."""
    metrics = {}
    n_nodes = result.n_nodes
    n_years = result.n_years
    
    # Aggregate population dynamics
    initial = result.initial_total_pop
    final = result.final_total_pop
    
    # Min total pop across years
    total_pops = result.yearly_total_pop if result.yearly_total_pop is not None else np.array([initial])
    min_pop = int(np.min(total_pops))
    min_year = int(np.argmin(total_pops))
    
    metrics["pop_crash_pct"] = (
        (initial - min_pop) / initial * 100.0 if initial > 0 else 0.0
    )
    metrics["final_pop_frac"] = final / initial if initial > 0 else 0.0
    metrics["recovery"] = 1.0 if final > initial * 0.5 else 0.0
    metrics["extinction"] = 1.0 if final == 0 else 0.0
    metrics["time_to_nadir"] = float(min_year)
    
    # Per-node disease deaths
    dis_deaths = result.yearly_disease_deaths  # (n_nodes, n_years)
    if dis_deaths is not None:
        total_dd = int(np.sum(dis_deaths))
        metrics["total_disease_deaths"] = float(total_dd)
        
        # Peak mortality fraction (any node, any year)
        pops = result.yearly_pop  # (n_nodes, n_years)
        if pops is not None:
            peak_mort = 0.0
            for ni in range(n_nodes):
                for yr in range(1, n_years):
                    pop_start = pops[ni, yr-1] if pops[ni, yr-1] > 0 else 1
                    mort_frac = dis_deaths[ni, yr] / pop_start
                    peak_mort = max(peak_mort, mort_frac)
            metrics["peak_mortality"] = peak_mort
        else:
            metrics["peak_mortality"] = 0.0
    else:
        metrics["total_disease_deaths"] = 0.0
        metrics["peak_mortality"] = 0.0
    
    # Resistance evolution
    mr = result.yearly_mean_resistance  # (n_nodes, n_years)
    if mr is not None and n_years > SIM_DISEASE_YEAR:
        pre = mr[:, SIM_DISEASE_YEAR]  # just before epidemic wave
        post = mr[:, -1]
        shifts = post - pre
        metrics["resistance_shift_mean"] = float(np.mean(shifts))
        metrics["resistance_shift_max"] = float(np.max(shifts))
    else:
        metrics["resistance_shift_mean"] = 0.0
        metrics["resistance_shift_max"] = 0.0
    
    # Additive variance retention
    va = result.yearly_va  # (n_nodes, n_years)
    if va is not None and n_years > SIM_DISEASE_YEAR:
        pre_va = va[:, SIM_DISEASE_YEAR]
        post_va = va[:, -1]
        retentions = np.where(pre_va > 0, post_va / pre_va, 1.0)
        metrics["va_retention_mean"] = float(np.mean(retentions))
    else:
        metrics["va_retention_mean"] = 1.0
    
    # EF1A dynamics
    ef = result.yearly_ef1a_freq  # (n_nodes, n_years)
    if ef is not None and n_years > 1:
        shifts = np.abs(ef[:, -1] - ef[:, 0])
        metrics["ef1a_shift_mean"] = float(np.mean(shifts))
    else:
        metrics["ef1a_shift_mean"] = 0.0
    
    # Node-level outcomes
    pops = result.yearly_pop
    if pops is not None:
        final_pops = pops[:, -1]
        metrics["n_extinct_nodes"] = float(np.sum(final_pops == 0))
        
        # North-south gradient: compare mortality at Sitka (0) vs Monterey (last)
        if n_nodes >= 2 and dis_deaths is not None:
            K = result.node_K if result.node_K is not None else np.ones(n_nodes)
            sitka_mort = np.sum(dis_deaths[0]) / max(K[0], 1)
            monterey_mort = np.sum(dis_deaths[-1]) / max(K[-1], 1)
            metrics["north_south_mortality_gradient"] = monterey_mort - sitka_mort
        else:
            metrics["north_south_mortality_gradient"] = 0.0
        
        # Fjord protection: Howe Sound (1, fjord) vs average of non-fjord nodes
        if n_nodes >= 3 and result.node_K is not None:
            hs_survival = final_pops[1] / max(result.node_K[1], 1)
            # Average of non-fjord nodes (Sitka=0, Monterey=2)
            non_fjord_survival = np.mean([
                final_pops[i] / max(result.node_K[i], 1)
                for i in range(n_nodes) if i != 1
            ])
            metrics["fjord_protection_effect"] = hs_survival - non_fjord_survival
        else:
            metrics["fjord_protection_effect"] = 0.0
    else:
        metrics["n_extinct_nodes"] = 0.0
        metrics["north_south_mortality_gradient"] = 0.0
        metrics["fjord_protection_effect"] = 0.0
    
    # ── SA Round 2: Pathogen evolution metrics ──────────────────────
    mv = getattr(result, 'yearly_mean_virulence', None)
    if mv is not None and mv.shape[-1] > 0:
        # mean_final_virulence: average across nodes at final year
        if mv.ndim == 2:  # (n_nodes, n_years)
            final_v = mv[:, -1]
            # Filter nodes that still have infected agents (v > 0)
            valid = final_v[final_v > 0]
            metrics["mean_final_virulence"] = float(np.mean(valid)) if len(valid) > 0 else 0.0
            # virulence_shift: change from initial
            init_v = mv[:, SIM_DISEASE_YEAR] if mv.shape[1] > SIM_DISEASE_YEAR else mv[:, 0]
            valid_init = init_v[init_v > 0]
            v_init_mean = float(np.mean(valid_init)) if len(valid_init) > 0 else 0.5
            v_final_mean = metrics["mean_final_virulence"]
            metrics["virulence_shift"] = v_final_mean - v_init_mean
        elif mv.ndim == 1:  # (n_years,)
            metrics["mean_final_virulence"] = float(mv[-1]) if mv[-1] > 0 else 0.0
            v_init_val = float(mv[SIM_DISEASE_YEAR]) if len(mv) > SIM_DISEASE_YEAR and mv[SIM_DISEASE_YEAR] > 0 else 0.5
            metrics["virulence_shift"] = metrics["mean_final_virulence"] - v_init_val
        else:
            metrics["mean_final_virulence"] = 0.0
            metrics["virulence_shift"] = 0.0
    else:
        metrics["mean_final_virulence"] = 0.0
        metrics["virulence_shift"] = 0.0
    
    # disease_death_fraction: fraction of all deaths caused by disease
    total_dd = metrics.get("total_disease_deaths", 0.0)
    # Fix: spatial result has yearly_* arrays, not total_* scalars
    ynd = getattr(result, 'yearly_natural_deaths', None)
    ysd = getattr(result, 'yearly_senescence_deaths', None)
    total_natural = float(np.sum(ynd)) if ynd is not None else 0.0
    total_senescence = float(np.sum(ysd)) if ysd is not None else 0.0
    all_deaths = total_dd + total_natural + total_senescence
    metrics["disease_death_fraction"] = total_dd / all_deaths if all_deaths > 0 else 0.0
    
    # ── SA Round 3: spawning + composite metrics ───────────────────
    recruits = result.yearly_recruits  # (n_nodes, n_years)
    node_K = result.node_K if result.node_K is not None else np.ones(n_nodes) * 5000
    
    # mean_recruitment_rate: mean annual recruits/K over pre-disease years
    if recruits is not None and SIM_DISEASE_YEAR > 0:
        if recruits.ndim == 2:  # (n_nodes, n_years)
            pre_disease_recruits = recruits[:, :SIM_DISEASE_YEAR]
            annual_rate = np.sum(pre_disease_recruits, axis=0) / np.sum(node_K)
        else:  # (n_years,)
            pre_disease_recruits = recruits[:SIM_DISEASE_YEAR]
            annual_rate = pre_disease_recruits / np.sum(node_K)
        metrics["mean_recruitment_rate"] = float(np.mean(annual_rate)) if len(annual_rate) > 0 else 0.0
    else:
        metrics["mean_recruitment_rate"] = 0.0
    
    # spawning_participation: mean annual recruits/K across ALL years (not just pre-disease)
    # Original binary metric was degenerate (always 1.0 over 3 pre-disease years).
    # This continuous version captures how recruitment rate changes through the epidemic.
    if recruits is not None:
        if recruits.ndim == 2:
            annual_total = np.sum(recruits, axis=0)  # sum across nodes per year
        else:
            annual_total = recruits
        metrics["spawning_participation"] = float(np.mean(annual_total / np.sum(node_K)))
    else:
        metrics["spawning_participation"] = 0.0
    
    # evolutionary_rescue_index: final_pop_frac × resistance_shift_mean
    metrics["evolutionary_rescue_index"] = (
        metrics.get("final_pop_frac", 0.0) * metrics.get("resistance_shift_mean", 0.0)
    )
    
    return metrics


def run_single_spatial(args):
    """Run a single spatial sensitivity analysis simulation.
    
    Args:
        args: tuple (run_index, sample_row, param_names, base_seed)
    
    Returns:
        dict with run_index, metrics, runtime, error
    """
    run_index, sample_row, param_names, base_seed = args
    run_index = int(run_index)
    base_seed = int(base_seed)
    
    t0 = time.time()
    
    try:
        # Build config with overrides
        config = default_config()
        overrides = sample_to_config_overrides(sample_row, param_names)
        apply_overrides_to_config(config, overrides)
        
        # Enable movement with 1 substep/day for SA
        # Full 24 substeps too expensive (20× overhead); 1 substep captures
        # spatial mixing and aggregation effects without hourly granularity
        config.movement.enabled = True
        config.movement.spatial_transmission = True
        config.movement.substeps_per_day = 1
        
        # Enable pathogen evolution for SA Round 2
        config.pathogen_evolution.enabled = True
        
        # Enable Beta genetics init
        config.genetics.q_init_mode = "beta"
        
        # Validate
        validate_config(config)
        
        # Unique seed per run
        seed = base_seed + run_index
        
        # Build fresh network for each run
        network = build_5node_network(rng_seed=seed + 1000000)
        
        # Run spatial simulation
        result = run_spatial_simulation(
            network=network,
            n_years=SIM_YEARS,
            disease_year=SIM_DISEASE_YEAR,
            initial_infected_per_node=SIM_INFECTED_PER_NODE,
            seed=seed,
            config=config,
        )
        
        metrics = extract_spatial_metrics(result)
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


if __name__ == "__main__":
    # Quick benchmark
    import time
    from scripts.sensitivity.param_spec import get_param_names, get_salib_problem
    from SALib.sample import morris as morris_sample
    
    names = get_param_names()
    problem = get_salib_problem()
    X = morris_sample.sample(problem, N=1, seed=12345)
    
    print(f"Benchmarking 3-node spatial sim (K=5000 each, 15K total)...")
    
    # Default params
    r = run_single_spatial((0, X[0], names, 99999))
    print(f"Run 1: {r['runtime']:.1f}s, error={r['error']}")
    if r['error'] is None:
        for m, v in r['metrics'].items():
            print(f"  {m}: {v:.4f}" if isinstance(v, float) else f"  {m}: {v}")
    
    # Another sample point
    r2 = run_single_spatial((1, X[1] if len(X) > 1 else X[0], names, 99999))
    print(f"\nRun 2: {r2['runtime']:.1f}s, error={r2['error']}")
