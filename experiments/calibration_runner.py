#!/usr/bin/env python3
"""Agentic calibration runner for SSWD-EvoEpi model.

Runs the full 896-node network with satellite SST, measures regional recovery
fractions, and scores against Willem's expert-informed targets.

Usage:
    python3 calibration_runner.py --config params.json --seed 42 --output results/calibration/round_00/
    
For multi-seed runs:
    python3 calibration_runner.py --config params.json --seeds 42,123,999 --output results/calibration/round_00/
"""

import argparse
import json
import sys
import os
import time
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional

# Add project root to path
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from sswd_evoepi.config import SimulationConfig, default_config, validate_config
from sswd_evoepi.spatial import NodeDefinition, build_network
from sswd_evoepi.model import run_spatial_simulation, SpatialSimResult
from sswd_evoepi.metrics import RECOVERY_TARGETS, ARRIVAL_TARGETS, rmsle, rmse_log, score_regions

# ── Paths ─────────────────────────────────────────────────────────────
SITES_PATH = 'data/nodes/all_sites.json'
DISTANCE_MATRIX_PATH = 'results/overwater/distance_matrix.npz'


def _estimate_mean_sst(lat: float) -> float:
    """Latitude → rough mean SST (°C) for initialization."""
    return max(4.0, min(20.0, 25.0 - 0.35 * lat))

def _estimate_sst_amplitude(lat: float) -> float:
    """Latitude → SST annual amplitude (°C)."""
    return max(2.0, min(6.0, 0.08 * lat - 1.0))


def load_sites() -> List[dict]:
    """Load site list."""
    with open(PROJECT_ROOT / SITES_PATH) as f:
        return json.load(f)


def build_node_defs(sites: List[dict], K: int = 5000, K_cv: float = 0.0,
                    seed: int = 42, n_connectivity: float = 0.3,
                    phi_open: float = 0.8, phi_fjord: float = 0.03,
                    K_ref: int = 5000) -> List[NodeDefinition]:
    """Build NodeDefinition list from site JSON.

    Args:
        sites: List of site dicts with latitude, longitude, region, name.
        K: Base carrying capacity per node.
        K_cv: Coefficient of variation for per-node K (lognormal).
            0 = uniform K (all nodes get K). >0 = lognormal variability
            with E[K] = K exactly (bias-corrected parameterization).
        seed: RNG seed for reproducible K draws.
    """
    n_sites = len(sites)

    # Generate per-node K values
    if K_cv > 0:
        rng = np.random.default_rng(seed + 9999)  # offset to avoid collision with sim RNG
        sigma = np.sqrt(np.log(1 + K_cv**2))
        mu = np.log(K) - sigma**2 / 2
        K_values = rng.lognormal(mu, sigma, size=n_sites)
        K_values = np.maximum(K_values, 100).astype(int)  # floor at 100
    else:
        K_values = np.full(n_sites, K, dtype=int)

    # Load per-site enclosedness metrics (raw — φ computed at runtime)
    encl_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'nodes', 'site_enclosedness.json')
    encl_by_name = {}
    if os.path.exists(encl_path):
        with open(encl_path) as ef:
            encl_data = json.load(ef)
        for entry in encl_data:
            encl_by_name[entry['name']] = entry
        print(f"Loaded enclosedness data for {len(encl_by_name)} sites")
    else:
        print("WARNING: site_enclosedness.json not found — using uniform φ=0.5")

    # Connectivity parameters (phi_open/phi_fjord now passed as arguments)
    n_conn = n_connectivity  # nonlinearity exponent (1=linear, >1=sharper contrast)

    node_defs = []
    for i, site in enumerate(sites):
        lat = float(site['latitude'])
        lon = float(site['longitude'])
        region = site.get('region', '?')
        site_name = site['name']

        # Compute flushing rate from enclosedness metrics at runtime
        encl = encl_by_name.get(site_name)
        if encl:
            enclosedness = encl['enclosedness_combined']
            depth_norm = encl.get('fjord_depth_norm', 0.0)
            # Effective enclosedness: fjord depth amplifies (mouth=50%, head=100%)
            effective = enclosedness * (0.5 + 0.5 * depth_norm)
            effective = min(effective, 1.0)
            # Apply connectivity nonlinearity
            effective_nl = effective ** n_conn
            phi = phi_open * (1.0 - effective_nl) + phi_fjord * effective_nl
            is_fjord = enclosedness > 0.6
        else:
            is_fjord = region in ('BC-N', 'AK-PWS') and lat > 50.0
            phi = 0.03 if is_fjord else 0.5
            effective_nl = 1.0 if is_fjord else 0.0

        # Density-invariant K scaling: scale habitat area so that population
        # density (individuals/m²) stays the same as at K_ref.
        _area_scale = int(K_values[i]) / K_ref if K_ref > 0 else 1.0
        _habitat_area = 25000.0 * _area_scale

        nd = NodeDefinition(
            node_id=i,
            name=site_name,
            lat=lat,
            lon=lon,
            subregion=region,
            habitat_area=_habitat_area,
            carrying_capacity=int(K_values[i]),
            is_fjord=is_fjord,
            sill_depth=30.0 if is_fjord else np.inf,
            flushing_rate=phi,
            mean_sst=_estimate_mean_sst(lat),
            sst_amplitude=_estimate_sst_amplitude(lat),
            sst_trend=0.02,
            salinity=22.0 if is_fjord else 32.0,
            fjord_depth_norm=depth_norm if encl else 0.0,
            depth_range=(5.0, 60.0),
        )
        # Store effective enclosedness for continuous alpha_self
        nd._effective_enclosedness = effective_nl if encl else 0.0
        node_defs.append(nd)
    return node_defs


def build_full_network(K: int = 5000, seed: int = 42, D_L: float = 400.0,
                       D_P: float = 15.0, D_P_max_range: float = None,
                       K_cv: float = 0.0, n_connectivity: float = 0.3,
                       alpha_self_open: float = 0.05,
                       alpha_self_fjord: float = 0.50,
                       r_total: float = 0.02,
                       phi_open: float = 0.8,
                       phi_fjord: float = 0.03,
                       K_ref: int = 5000):
    """Build 896-node network with overwater distances.

    alpha_self is computed continuously from enclosedness:
        alpha_self[i] = alpha_self_open + (alpha_self_fjord - alpha_self_open) * effective_enclosedness[i]
    This means enclosed sites retain more larvae locally (higher self-recruitment).

    When K != K_ref, habitat_area is scaled by K/K_ref so that population
    density (ind/m²) is preserved (density-invariant K scaling).
    """
    sites = load_sites()
    node_defs = build_node_defs(sites, K=K, K_cv=K_cv, seed=seed,
                                n_connectivity=n_connectivity,
                                phi_open=phi_open, phi_fjord=phi_fjord,
                                K_ref=K_ref)
    npz_path = str(PROJECT_ROOT / DISTANCE_MATRIX_PATH)

    # Compute continuous alpha_self from enclosedness
    alpha_self = np.array([
        alpha_self_open + (alpha_self_fjord - alpha_self_open) * getattr(nd, '_effective_enclosedness', 0.0)
        for nd in node_defs
    ], dtype=np.float64)

    network = build_network(
        node_defs,
        D_L=D_L,
        D_P=D_P,
        D_P_max_range=D_P_max_range,
        r_total=r_total,
        seed=seed,
        overwater_npz=npz_path,
        alpha_self=alpha_self,
    )
    return sites, network


def apply_param_overrides(config: SimulationConfig, overrides: dict):
    """Apply dotted parameter overrides to config."""
    for key, value in overrides.items():
        parts = key.split('.')
        obj = config
        for part in parts[:-1]:
            obj = getattr(obj, part)
        setattr(obj, parts[-1], value)


def compute_regional_recovery(result: SpatialSimResult, sites: List[dict]) -> Tuple[dict, dict]:
    """Compute recovery fraction per region.
    
    Recovery = final_population / peak_pre_crash_population.
    Pre-crash = max population in first 3 simulation years.
    """
    yearly_pop = result.yearly_pop  # shape: (n_nodes, n_years)
    n_years = yearly_pop.shape[1]
    
    # Group nodes by region
    region_nodes: Dict[str, List[int]] = {}
    for i, site in enumerate(sites):
        r = site.get('region', '?')
        region_nodes.setdefault(r, []).append(i)
    
    region_recovery = {}
    region_details = {}
    
    for region, idxs in sorted(region_nodes.items()):
        # Sum populations across nodes in region per year
        region_pop = [int(yearly_pop[idxs, y].sum()) for y in range(n_years)]
        
        # Peak = max in first 2 years (spinup baseline before disease bites)
        baseline_years = min(2, n_years)
        peak_pop = max(region_pop[:baseline_years]) if baseline_years > 0 else region_pop[0]
        final_pop = region_pop[-1]
        
        recovery = final_pop / peak_pop if peak_pop > 0 else 0.0
        
        region_recovery[region] = recovery
        # Per-year recruits and disease deaths for diagnostics
        yearly_recruits = result.yearly_recruits
        yearly_dd = result.yearly_disease_deaths
        region_recruits = [int(yearly_recruits[idxs, y].sum()) for y in range(n_years)] if yearly_recruits is not None else []
        region_dd = [int(yearly_dd[idxs, y].sum()) for y in range(n_years)] if yearly_dd is not None else []
        
        # Genetics: population-weighted mean traits per region per year
        def _weighted_mean_trait(trait_arr, pop_arr, idxs, n_years):
            """Compute population-weighted mean of a trait across region nodes."""
            if trait_arr is None:
                return []
            out = []
            for y in range(n_years):
                pops = pop_arr[idxs, y].astype(float)
                traits = trait_arr[idxs, y]
                total_pop = pops.sum()
                if total_pop > 0:
                    out.append(float((traits * pops).sum() / total_pop))
                else:
                    out.append(0.0)
            return out
        
        region_mean_r = _weighted_mean_trait(result.yearly_mean_resistance, yearly_pop, idxs, n_years)
        region_mean_t = _weighted_mean_trait(result.yearly_mean_tolerance, yearly_pop, idxs, n_years)
        region_mean_c = _weighted_mean_trait(result.yearly_mean_recovery, yearly_pop, idxs, n_years)
        
        # Additive variance (mean across nodes, not pop-weighted — V_A is per-node)
        def _mean_va(va_arr, idxs, n_years):
            if va_arr is None:
                return []
            return [float(va_arr[idxs, y].mean()) for y in range(n_years)]
        
        region_va_r = _mean_va(result.yearly_va, idxs, n_years)
        region_va_t = _mean_va(result.yearly_va_tolerance, idxs, n_years)
        region_va_c = _mean_va(result.yearly_va_recovery, idxs, n_years)
        
        region_details[region] = {
            'n_nodes': len(idxs),
            'peak_pop': int(peak_pop),
            'final_pop': int(final_pop),
            'recovery_frac': float(recovery),
            'crash_pct': float(100.0 * (1 - final_pop / peak_pop)) if peak_pop > 0 else 100.0,
            'yearly_totals': [int(p) for p in region_pop],
            'yearly_recruits': region_recruits,
            'yearly_disease_deaths': region_dd,
            'yearly_mean_resistance': region_mean_r,
            'yearly_mean_tolerance': region_mean_t,
            'yearly_mean_recovery': region_mean_c,
            'yearly_va_resistance': region_va_r,
            'yearly_va_tolerance': region_va_t,
            'yearly_va_recovery': region_va_c,
        }

        # Pathogen thermal adaptation: per-region mean final T_vbnc
        if result.T_vbnc_local is not None:
            region_details[region]['final_mean_T_vbnc'] = float(
                np.mean(result.T_vbnc_local[idxs])
            )

        # Community virulence: per-region mean final v_local
        if result.v_local is not None:
            region_details[region]['final_mean_v_local'] = float(
                np.mean(result.v_local[idxs])
            )
    
    return region_recovery, region_details


def score_against_targets(region_recovery: dict) -> dict:
    """Score using log-space RMSE (targets span 3 orders of magnitude).

    Delegates to sswd_evoepi.metrics.{rmsle, score_regions} and
    reformats the output into the JSON-compatible dict expected by
    downstream tools (format_results, viz, results.py).
    """
    per_region = score_regions(region_recovery)
    overall_rmse = rmsle(region_recovery)

    # Build per-region dict with extra fields for JSON compatibility
    scores = {}
    for region, info in per_region.items():
        scores[region] = {
            'target': info['target'],
            'actual': info['actual'],
            'target_pct': float(info['target'] * 100),
            'actual_pct': float(info['actual'] * 100),
            'log_error': info['log_error'],
            'log_sq_error': float(info['log_error'] ** 2),
            'within_2x': info['within_2x'],
            'within_5x': info['within_5x'],
        }

    return {
        'per_region': scores,
        'rmsle': float(overall_rmse),
        'rmse_log': float(overall_rmse),  # backward compat
        'within_2x': sum(1 for s in scores.values() if s['within_2x']),
        'within_5x': sum(1 for s in scores.values() if s['within_5x']),
        'n_targets': len(scores),
    }


def score_arrival_timing(result, sites: list, disease_year: int = 1) -> dict:
    """Score wavefront arrival timing against observed data.

    Args:
        result: SpatialSimResult with disease_arrival_day array.
        sites: List of site dicts with 'region' field.
        disease_year: Simulation year when disease was seeded.

    Returns:
        Dict with per-region timing, MAE, and individual scores.
    """
    if result.disease_arrival_day is None:
        return {'mae_months': float('inf'), 'per_region': {}, 'note': 'No wavefront data'}

    origin_day = disease_year * 365  # day disease was seeded

    # Build region → node indices
    region_nodes = {}
    for i, site in enumerate(sites):
        r = site.get('region', '?')
        region_nodes.setdefault(r, []).append(i)

    region_arrival = {}
    total_abs_error = 0.0
    n_targets = 0

    for region, target_months in ARRIVAL_TARGETS.items():
        idxs = region_nodes.get(region, [])
        if not idxs:
            continue

        # First arrival in region (earliest node)
        arrival_days = result.disease_arrival_day[idxs]
        reached = arrival_days[arrival_days >= 0]

        if len(reached) == 0:
            actual_months = float('inf')
        else:
            first_day = int(reached.min())
            actual_months = (first_day - origin_day) / 30.44  # avg days per month

        error = actual_months - target_months if actual_months != float('inf') else 60.0  # penalty

        region_arrival[region] = {
            'target_months': target_months,
            'actual_months': round(actual_months, 1) if actual_months != float('inf') else None,
            'error_months': round(error, 1),
            'reached': actual_months != float('inf'),
        }

        total_abs_error += abs(error)
        n_targets += 1

    mae = total_abs_error / n_targets if n_targets > 0 else float('inf')

    return {
        'mae_months': round(mae, 2),
        'n_targets': n_targets,
        'per_region': region_arrival,
    }


class EarlyStopChecker:
    """Year-level callback for early stopping and progress logging.

    Aborts clearly-failing calibration runs to save compute:
      - Year 3: total extinction or no epidemic (pop > 95% of initial)
      - Year 5: no regional gradient developing
    """

    def __init__(self, network, sites: List[dict], round_id: str, seed: int,
                 checkpoint_dir: Optional[Path] = None):
        self.network = network
        self.sites = sites
        self.round_id = round_id
        self.seed = seed
        self.early_stop_reason: Optional[str] = None
        self.checkpoint_dir = checkpoint_dir
        # Per-node peak population in first 2 years (for recovery calc)
        self._node_peak_pop: Optional[np.ndarray] = None
        self._initial_total_pop: Optional[int] = None
        # Region → node index mapping (built lazily)
        self._region_nodes: Optional[Dict[str, list]] = None

    def _build_region_map(self):
        if self._region_nodes is None:
            self._region_nodes = {}
            for i, site in enumerate(self.sites):
                r = site.get('region', '?')
                self._region_nodes.setdefault(r, []).append(i)

    def __call__(self, year: int, n_years: int):
        """Called after each simulated year.  Returns a stop reason string or None."""
        node_pops = np.array([n.n_alive for n in self.network.nodes])
        total_pop = int(node_pops.sum())
        alive_nodes = int(np.sum(node_pops > 0))

        # ── Track baselines ──────────────────────────────────────────
        if self._initial_total_pop is None:
            self._initial_total_pop = total_pop
            self._node_peak_pop = node_pops.copy()
        elif year <= 1:
            self._node_peak_pop = np.maximum(self._node_peak_pop, node_pops)

        # ── Progress logging ─────────────────────────────────────────
        print(f"  [{self.round_id}/seed_{self.seed}] "
              f"Year {year + 1}/{n_years}: "
              f"pop={total_pop:,}, alive_nodes={alive_nodes}",
              flush=True)

        # ── Yearly checkpoint save ───────────────────────────────────
        if self.checkpoint_dir is not None:
            self._build_region_map()
            region_snapshot = {}
            for region, idxs in self._region_nodes.items():
                peak = float(self._node_peak_pop[idxs].sum()) if self._node_peak_pop is not None else 0
                current = float(node_pops[idxs].sum())
                region_snapshot[region] = {
                    'pop': int(current),
                    'recovery': current / peak if peak > 0 else 0.0,
                }
            checkpoint = {
                'round': self.round_id,
                'seed': self.seed,
                'year': year + 1,
                'n_years': n_years,
                'total_pop': total_pop,
                'alive_nodes': alive_nodes,
                'initial_pop': self._initial_total_pop,
                'pop_frac': total_pop / self._initial_total_pop if self._initial_total_pop else 0,
                'regions': region_snapshot,
            }
            cp_file = self.checkpoint_dir / f'checkpoint_seed{self.seed}.json'
            try:
                with open(cp_file, 'w') as f:
                    json.dump(checkpoint, f, indent=2)
            except Exception:
                pass  # Don't crash on checkpoint write failure

        # ── Year 3 checks (0-indexed year == 2) ─────────────────────
        if year == 2:
            if total_pop == 0:
                self.early_stop_reason = "Total extinction by year 3"
                return self.early_stop_reason
            if self._initial_total_pop and total_pop > 0.95 * self._initial_total_pop:
                self.early_stop_reason = (
                    f"No epidemic by year 3 (pop {total_pop:,} > 95% "
                    f"of initial {self._initial_total_pop:,})"
                )
                return self.early_stop_reason

        # ── Year 5 checks (0-indexed year == 4) ─────────────────────
        if year == 4 and self._node_peak_pop is not None:
            self._build_region_map()
            recoveries = []
            for region, idxs in self._region_nodes.items():
                peak = float(self._node_peak_pop[idxs].sum())
                current = float(node_pops[idxs].sum())
                frac = current / peak if peak > 0 else 0.0
                recoveries.append(frac)

            if all(r < 0.01 for r in recoveries):
                self.early_stop_reason = "All regions < 1% recovery at year 5"
                return self.early_stop_reason
            if all(r > 0.50 for r in recoveries):
                self.early_stop_reason = "All regions > 50% recovery at year 5"
                return self.early_stop_reason

            spread = max(recoveries) - min(recoveries)
            if spread < 0.02:
                self.early_stop_reason = (
                    f"No regional variation at year 5 (spread={spread:.4f})"
                )
                return self.early_stop_reason

        return None  # continue simulation


def get_sentinel_eligible_nodes(sites: List[dict], site_filter: str = 'rocky',
                                 custom_node_ids: List[int] = None) -> List[int]:
    """Return list of node IDs eligible for sentinel agents.

    Args:
        sites: List of site dicts with 'habitat' field.
        site_filter: 'rocky' (rock/kelp/reef), 'all', or 'custom'.
        custom_node_ids: Node IDs when site_filter='custom'.

    Returns:
        List of node indices eligible for sentinels.
    """
    if site_filter == 'all':
        return list(range(len(sites)))
    elif site_filter == 'custom':
        return list(custom_node_ids or [])
    else:  # 'rocky' (default)
        eligible = []
        for i, site in enumerate(sites):
            habitat = site.get('habitat', '').lower()
            if any(x in habitat for x in ['rock', 'kelp', 'reef']):
                eligible.append(i)
        return eligible


def run_single(config: SimulationConfig, sites: List[dict], network, seed: int,
               disease_year: int = 3, n_years: int = 14,
               round_id: str = "XX", output_dir: Optional[Path] = None) -> dict:
    """Run one simulation and return structured results.

    Uses year-level early stopping to abort runs that are clearly failing,
    saving up to ~60% of compute on bad parameter sets.
    """
    config.simulation.seed = seed

    checker = EarlyStopChecker(network, sites, round_id, seed,
                               checkpoint_dir=output_dir)

    # Lightweight monthly recorder for wavefront visualization
    from sswd_evoepi.monthly_recorder import MonthlyRecorder
    monthly_rec = MonthlyRecorder(
        n_nodes=network.n_nodes,
        sst_start_year=getattr(config, '_sst_start_year', 2012),
    )

    # Compute sentinel eligible nodes if enabled
    _sentinel_node_ids = None
    if config.sentinel.enabled:
        _sentinel_node_ids = get_sentinel_eligible_nodes(
            sites,
            config.sentinel.site_filter,
            config.sentinel.custom_node_ids,
        )

    t0 = time.time()
    result = run_spatial_simulation(
        network=network,
        n_years=n_years,
        disease_year=disease_year,
        seed=seed,
        config=config,
        progress_callback=checker,
        monthly_recorder=monthly_rec,
        sentinel_node_ids=_sentinel_node_ids,
    )
    elapsed = time.time() - t0

    # ── Handle early stop ────────────────────────────────────────────
    if checker.early_stop_reason:
        print(f"  ⚡ EARLY STOP: {checker.early_stop_reason} ({elapsed:.0f}s)")
        return {
            'seed': seed,
            'wall_time_seconds': float(elapsed),
            'scoring': {
                'rmsle': float('inf'),
                'rmse_log': float('inf'),  # backward compat
                'within_2x': 0,
                'within_5x': 0,
                'n_targets': len(RECOVERY_TARGETS),
                'per_region': {},
            },
            'region_details': {},
            'region_recovery': {},
            'overall': {
                'pop_crash_pct': 0.0,
                'final_pop_frac': 0.0,
            },
            'early_stop': checker.early_stop_reason,
        }

    # ── Save monthly snapshot data for wavefront visualization ────
    if monthly_rec.n_frames > 0 and output_dir is not None:
        snap_path = Path(output_dir) / f'monthly_seed{seed}.npz'
        site_lats = np.array([s.get('latitude', s.get('lat', 0)) for s in sites])
        site_lons = np.array([s.get('longitude', s.get('lon', 0)) for s in sites])
        site_names = [s.get('name', f'node_{i}') for i, s in enumerate(sites)]
        monthly_rec.save(str(snap_path), site_lats, site_lons, site_names,
                         K=network.nodes[0].definition.carrying_capacity if network.nodes else 5000)
        print(f"  Monthly snapshots: {monthly_rec.n_frames} frames → {snap_path.name}")

    # ── Normal scoring ───────────────────────────────────────────────
    region_recovery, region_details = compute_regional_recovery(result, sites)
    scoring = score_against_targets(region_recovery)

    # Arrival timing (wavefront mode)
    arrival_scoring = score_arrival_timing(result, sites, disease_year=disease_year)

    return {
        'seed': seed,
        'wall_time_seconds': float(elapsed),
        'scoring': scoring,
        'arrival_timing': arrival_scoring,
        'region_details': region_details,
        'region_recovery': {k: float(v) for k, v in region_recovery.items()},
        'overall': {
            'pop_crash_pct': float(100.0 * (1.0 - result.final_total_pop / max(result.initial_total_pop, 1))),
            'final_pop_frac': float(result.final_total_pop / max(result.initial_total_pop, 1)),
        },
    }


def print_summary(results: List[dict], param_overrides: dict):
    """Pretty-print calibration results."""
    print(f"\n{'='*75}")
    print(f"CALIBRATION RESULTS ({len(results)} seed{'s' if len(results)>1 else ''})")
    print(f"{'='*75}")
    
    if param_overrides:
        print(f"\nParameter overrides:")
        for k, v in sorted(param_overrides.items()):
            print(f"  {k} = {v}")
    else:
        print(f"\nUsing default parameters (baseline)")
    
    # Filter out early-stopped results for scoring summary
    scored_results = [r for r in results if 'early_stop' not in r]
    early_stopped = [r for r in results if 'early_stop' in r]

    if early_stopped:
        print(f"\nEarly stopped: {len(early_stopped)}/{len(results)} seeds")
        for r in early_stopped:
            print(f"  seed {r['seed']}: {r['early_stop']}")

    # Average across scored seeds
    avg_scoring = {}
    for region in RECOVERY_TARGETS:
        actuals = [r['scoring']['per_region'][region]['actual']
                   for r in scored_results
                   if region in r['scoring'].get('per_region', {})]
        avg_scoring[region] = {
            'target': RECOVERY_TARGETS[region],
            'mean': float(np.mean(actuals)) if actuals else 0.0,
            'std': float(np.std(actuals)) if len(actuals) > 1 else 0.0,
        }
    
    avg_rmse = float(np.mean([r['scoring'].get('rmsle', r['scoring'].get('rmse_log', float('inf'))) for r in results]))
    avg_crash = float(np.mean([r['overall']['pop_crash_pct'] for r in results]))
    
    print(f"\nOverall crash: {avg_crash:.1f}%")
    print(f"RMSLE:   {avg_rmse:.3f}  (0 = perfect, <0.3 = within 2×)")
    
    print(f"\n{'Region':<10s} {'Target':>10s} {'Actual':>12s} {'LogErr':>10s} {'Grade'}")
    print(f"{'-'*50}")
    for region in RECOVERY_TARGETS:
        s = avg_scoring[region]
        log_err = np.log10(max(s['mean'], 1e-6)) - np.log10(max(s['target'], 1e-6))
        if abs(log_err) < 0.301:
            grade = '✓ within 2×'
        elif abs(log_err) < 0.699:
            grade = '~ within 5×'
        else:
            grade = '✗ off'
        
        mean_str = f"{s['mean']*100:.2f}%"
        if s['std'] > 0:
            mean_str += f" ±{s['std']*100:.2f}%"
        
        print(f"{region:<10s} {s['target']*100:>9.2f}% {mean_str:>12s} {log_err:>+9.2f}  {grade}")
    
    # Unconstrained regions
    print(f"\nUnconstrained regions (mean across seeds):")
    all_regions = set()
    for r in results:
        all_regions.update(r['region_recovery'].keys())
    for region in sorted(all_regions):
        if region not in RECOVERY_TARGETS:
            vals = [r['region_recovery'].get(region, 0) for r in results]
            mean_val = np.mean(vals)
            print(f"  {region:<10s}: {mean_val*100:.2f}%")
    
    times = [r['wall_time_seconds'] for r in results]
    print(f"\nWall time: {np.mean(times):.0f}s mean ({np.sum(times):.0f}s total)")


def main():
    parser = argparse.ArgumentParser(description='Agentic calibration runner')
    parser.add_argument('--config', type=str, required=True, help='Path to param overrides JSON')
    parser.add_argument('--seed', type=int, default=None, help='Single seed')
    parser.add_argument('--seeds', type=str, default=None, help='Comma-separated seeds')
    parser.add_argument('--output', type=str, required=True, help='Output directory')
    parser.add_argument('--K', type=int, default=5000, help='Carrying capacity per node')
    parser.add_argument('--K-ref', type=int, default=5000,
                        help='Reference K for density-invariant scaling. When K != K_ref, '
                             'habitat area and shedding are scaled so population density '
                             'and disease dynamics match K_ref behavior. Default 5000.')
    parser.add_argument('--K-cv', type=float, default=0.0, help='CV for per-node K variability (lognormal, 0=uniform)')
    parser.add_argument('--years', type=int, default=13, help='Simulation years (1yr spinup + 12yr 2013-2024)')
    parser.add_argument('--disease-year', type=int, default=1, help='Year disease starts (1 = 2013)')
    parser.add_argument('--D_L', type=float, default=400.0, help='Larval dispersal scale (km)')
    args = parser.parse_args()
    
    # Parse seeds
    if args.seeds:
        seeds = [int(s) for s in args.seeds.split(',')]
    elif args.seed is not None:
        seeds = [args.seed]
    else:
        seeds = [42]
    
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load config file (may be nested {param_overrides: {...}, K: ...} or flat)
    with open(args.config) as f:
        config_data = json.load(f)
    
    # Unwrap nested format if present
    if 'param_overrides' in config_data and isinstance(config_data['param_overrides'], dict):
        param_overrides = config_data['param_overrides']
    else:
        param_overrides = config_data

    # Extract sentinel config if present
    sentinel_config = config_data.get('sentinel', None)
    
    # Extract K_cv from CLI or param_overrides
    K_cv = args.K_cv
    if 'population.K_cv' in param_overrides:
        K_cv = float(param_overrides.pop('population.K_cv'))

    K_ref = args.K_ref
    density_scale = K_ref / args.K if args.K > 0 else 1.0

    print(f"=== Agentic Calibration Runner ===")
    print(f"Seeds: {seeds}")
    print(f"K: {args.K}, K_ref: {K_ref}, density_scale: {density_scale:.2f}, K_cv: {K_cv}")
    print(f"Years: {args.years}, Disease year: {args.disease_year}")
    print(f"Overrides: {len(param_overrides)} params")
    
    # Extract spatial params from overrides (needed before network build)
    D_P = param_overrides.pop('spatial.D_P', 15.0)
    D_P_max_range = param_overrides.pop('spatial.D_P_max_range', None)
    D_L = param_overrides.pop('spatial.D_L', args.D_L)
    n_connectivity = float(param_overrides.pop('spatial.n_connectivity', 0.3))
    alpha_self_open = float(param_overrides.pop('spatial.alpha_self_open', 0.05))
    alpha_self_fjord = float(param_overrides.pop('spatial.alpha_self_fjord', 0.50))
    r_total = float(param_overrides.pop('spatial.r_total', 0.02))
    phi_open = float(param_overrides.pop('spatial.phi_open', 0.8))
    phi_fjord = float(param_overrides.pop('spatial.phi_fjord', 0.03))
    # Store back so they appear in output metadata
    param_overrides['spatial.n_connectivity'] = n_connectivity
    param_overrides['spatial.alpha_self_open'] = alpha_self_open
    param_overrides['spatial.alpha_self_fjord'] = alpha_self_fjord
    param_overrides['spatial.r_total'] = r_total
    param_overrides['spatial.phi_open'] = phi_open
    param_overrides['spatial.phi_fjord'] = phi_fjord
    
    # Build network (once, shared across seeds)
    print(f"\nBuilding {len(load_sites())}-node network (D_P={D_P}, max_range={D_P_max_range}, n_conn={n_connectivity}, α_self=[{alpha_self_open:.2f},{alpha_self_fjord:.2f}], r_total={r_total}, φ=[{phi_open:.2f},{phi_fjord:.3f}])...")
    t0 = time.time()
    sites, network = build_full_network(K=args.K, seed=seeds[0], D_L=D_L,
                                         D_P=D_P, D_P_max_range=D_P_max_range,
                                         K_cv=K_cv, n_connectivity=n_connectivity,
                                         alpha_self_open=alpha_self_open,
                                         alpha_self_fjord=alpha_self_fjord,
                                         r_total=r_total,
                                         phi_open=phi_open,
                                         phi_fjord=phi_fjord,
                                         K_ref=K_ref)
    print(f"  Network built in {time.time()-t0:.1f}s ({len(sites)} nodes)")
    
    # Build config
    config = default_config()
    config.simulation.sst_source = 'monthly'
    config.simulation.sst_data_dir = str(PROJECT_ROOT / 'data' / 'sst' / 'site_sst')
    config.simulation.sst_start_year = 2012  # align SST with simulation calendar
    apply_param_overrides(config, param_overrides)

    # Apply sentinel config if present
    if sentinel_config and isinstance(sentinel_config, dict):
        from sswd_evoepi.config import SentinelSection
        config.sentinel = SentinelSection(**{
            k: v for k, v in sentinel_config.items()
            if k in {'enabled', 'n_per_site', 'shedding_fraction', 'site_filter', 'custom_node_ids'}
        })

    # Density-invariant K scaling: set disease shedding scale factor
    config.disease.density_scale = density_scale
    if density_scale != 1.0:
        print(f"  Density-invariant scaling: shedding ×{density_scale:.1f}, "
              f"habitat area ×{args.K/K_ref:.2f}")
    
    # Run simulations
    round_id = output_dir.name  # e.g. "round_00"
    results = []
    for i, seed in enumerate(seeds):
        print(f"\n--- Seed {seed} ({i+1}/{len(seeds)}) ---")
        result = run_single(config, sites, network, seed,
                           disease_year=args.disease_year, n_years=args.years,
                           round_id=round_id, output_dir=output_dir)
        results.append(result)
        
        # Save individual result
        with open(output_dir / f'result_seed{seed}.json', 'w') as f:
            json.dump(result, f, indent=2)
        print(f"  RMSLE = {result['scoring'].get('rmsle', result['scoring'].get('rmse_log', float('inf'))):.3f}, "
              f"crash = {result['overall']['pop_crash_pct']:.1f}%, "
              f"time = {result['wall_time_seconds']:.0f}s")
    
    # Print summary
    print_summary(results, param_overrides)
    
    # Save combined results
    combined = {
        'param_overrides': param_overrides,
        'K': args.K,
        'K_ref': K_ref,
        'density_scale': density_scale,
        'K_cv': K_cv,
        'years': args.years,
        'disease_year': args.disease_year,
        'seeds': seeds,
        'results': results,
        'mean_rmsle': float(np.mean([r['scoring'].get('rmsle', r['scoring'].get('rmse_log', float('inf'))) for r in results])),
        'mean_rmse_log': float(np.mean([r['scoring'].get('rmsle', r['scoring'].get('rmse_log', float('inf'))) for r in results])),  # backward compat
    }
    with open(output_dir / 'combined_results.json', 'w') as f:
        json.dump(combined, f, indent=2)
    
    print(f"\nResults saved to {output_dir}/")
    return combined['mean_rmse_log']


if __name__ == '__main__':
    main()
