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

# ── Paths ─────────────────────────────────────────────────────────────
SITES_PATH = 'data/nodes/all_sites.json'
DISTANCE_MATRIX_PATH = 'results/overwater/distance_matrix.npz'

# ── Calibration targets: recovery fraction by ~2024 (11yr post-crash) ──
RECOVERY_TARGETS = {
    'AK-PWS': 0.50,   # Prince William Sound — 50%
    'AK-FN': 0.50,    # North inshore SE Alaska — 50%
    'AK-FS': 0.20,    # South inshore SE Alaska — 20%
    'BC-N': 0.20,     # Northern BC — 20%
    'SS-S': 0.05,     # Puget Sound — 5%
    'JDF': 0.02,      # Juan de Fuca — 2%
    'OR': 0.0025,     # Oregon outer coast — 0.25%
    'CA-N': 0.001,    # Northern California — 0.1%
}


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


def build_node_defs(sites: List[dict], K: int = 5000) -> List[NodeDefinition]:
    """Build NodeDefinition list from site JSON."""
    node_defs = []
    for i, site in enumerate(sites):
        lat = float(site['latitude'])
        lon = float(site['longitude'])
        region = site.get('region', '?')
        is_fjord = region in ('BC-N', 'AK-PWS') and lat > 50.0
        
        nd = NodeDefinition(
            node_id=i,
            name=site['name'],
            lat=lat,
            lon=lon,
            subregion=region,
            habitat_area=25000.0,
            carrying_capacity=K,
            is_fjord=is_fjord,
            sill_depth=30.0 if is_fjord else np.inf,
            flushing_rate=0.03 if is_fjord else 0.5,
            mean_sst=_estimate_mean_sst(lat),
            sst_amplitude=_estimate_sst_amplitude(lat),
            sst_trend=0.02,
            salinity=22.0 if is_fjord else 32.0,
            depth_range=(5.0, 60.0),
        )
        node_defs.append(nd)
    return node_defs


def build_full_network(K: int = 5000, seed: int = 42, D_L: float = 400.0, D_P: float = 15.0):
    """Build 896-node network with overwater distances."""
    sites = load_sites()
    node_defs = build_node_defs(sites, K=K)
    npz_path = str(PROJECT_ROOT / DISTANCE_MATRIX_PATH)
    
    network = build_network(
        node_defs,
        D_L=D_L,
        D_P=D_P,
        seed=seed,
        overwater_npz=npz_path,
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
        
        region_details[region] = {
            'n_nodes': len(idxs),
            'peak_pop': int(peak_pop),
            'final_pop': int(final_pop),
            'recovery_frac': float(recovery),
            'crash_pct': float(100.0 * (1 - final_pop / peak_pop)) if peak_pop > 0 else 100.0,
            'yearly_totals': [int(p) for p in region_pop],
            'yearly_recruits': region_recruits,
            'yearly_disease_deaths': region_dd,
        }
    
    return region_recovery, region_details


def score_against_targets(region_recovery: dict) -> dict:
    """Score using log-space RMSE (targets span 3 orders of magnitude)."""
    scores = {}
    total_log_sq_error = 0.0
    n_targets = 0
    eps = 1e-6
    
    for region, target in RECOVERY_TARGETS.items():
        actual = region_recovery.get(region, 0.0)
        
        log_target = np.log10(max(target, eps))
        log_actual = np.log10(max(actual, eps))
        log_error = log_actual - log_target
        
        scores[region] = {
            'target': float(target),
            'actual': float(actual),
            'target_pct': float(target * 100),
            'actual_pct': float(actual * 100),
            'log_error': float(log_error),
            'log_sq_error': float(log_error ** 2),
            'within_2x': bool(abs(log_error) < 0.301),  # log10(2)
            'within_5x': bool(abs(log_error) < 0.699),  # log10(5)
        }
        
        total_log_sq_error += log_error ** 2
        n_targets += 1
    
    rmse_log = np.sqrt(total_log_sq_error / n_targets) if n_targets > 0 else float('inf')
    within_2x = sum(1 for s in scores.values() if s['within_2x'])
    within_5x = sum(1 for s in scores.values() if s['within_5x'])
    
    return {
        'per_region': scores,
        'rmse_log': float(rmse_log),
        'within_2x': within_2x,
        'within_5x': within_5x,
        'n_targets': n_targets,
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

    t0 = time.time()
    result = run_spatial_simulation(
        network=network,
        n_years=n_years,
        disease_year=disease_year,
        seed=seed,
        config=config,
        progress_callback=checker,
    )
    elapsed = time.time() - t0

    # ── Handle early stop ────────────────────────────────────────────
    if checker.early_stop_reason:
        print(f"  ⚡ EARLY STOP: {checker.early_stop_reason} ({elapsed:.0f}s)")
        return {
            'seed': seed,
            'wall_time_seconds': float(elapsed),
            'scoring': {
                'rmse_log': float('inf'),
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

    # ── Normal scoring ───────────────────────────────────────────────
    region_recovery, region_details = compute_regional_recovery(result, sites)
    scoring = score_against_targets(region_recovery)

    return {
        'seed': seed,
        'wall_time_seconds': float(elapsed),
        'scoring': scoring,
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
    
    avg_rmse = float(np.mean([r['scoring']['rmse_log'] for r in results]))
    avg_crash = float(np.mean([r['overall']['pop_crash_pct'] for r in results]))
    
    print(f"\nOverall crash: {avg_crash:.1f}%")
    print(f"RMSE(log10):   {avg_rmse:.3f}  (0 = perfect, <0.3 = within 2×)")
    
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
    
    # Load param overrides
    with open(args.config) as f:
        param_overrides = json.load(f)
    
    print(f"=== Agentic Calibration Runner ===")
    print(f"Seeds: {seeds}")
    print(f"K: {args.K}, Years: {args.years}, Disease year: {args.disease_year}")
    print(f"Overrides: {len(param_overrides)} params")
    
    # Build network (once, shared across seeds)
    print(f"\nBuilding {len(load_sites())}-node network...")
    t0 = time.time()
    sites, network = build_full_network(K=args.K, seed=seeds[0], D_L=args.D_L)
    print(f"  Network built in {time.time()-t0:.1f}s ({len(sites)} nodes)")
    
    # Build config
    config = default_config()
    config.simulation.sst_source = 'monthly'
    config.simulation.sst_data_dir = str(PROJECT_ROOT / 'data' / 'sst' / 'site_sst')
    config.simulation.sst_start_year = 2012  # align SST with simulation calendar
    apply_param_overrides(config, param_overrides)
    
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
        print(f"  RMSE(log10) = {result['scoring']['rmse_log']:.3f}, "
              f"crash = {result['overall']['pop_crash_pct']:.1f}%, "
              f"time = {result['wall_time_seconds']:.0f}s")
    
    # Print summary
    print_summary(results, param_overrides)
    
    # Save combined results
    combined = {
        'param_overrides': param_overrides,
        'K': args.K,
        'years': args.years,
        'disease_year': args.disease_year,
        'seeds': seeds,
        'results': results,
        'mean_rmse_log': float(np.mean([r['scoring']['rmse_log'] for r in results])),
    }
    with open(output_dir / 'combined_results.json', 'w') as f:
        json.dump(combined, f, indent=2)
    
    print(f"\nResults saved to {output_dir}/")
    return combined['mean_rmse_log']


if __name__ == '__main__':
    main()
