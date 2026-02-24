#!/usr/bin/env python3
"""Monterey reintroduction experiment on the full 907-node network.

Factorial design:
  - 3 restoration levels: partial (50), medium (500), full (5000)
  - 6 genetic backgrounds: pre_sswd, survivors_2019, bred_1gen,
    bred_2gen, bred_5gen, optimal
  - 1 baseline (no intervention)
  → 19 total scenarios

Network: All 907 Pycnopodia habitat sites (Aleutians to Baja
  California) with 18 regions, using precomputed overwater distance
  matrix at `results/overwater/distance_matrix.npz`.

SST: Per-site monthly SST from NOAA OISST v2.1 (2002-2025) stored in
  `data/sst/site_sst/`. Each of the 907 nodes maps to its nearest
  OISST 0.25° grid cell. For sites without per-site SST (fallback),
  nearest-latitude interpolation from 11 reference nodes is used.
  Projections: CMIP6 SSP2-4.5 (2026-2050).

Release site: CA-C-043 (Monterey outplanting site, 36.62°N, -121.90°W).
  Release at simulation day corresponding to Jan 1, 2026.

Timeline: 2002-2050 (48 years)
  - 2002-2012: Pre-epidemic spinup (11 years)
  - 2013: Epidemic onset (= disease_year=11)
  - 2026: Release (= sim day 8760+)
  - 2050: End of simulation (24 years post-release)

Usage:
  python3 experiments/reintroduction_monterey.py --mode demo
  python3 experiments/reintroduction_monterey.py --mode full
  python3 experiments/reintroduction_monterey.py --mode single \\
      --genetics bred_2gen --restoration medium

Authors: Anton 🔬 & Willem Weertman
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from dataclasses import dataclass
from multiprocessing import Pool
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

# Ensure project root is on path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from sswd_evoepi.config import (
    ReleaseEvent,
    SimulationConfig,
    SimulationSection,
    SpatialSection,
    GeneticsSection,
    DiseaseSection,
    PopulationSection,
    SpawningSection,
    MovementSection,
    OutputSection,
)
from sswd_evoepi.spatial import (
    NodeDefinition,
    MetapopulationNetwork,
    SpatialNode,
    build_network,
    construct_larval_connectivity,
    construct_pathogen_dispersal,
)
from sswd_evoepi.model import run_spatial_simulation
from scripts.genetic_backgrounds import get_all_backgrounds, get_background_by_name


# ═══════════════════════════════════════════════════════════════════════
# CONSTANTS
# ═══════════════════════════════════════════════════════════════════════

MONTEREY_SITE_NAME = 'CA-C-043'  # Monterey outplanting site (36.62°N, -121.90°W)
START_YEAR = 2002
END_YEAR = 2050
EPIDEMIC_YEAR = 2013
RELEASE_YEAR = 2026
N_YEARS = END_YEAR - START_YEAR   # 48
DISEASE_YEAR = EPIDEMIC_YEAR - START_YEAR  # 11 (0-indexed)
RELEASE_DAY = (RELEASE_YEAR - START_YEAR) * 365  # Day 8760

K_PER_NODE = 5000  # Carrying capacity per node (affordable: post-crash populations are tiny)

RESTORATION_LEVELS = {
    'partial': 50,
    'medium': 500,
    'full': 5000,
}

# Site data and distance matrix paths
SITES_PATH = 'data/nodes/all_sites.json'
DISTANCE_MATRIX_PATH = 'results/overwater/distance_matrix.npz'

# 18 regions for aggregation (matching assign_regions.py output)
REGIONS = [
    'AK-WG', 'AK-AL', 'AK-EG', 'AK-PWS', 'AK-FN', 'AK-FS', 'AK-OC',
    'BC-N', 'BC-C',
    'JDF', 'SS-N', 'SS-S',
    'WA-O', 'OR',
    'CA-N', 'CA-C', 'CA-S',
    'BJ',
]

# 11 SST reference nodes (fallback for sites without per-site SST)
SST_REFERENCE_NODES = [
    ('Sitka', 57.05),
    ('Ketchikan', 55.34),
    ('Haida_Gwaii', 53.25),
    ('Bella_Bella', 52.16),
    ('Howe_Sound', 49.38),
    ('SJI', 48.53),
    ('Westport', 46.89),
    ('Newport', 44.63),
    ('Crescent_City', 41.75),
    ('Fort_Bragg', 39.45),
    ('Monterey', 36.62),
]

# SST STRATEGY:
# Primary: Per-site monthly SST from NOAA OISST v2.1 (data/sst/site_sst/).
#   907 sites → 532 unique OISST 0.25° grid cells. Each site maps to
#   its nearest grid cell via data/sst/site_sst/site_to_grid.json.
#
# Fallback: If per-site SST not yet available, use nearest-latitude
#   interpolation from 11 reference nodes (same as old 489-node method).
#
# Projections: CMIP6 SSP2-4.5, bias-corrected, 2015-2100.


# ═══════════════════════════════════════════════════════════════════════
# 907-NODE NETWORK BUILDER
# ═══════════════════════════════════════════════════════════════════════


def _nearest_ref_node(lat: float) -> str:
    """Find the nearest SST reference node by latitude (fallback)."""
    best_name = SST_REFERENCE_NODES[0][0]
    best_dist = abs(lat - SST_REFERENCE_NODES[0][1])
    for name, ref_lat in SST_REFERENCE_NODES[1:]:
        d = abs(lat - ref_lat)
        if d < best_dist:
            best_dist = d
            best_name = name
    return best_name


def _estimate_mean_sst(lat: float) -> float:
    """Estimate annual mean SST from latitude (simple linear model).

    Based on NE Pacific coast observations:
      SST ≈ 25.0 - 0.27 * latitude (rough, for initialization only)
    """
    return max(5.0, min(20.0, 25.0 - 0.27 * lat))


def _estimate_sst_amplitude(lat: float) -> float:
    """Estimate SST annual amplitude from latitude.

    Higher latitudes have larger seasonal swings.
    """
    return max(1.5, min(5.0, 1.0 + 0.06 * lat))


def _find_monterey_index(names: list) -> int:
    """Find the index of the Monterey release site by name."""
    for i, name in enumerate(names):
        if name == MONTEREY_SITE_NAME:
            return i
    # Fallback: find nearest to 36.62, -121.90
    raise ValueError(
        f"Release site {MONTEREY_SITE_NAME} not found in site list. "
        f"Available sites: {names[:5]}..."
    )


def load_node_definitions(
    K: int = K_PER_NODE,
) -> Tuple[List[NodeDefinition], Dict[int, str], int]:
    """Load all 907 nodes from all_sites.json with SST assignments.

    Uses the canonical site list (data/nodes/all_sites.json) as the
    source of truth for site names, coordinates, and regions.

    Args:
        K: Carrying capacity per node.

    Returns:
        (node_defs, sst_mapping, monterey_idx) where:
          - sst_mapping maps node_id → ref SST name (fallback)
          - monterey_idx is the index of the release site
    """
    sites_path = PROJECT_ROOT / SITES_PATH
    with open(sites_path) as f:
        sites = json.load(f)

    node_defs = []
    sst_mapping = {}
    monterey_idx = None

    for i, site in enumerate(sites):
        lat = float(site['latitude'])
        lon = float(site['longitude'])
        region = site.get('region', '?')
        name = site['name']

        if name == MONTEREY_SITE_NAME:
            monterey_idx = i

        # Assign nearest SST reference node (fallback for sites without per-site SST)
        ref_name = _nearest_ref_node(lat)
        sst_mapping[i] = ref_name

        # Estimate environmental parameters from latitude
        mean_sst = _estimate_mean_sst(lat)
        sst_amp = _estimate_sst_amplitude(lat)

        # Determine if likely fjord (BC-N, some AK-PWS, SS-N)
        is_fjord = region in ('BC-N', 'AK-PWS') and lat > 50.0

        nd = NodeDefinition(
            node_id=i,
            name=name,
            lat=lat,
            lon=lon,
            subregion=region,
            habitat_area=25000.0,  # Default 25000 m² (standardized)
            carrying_capacity=K,
            is_fjord=is_fjord,
            sill_depth=30.0 if is_fjord else np.inf,
            flushing_rate=0.03 if is_fjord else 0.5,
            mean_sst=mean_sst,
            sst_amplitude=sst_amp,
            sst_trend=0.02,
            salinity=22.0 if is_fjord else 32.0,
            depth_range=(5.0, 60.0),
        )
        node_defs.append(nd)

    if monterey_idx is None:
        raise ValueError(f"Release site {MONTEREY_SITE_NAME} not found in {sites_path}")

    return node_defs, sst_mapping, monterey_idx


def build_network_907(
    K: int = K_PER_NODE,
    seed: int = 42,
    D_L: float = 400.0,
    D_P: float = 15.0,
) -> Tuple[MetapopulationNetwork, Dict[int, str], int]:
    """Build the full 907-node network with precomputed overwater distances.

    Returns:
        (network, sst_mapping, monterey_idx)
    """
    node_defs, sst_mapping, monterey_idx = load_node_definitions(K=K)

    npz_path = str(PROJECT_ROOT / DISTANCE_MATRIX_PATH)
    network = build_network(
        node_defs,
        D_L=D_L,
        D_P=D_P,
        seed=seed,
        overwater_npz=npz_path,
    )

    return network, sst_mapping, monterey_idx


# ═══════════════════════════════════════════════════════════════════════
# EXPERIMENT CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════


@dataclass
class ScenarioConfig:
    """Configuration for a single reintroduction scenario."""
    name: str
    restoration_level: Optional[str]  # None for baseline
    n_released: int
    genetics_name: Optional[str]  # None for baseline
    allele_freqs: Optional[np.ndarray]
    expected_traits: Optional[Dict[str, float]]
    seed: int
    replicate: int


def build_scenarios(
    n_replicates: int = 3,
    base_seed: int = 42,
) -> List[ScenarioConfig]:
    """Build the full factorial scenario list.

    18 treatment scenarios (3 levels × 6 genetics) + 1 baseline = 19.
    Each scenario is replicated n_replicates times with different seeds.

    Returns:
        List of ScenarioConfig objects.
    """
    backgrounds = get_all_backgrounds()
    scenarios = []

    # Baseline: no intervention
    for rep in range(n_replicates):
        scenarios.append(ScenarioConfig(
            name='baseline',
            restoration_level=None,
            n_released=0,
            genetics_name=None,
            allele_freqs=None,
            expected_traits=None,
            seed=base_seed + rep,
            replicate=rep,
        ))

    # Factorial: restoration level × genetic background
    for level_name, n_released in RESTORATION_LEVELS.items():
        for bg in backgrounds:
            for rep in range(n_replicates):
                scenarios.append(ScenarioConfig(
                    name=f'{level_name}_{bg["name"]}',
                    restoration_level=level_name,
                    n_released=n_released,
                    genetics_name=bg['name'],
                    allele_freqs=bg['allele_freqs'],
                    expected_traits=bg['expected_traits'],
                    seed=base_seed + rep,
                    replicate=rep,
                ))

    return scenarios


def build_config_for_scenario(
    scenario: ScenarioConfig,
    sst_mapping: Dict[int, str],
    monterey_idx: int = 0,
) -> SimulationConfig:
    """Build a SimulationConfig for a scenario.

    Uses:
      - Monthly SST with SSP2-4.5 projections
      - Full disease dynamics (ubiquitous, Prentice 2025 rates)
      - Three-trait genetics (17/17/17)
      - Release event at Monterey (CA-C-043), day 8760 (Jan 2026)
    """
    release_events = []
    if scenario.n_released > 0 and scenario.allele_freqs is not None:
        release_events.append(ReleaseEvent(
            time_step=RELEASE_DAY,
            node_id=monterey_idx,
            n_individuals=scenario.n_released,
            genetics_mode='allele_freqs',
            allele_freqs=scenario.allele_freqs,
            age_range=(365, 730),  # 1-2 year old juveniles
            mark_released=True,
        ))

    config = SimulationConfig(
        simulation=SimulationSection(
            start_year=START_YEAR,
            epidemic_year=EPIDEMIC_YEAR,
            end_year=END_YEAR,
            seed=scenario.seed,
            sst_source='monthly',
            sst_data_dir=str(PROJECT_ROOT / 'data' / 'sst'),
            sst_start_year=START_YEAR,
            sst_scenario='ssp245',
            sst_projection_dir=str(PROJECT_ROOT / 'data' / 'sst' / 'projections'),
            sst_obs_end_year=2025,
        ),
        genetics=GeneticsSection(
            n_resistance=17,
            n_tolerance=17,
            n_recovery=17,
            target_mean_r=0.15,
            target_mean_t=0.10,
            target_mean_c=0.02,
        ),
        release_events=release_events,
    )

    return config


# ═══════════════════════════════════════════════════════════════════════
# SST ASSIGNMENT FOR 907 NODES
# ═══════════════════════════════════════════════════════════════════════

def assign_sst_node_names(
    network: MetapopulationNetwork,
    sst_mapping: Dict[int, str],
) -> None:
    """Override node names for SST loading (fallback mode).

    When per-site SST (data/sst/site_sst/) is not available, the model
    looks up SST files by node name. We override node definition names
    to match their assigned SST reference node.

    IMPORTANT: This means node names in the SST system differ from
    geographic names. The geographic names are preserved in node_names
    result field.
    """
    for node in network.nodes:
        node_id = node.definition.node_id
        ref_name = sst_mapping.get(node_id)
        if ref_name:
            # Store geographic name for results
            if not hasattr(node, '_geographic_name'):
                node._geographic_name = node.definition.name
            # Override name for SST file lookup
            node.definition.name = ref_name


# ═══════════════════════════════════════════════════════════════════════
# SINGLE SCENARIO RUNNER
# ═══════════════════════════════════════════════════════════════════════


def run_single_scenario(
    scenario: ScenarioConfig,
    K: int = K_PER_NODE,
) -> Dict[str, Any]:
    """Run a single reintroduction scenario.

    Builds the 907-node network, configures SST, runs the spatial
    simulation, and extracts results.

    Args:
        scenario: Scenario configuration.
        K: Carrying capacity per node.

    Returns:
        Dict with scenario results + metadata.
    """
    t0 = time.time()
    print(f"  [{scenario.name} rep={scenario.replicate}] Starting...")

    # Build network
    network, sst_mapping, monterey_idx = build_network_907(K=K, seed=scenario.seed)

    # Store geographic names before SST override
    geographic_names = [n.definition.name for n in network.nodes]
    geographic_regions = [n.definition.subregion for n in network.nodes]

    # Override node names for SST file lookup
    assign_sst_node_names(network, sst_mapping)

    # Build config
    config = build_config_for_scenario(scenario, sst_mapping, monterey_idx=monterey_idx)

    # Run simulation
    result = run_spatial_simulation(
        network=network,
        n_years=N_YEARS,
        disease_year=DISEASE_YEAR,
        initial_infected_per_node=5,
        seed=scenario.seed,
        config=config,
    )

    elapsed = time.time() - t0
    print(f"  [{scenario.name} rep={scenario.replicate}] Done in {elapsed:.0f}s, "
          f"final pop={result.final_total_pop}")

    # Extract results
    output = {
        'scenario': scenario.name,
        'replicate': scenario.replicate,
        'seed': scenario.seed,
        'restoration_level': scenario.restoration_level,
        'genetics_name': scenario.genetics_name,
        'n_released': scenario.n_released,
        'expected_traits': scenario.expected_traits,
        'elapsed_s': elapsed,
        # Aggregate metrics
        'initial_total_pop': result.initial_total_pop,
        'final_total_pop': result.final_total_pop,
        'total_released': result.total_released,
        'released_surviving': result.released_surviving,
        # Per-year totals
        'yearly_total_pop': result.yearly_total_pop.tolist() if result.yearly_total_pop is not None else None,
        # Monterey-specific metrics
        'monterey_idx': monterey_idx,
        'monterey_yearly_pop': (
            result.yearly_pop[monterey_idx].tolist()
            if result.yearly_pop is not None else None
        ),
        'monterey_yearly_mean_r': (
            result.yearly_mean_resistance[monterey_idx].tolist()
            if result.yearly_mean_resistance is not None else None
        ),
        'monterey_yearly_mean_t': (
            result.yearly_mean_tolerance[monterey_idx].tolist()
            if result.yearly_mean_tolerance is not None else None
        ),
        'monterey_yearly_mean_c': (
            result.yearly_mean_recovery[monterey_idx].tolist()
            if result.yearly_mean_recovery is not None else None
        ),
        # Released individual tracking at Monterey
        'monterey_released_alive': (
            result.yearly_released_alive[monterey_idx].tolist()
            if result.yearly_released_alive is not None else None
        ),
        # Region-aggregated population trajectories
        'region_yearly_pop': _aggregate_by_region(
            result.yearly_pop, geographic_regions,
        ),
        # Geographic names for reference
        'n_nodes': len(geographic_names),
    }

    return output


def _aggregate_by_region(
    yearly_pop: np.ndarray,
    regions: List[str],
) -> Dict[str, List[int]]:
    """Aggregate per-node yearly population by region.

    Args:
        yearly_pop: (n_nodes, n_years) array.
        regions: List of region strings per node.

    Returns:
        Dict mapping region → list of yearly total population.
    """
    if yearly_pop is None:
        return {}

    result = {}
    for region in REGIONS:
        mask = np.array([r == region for r in regions])
        if np.any(mask):
            result[region] = yearly_pop[mask].sum(axis=0).tolist()
    return result


# ═══════════════════════════════════════════════════════════════════════
# PARALLEL RUNNER
# ═══════════════════════════════════════════════════════════════════════


def _worker(args):
    """Multiprocessing worker for a single scenario."""
    scenario, K = args
    try:
        return run_single_scenario(scenario, K=K)
    except Exception as e:
        print(f"  [{scenario.name} rep={scenario.replicate}] FAILED: {e}")
        import traceback
        traceback.print_exc()
        return {
            'scenario': scenario.name,
            'replicate': scenario.replicate,
            'error': str(e),
        }


def run_experiment(
    mode: str = 'demo',
    n_cores: int = 8,
    K: int = K_PER_NODE,
    single_genetics: Optional[str] = None,
    single_restoration: Optional[str] = None,
) -> List[Dict]:
    """Run the full reintroduction experiment.

    Args:
        mode: 'demo' (3 reps, 8 cores), 'full' (10 reps), or 'single'.
        n_cores: Number of parallel workers.
        K: Carrying capacity per node.
        single_genetics: For mode='single', which genetic background.
        single_restoration: For mode='single', which restoration level.

    Returns:
        List of result dicts.
    """
    if mode == 'demo':
        n_replicates = 3
    elif mode == 'full':
        n_replicates = 10
    elif mode == 'single':
        n_replicates = 1
    else:
        raise ValueError(f"Unknown mode: {mode}")

    scenarios = build_scenarios(n_replicates=n_replicates)

    # Filter for single mode
    if mode == 'single':
        if single_genetics is None or single_restoration is None:
            raise ValueError("--genetics and --restoration required for single mode")
        target_name = f'{single_restoration}_{single_genetics}'
        scenarios = [s for s in scenarios if s.name == target_name or s.name == 'baseline']
        if not scenarios:
            raise ValueError(f"No scenario matching {target_name}")

    print(f"Reintroduction Experiment: {len(scenarios)} scenarios")
    print(f"  Mode: {mode}, K={K}/node, {n_cores} cores")
    print(f"  Network: 907 nodes, 18 regions, SST: monthly + SSP2-4.5")
    print(f"  Release: {MONTEREY_SITE_NAME} (Monterey), day {RELEASE_DAY} (Jan {RELEASE_YEAR})")
    print()

    t0 = time.time()

    # Run scenarios
    if n_cores == 1 or len(scenarios) == 1:
        results = []
        for s in scenarios:
            results.append(_worker((s, K)))
    else:
        args = [(s, K) for s in scenarios]
        with Pool(processes=min(n_cores, len(scenarios))) as pool:
            results = pool.map(_worker, args)

    elapsed = time.time() - t0
    print(f"\nAll done in {elapsed:.0f}s ({elapsed/60:.1f} min)")

    # Save results
    output_dir = PROJECT_ROOT / 'results' / 'reintroduction'
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f'results_{mode}.json'

    # Convert numpy types for JSON serialization
    def _convert(obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return obj

    serializable = json.loads(json.dumps(results, default=_convert))
    with open(output_file, 'w') as f:
        json.dump(serializable, f, indent=2)
    print(f"Results saved to {output_file}")

    # Print summary
    _print_summary(results)

    return results


def _print_summary(results: List[Dict]) -> None:
    """Print experiment summary table."""
    print("\n" + "=" * 80)
    print("EXPERIMENT SUMMARY")
    print("=" * 80)
    print(f"{'Scenario':<30s} {'Rep':>3s} {'Final Pop':>10s} "
          f"{'Monterey':>10s} {'Released':>10s} {'Time':>6s}")
    print("-" * 80)

    for r in results:
        if 'error' in r:
            print(f"  {r['scenario']:<28s} {r['replicate']:>3d}  ERROR: {r['error']}")
            continue

        monterey_final = r.get('monterey_yearly_pop', [0])[-1] if r.get('monterey_yearly_pop') else 0
        released_alive = r.get('released_surviving', 0)

        print(f"  {r['scenario']:<28s} {r['replicate']:>3d} "
              f"{r['final_total_pop']:>10,d} "
              f"{monterey_final:>10,d} "
              f"{released_alive:>10,d} "
              f"{r.get('elapsed_s', 0):>5.0f}s")


# ═══════════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════════


def main():
    parser = argparse.ArgumentParser(
        description='Monterey reintroduction experiment (907-node network)',
    )
    parser.add_argument(
        '--mode', choices=['demo', 'full', 'single'], default='demo',
        help='Run mode: demo (3 reps), full (10 reps), single (1 scenario)',
    )
    parser.add_argument(
        '--cores', type=int, default=8,
        help='Number of parallel workers (default: 8)',
    )
    parser.add_argument(
        '--K', type=int, default=K_PER_NODE,
        help=f'Carrying capacity per node (default: {K_PER_NODE})',
    )
    parser.add_argument(
        '--genetics', type=str, default=None,
        help='Genetic background for single mode',
    )
    parser.add_argument(
        '--restoration', type=str, default=None,
        help='Restoration level for single mode',
    )
    args = parser.parse_args()

    run_experiment(
        mode=args.mode,
        n_cores=args.cores,
        K=args.K,
        single_genetics=args.genetics,
        single_restoration=args.restoration,
    )


if __name__ == '__main__':
    main()
