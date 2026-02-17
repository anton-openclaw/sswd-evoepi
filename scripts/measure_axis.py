#!/usr/bin/env python3
"""Modular scaling measurement — run one axis at a time.

Usage:
    python3 scripts/measure_axis.py single_N
    python3 scripts/measure_axis.py single_T
    python3 scripts/measure_axis.py NxT_matrix
    python3 scripts/measure_axis.py spawning_overhead
    python3 scripts/measure_axis.py spatial_K
    python3 scripts/measure_axis.py KxN_stress
    python3 scripts/measure_axis.py parallel
    python3 scripts/measure_axis.py disease_spatial

Each axis saves to results/performance/axis_<name>.json
"""
import sys, os, time, json, traceback, signal
sys.path.insert(0, '.')

import numpy as np
from sswd_evoepi.model import run_coupled_simulation, run_spatial_simulation
from sswd_evoepi.config import default_config
from sswd_evoepi.perf import PerfMonitor

os.makedirs('results/performance', exist_ok=True)


def get_memory_mb():
    try:
        with open('/proc/self/status') as f:
            for line in f:
                if line.startswith('VmRSS:'):
                    return int(line.split()[1]) / 1024
    except:
        return -1


def single_run(n, years, spawning=True, disease_year=999, seed=42):
    """Single-node simulation run with timing and memory."""
    cfg = default_config()
    if not spawning:
        cfg.spawning = None
    perf = PerfMonitor(enabled=True)
    mem_before = get_memory_mb()
    t0 = time.perf_counter()
    try:
        r = run_coupled_simulation(
            n_individuals=n, carrying_capacity=n,
            habitat_area=max(10000, n * 200),
            T_celsius=14.0, salinity=30.0, phi_k=0.02,
            n_years=years, disease_year=disease_year,
            initial_infected=min(5, max(1, n // 10)),
            seed=seed, config=cfg, perf=perf,
        )
        wall = time.perf_counter() - t0
        mem_after = get_memory_mb()
        return {
            'wall_s': round(wall, 4), 'final_pop': r.final_pop,
            'mem_delta_mb': round(mem_after - mem_before, 1),
            'mem_total_mb': round(mem_after, 1),
            'components': perf.summary(), 'ok': True,
        }
    except Exception as e:
        return {'wall_s': round(time.perf_counter() - t0, 4),
                'error': str(e)[:300], 'ok': False}


def spatial_run(n_nodes, K_per_node, years, disease_year=None, seed=42):
    """Spatial simulation with custom node count."""
    from sswd_evoepi.spatial import SpatialNode, NodeDefinition, MetapopulationNetwork
    from scipy.spatial.distance import cdist

    cfg = default_config()
    cfg.spawning = None
    cfg.movement.enabled = False

    nodes_def = []
    for i in range(n_nodes):
        nodes_def.append(NodeDefinition(
            node_id=i, name=f"N{i}",
            latitude=48.0 + i * (10.0 / max(n_nodes, 1)),
            longitude=-123.0 - i * (5.0 / max(n_nodes, 1)),
            carrying_capacity=K_per_node,
            habitat_area=max(10000, K_per_node * 200),
            mean_sst=12.0 + (i / max(n_nodes - 1, 1)) * 4.0,
            sst_amplitude=4.0, salinity=30.0,
            flushing_rate=0.02, is_fjord=(i == 1),
        ))

    positions = np.array([(nd.latitude, nd.longitude) for nd in nodes_def])
    dists = cdist(positions, positions) * 111.0  # approx km

    scale_c = 500.0
    C = np.exp(-dists / scale_c)
    np.fill_diagonal(C, 0)
    rs = C.sum(axis=1, keepdims=True)
    rs[rs == 0] = 1
    C = C / rs

    D = np.exp(-dists / 100.0) * 0.01
    np.fill_diagonal(D, 0)

    network = MetapopulationNetwork(nodes=[], C=C, D=D, n_nodes=n_nodes)
    for nd in nodes_def:
        network.nodes.append(SpatialNode(definition=nd))

    mem_before = get_memory_mb()
    t0 = time.perf_counter()
    try:
        r = run_spatial_simulation(
            network=network, n_years=years,
            disease_year=disease_year, seed=seed, config=cfg,
        )
        wall = time.perf_counter() - t0
        mem_after = get_memory_mb()
        return {
            'wall_s': round(wall, 4), 'final_pop': r.final_total_pop,
            'total_agents': n_nodes * K_per_node,
            'mem_delta_mb': round(mem_after - mem_before, 1),
            'mem_total_mb': round(mem_after, 1), 'ok': True,
        }
    except Exception as e:
        wall = time.perf_counter() - t0
        return {
            'wall_s': round(wall, 4), 'error': str(e)[:300],
            'total_agents': n_nodes * K_per_node, 'ok': False,
            'tb': traceback.format_exc()[-500:],
        }


def save(name, data):
    path = f'results/performance/axis_{name}.json'
    with open(path, 'w') as f:
        json.dump(data, f, indent=2, default=str)
    print(f"\nSaved {len(data)} points to {path}")


# ═══════════════════════════════════════════════════════════════════
# AXIS IMPLEMENTATIONS
# ═══════════════════════════════════════════════════════════════════

def axis_single_N():
    """Sweep N from 25 to 10000, 5yr, spawning on."""
    print("=== AXIS: Single-Node Population Size (N) ===")
    print(f"{'N':>8} {'Time':>10} {'Pop':>8} {'Mem(MB)':>8} {'Top':>25}")
    results = []
    for n in [25, 50, 100, 200, 500, 1000, 2000, 5000, 10000]:
        r = single_run(n, years=5, spawning=True)
        top = ''
        if r['ok'] and 'components' in r:
            c = r['components']
            t = max([(k, v) for k, v in c.items() if k != '_total_s'],
                    key=lambda x: x[1].get('total_s', 0), default=('?', {'pct': 0}))
            top = f"{t[0]} ({t[1]['pct']:.0f}%)"
        pop = r.get('final_pop', 'ERR')
        print(f"{n:>8} {r['wall_s']:>9.3f}s {pop:>8} {r.get('mem_total_mb', -1):>7.0f} {top:>25}")
        results.append({'n': n, **r})
        sys.stdout.flush()
    save('single_N', results)


def axis_single_T():
    """Sweep T from 1 to 200yr, 200 agents, spawning on."""
    print("=== AXIS: Single-Node Duration (T) ===")
    print(f"{'T':>8} {'Time':>10} {'Pop':>8} {'Mem(MB)':>8}")
    results = []
    for t in [1, 2, 5, 10, 20, 50, 100, 200]:
        r = single_run(200, years=t, spawning=True)
        pop = r.get('final_pop', 'ERR')
        print(f"{t:>8} {r['wall_s']:>9.3f}s {pop:>8} {r.get('mem_total_mb', -1):>7.0f}")
        results.append({'years': t, **r})
        sys.stdout.flush()
    save('single_T', results)


def axis_NxT_matrix():
    """N×T interaction matrix."""
    print("=== AXIS: N×T Interaction Matrix ===")
    ns = [50, 200, 500, 1000, 2000, 5000]
    ts = [5, 10, 20, 50]

    print(f"{'':>8}", end='')
    for t in ts:
        print(f"{'T='+str(t):>10}", end='')
    print()

    results = []
    for n in ns:
        print(f"{'N='+str(n):>8}", end='')
        for t in ts:
            r = single_run(n, years=t, spawning=True)
            val = f"{r['wall_s']:.2f}s" if r['ok'] else "FAIL"
            print(f"{val:>10}", end='')
            results.append({'n': n, 'years': t, **r})
            sys.stdout.flush()
        print()
    save('NxT_matrix', results)


def axis_spawning_overhead():
    """Spawning on vs off at multiple scales."""
    print("=== AXIS: Spawning Overhead ===")
    print(f"{'N':>8} {'No Spawn':>10} {'Spawn':>10} {'Overhead':>10} {'Mem':>8}")
    results = []
    for n in [50, 100, 200, 500, 1000, 2000, 5000, 10000]:
        r_off = single_run(n, years=10, spawning=False)
        r_on = single_run(n, years=10, spawning=True)
        oh = r_on['wall_s'] / r_off['wall_s'] if r_off['wall_s'] > 0 else -1
        print(f"{n:>8} {r_off['wall_s']:>9.3f}s {r_on['wall_s']:>9.3f}s {oh:>9.1f}× {r_on.get('mem_total_mb', -1):>7.0f}")
        results.append({
            'n': n, 'no_spawn_s': r_off['wall_s'], 'spawn_s': r_on['wall_s'],
            'overhead': round(oh, 2),
            'mem_no_spawn': r_off.get('mem_total_mb', -1),
            'mem_spawn': r_on.get('mem_total_mb', -1),
        })
        sys.stdout.flush()
    save('spawning_overhead', results)


def axis_spatial_K():
    """Spatial node count scaling — conservative, with fallback."""
    print("=== AXIS: Spatial Node Count (K) ===")
    print(f"{'K':>8} {'N/node':>8} {'Total':>8} {'Time':>10} {'Pop':>8} {'Mem':>8}")
    results = []
    K_per_node = 200
    for k in [2, 3, 5, 10, 20, 30, 50, 75, 100, 150]:
        total = k * K_per_node
        r = spatial_run(n_nodes=k, K_per_node=K_per_node, years=5)
        pop = r.get('final_pop', 'ERR')
        print(f"{k:>8} {K_per_node:>8} {total:>8} {r['wall_s']:>9.3f}s {pop:>8} {r.get('mem_total_mb', -1):>7.0f}")
        results.append({'k': k, 'K_per_node': K_per_node, 'total': total, **r})
        sys.stdout.flush()
        if not r['ok']:
            print(f"  ERROR: {r.get('error', '?')[:100]}")
            # Don't give up — try next point but note the failure
        if r.get('mem_total_mb', 0) > 28000:
            print("  Approaching 31GB memory limit — stopping K sweep.")
            break
    save('spatial_K', results)


def axis_KxN_stress():
    """K×N combined stress test — push total agent count."""
    print("=== AXIS: K×N Combined Stress Test (5yr) ===")
    print(f"{'K':>5} {'N/node':>8} {'Total':>8} {'Time':>10} {'Pop':>8} {'Mem':>8}")
    results = []
    combos = [
        (5, 100), (5, 200), (5, 500), (5, 1000),
        (10, 100), (10, 200), (10, 500), (10, 1000),
        (20, 100), (20, 200), (20, 500),
        (50, 100), (50, 200), (50, 500),
        (100, 100), (100, 200), (100, 500),
        (150, 100), (150, 200), (150, 500),
    ]
    for k, n in combos:
        total = k * n
        r = spatial_run(n_nodes=k, K_per_node=n, years=5)
        pop = r.get('final_pop', 'ERR')
        print(f"{k:>5} {n:>8} {total:>8} {r['wall_s']:>9.3f}s {pop:>8} {r.get('mem_total_mb', -1):>7.0f}")
        results.append({'k': k, 'n_per_node': n, 'total': total, **r})
        sys.stdout.flush()
        if not r['ok']:
            print(f"  ERROR: {r.get('error', '?')[:100]}")
        if r.get('mem_total_mb', 0) > 28000:
            print("  Approaching memory limit — stopping.")
            break
    save('KxN_stress', results)


def axis_parallel():
    """Parallel replicate throughput."""
    import multiprocessing as mp

    def _worker(args):
        n, years, seed = args
        return single_run(n, years, spawning=True, seed=seed)

    print("=== AXIS: Parallel Replicates (200 agents, 10yr) ===")
    print(f"{'Procs':>8} {'Runs':>8} {'Wall':>10} {'Throughput':>12} {'Efficiency':>12}")
    results = []

    # Measure single-run baseline
    t0 = time.perf_counter()
    _worker((200, 10, 0))
    single_time = time.perf_counter() - t0

    for n_proc in [1, 2, 4, 8, 16]:
        n_runs = max(n_proc * 2, 4)  # at least 4 runs for stable timing
        args = [(200, 10, seed) for seed in range(n_runs)]
        t0 = time.perf_counter()
        with mp.Pool(n_proc) as pool:
            pool.map(_worker, args)
        wall = time.perf_counter() - t0
        throughput = n_runs / wall
        ideal_throughput = n_proc / single_time
        efficiency = throughput / ideal_throughput * 100 if ideal_throughput > 0 else 0
        print(f"{n_proc:>8} {n_runs:>8} {wall:>9.2f}s {throughput:>10.1f} r/s {efficiency:>10.0f}%")
        results.append({
            'n_procs': n_proc, 'n_runs': n_runs,
            'wall_s': round(wall, 3), 'throughput': round(throughput, 2),
            'efficiency_pct': round(efficiency, 1),
            'single_run_s': round(single_time, 3),
        })
        sys.stdout.flush()
    save('parallel', results)


def axis_disease_spatial():
    """Disease impact on spatial simulations."""
    print("=== AXIS: Disease Impact — Spatial ===")
    print(f"{'Config':>30} {'Time':>10} {'Pop':>8} {'Mem':>8}")
    results = []

    configs = [
        (5, 200, 10, None, "5n×200, 10yr, no disease"),
        (5, 200, 10, 3, "5n×200, 10yr, disease yr3"),
        (5, 500, 10, None, "5n×500, 10yr, no disease"),
        (5, 500, 10, 3, "5n×500, 10yr, disease yr3"),
        (5, 1000, 10, None, "5n×1000, 10yr, no disease"),
        (5, 1000, 10, 3, "5n×1000, 10yr, disease yr3"),
        (20, 200, 10, None, "20n×200, 10yr, no disease"),
        (20, 200, 10, 3, "20n×200, 10yr, disease yr3"),
        (50, 200, 10, None, "50n×200, 10yr, no disease"),
        (50, 200, 10, 3, "50n×200, 10yr, disease yr3"),
    ]
    for k, n, yrs, dy, label in configs:
        r = spatial_run(k, n, yrs, disease_year=dy)
        pop = r.get('final_pop', 'ERR')
        print(f"{label:>30} {r['wall_s']:>9.3f}s {pop:>8} {r.get('mem_total_mb', -1):>7.0f}")
        results.append({'k': k, 'n': n, 'years': yrs, 'disease_year': dy, 'label': label, **r})
        sys.stdout.flush()
    save('disease_spatial', results)


AXES = {
    'single_N': axis_single_N,
    'single_T': axis_single_T,
    'NxT_matrix': axis_NxT_matrix,
    'spawning_overhead': axis_spawning_overhead,
    'spatial_K': axis_spatial_K,
    'KxN_stress': axis_KxN_stress,
    'parallel': axis_parallel,
    'disease_spatial': axis_disease_spatial,
}

if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1] not in AXES:
        print(f"Usage: {sys.argv[0]} <axis>")
        print(f"Available: {', '.join(AXES.keys())}")
        sys.exit(1)
    axis = sys.argv[1]
    print(f"\nStarting axis measurement: {axis}")
    print(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Memory: {get_memory_mb():.0f} MB\n")
    AXES[axis]()
    print(f"\nFinished: {time.strftime('%Y-%m-%d %H:%M:%S')}")
