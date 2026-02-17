#!/usr/bin/env python3
"""Exhaustive scaling analysis across all model axes.

Measures wall time, memory, and component breakdown across:
  Axis 1: N (agents per node) — single node
  Axis 2: T (simulation years) — single node
  Axis 3: K (number of spatial nodes)
  Axis 4: N×T interaction matrix
  Axis 5: Spawning on/off overhead at scale
  Axis 6: Disease impact at scale
  Axis 7: K×N combined (total system size)
  Axis 8: Parallel replicates (multiprocessing)

Each run has a per-config timeout to prevent hangs.
"""
import sys, os, time, json, signal, traceback
sys.path.insert(0, '.')

import numpy as np
import multiprocessing as mp
from sswd_evoepi.model import run_coupled_simulation, run_spatial_simulation
from sswd_evoepi.config import default_config, SimulationConfig
from sswd_evoepi.perf import PerfMonitor
from sswd_evoepi.spatial import make_5node_network, build_network


# ═══════════════════════════════════════════════════════════════════
# HELPERS
# ═══════════════════════════════════════════════════════════════════

def get_memory_mb():
    """Get current process RSS in MB."""
    try:
        with open('/proc/self/status') as f:
            for line in f:
                if line.startswith('VmRSS:'):
                    return int(line.split()[1]) / 1024  # KB to MB
    except:
        return -1


def single_node_run(n, years, spawning=True, disease_year=999, seed=42, timeout_s=120):
    """Run single-node simulation with timeout and memory tracking."""
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
            'wall_s': round(wall, 4),
            'final_pop': r.final_pop,
            'mem_mb': round(mem_after - mem_before, 1) if mem_before > 0 else -1,
            'mem_total_mb': round(mem_after, 1),
            'components': perf.summary(),
            'ok': True,
        }
    except Exception as e:
        wall = time.perf_counter() - t0
        return {'wall_s': round(wall, 4), 'error': str(e)[:200], 'ok': False}


def spatial_run(n_nodes, K_per_node, years, disease_year=None, seed=42, movement=False):
    """Run spatial simulation with variable node count."""
    cfg = default_config()
    cfg.spawning = None  # spatial sim uses pulse reproduction
    cfg.movement.enabled = movement
    
    # Build a network with n_nodes, each at K_per_node capacity
    from sswd_evoepi.spatial import NodeDefinition, MetapopulationNetwork
    
    nodes_def = []
    for i in range(n_nodes):
        nodes_def.append(NodeDefinition(
            node_id=i,
            name=f"Node_{i}",
            latitude=48.0 + i * 0.5,  # spread along coast
            longitude=-123.0 - i * 0.3,
            carrying_capacity=K_per_node,
            habitat_area=max(10000, K_per_node * 200),
            mean_sst=12.0 + i * 0.5,
            sst_amplitude=4.0,
            salinity=30.0,
            flushing_rate=0.02,
            is_fjord=(i == 1),
        ))
    
    # Build connectivity matrices (exponential kernel)
    positions = np.array([(nd.latitude, nd.longitude) for nd in nodes_def])
    # Simple distance-based C and D matrices
    from scipy.spatial.distance import cdist
    dists = cdist(positions, positions) * 111.0  # approx km
    
    # Connectivity: C_ij = exp(-d_ij / scale) (larval transport)
    scale_c = 200.0  # km
    C = np.exp(-dists / scale_c)
    np.fill_diagonal(C, 0)
    # Normalize rows
    row_sums = C.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    C = C / row_sums
    
    # Pathogen dispersal: D_ij = smaller scale
    scale_d = 50.0
    D = np.exp(-dists / scale_d) * 0.01
    np.fill_diagonal(D, 0)
    
    network = MetapopulationNetwork(
        nodes=[],
        C=C,
        D=D,
        n_nodes=n_nodes,
    )
    
    # Build using the spatial module's infrastructure
    from sswd_evoepi.spatial import SpatialNode
    for nd in nodes_def:
        sn = SpatialNode(definition=nd)
        network.nodes.append(sn)
    
    mem_before = get_memory_mb()
    t0 = time.perf_counter()
    try:
        r = run_spatial_simulation(
            network=network,
            n_years=years,
            disease_year=disease_year,
            seed=seed,
            config=cfg,
        )
        wall = time.perf_counter() - t0
        mem_after = get_memory_mb()
        
        return {
            'wall_s': round(wall, 4),
            'final_pop': r.final_total_pop,
            'total_agents': n_nodes * K_per_node,
            'mem_mb': round(mem_after - mem_before, 1) if mem_before > 0 else -1,
            'mem_total_mb': round(mem_after, 1),
            'ok': True,
        }
    except Exception as e:
        wall = time.perf_counter() - t0
        return {
            'wall_s': round(wall, 4),
            'error': str(e)[:200],
            'ok': False,
            'traceback': traceback.format_exc()[-500:],
        }


def parallel_worker(args):
    """Worker for parallel replicate testing."""
    n, years, seed = args
    return single_node_run(n, years, spawning=True, seed=seed)


# ═══════════════════════════════════════════════════════════════════
# MAIN MEASUREMENT SEQUENCE
# ═══════════════════════════════════════════════════════════════════

def main():
    all_results = {}
    
    # ── AXIS 1: Population size (extended) ────────────────────────
    print("\n" + "="*60)
    print(" AXIS 1: Population Size (N) — 5yr, spawning enabled")
    print("="*60)
    print(f"{'N':>8} {'Time':>10} {'Pop':>8} {'Mem(MB)':>8} {'Top':>20}")
    
    axis1 = []
    for n in [25, 50, 100, 200, 500, 1000, 2000, 5000, 10000]:
        r = single_node_run(n, years=5, spawning=True)
        top = 'N/A'
        if r['ok'] and 'components' in r:
            comps = r['components']
            top_item = max(
                [(k, v) for k, v in comps.items() if k != '_total_s'],
                key=lambda x: x[1].get('total_s', 0), default=('?', {'pct': 0})
            )
            top = f"{top_item[0]} ({top_item[1]['pct']:.0f}%)"
        status = f"{r.get('final_pop', 'ERR'):>8}" if r['ok'] else "    ERR"
        print(f"{n:>8} {r['wall_s']:>9.3f}s {status} {r.get('mem_total_mb', -1):>7.0f} {top:>20}")
        axis1.append({'n': n, **r})
        sys.stdout.flush()
    all_results['axis_N'] = axis1
    
    # ── AXIS 2: Duration (extended) ───────────────────────────────
    print("\n" + "="*60)
    print(" AXIS 2: Duration (T years) — 200 agents, spawning enabled")
    print("="*60)
    print(f"{'T':>8} {'Time':>10} {'Pop':>8} {'Mem(MB)':>8}")
    
    axis2 = []
    for t in [1, 2, 5, 10, 20, 50, 100, 200]:
        r = single_node_run(200, years=t, spawning=True, timeout_s=300)
        status = f"{r.get('final_pop', 'ERR'):>8}" if r['ok'] else "    ERR"
        print(f"{t:>8} {r['wall_s']:>9.3f}s {status} {r.get('mem_total_mb', -1):>7.0f}")
        axis2.append({'years': t, **r})
        sys.stdout.flush()
    all_results['axis_T'] = axis2
    
    # ── AXIS 3: Spatial nodes (K) ─────────────────────────────────
    print("\n" + "="*60)
    print(" AXIS 3: Spatial Nodes (K) — 200 agents/node, 5yr, no disease")
    print("="*60)
    print(f"{'K':>8} {'Total N':>8} {'Time':>10} {'Pop':>8} {'Mem(MB)':>8}")
    
    axis3 = []
    for k in [2, 3, 5, 10, 20, 50, 100, 150]:
        r = spatial_run(n_nodes=k, K_per_node=200, years=5)
        total = k * 200
        status = f"{r.get('final_pop', 'ERR'):>8}" if r['ok'] else "    ERR"
        print(f"{k:>8} {total:>8} {r['wall_s']:>9.3f}s {status} {r.get('mem_total_mb', -1):>7.0f}")
        axis3.append({'k': k, 'total_agents': total, **r})
        sys.stdout.flush()
        if not r['ok']:
            print(f"         Error: {r.get('error', 'unknown')[:80]}")
            break
    all_results['axis_K'] = axis3
    
    # ── AXIS 4: N×T interaction matrix ────────────────────────────
    print("\n" + "="*60)
    print(" AXIS 4: N×T Interaction Matrix (spawning enabled)")
    print("="*60)
    ns_matrix = [50, 200, 500, 1000, 2000]
    ts_matrix = [5, 10, 20, 50]
    
    print(f"{'':>8}", end='')
    for t in ts_matrix:
        print(f"{'T='+str(t):>10}", end='')
    print()
    
    axis4 = []
    for n in ns_matrix:
        print(f"{'N='+str(n):>8}", end='')
        for t in ts_matrix:
            r = single_node_run(n, years=t, spawning=True, timeout_s=120)
            val = f"{r['wall_s']:.2f}s" if r['ok'] else "ERR"
            print(f"{val:>10}", end='')
            axis4.append({'n': n, 'years': t, **r})
            sys.stdout.flush()
        print()
    all_results['axis_NxT'] = axis4
    
    # ── AXIS 5: Spawning overhead at scale ────────────────────────
    print("\n" + "="*60)
    print(" AXIS 5: Spawning Overhead (10yr)")
    print("="*60)
    print(f"{'N':>8} {'No Spawn':>10} {'Spawn':>10} {'Overhead':>10}")
    
    axis5 = []
    for n in [50, 200, 500, 1000, 2000, 5000]:
        r_off = single_node_run(n, years=10, spawning=False)
        r_on = single_node_run(n, years=10, spawning=True)
        oh = r_on['wall_s'] / r_off['wall_s'] if r_off['wall_s'] > 0 else -1
        print(f"{n:>8} {r_off['wall_s']:>9.3f}s {r_on['wall_s']:>9.3f}s {oh:>9.1f}×")
        axis5.append({'n': n, 'no_spawn_s': r_off['wall_s'], 'spawn_s': r_on['wall_s'], 'overhead': round(oh, 2)})
        sys.stdout.flush()
    all_results['axis_spawning'] = axis5
    
    # ── AXIS 6: K×N combined (total system stress) ────────────────
    print("\n" + "="*60)
    print(" AXIS 6: K×N Combined Stress Test (5yr, no disease)")
    print("="*60)
    print(f"{'K':>5} {'N/node':>8} {'Total':>8} {'Time':>10} {'Mem(MB)':>8}")
    
    axis6 = []
    combos = [
        (5, 200), (5, 500), (5, 1000),
        (10, 200), (10, 500),
        (20, 200), (20, 500),
        (50, 200), (50, 500),
        (100, 200), (100, 500),
        (150, 200), (150, 500),
    ]
    for k, n in combos:
        total = k * n
        r = spatial_run(n_nodes=k, K_per_node=n, years=5)
        status = f"{r.get('final_pop', 'ERR'):>8}" if r['ok'] else "    ERR"
        print(f"{k:>5} {n:>8} {total:>8} {r['wall_s']:>9.3f}s {r.get('mem_total_mb', -1):>7.0f}")
        axis6.append({'k': k, 'n_per_node': n, 'total': total, **r})
        sys.stdout.flush()
        if not r['ok']:
            print(f"      Error: {r.get('error', 'unknown')[:80]}")
            if r.get('mem_total_mb', 0) > 25000:  # approaching 31GB limit
                print("      Approaching memory limit, stopping K×N sweep")
                break
    all_results['axis_KxN'] = axis6
    
    # ── AXIS 7: Parallel replicates ───────────────────────────────
    print("\n" + "="*60)
    print(" AXIS 7: Parallel Replicates (200 agents, 10yr)")
    print("="*60)
    print(f"{'Procs':>8} {'Total':>8} {'Wall':>10} {'Throughput':>12}")
    
    axis7 = []
    base_args = [(200, 10, seed) for seed in range(16)]
    
    for n_proc in [1, 2, 4, 8, 16]:
        t0 = time.perf_counter()
        with mp.Pool(n_proc) as pool:
            results = pool.map(parallel_worker, base_args[:n_proc * 2])
        wall = time.perf_counter() - t0
        n_runs = n_proc * 2
        throughput = n_runs / wall
        print(f"{n_proc:>8} {n_runs:>8} {wall:>9.2f}s {throughput:>10.1f} runs/s")
        axis7.append({'n_procs': n_proc, 'n_runs': n_runs, 'wall_s': round(wall, 3), 
                       'throughput': round(throughput, 2)})
        sys.stdout.flush()
    all_results['axis_parallel'] = axis7
    
    # ── AXIS 8: Disease impact at scale ───────────────────────────
    print("\n" + "="*60)
    print(" AXIS 8: Disease Impact — Spatial (5 nodes, 10yr)")
    print("="*60)
    print(f"{'K/node':>8} {'No Dis':>10} {'Disease':>10} {'Ratio':>8}")
    
    axis8 = []
    for K in [200, 500, 1000]:
        r_off = spatial_run(5, K, years=10, disease_year=None)
        r_on = spatial_run(5, K, years=10, disease_year=3)
        ratio = r_on['wall_s'] / r_off['wall_s'] if r_off['wall_s'] > 0 else -1
        print(f"{K:>8} {r_off['wall_s']:>9.3f}s {r_on['wall_s']:>9.3f}s {ratio:>7.2f}×")
        axis8.append({'K_per_node': K, 'no_disease_s': r_off['wall_s'], 
                       'disease_s': r_on['wall_s'], 'ratio': round(ratio, 2)})
        sys.stdout.flush()
    all_results['axis_disease_spatial'] = axis8
    
    # ── Save all results ──────────────────────────────────────────
    os.makedirs('results/performance', exist_ok=True)
    with open('results/performance/exhaustive_scaling_data.json', 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    
    # ── Summary ───────────────────────────────────────────────────
    print("\n" + "="*60)
    print(" EXHAUSTIVE SCALING SUMMARY")
    print("="*60)
    
    # Fit power laws
    n_data = [(r['n'], r['wall_s']) for r in axis1 if r.get('ok') and r['wall_s'] > 0.01]
    if len(n_data) >= 3:
        ns, ts = zip(*n_data)
        slope_n = np.polyfit(np.log(ns), np.log(ts), 1)[0]
        print(f"  N scaling (single node): O(N^{slope_n:.2f})")
    
    t_data = [(r['years'], r['wall_s']) for r in axis2 if r.get('ok') and r['wall_s'] > 0.01]
    if len(t_data) >= 3:
        ys, ts = zip(*t_data)
        slope_t = np.polyfit(np.log(ys), np.log(ts), 1)[0]
        print(f"  T scaling (single node): O(T^{slope_t:.2f})")
    
    k_data = [(r['k'], r['wall_s']) for r in axis3 if r.get('ok') and r['wall_s'] > 0.01]
    if len(k_data) >= 3:
        ks, ts = zip(*k_data)
        slope_k = np.polyfit(np.log(ks), np.log(ts), 1)[0]
        print(f"  K scaling (spatial): O(K^{slope_k:.2f})")
    
    # Machine limits
    max_single = max((r['n'] for r in axis1 if r.get('ok')), default=0)
    max_spatial_total = max((r['total'] for r in axis6 if r.get('ok')), default=0)
    print(f"\n  Machine limits:")
    print(f"    Max single-node N tested: {max_single}")
    print(f"    Max spatial total agents: {max_spatial_total}")
    
    if axis7:
        best_tp = max(r['throughput'] for r in axis7)
        best_procs = next(r['n_procs'] for r in axis7 if r['throughput'] == best_tp)
        print(f"    Best parallel throughput: {best_tp:.1f} runs/s ({best_procs} processes)")
    
    print(f"\n  Data saved: results/performance/exhaustive_scaling_data.json")


if __name__ == '__main__':
    main()
