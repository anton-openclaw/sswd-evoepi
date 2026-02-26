#!/usr/bin/env python3
"""Benchmark Numba JIT movement kernel vs pure NumPy baseline.

Compares wall-clock time for the movement substep loop under:
  1. Pure NumPy (existing vectorized code)
  2. Numba JIT (serial, @njit)
  3. Numba JIT + parallel (@njit(parallel=True) with prange)

Usage:
    python scripts/benchmark_numba.py [--agents N] [--substeps S] [--iters I]
"""

from __future__ import annotations

import argparse
import sys
import time
import numpy as np
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))


# ═══════════════════════════════════════════════════════════════════════
# Pure NumPy baseline (copied from movement.py, no imports needed)
# ═══════════════════════════════════════════════════════════════════════

TWO_PI = 2.0 * np.pi


def _reflect_numpy(x: np.ndarray, side: float) -> np.ndarray:
    period = 2.0 * side
    x = x % period
    x = np.where(x > side, period - x, x)
    return x


def movement_numpy(headings, x, y, speed_dt, all_turns, habitat_side, substeps):
    """Pure NumPy movement (baseline)."""
    for s in range(substeps):
        headings = (headings + all_turns[s]) % TWO_PI
        dx = speed_dt * np.cos(headings)
        dy = speed_dt * np.sin(headings)
        x = _reflect_numpy(x + dx, habitat_side)
        y = _reflect_numpy(y + dy, habitat_side)
    return headings, x, y


# ═══════════════════════════════════════════════════════════════════════
# Numba JIT kernels
# ═══════════════════════════════════════════════════════════════════════

try:
    import numba
    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False

if HAS_NUMBA:
    @numba.njit(cache=True)
    def movement_numba_serial(headings, x, y, speed_dt, all_turns,
                              habitat_side, n_substeps):
        TWO_PI = 2.0 * np.pi
        period = 2.0 * habitat_side
        n = len(headings)
        for s in range(n_substeps):
            for i in range(n):
                headings[i] = (headings[i] + all_turns[s, i]) % TWO_PI
                dx = speed_dt[i] * np.cos(headings[i])
                dy = speed_dt[i] * np.sin(headings[i])
                new_x = (x[i] + dx) % period
                if new_x < 0.0:
                    new_x += period
                if new_x > habitat_side:
                    new_x = period - new_x
                x[i] = new_x
                new_y = (y[i] + dy) % period
                if new_y < 0.0:
                    new_y += period
                if new_y > habitat_side:
                    new_y = period - new_y
                y[i] = new_y
        return headings, x, y

    @numba.njit(parallel=True, cache=True)
    def movement_numba_parallel(headings, x, y, speed_dt, all_turns,
                                habitat_side, n_substeps):
        TWO_PI = 2.0 * np.pi
        period = 2.0 * habitat_side
        n = len(headings)
        for s in range(n_substeps):
            for i in numba.prange(n):
                headings[i] = (headings[i] + all_turns[s, i]) % TWO_PI
                dx = speed_dt[i] * np.cos(headings[i])
                dy = speed_dt[i] * np.sin(headings[i])
                new_x = (x[i] + dx) % period
                if new_x < 0.0:
                    new_x += period
                if new_x > habitat_side:
                    new_x = period - new_x
                x[i] = new_x
                new_y = (y[i] + dy) % period
                if new_y < 0.0:
                    new_y += period
                if new_y > habitat_side:
                    new_y = period - new_y
                y[i] = new_y
        return headings, x, y


# ═══════════════════════════════════════════════════════════════════════
# Benchmark harness
# ═══════════════════════════════════════════════════════════════════════

def make_test_data(n_agents: int, substeps: int, rng: np.random.Generator):
    """Generate synthetic test data matching real simulation layout."""
    headings = rng.uniform(0, TWO_PI, n_agents).astype(np.float64)
    x = rng.uniform(0, 500.0, n_agents).astype(np.float64)
    y = rng.uniform(0, 500.0, n_agents).astype(np.float64)
    speed_dt = rng.uniform(1.0, 60.0, n_agents).astype(np.float64)
    all_turns = rng.normal(0, 0.5, (substeps, n_agents)).astype(np.float64)
    habitat_side = 500.0
    return headings, x, y, speed_dt, all_turns, habitat_side


def benchmark_one(fn, label, headings, x, y, speed_dt, all_turns,
                  habitat_side, substeps, n_iters):
    """Benchmark a single movement function."""
    times = []
    for _ in range(n_iters):
        # Copy inputs (functions modify in-place for Numba)
        h = headings.copy()
        xc = x.copy()
        yc = y.copy()
        t0 = time.perf_counter()
        fn(h, xc, yc, speed_dt, all_turns, habitat_side, substeps)
        t1 = time.perf_counter()
        times.append(t1 - t0)
    times = np.array(times)
    return {
        'label': label,
        'mean': times.mean(),
        'std': times.std(),
        'min': times.min(),
        'median': np.median(times),
        'times': times,
    }


def run_benchmarks(n_agents: int, substeps: int, n_iters: int, n_nodes: int):
    """Run all benchmark variants."""
    rng = np.random.default_rng(42)

    print(f"╔══════════════════════════════════════════════════════════════╗")
    print(f"║  SSWD-EvoEpi Movement Kernel Benchmark                     ║")
    print(f"╠══════════════════════════════════════════════════════════════╣")
    print(f"║  Agents/node: {n_agents:<6d}  Substeps: {substeps:<4d}  Iterations: {n_iters:<4d} ║")
    print(f"║  Simulated nodes: {n_nodes:<6d}                                    ║")
    print(f"╚══════════════════════════════════════════════════════════════╝")
    print()

    headings, x, y, speed_dt, all_turns, habitat_side = make_test_data(
        n_agents, substeps, rng
    )

    results = []

    # 1. Pure NumPy baseline
    res = benchmark_one(
        movement_numpy, "Pure NumPy",
        headings, x, y, speed_dt, all_turns, habitat_side, substeps, n_iters,
    )
    results.append(res)

    if HAS_NUMBA:
        # 2. Numba JIT (warmup first)
        print("  Warming up Numba JIT (serial)...", end=" ", flush=True)
        h, xc, yc = headings.copy(), x.copy(), y.copy()
        t0 = time.perf_counter()
        movement_numba_serial(h, xc, yc, speed_dt.copy(), all_turns.copy(),
                              habitat_side, substeps)
        jit_warmup = time.perf_counter() - t0
        print(f"done ({jit_warmup:.2f}s JIT compile)")

        res = benchmark_one(
            movement_numba_serial, "Numba JIT (serial)",
            headings, x, y, speed_dt, all_turns, habitat_side, substeps, n_iters,
        )
        results.append(res)

        # 3. Numba JIT parallel (warmup first)
        print("  Warming up Numba JIT (parallel)...", end=" ", flush=True)
        h, xc, yc = headings.copy(), x.copy(), y.copy()
        t0 = time.perf_counter()
        movement_numba_parallel(h, xc, yc, speed_dt.copy(), all_turns.copy(),
                                habitat_side, substeps)
        jit_warmup_p = time.perf_counter() - t0
        print(f"done ({jit_warmup_p:.2f}s JIT compile)")

        res = benchmark_one(
            movement_numba_parallel, "Numba JIT (parallel)",
            headings, x, y, speed_dt, all_turns, habitat_side, substeps, n_iters,
        )
        results.append(res)
    else:
        print("  ⚠ Numba not available, skipping JIT benchmarks")

    print()

    # ── Results table ────────────────────────────────────────────────
    baseline = results[0]['median']
    print(f"{'Method':<25s} {'Median (ms)':>12s} {'Min (ms)':>10s} {'Speedup':>8s}")
    print("─" * 60)
    for res in results:
        speedup = baseline / res['median'] if res['median'] > 0 else float('inf')
        print(f"  {res['label']:<23s} {res['median']*1000:>10.3f}   {res['min']*1000:>8.3f}   {speedup:>6.2f}×")

    # ── Per-day estimate ──────────────────────────────────────────────
    print()
    print("── Per-day wall time estimate (all nodes) ──")
    for res in results:
        per_day_s = res['median'] * n_nodes
        print(f"  {res['label']:<23s}  {per_day_s:.4f}s/day  "
              f"({per_day_s * 4745:.1f}s for 13-year run)")

    # ── Numerical agreement check ────────────────────────────────────
    if HAS_NUMBA and len(results) >= 2:
        print()
        print("── Numerical agreement ──")
        # NumPy returns new arrays; Numba modifies in-place
        h_np, x_np, y_np = movement_numpy(
            headings.copy(), x.copy(), y.copy(),
            speed_dt, all_turns, habitat_side, substeps,
        )

        h_jit, x_jit, y_jit = headings.copy(), x.copy(), y.copy()
        movement_numba_serial(h_jit, x_jit, y_jit, speed_dt.copy(),
                              all_turns.copy(), habitat_side, substeps)

        dh = np.max(np.abs(h_np - h_jit))
        dx = np.max(np.abs(x_np - x_jit))
        dy = np.max(np.abs(y_np - y_jit))
        print(f"  Max |Δheading|: {dh:.2e}")
        print(f"  Max |Δx|:       {dx:.2e}")
        print(f"  Max |Δy|:       {dy:.2e}")
        ok = dh < 1e-10 and dx < 1e-10 and dy < 1e-10
        print(f"  Agreement: {'✓ PASS' if ok else '✗ FAIL (float precision differences)'}")

    return results


def main():
    parser = argparse.ArgumentParser(description="Benchmark Numba movement kernel")
    parser.add_argument("--agents", type=int, default=5000,
                        help="Agents per node (default: 5000)")
    parser.add_argument("--substeps", type=int, default=24,
                        help="Substeps per day (default: 24)")
    parser.add_argument("--iters", type=int, default=50,
                        help="Benchmark iterations (default: 50)")
    parser.add_argument("--nodes", type=int, default=896,
                        help="Simulated node count for per-day estimates (default: 896)")
    args = parser.parse_args()

    run_benchmarks(args.agents, args.substeps, args.iters, args.nodes)


if __name__ == "__main__":
    main()
