#!/usr/bin/env python3
"""Quick debug run: 15 years, BH fix verified, individual snapshots + wildfire.

- Years 0-4: spinup (should be STABLE now with BH fix)
- Year 5: disease
- Years 5-14: epidemic
- Individual snapshots recorded daily during epidemic for wildfire viz
"""

import sys
import time
from pathlib import Path

import numpy as np

project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

from sswd_evoepi.config import default_config
from sswd_evoepi.model import run_spatial_simulation, make_effect_sizes
from sswd_evoepi.snapshots import SnapshotRecorder
from sswd_evoepi.spatial import make_5node_network

N_YEARS = 15
DISEASE_YEAR = 5
SEED = 42
OUTPUT_DIR = project_root / "results" / "debug_bh_fix"


def main():
    print("=" * 60)
    print("DEBUG RUN: BH fix + individual snapshots + wildfire")
    print("=" * 60)

    config = default_config()
    network = make_5node_network(seed=SEED)

    for node in network.nodes:
        nd = node.definition
        print(f"  [{nd.node_id}] {nd.name}: K={nd.carrying_capacity}")

    # Snapshot recorder: daily during epidemic (years 4-14 = days 1460-5110)
    snap = SnapshotRecorder(
        enabled=True,
        interval_days=7,  # weekly to keep memory reasonable
        start_day=DISEASE_YEAR * 365 - 365,  # 1 year before disease
        end_day=N_YEARS * 365,
    )

    t0 = time.time()

    def progress(year, total):
        elapsed = time.time() - t0
        print(f"  Year {year:3d}/{total}  [{elapsed:.0f}s]", flush=True)

    result = run_spatial_simulation(
        network=network,
        n_years=N_YEARS,
        disease_year=DISEASE_YEAR,
        initial_infected_per_node=5,
        seed=SEED,
        config=config,
        progress_callback=progress,
        snapshot_recorder=snap,
    )

    elapsed = time.time() - t0
    print(f"\nDone in {elapsed:.1f}s")

    # ── Check pre-epidemic stability ──────────────────────────
    print("\n─── Pre-Epidemic Population Check ───")
    names = result.node_names
    K_vals = result.node_K
    for i in range(result.n_nodes):
        pops = result.yearly_pop[i, :DISEASE_YEAR]
        print(f"  {names[i]:25s}  K={K_vals[i]:4d}  "
              f"pop=[{', '.join(str(int(p)) for p in pops)}]")

    print(f"\n  Total: {[int(result.yearly_total_pop[y]) for y in range(DISEASE_YEAR)]}")

    # Check stability: is year 4 pop within 10% of year 0?
    y0 = int(result.yearly_total_pop[0])
    y4 = int(result.yearly_total_pop[DISEASE_YEAR - 1])
    drift = abs(y4 - y0) / y0 * 100
    stable = drift < 10
    print(f"  Drift from year 0→{DISEASE_YEAR-1}: {drift:.1f}% {'✓ STABLE' if stable else '✗ UNSTABLE'}")

    # Post-epidemic
    print("\n─── Year-by-Year Total Population ───")
    for yr in range(N_YEARS):
        marker = " ← disease" if yr == DISEASE_YEAR else ""
        print(f"  Year {yr:3d}: {result.yearly_total_pop[yr]:6d}{marker}")

    # Save results
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        OUTPUT_DIR / "simulation_data.npz",
        yearly_pop=result.yearly_pop,
        yearly_total_pop=result.yearly_total_pop,
        yearly_disease_deaths=result.yearly_disease_deaths,
        yearly_mean_resistance=result.yearly_mean_resistance,
        node_K=result.node_K,
    )

    # Save snapshots
    snap_path = OUTPUT_DIR / "snapshots.npz"
    snap.save(str(snap_path))
    n_snaps = len(snap.snapshots)
    mem_mb = snap.memory_estimate_mb()
    print(f"\nSnapshots: {n_snaps} captured, {mem_mb:.1f} MB")

    # ── Wildfire visualization ────────────────────────────────
    print("\nGenerating wildfire animation...")
    try:
        from viz.wildfire import render_wildfire, render_wildfire_static

        # Compute habitat side lengths from areas
        hab_sides = {}
        for node in network.nodes:
            nd = node.definition
            hab_sides[nd.node_id] = np.sqrt(nd.habitat_area)

        render_wildfire(
            snap,
            output_path=str(OUTPUT_DIR / "wildfire.gif"),
            node_names=names,
            fps=8,
            dpi=80,
            max_frames=200,
            habitat_sides=hab_sides,
        )

        render_wildfire_static(
            snap,
            output_dir=str(OUTPUT_DIR / "wildfire_stills"),
            node_names=names,
            habitat_sides=hab_sides,
        )
    except Exception as e:
        print(f"Wildfire viz error: {e}")
        import traceback
        traceback.print_exc()

    print(f"\nAll output in {OUTPUT_DIR}/")
    return 0 if stable else 1


if __name__ == "__main__":
    sys.exit(main())
