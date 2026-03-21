#!/usr/bin/env python3
"""Generate spawning dynamics visualizations + animated GIFs (Phase 13).

Produces:
  PART A (static PNGs):
    1. spawning_event_profile.png
    2. spawning_participation.png
    3. readiness_cascade.png
    4. spawning_density_scatter.png
    5. spawning_before_after.png

  PART B (animated GIFs):
    1. spawning_mass_event_highdensity.gif
    2. spawning_mass_event_lowdensity.gif
    3. spawning_param_comparison.gif
    4. spawning_full_season.gif

Usage:
    python3 scripts/generate_spawning_viz.py

Author: Anton 🔬
"""

import sys
import os
import time
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

import numpy as np

from sswd_evoepi.viz.spawning_viz import (
    run_simulation_with_snapshots,
    plot_spawning_event_profile,
    plot_spawning_participation,
    plot_readiness_cascade,
    plot_spawning_vs_density,
    plot_spawning_before_after,
    create_spawning_animation,
    create_comparison_animation,
)

# Output directories
OUT_DIR = project_root / 'results' / 'continuous_settlement' / 'viz' / 'spawning'
ANIM_DIR = OUT_DIR / 'animations'
OUT_DIR.mkdir(parents=True, exist_ok=True)
ANIM_DIR.mkdir(parents=True, exist_ok=True)


def timer(label):
    """Simple timer context manager."""
    class Timer:
        def __enter__(self):
            self.start = time.time()
            print(f"\n{'='*60}")
            print(f"  {label}")
            print(f"{'='*60}")
            return self
        def __exit__(self, *args):
            elapsed = time.time() - self.start
            print(f"  ✓ {label} — {elapsed:.1f}s")
    return Timer()


def main():
    print("=" * 70)
    print("  SSWD-EvoEpi Phase 13: Spawning Dynamics Visualization")
    print("=" * 70)
    t0 = time.time()

    hab_side = np.sqrt(10000.0)  # 100m side

    # ── SIM 1: High density, new params (year 1, pre-disease) ──────
    # Full season for event profile, readiness cascade, and full-season GIF
    with timer("Sim 1: High density, new params, full season (year 1)"):
        snaps_high_full, result_high = run_simulation_with_snapshots(
            n_individuals=500,
            carrying_capacity=500,
            habitat_area=10000.0,
            T_celsius=10.5,
            n_years=3,
            seed=42,
            snapshot_start_doy=None,  # entire season
            snapshot_end_doy=None,
            snapshot_year=1,  # year 1 (pre-disease, high density)
            disease_year=None,  # no disease
        )
        print(f"  Captured {len(snaps_high_full)} daily snapshots")

    # ── SIM 2: Low density, new params (simulate post-crash) ──────
    with timer("Sim 2: Low density (post-crash), new params, full season"):
        snaps_low_full, result_low = run_simulation_with_snapshots(
            n_individuals=50,
            carrying_capacity=500,
            habitat_area=10000.0,
            T_celsius=10.5,
            n_years=3,
            seed=42,
            snapshot_start_doy=None,
            snapshot_end_doy=None,
            snapshot_year=1,
            disease_year=None,
        )
        print(f"  Captured {len(snaps_low_full)} daily snapshots")

    # ── SIM 3: High density, OLD params ───────────────────────────
    with timer("Sim 3: High density, OLD params (pre-overhaul), full season"):
        old_overrides = {
            'female_max_bouts': 1,       # single spawn
            'male_refractory_days': 21,   # 21-day refractory
            'induction_male_to_female': 0.30,  # weaker induction
            'induction_female_to_male': 0.50,  # weaker induction
            'readiness_induction_prob': 0.0,   # no readiness induction
            'gravity_enabled': False,          # no gravity
        }
        snaps_old_full, result_old = run_simulation_with_snapshots(
            n_individuals=500,
            carrying_capacity=500,
            habitat_area=10000.0,
            T_celsius=10.5,
            n_years=3,
            seed=42,  # same seed for comparison
            snapshot_start_doy=None,
            snapshot_end_doy=None,
            snapshot_year=1,
            disease_year=None,
            spawning_overrides=old_overrides,
        )
        print(f"  Captured {len(snaps_old_full)} daily snapshots")

    # ── Find peak spawning window for focused GIFs ────────────────
    if snaps_high_full:
        peak_idx = max(range(len(snaps_high_full)),
                       key=lambda i: snaps_high_full[i].n_spawning_today)
        peak_doy = snaps_high_full[peak_idx].doy
        print(f"\n  Peak spawning at DOY {peak_doy} "
              f"({snaps_high_full[peak_idx].n_spawning_today} spawners)")

        # Extract ±30 day window around peak
        window_start = peak_doy - 30
        window_end = peak_doy + 30

        snaps_high_window = [s for s in snaps_high_full
                             if _doy_distance(s.doy, peak_doy) <= 30]
        snaps_low_window = [s for s in snaps_low_full
                            if _doy_distance(s.doy, peak_doy) <= 30]
        snaps_old_window = [s for s in snaps_old_full
                            if _doy_distance(s.doy, peak_doy) <= 30]
    else:
        snaps_high_window = snaps_high_full
        snaps_low_window = snaps_low_full
        snaps_old_window = snaps_old_full

    # ════════════════════════════════════════════════════════════════
    # PART A: STATIC VISUALIZATIONS
    # ════════════════════════════════════════════════════════════════

    # A1. Mass Spawning Event Profile
    with timer("A1: Spawning Event Profile"):
        plot_spawning_event_profile(
            snaps_high_full,
            save_path=str(OUT_DIR / 'spawning_event_profile.png'),
        )

    # A2. Spawning Participation Rate
    with timer("A2: Spawning Participation Rate"):
        # Collect end-of-season snapshots per year from additional sims
        # For simplicity, use the full-season data we have
        snaps_by_year = {}

        # Year 1 from high-density sim
        snaps_by_year[1] = snaps_high_full

        # Run year 2 (with disease) for comparison
        snaps_disease_yr2, _ = run_simulation_with_snapshots(
            n_individuals=500,
            carrying_capacity=500,
            habitat_area=10000.0,
            T_celsius=10.5,
            n_years=4,
            seed=42,
            snapshot_start_doy=None,
            snapshot_end_doy=None,
            snapshot_year=2,
            disease_year=2,  # disease hits in year 2
        )
        snaps_by_year[2] = snaps_disease_yr2

        # Year 0 (early, population still equilibrating)
        snaps_yr0, _ = run_simulation_with_snapshots(
            n_individuals=500,
            carrying_capacity=500,
            habitat_area=10000.0,
            T_celsius=10.5,
            n_years=2,
            seed=42,
            snapshot_start_doy=None,
            snapshot_end_doy=None,
            snapshot_year=0,
            disease_year=None,
        )
        snaps_by_year[0] = snaps_yr0

        plot_spawning_participation(
            snaps_by_year,
            n_individuals=500,
            save_path=str(OUT_DIR / 'spawning_participation.png'),
        )

    # A3. Readiness Induction Cascade Diagram
    with timer("A3: Readiness Cascade Diagram"):
        plot_readiness_cascade(
            snaps_high_full,
            save_path=str(OUT_DIR / 'readiness_cascade.png'),
        )

    # A4. Spawning Intensity vs Density Scatter
    with timer("A4: Spawning Density Scatter"):
        plot_spawning_vs_density(
            [
                ('High density (500/K=500)', snaps_high_full, 500),
                ('Low density (50/K=500)', snaps_low_full, 500),
            ],
            save_path=str(OUT_DIR / 'spawning_density_scatter.png'),
        )

    # A5. Before/After Spawning Overhaul
    with timer("A5: Before/After Comparison"):
        plot_spawning_before_after(
            snaps_old_full,
            snaps_high_full,
            save_path=str(OUT_DIR / 'spawning_before_after.png'),
        )

    # ════════════════════════════════════════════════════════════════
    # PART B: ANIMATED GIFS
    # ════════════════════════════════════════════════════════════════

    # GIF 1: Mass spawning event at high density (±30 days around peak)
    with timer("GIF 1: Mass spawning event — high density"):
        create_spawning_animation(
            snaps_high_window,
            hab_side=hab_side,
            fps=4,
            title_prefix='HIGH DENSITY  |  ',
            save_path=str(ANIM_DIR / 'spawning_mass_event_highdensity.gif'),
            dpi=80,
        )

    # GIF 2: Same window at low density
    with timer("GIF 2: Mass spawning event — low density"):
        create_spawning_animation(
            snaps_low_window,
            hab_side=hab_side,
            fps=4,
            title_prefix='LOW DENSITY  |  ',
            save_path=str(ANIM_DIR / 'spawning_mass_event_lowdensity.gif'),
            dpi=80,
        )

    # GIF 3: Parameter comparison (old vs new, windowed)
    with timer("GIF 3: Parameter comparison animation"):
        create_comparison_animation(
            snaps_old_window,
            snaps_high_window,
            hab_side=hab_side,
            fps=4,
            label_a='Old Params',
            label_b='Phase 12 Params',
            max_female_bouts_a=1,
            max_male_bouts_a=3,
            max_female_bouts_b=2,
            max_male_bouts_b=3,
            save_path=str(ANIM_DIR / 'spawning_param_comparison.gif'),
            dpi=80,
        )

    # GIF 4: Full season overview at 8 fps
    with timer("GIF 4: Full season overview"):
        create_spawning_animation(
            snaps_high_full,
            hab_side=hab_side,
            fps=8,
            title_prefix='FULL SEASON  |  ',
            save_path=str(ANIM_DIR / 'spawning_full_season.gif'),
            dpi=80,
        )

    # ── Summary ───────────────────────────────────────────────────
    total_time = time.time() - t0
    print(f"\n{'='*70}")
    print(f"  Phase 13 COMPLETE — {total_time:.1f}s total")
    print(f"{'='*70}")

    # List outputs
    print(f"\n  Static PNGs:")
    for f in sorted(OUT_DIR.glob('*.png')):
        size_kb = f.stat().st_size / 1024
        print(f"    {f.name} ({size_kb:.0f} KB)")

    print(f"\n  Animated GIFs:")
    for f in sorted(ANIM_DIR.glob('*.gif')):
        size_mb = f.stat().st_size / (1024 * 1024)
        print(f"    {f.name} ({size_mb:.1f} MB)")

    print(f"\n  Output directory: {OUT_DIR}")
    print(f"  Animation directory: {ANIM_DIR}")


def _doy_distance(doy1: int, doy2: int) -> int:
    """Circular distance between two days of year."""
    direct = abs(doy1 - doy2)
    wrapped = 365 - direct
    return min(direct, wrapped)


if __name__ == '__main__':
    main()
