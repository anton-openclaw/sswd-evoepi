#!/usr/bin/env python3
"""Profile the spatial simulation to find where time is spent.

Runs a 3-node 20yr spatial sim (Sitka/Howe Sound/Monterey) with both
cProfile and manual per-component timers via monkey-patching.

The approach: wrap key functions called inside run_spatial_simulation's
daily loop with timing decorators that accumulate into a global dict.
This avoids modifying model.py while giving us per-component breakdowns.

Output: results/continuous_mortality/profile_output.txt
"""

import cProfile
import io
import os
import pstats
import sys
import time
from collections import defaultdict
from functools import wraps
from pathlib import Path

# Ensure project root is on path
project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root))
os.chdir(project_root)

import numpy as np

from sswd_evoepi.config import SimulationConfig, default_config
from sswd_evoepi.spatial import NodeDefinition, build_network


# ═══════════════════════════════════════════════════════════════════════
# TIMER INFRASTRUCTURE
# ═══════════════════════════════════════════════════════════════════════

class ComponentTimers:
    """Global accumulator for per-component wall-clock times."""
    def __init__(self):
        self.totals = defaultdict(float)
        self.counts = defaultdict(int)

    def record(self, name: str, elapsed: float):
        self.totals[name] += elapsed
        self.counts[name] += 1

    def report(self, overall_seconds: float, n_days: int) -> str:
        lines = []
        lines.append("")
        lines.append("=" * 72)
        lines.append(" Component Breakdown (monkey-patched timers)")
        lines.append("=" * 72)
        lines.append(
            f"{'Component':<28} {'Total (s)':>10} {'Per-day (ms)':>13} "
            f"{'% of total':>10} {'Calls':>8}"
        )
        lines.append(
            f"{'─' * 28} {'─' * 10} {'─' * 13} {'─' * 10} {'─' * 8}"
        )

        # Sort by total time descending
        sorted_items = sorted(
            self.totals.items(), key=lambda x: -x[1]
        )
        tracked_total = sum(self.totals.values())

        for name, total_s in sorted_items:
            per_day_ms = (total_s / n_days * 1000) if n_days > 0 else 0
            pct = (total_s / overall_seconds * 100) if overall_seconds > 0 else 0
            calls = self.counts[name]
            lines.append(
                f"{name:<28} {total_s:>10.3f} {per_day_ms:>13.3f} "
                f"{pct:>9.1f}% {calls:>8}"
            )

        lines.append(f"{'─' * 28} {'─' * 10} {'─' * 13} {'─' * 10} {'─' * 8}")
        lines.append(
            f"{'TRACKED TOTAL':<28} {tracked_total:>10.3f} "
            f"{(tracked_total / n_days * 1000) if n_days > 0 else 0:>13.3f} "
            f"{(tracked_total / overall_seconds * 100) if overall_seconds > 0 else 0:>9.1f}%"
        )
        untracked = overall_seconds - tracked_total
        lines.append(
            f"{'UNTRACKED (overhead)':<28} {untracked:>10.3f} "
            f"{(untracked / n_days * 1000) if n_days > 0 else 0:>13.3f} "
            f"{(untracked / overall_seconds * 100) if overall_seconds > 0 else 0:>9.1f}%"
        )
        lines.append(
            f"{'WALL CLOCK TOTAL':<28} {overall_seconds:>10.3f}"
        )
        lines.append("=" * 72)
        return "\n".join(lines)


TIMERS = ComponentTimers()


def timed_wrap(original_fn, component_name):
    """Wrap a function to accumulate its wall-clock time."""
    @wraps(original_fn)
    def wrapper(*args, **kwargs):
        t0 = time.perf_counter()
        result = original_fn(*args, **kwargs)
        TIMERS.record(component_name, time.perf_counter() - t0)
        return result
    return wrapper


# ═══════════════════════════════════════════════════════════════════════
# MONKEY-PATCH KEY FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════

def apply_patches():
    """Wrap functions called in the spatial sim's daily/annual loop."""
    import sswd_evoepi.model as model_mod
    import sswd_evoepi.disease as disease_mod
    import sswd_evoepi.spatial as spatial_mod
    import sswd_evoepi.reproduction as repro_mod
    import sswd_evoepi.spawning as spawning_mod
    import sswd_evoepi.environment as env_mod
    import sswd_evoepi.genetics as genetics_mod

    # Movement — called per-node per-day
    try:
        from sswd_evoepi import movement as mov_mod
        mov_mod.daily_movement = timed_wrap(
            mov_mod.daily_movement, "movement"
        )
    except ImportError:
        pass

    # Disease — called per-node per-day
    disease_mod.daily_disease_update = timed_wrap(
        disease_mod.daily_disease_update, "disease"
    )

    # Pathogen dispersal — called once per day
    spatial_mod.pathogen_dispersal_step = timed_wrap(
        spatial_mod.pathogen_dispersal_step, "pathogen_dispersal"
    )

    # Settlement (continuous) — called per-node per-day
    model_mod.settle_daily_cohorts = timed_wrap(
        model_mod.settle_daily_cohorts, "settlement"
    )

    # Spawning — called per-node per-day during season
    spawning_mod.spawning_step = timed_wrap(
        spawning_mod.spawning_step, "spawning"
    )

    # Daily mortality — called per-node per-day
    model_mod.daily_natural_mortality = timed_wrap(
        model_mod.daily_natural_mortality, "daily_mortality"
    )

    # Daily growth — called per-node per-day
    model_mod.daily_growth_and_aging = timed_wrap(
        model_mod.daily_growth_and_aging, "daily_growth"
    )

    # Environment updates — called inline (SST, flushing, salinity)
    # These are called inline via sst_with_trend / seasonal_flushing,
    # wrap the imported versions in model's namespace
    env_mod.sst_with_trend = timed_wrap(
        env_mod.sst_with_trend, "environment"
    )
    env_mod.seasonal_flushing = timed_wrap(
        env_mod.seasonal_flushing, "environment"
    )

    # Genetics (annual) — compute_allele_frequencies, compute_additive_variance
    genetics_mod.compute_allele_frequencies = timed_wrap(
        genetics_mod.compute_allele_frequencies, "genetics_annual"
    )
    genetics_mod.compute_additive_variance = timed_wrap(
        genetics_mod.compute_additive_variance, "genetics_annual"
    )

    # Larval dispersal (annual) — distribute_larvae
    spatial_mod.distribute_larvae = timed_wrap(
        spatial_mod.distribute_larvae, "larval_dispersal"
    )


# ═══════════════════════════════════════════════════════════════════════
# 3-NODE NETWORK SETUP
# ═══════════════════════════════════════════════════════════════════════

def make_3node_network(seed: int = 42):
    """Create a 3-node test network: Sitka / Howe Sound / Monterey."""
    node_defs = [
        NodeDefinition(
            node_id=0,
            name="Sitka",
            lat=57.05, lon=-135.33,
            subregion="AK-SE",
            habitat_area=500_000.0,
            carrying_capacity=5000,
            is_fjord=False,
            mean_sst=8.5,
            sst_amplitude=3.5,
            sst_trend=0.02,
            salinity=32.0,
            flushing_rate=0.5,
        ),
        NodeDefinition(
            node_id=1,
            name="Howe Sound",
            lat=49.38, lon=-123.25,
            subregion="BC-SW",
            habitat_area=300_000.0,
            carrying_capacity=5000,
            is_fjord=True,
            sill_depth=40.0,
            mean_sst=10.0,
            sst_amplitude=3.0,
            sst_trend=0.02,
            salinity=28.0,
            flushing_rate=0.02,
        ),
        NodeDefinition(
            node_id=2,
            name="Monterey",
            lat=36.62, lon=-121.90,
            subregion="CA-C",
            habitat_area=400_000.0,
            carrying_capacity=5000,
            is_fjord=False,
            mean_sst=13.0,
            sst_amplitude=2.5,
            sst_trend=0.03,
            salinity=33.5,
            flushing_rate=0.5,
        ),
    ]
    return build_network(node_defs, seed=seed)


# ═══════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════

def main():
    output_lines = []
    def log(msg=""):
        print(msg)
        output_lines.append(msg)

    log("=" * 72)
    log(" SSWD-EvoEpi Spatial Simulation Profiling")
    log(f" {time.strftime('%Y-%m-%d %H:%M:%S')}")
    log("=" * 72)
    log()

    # Config
    n_years = 20
    n_days = n_years * 365
    seed = 42
    disease_year = 3

    cfg = default_config()
    # Ensure pathogen evolution is enabled
    if hasattr(cfg, 'pathogen_evolution'):
        cfg.pathogen_evolution.enabled = True

    log(f"Configuration:")
    log(f"  Nodes: 3 (Sitka 8.5°C / Howe Sound 10°C fjord / Monterey 13°C)")
    log(f"  K per node: 5000")
    log(f"  Years: {n_years} ({n_days} days)")
    log(f"  Disease year: {disease_year}")
    log(f"  Pathogen evolution: {getattr(cfg, 'pathogen_evolution', None) and cfg.pathogen_evolution.enabled}")
    log(f"  Seed: {seed}")
    log()

    # Build network
    log("Building 3-node network...")
    network = make_3node_network(seed=seed)
    log(f"  Nodes: {[n.name for n in network.nodes]}")
    log(f"  C matrix shape: {network.C.shape}")
    log(f"  D matrix shape: {network.D.shape}")
    log()

    # Apply monkey patches BEFORE importing run_spatial_simulation
    # (the function uses module-level references, so we patch the modules)
    log("Applying instrumentation patches...")
    apply_patches()
    log()

    # Import after patching — but run_spatial_simulation uses local imports
    # for some functions (sst_with_trend, seasonal_flushing, daily_movement).
    # Those are imported inside the function body, so we need to also
    # patch the modules those are imported FROM.
    from sswd_evoepi.model import run_spatial_simulation

    # ── RUN 1: Manual timers only (warm-up + timing) ─────────────
    log("─" * 72)
    log(" RUN 1: Instrumented run (monkey-patched timers)")
    log("─" * 72)

    # Rebuild network for fresh state
    network = make_3node_network(seed=seed)

    t_start = time.perf_counter()
    result = run_spatial_simulation(
        network=network,
        n_years=n_years,
        disease_year=disease_year,
        initial_infected_per_node=5,
        seed=seed,
        config=cfg,
    )
    t_total = time.perf_counter() - t_start

    log(f"\nTotal wall-clock time: {t_total:.2f}s")
    log(f"Final total pop: {result.final_total_pop}")
    log(f"Peak disease prevalence: {[f'{p:.3f}' for p in result.peak_disease_prevalence]}")

    # Print component breakdown
    report = TIMERS.report(t_total, n_days)
    log(report)

    # ── RUN 2: cProfile ──────────────────────────────────────────
    log()
    log("─" * 72)
    log(" RUN 2: cProfile (top 50 cumulative)")
    log("─" * 72)

    # Rebuild network for fresh state
    network2 = make_3node_network(seed=seed)
    TIMERS.totals.clear()
    TIMERS.counts.clear()

    pr = cProfile.Profile()
    pr.enable()
    result2 = run_spatial_simulation(
        network=network2,
        n_years=n_years,
        disease_year=disease_year,
        initial_infected_per_node=5,
        seed=seed,
        config=cfg,
    )
    pr.disable()

    s = io.StringIO()
    ps = pstats.Stats(pr, stream=s).sort_stats('cumulative')
    ps.print_stats(50)
    cprofile_output = s.getvalue()
    log(cprofile_output)

    # Also print top 50 by tottime
    log()
    log("─" * 72)
    log(" cProfile (top 50 by tottime — self time)")
    log("─" * 72)
    s2 = io.StringIO()
    ps2 = pstats.Stats(pr, stream=s2).sort_stats('tottime')
    ps2.print_stats(50)
    log(s2.getvalue())

    # ── Save output ──────────────────────────────────────────────
    output_path = project_root / "results" / "continuous_mortality" / "profile_output.txt"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.write("\n".join(output_lines))
    log(f"\nOutput saved to: {output_path}")


if __name__ == "__main__":
    main()
