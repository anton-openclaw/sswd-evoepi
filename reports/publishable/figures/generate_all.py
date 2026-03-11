#!/usr/bin/env python3
"""Master runner — generate all publishable figures and copy progress report figs."""
import sys, os, shutil, pathlib, importlib, time

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

# Figures to generate (module name, description)
FIGURE_SCRIPTS = [
    ('gen_fig_coastwide_population', 'Fig 1: Coast-wide population'),
    ('gen_fig_regional_recovery_bars', 'Fig 2: Regional recovery bars'),
    ('gen_fig_regional_trajectories', 'Fig 3: Regional trajectories'),
    ('gen_fig_climate_effect', 'Fig 4: Climate effect'),
    ('gen_fig_resistance_evolution', 'Fig 5: Resistance evolution'),
    ('gen_fig_genetic_diversity', 'Fig 6: Genetic diversity'),
    ('gen_fig_pathogen_traits', 'Fig 7: Pathogen traits'),
    ('gen_fig_evo_comparison', 'Fig 8: Evo comparison'),
    ('gen_fig_recovery_vs_latitude', 'Fig 9: Recovery vs latitude'),
    ('gen_fig_heatmap_population', 'Fig 10: Population heatmap'),
]

# Progress report figures to copy
PROGRESS_DIR = HERE.parent.parent / 'progress_report' / 'figures'
COPY_FIGS = [
    'fig_network_map.png',
    'fig_within_site_disease.png',
    'fig_disease_snapshots.png',
]


def main():
    t0 = time.time()
    print("=" * 60)
    print("Generating all publishable figures")
    print("=" * 60)

    ok, fail = 0, 0
    for mod_name, desc in FIGURE_SCRIPTS:
        print(f"\n▸ {desc} ({mod_name})")
        try:
            mod = importlib.import_module(mod_name)
            mod.main()
            ok += 1
        except Exception as e:
            print(f"  ✗ FAILED: {e}")
            import traceback; traceback.print_exc()
            fail += 1

    # Copy progress report figures
    print(f"\n{'=' * 60}")
    print("Copying progress report figures")
    print("=" * 60)
    for fname in COPY_FIGS:
        src = PROGRESS_DIR / fname
        dst = HERE / fname
        if src.exists():
            shutil.copy2(src, dst)
            print(f"  ✓ Copied {fname}")
        else:
            print(f"  ⚠ Not found: {src}")

    elapsed = time.time() - t0
    print(f"\n{'=' * 60}")
    print(f"Done: {ok} generated, {fail} failed, {elapsed:.1f}s elapsed")
    print("Output directory:", HERE)
    print("=" * 60)

    # List all figures
    pngs = sorted(HERE.glob('fig_*.png'))
    print(f"\nFigures ({len(pngs)}):")
    for p in pngs:
        sz = p.stat().st_size / 1024
        print(f"  {p.name:45s} {sz:7.0f} KB")

    return fail == 0


if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1)
