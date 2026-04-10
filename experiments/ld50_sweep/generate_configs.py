#!/usr/bin/env python3
"""Generate 256 LD50 grid search configs (4^4 parameter combinations).

Grid design (from LD50_GRID_SEARCH_PLAN.md):
  r_growth:    [0.08, 0.16, 0.32, 0.64]    — pathogen growth rate
  delta_clear: [0.25, 0.40, 0.55, 0.80]    — base immune clearance rate
  LD50_base:   [1e4, 5e4, 2.5e5, 1e6]      — load for 50% daily mortality
  sigma_load:  [50, 150, 450, 1350]         — shedding rate

Each config starts from W330 calibrated baseline and enables the
load-dependent disease model with swept parameters. Fixed LD50 params
(n_hill=3, p_death_max=0.15, etc.) use defaults from LoadDependentSection.

Naming: LD50_G{g}_C{c}_L{l}_S{s}.json  (indices 0-3)

The calibration_runner.py applies overrides via dotted-key traversal:
  "disease_load_dependent.r_growth" → config.disease_load_dependent.r_growth
This works out of the box — no runner patch needed.
"""

import json
import os
from itertools import product
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent.parent
BASELINE_FILE = REPO / "experiments" / "calibration" / "W330_config.json"
OUTPUT_DIR = Path(__file__).resolve().parent / "configs"

# ── Parameter grid ─────────────────────────────────────────────────────
PARAM_GRID = {
    "r_growth":    [0.08, 0.16, 0.32, 0.64],
    "delta_clear": [0.25, 0.40, 0.55, 0.80],
    "LD50_base":   [1e4, 5e4, 2.5e5, 1e6],
    "sigma_load":  [50.0, 150.0, 450.0, 1350.0],
}

# Map short names to index letters used in filename
INDEX_LETTERS = ["G", "C", "L", "S"]
PARAM_NAMES = list(PARAM_GRID.keys())
PARAM_LEVELS = [PARAM_GRID[name] for name in PARAM_NAMES]

# Defaults (for reference — the non-swept LD50 params use LoadDependentSection defaults)
FIXED_LD50_DEFAULTS = {
    "n_hill": 3.0,
    "p_death_max": 0.15,
    "L_max": 1e6,
    "Ea_growth": 5000.0,
    "L_clear": 10.0,
    "L_symp": 1e4,
    "L_ref": 1e5,
    "alpha_reinfect": 0.01,
    "L_init_base": 1e5,
}


def main():
    # Load W330 baseline param_overrides
    if not BASELINE_FILE.exists():
        print(f"ERROR: Baseline config not found: {BASELINE_FILE}")
        print("  Expected W330_config.json from calibration round W330")
        raise SystemExit(1)

    with open(BASELINE_FILE) as f:
        baseline = json.load(f)

    if "param_overrides" in baseline:
        base_overrides = baseline["param_overrides"]
    else:
        base_overrides = baseline

    print(f"Loaded W330 baseline: {len(base_overrides)} param overrides")
    print(f"Grid: {' x '.join(str(len(v)) for v in PARAM_LEVELS)} = "
          f"{len(list(product(*PARAM_LEVELS)))} combinations")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    count = 0

    for indices in product(*[range(len(levels)) for levels in PARAM_LEVELS]):
        g, c, l, s = indices
        values = {
            name: PARAM_LEVELS[i][idx]
            for i, (name, idx) in enumerate(zip(PARAM_NAMES, indices))
        }

        scenario_name = f"LD50_G{g}_C{c}_L{l}_S{s}"

        # Build param_overrides: W330 baseline + load_dependent overrides
        param_overrides = dict(base_overrides)  # copy

        # Enable load-dependent disease model
        param_overrides["disease_load_dependent.enabled"] = True

        # Set swept parameters
        for param_name, param_value in values.items():
            param_overrides[f"disease_load_dependent.{param_name}"] = param_value

        config = {
            "param_overrides": param_overrides,
            "scenario": {
                "name": scenario_name,
                "type": "ld50_sweep",
                "param_indices": {
                    INDEX_LETTERS[i]: idx for i, idx in enumerate(indices)
                },
                "param_values": {
                    name: value for name, value in values.items()
                },
                "param_names": PARAM_NAMES,
                "grid_levels": {
                    name: levels for name, levels in zip(PARAM_NAMES, PARAM_LEVELS)
                },
            },
        }

        out_path = OUTPUT_DIR / f"{scenario_name}.json"
        with open(out_path, "w") as f:
            json.dump(config, f, indent=2)

        count += 1

    print(f"\nGenerated {count} configs in {OUTPUT_DIR}/")
    print(f"  Naming: LD50_G{{0-3}}_C{{0-3}}_L{{0-3}}_S{{0-3}}.json")
    print(f"  Each config enables disease_load_dependent with swept params")
    print(f"  + W330 calibrated baseline ({len(base_overrides)} disease/spatial params)")

    # Print a sample config summary
    sample = "LD50_G0_C0_L0_S0"
    with open(OUTPUT_DIR / f"{sample}.json") as f:
        sample_cfg = json.load(f)
    print(f"\nSample ({sample}):")
    for name in PARAM_NAMES:
        key = f"disease_load_dependent.{name}"
        print(f"  {key} = {sample_cfg['param_overrides'][key]}")
    print(f"  disease_load_dependent.enabled = {sample_cfg['param_overrides']['disease_load_dependent.enabled']}")
    print(f"  Total param_overrides keys: {len(sample_cfg['param_overrides'])}")


if __name__ == "__main__":
    main()
