#!/usr/bin/env python3
"""Generate 40 reintroduction experiment configs.

Design:
  L1 (CA): Release at all CA nodes (CA-S, CA-C, CA-N)
  L2 (WA): Release at all WA-area nodes (JDF, WA-O, SS-S)
  R1-R5: resistance = 0.31, 0.15, 0.20, 0.50, 1.0
  D1-D4: n_individuals per site = 50, 200, 500, 1000
  Release day: year 13 = day 4745
  genetics_mode: trait_targets

Total: 2 locations × 5 resistance × 4 density = 40 configs
"""

import json
import os
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent.parent
SITES_FILE = REPO / "data" / "nodes" / "all_sites.json"
BASELINE_FILE = REPO / "experiments" / "forecast" / "F01_baseline.json"
OUTPUT_DIR = Path(__file__).resolve().parent / "configs"

# Location definitions
LOCATIONS = {
    "L1": {"label": "CA", "regions": {"CA-S", "CA-C", "CA-N"}},
    "L2": {"label": "WA", "regions": {"JDF", "WA-O", "SS-S"}},
}

# Resistance levels
RESISTANCES = {
    "R1": 0.31,
    "R2": 0.15,
    "R3": 0.20,
    "R4": 0.50,
    "R5": 1.0,
}

# Density (n_individuals per site)
DENSITIES = {
    "D1": 50,
    "D2": 200,
    "D3": 500,
    "D4": 1000,
}

RELEASE_DAY = 4745  # year 13 = 13 * 365


def main():
    # Load sites
    with open(SITES_FILE) as f:
        sites = json.load(f)

    # Build region -> node_id mapping
    location_nodes = {}
    for loc_key, loc_info in LOCATIONS.items():
        node_ids = []
        for i, site in enumerate(sites):
            if site.get("region") in loc_info["regions"]:
                node_ids.append(i)
        location_nodes[loc_key] = sorted(node_ids)
        print(f"{loc_key} ({loc_info['label']}): {len(node_ids)} nodes")

    # Load baseline param_overrides
    with open(BASELINE_FILE) as f:
        baseline = json.load(f)
    param_overrides = baseline["param_overrides"]

    # Generate configs
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    count = 0

    for loc_key in sorted(LOCATIONS.keys()):
        node_ids = location_nodes[loc_key]
        for res_key in sorted(RESISTANCES.keys()):
            resistance = RESISTANCES[res_key]
            for den_key in sorted(DENSITIES.keys()):
                n_individuals = DENSITIES[den_key]

                scenario_name = f"REINTRO_{loc_key}_{res_key}_{den_key}"

                # Build release events - one per node
                release_events = []
                for node_id in node_ids:
                    release_events.append({
                        "time_step": RELEASE_DAY,
                        "node_id": node_id,
                        "n_individuals": n_individuals,
                        "genetics_mode": "trait_targets",
                        "trait_targets": {"resistance": resistance},
                    })

                config = {
                    "param_overrides": param_overrides,
                    "release_events": release_events,
                }

                out_path = OUTPUT_DIR / f"{scenario_name}.json"
                with open(out_path, "w") as f:
                    json.dump(config, f, indent=2)

                count += 1
                print(f"  {scenario_name}: {len(release_events)} release events, "
                      f"resistance={resistance}, n={n_individuals}")

    print(f"\nGenerated {count} configs in {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
