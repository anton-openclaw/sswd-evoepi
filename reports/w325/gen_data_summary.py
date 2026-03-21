#!/usr/bin/env python3
"""Generate data_summary.json from W325 result."""
import json
import numpy as np

REPO = "/home/starbot/.openclaw/workspace/sswd-evoepi"
RESULT = f"{REPO}/results/k1000_scaled_sweep/W325/result_seed42.json"
CONFIG = f"{REPO}/experiments/calibration/W325_config.json"
NPZ = f"{REPO}/results/k1000_scaled_sweep/W325/monthly_seed42.npz"
OUT = f"{REPO}/reports/w325/data_summary.json"

COASTLINE_ORDER = [
    "BJ", "CA-S", "CA-C", "CA-N", "OR", "WA-O", "JDF", "SS-S", "SS-N",
    "BC-C", "BC-N", "AK-FS", "AK-FN", "AK-OC", "AK-PWS", "AK-EG", "AK-WG", "AK-AL"
]

with open(RESULT) as f:
    result = json.load(f)
with open(CONFIG) as f:
    config = json.load(f)

# Load NPZ for population dynamics summary
data = np.load(NPZ, allow_pickle=True)
populations = data["populations"]  # (159, 896)
infected = data["infected"]
sim_days = data["sim_days"]
site_names = data["site_names"]
site_lats = data["site_lats"]

# Map sites to regions — strip trailing numeric site ID
# BJ-001 → BJ, AK-PWS-001 → AK-PWS
def _site_to_region(name):
    parts = str(name).split("-")
    if len(parts) >= 2 and parts[-1].isdigit():
        return "-".join(parts[:-1])
    return "-".join(parts[:2])

regions_per_site = [_site_to_region(n) for n in site_names]

# Compute region-level mean latitude
region_lats = {}
for r, lat in zip(regions_per_site, site_lats):
    region_lats.setdefault(r, []).append(lat)
region_lats = {r: float(np.mean(v)) for r, v in region_lats.items()}

# Population dynamics summary per region
years = 2012 + sim_days / 365.25
region_pop_ts = {}
region_inf_ts = {}
for reg in COASTLINE_ORDER:
    mask = np.array([r == reg for r in regions_per_site])
    region_pop_ts[reg] = populations[:, mask].sum(axis=1).tolist()
    region_inf_ts[reg] = infected[:, mask].sum(axis=1).tolist()

# Build summary
summary = {
    "run_id": "W325",
    "seed": 42,
    "wall_time_seconds": result["wall_time_seconds"],
    "coastline_order": COASTLINE_ORDER,
    "config": config["param_overrides"],
    "scoring": {
        "rmsle": result["scoring"]["rmsle"],
        "within_2x": result["scoring"]["within_2x"],
        "within_5x": result["scoring"]["within_5x"],
        "n_targets": result["scoring"]["n_targets"],
        "per_region": result["scoring"]["per_region"],
    },
    "region_recovery": result["region_recovery"],
    "region_mean_latitude": {r: region_lats.get(r, None) for r in COASTLINE_ORDER},
    "arrival_timing": result["arrival_timing"],
    "overall": result["overall"],
    "evolution": {},
    "pathogen_evolution": {},
    "population_dynamics": {},
}

rd = result["region_details"]
for reg in COASTLINE_ORDER:
    if reg in rd:
        d = rd[reg]
        summary["evolution"][reg] = {
            "yearly_mean_resistance": d["yearly_mean_resistance"],
            "yearly_mean_tolerance": d["yearly_mean_tolerance"],
            "yearly_mean_recovery": d["yearly_mean_recovery"],
            "yearly_va_resistance": d["yearly_va_resistance"],
            "yearly_va_tolerance": d["yearly_va_tolerance"],
            "yearly_va_recovery": d["yearly_va_recovery"],
        }
        summary["pathogen_evolution"][reg] = {
            "final_mean_T_vbnc": d["final_mean_T_vbnc"],
            "final_mean_v_local": d["final_mean_v_local"],
        }
        summary["population_dynamics"][reg] = {
            "n_nodes": d["n_nodes"],
            "peak_pop": d["peak_pop"],
            "final_pop": d["final_pop"],
            "recovery_frac": d["recovery_frac"],
            "crash_pct": d["crash_pct"],
            "yearly_totals": d["yearly_totals"],
            "yearly_recruits": d["yearly_recruits"],
            "yearly_disease_deaths": d["yearly_disease_deaths"],
        }

with open(OUT, "w") as f:
    json.dump(summary, f, indent=2)

print(f"Wrote {OUT}")
print(f"  Regions: {len(summary['evolution'])}")
print(f"  RMSLE: {summary['scoring']['rmsle']:.4f}")
