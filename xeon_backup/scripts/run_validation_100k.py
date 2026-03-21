#!/usr/bin/env python3
"""100K validation run — 5 nodes, K=100K each, 20 years."""
import sys, time, json, os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sswd_evoepi.config import default_config
from sswd_evoepi.model import run_spatial_simulation
from sswd_evoepi.spatial import NodeDefinition, SpatialNode, MetapopulationNetwork

K = 100_000
SEED = 42
N_YEARS = 20
DISEASE_YEAR = 3
OUTDIR = "results/validation_100k"

os.makedirs(OUTDIR, exist_ok=True)

cfg = default_config()
cfg.movement.enabled = True
cfg.movement.substeps_per_day = 1

node_specs = [
    ("Sitka",      57.05, -135.33, 8.5,  3.0, 30.0, False, 0.50),
    ("Howe_Sound", 49.38, -123.24, 9.5,  4.0, 28.0, True,  0.15),
    ("SJI",        48.53, -123.01, 9.0,  3.5, 30.0, False, 0.50),
    ("Newport",    44.63, -124.05, 10.0, 3.5, 33.0, False, 0.50),
    ("Monterey",   36.62, -121.90, 12.5, 2.5, 33.5, False, 0.50),
]

node_defs = []
for i, (name, lat, lon, sst, amp, sal, fjord, flush) in enumerate(node_specs):
    kwargs = dict(
        node_id=i, name=name, lat=lat, lon=lon,
        subregion=name[:2].upper(), habitat_area=666667.0,
        carrying_capacity=K, mean_sst=sst, sst_amplitude=amp,
        salinity=sal, sst_trend=0.02, flushing_rate=flush, is_fjord=fjord,
    )
    if fjord:
        kwargs["sill_depth"] = 30.0
    node_defs.append(NodeDefinition(**kwargs))

# Build network
n = len(node_defs)
D = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        if i != j:
            dlat = np.radians(node_defs[i].lat - node_defs[j].lat)
            dlon = np.radians(node_defs[i].lon - node_defs[j].lon)
            a = (np.sin(dlat/2)**2 +
                 np.cos(np.radians(node_defs[i].lat)) *
                 np.cos(np.radians(node_defs[j].lat)) *
                 np.sin(dlon/2)**2)
            D[i, j] = 6371 * 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

nodes = [SpatialNode(definition=nd) for nd in node_defs]
rng = np.random.default_rng(SEED)
network = MetapopulationNetwork(nodes=nodes, C=np.zeros((n,n)), D=np.zeros((n,n)), distances=D, rng=rng)

print(f"=== 100K Validation (vectorized disease) ===")
print(f"K={K:,}, {n} nodes, {N_YEARS}yr, seed={SEED}")
print(f"Total agents: ~{K*n:,}")
print(f"Config: rho_rec={cfg.disease.rho_rec}, target_mean_c={cfg.genetics.target_mean_c}")
sys.stdout.flush()

t0 = time.perf_counter()
r = run_spatial_simulation(
    network=network,
    n_years=N_YEARS,
    disease_year=DISEASE_YEAR,
    initial_infected_per_node=10,
    seed=SEED,
    config=cfg,
)
elapsed = time.perf_counter() - t0

# Extract results
names = [nd.name for nd in node_defs]
results = {"elapsed_s": round(elapsed, 1), "seed": SEED, "K": K, "n_nodes": n}
for k_idx, name in enumerate(names):
    final_pop = int(r.yearly_pop[k_idx, -1])
    crash = 1.0 - final_pop / K
    total_deaths = int(np.sum(r.yearly_disease_deaths[k_idx, :]))
    total_recs = int(np.sum(r.yearly_recoveries[k_idx, :]))
    mean_r_final = float(r.yearly_mean_resistance[k_idx, -1])
    mean_t_final = float(r.yearly_mean_tolerance[k_idx, -1])
    mean_c_final = float(r.yearly_mean_recovery[k_idx, -1])
    results[name] = {
        "final_pop": final_pop, "crash_pct": round(crash*100, 1),
        "disease_deaths": total_deaths, "recoveries": total_recs,
        "mean_r": round(mean_r_final, 4), "mean_t": round(mean_t_final, 4),
        "mean_c": round(mean_c_final, 4),
    }
    print(f"  {name}: crash={crash:.1%}, deaths={total_deaths:,}, recs={total_recs}, "
          f"Δr={mean_r_final - 0.15:+.4f}, Δt={mean_t_final - 0.10:+.4f}, Δc={mean_c_final - 0.02:+.4f}")
    sys.stdout.flush()

total_alive = sum(results[n]["final_pop"] for n in names)
total_crash = 1 - total_alive / (K * len(names))
total_recs = sum(results[n]["recoveries"] for n in names)
print(f"\nOverall: {elapsed:.1f}s, crash={total_crash:.1%}, total_recoveries={total_recs}")

with open(f"{OUTDIR}/results.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"Saved to {OUTDIR}/results.json")
