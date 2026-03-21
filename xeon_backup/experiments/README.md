# Reintroduction Experiment: Monterey Outplanting

## Motivation

*Pycnopodia helianthoides* (sunflower sea star) populations crashed >90% across their range due to Sea Star Wasting Disease (SSWD) starting in 2013. In December 2025, the first successful captive-bred outplanting placed 47/48 juveniles at a Monterey site with 100% 4-week survival.

This experiment uses our coupled eco-evolutionary epidemiological ABM to evaluate: **What combination of restoration scale and genetic composition maximizes the probability of population recovery?**

## Network

**Full 907-node network** covering the entire Pycnopodia range from the Aleutian Islands (61°N) to Baja California (27°N). This is NOT the 11-node stepping-stone used for sensitivity analysis — that was designed for SA efficiency with uniform spacing. The full network captures:

- Realistic population connectivity (overwater distance-based dispersal)
- Regional population structure (18 regions: AK-WG, AK-AL, AK-EG, AK-PWS, AK-FN, AK-FS, AK-OC, BC-N, BC-C, JDF, SS-N, SS-S, WA-O, OR, CA-N, CA-C, CA-S, BJ)
- The actual Monterey outplanting site (CA-C-043: 36.62°N, -121.90°W)

Site data: `data/nodes/all_sites.json` (907 sites, systematic Region-NNN naming)
Distance matrix: `results/overwater/distance_matrix.npz` (907×907 precomputed overwater distances)
Carrying capacity: K = 5,000 per node (SSWD crashes populations to <5% of K within years, so effective agent counts are small)

## SST Data

**Primary**: Per-site monthly SST from NOAA OISST v2.1 (2002–2025). Each of the 907 sites maps to its nearest OISST 0.25° grid cell (532 unique cells). Data stored at `data/sst/site_sst/`.

**Fallback**: If per-site SST is not yet available, 11 reference nodes provide nearest-latitude interpolation (same as legacy method).

**Projections**: Bias-corrected CMIP6 SSP2-4.5 (GFDL-ESM4 + IPSL-CM6A-LR ensemble mean), 2026–2100.

## Factorial Design

### Restoration Levels
| Level | Individuals Released |
|---|---|
| Partial | 50 |
| Medium | 500 |
| Full | 5,000 |

### Genetic Backgrounds
| Background | Resistance | Tolerance | Recovery | Description |
|---|---|---|---|---|
| `pre_sswd` | 0.15 | 0.10 | 0.02 | Pre-2013 wild-type standing variation |
| `survivors_2019` | 0.18 | 0.14 | 0.03 | Post-epidemic survivors (~6 years selection) |
| `bred_1gen` | ~0.32 | ~0.19 | ~0.03 | 1 generation selective breeding |
| `bred_2gen` | ~0.44 | ~0.23 | ~0.03 | 2 generations selective breeding |
| `bred_5gen` | ~0.77 | ~0.41 | ~0.05 | 5 generations selective breeding |
| `optimal` | 1.00 | 0.10 | 0.02 | Theoretical max (all R-loci fixed) |

Selective breeding uses truncation selection (top 20%) with weighted index: 0.7 × resistance + 0.2 × tolerance + 0.1 × recovery.

### Total: 19 scenarios
- 3 restoration levels × 6 genetic backgrounds = 18 treatment scenarios
- 1 baseline (no intervention)
- Each scenario replicated (3 reps for demo, 10 for full)

## Release Protocol

- **Site**: CA-C-043 (Monterey outplanting site, 36.62°N, -121.90°W)
- **Timing**: Simulation day 8760 (January 1, 2026)
- **Age**: 1–2 year old juveniles (365–730 days)
- **Genetics mode**: `allele_freqs` (per-locus protective allele frequencies)
- **Tagging**: All released individuals marked with `Origin.CAPTIVE_BRED`

## Timeline

| Year | Event |
|---|---|
| 2002 | Simulation start (11 years pre-epidemic spinup) |
| 2013 | SSWD epidemic onset |
| 2019 | Approximate genetic state of survivors |
| 2026 | Captive-bred release at Monterey |
| 2050 | Simulation end (24 years post-release) |

## Tracked Outputs

Per scenario:
1. **Population trajectory** — yearly total population (all 907 nodes)
2. **Monterey population** — yearly pop at release site
3. **Released individual survival** — how many captive-bred survive per year
4. **Trait evolution at Monterey** — mean resistance, tolerance, recovery
5. **Region-aggregated population** — yearly pop summed by 18 regions
6. **Larval dispersal** — offspring spreading from Monterey to adjacent nodes

## How to Run

```bash
# Demo: 3 replicates, 8 cores
python3 experiments/reintroduction_monterey.py --mode demo

# Full: 10 replicates on Xeon
python3 experiments/reintroduction_monterey.py --mode full --cores 48

# Single scenario for testing
python3 experiments/reintroduction_monterey.py --mode single \
    --genetics bred_2gen --restoration medium

# Custom K (default 5000)
python3 experiments/reintroduction_monterey.py --mode demo --K 1000
```

## Expected Outputs

Results saved to `results/reintroduction/`:
- `results_demo.json` — demo run results
- `results_full.json` — full run results

Each result contains per-scenario: population trajectories, trait evolution, release survival, regional aggregates.

## Key Questions

1. **Minimum viable release size**: How many individuals are needed for detectable population effect?
2. **Genetic gain vs. effort**: Does 5 generations of breeding meaningfully outperform using 2019 survivors?
3. **Spillover**: Do released individuals' offspring disperse to repopulate adjacent nodes?
4. **Resistance threshold**: Is there a resistance level above which Monterey population persists?
5. **Optimal strategy**: What combination of n_released × genetic background maximizes range-wide recovery per dollar?
