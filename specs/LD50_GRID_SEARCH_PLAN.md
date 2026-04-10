# LD50 Grid Search Plan

## Objective
Explore the load-dependent disease parameter space with precomputed simulations,
then provide an interactive interface to explore outcomes.

## Compute Budget
- Machine: Ryzen 7 5700 (8c/16t), 32GB RAM
- Budget: 14 days
- K=1000: 2.2 hr/run, 714MB, 8 parallel → ~1,200 runs max
- Estimated: 768 runs → 8.8 days (leaves buffer for UI)

## Grid Design: 4 parameters × 4 levels × 3 seeds = 768 runs

### Parameters Swept (4)

1. **r_growth** — pathogen growth rate
   Levels: [0.08, 0.16, 0.32, 0.64]
   Rationale: from slow (immune system usually wins) to fast (overwhelms most hosts)
   Default: 0.16

2. **delta_clear** — base immune clearance rate  
   Levels: [0.25, 0.40, 0.55, 0.80]
   Rationale: from weak immunity to strong; must compare against r_growth
   Default: 0.55

3. **LD50_base** — pathogen load for 50% daily mortality
   Levels: [1e4, 5e4, 2.5e5, 1e6]
   Rationale: from lethal at low load to very tolerant of high load
   Default: 5e4

4. **sigma_load** — shedding rate (transmission feedback)
   Levels: [50, 150, 450, 1350]
   Rationale: log-spaced; controls epidemic speed and R0
   Default: 150

### Fixed Parameters
- All other LD50 params at defaults (n_hill=3, p_death_max=0.15, etc.)
- W330 calibrated baseline (disease params, spatial, SST)
- K=1000, K_ref=5000
- Years=38 (2012-2050), disease_year=1
- 3 seeds per combination: 42, 137, 256

### Why These 4 Parameters
- r_growth / delta_clear **ratio** controls whether infection is clearable.
  Sweeping both independently reveals the transition boundary.
- LD50_base controls **lethality** — same infection dynamics but different death toll
- sigma_load controls **transmission** — same within-host dynamics but different spread

Together they span the three axes of disease ecology:
  Growth/Clearance → individual outcome
  LD50 → population mortality  
  Shedding → epidemic dynamics

### Output per Run
- result_seed{N}.json (scoring, region recovery, yearly metrics)
- monthly_seed{N}.npz (monthly population + infection + load per node)

### Naming Convention
```
results/ld50_sweep/G{g}_C{c}_L{l}_S{s}/result_seed{N}.json
```
Where g,c,l,s are parameter indices 0-3.

## Interactive Interface

After simulations complete, build a local web dashboard (Streamlit or Panel):

- **Parameter sliders** for all 4 grid dimensions
- **Region recovery heatmap** showing all 18 regions
- **Population trajectories** (time series per region)
- **Infection curves** (prevalence over time)
- **Mean pathogen load** over time
- **Resistance evolution** over time
- Seed toggle (show individual seeds or mean±SD)
- Nearest-neighbor interpolation for slider positions between grid points

All data precomputed — the UI just loads and displays.

## Schedule
- Day 0: Generate configs + launch script, start running
- Days 1-9: Simulations running (8 parallel, ~96 per day)
- Day 9: Extract results into dashboard data format
- Days 10-12: Build interactive dashboard
- Days 13-14: Buffer / refinement
