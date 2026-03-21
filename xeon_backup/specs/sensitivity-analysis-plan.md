# SSWD-EvoEpi Sensitivity Analysis Plan

**Date:** 2026-02-16
**Author:** Anton 🔬
**Status:** DRAFT — awaiting Willem's approval

---

## 1. Objective

Identify which of the model's uncertain parameters most strongly influence key ecological and evolutionary outcomes. This will:
1. Prioritize future empirical work (measure what matters most)
2. Identify parameters that can be safely approximated
3. Reveal parameter interactions driving emergent behavior
4. Build confidence intervals on model predictions

## 2. Method: Sobol Sensitivity Analysis

**Why Sobol over simpler methods (Morris/OAT)?**
- Decomposes variance into first-order (S1) and total-order (ST) effects
- Captures interaction effects (ST - S1 = interactions)
- Global method — explores full parameter space, not just local perturbations
- Gold standard for computationally affordable ABMs

**Sampling:** Saltelli's extension of Sobol sequences (SALib library)
- For k parameters: N × (2k + 2) model evaluations needed
- With k=25 parameters and N=512 base samples: 512 × 52 = **26,624 runs**
- At ~5s/run (500 agents × 20yr): **~37 hours sequential, ~5 hours on 8 cores**

**BUT**: This is aggressive. A phased approach is smarter.

## 3. Phased Approach

### Phase 1: Morris Screening (2-3 hours)
**Goal:** Cheap pre-screen to identify the ~10-15 parameters that actually matter, eliminating the rest.

- Morris method (Elementary Effects): N × (k+1) = ~400 runs for k=25, N=15
- ~30 minutes on 8 cores
- Output: μ* (mean absolute effect) and σ (interaction indicator) for each parameter
- **Decision point:** Parameters with μ* < threshold get fixed at defaults

### Phase 2: Sobol Analysis on Screened Parameters (~10-15 params)
- k_screened parameters × N=1024 base samples
- E.g., k=12: 1024 × 26 = 26,624 runs → ~5 hours on 8 cores
- First-order (S1) and total-order (ST) indices with confidence intervals
- Bootstrap CIs (1000 resamples)

### Phase 3: Deep Dives (if needed)
- 2D interaction maps for top parameter pairs
- Scatter plots for non-linear relationships
- Threshold/bifurcation detection

## 4. The 25 Uncertain Parameters

Mapped from config.py with ranges justified by literature (parameters-extracted.md).

### Disease Module (10 parameters)
| # | Config Path | Default | Low | High | Basis |
|---|------------|---------|-----|------|-------|
| 1 | `disease.a_exposure` | 0.75 | 0.30 | 1.50 | Lupo 2020 analogue (0.61-0.91); ★☆☆ |
| 2 | `disease.K_half` | 87000 | 20000 | 200000 | Lupo 2020 (54K-120K); order-of-magnitude uncertainty |
| 3 | `disease.sigma_1_eff` | 5.0 | 1.0 | 25.0 | Field-effective shedding; no direct measurement |
| 4 | `disease.sigma_2_eff` | 50.0 | 10.0 | 250.0 | Field-effective shedding; I2 >> I1 assumed |
| 5 | `disease.sigma_D` | 15.0 | 3.0 | 75.0 | Saprophytic burst; CE-6 corrected; no empirical anchor |
| 6 | `disease.rho_rec` | 0.05 | 0.0 | 0.20 | Recovery rate: ZERO empirical basis (★☆☆) |
| 7 | `disease.mu_EI1_ref` | 0.57 | 0.20 | 1.00 | E→I1 rate; Prentice 2025 (3-7d range) |
| 8 | `disease.mu_I2D_ref` | 0.173 | 0.08 | 0.35 | I2→D rate; ~7d symptomatic (Prentice 2025) |
| 9 | `disease.P_env_max` | 500.0 | 50.0 | 5000.0 | Background Vibrio input; no measurement |
| 10 | `disease.T_ref` | 20.0 | 17.0 | 23.0 | V. pectenicida T_opt; Lambert 1998 (★★☆) |

### Population Module (7 parameters)
| # | Config Path | Default | Low | High | Basis |
|---|------------|---------|-----|------|-------|
| 11 | `population.F0` | 1e7 | 1e6 | 1e8 | Fecundity: 8-114M range (★☆☆) |
| 12 | `population.gamma_fert` | 4.5 | 1.0 | 10.0 | Fertilization kinetics; theoretical (★☆☆) |
| 13 | `population.settler_survival` | 0.03 | 0.005 | 0.10 | B-H s0; Pisaster proxy <3% (★☆☆) |
| 14 | `population.alpha_srs` | 1.35 | 1.0 | 1.8 | SRS Pareto shape; Árnason 2023 (★★☆) |
| 15 | `population.senescence_age` | 50.0 | 20.0 | 80.0 | Unknown longevity (★☆☆) |
| 16 | `population.k_growth` | 0.08 | 0.03 | 0.15 | VB growth rate; anecdotal (★☆☆) |
| 17 | `population.L_min_repro` | 400.0 | 200.0 | 500.0 | Min reproductive size; estimated (★☆☆) |

### Genetics Module (3 parameters)
| # | Config Path | Default | Low | High | Basis |
|---|------------|---------|-----|------|-------|
| 18 | `genetics.n_additive` | 51 | 10 | 51 | Loci count: 51 (R&R) but may overcount (★★☆) |
| 19 | `genetics.s_het` | 0.19 | 0.05 | 0.35 | EF1A heterozygote advantage; Wares 2016 (★★★ but uncertain) |
| 20 | `genetics.q_ef1a_init` | 0.24 | 0.10 | 0.50 | Initial EF1A freq; geographic variation 0-0.88 (★★★) |

### Spawning Module (3 parameters)
| # | Config Path | Default | Low | High | Basis |
|---|------------|---------|-----|------|-------|
| 21 | `spawning.p_spontaneous_female` | 0.012 | 0.005 | 0.025 | Calibrated but uncertain (★★☆) |
| 22 | `spawning.induction_female_to_male` | 0.80 | 0.40 | 0.95 | Willem's estimate (★★☆) |
| 23 | `spawning.susceptibility_multiplier` | 2.0 | 1.0 | 4.0 | Post-spawning immunosuppression; field obs (★☆☆) |

### Environmental (2 parameters)
| # | Config Path | Default | Low | High | Basis |
|---|------------|---------|-----|------|-------|
| 24 | `disease.T_vbnc` | 12.0 | 8.0 | 15.0 | VBNC midpoint; Erken 2013 (★★☆) |
| 25 | `disease.s_min` | 10.0 | 5.0 | 15.0 | Salinity minimum for Vibrio (★★☆) |

### Special Handling
- **`genetics.n_additive`** (#18): Integer parameter — sample from {10, 20, 30, 40, 51} discrete set
- **`disease.rho_rec`** (#6): Include 0.0 in range (no recovery scenario)
- **`population.F0`** (#11): Log-uniform sampling (spans 2 orders of magnitude)
- **Constraint**: n_loci = n_additive + 1 (automatically enforced)

## 5. Output Metrics (Response Variables)

For each run, extract these from `CoupledSimResult`:

| Metric | Description | Ecological Relevance |
|--------|-------------|---------------------|
| `pop_crash_pct` | (initial - min_pop) / initial × 100 | Severity of epidemic |
| `final_pop_frac` | final_pop / initial_pop | Long-term persistence |
| `recovery` | final_pop > 0.5 × initial_pop (bool→0/1) | Whether population recovers |
| `peak_mortality` | peak_mortality_fraction | Epidemic intensity |
| `resistance_shift` | mean_resistance[final] - mean_resistance[pre-epidemic] | Selection signal |
| `va_retention` | va[final] / va[pre-epidemic] | Evolutionary potential |
| `ef1a_shift` | abs(ef1a[final] - ef1a[0]) | Overdominant locus dynamics |
| `time_to_nadir` | year of minimum population | Speed of crash |
| `total_disease_deaths` | cumulative disease mortality | Cumulative impact |
| `extinction` | final_pop == 0 (bool→0/1) | Complete loss |

## 6. Simulation Configuration

- **N:** 200 agents (faster than 500; captures dynamics; scales O(N^0.62))
- **T:** 20 years (3yr spinup + disease at yr 3 + 17yr recovery window)
- **Temperature:** 15°C constant (mid-range, above VBNC threshold)
- **Salinity:** 30 psu (full marine)
- **Seed:** Unique per run (base_seed + run_index)
- **Spawning:** Enabled (new system — this IS the model we're analyzing)
- **Estimated runtime:** ~2.5s/run at N=200

## 7. Cron Job Architecture

### Job 0: Setup & Validation (`sa-0-setup`)
- Install dependencies (SALib already done)
- Build and test the sensitivity analysis runner script
- Generate Saltelli sample matrices
- Validate 5 edge-case configs (extreme corners of parameter space)
- **Runtime:** ~5 min
- **Schedule:** Immediate (one-shot)

### Job 1: Morris Screening (`sa-1-morris`)
- Run Morris method (15 trajectories × 26 = 390 runs)
- ~390 × 2.5s = 16 min sequential, ~2 min on 8 cores
- Compute μ* and σ for all 25 parameters
- Identify "important" subset (μ* > threshold)
- Save results + screening plot
- **Depends on:** Job 0
- **Schedule:** 15 min after Job 0

### Job 2: Sobol Sampling (`sa-2-sobol-sample`)
- Generate Saltelli samples for screened parameters only
- Validate sample matrix (parameter bounds, constraints)
- Split into batches for parallel execution
- **Depends on:** Job 1 (needs screened parameter list)
- **Schedule:** After Job 1

### Job 3a-3d: Sobol Execution (4 parallel batch jobs)
- `sa-3a-sobol-batch1` through `sa-3d-sobol-batch4`
- Each runs ~25% of the Saltelli samples
- Results saved as numpy arrays
- Checkpointing: save intermediate results every 500 runs
- **Depends on:** Job 2
- **Runtime estimate:** ~1.5 hours each on 8 cores

### Job 4: Sobol Analysis (`sa-4-analysis`)
- Collect all batch results
- Compute S1, ST indices with bootstrap CIs
- Generate:
  - Bar charts (S1 vs ST for each metric)
  - Heatmap (parameter × metric importance matrix)
  - Scatter plots (top 5 parameters vs each metric)
  - Interaction matrix (ST - S1 for each parameter)
- **Depends on:** Jobs 3a-3d all complete
- **Schedule:** After last batch completes

### Job 5: Report & Email (`sa-5-report`)
- Compile markdown report with embedded figures
- Email to Willem with key findings
- Commit all results to git
- **Depends on:** Job 4
- **Schedule:** Morning (target ~6-8 AM PST)

### Fault Tolerance
- Each batch saves progress to disk every 500 runs
- If a batch crashes, it can resume from last checkpoint
- Individual run failures logged but don't stop the batch (NaN for that row)
- All data in `results/sensitivity/` directory

## 8. Risk Assessment

| Risk | Mitigation |
|------|-----------|
| Parameter combos crash the model | Try/except per run; log failures; NaN fill |
| Extreme params → numerical instability | Validate bounds; cap populations at 10× K |
| Too many runs for overnight | Morris screening reduces Sobol parameter count |
| Memory issues with parallel runs | Cap at 8 parallel (measured sweet spot) |
| N=200 too small for reliable dynamics | Run 10 seed replicates for Morris validation |
| `n_additive` changes break model | Recompute n_loci = n_additive + 1; regenerate effect sizes |
| Some params interact in unexpected ways | Sobol ST captures this; Phase 3 deep dives if needed |

## 9. Timeline

| Time (PST) | Job | Duration |
|------------|-----|----------|
| ~12:00 AM | sa-0-setup | 5 min |
| ~12:10 AM | sa-1-morris | 15-30 min |
| ~12:45 AM | sa-2-sobol-sample | 5 min |
| ~1:00 AM | sa-3a through sa-3d (parallel) | ~1.5-2 hrs each |
| ~3:30 AM | sa-4-analysis | 15 min |
| ~4:00 AM | sa-5-report + email | 10 min |
| **~4:30 AM** | **All complete** | |

Conservative estimate: done by 6 AM even with retries.

## 10. Deliverables

1. **Morris screening plot** — μ* vs σ for all 25 parameters
2. **Sobol indices table** — S1, ST, CI for each parameter × metric
3. **Importance heatmap** — which parameters matter for which outcomes
4. **Scatter plots** — non-linear relationships for top parameters
5. **Interaction matrix** — which parameters interact
6. **Written report** — key findings, recommendations for empirical priorities
7. **All raw data** — JSON/NPZ for reproducibility
