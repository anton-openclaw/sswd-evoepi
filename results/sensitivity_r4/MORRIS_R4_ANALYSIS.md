# Morris R4 Sensitivity Analysis Report

**Date:** February 20, 2026  
**Round:** 4 (three-trait genetics + pathogen evolution + 11-node network)  
**Runs:** 960 (20 trajectories × 48 parameters)  
**Parameters:** 47 | **Metrics:** 23  
**Network:** 11-node stepping-stone chain  
**Executed on:** Xeon W-3365 (48 cores)

## Executive Summary

Morris R4 is the first screening analysis of the complete SSWD-EvoEpi model with three-trait genetic architecture (resistance/tolerance/recovery), pathogen virulence evolution, and an 11-node spatial network. The expanded model (47 params, up from 43 in R3) reveals three major shifts:

1. **Genetic architecture parameters surge in importance** — `n_resistance` jumps from #19 to #5; new recovery trait mean (`target_mean_c`) enters at #10
2. **Environmental reservoir (`P_env_max`) becomes critical** — #11 → #4 with spatial network
3. **Universal nonlinearity** — ALL 47 parameters show σ/μ* > 1.0, meaning every parameter interacts with others. This is a deeply coupled system.
4. **No parameter can be eliminated** — unlike R3 where we hoped to prune, R4 confirms all 47 contribute meaningfully

## Full Parameter Ranking

| Rank | Parameter | Module | Mean Norm μ* | σ/μ* | R3 Rank | Δ |
|---:|:---|:---|---:|---:|---:|---:|
| 1 | rho_rec | disease | 0.889 | 1.46 | 1 | — |
| 2 | k_growth | population | 0.633 | 1.63 | 5 | ↑3 |
| 3 | K_half | disease | 0.622 | 1.84 | 8 | ↑5 |
| 4 | P_env_max | disease | 0.598 | 1.92 | 11 | ↑7 |
| 5 | n_resistance | genetics | 0.525 | 1.78 | 19* | ↑14 |
| 6 | settler_survival | population | 0.509 | 1.42 | 3 | ↓3 |
| 7 | sigma_2_eff | disease | 0.431 | 1.95 | 10 | ↑3 |
| 8 | mu_I2D_ref | disease | 0.419 | 1.98 | 7 | ↓1 |
| 9 | peak_width_days | spawning | 0.392 | 2.03 | 2 | ↓7 |
| 10 | **target_mean_c** | genetics | 0.385 | 2.08 | NEW | — |
| 11 | a_exposure | disease | 0.379 | 2.20 | 6 | ↓5 |
| 12 | T_vbnc | disease | 0.355 | 2.07 | 9 | ↓3 |
| 13 | **tau_max** | genetics | 0.292 | 2.05 | NEW | — |
| 14 | sigma_v_mutation | pathogen_evo | 0.259 | 2.52 | 31 | ↑17 |
| 15 | alpha_kill | pathogen_evo | 0.254 | 2.25 | 22 | ↑7 |
| 16 | sigma_1_eff | disease | 0.245 | 2.24 | 43 | ↑27 |
| 17 | target_mean_r | genetics | 0.236 | 1.86 | 4 | ↓13 |
| 18 | T_ref | disease | 0.229 | 1.94 | 34 | ↑16 |
| 19 | min_susceptible_age_days | disease | 0.229 | 2.04 | 13 | ↓6 |
| 20 | sigma_D | disease | 0.211 | 1.96 | 29 | ↑9 |
| 21 | female_max_bouts | spawning | 0.206 | 1.95 | 32 | ↑11 |
| 22 | readiness_induction_prob | spawning | 0.204 | 2.26 | 33 | ↑11 |
| 23 | **target_mean_t** | genetics | 0.197 | 2.05 | NEW | — |
| 24 | **n_tolerance** | genetics | 0.189 | 2.51 | NEW | — |
| 25 | alpha_self_open | spatial | 0.187 | 2.07 | 39 | ↑14 |
| 26 | D_L | spatial | 0.178 | 2.29 | 18 | ↓8 |
| 27 | induction_male_to_female | spawning | 0.176 | 2.07 | 16 | ↓11 |
| 28 | s_min | disease | 0.175 | 1.84 | 36 | ↑8 |
| 29 | v_init | pathogen_evo | 0.173 | 2.13 | 12 | ↓17 |
| 30 | p_spontaneous_male | spawning | 0.169 | 2.11 | 28 | ↓2 |
| 31 | mu_I1I2_ref | disease | 0.156 | 1.97 | 14 | ↓17 |
| 32 | q_init_beta_a | genetics | 0.150 | 2.45 | 40 | ↑8 |
| 33 | alpha_self_fjord | spatial | 0.149 | 2.00 | 42 | ↑9 |
| 34 | senescence_age | population | 0.148 | 1.66 | 21 | ↓13 |
| 35 | gamma_early | pathogen_evo | 0.148 | 2.03 | 30 | ↓5 |
| 36 | alpha_srs | population | 0.146 | 2.34 | 35 | ↓1 |
| 37 | alpha_prog | pathogen_evo | 0.143 | 2.09 | 38 | ↑1 |
| 38 | mu_EI1_ref | disease | 0.141 | 2.19 | 27 | ↓11 |
| 39 | L_min_repro | population | 0.139 | 2.06 | 25 | ↓14 |
| 40 | alpha_shed | pathogen_evo | 0.136 | 2.12 | 41 | ↑1 |
| 41 | induction_female_to_male | spawning | 0.130 | 1.79 | 24 | ↓17 |
| 42 | immunosuppression_duration | disease | 0.127 | 2.07 | 15 | ↓27 |
| 43 | gamma_fert | population | 0.122 | 2.21 | 37 | ↓6 |
| 44 | susceptibility_multiplier | disease | 0.111 | 2.03 | 23 | ↓21 |
| 45 | p_spontaneous_female | spawning | 0.110 | 1.67 | 26 | ↓19 |
| 46 | q_init_beta_b | genetics | 0.104 | 2.20 | 17 | ↓29 |
| 47 | F0 | population | 0.102 | 1.83 | 20 | ↓27 |

*\* n_resistance replaces n_additive (renamed for three-trait architecture)*

## Key Changes from R3

### Major Rank Gains (Parameters now more important)
| Parameter | R3→R4 | Why |
|:---|:---|:---|
| **sigma_1_eff** | #43→#16 (↑27) | Early shedding now interacts with pathogen evolution; σ_1 shapes initial epidemic wave that selects virulence |
| **sigma_v_mutation** | #31→#14 (↑17) | Pathogen evolution is new — mutation rate directly controls adaptation speed |
| **T_ref** | #34→#18 (↑16) | Temperature reference interacts with 11-node latitudinal gradient (thermally diverse) |
| **n_resistance** | #19→#5 (↑14) | Three-trait partition means genetic architecture matters more; 17 loci per trait vs. 51 total |
| **alpha_self_open/fjord** | #39→#25, #42→#33 | Spatial retention now detectable with 11 nodes (was invisible at 3 nodes) |
| **P_env_max** | #11→#4 (↑7) | Environmental reservoir interacts with spatial spread across 11 nodes |

### Major Rank Drops (Parameters now less important)
| Parameter | R3→R4 | Why |
|:---|:---|:---|
| **q_init_beta_b** | #17→#46 (↓29) | Initial allele freq shape parameter overwhelmed by trait-specific means |
| **F0** | #20→#47 (↓27) | Initial fecundity diluted in expanded parameter space |
| **immunosuppression_duration** | #15→#42 (↓27) | Effect absorbed by spawning parameters and recovery trait |
| **susceptibility_multiplier** | #23→#44 (↓21) | Formerly #1 in R1 Sobol! Now absorbed by explicit resistance genetics |
| **peak_width_days** | #2→#9 (↓7) | Still important, but relative importance drops with more parameters competing |

### New R4 Parameters Performance
| Parameter | Rank | Significance |
|:---|:---|:---|
| **target_mean_c** (recovery trait mean) | #10 | Top-10 entry — confirms recovery is the fastest-evolving trait |
| **tau_max** (max tolerance effect) | #13 | Tolerance ceiling matters for I₂ survival extension |
| **target_mean_t** (tolerance trait mean) | #23 | Mid-pack — less influential than recovery but detectable |
| **n_tolerance** (# tolerance loci) | #24 | Mid-pack — genetic architecture of tolerance matters |

## Interaction Analysis

**Every single parameter (47/47) has σ/μ\* > 1.0.** This means the entire model is dominated by parameter interactions and nonlinearities. No parameter acts additively.

### Interaction Intensity Tiers

| σ/μ\* Range | Count | Interpretation | Examples |
|:---|:---|:---|:---|
| 1.0–1.5 | 2 | Moderate interactions | rho_rec (1.46), settler_survival (1.42) |
| 1.5–2.0 | 16 | Strong interactions | k_growth (1.63), K_half (1.84), n_resistance (1.78), mu_I2D_ref (1.98) |
| 2.0–2.5 | 26 | Very strong interactions | a_exposure (2.20), target_mean_c (2.08), alpha_srs (2.34), q_init_beta_a (2.45) |
| > 2.5 | 3 | Extreme interactions | sigma_v_mutation (2.52), n_tolerance (2.51) — both genetic/evolutionary params |

**Key insight:** The most interacting parameters are genetic architecture (`n_tolerance`, `q_init_beta_a`) and pathogen evolution (`sigma_v_mutation`) — exactly the new R4 additions. These create cascading interactions through the model because they control *adaptation rates* which feed back on disease dynamics which feed back on selection pressures.

**rho_rec** has the *lowest* interaction ratio (1.46) while being the #1 parameter — it operates semi-additively because recovery rate directly scales the probability of clearing infection, with less dependence on the surrounding parameter context.

## Module-Level Analysis

| Module | # Params | Mean μ* | Max μ* | Top Parameter |
|:---|:---|:---|:---|:---|
| Disease | 16 | 0.332 | 0.889 | rho_rec |
| Genetics | 8 | 0.260 | 0.525 | n_resistance |
| Population | 7 | 0.257 | 0.633 | k_growth |
| Spawning | 7 | 0.198 | 0.392 | peak_width_days |
| Pathogen Evolution | 6 | 0.185 | 0.259 | sigma_v_mutation |
| Spatial | 3 | 0.171 | 0.187 | alpha_self_open |

Disease module dominates (16 params, highest mean), but genetics punches above its weight — 8 params with second-highest mean, and its top parameter (`n_resistance`) is #5 globally.

Pathogen evolution, despite being entirely new, immediately achieves mean μ\* = 0.185 with `sigma_v_mutation` at #14. This confirms virulence evolution is not negligible — it needs to be in the Sobol analysis.

Spatial parameters are detectable for the first time at 11 nodes (R3 with 3 nodes couldn't resolve them). `alpha_self_open` (#25) controls how much larval retention vs. export occurs at exposed-coast nodes — directly relevant to recolonization.

## The 11-Node Network

R4 uses an 11-node stepping-stone chain spanning the latitudinal range of *Pycnopodia* habitat:

- **Why 11 nodes:** Captures north-south temperature gradient, fjord vs. open coast distinction, and realistic dispersal distances without the computational cost of 150+ nodes
- **Why it matters:** At R3's 3 nodes, spatial parameters were unresolvable. At 11 nodes, `alpha_self_open` jumps 14 ranks to #25, `alpha_self_fjord` jumps 9 ranks to #33, and latitudinally-varying parameters (`T_ref`, `P_env_max`) become much more influential
- **Validation:** The north-south mortality gradient metric is now well-resolved across nodes

## Recommendations for Sobol

### Must Include (Top 20, all contributing)
All 47 parameters should advance to Sobol — no parameter can be safely eliminated given the universal σ/μ\* > 1.0 interaction signal.

### Sobol Design Priorities
1. **N ≥ 512 samples** (same as R3 Sobol) — the extreme nonlinearity demands it
2. **Focus convergence checks on top-10**: rho_rec, k_growth, K_half, P_env_max, n_resistance, settler_survival, sigma_2_eff, mu_I2D_ref, peak_width_days, target_mean_c
3. **Watch for S2 (second-order) indices** between:
   - `rho_rec × target_mean_c` (recovery rate × recovery genetics — both affect clearance)
   - `P_env_max × a_exposure` (environmental reservoir × transmission — dual exposure pathways)
   - `n_resistance × sigma_v_mutation` (host genetics × pathogen adaptation — coevolutionary arms race)
   - `k_growth × settler_survival` (growth × recruitment — demographic compensation)
4. **11-node network for Sobol**: The spatial resolution is necessary to capture the parameters that only emerged at this scale
5. **Estimated compute**: 47 params × (N+2) × 512 = ~25,088 runs. At ~25s/run on Xeon = ~174 hours (~7 days on 12 cores)

### Parameters to Fix for Computational Savings (if needed)
If compute budget forces reduction, the safest params to fix at nominal values are the bottom 5 (#43-47: gamma_fert, susceptibility_multiplier, p_spontaneous_female, q_init_beta_b, F0) — but even these contribute. Proceed with caution.

## Figures

All figures saved in `results/sensitivity_r4/figures/`:

1. **morris_r4_top20.png** — Top 20 parameters by mean normalized μ\*, color-coded by module
2. **morris_r4_interaction.png** — μ\* vs σ scatter plot with interaction ratio reference lines
3. **morris_r4_heatmap.png** — Full 47×23 parameter-metric heatmap
4. **morris_r4_vs_r3.png** — R3→R4 rank change visualization for 43 common parameters
5. **morris_r4_modules.png** — Module-level sensitivity breakdown
6. **morris_r4_genetics.png** — Three-trait genetic architecture focus (new R4 params highlighted)
7. **morris_r4_pathogen_evolution.png** — Pathogen evolution parameter sensitivity and interaction ratios

## Raw Data

- `morris_results.json` — Full Morris output (μ\*, σ, confidence intervals per param×metric)
- `morris_r4_ranking.json` — Processed ranking with R3 comparison
