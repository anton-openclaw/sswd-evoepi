# SA Parameter Selection: Fix vs. Sweep for Full-Model Sensitivity Analysis

**SSWD-EvoEpi Model — 896-Node Full Spatial Network**
**Date:** March 7, 2026
**Authors:** Anton Star & Willem Weertman

---

## Executive Summary

We recommend **13 free (SWEEP) parameters** and **47 fixed parameters** for the full-model Morris + Sobol sensitivity analysis. This reduction from 60 total parameters to 13 makes the analysis computationally feasible on the Xeon W-3365 cluster (~5 hours for Morris, ~4 days for Sobol).

The 13 SWEEP parameters were selected because they either:
- Ranked in the top tier of SA R4 Sobol sensitivity (ST > 0.03)
- Are new parameters added since SA R4 with unknown sensitivity and high expected influence
- Control biological outcomes central to the six target questions

**Critical insight:** The model has fundamentally changed since SA R4. The addition of dynamic P_env (alpha_env, delta_env, P_env_floor), wavefront spread (CDT, wavefront_D_P), and enclosedness-based flushing (n_connectivity) means the old SA R4 rankings are necessary but not sufficient. The new parameters likely dominate the full-model dynamics — alpha_env alone drove more calibration variation than any original parameter.

---

## 1. SA R4 Results Summary

SA R4 evaluated 47 parameters on an 11-node stepping-stone network using Morris screening (r=50 trajectories, 2,400 runs) and Sobol analysis (N=512, 25,088 runs). Key findings:

### Top 10 by Sobol ST (pop_crash_pct)

| Rank | Parameter | S1 | ST | Module |
|------|-----------|----:|-----:|--------|
| 1 | K_half | 0.182 | 0.456 | Disease |
| 2 | a_exposure | 0.112 | 0.337 | Disease |
| 3 | P_env_max | 0.049 | 0.251 | Disease |
| 4 | sigma_2_eff | 0.041 | 0.232 | Disease |
| 5 | sigma_D | 0.075 | 0.141 | Disease |
| 6 | T_vbnc | −0.015 | 0.040 | Disease |
| 7 | k_growth | −0.009 | 0.035 | Population |
| 8 | peak_width_days | 0.018 | 0.035 | Spawning |
| 9 | settler_survival | 0.010 | 0.033 | Population |
| 10 | target_mean_r | 0.008 | 0.021 | Genetics |

### Bottom quartile (all ST < 0.005)

Parameters ranked 16–47 all had ST < 0.006 for pop_crash_pct. These are safe to fix.

### Morris vs. Sobol discrepancies

- **rho_rec**: Morris #1 (μ\*=0.889) → Sobol #32 (ST=0.004). Morris inflated by interactions.
- **sigma_D**: Morris #20 (μ\*=0.211) → Sobol #5 (ST=0.141). Nonlinear feedback underestimated by OAT.
- **P_env_max**: Morris #4 (μ\*=0.598) → Sobol #3 (ST=0.251). Consistent.

---

## 2. Model Changes Since SA R4

The full model differs from the SA R4 configuration in five major ways:

1. **Dynamic P_env** (replaces static P_env_max): Environmental pathogen pool is now driven by host shedding (alpha_env), natural decay (delta_env), and a community floor (P_env_floor). P_env_max is no longer used.

2. **Wavefront disease spread**: Disease originates in CA-S and spreads north via cumulative dose accumulation (CDT) and long-range pathogen dispersal (wavefront_D_P).

3. **Enclosedness-based flushing**: Per-site flushing rates computed from GIS-derived enclosedness metrics, mapped through a nonlinear function controlled by n_connectivity.

4. **Pathogen thermal adaptation**: T_vbnc evolves per-site toward local SST, bounded by T_vbnc_min = 9°C (literature-constrained biophysical floor).

5. **Community virulence evolution**: Pathogen virulence (v_local) evolves per-site based on host density and temperature, creating an emergent north-south virulence gradient.

These changes mean **P_env_max (SA R4 #3, ST=0.251) is structurally obsolete** — its role is entirely absorbed by alpha_env + delta_env + P_env_floor. Similarly, sigma_D's role is partially subsumed by alpha_env, which controls the fraction of ALL shedding entering the environmental pool.

---

## 3. SWEEP Parameters (13)

### 3.1 Disease Transmission Core (3 parameters)

These dominated SA R4 and will remain influential regardless of model changes.

| # | Parameter | Default | Sweep Range | Scale | Reason |
|---|-----------|---------|-------------|-------|--------|
| 1 | **K_half** | 87,000 | 20,000–200,000 | log | (E) ST=0.456, #1 Sobol R4 |
| 2 | **a_exposure** | 0.75 | 0.30–1.50 | linear | (E) ST=0.337, #2 Sobol R4 |
| 3 | **sigma_2_eff** | 50.0 | 10–250 | log | (E) ST=0.232, #4 Sobol R4 |

**Rationale:** K_half, a_exposure, and sigma_2_eff form the irreducible disease transmission kernel. They jointly define the dose-response curve and shedding feedback loop. SA R4 showed they account for >70% of total-order sensitivity for pop_crash_pct. No amount of model restructuring changes their fundamental role.

### 3.2 Dynamic P_env System (3 parameters)

New parameters that replaced P_env_max. Unknown sensitivity but expected to be dominant.

| # | Parameter | Default | Sweep Range | Scale | Reason |
|---|-----------|---------|-------------|-------|--------|
| 4 | **alpha_env** | 0.10 | 0.01–0.50 | log | (F,G,H) THE gradient driver |
| 5 | **delta_env** | 0.05 | 0.01–0.20 | log | (F,G,H) Pool persistence |
| 6 | **P_env_floor** | 50.0 | 5–500 | log | (F,G,H) Chronic disease pressure |

**Rationale:** alpha_env is the single most important calibration parameter discovered in W01–W134. It controls how much host shedding amplifies the environmental pathogen pool — the mechanism that creates the north-south recovery gradient. When alpha_env is too low, all regions recover equally. When too high, all regions crash to extinction. delta_env determines pathogen persistence (half-life 3.5–70 days across sweep range). P_env_floor sets the minimum chronic disease pressure from the multi-species vibrio community, determining whether southern populations can ever escape disease suppression.

### 3.3 Wavefront Dynamics (2 parameters)

New parameters controlling spatial disease spread timing.

| # | Parameter | Default | Sweep Range | Scale | Reason |
|---|-----------|---------|-------------|-------|--------|
| 7 | **CDT** | 10.0 | 1–100 | log | (F,H) Wavefront activation threshold |
| 8 | **wavefront_D_P** | 150.0 | 50–500 | log | (F,H) Long-range pathogen dispersal |

**Rationale:** CDT (cumulative dose threshold) and wavefront_D_P jointly control crash timing — one of the six biological questions. Higher CDT = more accumulated pathogen dose needed before a node activates = slower wavefront. Larger wavefront_D_P = further pathogen reaches per timestep = faster wavefront. These two parameters directly determine whether the modeled 3-year south-to-north wavefront matches observations (Gravem et al. 2021).

### 3.4 Spatial Structure (1 parameter)

New parameter controlling enclosedness-flushing mapping.

| # | Parameter | Default | Sweep Range | Scale | Reason |
|---|-----------|---------|-------------|-------|--------|
| 9 | **n_connectivity** | 1.0 | 0.3–3.0 | linear | (F,G,H) Flushing nonlinearity |

**Rationale:** n_connectivity controls whether the flushing rate drops sharply at moderate enclosedness (n < 1) or mostly at high enclosedness (n > 1). The flushing literature (Liu et al. 2019, Khangaonkar et al. 2022) suggests n ≈ 0.5–0.7 is most physically realistic, but this is uncertain. Flushing directly modulates pathogen retention time, which affects recovery in enclosed sites like PWS and Hood Canal. This parameter has no SA R4 data because it didn't exist.

### 3.5 Temperature & Demography (4 parameters)

SA R4 second-tier parameters plus the key evolutionary architecture parameter.

| # | Parameter | Default | Sweep Range | Scale | Reason |
|---|-----------|---------|-------------|-------|--------|
| 10 | **T_vbnc** | 12.0 | 8–15 | linear | (E,H) ST=0.040, seasonal disease window |
| 11 | **k_growth** | 0.08 | 0.03–0.15 | linear | (E) ST=0.035, recovery speed |
| 12 | **settler_survival** | 0.001 | 0.005–0.10 | log | (E) ST=0.033, recruitment capacity |
| 13 | **n_resistance** | 17 | 5–30 | linear | (E,H) ST=0.020 (pop), ST=0.275 (evol) |

**Rationale:** T_vbnc determines the seasonal on/off switch for disease — critical for the La Niña rebound question. k_growth controls how fast individuals reach reproductive size after a crash — directly impacts recovery timescale. settler_survival controls recruitment success — the demographic bottleneck for crashed populations. n_resistance controls genetic variance and thus evolutionary speed — essential for the evolutionary dynamics question. Despite its modest ST=0.020 for pop_crash_pct, n_resistance has ST=0.275 for resistance_shift_mean, making it the single most important parameter for evolutionary outcomes.

---

## 4. FIX Parameters (47)

### 4.1 Disease Module — Fixed (13 parameters)

| Parameter | Fixed Value | Reason Code | Justification |
|-----------|-------------|-------------|---------------|
| sigma_1_eff | 5.0 | A | ST=0.016, below threshold |
| sigma_D | 15.0 | A | ST=0.141 in R4, but alpha_env now subsumes role in full model; fix to maintain budget |
| rho_rec | 0.05 | A | ST=0.004; Morris #1 was inflated by interactions |
| mu_EI1_ref | 0.233 | A,B | ST=0.004; Prentice et al. 2025 constraints |
| mu_I1I2_ref | 0.434 | A,B | ST=0.009; temperature-corrected from experiments |
| mu_I2D_ref | 0.563 | A,B | ST=0.014; borderline, but progression rates tightly coupled |
| P_env_max | 500.0 | C | Structurally obsolete — replaced by dynamic P_env system |
| T_ref | 20.0 | A,B | ST=0.004; well-constrained by V. pectenicida literature |
| s_min | 10.0 | A,B | ST=0.004; marine Vibrio salinity requirements |
| susceptibility_multiplier | 2.0 | A | ST=0.003; superseded by explicit genetic resistance |
| immunosuppression_duration | 28 | A | ST=0.005 |
| min_susceptible_age_days | 0 | A | ST=0.010; set to 0 (conservative, backward compatible) |
| k_vbnc | 1.0 | C | Design choice — standard sigmoid steepness |

### 4.2 Disease Module — New Parameters Fixed (4 parameters)

| Parameter | Fixed Value | Reason Code | Justification |
|-----------|-------------|-------------|---------------|
| T_vbnc_min | 9.0 | B | Literature-constrained biophysical floor (Vibrio VBNC <10°C; see specs/vibrio_thermal_limits_review.md) |
| pathogen_adapt_rate (s0) | 0.001 | B,G | Calibration-narrowed; interacts with T_vbnc which is swept |
| v_adapt_rate | 0.001 | G | Community virulence drift rate; secondary to alpha_env for gradient |
| K_cv | 0.0 | C,D | Carrying capacity variability; adds stochasticity but not a key gradient driver |

### 4.3 Spatial Module — Fixed (5 parameters)

| Parameter | Fixed Value | Reason Code | Justification |
|-----------|-------------|-------------|---------------|
| D_L | 400.0 | A | ST=0.003; larval dispersal scale well-constrained by PLD × current speed |
| D_P | 15.0 | A,B | Standard pathogen dispersal; wavefront_D_P handles long-range separately |
| alpha_self_fjord | 0.30 | A | ST=0.004; estuarine self-retention fraction |
| alpha_self_open | 0.10 | A | ST=0.003; open coast self-retention |
| phi_open | 0.80 | C | Boundary condition for enclosedness mapping |
| phi_fjord | 0.03 | C | Boundary condition for enclosedness mapping |

### 4.4 Population Module — Fixed (7 parameters)

| Parameter | Fixed Value | Reason Code | Justification |
|-----------|-------------|-------------|---------------|
| F0 | 1.0×10⁷ | A | ST=0.004; fecundity order of magnitude |
| gamma_fert | 4.5 | A | ST=0.004; fertilization kinetics |
| alpha_srs | 1.35 | A | ST=0.004; size-recruitment survival shape |
| senescence_age | 50.0 | A | ST=0.003; echinoderm negligible senescence |
| L_min_repro | 400.0 | A | ST=0.002; lowest-ranked R4 parameter |
| L_inf | 1000.0 | B | Well-constrained by morphometric data |
| fecundity_exp | 2.5 | B | Allometric scaling |

### 4.5 Genetics Module — Fixed (7 parameters)

| Parameter | Fixed Value | Reason Code | Justification |
|-----------|-------------|-------------|---------------|
| n_tolerance | 17 | A | ST=0.004; partition of 51 loci |
| target_mean_r | 0.15 | A | ST=0.021; moderate sensitivity but degenerate with n_resistance |
| target_mean_t | 0.10 | A | ST=0.004 |
| target_mean_c | 0.02 | A | ST=0.003 |
| tau_max | 0.85 | A | ST=0.005 |
| q_init_beta_a | 2.0 | A | ST=0.004; allele frequency distribution shape |
| q_init_beta_b | 8.0 | A | ST=0.003 |

### 4.6 Spawning Module — Fixed (7 parameters)

| Parameter | Fixed Value | Reason Code | Justification |
|-----------|-------------|-------------|---------------|
| peak_width_days | 60.0 | A | ST=0.035; borderline but spawning timing is secondary |
| p_spontaneous_female | 0.012 | A | ST=0.004 |
| p_spontaneous_male | 0.0125 | A | ST=0.004 |
| female_max_bouts | 2 | A | ST=0.004 |
| induction_female_to_male | 0.80 | A | ST=0.004 |
| induction_male_to_female | 0.60 | A | ST=0.004 |
| readiness_induction_prob | 0.50 | A | ST=0.005 |

### 4.7 Pathogen Evolution Module — Fixed (6 parameters)

| Parameter | Fixed Value | Reason Code | Justification |
|-----------|-------------|-------------|---------------|
| alpha_kill | 2.0 | A,B | ST=0.006; Anderson-May trade-off theory constrains |
| alpha_shed | 1.5 | A,B | ST=0.005; convex trade-off requirement |
| alpha_prog | 1.0 | A,B | ST=0.004; linear scaling assumption |
| gamma_early | 0.3 | A | ST=0.004; asymptomatic transmission fraction |
| sigma_v_mutation | 0.02 | A | ST=0.005; mutation step size |
| v_init | 0.5 | A,B | ST=0.005; ancestral virulence at host-shift |

---

## 5. Borderline Decisions

Three parameters warrant discussion:

### sigma_D (FIXED at 15.0 — borderline SWEEP)
SA R4 ranked sigma_D #5 (ST=0.141). In the old static P_env model, dead-animal shedding created a critical positive feedback loop. With dynamic P_env, alpha_env now controls the fraction of ALL shedding (including from deaths) entering the environmental pool. sigma_D and alpha_env are partially degenerate — both increase pathogen loading from hosts. Sweeping both would create identifiability problems. We fix sigma_D and sweep alpha_env because alpha_env is the more mechanistically fundamental parameter in the new model.

### peak_width_days (FIXED at 60.0 — borderline SWEEP)
SA R4 ranked peak_width_days #8 (ST=0.035). However, it primarily affects spawning synchrony and recruitment timing — it doesn't directly control any of the six biological questions. In calibration, changing peak_width_days had minimal impact on the recovery gradient compared to alpha_env or n_connectivity. Fix at 60.0 (literature-supported 4-month effective spawning season).

### target_mean_r (FIXED at 0.15 — borderline SWEEP)
SA R4 ranked target_mean_r #10 (ST=0.021). It sets the initial population mean resistance, determining how far the population must evolve. However, it's partially degenerate with n_resistance (both control initial genetic variance and evolutionary potential). We sweep n_resistance because it controls the mechanistic architecture, and fix target_mean_r at the calibrated value.

---

## 6. Biological Questions Coverage

| Biological Question | Key SWEEP Parameters | Coverage |
|---------------------|---------------------|----------|
| **Recovery gradient** (AK 50% vs CA-S 0.1%) | alpha_env, delta_env, P_env_floor, n_connectivity, T_vbnc | ✓ Excellent |
| **Crash timing** (wavefront speed) | CDT, wavefront_D_P, K_half, a_exposure | ✓ Excellent |
| **Evolutionary dynamics** (host trait evolution) | n_resistance, k_growth, settler_survival | ✓ Good |
| **Pathogen adaptation** (virulence/thermal evolution) | T_vbnc, alpha_env, delta_env | ○ Moderate — v_adapt_rate fixed |
| **Reintroduction outcomes** (captive-bred persistence) | settler_survival, alpha_env, K_half | ✓ Good |
| **La Niña rebound** (year-9 recovery spike) | T_vbnc, alpha_env, delta_env, k_growth | ✓ Good |

The 13-parameter set covers all six questions. Pathogen adaptation coverage is moderate because v_adapt_rate and pathogen_adapt_rate are fixed, but T_vbnc (the ancestral threshold) and alpha_env (which drives environmental pool and thus selection pressure) capture the key axes.

---

## 7. Computational Cost Estimates

**Platform:** Intel Xeon W-3365, 64 cores (128 threads), 503 GB RAM
**Per-run time:** ~35 min average (range 30–45 min depending on early stopping)

### Phase 1: Morris Screening (verification)

| Parameter | Value |
|-----------|-------|
| Free parameters (k) | 13 |
| Trajectories (r) | 20 |
| Total runs | r × (k+1) = 280 |
| Serial compute | 163 hours |
| Parallel (48 workers) | **~5 hours** |
| Purpose | Verify new parameters rank as expected; confirm old insensitive params remain insensitive |

### Phase 2: Sobol Analysis

| Parameter | Value |
|-----------|-------|
| Free parameters (k) | 13 |
| Base samples (N) | 256 |
| Total runs | N × (2k+2) = 7,168 |
| Serial compute | 4,181 hours |
| Parallel (48 workers) | **~4 days** |
| Purpose | Variance decomposition, S1/ST indices, interaction structure |

### Combined

| | Runs | Wall Time |
|--|------|-----------|
| Morris | 280 | 5 hours |
| Sobol | 7,168 | 4 days |
| **Total** | **7,448** | **~5 days** |

### Alternative: Larger Sobol (higher confidence)

With N=512: 14,336 runs, ~8 days. Recommended if initial N=256 shows wide confidence intervals on top parameters.

---

## 8. Summary Classification Table

### SWEEP Parameters

| # | Parameter | Module | Default | Min | Max | Scale | Reason | SA R4 ST |
|---|-----------|--------|---------|-----|-----|-------|--------|----------|
| 1 | alpha_env | Disease (new) | 0.10 | 0.01 | 0.50 | log | F,G,H | — |
| 2 | delta_env | Disease (new) | 0.05 | 0.01 | 0.20 | log | F,G,H | — |
| 3 | P_env_floor | Disease (new) | 50.0 | 5 | 500 | log | F,G,H | — |
| 4 | K_half | Disease | 87,000 | 20,000 | 200,000 | log | E | 0.456 |
| 5 | a_exposure | Disease | 0.75 | 0.30 | 1.50 | linear | E | 0.337 |
| 6 | sigma_2_eff | Disease | 50.0 | 10 | 250 | log | E | 0.232 |
| 7 | CDT | Disease (new) | 10.0 | 1 | 100 | log | F,H | — |
| 8 | wavefront_D_P | Disease (new) | 150.0 | 50 | 500 | log | F,H | — |
| 9 | n_connectivity | Spatial (new) | 1.0 | 0.3 | 3.0 | linear | F,G,H | — |
| 10 | T_vbnc | Disease | 12.0 | 8 | 15 | linear | E,H | 0.040 |
| 11 | k_growth | Population | 0.08 | 0.03 | 0.15 | linear | E | 0.035 |
| 12 | settler_survival | Population | 0.001 | 0.005 | 0.10 | log | E | 0.033 |
| 13 | n_resistance | Genetics | 17 | 5 | 30 | linear | E,H | 0.020 |

### FIX Parameters (complete list)

| Parameter | Module | Fixed Value | Reason | SA R4 ST |
|-----------|--------|-------------|--------|----------|
| sigma_1_eff | Disease | 5.0 | A | 0.016 |
| sigma_D | Disease | 15.0 | A* | 0.141 |
| rho_rec | Disease | 0.05 | A | 0.004 |
| mu_EI1_ref | Disease | 0.233 | A,B | 0.004 |
| mu_I1I2_ref | Disease | 0.434 | A,B | 0.009 |
| mu_I2D_ref | Disease | 0.563 | A,B | 0.014 |
| P_env_max | Disease | 500.0 | C | 0.251† |
| T_ref | Disease | 20.0 | A,B | 0.004 |
| s_min | Disease | 10.0 | A,B | 0.004 |
| susceptibility_multiplier | Disease | 2.0 | A | 0.003 |
| immunosuppression_duration | Disease | 28 | A | 0.005 |
| min_susceptible_age_days | Disease | 0 | A | 0.010 |
| k_vbnc | Disease (new) | 1.0 | C | — |
| T_vbnc_min | Disease (new) | 9.0 | B | — |
| pathogen_adapt_rate | Disease (new) | 0.001 | B,G | — |
| v_adapt_rate | Disease (new) | 0.001 | G | — |
| K_cv | Population (new) | 0.0 | C,D | — |
| D_L | Spatial | 400.0 | A | 0.003 |
| D_P | Spatial | 15.0 | A,B | — |
| alpha_self_fjord | Spatial | 0.30 | A | 0.004 |
| alpha_self_open | Spatial | 0.10 | A | 0.003 |
| phi_open | Spatial (new) | 0.80 | C | — |
| phi_fjord | Spatial (new) | 0.03 | C | — |
| F0 | Population | 1.0×10⁷ | A | 0.004 |
| gamma_fert | Population | 4.5 | A | 0.004 |
| alpha_srs | Population | 1.35 | A | 0.004 |
| senescence_age | Population | 50.0 | A | 0.003 |
| L_min_repro | Population | 400.0 | A | 0.002 |
| settler_survival | — | — | — | — |
| n_tolerance | Genetics | 17 | A | 0.004 |
| target_mean_r | Genetics | 0.15 | A | 0.021 |
| target_mean_t | Genetics | 0.10 | A | 0.004 |
| target_mean_c | Genetics | 0.02 | A | 0.003 |
| tau_max | Genetics | 0.85 | A | 0.005 |
| q_init_beta_a | Genetics | 2.0 | A | 0.004 |
| q_init_beta_b | Genetics | 8.0 | A | 0.003 |
| peak_width_days | Spawning | 60.0 | A* | 0.035 |
| p_spontaneous_female | Spawning | 0.012 | A | 0.004 |
| p_spontaneous_male | Spawning | 0.0125 | A | 0.004 |
| female_max_bouts | Spawning | 2 | A | 0.004 |
| induction_female_to_male | Spawning | 0.80 | A | 0.004 |
| induction_male_to_female | Spawning | 0.60 | A | 0.004 |
| readiness_induction_prob | Spawning | 0.50 | A | 0.005 |
| alpha_kill | Path. Evol. | 2.0 | A,B | 0.006 |
| alpha_shed | Path. Evol. | 1.5 | A,B | 0.005 |
| alpha_prog | Path. Evol. | 1.0 | A,B | 0.004 |
| gamma_early | Path. Evol. | 0.3 | A | 0.004 |
| sigma_v_mutation | Path. Evol. | 0.02 | A | 0.005 |
| v_init | Path. Evol. | 0.5 | A,B | 0.005 |

\* Borderline — see §5 for discussion.
† P_env_max had ST=0.251 in R4 but is structurally obsolete in the full model (replaced by dynamic P_env).

### Reason Codes

| Code | Meaning |
|------|---------|
| A | SA R4 showed insensitive (low μ\*, low ST) |
| B | Well-constrained by literature |
| C | Structural/design choice |
| D | Computationally necessary |
| E | SA R4 showed highly sensitive |
| F | New parameter, sensitivity unknown |
| G | Poorly constrained by literature |
| H | Controls key biological outcome |

---

## 9. Recommended SA Output Metrics

For the full-model SA, we recommend 8 summary statistics that map to the biological questions:

| Metric | Biological Question | Definition |
|--------|---------------------|------------|
| recovery_gradient | Recovery gradient | AK-PWS recovery − CA-S recovery |
| rmse_log_recovery | Recovery gradient | Log-space RMSE against 8 regional targets |
| wavefront_mae | Crash timing | MAE of regional disease arrival vs. observed |
| final_pop_frac | Overall crash severity | Total population at year 13 / peak population |
| resistance_shift_mean | Evolutionary dynamics | Δr̄ from initial to final |
| v_local_gradient | Pathogen adaptation | Mean v_local(AK) − mean v_local(CA) |
| reintro_persistence | Reintroduction outcomes | Pop at released site 5 years post-release / released |
| la_nina_rebound | La Niña rebound | Pop at year 9 / pop at year 7 (rebound ratio) |

---

*This document is a decision record. The recommendations are opinionated and based on SA R4 data, 134 calibration rounds, and model structural analysis. Revise after Phase 1 Morris results if the full-model sensitivity landscape differs substantially from expectations.*
