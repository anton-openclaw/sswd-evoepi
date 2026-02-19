# SA Round 3 — Morris Screening Analysis

**Date:** 2026-02-19
**Runs:** 850/880 successful (97%), 30 errors (extreme param combos → crashes)
**Config:** 3-node (Sitka/Howe/Monterey), K=5000, 20yr, PE on, 1 substep/day movement

---

## Global Parameter Ranking

Mean normalized μ* across all 20 metrics. Higher = more influential.

| Rank | Parameter | μ*_norm | Category | Change from R1 |
|------|-----------|---------|----------|----------------|
| 1 | **disease.rho_rec** | 0.666 | Disease | ↑ was #14 in Sobol R1! |
| 2 | **spawning.peak_width_days** | 0.544 | Spawning | 🆕 NEW |
| 3 | **population.settler_survival** | 0.490 | Population | ↑ was #19 in Sobol R1 |
| 4 | **disease.a_exposure** | 0.488 | Disease | → was #3 |
| 5 | **genetics.target_mean_r** | 0.449 | Genetics | 🆕 NEW |
| 6 | population.k_growth | 0.378 | Population | ↑ was #13 |
| 7 | disease.sigma_D | 0.315 | Disease | ↑ |
| 8 | population.L_min_repro | 0.309 | Population | → was #10 |
| 9 | **disease.mu_I2D_ref** | 0.308 | Disease | ↓ was #1 in Sobol R1! |
| 10 | disease.P_env_max | 0.299 | Disease | ↑ was #23 |
| 11 | population.F0 | 0.299 | Population | ↑ |
| 12 | disease.K_half | 0.290 | Disease | ↑ |
| 13 | disease.sigma_2_eff | 0.280 | Disease | ↓ was #4 |
| 14 | disease.mu_EI1_ref | 0.276 | Disease | → was #6 |
| 15 | spawning.p_spontaneous_male | 0.269 | Spawning | 🆕 NEW |
| ... | ... | ... | ... | ... |
| 20 | **pathogen_evolution.alpha_kill** | 0.249 | PE | 🆕 NEW |
| 26 | pathogen_evolution.sigma_v_mutation | 0.231 | PE | 🆕 NEW |
| 33 | pathogen_evolution.alpha_prog | 0.207 | PE | 🆕 NEW |
| 37 | pathogen_evolution.v_init | 0.189 | PE | 🆕 NEW |
| 42 | pathogen_evolution.alpha_shed | 0.176 | PE | 🆕 NEW |
| 43 | population.gamma_fert | 0.174 | Population | → (last) |

---

## Category-Specific Top Drivers

### Population Outcomes (crash, recovery, extinction, mortality)
1. **disease.rho_rec** (0.688) — recovery rate dominates population fate
2. **population.settler_survival** (0.611) — demographic baseline
3. **spawning.peak_width_days** (0.550) — ⚠️ spawning timing matters!
4. disease.a_exposure (0.501)
5. genetics.target_mean_r (0.412)

### Evolutionary Outcomes (resistance shift, Va retention, EF1A, rescue index)
1. **disease.rho_rec** (0.691) — recovery window → selection window
2. **genetics.target_mean_r** (0.662) — starting genetic diversity
3. disease.a_exposure (0.528)
4. population.settler_survival (0.453)
5. disease.sigma_D (0.422)

### Spatial Patterns (N-S gradient, fjord protection, node extinction)
1. **population.settler_survival** (0.704) — demographic resilience
2. disease.rho_rec (0.672)
3. **spawning.peak_width_days** (0.615) — season timing × latitude interaction
4. population.k_growth (0.352)
5. disease.a_exposure (0.315)

### Pathogen Evolution (virulence, virulence shift)
1. **spawning.peak_width_days** (1.000!) — virulence evolution is timing-dependent
2. disease.rho_rec (0.827)
3. disease.a_exposure (0.827)
4. population.L_min_repro (0.758)
5. disease.sigma_D (0.691)

### Spawning/Recruitment (participation, recruitment rate)
1. **population.k_growth** (1.000) — growth rate → reproductive maturity
2. population.senescence_age (0.133) — distant second
3. spawning.peak_width_days (0.071)

---

## Key Surprises

### 1. rho_rec went from Sobol R1 #14 to Morris R3 #1
In Round 1 (23 params, no PE, annual mortality), rho_rec dropped from Morris #1 to Sobol #14. Now it's back at #1 in Morris. This could mean:
- The new model features (continuous mortality, spawning overhaul) amplify recovery's importance
- Or Morris is again overweighting it due to extreme perturbation (will Sobol agree?)
- **Key question for Sobol: does rho_rec stay #1 or drop again?**

### 2. spawning.peak_width_days is a top-5 driver and #1 for virulence evolution
A spawning timing parameter driving pathogen virulence evolution was unexpected. Mechanism: narrow spawning peak → synchronized settlement → age-structured epidemic vulnerability. Wide peak → diffuse settlement → disease pressure spread over time. This interacts with PE because epidemic intensity determines selection pressure on pathogen virulence.

### 3. mu_I2D_ref dropped from Sobol R1 #1 to Morris R3 #9
The former king of the sensitivity landscape. Could mean continuous mortality dilutes its impact (deaths spread daily, not lumped). Or Morris just ranks it differently (marginal vs variance). **Sobol will tell us the truth.**

### 4. target_mean_r at #5 — initial genetic diversity is critical
This was new in Round 2 (never Sobol-tested). Starting resistance level determines whether selection can operate. At target_mean_r=0.05, nobody has enough resistance to benefit from recovery; at 0.30, there's a meaningful high-resistance tail.

### 5. PE params are mid-pack, not dominant
alpha_kill is the highest PE param at #20. The pathogen's evolutionary machinery isn't driving broad model behavior — host dynamics dominate. But PE params DO matter for virulence-specific metrics.

### 6. Zero params eliminated (0/43 below 5% threshold)
Every parameter influences at least one metric. Round 1 had the same result (0/23 eliminated). **We're running all 43 through Sobol.**

---

## Sobol Recommendations

### Parameter list: ALL 43 (no elimination)
No parameter falls below the 5% threshold. Morris is a screening tool — it can't eliminate what it can't distinguish. All 43 go to Sobol.

### Range adjustments: None recommended
Ranges look appropriate — 97% success rate suggests no extreme-parameter crashes are systematic.

### Sobol sizing: N=512, 43 params
- Total Saltelli samples: 2 × 43 × (512 + 1) = 44,118 runs
- At 77s/run, 8 cores: **~117 hours (4.9 days)**
- With checkpoint every 500 runs → resumable

### Priority calibration candidates (from Morris)
1. **disease.rho_rec** — #1 overall, zero empirical basis
2. **spawning.peak_width_days** — #2, moderate field data available
3. **population.settler_survival** — #3, very poorly constrained
4. **disease.a_exposure** — #4, lab-only estimates
5. **genetics.target_mean_r** — #5, bounded by Schiebelhut but wide range

---

## Comparison: Round 1 Morris vs Round 3 Morris

| R1 Rank | R1 Parameter | R3 Rank | Change |
|---------|-------------|---------|--------|
| 1 | settler_survival | 3 | ↓2 |
| 2 | rho_rec | 1 | ↑1 |
| 3 | mu_I2D_ref | 9 | ↓6 |
| 4 | F0 | 11 | ↓7 |
| 5 | a_exposure | 4 | ↑1 |
| — | peak_width_days | 2 | 🆕 |
| — | target_mean_r | 5 | 🆕 |

The landscape shifted: disease rate params lost ground, reproduction/genetics gained.
