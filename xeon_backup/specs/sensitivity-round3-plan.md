# Sensitivity Analysis Round 3 — Full Plan

**Date:** 2026-02-19
**Authors:** Willem Weertman & Anton 🔬
**Status:** READY FOR REVIEW

---

## Context

The model has undergone 7 major upgrades since SA Round 1 (Feb 16-17):

| Feature | Impact on SA |
|---------|-------------|
| **Pathogen evolution** | 6 new params, 3 new metrics, co-evolutionary dynamics |
| **Continuous settlement** | Smooth recruitment, PLD-dependent timing |
| **Spawning overhaul** | Multi-bout, cascade induction, readiness induction, immunosuppression |
| **Juvenile immunity** | 1 new param, disease-free settlement window |
| **Beta genetics init** | 3 new params, per-locus allele frequency variation |
| **Continuous mortality** | Daily death/growth, no year-boundary artifacts |
| **Configurable larval retention** | 2 new params, fjord-specific self-recruitment |

SA Round 2 was killed twice (virulence recording bug, then continuous settlement made results invalid). **Round 3 is the first complete SA on the stabilized model.**

---

## Parameters: 43 total

### Existing 39 (from param_spec.py)

**Disease (13):** a_exposure, K_half, sigma_1_eff, sigma_2_eff, sigma_D, rho_rec, mu_EI1_ref, mu_I1I2_ref, mu_I2D_ref, P_env_max, T_ref, T_vbnc, s_min

**Disease-spawning coupling (2):** susceptibility_multiplier, immunosuppression_duration

**Population (7):** F0, gamma_fert, settler_survival, alpha_srs, senescence_age, k_growth, L_min_repro

**Genetics (4):** n_additive, target_mean_r, q_init_beta_a, q_init_beta_b

**Spawning (3):** p_spontaneous_female, induction_female_to_male, induction_male_to_female

**Spatial (3):** D_L, alpha_self_fjord, alpha_self_open

**Pathogen evolution (6):** alpha_kill, alpha_shed, alpha_prog, gamma_early, sigma_v_mutation, v_init

**Juvenile immunity (1):** min_susceptible_age_days

### New params to add: 4

| Parameter | Range | Dist | Rationale |
|-----------|-------|------|-----------|
| `spawning.p_spontaneous_male` | [0.005, 0.025] | uniform | Separate from female rate; controls male multi-bout frequency |
| `spawning.peak_width_days` | [30, 90] | uniform | Controls spawning season shape — tight peak vs diffuse. Major lever on when larvae arrive relative to disease. |
| `spawning.readiness_induction_prob` | [0.1, 0.8] | uniform | Social spawning trigger strength — scales cascade participation |
| `spawning.female_max_bouts` | [1, 3] | discrete [1, 2, 3] | Female reproductive investment — 1 bout (conservative) vs 3 (maximum) |

**Total: 43 parameters**

### Parameters NOT included (and why)

| Parameter | Reason for exclusion |
|-----------|---------------------|
| `male_max_bouts` | Co-varies with female_max_bouts; male multi-bout already captured by p_spontaneous_male |
| `male_refractory_days` | Currently 0 (calibrated); low uncertainty |
| `cascade_window` | 3 days is well-constrained by chemical cue persistence biology |
| `cascade_radius` | Tightly coupled to habitat_area; already implicitly varied by D_L |
| `gravity_strength/range` | Movement already dominant — gravity is a fine detail within it |
| `season_start_doy / season_end_doy` | Well-constrained by field observations (Nov-Jul) |
| `EF1A s_het / q_ef1a_init` | Pisaster finding, not Pycnopodia. Structural comparison, not SA. |

---

## Output Metrics: 20

### Existing 17

**Population (7):** pop_crash_pct, final_pop_frac, recovery, extinction, peak_mortality, time_to_nadir, total_disease_deaths

**Evolution (4):** resistance_shift_mean, resistance_shift_max, va_retention_mean, ef1a_shift_mean

**Spatial (3):** n_extinct_nodes, north_south_mortality_gradient, fjord_protection_effect

**Pathogen evolution (3):** mean_final_virulence, virulence_shift, disease_death_fraction

### New metrics: 3

| Metric | Definition | Rationale |
|--------|-----------|-----------|
| `spawning_participation` | Fraction of adult females that spawn at least once per year (averaged over pre-disease years) | Validates spawning overhaul is functioning; sensitive to p_spontaneous_female, peak_width_days, readiness_induction |
| `mean_recruitment_rate` | Mean annual recruits / K (averaged over pre-disease years) | Integrates spawning + settlement + juvenile survival; sensitive to spawning, retention, juvenile immunity |
| `evolutionary_rescue_index` | final_pop_frac × resistance_shift_mean | Composite: captures whether resistance evolution actually translates to population persistence (not just allele frequency change in tiny remnants) |

---

## Spatial Configuration

**3 nodes** (same as Round 1 — proven, interpretable):
- **Sitka** (57°N, 8.5°C, open coast) — cold, high-latitude
- **Howe Sound** (49.5°N, 10°C, fjord) — refugium test
- **Monterey** (36.6°N, 13°C, open coast) — warm, southern range

**K = 5000 per node** (15,000 total). Field densities of 0.1–0.2 ind/m².

**Movement:** 1 substep/day (captures spatial mixing; 24 substeps = 2× runtime for marginal SA benefit).

**Disease year:** 3 (years 0-2 = spinup to demographic equilibrium).

**Simulation:** 20 years.

---

## Execution Plan

### Phase 1: Infrastructure update (~2 hours, 1 cron job)

1. Add 4 new params to `param_spec.py` with ranges and `sample_to_config_overrides` mappings
2. Add 3 new metrics to `spatial_runner.py`
3. Fix the existing `run_single_spatial` TypeError (str + int from Feb 18)
4. Smoke test: 3 runs with random samples, verify all 20 metrics populate
5. Commit + push

### Phase 2: Morris screening (~3 hours compute, 1 cron job)

- **43 params × 20 trajectories + 20 = 880 runs**
- 8 cores, ~77s/run → ~880 × 77 / 8 = **2.4 hours**
- Purpose: eliminate low-influence params for Sobol (goal: reduce to 25-30)
- Threshold: μ* < 5% of max μ* for ANY metric → eliminated
- Output: `results/sensitivity_r3/morris_screening.json`

### Phase 3: Morris analysis + email (30min, 1 cron job) → DECISION POINT

1. Which params survived, which eliminated
2. Morris rankings per metric category (population, evolution, spatial, PE)
3. Comparison with Round 1 Morris (what changed?)
4. Flag surprises (new params ranking high, old top-5 dropping)
5. Recommend final param list + any range adjustments for Sobol
6. Email full results to Willem
7. **WAIT for Willem's go/no-go before launching Sobol**

### Phase 4: Sobol N=512 (3-5 days compute, 1 cron job)

Based on Morris results:

| If Morris cuts to... | Runs | Wall time (8 cores) |
|---------------------|------|---------------------|
| 30 params | 31,744 | 84h (3.5 days) |
| 35 params | 36,864 | 98h (4 days) |
| 40 params | 41,984 | 112h (4.7 days) |

**N=512** for tighter Sobol indices and reliable interaction estimates.
**Bootstrap CIs** (1000 resamples) on all indices.

**Execution:**
- Run with checkpoint every 500 runs
- `setsid` + `disown` for process persistence (lesson from Round 1)
- Monitor progress via file-based checkpoints

### Phase 5: Sobol analysis (~1 hour, 1 cron job)

1. Compute Sobol S1, ST with 95% bootstrap CIs for all surviving params × 20 metrics
2. Round 1 vs Round 3 comparison (how did new features change the landscape?)
3. Identify top drivers per metric category
4. Flag interaction structure changes
5. Do PE params, spawning params, juvenile immunity matter?

### Phase 6: Report + visualizations (~1 hour, 1 cron job)

1. Full visualization library figures:
   - Parameter ranking bars (mean ST)
   - S1 vs ST heatmap (interaction detection)
   - Round 1 vs Round 3 comparison panels
   - Metric-specific driver plots
   - New-param-specific breakdowns
2. Comprehensive report: `results/sensitivity_r3/report/SA_ROUND3_REPORT.md`
3. Email report to Willem
4. WhatsApp summary

---

## Key Scientific Questions for Round 3

1. **Do the Round 1 top 5 still dominate?** (mu_I2D_ref, susceptibility_multiplier, a_exposure, sigma_2_eff, n_additive)
2. **Does pathogen evolution change the sensitivity landscape?** If PE params rank high, co-evolution fundamentally alters model behavior. If they rank low, host-only dynamics dominate.
3. **How important is spawning biology?** With the full overhaul, do spawning params drive population outcomes more than Round 1 suggested?
4. **Does juvenile immunity matter?** min_susceptible_age_days could be a major lever for recovery if it protects settlers from epidemic exposure.
5. **Do continuous mortality + settlement change the interaction structure?** Round 1 found massive interactions (ST >> S1). Are these preserved, amplified, or reduced with smoother dynamics?
6. **What's the evolutionary rescue index sensitive to?** This composite metric is the bottom-line question: can evolution save populations?

---

## Schedule

| Phase | Start | Duration | Dependencies |
|-------|-------|----------|-------------|
| **1: Infrastructure** | Now | 2h | None |
| **2: Morris** | After Phase 1 | 3h | Phase 1 |
| **3: Morris analysis + email** | After Phase 2 | 30min | Phase 2 |
| ⏸️ | **Wait for Willem's approval** | — | Phase 3 |
| **4: Sobol N=512** | After approval | 84-112h | Phase 3 + approval |
| **5: Sobol analysis** | After Phase 4 | 1h | Phase 4 |
| **6: Report + viz** | After Phase 5 | 1h | Phase 5 |

**Total: ~90-118 hours compute** depending on Morris screening effectiveness.

**Earliest completion:** Monday Feb 23 evening (if Morris eliminates ≥10 params).
**Latest completion:** Wednesday Feb 25 (if minimal elimination).

---

## Risk Mitigation

| Risk | Mitigation |
|------|-----------|
| Process killed (OOM, crash) | Checkpoint every 500 runs; resume from last checkpoint |
| Network outage | `setsid` + `disown`; no dependency on WhatsApp/gateway |
| SA runner bug | Smoke test 3 runs before launch |
| Extreme params → instant crash | Return NaN metrics (already handled); track error rate |
| Round takes >3 days | Fall back to N=128 (26h); still scientifically useful |

---

## Differences from Round 1

| Aspect | Round 1 | Round 3 |
|--------|---------|---------|
| Parameters | 23 | 43 |
| Metrics | 14 | 20 |
| PE enabled | No | Yes |
| Settlement | Annual pulse | Continuous (PLD-based daily) |
| Spawning | Simple | Multi-bout, cascade, readiness induction |
| Mortality | Annual lump-sum | Daily continuous |
| Growth | Annual lump-sum | Daily VB |
| Genetics init | Uniform q=0.05 | Beta(a,b) per-locus |
| Juvenile immunity | None | Configurable refractory period |
| Movement substeps | 1/day | 1/day (same) |
| Per-run time | ~24s | ~77s |
| Sobol N | 256 | **512** |
| Bootstrap CIs | None | 1000 resamples |
