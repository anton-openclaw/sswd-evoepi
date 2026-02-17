# SSWD-EvoEpi Sensitivity Analysis Report
## February 17, 2026

### Executive Summary

We ran a two-phase global sensitivity analysis (Morris screening → Sobol variance decomposition) on the SSWD-EvoEpi model with 23 uncertain parameters. The analysis used a 3-node spatial network (Sitka, Howe Sound, Monterey; K=5000 each; ~15,000 total agents) with movement enabled at 1 substep/day over 20 simulated years.

**The single most important finding: disease progression speed (mu_I2D_ref) dominates model behavior more than any other parameter.** How fast late-stage wasting kills an infected star is the #1 driver of evolutionary outcomes, genetic diversity loss, and node extinction risk.

**Second key finding: massive parameter interactions.** Total-order indices (ST) far exceed first-order indices (S1) across most parameters and metrics, meaning the model's behavior is driven by *combinations* of parameters acting together — not by any single parameter in isolation. This has profound implications for calibration strategy.

---

### Methods

| Phase | Method | Runs | Cores | Runtime |
|-------|--------|------|-------|---------|
| 1 | Morris screening (20 trajectories) | 480 | 8 | 16 min |
| 2 | Sobol/Saltelli (N=256) | 12,288 | 8 | 10.5 hrs |

**Spatial configuration:** 3 nodes — Sitka (subarctic, 8.5°C), Howe Sound/Porteau Cove (fjord, 10°C), Monterey (warm temperate, 13°C). K=5000 per site based on Willem's field density estimates (0.1–0.2 ind/m²). Movement: 1 substep/day.

**14 output metrics tracked:** Population crash %, final population fraction, recovery (boolean), extinction, peak mortality, time to nadir, total disease deaths, mean/max resistance shift, additive variance retention, EF1A shift, nodes extinct, N-S mortality gradient, fjord protection effect.

**23 parameters varied** across disease (13), population (7), genetics (1), and spawning (2) modules. EF1A parameters excluded (Pisaster finding, not Pycnopodia).

---

### Phase 1: Morris Screening

**Result: All 23 parameters survived screening.** No dimensionality reduction was possible — every parameter influenced at least one metric at ≥5% of the maximum μ* for that metric.

**Morris top drivers:**
- **Population outcomes:** settler_survival (μ*=26.0), rho_rec (20.3), mu_I2D_ref (16.1)
- **Evolutionary outcomes:** rho_rec (μ*=0.100), a_exposure (0.034), mu_EI1_ref (0.026)

**Key lesson:** Morris and Sobol rankings **disagree** in important ways. Morris flagged settler_survival and rho_rec as the top population drivers; Sobol found susceptibility_multiplier and mu_I2D_ref. This is because Morris measures *marginal* effects (one-at-a-time perturbations from extreme values), while Sobol captures *variance* across the full parameter space including interactions. **Morris alone would have led us to focus on the wrong parameters.**

---

### Phase 2: Sobol Variance Decomposition

#### Global Parameter Ranking (Mean ST across all 14 metrics)

| Rank | Parameter | Mean ST | Category | Confidence |
|------|-----------|---------|----------|------------|
| 1 | **mu_I2D_ref** (I2→death rate) | 0.638 | Disease | ★☆☆ |
| 2 | **susceptibility_multiplier** | 0.540 | Disease | ★☆☆ |
| 3 | **a_exposure** (exposure rate) | 0.473 | Disease | ★☆☆ |
| 4 | **sigma_2_eff** (late shedding) | 0.456 | Disease | ★☆☆ |
| 5 | **n_additive** (# resistance loci) | 0.431 | Genetics | ★★☆ |
| 6 | mu_EI1_ref (E→I1 rate) | 0.403 | Disease | ★★☆ |
| 7 | K_half (half-saturation) | 0.387 | Disease | ★☆☆ |
| 8 | T_vbnc (VBNC threshold) | 0.385 | Disease | ★★☆ |
| 9 | s_min (min susceptibility) | 0.384 | Disease | ★☆☆ |
| 10 | L_min_repro (min repro size) | 0.381 | Population | ★★★ |
| ... | ... | ... | ... | ... |
| 14 | **rho_rec** (recovery rate) | 0.370 | Disease | ★☆☆ |
| 19 | settler_survival | 0.339 | Population | ★☆☆ |
| 23 | P_env_max | 0.293 | Disease | ★☆☆ |

#### Metric-Specific Top Drivers

**Population outcomes (crash, recovery, extinction):**
- `susceptibility_multiplier` is #1 for whether populations crash and whether they recover
- `mu_I2D_ref` is #1 for how many nodes go extinct
- `senescence_age` surprisingly important (#3 for crash/recovery) — background mortality baseline matters

**Evolutionary outcomes (resistance shift, Va retention):**
- `mu_I2D_ref` dominates (#1 for resistance shift, Va retention, EF1A shift)
- `n_additive` is #2 for mean resistance shift — genetic architecture matters hugely
- `a_exposure` is #2 for Va retention — exposure rate determines selection intensity

**Spatial patterns (N-S gradient, fjord protection):**
- `a_exposure` is #1 for fjord protection effect
- `susceptibility_multiplier` is #1 for N-S mortality gradient
- These are the most interaction-dominated metrics (ST >> S1 everywhere)

---

### Key Scientific Insights

#### 1. Speed of Death > Probability of Recovery

rho_rec (recovery rate), which Morris flagged as the #1 evolutionary driver, dropped to #14 in the Sobol ranking. Meanwhile mu_I2D_ref (how fast late-stage wasting kills) rose to #1 overall. 

**Interpretation:** The rate of disease progression through the lethal I2 stage sets the *time window* for recovery to operate. A fast-killing disease leaves no opportunity for immune clearance, regardless of the recovery rate. This is consistent with lab observations of rapid progression once wasting symptoms appear (days to weeks, not months).

#### 2. Massive Parameter Interactions

For most metrics, ST >> S1. Examples:
- **Extinction:** sigma_2_eff ST=1.51, S1≈0 → extinction is *entirely* driven by interactions
- **Fjord protection:** a_exposure ST=0.96, S1=-0.12 → negative S1 means the parameter's effect *reverses sign* depending on other parameters
- **Recovery:** susceptibility_multiplier ST=0.96, S1=0.38 → 60% of its influence comes through interactions

**Implication:** The model cannot be calibrated one parameter at a time. Parameters form coupled feedback loops (shedding × exposure × death rate = epidemic intensity; resistance × recovery × reproduction = evolutionary potential). Calibration must be done jointly, ideally with ABC or MCMC.

#### 3. Genetic Architecture is a Structural Choice, Not a Detail

n_additive (number of resistance loci) ranked #5 globally and #2 for mean resistance shift. This isn't a parameter we can calibrate from data — it's a model design choice with major consequences. The planned structural comparison (10 vs 20 vs 51 loci) is essential.

#### 4. Population Parameters Matter More Than Expected

senescence_age (#17 globally) is #3 for population crash and recovery. L_min_repro (#10) and k_growth (#13) are in the top 15. Background demography sets the stage that disease dynamics play out on — faster-growing populations with earlier reproduction buffer disease impacts better.

#### 5. Spawning Parameters Have Targeted Effects

p_spontaneous_female is #1 for time_to_nadir (controls epidemic speed via population turnover timing) but low for other metrics. induction_female_to_male is #2 for node extinction. Spawning mechanics don't drive the big-picture outcomes but matter for specific spatial and temporal patterns.

---

### Implications for Next Steps

1. **Priority calibration targets:** mu_I2D_ref, susceptibility_multiplier, a_exposure, sigma_2_eff — these 4 disease parameters drive the most variance and have the lowest empirical confidence (★☆☆)

2. **Structural comparisons needed:**
   - Genetic architecture: 10 vs 20 vs 51 additive loci
   - Beta vs uniform allele frequency initialization
   - With/without EF1A overdominant locus
   - **Pathogen evolution** (co-evolutionary dynamics)

3. **Joint calibration:** Given the massive interactions, we need MCMC/ABC calibration against empirical targets (Hamilton 2021 decline rates, Schiebelhut allele frequency shifts, monitoring data), not parameter-by-parameter tuning

4. **New parameters to add to SA:** 
   - `alpha_self_fjord` / `alpha_self_open` (larval retention — just made configurable)
   - `q_init_beta_a`, `q_init_beta_b` (allele frequency distribution shape)
   - `target_mean_r` (initial mean resistance)

---

### Model Changes Made Today (Feb 17)

1. **Beta-distributed allele frequencies** — Each locus now gets an independent random starting frequency drawn from Beta(2, 8), producing realistic among-locus variation. Old uniform-q mode available as fallback.

2. **Cause-of-death tracking** — New `DeathCause` enum (DISEASE, NATURAL, SENESCENCE) stamped on every agent at death. Per-site per-year tracking in spatial sim.

3. **Configurable larval retention** — `alpha_self_fjord` (0.30) and `alpha_self_open` (0.10) now tunable in config.

4. **target_mean_r raised to 0.15** — Supported by Schiebelhut data requiring pre-outbreak q≈0.15–0.25 at immune loci to produce observed Δq=0.08–0.15.

5. **Found hidden bug:** Old q_init=0.05 didn't account for EF1A's ~0.058 contribution to mean r, so actual initial resistance was ~0.11, not 0.05.

---

### Files

- Morris results: `results/sensitivity/morris_screening.json`
- Sobol indices: `results/sensitivity/sobol_indices.json`
- Raw Sobol data: `results/sensitivity/sobol_batch0_results.npz` (12,288 × 14)
- Plots: `results/sensitivity/report/*.png`
- This report: `results/sensitivity/report/SENSITIVITY_ANALYSIS_REPORT.md`

---

*Analysis by Anton 🔬 | SSWD-EvoEpi v0.2 | 12,768 total simulation runs | ~11 hours compute on 8× AMD Ryzen 7 5700*
