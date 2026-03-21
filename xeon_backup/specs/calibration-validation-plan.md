# SSWD-EvoEpi Calibration & Validation Plan

**Date:** 2026-02-19  
**Authors:** Anton (AI) & Willem Weertman  
**Status:** Draft — awaiting review

---

## 1. Overview

This document formalizes the workflow from sensitivity analysis through calibration, convergence validation, and production runs. The core principle: **calibrate cheap, validate expensive**.

| Phase | Purpose | K per node | Runs | Est. Time (8 cores) |
|-------|---------|-----------|------|---------------------|
| 1. Sensitivity Analysis | Identify influential parameters | 5,000 | ~45,000 | ~4 days |
| 2. Calibration (ABC-SMC) | Fit parameters to empirical data | 5,000 | 10,000–50,000 | 1–5 days |
| 3. Convergence Validation | Verify N-independence | 1K–100K | ~350 | ~12 hours |
| 4. Scale Correction | Fix N-dependent parameters | 50,000 | ~200 | ~2 days |
| 5. Production Scenarios | Paper-quality results | 50,000–100,000 | ~500 | 2–5 days |

**Total compute:** ~2–3 weeks at 8 cores. Scales linearly with more cores.

---

## 2. Phase 1: Sensitivity Analysis (IN PROGRESS)

**Status:** Sobol N=512 running, ETA ~Feb 24–25.

### Outputs
- First-order (S₁) and total-order (S_T) Sobol indices with bootstrap CIs
- Parameter ranking by metric category
- Identification of interaction structure (S_T − S₁)
- Parameters to focus calibration on (top 10–15 by S_T)

### Decision Point
After Phase 1, we select the **calibration parameter set**: parameters with S_T > 0.05 for at least one target metric. Expected: 15–25 of the 43 parameters. The rest are fixed at literature values or midpoints.

---

## 3. Phase 2: Calibration via ABC-SMC

### 3.1 Method

**Approximate Bayesian Computation with Sequential Monte Carlo** (ABC-SMC, Toni et al. 2009). Chosen because:
- No likelihood function needed (our model has no closed-form likelihood)
- Handles multimodal posteriors
- Naturally produces uncertainty estimates
- Well-established for ecological IBMs

### 3.2 Calibration Targets

Empirical data to fit against:

| Target | Source | Value | Metric in Model |
|--------|--------|-------|-----------------|
| Population crash (Pycnopodia) | Harvell et al. 2019, Montecino-Latorre et al. 2016 | 80–99% decline | `pop_crash_pct` |
| North–south mortality gradient | Hamilton et al. 2021 | Southern populations hit harder | `north_south_mortality_gradient` > 0 |
| Fjord/deep-water refugia | Hamilton et al. 2021 | Higher survival in semi-enclosed waters | `fjord_protection_effect` > 0 |
| Allele frequency shift at immune loci | Schiebelhut et al. 2018 | Δq = 0.08–0.15 | `resistance_shift_mean` |
| Timeline: onset to nadir | Montecino-Latorre et al. 2016 | 2–5 years post-introduction | `time_to_nadir` |
| Wasting progression time | Hewson et al. 2014, Kohl et al. 2016 | ~7–21 days from symptoms to death | (implicit in disease rates) |
| Outplanting survival (if available) | Sunflower Star Lab 2025 | 47/48 survived 4 weeks | (conservation module, Phase 5) |

### 3.3 Summary Statistics

Distance function between observed and simulated:

```
d(sim, obs) = Σ_i w_i × |S_i(sim) - S_i(obs)| / σ_i
```

Where:
- S_i = summary statistic (e.g., pop_crash_pct, resistance_shift_mean)
- σ_i = normalization (range or empirical SD)
- w_i = weight (higher for better-constrained targets)

Proposed weights:

| Statistic | Weight | Justification |
|-----------|--------|---------------|
| pop_crash_pct | 1.0 | Well-documented, narrow range |
| resistance_shift_mean | 1.0 | Schiebelhut quantitative data |
| time_to_nadir | 0.5 | Approximate (2–5 yr window) |
| north_south_gradient (sign) | 0.5 | Qualitative constraint |
| fjord_protection (sign) | 0.5 | Qualitative constraint |

### 3.4 ABC-SMC Protocol

1. **Prior:** Uniform over SA ranges (same as Sobol bounds) for calibration parameters; fixed for non-influential parameters
2. **Particles:** N_particles = 1,000
3. **Populations:** T = 5–8 sequential populations (adaptive ε schedule)
4. **Acceptance threshold:** Start at ε₁ = 75th percentile of prior-predictive distances, shrink by 50% each population until ε stabilizes
5. **Perturbation kernel:** Component-wise uniform with adaptive width (Beaumont et al. 2009)
6. **Stopping:** When acceptance rate < 1% or ε change < 5% between populations
7. **Seeds per evaluation:** 3 (median distance across seeds to reduce stochastic noise)

### 3.5 Computational Cost

- Per evaluation: 3 seeds × ~77s = ~230s
- Per population: 1,000 accepted particles / acceptance_rate
  - Population 1 (loose ε): ~3,000 evaluations → ~8 hours on 12 cores
  - Population 5 (tight ε): ~10,000 evaluations → ~27 hours on 12 cores
- **Total estimate: 10,000–50,000 evaluations → 1–5 days on Ryzen, 6–24 hours on Xeon**

### 3.6 Outputs

- Posterior distributions for each calibrated parameter
- Joint posterior (correlation structure between parameters)
- Posterior predictive check: simulate 1,000 parameter sets from posterior, compare to observations
- Parameter identifiability assessment: which parameters are well-constrained vs. still diffuse?

### 3.7 Software

Use `pyABC` (Klinger et al. 2018) — mature Python ABC-SMC library with:
- Built-in parallelization (multiprocessing, Redis, Dask)
- Adaptive ε scheduling
- Visualization tools
- SQLite database for results

---

## 4. Phase 3: Convergence Validation

### 4.1 Purpose

Verify that calibrated parameters produce consistent dynamics at larger population sizes. Tests the assumption that per-capita rates are N-independent.

### 4.2 Protocol

Fix all parameters at posterior median from Phase 2. Run:

| K per node | Total agents (3 nodes) | Seeds | Metrics tracked |
|-----------|----------------------|-------|----------------|
| 1,000 | 3,000 | 10 | All 20 SA metrics |
| 2,500 | 7,500 | 10 | All 20 |
| 5,000 | 15,000 | 10 | All 20 |
| 10,000 | 30,000 | 10 | All 20 |
| 25,000 | 75,000 | 10 | All 20 |
| 50,000 | 150,000 | 10 | All 20 |
| 100,000 | 300,000 | 10 | All 20 |

**Total: 70 runs** (feasible in ~12 hours on Ryzen including the large-K runs).

### 4.3 Convergence Criteria

For each metric, compute across the 10 seeds at each K:
- **Mean:** Should stabilize (< 5% change between K doublings)
- **CV (coefficient of variation):** Should decrease monotonically
- **CV threshold:** CV < 5% for population metrics, CV < 10% for evolutionary metrics

### 4.4 Expected Results

| Metric category | Expected convergence K | Rationale |
|----------------|----------------------|-----------|
| Population (crash, recovery) | ~5,000–10,000 | Demographic noise ~ 1/√N |
| Spatial (gradient, fjord) | ~10,000–25,000 | Per-node stochasticity |
| Evolutionary (Δr̄, Va) | ~25,000–50,000 | Ne = K × 10⁻³, need Ne×s > 1 |
| Pathogen evolution (Δv) | ~10,000–25,000 | Pathogen has larger Ne (no SRS) |

### 4.5 Decision Point

- If all metrics converge by K = 50,000: **use K = 50,000 for production**
- If evolutionary metrics need K = 100,000: **use K = 100,000 for production**
- If convergence requires K > 100,000: reconsider SRS implementation (Ne/N ratio may be too extreme)

---

## 5. Phase 4: Scale Correction

### 5.1 Purpose

If Phase 3 reveals systematic bias (e.g., crash is 95% at K=5K but 92% at K=50K), identify which parameters need adjustment and apply a correction.

### 5.2 Method

1. Compare posterior predictive distributions at K=5K vs K=50K for each calibration target
2. For each target with significant bias (|Δmean| > 2σ):
   - Identify the most influential parameter for that metric (from SA)
   - Adjust that parameter to minimize bias
3. Re-validate: confirm correction doesn't break other metrics

### 5.3 Expected Scale-Dependent Parameters

| Parameter | Why N-dependent | Expected direction |
|-----------|----------------|-------------------|
| a_exposure | Epidemic stochasticity at small N biases toward stronger exposure | Decrease slightly at large N |
| settler_survival | Beverton-Holt density-dependence relative to K | May need rescaling |
| rho_rec | Recovery interacts with stochastic fadeout | Minimal correction expected |

### 5.4 Computational Cost

~200 runs at K=50K (~2 days on Ryzen).

---

## 6. Phase 5: Production Scenarios

### 6.1 Scenario Design

Using calibrated + scale-corrected parameters at K = 50,000–100,000 per node.

#### 6.1.1 Baseline Scenarios

| Scenario | Nodes | Config | Seeds | Purpose |
|----------|-------|--------|-------|---------|
| A. No intervention | 489 | Full coastline, disease at yr 3 | 20 | Reference trajectory |
| B. Status quo | 489 | Same as A, 100-year horizon | 20 | Long-term extinction/recovery |
| C. Natural recovery only | 150 | Reduced coast, no captive breeding | 20 | Can evolution alone save them? |

#### 6.1.2 Conservation Scenarios

| Scenario | Config | Seeds | Question |
|----------|--------|-------|---------|
| D. Captive outplanting (Monterey protocol) | 50 captive-bred juveniles/yr at 5 sites | 20 | Does current outplanting scale matter? |
| E. Aggressive outplanting | 500/yr at 20 sites | 20 | What scale is needed? |
| F. Genetic rescue | Outplant from high-resistance broodstock | 20 | Does selective breeding accelerate recovery? |
| G. Assisted gene flow | Move individuals between refugia | 20 | Does connectivity management help? |
| H. Pathogen management | Reduce environmental Vibrio at key sites | 20 | Is environmental intervention viable? |

#### 6.1.3 Sensitivity Scenarios

| Scenario | Config | Seeds | Question |
|----------|--------|-------|---------|
| I. Genetic architecture | 10 vs 20 vs 51 loci | 10 each | Does architecture change qualitative outcomes? |
| J. Climate warming (+2°C) | SST shifted | 20 | How does warming affect disease-evolution dynamics? |
| K. No pathogen evolution | PE disabled | 20 | How much does co-evolution matter? |
| L. No SRS | α_srs = ∞ (standard reproduction) | 20 | Is sweepstake reproduction important? |

### 6.2 Computational Cost

| Group | Runs | K | Est. time (8 cores) |
|-------|------|---|---------------------|
| Baselines (A–C) | 60 | 100K | ~2 days |
| Conservation (D–H) | 100 | 50K | ~2 days |
| Sensitivity (I–L) | 90 | 50K | ~1.5 days |
| **Total** | **250** | — | **~5 days** |

### 6.3 Outputs

For each scenario:
- Population trajectory (mean ± 95% CI across seeds)
- Resistance evolution trajectory
- Virulence evolution trajectory
- Spatial maps (per-node survival at year 20, 50, 100)
- Extinction probability
- Time to recovery (if applicable)
- Comparison to baseline (effect size + significance)

---

## 7. Timeline

| Week | Phase | Status |
|------|-------|--------|
| Feb 19–25 | Phase 1: Sobol SA | IN PROGRESS |
| Feb 25–26 | Phase 1: Sobol analysis + report | — |
| Feb 26 – Mar 2 | Phase 2: ABC-SMC calibration | — |
| Mar 2–3 | Phase 3: Convergence validation | — |
| Mar 3–4 | Phase 4: Scale correction (if needed) | — |
| Mar 4–8 | Phase 5: Production scenarios | — |
| Mar 8–14 | Analysis, visualization, writing | — |

---

## 8. Software Dependencies

| Package | Purpose | Status |
|---------|---------|--------|
| SALib | Sensitivity analysis (Morris, Sobol) | ✅ Installed |
| pyABC | ABC-SMC calibration | ❌ Needs install |
| matplotlib/plotly | Visualization | ✅ Installed |
| NumPy/SciPy | Core computation | ✅ Installed |

---

## 9. Risk Register

| Risk | Probability | Impact | Mitigation |
|------|------------|--------|-----------|
| Sobol run crashes/hangs | Low | 1–2 day delay | Checkpointing every 500 runs |
| ABC doesn't converge | Medium | 1 week delay | Start with fewer calibration targets, widen priors |
| Scale correction fails (strong N-dependence) | Low | Conceptual problem | Calibrate directly at large K with more compute |
| Power outage during long run | Medium | Hours lost | Checkpoint + auto-restart |
| Posterior is multimodal | Medium | Interpretation difficulty | Report all modes, run scenarios from each |
| Key parameter unidentifiable | High (for rho_rec) | Scientific limitation | Report explicitly, test sensitivity of conclusions to rho_rec value |
