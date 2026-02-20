# SA Round 3 — Morris Screening Analysis

**Date:** 2026-02-19 (v2.1 — corrected re-run, degenerate metrics fixed)  
**Runs:** 880/880 successful (100%), 0 errors  
**Config:** 3-node spatial (Sitka / Howe Sound / Monterey), K = 5,000 each, 20-year simulation, pathogen evolution enabled, 1 movement substep/day  
**Wall time:** 110.9 min on 8 cores

---

## 1. Variable Definitions

### 1.1 Parameters (43 total)

All 43 parameters tested, organized by module. Each is varied independently across its sampling range.

#### Disease Module (15 parameters)

| Symbol in Model | Name | Description | Range | Units | Distribution | Empirical Confidence |
|----------------|------|-------------|-------|-------|-------------|---------------------|
| `a_exposure` | Exposure rate | Rate at which susceptible individuals transition to exposed (S→E) per unit pathogen concentration. Controls epidemic speed. | [0.30, 1.50] | d⁻¹ | Uniform | ★☆☆ Lab only |
| `K_half` | Half-infective dose | Pathogen concentration at which probability of infection = 50% (Michaelis-Menten saturation). | [2×10⁴, 2×10⁵] | bacteria mL⁻¹ | Log-uniform | ★☆☆ |
| `sigma_1_eff` | I₁ shedding rate | Field-effective pathogen shedding rate from early-infected (I₁) hosts. | [1.0, 25.0] | bacteria mL⁻¹ d⁻¹ | Log-uniform | ★☆☆ |
| `sigma_2_eff` | I₂ shedding rate | Field-effective pathogen shedding rate from late-infected (I₂) hosts. Late stage sheds 10–250× more. | [10.0, 250.0] | bacteria mL⁻¹ d⁻¹ | Log-uniform | ★☆☆ |
| `sigma_D` | Saprophytic burst | Pathogen release from decaying carcasses. Creates environmental reservoir after host death. | [3.0, 75.0] | bacteria mL⁻¹ d⁻¹ | Log-uniform | ★☆☆ |
| `rho_rec` | Recovery rate | Daily probability of clearing infection and transitioning to recovered (R). **No empirical basis.** | [0.0, 0.20] | d⁻¹ | Uniform | ★☆☆ None |
| `mu_EI1_ref` | E→I₁ progression | Rate of latent-to-early-infected transition at reference temperature T_ref. | [0.20, 1.00] | d⁻¹ | Uniform | ★★☆ |
| `mu_I1I2_ref` | I₁→I₂ progression | Rate of early-to-late infection transition at T_ref. | [0.15, 0.80] | d⁻¹ | Uniform | ★★☆ |
| `mu_I2D_ref` | I₂→Death rate | Disease-induced mortality rate from late infection at T_ref. | [0.08, 0.35] | d⁻¹ | Uniform | ★★☆ |
| `P_env_max` | Background Vibrio | Environmental V. pectenicida input rate (independent of shedding). Represents waterborne reservoir. | [50, 5,000] | bacteria mL⁻¹ d⁻¹ | Log-uniform | ★☆☆ |
| `T_ref` | Pathogen T_opt | Optimal temperature for V. pectenicida activity (Arrhenius reference). | [17.0, 23.0] | °C | Uniform | ★★☆ |
| `susceptibility_multiplier` | Immunosuppression factor | Multiplicative increase in disease susceptibility during post-spawning immunosuppression window. | [1.0, 4.0] | dimensionless | Uniform | ★☆☆ |
| `T_vbnc` | VBNC midpoint | Temperature below which V. pectenicida enters viable-but-non-culturable (VBNC) state. | [8.0, 15.0] | °C | Uniform | ★★☆ |
| `s_min` | Salinity minimum | Minimum salinity for Vibrio survival. Fjords with lower salinity get partial protection. | [5.0, 15.0] | PSU | Uniform | ★★☆ |
| `min_susceptible_age_days` | Juvenile immunity | Days post-settlement before juveniles become susceptible to infection. | [0, 180] | days | Uniform | ★☆☆ |
| `immunosuppression_duration` | Immunosuppression window | Duration of elevated disease susceptibility after spawning. | [7, 56] | days | Uniform | ★★☆ |

#### Population Module (7 parameters)

| Symbol in Model | Name | Description | Range | Units | Distribution | Empirical Confidence |
|----------------|------|-------------|-------|-------|-------------|---------------------|
| `F0` | Reference fecundity | Egg production per spawning bout for a reference-size female. | [10⁶, 10⁸] | eggs | Log-uniform | ★☆☆ |
| `gamma_fert` | Fertilization kinetics | Shape parameter controlling sperm-limitation of fertilization at low male density. | [1.0, 10.0] | dimensionless | Uniform | ★☆☆ |
| `settler_survival` | Settler survival (s₀) | Beverton-Holt baseline fraction of settling larvae that survive to recruitment. | [0.005, 0.10] | fraction | Log-uniform | ★☆☆ |
| `alpha_srs` | SRS Pareto shape | Sweepstakes Reproductive Success shape parameter (α). Lower α = more reproductive variance. Hedgecock (2023): α ≈ 1.35. | [1.0, 1.8] | dimensionless | Uniform | ★★☆ |
| `senescence_age` | Senescence onset | Age at which daily natural mortality begins increasing exponentially (Gompertz). | [20.0, 80.0] | years | Uniform | ★☆☆ |
| `k_growth` | VB growth rate | Von Bertalanffy growth rate constant. Faster growth → earlier reproductive maturity. | [0.03, 0.15] | yr⁻¹ | Uniform | ★☆☆ |
| `L_min_repro` | Min reproductive size | Minimum arm disc radius for spawning eligibility. | [200, 500] | mm | Uniform | ★☆☆ |

#### Genetics Module (4 parameters)

| Symbol in Model | Name | Description | Range | Units | Distribution | Empirical Confidence |
|----------------|------|-------------|-------|-------|-------------|---------------------|
| `n_additive` | Additive loci | Number of additive resistance loci (n_loci = n_additive + 1 for EF1A). | {10, 20, 30, 40, 51} | count | Discrete | ★★☆ |
| `target_mean_r` | Initial mean resistance | Population-mean innate resistance (r_i) at simulation start. Pre-outbreak standing variation. | [0.05, 0.30] | dimensionless | Uniform | ★☆☆ |
| `q_init_beta_a` | Beta shape a | Shape parameter 'a' of Beta(a,b) distribution for per-locus initial allele frequencies. | [1.0, 5.0] | dimensionless | Uniform | ★☆☆ |
| `q_init_beta_b` | Beta shape b | Shape parameter 'b' of Beta(a,b) distribution for per-locus initial allele frequencies. | [3.0, 15.0] | dimensionless | Uniform | ★☆☆ |

#### Spawning Module (6 parameters)

| Symbol in Model | Name | Description | Range | Units | Distribution | Empirical Confidence |
|----------------|------|-------------|-------|-------|-------------|---------------------|
| `p_spontaneous_female` | Female spontaneous spawn | Daily probability a ripe female initiates spawning spontaneously (without induction). | [0.005, 0.025] | d⁻¹ | Uniform | ★★☆ |
| `p_spontaneous_male` | Male spontaneous spawn | Daily probability a ripe male spawns spontaneously. | [0.005, 0.025] | d⁻¹ | Uniform | ★★☆ |
| `induction_female_to_male` | F→M cascade induction (κ_fm) | Probability that a spawning female triggers spawning in a nearby ripe male. Asymmetric: females strongly induce males. | [0.40, 0.95] | probability | Uniform | ★★☆ |
| `induction_male_to_female` | M→F cascade induction (κ_mf) | Probability that a spawning male induces a nearby ripe female. Weaker than F→M. | [0.10, 0.60] | probability | Uniform | ★★☆ |
| `peak_width_days` | Spawning peak width (σ_spawn) | Standard deviation of the Gaussian seasonal spawning window. Narrow = synchronized pulse; wide = diffuse season. | [30, 90] | days | Uniform | ★★☆ |
| `readiness_induction_prob` | Social readiness induction | Probability of socially-induced spawning readiness (nearby spawner accelerates gonad maturation). | [0.10, 0.80] | probability | Uniform | ★☆☆ |
| `female_max_bouts` | Female max bouts | Maximum spawning bouts per female per season. Males can spawn 2–3× per season regardless. | {1, 2, 3} | count | Discrete | ★★☆ |

#### Pathogen Evolution Module (6 parameters)

| Symbol in Model | Name | Description | Range | Units | Distribution | Empirical Confidence |
|----------------|------|-------------|-------|-------|-------------|---------------------|
| `alpha_kill` | Kill rate exponent | How virulence phenotype (v) scales death rate: μ_I2D ∝ v^α_kill. Higher = steeper virulence-mortality relationship. | [1.0, 4.0] | dimensionless | Uniform | ★☆☆ |
| `alpha_shed` | Shedding exponent | How virulence scales pathogen shedding: σ ∝ v^α_shed. The virulence-transmission tradeoff mechanism. | [0.5, 3.0] | dimensionless | Uniform | ★☆☆ |
| `alpha_prog` | Progression exponent | How virulence scales I₁→I₂ progression rate. | [0.5, 2.0] | dimensionless | Uniform | ★☆☆ |
| `gamma_early` | Early shedding attenuation | Fraction of full shedding during I₁ stage (early infection). 0 = no early shedding, 1 = full. | [0.0, 1.0] | fraction | Uniform | ★☆☆ |
| `sigma_v_mutation` | Mutation step size | Standard deviation of virulence mutation during transmission (normal perturbation). | [0.005, 0.10] | dimensionless | Log-uniform | ★☆☆ |
| `v_init` | Initial virulence | Pathogen virulence phenotype at simulation start (v ∈ [0,1] scale). | [0.2, 0.8] | dimensionless | Uniform | ★☆☆ |

#### Spatial Module (3 parameters)

| Symbol in Model | Name | Description | Range | Units | Distribution | Empirical Confidence |
|----------------|------|-------------|-------|-------|-------------|---------------------|
| `D_L` | Larval dispersal scale | Characteristic dispersal distance for larval transport kernel (exponential decay). | [100, 1,000] | km | Log-uniform | ★☆☆ |
| `alpha_self_fjord` | Fjord self-recruitment | Fraction of larvae retained locally in fjord habitats (higher due to sill-restricted circulation). | [0.10, 0.50] | fraction | Uniform | ★☆☆ |
| `alpha_self_open` | Open coast self-recruitment | Fraction of larvae retained locally at open coast sites. | [0.02, 0.20] | fraction | Uniform | ★☆☆ |

### 1.2 Metrics (20 total)

| Metric | Description | Range | Interpretation |
|--------|-------------|-------|---------------|
| `pop_crash_pct` | Maximum decline from initial population (%) | [0, 100] | Higher = more severe epidemic |
| `final_pop_frac` | Final population / initial population | [0, ∞) | < 1 = net decline; > 1 = recovery overshoot |
| `recovery` | Binary: did population recover to > 50% of initial? | {0, 1} | 1 = demographic recovery |
| `extinction` | Binary: did total population reach 0? | {0, 1} | 1 = metapopulation extinction |
| `peak_mortality` | Highest single-year mortality fraction at any node | [0, 1] | Intensity of epidemic peak |
| `time_to_nadir` | Year of minimum population | [0, 20] | Speed of population crash |
| `total_disease_deaths` | Cumulative disease-caused deaths across all nodes | [0, ∞) | Total epidemic toll |
| `resistance_shift_mean` | Mean Δr̄ (change in mean resistance) across nodes, post-disease vs pre-disease | (-∞, ∞) | Positive = adaptive evolution |
| `resistance_shift_max` | Maximum Δr̄ at any single node | (-∞, ∞) | Strongest local selection signal |
| `va_retention_mean` | Mean V_A(post) / V_A(pre) across nodes | [0, ∞) | < 1 = genetic variation lost; selection exhausting standing variation |
| `ef1a_shift_mean` | Mean absolute change in EF1A (overdominant) locus frequency | [0, 1] | Overdominant selection signal |
| `n_extinct_nodes` | Number of nodes with zero final population | [0, 3] | Spatial extent of extinction |
| `north_south_mortality_gradient` | (Monterey cumulative disease deaths / K) − (Sitka / K) | (-∞, ∞) | Positive = south hit harder (expected from temperature gradient) |
| `fjord_protection_effect` | (Howe Sound final / K) − mean(non-fjord final / K) | (-∞, ∞) | Positive = fjord confers protection |
| `mean_final_virulence` | Mean pathogen virulence (v) across nodes at year 20 | [0, 1] | Evolved virulence level |
| `virulence_shift` | Final mean virulence − initial virulence | (-1, 1) | Direction of pathogen evolution |
| `disease_death_fraction` | Disease deaths / (disease + natural + senescence deaths) | [0, 1] | Relative importance of disease as mortality source |
| `spawning_participation` | Fraction of pre-disease years with nonzero recruitment | [0, 1] | Spawning system functionality |
| `mean_recruitment_rate` | Mean annual recruits / total K during pre-disease years | [0, ∞) | Demographic input rate |
| `evolutionary_rescue_index` | final_pop_frac × resistance_shift_mean | (-∞, ∞) | Composite: populations that both survived AND evolved resistance |

### 1.3 Methods

**Morris screening** (Elementary Effects Method): One-at-a-time perturbation of each parameter across r = 20 trajectories through the 43-dimensional space. For each parameter and metric, we compute:
- **μ\*** (mean absolute elementary effect): overall importance — how much does this parameter move this metric on average?
- **σ** (standard deviation of elementary effects): interaction/nonlinearity — does the parameter's effect depend on what other parameters are doing?
- When **σ > μ\***: the parameter's effect is dominated by interactions with other parameters, not its direct (additive) effect.

**Global ranking**: For each metric, μ\* values are normalized by the maximum μ\* for that metric (so the most important parameter for each metric = 1.0). The global rank is the mean of these normalized scores across all 20 metrics.

**Elimination threshold**: Parameters with global μ\*_norm < 5% of the top-ranked parameter are candidates for elimination from Sobol analysis.

---

## 2. Global Parameter Ranking

Mean normalized μ* across all 20 metrics. Higher = more influential. (v2.1 — corrected with non-degenerate disease_death_fraction and spawning_participation)

| Rank | Parameter | μ\*_norm | Mean Rank | Top5 in | Category | vs. Sobol R1 |
|------|-----------|---------|-----------|---------|----------|--------------|
| 1 | **rho_rec** | 0.642 | 4.7 | 15/20 | Disease | ↑ was #14 |
| 2 | **peak_width_days** | 0.573 | 9.3 | 13/20 | Spawning | 🆕 |
| 3 | **settler_survival** | 0.571 | 7.5 | 13/20 | Population | ↑ was #19 |
| 4 | **target_mean_r** | 0.432 | 11.2 | 8/20 | Genetics | 🆕 |
| 5 | **k_growth** | 0.402 | 11.2 | 5/20 | Population | ↑ was #13 |
| 6 | **a_exposure** | 0.394 | 8.0 | 10/20 | Disease | → was #3 |
| 7 | mu_I2D_ref | 0.322 | 13.1 | 4/20 | Disease | ↓ was #1 |
| 8 | K_half | 0.303 | 15.6 | 3/20 | Disease | ↑ |
| 9 | T_vbnc | 0.298 | 16.2 | 3/20 | Disease | — |
| 10 | sigma_2_eff | 0.298 | 13.5 | 6/20 | Disease | ↓ was #4 |
| 11 | P_env_max | 0.287 | 17.6 | 3/20 | Disease | ↑ was #23 |
| 12 | v_init | 0.263 | 17.4 | 0/20 | Path. Evo | 🆕 |
| 13 | min_susceptible_age_days | 0.249 | 20.3 | 1/20 | Disease | 🆕 |
| 14 | mu_I1I2_ref | 0.241 | 23.1 | 2/20 | Disease | 🆕 |
| 15 | immunosuppression_duration | 0.231 | 24.9 | 2/20 | Disease | 🆕 |
| 16 | induction_male_to_female | 0.227 | 20.2 | 0/20 | Spawning | 🆕 |
| 17 | q_init_beta_b | 0.225 | 26.1 | 1/20 | Genetics | 🆕 |
| 18 | D_L | 0.222 | 23.0 | 0/20 | Spatial | — |
| 19 | n_additive | 0.221 | 23.2 | 1/20 | Genetics | ↓ was #5 |
| 20 | F0 | 0.220 | 23.1 | 2/20 | Population | ↑ |
| 21 | senescence_age | 0.220 | 21.4 | 1/20 | Population | — |
| 22 | alpha_kill | 0.219 | 25.1 | 0/20 | Path. Evo | 🆕 |
| 23 | susceptibility_multiplier | 0.215 | 22.2 | 1/20 | Disease | ↓ |
| 24 | induction_female_to_male | 0.214 | 23.4 | 2/20 | Spawning | 🆕 |
| 25 | L_min_repro | 0.214 | 24.2 | 2/20 | Population | → was #10 |
| 26 | p_spontaneous_female | 0.212 | 19.8 | 0/20 | Spawning | — |
| 27 | mu_EI1_ref | 0.209 | 20.1 | 0/20 | Disease | → was #6 |
| 28 | p_spontaneous_male | 0.198 | 27.9 | 0/20 | Spawning | 🆕 |
| 29 | sigma_D | 0.196 | 24.1 | 0/20 | Disease | ↑ |
| 30 | gamma_early | 0.193 | 28.8 | 0/20 | Path. Evo | 🆕 |
| 31 | sigma_v_mutation | 0.192 | 27.6 | 0/20 | Path. Evo | 🆕 |
| 32 | female_max_bouts | 0.187 | 27.6 | 0/20 | Spawning | 🆕 |
| 33 | readiness_induction_prob | 0.180 | 30.1 | 0/20 | Spawning | 🆕 |
| 34 | T_ref | 0.176 | 26.1 | 0/20 | Disease | — |
| 35 | alpha_srs | 0.169 | 29.2 | 0/20 | Population | — |
| 36 | s_min | 0.168 | 30.5 | 0/20 | Disease | — |
| 37 | gamma_fert | 0.167 | 27.4 | 0/20 | Population | — |
| 38 | alpha_prog | 0.164 | 32.4 | 1/20 | Path. Evo | 🆕 |
| 39 | alpha_self_open | 0.162 | 30.9 | 0/20 | Spatial | 🆕 |
| 40 | q_init_beta_a | 0.161 | 28.1 | 0/20 | Genetics | 🆕 |
| 41 | alpha_shed | 0.160 | 30.9 | 0/20 | Path. Evo | 🆕 |
| 42 | alpha_self_fjord | 0.149 | 31.6 | 0/20 | Spatial | 🆕 |
| 43 | sigma_1_eff | 0.146 | 27.3 | 1/20 | Disease | — (last) |

**Elimination threshold (5% of max):** 0.032. **All 43 parameters above threshold (minimum: sigma_1_eff at 0.146) — zero eliminated.**

---

## 3. Category-Specific Top Drivers

### 3.1 Population Outcomes (crash, recovery, extinction, mortality, time-to-nadir, total deaths)
| Rank | Parameter | Normalized μ\* |
|------|-----------|---------------|
| 1 | rho_rec | 0.688 |
| 2 | settler_survival | 0.611 |
| 3 | peak_width_days | 0.550 |
| 4 | a_exposure | 0.501 |
| 5 | target_mean_r | 0.412 |

**Interpretation:** Recovery rate dominates population fate — even more than disease exposure. This makes biological sense: in a chronic epidemic, what matters most is not how fast you get sick but whether you can clear the infection. settler_survival (juvenile demographic input) is the second lever — it controls whether the population can replace disease losses faster than they accumulate.

### 3.2 Evolutionary Outcomes (resistance shift, V_A retention, EF1A dynamics, rescue index)
| Rank | Parameter | Normalized μ\* |
|------|-----------|---------------|
| 1 | rho_rec | 0.691 |
| 2 | target_mean_r | 0.662 |
| 3 | a_exposure | 0.528 |
| 4 | settler_survival | 0.453 |
| 5 | sigma_D | 0.422 |

**Interpretation:** rho_rec is #1 because recovery creates the selection window — only individuals who survive long enough to recover can pass on resistance alleles. target_mean_r matters because it sets how much standing variation selection has to work with: at r̄ = 0.05, almost nobody has enough resistance for recovery to be selectively advantageous; at r̄ = 0.25, there's a meaningful high-resistance tail.

### 3.3 Spatial Patterns (N–S gradient, fjord protection, node extinction)
| Rank | Parameter | Normalized μ\* |
|------|-----------|---------------|
| 1 | settler_survival | 0.704 |
| 2 | rho_rec | 0.672 |
| 3 | peak_width_days | 0.615 |
| 4 | k_growth | 0.352 |
| 5 | a_exposure | 0.315 |

**Interpretation:** settler_survival determines which nodes can sustain themselves demographically — nodes below the recruitment threshold go extinct regardless of disease. peak_width_days interacts with latitude because spawning seasonality is temperature-modulated: narrow peak × warm site = different epidemic-recruitment phasing than narrow peak × cold site.

### 3.4 Pathogen Evolution (final virulence, virulence shift)
| Rank | Parameter | Normalized μ\* |
|------|-----------|---------------|
| 1 | peak_width_days | 1.000 |
| 2 | rho_rec | 0.827 |
| 3 | a_exposure | 0.827 |
| 4 | L_min_repro | 0.758 |
| 5 | sigma_D | 0.691 |

**Interpretation:** peak_width_days is the absolute #1 driver of virulence evolution — this was the biggest surprise. Mechanism: narrow spawning peak → synchronized larval settlement → large age-cohort of naive susceptibles entering the population simultaneously → intense epidemic wave → strong selection pressure on pathogen virulence. Wide peak diffuses settlement → disease encounters spread over time → weaker selection on virulence. Note that alpha_kill (#20 globally) doesn't even appear in this top 5 — the shape of the host population's temporal dynamics matters more than the pathogen's own genetic architecture.

### 3.5 Spawning/Recruitment (participation rate, recruitment rate, disease death fraction)
| Rank | Parameter | Normalized μ\* |
|------|-----------|---------------|
| 1 | settler_survival | 0.714 |
| 2 | k_growth | 0.703 |
| 3 | peak_width_days | 0.632 |
| 4 | rho_rec | 0.410 |
| 5 | target_mean_r | 0.372 |

**Interpretation:** [v2.1 CORRECTED — previously showed k_growth dominating at 1.000 with a 7.5× gap to #2, which was an artifact of disease_death_fraction and spawning_participation being degenerate (zero variance). With the corrected re-run, these metrics now have healthy variance and the ranking is more balanced.] settler_survival and k_growth jointly dominate — one controls demographic input from the sea (larval settlement success), the other controls how fast individuals reach reproductive size. peak_width_days matters because spawning synchrony determines whether larval settlement is a single pulse or diffuse season, with cascading effects on recruitment timing relative to epidemic waves. Notably, disease parameters (rho_rec, target_mean_r) now appear in this category's top 5 — disease-mediated adult mortality feeds back into spawning participation (fewer adults = fewer spawners = lower participation).

---

## 4. Absolute Effect Sizes

The normalized rankings (Section 2) show *relative* importance — which parameters matter more than others. This section reports **raw μ\*** values in each metric's native units, showing how much each parameter actually moves the outcome when perturbed across its range.

**How to read μ\*:** The mean absolute elementary effect. If μ\* = 9.5 for pop_crash_pct, that means perturbing this parameter (one-at-a-time, across Morris trajectories) changes population crash by ±9.5 percentage points on average.

### 4.1 Population Crash (%)

Baseline crash in the default parameterization is ~95–99%. μ\* is in percentage-point units.

| Rank | Parameter | μ\* (pp) | σ | 95% CI | σ/μ\* |
|------|-----------|---------|---|--------|-------|
| 1 | a_exposure | 10.27 | 22.14 | ±8.93 | 2.16 |
| 2 | P_env_max | 9.52 | 19.32 | ±8.51 | 2.03 |
| 3 | rho_rec | 9.50 | 8.99 | ±4.22 | 0.95 |
| 4 | target_mean_r | 9.41 | 8.13 | ±3.29 | 0.86 |
| 5 | sigma_2_eff | 9.30 | — | — | — |
| 6 | K_half | 8.33 | — | — | — |
| 7 | settler_survival | 7.66 | 10.82 | ±3.48 | 1.41 |
| 8 | mu_I2D_ref | 6.29 | 7.92 | ±3.07 | 1.26 |
| 9 | p_spontaneous_male | 5.40 | — | — | — |
| 10 | T_vbnc | 4.99 | — | — | — |

**Interpretation:** The top parameters each shift crash magnitude by 5–10 percentage points — substantial given that baseline crash is ~96%. Note a_exposure has μ\* = 10.3 but σ/μ\* = 2.16 — its effect is *highly* interaction-dependent (the crash impact of exposure rate depends strongly on what recovery rate, initial resistance, etc. are doing). In contrast, rho_rec and target_mean_r have σ/μ\* ≈ 0.9 — their effects are more consistent (closer to additive).

### 4.2 Total Disease Deaths (individuals)

Total across all 3 nodes × 20 years. Baseline total K = 15,000.

| Rank | Parameter | μ\* (deaths) | σ | σ/μ\* |
|------|-----------|-------------|---|-------|
| 1 | settler_survival | 54,409 | — | — |
| 2 | peak_width_days | 30,943 | — | — |
| 3 | k_growth | 23,595 | — | — |
| 4 | mu_I2D_ref | 15,062 | — | — |
| 5 | rho_rec | 11,800 | — | — |
| 6 | K_half | 11,387 | — | — |
| 7 | a_exposure | 10,756 | — | — |

**Interpretation:** settler_survival dominates total deaths (μ\* = 54,409 — over 3× total K!) because higher settler survival → larger standing population → more individuals available to die from disease. This is a demographic amplification effect: more recruitment sustains a larger host population that disease then kills. peak_width_days at 30,943 deaths reflects the same amplification through synchronized settlement creating vulnerable age cohorts.

### 4.3 Mean Resistance Shift (Δr̄)

Change in population-mean innate resistance from pre-epidemic (year 3) to final year. Schiebelhut observed Δq = 0.08–0.15 at individual loci in wild populations.

| Rank | Parameter | μ\* (Δr̄) | σ | σ/μ\* |
|------|-----------|----------|---|-------|
| 1 | rho_rec | 0.073 | 0.125 | 1.72 |
| 2 | target_mean_r | 0.048 | 0.069 | 1.45 |
| 3 | a_exposure | 0.034 | 0.051 | 1.49 |
| 4 | P_env_max | 0.028 | 0.058 | 2.07 |
| 5 | K_half | 0.027 | — | — |
| 6 | L_min_repro | 0.024 | 0.047 | 2.01 |
| 7 | peak_width_days | 0.023 | 0.035 | 1.49 |

**Interpretation:** rho_rec dominates evolutionary response with μ\* = 0.073 — meaning perturbing recovery rate shifts mean resistance by ±0.07 on average. This is in the same order of magnitude as Schiebelhut's observed Δq, suggesting recovery rate alone can determine whether the model produces realistic selection signals. But σ/μ\* = 1.72 means this effect is strongly interaction-dependent. target_mean_r at 0.048 confirms that standing variation is the raw material for selection — you can't evolve resistance you don't have.

### 4.4 Final Virulence (v, dimensionless 0–1 scale)

Pathogen virulence at year 20, starting from v_init ∈ [0.2, 0.8].

| Rank | Parameter | μ\* (Δv) | σ | σ/μ\* |
|------|-----------|----------|---|-------|
| 1 | peak_width_days | 0.167 | 0.345 | 2.07 |
| 2 | rho_rec | 0.138 | 0.290 | 2.10 |
| 3 | a_exposure | 0.138 | 0.300 | 2.17 |
| 4 | L_min_repro | 0.127 | 0.294 | 2.32 |
| 5 | sigma_D | 0.116 | 0.264 | 2.28 |

**Interpretation:** All top parameters show σ/μ\* > 2.0 for virulence — virulence evolution is **entirely interaction-driven**. No single parameter has a consistent, additive effect on evolved virulence. This makes biological sense: virulence evolution depends on the joint epidemiological context (host density × transmission × recovery × mortality), not any single axis. peak_width_days at μ\* = 0.167 means spawning timing can shift evolved virulence by ±0.17 on a 0–1 scale — a massive effect.

### 4.5 Fjord Protection Effect (Δ survival fraction)

Howe Sound final survival minus mean non-fjord survival. Positive = fjord is protective.

| Rank | Parameter | μ\* | σ | σ/μ\* |
|------|-----------|-----|---|-------|
| 1 | settler_survival | 0.243 | 0.305 | 1.25 |
| 2 | peak_width_days | 0.199 | 0.277 | 1.39 |
| 3 | rho_rec | 0.160 | 0.241 | 1.51 |
| 4 | mu_I2D_ref | 0.123 | 0.241 | 1.96 |
| 5 | F0 | 0.133 | — | — |

**Interpretation:** settler_survival has the largest effect on fjord protection (μ\* = 0.24 — nearly a quarter of the entire survival scale). The fjord advantage depends on whether recruitment can keep up with losses: high settler survival amplifies the self-recruitment advantage fjords have due to sill-restricted larval retention.

### 4.6 Peak Mortality Rate (fraction of node population in worst year)

| Rank | Parameter | μ\* | σ | σ/μ\* |
|------|-----------|-----|---|-------|
| 1 | peak_width_days | 33.81 | — | — |
| 2 | settler_survival | 20.81 | — | — |
| 3 | T_vbnc | 11.43 | — | — |
| 4 | a_exposure | 11.16 | — | — |
| 5 | mu_I2D_ref | 10.59 | — | — |

**Interpretation:** peak_width_days has by far the largest absolute effect (μ\* = 33.8, in percentage-point units). Narrow spawning peak creates a synchronized pulse of naive juveniles → epidemic spike. Note: values >1.0 (>100%) are possible because this metric measures disease deaths / start-of-year population, and within-year recruitment can inflate the denominator.

### 4.7 Previously-Degenerate Metrics — Now Resolved

In the initial Morris run, disease_death_fraction and spawning_participation showed zero variance (always ≈ 1.0). The corrected re-run (880/880 runs, 0 errors) resolved this:

#### disease_death_fraction (std = 0.137)

Fraction of total deaths caused by disease (vs natural + senescence). No longer saturated at 1.0.

| Rank | Parameter | μ\* | σ | σ/μ\* |
|------|-----------|-----|---|-------|
| 1 | k_growth | 0.238 | 0.118 | 0.50 |
| 2 | rho_rec | 0.072 | 0.079 | 1.10 |
| 3 | settler_survival | 0.069 | 0.059 | 0.86 |
| 4 | min_susceptible_age_days | 0.065 | 0.106 | 1.63 |
| 5 | a_exposure | 0.054 | 0.072 | 1.32 |

**Interpretation:** k_growth dominates with σ/μ\* = 0.50 (nearly additive effect!) — faster growth means individuals reach larger sizes where natural and senescence mortality apply, diluting disease's share. The low σ/μ\* is notable: this is one of the few near-additive effects in the entire analysis. min_susceptible_age_days at #4 makes mechanistic sense: longer juvenile immunity delays disease exposure, reducing disease deaths relative to other causes.

#### spawning_participation (std = 0.196)

Fraction of pre-disease years with nonzero recruitment. Now shows meaningful variation.

| Rank | Parameter | μ\* | σ | σ/μ\* |
|------|-----------|-----|---|-------|
| 1 | settler_survival | 0.323 | 0.172 | 0.53 |
| 2 | peak_width_days | 0.265 | 0.376 | 1.42 |
| 3 | rho_rec | 0.093 | 0.119 | 1.28 |
| 4 | k_growth | 0.073 | 0.105 | 1.45 |
| 5 | target_mean_r | 0.071 | 0.097 | 1.35 |

**Interpretation:** settler_survival has the largest and most additive effect (σ/μ\* = 0.53) — below a threshold settler survival, recruitment drops to zero regardless of other conditions. peak_width_days at #2 with σ/μ\* = 1.42 reflects a nonlinear interaction: narrow peaks concentrate spawning into a window that either succeeds or fails entirely, while wide peaks spread the bet. The appearance of disease parameters (rho_rec, target_mean_r) confirms feedback: disease kills adults, reducing spawner counts, which lowers participation.

**Key finding:** Both metrics' σ/μ\* ratios are lower than the model average, meaning direct (additive) effects are stronger here than for most other metrics. This makes sense — these are demographic ratios that respond more mechanically to parameter changes than, say, evolutionary rescue or virulence evolution.

### 4.8 σ/μ\* Summary — Where Interactions Dominate

| Parameter | pop_crash | resist_shift | virulence | fjord_protect | Mean σ/μ\* |
|-----------|-----------|-------------|-----------|---------------|-----------|
| rho_rec | 0.9 | 1.7 | 2.1 | 1.5 | 1.6 |
| peak_width_days | 1.9 | 1.5 | 2.1 | 1.4 | 1.7 |
| settler_survival | 1.4 | 1.6 | 2.6 | 1.3 | 1.7 |
| a_exposure | **2.2** | 1.5 | **2.2** | 1.7 | 1.9 |
| target_mean_r | 0.9 | 1.5 | **2.6** | 1.4 | 1.6 |

**Key pattern:** σ/μ\* > 1 everywhere — interactions dominate for every top parameter on every metric. For virulence evolution, σ/μ\* > 2 universally. This confirms the Round 1 finding that **Sobol total-order indices (S_T) will be much larger than first-order (S₁)**. Morris μ\* alone cannot be trusted for ranking — it conflates direct and interaction effects. The Sobol decomposition is essential.

---

## 5. Key Findings

### 4.1 rho_rec: From Nobody to King (Sobol R1 #14 → Morris R3 #1)

In Round 1 (23 parameters, annual mortality steps, no pathogen evolution), rho_rec was Morris #1 but dropped to Sobol #14 — we dismissed it as a Morris artifact. Now it's Morris #1 again. What changed?

**Continuous daily mortality.** Under the old annual-timestep model, mortality was applied as a single annual purge. Recovery rate operating on a daily timescale was poorly resolved — a fast daily recovery gets "washed out" when the yearly accounting collapses everything into one event. With continuous daily mortality and daily recovery, each day's survival depends on the balance between recovery and death. The compounding effect of ρ_rec = 0.05 d⁻¹ over 20 days (≈ 64% cleared) vs ρ_rec = 0.01 (≈ 18% cleared) is now properly resolved.

**Key question for Sobol: Does rho_rec hold at #1 when we account for interactions, or does it drop again?** Morris overweights parameters with large marginal effects. If rho_rec's effect is mostly additive (σ ≈ μ\*), it will hold. If its effect is heavily interaction-dependent (σ >> μ\*), Sobol may downweight it as in Round 1.

### 4.2 peak_width_days: Spawning Timing Drives Virulence Evolution

This parameter was completely invisible before the spawning overhaul. Its emergence as #2 globally and **#1 for virulence evolution** demonstrates that the spawning system isn't just biological realism — it exposed a previously hidden sensitivity axis.

The mechanism chain: σ_spawn → settlement synchrony → age-cohort structure → epidemic pulse intensity → selection gradient on virulence. This is an indirect, multi-step pathway that only appears when all components (spawning model, continuous settlement, juvenile immunity, pathogen evolution) are coupled.

### 4.3 mu_I2D_ref Dethroned (Sobol R1 #1 → Morris R3 #9)

The former dominant parameter. Its decline likely reflects two changes:
1. **Continuous mortality** distributes deaths across days rather than concentrating them at year-end, reducing the marginal impact of the death rate.
2. **Expanded parameter space** (43 vs 23) means more parameters compete for variance, mechanically reducing any individual parameter's share.

### 4.4 Pathogen Evolution Parameters Are Secondary (#20–42)

alpha_kill is the most influential PE parameter at global rank #20. The pathogen's evolutionary machinery (mutation rate, scaling exponents, initial virulence) doesn't drive broad model behavior — host dynamics set the stage. But PE parameters **do** matter specifically for virulence metrics: they determine how fast and how far virulence evolves, given the host-imposed selection environment.

### 4.5 Zero Parameters Eliminated

All 43 parameters have μ\*_norm > 0.17 (minimum: gamma_fert at 0.174), well above the 5% elimination threshold of 0.033. Every parameter influences at least one metric meaningfully. Round 1 had the same result (0/23 eliminated).

---

## 6. Round 1 vs Round 3 Comparison

| Sobol R1 Rank | Parameter | Morris R3 Rank | Shift | Explanation |
|------|-----------|------|-------|-------------|
| 1 | mu_I2D_ref | 9 | ↓8 | Continuous mortality dilutes death-rate impact |
| 2 | susceptibility_multiplier | 18 | ↓16 | Expanded param space, spawning overhaul |
| 3 | a_exposure | 4 | ↑1 | Stable across rounds |
| 4 | sigma_2_eff | 13 | ↓9 | Continuous mortality effect |
| 5 | n_additive | 21 | ↓16 | Genetic architecture matters less with new init |
| 14 | rho_rec | 1 | ↑13 | Continuous mortality amplifies daily recovery |
| 19 | settler_survival | 3 | ↑16 | Demographics gain importance with daily resolution |
| — | peak_width_days | 2 | 🆕 | Spawning overhaul exposed hidden axis |
| — | target_mean_r | 5 | 🆕 | Genetic initialization matters |
| 23 | P_env_max | 10 | ↑13 | Environmental reservoir more important in continuous model |

**Summary:** The landscape shifted from disease-rate-dominated (Round 1: mu_I2D_ref, sigma_2_eff, susceptibility_multiplier in top 5) to recovery-and-demography-dominated (Round 3: rho_rec, peak_width_days, settler_survival). This reflects the model maturation: continuous daily dynamics, spawning overhaul, and pathogen evolution collectively rebalanced parameter importance.

---

## 7. Sobol Plan

- **Parameters:** All 43 (no elimination)
- **Sample size:** N = 512 (Saltelli sampling → 45,056 total runs)
- **Cores:** 8 (parallel pool)
- **Estimated time:** ~117 hours (4.9 days) at 77s/run mean
- **Checkpoint:** Every 500 runs (resumable after interruption)
- **Analysis:** S₁ (first-order), S_T (total-order), bootstrap CIs (1,000 resamples)
- **Completion marker:** `results/sensitivity_r3/sobol_complete.json`

### Priority Calibration Targets (post-Sobol)
1. **rho_rec** — #1, ★☆☆ confidence, zero empirical data
2. **peak_width_days** — #2, ★★☆ confidence, some phenological field data
3. **settler_survival** — #3, ★☆☆ confidence, typical marine invertebrate unknown
4. **a_exposure** — #4, ★☆☆ confidence, lab dose-response needed
5. **target_mean_r** — #5, ★☆☆ confidence, bounded by Schiebelhut Δq observations
