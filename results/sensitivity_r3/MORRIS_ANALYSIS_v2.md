# SA Round 3 — Morris Screening Analysis

**Date:** 2026-02-19  
**Runs:** 850/880 successful (97%), 30 errors (extreme parameter combinations → crashes)  
**Config:** 3-node spatial (Sitka / Howe Sound / Monterey), K = 5,000 each, 20-year simulation, pathogen evolution enabled, 1 movement substep/day

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

Mean normalized μ* across all 20 metrics. Higher = more influential.

| Rank | Parameter | μ\*_norm | Category | vs. Sobol R1 |
|------|-----------|---------|----------|--------------|
| 1 | **rho_rec** | 0.666 | Disease | ↑ was #14 |
| 2 | **peak_width_days** | 0.544 | Spawning | 🆕 |
| 3 | **settler_survival** | 0.490 | Population | ↑ was #19 |
| 4 | **a_exposure** | 0.488 | Disease | → was #3 |
| 5 | **target_mean_r** | 0.449 | Genetics | 🆕 |
| 6 | k_growth | 0.378 | Population | ↑ was #13 |
| 7 | sigma_D | 0.315 | Disease | ↑ |
| 8 | L_min_repro | 0.309 | Population | → was #10 |
| 9 | mu_I2D_ref | 0.308 | Disease | ↓ was #1 |
| 10 | P_env_max | 0.299 | Disease | ↑ was #23 |
| 11 | F0 | 0.299 | Population | ↑ |
| 12 | K_half | 0.290 | Disease | ↑ |
| 13 | sigma_2_eff | 0.280 | Disease | ↓ was #4 |
| 14 | mu_EI1_ref | 0.276 | Disease | → was #6 |
| 15 | p_spontaneous_male | 0.269 | Spawning | 🆕 |
| 16 | q_init_beta_a | 0.264 | Genetics | 🆕 |
| 17 | induction_female_to_male | 0.262 | Spawning | 🆕 |
| 18 | susceptibility_multiplier | 0.250 | Disease | ↓ |
| 19 | T_vbnc | 0.250 | Disease | — |
| 20 | alpha_kill | 0.249 | Path. Evo | 🆕 |
| 21 | n_additive | 0.243 | Genetics | ↓ was #5 |
| 22 | min_susceptible_age_days | 0.242 | Disease | 🆕 |
| 23 | T_ref | 0.235 | Disease | — |
| 24 | senescence_age | 0.233 | Population | — |
| 25 | readiness_induction_prob | 0.232 | Spawning | 🆕 |
| 26 | sigma_v_mutation | 0.231 | Path. Evo | 🆕 |
| 27 | immunosuppression_duration | 0.229 | Disease | 🆕 |
| 28 | induction_male_to_female | 0.229 | Spawning | 🆕 |
| 29 | female_max_bouts | 0.221 | Spawning | 🆕 |
| 30 | gamma_early | 0.213 | Path. Evo | 🆕 |
| 31 | mu_I1I2_ref | 0.213 | Disease | 🆕 |
| 32 | s_min | 0.208 | Disease | — |
| 33 | alpha_prog | 0.207 | Path. Evo | 🆕 |
| 34 | D_L | 0.202 | Spatial | — |
| 35 | p_spontaneous_female | 0.199 | Spawning | — |
| 36 | alpha_self_open | 0.197 | Spatial | 🆕 |
| 37 | v_init | 0.189 | Path. Evo | 🆕 |
| 38 | q_init_beta_b | 0.187 | Genetics | 🆕 |
| 39 | alpha_srs | 0.185 | Population | — |
| 40 | alpha_self_fjord | 0.183 | Spatial | 🆕 |
| 41 | sigma_1_eff | 0.182 | Disease | — |
| 42 | alpha_shed | 0.176 | Path. Evo | 🆕 |
| 43 | gamma_fert | 0.174 | Population | — (last) |

**Elimination threshold (5% of max):** 0.033. **All 43 parameters above threshold — zero eliminated.**

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
| 1 | k_growth | 1.000 |
| 2 | senescence_age | 0.133 |
| 3 | peak_width_days | 0.071 |

**Interpretation:** Growth rate dominates because it determines when individuals reach reproductive size (L_min_repro). Pre-disease recruitment is almost entirely a function of demographics, not disease or genetics — as expected, since disease is introduced at year 3.

---

## 4. Key Findings

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

## 5. Round 1 vs Round 3 Comparison

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

## 6. Sobol Plan

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
