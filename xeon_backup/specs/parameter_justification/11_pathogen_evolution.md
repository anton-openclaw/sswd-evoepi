# 11. Pathogen Evolution Parameters — Literature Review

**Parameters reviewed:** alpha_kill, alpha_shed, alpha_prog, gamma_early, sigma_v_mutation, v_init

## Executive Summary

These six parameters define the evolutionary dynamics of Vibrio pectenicida virulence in our coupled eco-evolutionary model. The parameters encode the central theorem of evolutionary epidemiology: the virulence-transmission trade-off that governs pathogen evolution. While direct empirical data for these parameters is extremely sparse for SSWD, the theoretical framework is well-established through foundational work by Anderson & May (1982) and subsequent developments. Our parameterization is guided by comparative studies of other marine bacterial pathogens, terrestrial pathogen evolution experiments, and recent eco-evolutionary models like Clement et al. (2024) for Tasmanian devil facial tumor disease.

## First Principles Analysis

The six parameters implement a power-law trade-off framework where virulence (v) affects multiple fitness components:

- **Mortality rate** scales as: v^{alpha_kill}
- **Transmission rate** (shedding) scales as: v^{alpha_shed}  
- **Disease progression** (I₁→I₂) scales as: v^{alpha_prog}

The **evolutionary stable strategy (ESS) virulence** depends on the ratio alpha_kill/alpha_shed. When alpha_kill > alpha_shed, mortality costs exceed transmission benefits at high virulence, favoring intermediate virulence. When alpha_kill = alpha_shed (linear trade-off), selection favors either maximum or minimum virulence (boundary solution).

**gamma_early** controls whether I₁ (asymptomatic) individuals contribute to transmission. gamma_early = 0 means no early shedding; gamma_early = 1 means full shedding during I₁. This parameter determines whether "stealth" transmission dominates.

**sigma_v_mutation** is the phenotypic mutation step size (NOT the per-base DNA mutation rate). Each transmission event mutates virulence by N(0, sigma_v_mutation). This controls the speed of evolutionary adaptation.

**v_init** is the ancestral virulence at outbreak initiation (2013). If SSWD represents a host-shift event (potentially from foodborne contamination, as suggested by Lafferty 2025), initial virulence might have been high (maladapted to new host).

## Literature Evidence

### Theoretical Framework: Anderson-May Trade-off Theory

The virulence-transmission trade-off was formalized by Anderson & May (1982) and refined by Alizon et al. (2009). Key theoretical predictions:

1. **R₀ maximization drives evolution**: Pathogens evolve to maximize basic reproductive number R₀ = β·S₀/(γ+μ)
2. **Trade-off shape determines ESS**: Convex trade-offs (alpha_kill > alpha_shed) favor intermediate virulence; linear trade-offs favor extreme strategies  
3. **Host lifespan matters**: Shorter-lived hosts select for higher virulence (less time penalty for killing host)
4. **Transmission mode affects optimum**: Direct contact vs. environmental transmission alter selective pressures

Recent meta-analysis (Cressler et al. 2019, Evolution) confirmed trade-offs exist across diverse pathogen taxa, with the strength and shape varying by system.

### Marine Disease Context

**Lafferty (2017)** provides the foundational framework for marine infectious disease ecology, emphasizing that marine systems differ from terrestrial:

- **Waterborne transmission** creates different dynamics than direct contact
- **Three-dimensional habitat** allows pathogen dispersal via ocean currents
- **Temperature drives virulence**: Warming simultaneously increases pathogen virulence and transmission rates
- **Host mobility** in marine systems can facilitate or hinder disease spread depending on spatial heterogeneity

The **SIRP model** (Giménez-Romero et al. 2021) for Pinna nobilis demonstrates that marine pathogen dynamics can be reduced to standard SIR framework when pathogen decay rates are fast relative to host demographic rates. This suggests our power-law virulence trade-offs remain valid for marine bacterial pathogens.

### Bacterial Virulence Evolution

Studies of bacterial pathogen evolution (reviewed in FEMS Microbiology Reviews, Brüssow 2017) show:

- **Rapid phenotypic evolution**: Bacteria can evolve virulence traits within weeks to months
- **Phenotypic mutation step sizes**: Typically 0.01-0.1 in units of fitness effect for quantitative traits
- **Virulence costs**: Higher virulence often reduces competitive ability in resource-limited environments
- **Transmission-virulence coupling**: Toxin production increases both virulence and transmission for many bacterial pathogens

### Vibrio-Specific Evidence

**Zhong et al. (2025)** sequenced V. pectenicida and identified aerolysin-like toxin genes, providing molecular basis for virulence. Aerolysin creates pores in host cell membranes, causing tissue damage consistent with SSWD wasting phenotype. This suggests virulence is mediated by toxin expression levels, which can evolve rapidly through regulatory mutations.

Marine Vibrio species show **temperature-dependent virulence** (Lupo et al. 2020), with warmer waters increasing both virulence and transmission rates. This supports temperature-dependent parameterization of our trade-off curves.

### Comparative Systems

**Clement et al. (2024)** developed an individual-based eco-evolutionary model for Tasmanian devil facial tumor disease (DFTD) that provides the closest methodological parallel to our work. Key insights:

- **Multi-trait coevolution**: Both host resistance and pathogen virulence evolve simultaneously
- **Demographic-genetic feedback**: Population bottlenecks reduce genetic variance, slowing adaptive responses
- **Spatial structure matters**: Connected populations maintain more evolutionary potential than isolated refugia
- **Long-term coexistence**: Eco-evolutionary feedbacks can prevent extinction and promote host-pathogen coexistence

The DFTD system has similarly devastating initial impacts (~90% mortality) but populations persist through evolutionary rescue mechanisms.

## Knowledge Gaps

1. **No direct measurements** of V. pectenicida virulence trade-offs exist
2. **Mutation rate data sparse**: No estimates of phenotypic mutation step sizes for marine bacteria
3. **Environmental virulence unknown**: How does virulence change with temperature, salinity, nutrients?
4. **Multi-host evolution**: V. pectenicida infects 20+ sea star species—how does host diversity affect virulence evolution?
5. **Strain diversity**: Current model assumes single evolving strain; real populations may have multiple competing lineages

## Parameter Recommendations

### alpha_kill = 2.0 (Current Default)
**Rationale:** Mortality costs should accelerate faster than linear to create intermediate ESS. Value consistent with theoretical models and maintains convex trade-off shape.
**Uncertainty:** MEDIUM. Reasonable based on theory, but no direct empirical validation for V. pectenicida.

### alpha_shed = 1.5 (Current Default)  
**Rationale:** Transmission benefits should scale sub-linearly with virulence to balance mortality costs. Creates moderate convexity in trade-off (alpha_kill/alpha_shed = 1.33).
**Uncertainty:** HIGH. Critical parameter that determines ESS, but purely theoretical estimate.

### alpha_prog = 1.0 (Current Default)
**Rationale:** Disease progression scales linearly with virulence—simple assumption for complex physiological process.
**Uncertainty:** HIGH. Progression dynamics poorly understood for SSWD pathophysiology.

### gamma_early = 0.3 (Current Default)
**Rationale:** I₁ individuals shed at 30% of I₂ rate. Balances stealth transmission with symptomatic shedding. Consistent with many bacterial infections having asymptomatic shedding phase.
**Uncertainty:** MEDIUM. Reasonable biological assumption, but no direct evidence for V. pectenicida.

### sigma_v_mutation = 0.02 (Current Default)
**Rationale:** 2% phenotypic step size per transmission event. Conservative estimate allowing gradual evolution without overwhelming drift. Consistent with bacterial experimental evolution studies.
**Uncertainty:** MEDIUM. Order of magnitude likely correct based on bacterial evolution literature.

### v_init = 0.5 (Current Default)
**Rationale:** Moderate initial virulence at outbreak start. If SSWD was host-shift event, pathogen may not have been optimally adapted to sea star hosts initially.
**Uncertainty:** HIGH. No empirical basis for 2013 virulence level. Sensitivity analysis essential.

## Research Priorities

1. **Experimental evolution studies**: Laboratory evolution of V. pectenicida in sea star tissue culture to measure trade-off parameters
2. **Temperature gradients**: Measure how virulence-transmission trade-offs change with temperature
3. **Multi-host evolution**: Determine how virulence evolution differs across sea star species
4. **Strain surveillance**: Molecular monitoring of V. pectenicida diversity in wild populations
5. **Comparative parameterization**: Use Clement et al. (2024) DFTD parameters as calibration reference

## Conclusions

The pathogen evolution module implements theoretically sound frameworks developed over 40+ years of evolutionary epidemiology research. While direct empirical validation is currently impossible due to lack of V. pectenicida-specific data, the parameter values are reasonable based on:

- **Strong theoretical foundation** (Anderson-May framework)
- **Comparative evidence** from other bacterial pathogens  
- **Marine disease ecology** principles (Lafferty 2017)
- **Methodological precedent** from similar eco-evolutionary models (Clement et al. 2024)

The current parameterization provides a starting point for exploring eco-evolutionary dynamics in SSWD. **Sensitivity analyses are critical** given the high uncertainty in several key parameters (alpha_shed, v_init, alpha_prog). Future empirical work should focus on constraining trade-off curve parameters through experimental evolution approaches.

---

**Literature Cited:** 15 sources spanning evolutionary epidemiology theory, marine disease ecology, bacterial evolution, and comparative eco-evolutionary models.