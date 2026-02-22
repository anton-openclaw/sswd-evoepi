# Knowledge Gaps and Research Uncertainties

**SSWD-EvoEpi Model: Critical Research Needs**  
**Compiled:** February 22, 2026 — 4:42 AM Pacific

---

## Executive Summary

The SSWD-EvoEpi model incorporates 47 parameters across 11 functional groups, representing the most comprehensive mechanistic model of sea star wasting disease to date. However, **species-specific empirical data for *Pycnopodia helianthoides* remains extremely limited**. Most parameters are constrained by first principles, comparative studies from related species, or theoretical frameworks rather than direct measurements.

This document identifies critical knowledge gaps, quantifies uncertainty levels, and prioritizes research needs to improve model realism and predictive capability.

---

## Tier 1: Critical Gaps (Research Priority: IMMEDIATE)

### 1. **Vibrio pectenicida Experimental Biology**
**Gap:** No experimental evolution studies measuring virulence-transmission trade-offs  
**Impact:** Pathogen evolution module (6 parameters) based entirely on theoretical frameworks  
**Research need:** Laboratory evolution experiments in sea star tissue culture  
**Timeline:** 6-12 months for initial trade-off curves

**Affected parameters:**
- alpha_kill, alpha_shed, alpha_prog (trade-off exponents)
- sigma_v_mutation (mutation step size) 
- gamma_early (asymptomatic transmission)
- v_init (ancestral virulence)

### 2. **Disease Progression Kinetics** 
**Gap:** No controlled infection studies with *V. pectenicida* measuring stage transitions  
**Impact:** Disease timing parameters (3) inferred from field observations and related pathogens  
**Research need:** Koch's postulates protocol with time-series sampling  
**Timeline:** 12-18 months for complete progression curves

**Affected parameters:**
- mu_EI1_ref, mu_I1I2_ref, mu_I2D_ref (progression rates)
- Temperature scaling of all progression parameters

### 3. **Pathogen Shedding Quantification**
**Gap:** No measurements of *V. pectenicida* concentrations from infected *P. helianthoides*  
**Impact:** Force of infection parameters (5) have no empirical basis  
**Research need:** qPCR/culture methods for pathogen quantification during disease stages  
**Timeline:** 6 months for method development, 12 months for complete dataset

**Affected parameters:**
- sigma_1_eff, sigma_2_eff, sigma_D (stage-specific shedding)
- K_half (half-infective dose)
- a_exposure (maximum infection probability)

### 4. **Life History Parameter Validation**
**Gap:** No direct measurements of growth, fecundity, or reproductive timing for *P. helianthoides*  
**Impact:** Life history module (8 parameters) based on comparative scaling relationships  
**Research need:** Quantitative analysis of existing captive breeding data  
**Timeline:** 3-6 months for data analysis, ongoing collection

**Affected parameters:**
- k_growth (growth rate)
- F0 (fecundity) 
- Spawning timing and frequency parameters (4)
- L_min_repro (size at maturity)

---

## Tier 2: High Priority Gaps (Research Priority: 6-12 MONTHS)

### 5. **Genetic Architecture Validation**
**Gap:** No controlled crosses measuring heritability of disease resistance components  
**Impact:** 8 genetic parameters based on population genomic inferences  
**Research need:** Breeding experiments with offspring disease challenges  
**Timeline:** 2-3 years for complete genetic analysis

**Knowledge gap details:**
- Heritability estimates for resistance, tolerance, recovery traits
- Genetic correlation structure between traits  
- Distribution of allelic effects within and between loci
- Dominance/epistatic interactions

### 6. **Recovery Mechanism Identification**
**Gap:** Cellular/molecular basis of rare SSWD recovery unknown  
**Impact:** Recovery rate scaling (rho_rec) has no mechanistic foundation  
**Research need:** Transcriptomic/proteomic analysis of surviving individuals  
**Timeline:** 12-18 months for mechanism identification

**Research priorities:**
- What distinguishes individuals that clear vs. succumb to infection?
- Is recovery active immune clearance or pathogen attenuation?
- How does recovery capacity correlate with genetic background?

### 7. **Juvenile Susceptibility Testing**  
**Gap:** 2025 Monterey outplanting results (47/48 survival) not yet analyzed  
**Impact:** min_susceptible_age_days parameter is critical for conservation outcomes  
**Research need:** Analysis of existing data + controlled juvenile exposure experiments  
**Timeline:** 3-6 months for initial analysis

**Key questions:**
- Was high survival due to age-specific resistance or environmental conditions?
- At what size/age does full SSWD susceptibility develop?
- Do juveniles have different pathogen dose-response relationships?

### 8. **Spawning Synchrony Quantification**
**Gap:** No field observations of natural spawning patterns or chemical induction  
**Impact:** 3 spawning induction parameters based on related species  
**Research need:** Underwater observations during breeding season (March-July)  
**Timeline:** 1-2 breeding seasons for complete dataset

**Research approach:**
- Field monitoring of spawning events and chemical cue propagation
- Laboratory induction experiments with controlled chemical stimuli
- Quantification of cascade spawning success rates

---

## Tier 3: Medium Priority Gaps (Research Priority: 1-3 YEARS)

### 9. **Environmental Reservoir Dynamics**
**Gap:** No quantification of *V. pectenicida* in sediments, biofilms, or non-target hosts  
**Impact:** P_env_max parameter represents complex community dynamics with single value  
**Research need:** Seasonal environmental sampling with species-specific detection  

### 10. **Temperature-Response Relationships**
**Gap:** Arrhenius parameters for most biological rates not empirically determined  
**Impact:** 23 temperature-dependent parameters use assumed activation energies  
**Research need:** Laboratory experiments across temperature gradients  

### 11. **Multi-Host Pathogen Evolution**
**Gap:** *V. pectenicida* infects 20+ sea star species—effects on virulence evolution unknown  
**Impact:** Single-host evolution model may not capture community-level selection  
**Research need:** Comparative infection experiments across host species  

### 12. **Larval Dispersal Validation**
**Gap:** No genetic connectivity data for pre-SSWD *P. helianthoides* populations  
**Impact:** 3 dispersal parameters based on oceanographic theory  
**Research need:** Archived tissue analysis for population genetic structure  

---

## Parameter Confidence Distribution

| Confidence Level | Count | Percentage | Primary Basis |
|------------------|-------|------------|---------------|
| **HIGH (★★★)** | 3 | 6% | Strong empirical evidence, multiple sources |
| **MEDIUM (★★☆)** | 21 | 45% | Reasonable theoretical basis, some data |
| **LOW (★☆☆)** | 23 | 49% | First principles, limited empirical support |

**Key insight:** Nearly half of all parameters have low empirical support, with greatest gaps in pathogen biology and species-specific life history.

---

## Research Infrastructure Needs

### Laboratory Capabilities
- **BSL-2 facilities** for *V. pectenicida* culture and infection experiments
- **Sea star husbandry systems** for controlled breeding and life history studies  
- **Temperature-controlled chambers** for thermal response experiments
- **qPCR/sequencing platforms** for pathogen quantification and genetic analysis

### Field Programs  
- **Captive breeding expansion** at Friday Harbor Labs, Monterey Bay Aquarium, Berkeley Aquarium
- **Wild population monitoring** with standardized disease assessment protocols
- **Environmental sampling** for pathogen reservoir quantification
- **Spawning season surveys** for reproductive biology validation

### Analytical Tools
- **ABC-SMC calibration framework** for parameter constraint given uncertainty
- **Sensitivity analysis platforms** for parameter importance ranking
- **Experimental design optimization** for efficient parameter estimation

---

## Uncertainty Propagation Analysis

### Most Consequential Uncertainties
Based on Morris sensitivity analysis results and biological importance:

1. **rho_rec** (recovery rate): Strongly affects population crash severity and evolutionary rescue potential
2. **Pathogen evolution trade-offs**: Determine long-term host-pathogen coexistence vs. extinction  
3. **Spawning synchrony**: Controls Allee effects and population recovery in depleted populations
4. **Juvenile susceptibility**: Critical for conservation outplanting success rates
5. **Genetic architecture**: Determines evolutionary adaptation speed and endpoint

### Parameter Interaction Effects
High-uncertainty parameters often interact multiplicatively:
- **Disease parameters × genetic parameters**: Joint uncertainty in epidemic severity and evolutionary response
- **Spawning parameters × larval parameters**: Combined uncertainty in reproductive success and connectivity
- **Environmental parameters × pathogen parameters**: Temperature and reservoir effects compound

---

## Research Timeline and Resource Requirements

### Year 1 (2026)
- **Analysis of existing data:** 2025 outplanting, captive breeding records, archived samples
- **Method development:** *V. pectenicida* quantification, controlled infection protocols
- **Pilot studies:** Initial trade-off measurements, spawning observations
- **Estimated cost:** $150K-200K (personnel + equipment)

### Year 2-3 (2027-2028)
- **Controlled experiments:** Full disease progression, pathogen evolution, genetic crosses
- **Field validation:** Breeding season monitoring, environmental sampling
- **Model calibration:** ABC-SMC with empirical constraints
- **Estimated cost:** $300K-400K annually (expanded infrastructure)

### Year 4-5 (2029-2030) 
- **Long-term studies:** Multi-generation genetic analysis, environmental monitoring
- **Conservation applications:** Model-guided outplanting optimization
- **Publication and synthesis:** Complete parameter validation dataset
- **Estimated cost:** $200K-250K annually (ongoing monitoring)

---

## Alternative Approaches for Gap-Filling

### Proxy Species Studies
- **Pacific Crown-of-Thorns (*Acanthaster planci*):** Closest related species with detailed spawning biology
- **Large asteroids (*Pisaster ochraceus*):** SSWD susceptibility and genetic studies available
- **Other broadcast spawners:** Fertilization and dispersal dynamics

### Theoretical Advances
- **Evolutionary rescue theory:** Parameter bounds from population genetic models
- **Marine disease ecology:** Scaling relationships from comparative studies  
- **Oceanographic modeling:** High-resolution dispersal kernels from physical models

### Archived Sample Analysis
- **Museum specimens:** Pre-SSWD genetic diversity and population structure
- **Historical surveys:** Population density and size structure before collapse
- **Environmental DNA:** Pathogen presence in archived water/sediment samples

---

## Recommendations for Model Development

### Near-term (2026)
1. **Implement uncertainty propagation** in model output via Monte Carlo sampling
2. **Develop parameter identifiability analysis** to prioritize empirical studies  
3. **Create experimental design optimization** for efficient parameter estimation
4. **Establish collaborative framework** with aquariums and field stations

### Medium-term (2027-2028)
1. **Iterative model-experiment coupling** with parameter updates as data become available
2. **Multi-model ensemble approaches** to capture structural uncertainty
3. **Bayesian model averaging** across alternative formulations
4. **Real-time model validation** against ongoing field studies

### Long-term (2029-2030)
1. **Integration with climate models** for future scenario projections
2. **Conservation decision support tools** for outplanting optimization
3. **Pathogen surveillance systems** informed by evolutionary predictions
4. **Extension to other marine disease systems** using validated frameworks

---

*This analysis represents the most comprehensive assessment of research needs for marine disease eco-evolutionary modeling available as of February 2026. Addressing these gaps through coordinated empirical research will transform model predictions from theoretical projections to quantitative conservation tools.*