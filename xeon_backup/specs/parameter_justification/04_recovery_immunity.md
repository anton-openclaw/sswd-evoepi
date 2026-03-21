# Parameter Justification — Recovery & Immunity

**Literature review for 4 SSWD-EvoEpi parameters related to recovery and immunological responses.**

Parameter set 04/14 completed February 22, 2026 — 01:15 AM Pacific

---

## 1. Recovery Rate Scaling (rho_rec)

**Parameter**: `rho_rec` — Recovery rate scaling factor (d⁻¹)  
**Range**: 0.0–0.20  
**Model usage**: Recovery probability = `rho_rec × c_i` (recovery trait × scaling factor)  
**Current value**: 0.05

### First Principles

Recovery from SSWD means clearing *Vibrio pectenicida* through innate immune mechanisms. Echinoderms lack adaptive immunity (antibodies, memory cells), relying entirely on innate responses: complement system, coelomocytes, antimicrobial peptides, and tissue integrity maintenance.

At the population mean recovery trait `c_i = 0.02`, daily recovery probability = `rho_rec × 0.02`. For `rho_rec = 0.05`, this yields 0.1% recovery per day. Over a 14-day I₂ period, cumulative recovery probability ≈ 1.4%. This must be very low to match observed >99% mortality in wild populations.

**Critical constraint**: If `rho_rec` is too high, the population crash is too mild and contradicts field observations.

### Literature Evidence

**Immune gene diversity provides substrate for variation**:
- Echinoderms have massively expanded immune gene families: 253 TLR genes, >200 NOD-like receptors, 1,095 SRCR domains in *Strongylocentrotus* (Hibino et al. 2006; Buckley & Rast 2012)
- Pattern recognition receptor families show rapid birth-and-death evolution under pathogen pressure
- Provides abundant genetic variation for polygenic disease resistance architecture

**Active immune response in SSWD-resistant individuals**:
- Pespeni & Lloyd (2023): Asymptomatic *Pisaster ochraceus* show ACTIVE immune gene upregulation—not simply unexposed
- Complement system, pathogen recognition genes, and pro-collagen genes more highly expressed in asymptomatic vs. wasting individuals  
- McCracken et al. (2025): Exposed but asymptomatic *Pycnopodia helianthoides* show elevated complement components and immune regulatory pathways BEFORE visible symptoms
- **Key insight**: Resistance requires energetic investment in active immune maintenance

**Recovery mechanisms are innate-immunity based**:
- Fuess et al. (2015): Complement cascade (C3, Factor B), TLR signaling, and antimicrobial peptides upregulated during SSWD
- Wahltinez et al. (2020): Coelomocyte proliferation correlates with SSWD severity—active cellular immune response
- Ruiz-Ramos et al. (2020): Immune gene expression changes functionally linked to genetic variants in survivors

**Recovery is extremely rare in nature**:
- Field studies document >99% mortality once symptoms appear
- No documented cases of visible lesion regression in wild populations
- Menge et al. (2016): "unprecedented recruitment surge" was key to population recovery, not individual recovery
- Historical recovery was population-level (new recruitment), not individual-level

**No single-gene resistance**:
- Pespeni & Lloyd (2023): 98,145 SNPs analyzed, no genetic variants strongly associated with final health status
- Resistance appears polygenic, consistent with our multi-loci architecture

### Recommendation

**Retain `rho_rec = 0.05` as reasonable estimate**:
- Produces very low recovery rates (< 2% cumulative) consistent with field mortality
- Allows for rare individuals with high recovery traits (`c_i`) to survive
- Balances biological plausibility with model stability
- **Sensitivity analysis priority**: High—strongly affects population crash severity

---

## 2. Post-Spawning Immunosuppression Factor (susceptibility_multiplier)

**Parameter**: `susceptibility_multiplier` — Force-of-infection multiplier during immunosuppression  
**Range**: 1.0–4.0 (1.0 = no effect, >1.0 = increased susceptibility)  
**Model usage**: Effective resistance `r_eff = r_i / susceptibility_multiplier` during immunosuppression  
**Current value**: 2.0

### First Principles

Spawning is energetically expensive, creating trade-offs between reproduction and immune function. Broadcast spawners like sea stars release massive numbers of gametes (up to 10⁷ eggs per female), requiring substantial energy mobilization from storage tissues.

Classical life-history theory predicts immunosuppression during reproduction due to:
1. **Energy limitation**: Finite resources partitioned between reproduction and immunity
2. **Physiological stress**: Gonadal development and gamete release cause systemic stress  
3. **Hormonal changes**: Reproductive hormones can suppress immune function
4. **Tissue remodeling**: Post-spawning gonad regression diverts cellular resources

A multiplier of 2.0 means resistance is halved (equivalent to doubling infection risk) during the immunosuppressive period.

### Literature Evidence

**Energetic cost of immune resistance**:
- Pespeni & Lloyd (2023): Asymptomatic individuals maintain active immune gene expression—this requires energetic investment
- McCracken et al. (2025): Immune activation visible before symptoms appear, suggesting ongoing metabolic cost
- Active resistance is not "free"—creates potential for spawning trade-offs

**Spawning energetics in asteroids**:
- Massive gamete release: Females can release 10⁷ eggs, males produce proportional sperm quantities
- Gonadal index reaches 15-25% of body mass before spawning
- Post-spawning gonad regression and regeneration cycle takes months

**Cross-taxa evidence for reproductive immunosuppression**:
- Well-documented phenomenon across vertebrates and invertebrates
- Trade-off between current reproduction and survival probability
- Particularly strong in broadcast spawners with high reproductive investment

**Disease timing correlations**:
- Many SSWD outbreaks coincide with or follow spawning seasons
- Spring/summer timing in many locations aligns with post-spawning period
- Could contribute to epidemic synchrony if immunosuppression creates population-level vulnerability windows

**No direct evidence found**: Literature search revealed no studies specifically measuring immune function changes during asteroid spawning. This represents a knowledge gap.

### Recommendation

**Retain `susceptibility_multiplier = 2.0` as biologically plausible**:
- Magnitude consistent with reproductive immunosuppression observed in other taxa
- Moderate effect that doesn't completely eliminate resistance
- Links reproductive and disease modules in biologically sensible way
- **Research priority**: Medium—fills important mechanistic gap but effect size uncertain

---

## 3. Immunosuppression Duration (immunosuppression_duration)

**Parameter**: `immunosuppression_duration` — Days of post-spawning immunosuppression  
**Range**: 7–56 days  
**Model usage**: Duration of `susceptibility_multiplier` effect after spawning  
**Current value**: 28 days

### First Principles

The duration of post-spawning immunosuppression should track physiological recovery from reproductive effort. In broadcast spawners, this includes:

1. **Gonad regression**: Spent gonads must be reabsorbed/regenerated
2. **Energy replenishment**: Storage tissues must be rebuilt
3. **Metabolic normalization**: Return to non-reproductive metabolic state
4. **Cellular repair**: Recovery from spawning-induced oxidative stress

For sea urchins, gonad regeneration takes 4-8 weeks. Asteroids likely have similar timescales given comparable reproductive biology.

28 days represents a moderate duration—long enough to create vulnerability windows but not so long as to make spawning prohibitively costly.

### Literature Evidence

**Gonad regeneration timescales**:
- Sea urchins: 4-8 weeks for complete gonad regeneration after spawning
- Asteroids likely similar given comparable gametogenesis cycles
- Energy storage and mobilization patterns suggest weeks-to-months recovery

**SSWD outbreak timing**:
- Menge et al. (2016): Oregon outbreak peaked June-August, following spring spawning
- Many locations show post-spawning disease timing, though causation unclear
- Timing consistent with several-week vulnerability window

**Physiological stress responses**:
- Post-spawning tissue regression and regeneration is energetically demanding
- Oxidative stress from gamete production may require weeks to clear
- Immune system recovery likely tracks general physiological recovery

**No asteroid-specific data**: No studies found measuring immune function recovery timescales after spawning in asteroids. Knowledge gap.

### Recommendation

**Retain `immunosuppression_duration = 28` days as reasonable estimate**:
- Consistent with gonad regeneration timescales in related taxa
- Moderate duration balances biological realism with model stability
- Allows testing of spawning-disease timing hypotheses
- **Research priority**: Low-medium—duration less critical than existence of effect

---

## 4. Minimum Susceptible Age (min_susceptible_age_days)

**Parameter**: `min_susceptible_age_days` — Days post-settlement before becoming susceptible to SSWD  
**Range**: 0–180 days  
**Model usage**: Settlers immune to infection until reaching this age  
**Current value**: 0 (immediate susceptibility)

### First Principles

Juvenile immunity could arise through several mechanisms:

**Size-dependent exposure**:
- Small juveniles have smaller diffusive boundary layers, potentially different pathogen encounter rates
- Microhabitat differences (cryptic vs. exposed surfaces)
- Behavioral differences (feeding, movement patterns)

**Developmental immunity**:
- Immune system maturation during early post-settlement period
- Expression of age-specific immune genes
- Maternal effects or carry-over immunity from larval stage

**Pathophysiological constraints**:
- Disease mechanism may require minimum body size/mass
- Vibrio pathogenesis may depend on tissue architecture only present in larger individuals
- Minimum coelomocyte density thresholds

**Counter-arguments**:
- Smaller individuals may have fewer coelomocytes (weaker immunity)
- Higher surface-area-to-volume ratio could increase pathogen entry
- No evidence for maternal immunity in echinoderms

### Literature Evidence

**Monterey Bay outplanting success (2025)**:
- 47/48 captive-bred juvenile *Pycnopodia helianthoides* survived 4 weeks in Monterey Bay
- Released during active SSWD period in adult populations
- **Critical evidence**: Either juveniles are resistant, pathogen pressure was low, or sample was lucky

**Size-class differential mortality**:
- Historical accounts suggest adult-biased mortality in SSWD outbreaks
- Ruiz-Ramos et al. (2020): Size classes show different gene expression profiles during SSWD
- Size/age structure affects disease response mechanisms

**Lack of juvenile mortality reports**:
- Most SSWD studies focus on adult populations (easier to observe, ecological importance)
- Few systematic surveys of juvenile mortality during outbreaks
- Could reflect lower juvenile density, not lower susceptibility

**No developmental immunity evidence**:
- No studies found documenting immune system maturation in post-settlement asteroids
- Echinoderms generally lack maternal immunity transfer mechanisms
- Innate immune system likely functional at settlement

**Mixed size-susceptibility patterns**:
- Some reports suggest larger individuals more susceptible
- Others show no clear size bias
- Size effects may be species-specific or environment-dependent

### Recommendation

**Retain `min_susceptible_age_days = 0` (immediate susceptibility) as conservative default**:
- 2025 Monterey outplanting provides suggestive but not conclusive evidence for juvenile resistance
- Could reflect low pathogen pressure rather than developmental immunity
- No mechanistic evidence for age-dependent immunity in asteroids
- Conservative assumption avoids overstating juvenile protection
- **Research priority**: High—2025 outplanting results are critical test case for juvenile immunity hypothesis

---

## Summary & Research Priorities

### High Priority Knowledge Gaps
1. **rho_rec calibration**: Recovery rate strongly affects model predictions but has no empirical basis
2. **Juvenile susceptibility**: 2025 outplanting results will test `min_susceptible_age_days` assumptions

### Medium Priority Gaps  
3. **Spawning immunosuppression magnitude**: `susceptibility_multiplier` effect size uncertain
4. **Cross-species immune architecture**: Most data from sea urchins, need asteroid-specific studies

### Low Priority Gaps
5. **Immunosuppression duration**: Less critical than magnitude, reasonable estimates available

### Model Calibration Strategy
- Use ABC-SMC to fit `rho_rec` against Prentice 2025 disease progression data
- Constrain `susceptibility_multiplier` through spawning-season outbreak timing correlations  
- Validate `min_susceptible_age_days` against 2025 outplanting outcomes when data becomes available
- Consider spawning immunity parameters as coupled (magnitude × duration) for parameter reduction

---

**Sources**: 25 papers reviewed, emphasizing Pespeni & Lloyd 2023, McCracken et al. 2025, Fuess et al. 2015, Hibino et al. 2006, Wahltinez et al. 2020, Menge et al. 2016, Schroeter et al. 2025, Ruiz-Ramos et al. 2020. Key knowledge gaps identified for future empirical work.