# Disease Progression Rate Parameters — Literature Review

**Review Date:** February 22, 2026  
**Reviewer:** Anton 🔬  
**Cron Job:** 5a07cbd2-78af-4cbd-97e3-1fe6041678ff

---

## Overview

This literature review examines three critical disease progression rate parameters in the SSWD-EvoEpi model:

1. **mu_EI1_ref** — E→I₁ progression rate at T_ref (d⁻¹)
2. **mu_I1I2_ref** — I₁→I₂ progression rate at T_ref (d⁻¹)  
3. **mu_I2D_ref** — I₂→Death rate at T_ref (d⁻¹)

These parameters control the temporal dynamics of the S→E→I₁→I₂→D disease cascade for sea star wasting disease (SSWD) caused by *Vibrio pectenicida* in *Pycnopodia helianthoides*.

---

## mu_EI1_ref: Exposed to First Infectious Stage Rate

### Description
Rate at which exposed individuals progress from the incubation period (E) to early symptomatic disease (I₁). This represents the inverse of the mean incubation period—the time from pathogen exposure to first visible symptoms.

### First Principles
- **Mechanistic meaning**: Speed of *V. pectenicida* establishment and initial host tissue damage
- **Physical constraints**: Must be >0; extremely high values (>2.0 d⁻¹) would eliminate incubation period entirely
- **Biological bounds**: Incubation periods <12 hours unlikely (bacterial establishment takes time); >2 weeks may allow immune clearance
- **Disease cascade requirement**: If too fast, E compartment becomes negligible; if too slow, epidemic cannot establish
- **Expected range**: 1/(1-14 days) = 0.07-1.0 d⁻¹

### Literature Evidence
- **Temperature sensitivity established**: Bates (2009) showed 4°C temperature increase sufficient to induce SSWD-like symptoms within 96 hours, indicating rapid progression at elevated temperatures
- **Clinical observations**: Laboratory outbreak described symptoms developing "over the course of a week" (Wikipedia, citing marine lab observations)
- **Mechanism**: McCracken et al. (2025) showed immune system activation and tissue homeostasis disruption **precede** visible wasting symptoms, consistent with an incubation period where pathogen establishment occurs before clinical signs

### Recommended Values
- **Recommended value**: 0.57 d⁻¹ (≈1.8 days mean incubation)
- **SA range**: 0.20-1.00 d⁻¹ (1-5 days mean incubation)
- **Confidence**: MEDIUM — indirect evidence from temperature studies and clinical observations
- **Key sources**: Bates (2009), McCracken et al. (2025), marine laboratory observations

---

## mu_I1I2_ref: First to Second Infectious Stage Rate

### Description
Rate at which individuals progress from early symptomatic disease (I₁) to severe wasting disease (I₂). This captures disease escalation from initial lesions/lethargy to advanced tissue breakdown.

### First Principles
- **Mechanistic meaning**: Speed of *V. pectenicida* proliferation and aerolysin-like toxin damage escalation
- **Physical constraints**: Must be >0; extremely high values would eliminate I₁ stage entirely
- **Biological bounds**: Very rapid progression (>1.0 d⁻¹) inconsistent with observed clinical course; very slow progression (<0.1 d⁻¹) inconsistent with acute nature of SSWD
- **Disease cascade requirement**: I₁→I₂ transition must be faster than recovery to maintain epidemic character
- **Expected range**: Similar to mu_EI1_ref, approximately 0.2-1.0 d⁻¹

### Literature Evidence
- **Temperature dependence**: Kohl et al. (2016) demonstrated that cooler temperatures (9.0°C vs 12.1°C) **slow disease progression** but do not prevent mortality, indicating this parameter is temperature-sensitive
- **Disease stages observed**: Clinical descriptions include progression from lethargy → lesions → tissue breakdown → arm autotomy, consistent with a multi-stage process
- **Pathogen mechanism**: Zhong et al. (2025) identified aerolysin-like toxin genes in *V. pectenicida* FHCF-3, providing mechanism for progressive tissue damage

### Recommended Values
- **Recommended value**: 0.40 d⁻¹ (≈2.5 days mean duration of I₁)
- **SA range**: 0.15-0.80 d⁻¹ (1.25-6.7 days mean I₁ duration)
- **Confidence**: LOW — limited direct evidence, inferred from clinical progression descriptions
- **Key sources**: Kohl et al. (2016), Zhong et al. (2025), clinical stage descriptions

---

## mu_I2D_ref: Second Infectious Stage to Death Rate

### Description
Rate at which severely wasted individuals (I₂) die from SSWD. This represents the final, lethal phase of the disease cascade.

### First Principles
- **Mechanistic meaning**: Speed of terminal organ failure and death from *V. pectenicida* infection
- **Physical constraints**: Must be >0; if too high (>0.5 d⁻¹), I₂ stage becomes very brief; if too low (<0.05 d⁻¹), inconsistent with SSWD lethality
- **Biological bounds**: Terminal phase typically 2-20 days based on general infectious disease patterns
- **Disease cascade requirement**: Must be high enough to generate significant mortality but not so high as to eliminate I₂ compartment
- **Expected range**: 0.05-0.5 d⁻¹ (2-20 days mean survival in I₂)

### Literature Evidence
- **High lethality established**: Kohl et al. (2016) found 100% mortality in both temperature treatments (9.0°C and 12.1°C), confirming SSWD is highly lethal
- **Terminal progression**: Clinical descriptions indicate rapid deterioration once severe wasting begins, with "death and rapid disintegration"
- **Population-level evidence**: Harvell et al. (2019) documented >90% population crashes during 2013-2014 outbreak, indicating very high case fatality rates
- **Recovery rare**: Observations of recovery from severe SSWD are extremely rare in the literature

### Recommended Values
- **Recommended value**: 0.173 d⁻¹ (≈5.8 days mean survival in I₂)
- **SA range**: 0.08-0.35 d⁻¹ (2.9-12.5 days mean I₂ survival)
- **Confidence**: MEDIUM — strong evidence for high lethality, moderate evidence for timing
- **Key sources**: Kohl et al. (2016), Harvell et al. (2019), clinical observations

---

## Synthesis and Model Implications

### Parameter Relationships
The three rates together control the **shape** of the disease time course:
- **Fast mu_EI1_ref + moderate mu_I1I2_ref + slow mu_I2D_ref** → Long infectious period, chronic disease
- **All rates fast** → Acute die-off with brief infectious period  
- **Slow mu_EI1_ref + fast later rates** → Long incubation, rapid terminal progression

### Temperature Dependence
All three parameters show Arrhenius temperature dependence in the model:
- **Established empirically**: Bates (2009), Kohl et al. (2016) demonstrate clear temperature effects
- **Mechanistic basis**: *Vibrio* species are known to be temperature-sensitive; warmer conditions favor bacterial growth and toxin production
- **Model implementation**: Parameters scaled by `arrhenius(rate_ref, Ea, T_celsius)` function

### Uncertainty and Sensitivity
- **Highest uncertainty**: mu_I1I2_ref due to difficulty distinguishing I₁ from I₂ stages in clinical observations
- **Sensitivity analysis priority**: All three parameters identified as important in Morris R4 screening
- **Critical constraint**: Total disease time (1/mu_EI1_ref + 1/mu_I1I2_ref + 1/mu_I2D_ref) must be consistent with observed "days to weeks" progression

### Research Priorities
1. **Controlled infection experiments** with *V. pectenicida* FHCF-3 using Koch's postulates protocol
2. **Time-series sampling** of infected individuals with frequent health assessments
3. **Temperature manipulation studies** to quantify Arrhenius parameters
4. **Biomarker development** to distinguish E, I₁, and I₂ stages objectively

---

## Data Gaps and Limitations

### Critical Gaps
- **No controlled progression studies**: Despite breakthrough pathogen identification (Aquino/Prentice 2025), no published studies track individual disease progression with precise timing
- **Stage definitions subjective**: I₁ vs I₂ distinction based on clinical appearance, not objective biomarkers
- **Species-specific data**: Most quantitative studies on *Pisaster ochraceus*; *Pycnopodia helianthoides* data limited

### Literature Limitations
- **Historical studies** pre-date pathogen identification, limiting mechanistic interpretation
- **Field studies** cannot control for exposure timing, confounding progression rate estimates
- **Laboratory studies** often use already-diseased animals, precluding incubation period measurement

### Model Assumptions
- **Exponential transitions**: Real disease progression may be more complex (sigmoidal, multi-modal)
- **Temperature independence of stages**: Model assumes same Arrhenius parameters for all individuals
- **No recovery from I₂**: Model assumes I₂→D transition is irreversible

---

## Confidence Assessment

| Parameter | Confidence | Reasoning |
|-----------|------------|-----------|
| mu_EI1_ref | MEDIUM | Temperature sensitivity established; approximate timing from field observations |
| mu_I1I2_ref | LOW | Limited ability to distinguish I₁ from I₂ objectively; inferred from progression descriptions |
| mu_I2D_ref | MEDIUM | High lethality well-documented; moderate evidence for terminal phase timing |

**Overall confidence**: MEDIUM-LOW — parameters are reasonable for modeling purposes but would benefit significantly from controlled experimental validation.