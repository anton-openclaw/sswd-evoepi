# Master Parameter Justification Table

**SSWD-EvoEpi Model: Complete Parameter Documentation**  
**Compiled:** February 22, 2026 — 4:38 AM Pacific  
**47 Parameters across 11 Functional Groups**

---

## Table of Contents
1. [Disease Progression Parameters](#1-disease-progression-parameters)
2. [Pathogen Shedding & Dose-Response](#2-pathogen-shedding--dose-response)
3. [Environmental Pathogen Dynamics](#3-environmental-pathogen-dynamics)
4. [Recovery & Immunity](#4-recovery--immunity)
5. [Growth & Life History](#5-growth--life-history)
6. [Fecundity & Recruitment](#6-fecundity--recruitment)
7. [Genetic Architecture](#7-genetic-architecture)
8. [Spawning Timing](#8-spawning-timing)
9. [Spawning Induction](#9-spawning-induction)
10. [Larval Dispersal](#10-larval-dispersal)
11. [Pathogen Evolution](#11-pathogen-evolution)

---

## 1. Disease Progression Parameters

Controls the temporal dynamics of the S→E→I₁→I₂→D disease cascade for SSWD.

| Parameter | Description | Range | Default | Confidence | Key Literature |
|-----------|-------------|-------|---------|------------|----------------|
| **mu_EI1_ref** | E→I₁ progression rate at T_ref (d⁻¹) | 0.20–1.00 | 0.57 | MEDIUM | Bates (2009), McCracken et al. (2025) |
| **mu_I1I2_ref** | I₁→I₂ progression rate at T_ref (d⁻¹) | 0.15–0.80 | 0.40 | LOW | Kohl et al. (2016), Zhong et al. (2025) |
| **mu_I2D_ref** | I₂→Death rate at T_ref (d⁻¹) | 0.08–0.35 | 0.173 | MEDIUM | Kohl et al. (2016), Harvell et al. (2019) |

### Justification
Temperature sensitivity of disease progression established by Bates (2009) showing 4°C increase sufficient for SSWD-like symptoms within 96 hours. All parameters show Arrhenius temperature dependence. Total disease time (sum of inverse rates) must be consistent with observed "days to weeks" progression. High priority for controlled infection experiments with V. pectenicida.

---

## 2. Pathogen Shedding & Dose-Response

Controls force of infection following Michaelis-Menten kinetics for pathogen dose-response.

| Parameter | Description | Range | Default | Confidence | Key Literature |
|-----------|-------------|-------|---------|------------|----------------|
| **a_exposure** | Maximum daily infection probability (d⁻¹) | 0.30–1.50 | 0.75 | MODERATE | Lafferty (2017), Giménez-Romero (2021) |
| **K_half** | Half-infective dose (bact/mL) | 20,000–200,000 | 87,000 | MODERATE | Marine Vibrio literature |
| **sigma_1_eff** | I₁ shedding rate (field-effective) | 1.0–25.0 | 5.0 | WEAK-MODERATE | McCracken (2023, 2025) |
| **sigma_2_eff** | I₂ shedding rate (field-effective) | 10.0–250.0 | 50.0 | MODERATE | Work et al. (2021), Zhong (2025) |
| **sigma_D** | Saprophytic burst from dead (field-effective) | 3.0–75.0 | 15.0 | WEAK | General decomposition ecology |

### Justification  
Force of infection: λᵢ = a × P/(K_half + P) × (1 - r_eff) × S_sal × f_size(Lᵢ). K_half consistent with marine Vibrio pathogenesis literature (10⁵-10⁶ CFU/mL range). Shedding ratios (σ₂/σ₁ ≈ 10, σ_D/σ₁ ≈ 3) reflect disease pathology progression. Critical gap: quantitative V. pectenicida shedding measurements.

---

## 3. Environmental Pathogen Dynamics

Controls temperature-dependent Vibrio ecology and environmental reservoirs.

| Parameter | Description | Range | Default | Confidence | Key Literature |
|-----------|-------------|-------|---------|------------|----------------|
| **P_env_max** | Background environmental Vibrio input (bact/mL/d) | 50–5000 | 500.0 | MODERATE | Lupo et al. (2020), Aalto et al. (2020) |
| **T_ref** | V. pectenicida optimal temperature (°C) | 17–23 | 20.0 | HIGH | Eisenlord et al. (2016), Bates et al. (2009) |
| **T_vbnc** | VBNC midpoint temperature (°C) | 8–15 | 12.0 | MODERATE | Vibrio VBNC literature, Kohl et al. (2016) |
| **s_min** | Minimum salinity for Vibrio viability (psu) | 5–15 | 10.0 | MODERATE | Marine Vibrio requirements, UW (2014) |

### Justification
T_ref=20°C well-supported: above winter SST (4-8°C), achievable during warm periods (15-21°C), matches V. aquamarinus optimum, consistent with VBNC biology. VBNC transition creates seasonal ON/OFF switch. P_env represents multi-species pathogen maintenance without explicit modeling. Salinity creates spatial boundaries in estuarine systems.

---

## 4. Recovery & Immunity

Controls pathogen clearance and immunological responses, including spawning trade-offs.

| Parameter | Description | Range | Default | Confidence | Key Literature |
|-----------|-------------|-------|---------|------------|----------------|
| **rho_rec** | Recovery rate scaling factor (d⁻¹) | 0.0–0.20 | 0.05 | LOW | Field mortality rates, immune gene studies |
| **susceptibility_multiplier** | Post-spawning immunosuppression factor | 1.0–4.0 | 2.0 | LOW-MODERATE | Life-history theory |
| **immunosuppression_duration** | Duration of post-spawning vulnerability (days) | 7–56 | 28 | LOW-MODERATE | Gonad regeneration timescales |
| **min_susceptible_age_days** | Minimum age for SSWD susceptibility (days) | 0–180 | 0 | MEDIUM | 2025 Monterey outplanting |

### Justification
Recovery extremely rare in field (>99% mortality), but genetic basis exists (Pespeni & Lloyd 2023). Spawning immunosuppression reflects energy trade-offs in broadcast spawners. 2025 Monterey outplanting (47/48 juveniles survived) provides critical test of juvenile immunity hypothesis. Recovery probability = rho_rec × c_i must be very low to match field observations.

---

## 5. Growth & Life History  

Controls individual size-age relationships and demographic rates.

| Parameter | Description | Range | Default | Confidence | Key Literature |
|-----------|-------------|-------|---------|------------|----------------|
| **k_growth** | Von Bertalanffy growth rate (yr⁻¹) | 0.03–0.15 | 0.08 | LOW | Arctic brittle star comparison |
| **L_min_repro** | Minimum reproductive size (mm arm radius) | 200–500 | 400 | LOW | L_min_repro/L_inf ratios |
| **senescence_age** | Senescence onset age (yr) | 20–80 | 50 | LOW | Echinoderm negligible senescence |
| **alpha_srs** | Size-recruitment survival Pareto shape | 1.0–1.8 | 1.35 | MEDIUM | Size-selective mortality theory |

### Justification
Critical knowledge gaps: no direct growth or aging data for P. helianthoides. K/M ≈ 1.0 relationship suggests k ≈ 0.05-0.10 yr⁻¹. Echinoderms show negligible senescence (red sea urchin >100 years), challenging discrete senescence concept. Size-selective recruitment universal in marine systems. Priority: captive breeding growth measurements.

---

## 6. Fecundity & Recruitment

Controls reproductive potential and early life stage survival.

| Parameter | Description | Range | Default | Confidence | Key Literature |
|-----------|-------------|-------|---------|------------|----------------|
| **F0** | Reference fecundity (eggs/female) | 1e6–1e8 | 1e7 | LOW | Crown-of-thorns starfish comparison |
| **gamma_fert** | Fertilization kinetics parameter | 1.0–10.0 | 4.5 | LOW | Lundquist & Botsford (2004) |
| **settler_survival** | Beverton-Holt settler survival s0 | 0.005–0.10 | 0.03 | LOW | Morris SA ranking #6 |

### Justification
F0 based on body size scaling from other large echinoderms (1-100 million range). Fertilization Allee effects critical for crashed populations—gamma_fert controls density threshold steepness. settler_survival absorbs all larval mortality not explicitly modeled. Parameters interact multiplicatively; calibrate as coupled system against pre-SSWD equilibrium.

---

## 7. Genetic Architecture

Controls three-trait polygenic architecture (resistance/tolerance/recovery) and initialization.

| Parameter | Description | Range | Default | Confidence | Key Literature |
|-----------|-------------|-------|---------|------------|----------------|
| **n_resistance** | Number of resistance loci | 5–30 (discrete) | 17 | MODERATE | Schiebelhut et al. (2018) 51 loci |
| **n_tolerance** | Number of tolerance loci | 5–30 (discrete) | 17 | MODERATE | Pespeni & Lloyd (2023) collagen genes |
| **target_mean_r** | Initial population mean resistance | 0.05–0.30 | 0.15 | MODERATE | Population crash severity |
| **target_mean_t** | Initial population mean tolerance | 0.02–0.30 | 0.10 | MODERATE | Constitutive stress responses |
| **target_mean_c** | Initial population mean recovery | 0.02–0.25 | 0.02 | MODERATE | Pathogen-specific clearance |
| **tau_max** | Maximum tolerance mortality reduction | 0.3–0.95 | 0.85 | MODERATE | Tissue repair vs. pathogen damage |
| **q_init_beta_a** | Beta distribution shape α | 1.0–5.0 | 2.5 | LOW | Allele frequency distributions |
| **q_init_beta_b** | Beta distribution shape β | 3.0–15.0 | 8.0 | LOW | Right-skewed for low frequencies |

### Justification
Polygenic small-effect model strongly supported: Burton et al. (2022) and Pespeni & Lloyd (2023) found no major-effect loci in comprehensive screens. Resistance involves active immune gene expression. Three traits have different epidemiological consequences: resistance reduces transmission, tolerance creates silent spreaders, recovery removes infecteds. Initial values balance crash severity with evolutionary potential.

---

## 8. Spawning Timing

Controls seasonal reproductive dynamics and temporal clustering.

| Parameter | Description | Range | Default | Confidence | Key Literature |
|-----------|-------------|-------|---------|------------|----------------|
| **p_spontaneous_female** | Female daily spawning probability (d⁻¹) | 0.005–0.020 | 0.012 | MODERATE | March-July season, May-June peak |
| **p_spontaneous_male** | Male daily spawning probability (d⁻¹) | 0.008–0.020 | 0.0125 | MODERATE | Multiple bouts per season |
| **peak_width_days** | Seasonal peak standard deviation (days) | 30–90 | 60.0 | MODERATE-HIGH | Animal Diversity Web |
| **female_max_bouts** | Maximum female spawning bouts | 1–3 | 2 | LOW | Energetic constraints |

### Justification  
Animal Diversity Web reports spawning "between March and July" with "main peak in May and June". peak_width_days=60 gives ~4-month effective season (±2σ) with 2-month peak window. Values balance participation with synchrony for fertilization success. Critical gap: no quantitative spawning frequency data for P. helianthoides.

---

## 9. Spawning Induction

Controls cascade spawning through chemical cues for reproductive synchrony.

| Parameter | Description | Range | Default | Confidence | Key Literature |
|-----------|-------------|-------|---------|------------|----------------|
| **induction_female_to_male** | Female→male spawning induction probability | 0.50–0.90 | 0.80 | MODERATE | Crown-of-thorns starfish studies |
| **induction_male_to_female** | Male→female spawning induction probability | 0.30–0.70 | 0.60 | MODERATE | Risk-reward asymmetry |
| **readiness_induction_prob** | Social facilitation of gonadal maturation | 0.20–0.80 | 0.50 | LOW | Sea cucumber aggregation |

### Justification
Female→male induction stronger due to: expensive eggs create strong chemical signals, males have more to gain from responding to rare female spawning events. Male→female weaker due to higher female spawning costs requiring selectivity. Readiness induction operates over longer distances (300m vs 200m) and timescales. Critical for reproductive success in depleted populations.

---

## 10. Larval Dispersal

Controls larval connectivity and spatial population structure.

| Parameter | Description | Range | Default | Confidence | Key Literature |
|-----------|-------------|-------|---------|------------|----------------|
| **D_L** | Dispersal scale (km) | 200–600 | 400.0 | MODERATE | PLD × current speed scaling |
| **alpha_self_fjord** | Fjord self-recruitment fraction | 0.20–0.40 | 0.30 | MODERATE | Estuarine circulation |
| **alpha_self_open** | Open coast self-recruitment fraction | 0.05–0.15 | 0.10 | MODERATE | Longshore export |

### Justification
D_L based on PLD (14-70 days) × NE Pacific currents (5-20 cm/s) = ~63-544 km maximum transport. Empirical dispersal typically 10-30% of maximum due to eddies, retention. Fjord vs. open coast differentiation has strong physical basis: estuarine circulation vs. longshore export. 11-node network remains connected with stepping-stone structure preserved.

---

## 11. Pathogen Evolution

Controls virulence evolution through trade-off relationships and mutation dynamics.

| Parameter | Description | Range | Default | Confidence | Key Literature |
|-----------|-------------|-------|---------|------------|----------------|
| **alpha_kill** | Virulence-mortality scaling exponent | 1.5–3.0 | 2.0 | MEDIUM | Anderson-May trade-off theory |
| **alpha_shed** | Virulence-transmission scaling exponent | 1.0–2.5 | 1.5 | HIGH | Convex trade-off requirement |
| **alpha_prog** | Virulence-progression scaling exponent | 0.5–1.5 | 1.0 | HIGH | Linear assumption |
| **gamma_early** | I₁ relative shedding rate | 0.0–0.8 | 0.3 | MEDIUM | Asymptomatic transmission |
| **sigma_v_mutation** | Virulence mutation step size | 0.005–0.05 | 0.02 | MEDIUM | Bacterial evolution studies |
| **v_init** | Initial virulence at outbreak start | 0.2–0.8 | 0.5 | HIGH | Host-shift hypothesis |

### Justification
Implements Anderson-May virulence-transmission trade-off with ESS at intermediate virulence when alpha_kill > alpha_shed. Zhong et al. (2025) identified aerolysin-like toxin genes providing molecular basis. Closest methodological parallel: Clement et al. (2024) DFTD eco-evolutionary model. Critical gap: no V. pectenicida-specific trade-off measurements. Priority: experimental evolution studies.

---

## Research Priority Summary

### Critical Gaps (Immediate Priority)
1. **V. pectenicida experimental evolution** — trade-off parameter measurements
2. **Captive breeding quantification** — growth rates, fecundity, spawning frequencies
3. **2025 outplanting outcome analysis** — juvenile susceptibility validation
4. **Controlled infection experiments** — disease progression rates and shedding dynamics

### High Priority Gaps
5. **Temperature-dependent parameter scaling** — Arrhenius relationships for all biological rates
6. **Genetic architecture empirical validation** — controlled crosses for heritability estimates
7. **Larval dispersal validation** — genetic connectivity patterns in pre-SSWD populations
8. **Spawning synchrony quantification** — field observations of natural spawning patterns

### Medium Priority Gaps
9. **Life history parameter validation** — aging, growth, and reproduction in captivity
10. **Environmental reservoir quantification** — sediment and multi-species Vibrio sources
11. **Recovery mechanism identification** — what determines the rare successful clearance events?
12. **Multi-host virulence evolution** — how does host species diversity affect pathogen adaptation?

---

## Literature Database
**Total sources:** 103+ papers spanning marine disease ecology, echinoderm biology, evolutionary epidemiology, population genetics, and marine connectivity  
**Key foundational papers:** Schiebelhut et al. (2018), Pespeni & Lloyd (2023), Anderson & May (1982), Lafferty (2017), Clement et al. (2024)  
**Species-specific data:** Severely limited — most parameters based on comparative studies and first principles

---

*This master table represents the most comprehensive parameter justification available for the SSWD-EvoEpi model as of February 2026. Parameter confidence levels reflect the strength of available evidence, with many parameters requiring empirical validation through ongoing and planned experimental programs.*