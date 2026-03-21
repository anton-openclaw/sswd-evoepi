# Growth & Life History Parameters Literature Review

## Summary

This document reviews the literature basis for 4 growth and life history parameters in the SSWD-EvoEpi model:

1. **k_growth** — Von Bertalanffy growth rate (yr⁻¹), range: 0.03–0.15
2. **L_min_repro** — Min reproductive size (mm arm radius), range: 200–500
3. **senescence_age** — Senescence onset age (yr), range: 20–80
4. **alpha_srs** — Size-recruitment survival Pareto shape, range: 1.0–1.8

## Parameter 1: k_growth (Von Bertalanffy Growth Rate)

**Range:** 0.03–0.15 yr⁻¹  
**Current default:** 0.08 yr⁻¹  
**Confidence:** ★☆☆ (Low)

### First Principles
The von Bertalanffy growth constant k determines how rapidly an organism approaches its asymptotic size L_inf. At k=0.05 yr⁻¹, reaching 95% of L_inf takes ~60 years, while at k=0.15 yr⁻¹, only ~20 years are required. For *Pycnopodia*, the world's largest sea star (up to 650 mm arm radius), we expect relatively slow growth given the large body size and long lifespan typical of echinoderms. However, captive juveniles show rapid early growth, suggesting possible biphasic growth patterns.

### Literature Evidence
Direct growth rate data for *Pycnopodia helianthoides* are extremely limited. No von Bertalanffy parameters have been published specifically for this species.

**Comparative Evidence from Other Echinoderms:**
- **Arctic brittle star (*Ophiopleura borealis*):** k = 0.01–0.09 yr⁻¹ with 25-32 year lifespans (Frontiers Marine Science, 2025)
- **Mediterranean sea star (*Astropecten aranciacus*):** von Bertalanffy analysis applied successfully (Helgoland Marine Research, 2016)
- **General echinoderm pattern:** Growth is typically slow relative to other marine taxa

**Indirect Evidence:**
- *Pycnopodia* reaches maximum sizes of 40-65 cm arm radius (400-650 mm) based on field observations (Race Rocks Ecological Reserve)
- Captive breeding programs show juveniles can grow rapidly under optimal conditions (Hodin et al. 2021)
- The species' pre-SSWD populations were stable, indicating growth and mortality rates were balanced

**Constraint from Mortality-Growth Relationship:**
Comparative studies suggest von Bertalanffy growth constant K correlates with adult mortality rate M, with K/M ≈ 1.0 (American Naturalist, 1992). If adult mortality is ~0.05-0.10 yr⁻¹, then k should be similar.

### Recommendations
The current range (0.03-0.15 yr⁻¹) is reasonable but poorly constrained. The lower bound reflects the slow growth expected for large, long-lived echinoderms, while the upper bound allows for faster growth observed in captive juveniles. Priority should be given to growth studies in captive populations where age is known.

## Parameter 2: L_min_repro (Minimum Reproductive Size)

**Range:** 200–500 mm arm radius  
**Current default:** 400 mm arm radius  
**Confidence:** ★☆☆ (Low)

### First Principles
Reproductive maturity requires sufficient body mass to support gametogenesis while maintaining somatic functions. For broadcast spawners like *Pycnopodia*, eggs are energetically expensive and require substantial energy reserves. The L_min_repro/L_inf ratio typically ranges from 0.3-0.6 for marine invertebrates, indicating reproduction begins at 30-60% of maximum size.

### Literature Evidence
No direct measurements of size at sexual maturity exist for *Pycnopodia helianthoides*.

**Reproductive Biology Context:**
- Reproductive season: March-July (Animal Diversity Web)
- Broadcast spawning with external fertilization
- Juveniles begin with 5 arms, developing up to 24 arms as adults
- Maximum sizes: 400-650 mm arm radius in nature

**Indirect Size Estimates:**
Given L_inf ≈ 1000 mm (model parameter) and natural maximum sizes of 650 mm, a range of 200-500 mm represents 20-50% of model L_inf, which is consistent with general patterns for marine invertebrates.

**Captive Breeding Context:**
Current captive breeding programs (Hodin et al. 2021, AZA SAFE Program 2024) provide opportunities to directly observe size at first reproduction, but such data are not yet published.

### Recommendations
The 200-500 mm range is reasonable based on general biological principles, but lacks empirical validation. The current default of 400 mm (40% of L_inf) falls within expected ranges for marine broadcast spawners. Captive breeding programs should prioritize documenting size at sexual maturity.

## Parameter 3: senescence_age (Senescence Onset Age)

**Range:** 20–80 years  
**Current default:** 50 years  
**Confidence:** ★☆☆ (Low)

### First Principles
Echinoderms are famous for showing "negligible senescence" — no age-related increase in mortality or decline in physiological function. Many species exhibit indeterminate growth, lifelong reproduction, and extreme longevity. This fundamentally challenges the concept of a discrete "senescence age" for *Pycnopodia*.

### Literature Evidence
**Echinoderm Longevity and Senescence:**
- Red sea urchin (*Strongylocentrotus franciscanus*): >100 year lifespan with negligible senescence (Cell Reports, 2024)
- Both long- and short-lived sea urchin species demonstrate negligible senescence (Genes, 2020)
- Echinoderms maintain regenerative capacity throughout life
- No age-related increase in mortality or disease susceptibility documented

**Implications for Pycnopodia:**
- No direct aging studies exist for *Pycnopodia*
- No growth rings or other aging structures available (unlike fish otoliths)
- Pre-SSWD populations included large, presumably old individuals
- The concept of discrete senescence may not apply to this species

**Model Implementation Context:**
The model implements senescence as increased mortality beginning at senescence_age. For echinoderms, this may be biologically inappropriate, but some cutoff may be necessary for computational tractability.

### Recommendations
The biological evidence suggests echinoderms may not undergo classic senescence. The 20-80 year range encompasses uncertainty, but even the lower bound may be too conservative. Consider alternative formulations such as very gradual age-related mortality increases or eliminate senescence entirely, relying only on background mortality and disease.

## Parameter 4: alpha_srs (Size-Recruitment Survival Pareto Shape)

**Range:** 1.0–1.8  
**Current default:** 1.35  
**Confidence:** ★★☆ (Medium)

### First Principles
Size-selective mortality is universal in marine recruitment. Larger settlers have advantages including: (1) higher energy reserves for post-settlement survival, (2) reduced vulnerability to size-limited predators, (3) better competitive ability for space and resources, and (4) enhanced physiological buffering capacity. The Pareto shape parameter α>1 creates increasingly steep survival advantages with size.

### Literature Evidence
**Direct Evidence:**
No size-recruitment survival curves have been published for *Pycnopodia* settlers.

**Comparative Evidence:**
- Size-selective mortality well-documented in marine invertebrate recruitment
- Larger larvae/settlers consistently show higher survival rates across taxa
- Effect sizes vary by species, habitat, and predator assemblage

**Genetic Context:**
The model includes genetic effects on larval size through parental traits. This creates a realistic link between adult genetics, offspring size, and recruitment success, providing selective pressure for larger size.

**Interaction with Settlement Success:**
Parameter interacts with settler_survival (B-H s₀ = 0.03) to determine overall recruitment rates. The combination must produce realistic population replacement rates.

### Recommendations
The range 1.0-1.8 captures reasonable uncertainty around size-selective recruitment mortality. α=1.0 represents no size effect, while α=1.8 creates strong size advantages. The default α=1.35 represents moderate size selection, which is biologically plausible given the importance of size in marine recruitment.

## Cross-Parameter Constraints

### Growth-Mortality Balance
Pre-SSWD *Pycnopodia* populations were stable, requiring balance between:
- Growth rate (k_growth)
- Age at maturity (function of k_growth and L_min_repro)
- Adult mortality (related to senescence_age)
- Recruitment success (alpha_srs, settler_survival)

### Size Structure Effects
Moritsch (2018) demonstrated that functional recovery requires not just population recovery but size structure recovery. Smaller individuals provide less predation pressure, meaning population counts alone are insufficient metrics.

### Population Viability
The NOAA Status Review (Lowry et al. 2022) and associated Population Viability Analysis provide the most comprehensive demographic data available, though they focus on decline rates rather than growth parameters.

## Knowledge Gaps and Research Priorities

1. **Direct growth measurements** from captive populations where age is known
2. **Size at sexual maturity** from breeding programs
3. **Maximum lifespan estimates** through alternative aging methods
4. **Size-selective recruitment survival** from settlement studies
5. **Growth rate variation** across environmental conditions and life stages

## Literature Cited

- Frontiers in Marine Science (2025). Slow growth and high longevity characterize the common, large Arctic brittle star, Ophiopleura borealis.
- Cell Reports (2024). Genomic signatures of exceptional longevity and negligible aging in the long-lived red sea urchin.
- Genes (2020). Senescence and Longevity of Sea Urchins.
- Hodin, J. et al. (2021). Progress Toward Complete Life-Cycle Culturing of the Endangered Sunflower Star, Pycnopodia helianthoides. Biological Bulletin, 241(3), 243-258.
- Lowry, L. et al. (2022). ESA Status Review Report: Sunflower Sea Star. NMFS, 89 pp.
- Moritsch, M.M. (2018). Ecological causes and consequences of Sea Star Wasting Syndrome on the Pacific coast. PhD dissertation, UC Santa Cruz.
- Race Rocks Ecological Reserve. Pycnopodia heliathoides: Sunflower star.
- The American Naturalist (1992). Patterns of Survival, Growth, and Maturation in Snakes and Lizards.