# RESEARCH_GAPS.md — Targeted Literature Review (Phase 14)

**Date:** 2026-02-14  
**Author:** Anton 🔬  
**Purpose:** Fill knowledge gaps identified in REVIEW.md via targeted web searches

---

## 1. Searches Performed

15 targeted web searches across 5 categories:
1. SSWD biology & recovery status (3 searches)
2. Population genetics & SRS theory (3 searches)
3. Disease modeling methodology (3 searches)
4. Environmental forcing & heatwaves (2 searches)
5. Conservation & captive breeding (2 searches)
6. Sensitivity analysis methods (2 searches)

---

## 2. Key Findings

### 2.1 Prentice et al. 2025 — Etiology Confirmed

**Citation:** Prentice et al. (2025). "Vibrio pectenicida strain FHCF-3 is a causative agent of sea star wasting disease." *Nature Ecology & Evolution*.  
**DOI:** 10.1038/s41559-025-02797-2  
**Access:** PubMed 40760083 (abstract only; Nature paywall)

**Key findings (from press releases and abstract):**
- Koch's postulates fulfilled: pure V. pectenicida cultures injected into healthy Pycnopodia induced SSWD and mortality
- All exposed sea stars except one (lowest dose) died
- V. pectenicida dominant in coelomic fluid of wasting stars, nearly absent in healthy controls
- Genome (MRA.00287-25) reveals aerolysin-like toxin genes — potential mechanism for tissue destruction

**Implications for model:**
- Confirms V. pectenicida as THE causative agent (not a densovirus or multi-pathogen syndrome)
- Dose-dependent mortality supports our force-of-infection framework
- The single survivor at lowest dose is consistent with individual resistance variation
- Aerolysin-like toxins suggest tissue necrosis mechanism — compatible with our I₂→D progression

**Status:** ⬜ Scrape request (Nature paywall). Priority: HIGH.

### 2.2 Hewson et al. 2025 — Review Synthesis

**Citation:** Hewson, I., Johnson, M.R., & Reyes-Chavez, B. (2025). "Lessons Learned from the Sea Star Wasting Disease Investigation." *Annual Review of Marine Science*, 17(1), 257–279.  
**DOI:** 10.1146/annurev-marine-040623-082617  
**Access:** Annual Reviews paywall

**Key findings (from secondary sources):**
- V. pectenicida NOT detected in non-Pycnopodia species with SSWD — the disease may be multi-agent across species
- Previous SSaDV (densovirus) hypothesis definitively ruled out
- Temperature strongly correlates with SSWD outbreaks — consistent with our Arrhenius framework
- Environmental stressors (marine heatwaves, low oxygen) may predispose to disease

**Implications for model:**
- Our model is correctly species-specific for Pycnopodia
- Multi-species SSWD may involve different pathogens — our V. pectenicida focus is appropriate for sunflower stars
- Environmental stress interactions (hypoxia, MHW) could be future extensions

**Status:** ⬜ Scrape request (Annual Reviews paywall). Priority: HIGH.

### 2.3 Sunflower Star Conservation Progress — December 2025

**Source:** Sunflower Star Laboratory press release (2025-12-11)  
**URL:** https://www.sunflowerstarlab.org/news/sunflower-stars-restoration-research-california  
**Access:** ✅ Fully accessed

**Key findings:**
- First experimental outplanting of captive-reared juvenile Pycnopodia in Monterey, CA
- 47/48 (98%) survived 4 weeks in ocean conditions
- eDNA monitoring tools being validated (Stanford University)
- California Academy of Sciences plans first captive spawn in 2026
- Cryopreserved sperm banks established — hundreds of thousands of larvae, millions of sperm
- Wild sunflower stars spotted in California in summer 2025 (first since 2018!)
- SAFE program, PCOR Initiative, and Pycnopodia Recovery Working Group all active
- Prior trials conducted by Dr. Jason Hodin at Friday Harbor Labs (Willem's institution!)

**Implications for model:**
- Conservation module should model the actual Monterey outplanting protocol: 48 juveniles, 4-week holding
- 98% survival over 4 weeks provides empirical survival parameter for released juveniles
- Cryopreservation means genetic diversity can be maintained across breeding seasons
- eDNA could provide field Vibrio concentration data we currently lack
- Wild sightings in CA suggest very low-level persistence (our model shows N=1 at Monterey — consistent!)
- FHL connection means Willem can directly inform the conservation module design

**Status:** ✅ Saved summary. No PDF needed.

### 2.4 Clement et al. 2024 — Eco-Evolutionary IBM for Disease-Host Coevolution

**Citation:** Clement, D.T. et al. (2024). "Coevolution promotes the coexistence of Tasmanian devils and a fatal, transmissible cancer." *Evolution*, qpae143.  
**DOI:** 10.1093/evolut/qpae143  
**Access:** Oxford Academic paywall; PMC review article available (PMC12459819)

**Key findings (from review & abstract):**
- Individual-based eco-evolutionary model of devil-DFTD coevolution
- Parameterized with ~20 years of devil demography, DFTD epidemiology, and GWAS data
- High probability of devil-DFTD coexistence predicted, with greater devil recovery than purely ecological models
- **Coevolution is the key**: ecological-only models predicted extinction; adding host evolution changed the outcome
- Pathogen attenuation also contributed to coexistence

**Implications for model:**
- Directly validates our approach: eco-evolutionary IBM is the RIGHT framework for this question
- We should consider pathogen evolution (virulence attenuation) as a model extension
- Their GWAS-parameterization is what we should aim for once Pycnopodia GWAS data exists
- Code likely available on Dryad/GitHub — methodological reference for sensitivity analysis

**Status:** ⬜ Scrape request (Oxford Academic paywall). Priority: CRITICAL. Also check for code repository.

### 2.5 GitHub Code: siskavera/tasmanian-devil

**URL:** https://github.com/siskavera/tasmanian-devil  
**Description:** Spatially explicit, individual-based framework to model the spread of DFTD in Tasmania  
**Access:** ✅ Public repository

**Implications:** Reference implementation for spatial IBM disease modeling. Should clone and review for architectural patterns (spatial network, transmission, population dynamics).

**Status:** ⬜ Clone and review.

### 2.6 Eco-Evolutionary IBMs for Marine Adaptation

**Citation:** Cosandey-Godin et al. (2021). "Individual-based eco-evolutionary models for understanding adaptation in changing seas." *Proc R Soc B*, 288: 20212006.  
**PMC:** PMC8580472  
**Access:** ✅ Open access

**Key findings (from abstract):**
- Literature review of eco-evolutionary IBMs in marine systems
- Three themes: (1) genetic architecture effects on adaptation, (2) gene flow facilitating rapid adaptation, (3) multiple stressors
- Genetically explicit IBMs simulate polygenic traits using frameworks like SLiM
- Demonstrates that genetic architecture (number of loci, effect size distribution) significantly affects adaptive capacity predictions
- Provides framework for user exploration

**Implications for model:**
- Confirms that genetic architecture choice (51 vs fewer loci) matters for adaptive predictions — our concern (§3.2 of REVIEW.md) is well-founded
- Gene flow between nodes could accelerate resistance spread — our C matrix handles this
- Multiple stressors (disease + temperature) interactions are important — we have temperature-dependent disease but could add heatwave events
- SLiM is the dominant framework for this type of work — future versions could leverage SLiM for genetics component

**Status:** ✅ Accessed. Summary above.

### 2.7 Sweepstakes Reproductive Success — Ne/N Evidence

**Sources searched:**
- Árnason et al. 2023 (eLife): SRS via pervasive selective sweeps, α ≈ 1.35
- Hedgecock & Pudovkin 2011 review: Ne/N < 10⁻³ in some marine species
- Southern bluefin tuna (Science Advances): Ne/N > 0.1 (NOT SRS) — shows SRS is not universal
- NZ snapper study: SRS absent in protected population — fishing may amplify SRS signal

**Implications for model:**
- Ne/N ≈ 10⁻³ is an extreme scenario, not universal
- For Pycnopodia specifically: no empirical Ne/N estimate exists
- Our α_srs = 1.35 produces heavy SRS but should be tested across [1.0, 2.0]
- Sensitivity analysis should include α_srs as a key parameter
- The snapper study suggests environmental conditions affect SRS strength — annual variation could be important

### 2.8 Evolutionary Rescue Theory

**Sources:**
- Alexander et al. 2014 (Evol Appl): "Evolutionary rescue: linking theory for conservation and medicine"
- Carlson et al. 2014 (Phil Trans B): "Evolutionary rescue in a changing world"
- Bell & Gonzalez 2009 (Evol Appl): Foundational experimental work

**Key theoretical insights:**
- Evolutionary rescue is unlikely for small, genetically depauperate populations of organisms with long generation times
- Pycnopodia: generation time ~10 years, critically low N, SRS reduces Ne further
- Standing genetic variation is more likely to provide rescue than new mutations at μ=10⁻⁸
- SRS is a double-edged sword: reduces Ne (bad for genetic variation) but allows rapid frequency shifts (good for adaptation)
- The critical question is whether per-generation selection differential × initial V_A is sufficient to outpace demographic decline

**Implications for model:**
- Spontaneous evolutionary rescue in Pycnopodia may require 50–100+ generations (500–1000 years)
- This is too slow for conservation relevance — reinforces the importance of the conservation module
- Captive breeding + genetic rescue (importing diverse alleles from distant populations) could accelerate the process
- Our model can test this by comparing: (a) no intervention, (b) captive-bred release, (c) assisted gene flow

### 2.9 Pycnopodia Reference Genome

**Citation:** Schiebelhut et al. (2024). "A reference genome for ecological restoration of the sunflower sea star, Pycnopodia helianthoides." *Journal of Heredity*, esad054.  
**NCBI:** GCA_032158295.1  
**Access:** Journal of Heredity (partially open)

**Key findings:**
- First reference-quality genome for Pycnopodia
- Enables population genomic, comparative genomic analyses
- Will "underwrite" analyses of SSW and environmental stressors
- Foundation for future GWAS studies of SSWD resistance

**Implications:**
- Once GWAS for SSWD resistance is conducted using this genome, we can replace our 52-locus model with empirically-parameterized architecture
- This is the most impactful future dataset for our model

### 2.10 Vibrio Survival in Seawater

**Sources:**
- Kaspar & Tamplin 1993 (AEM): V. vulnificus survival vs temperature/salinity
- Biosca et al. 1999 (AEM): V. vulnificus biotype 2 long-term survival
- V. splendidus temperature/salinity responses (Lacoste 2003)

**Findings:**
- Vibrio survival is ENHANCED at warm temperatures (13–22°C) — our model's counter-intuitive faster decay at cold T is consistent with V. vulnificus data
- At cold temps (< 10°C), Vibrio enters VBNC state — our sigmoidal VBNC model is appropriate
- Salinity: Vibrio needs marine salinity (> 10 psu), dies in freshwater — our salinity modifier is correct
- Biological factors (grazing, competition) reduce Vibrio in unsterilized seawater — our model doesn't include this, but flushing rate serves as a proxy

**Limitation:** No survival/decay data specific to V. pectenicida in seawater. We're using V. vulnificus analogs.

---

## 3. Scrape Request List

Papers to request from Willem or attempt through Unpaywall/Sci-Hub:

| Priority | Citation | DOI | Why We Need It |
|----------|----------|-----|---------------|
| CRITICAL | Clement et al. 2024, Evolution | 10.1093/evolut/qpae143 | Closest methodological analog; GWAS-parameterized eco-evo IBM |
| HIGH | Prentice et al. 2025, Nat Ecol Evol | 10.1038/s41559-025-02797-2 | Koch's postulates, dose-response, V. pectenicida details |
| HIGH | Hewson et al. 2025, ARMS | 10.1146/annurev-marine-040623-082617 | Comprehensive SSWD review, multi-species evidence |
| MEDIUM | Prowse et al. 2016, Ecosphere | 10.1002/ecs2.1238 | Sensitivity analysis protocol for stochastic IBMs |
| MEDIUM | Schiebelhut et al. 2024, J Heredity | 10.1093/jhered/esad054 | Pycnopodia reference genome and pop gen implications |
| LOW | Lundquist & Botsford 2004, Ecol Appl | 10.1890/02-5325 | Original Allee effect model for broadcast spawners |

---

## 4. Code Repositories to Review

| Repo | URL | Relevance |
|------|-----|-----------|
| siskavera/tasmanian-devil | github.com/siskavera/tasmanian-devil | Spatial IBM for DFTD in Tasmania — architectural reference |
| Clement et al. 2024 code | Check Dryad / paper supplementary | GWAS-parameterized eco-evo IBM code |
| SLiM | messerlab.org/slim | Dominant eco-evo simulation framework; potential backend for genetics |

---

## 5. Key Research Gaps Summary

### Gaps that would MOST improve model confidence:

1. **Pycnopodia GWAS for SSWD resistance** — transforms genetic architecture from theoretical to empirical. Waiting on follow-up work using the 2024 reference genome.

2. **Field Vibrio concentration measurements** — would directly calibrate σ_eff, P_env, and ξ(T). Could use eDNA surveys at SSWD-affected sites.

3. **Clement 2024 methodology transfer** — their validation approach (comparison with 20 years of field data) is exactly what we need. Full paper access is critical.

4. **Multi-generational allele frequency data** — pre-SSWD vs survivor genetic comparisons in Pycnopodia would calibrate per-locus shift targets.

5. **Captive breeding release data** — the December 2025 outplanting provides the first empirical datapoint; longer-term survival and reproduction data will follow.

### Gaps that are interesting but lower priority:

6. Marine heatwave forcing data (available from NOAA Blobtracker)
7. Larval dispersal modeling via ROMS particle tracking
8. Urchin barren feedback dynamics
9. Pathogen evolution / virulence attenuation
10. Multi-species SSWD interactions
