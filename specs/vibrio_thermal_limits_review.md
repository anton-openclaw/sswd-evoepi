# Vibrio Thermal Adaptation Limits — Literature Review

**Date:** 2026-03-06  
**Purpose:** Constrain T_vbnc_min parameter for SSWD-EvoEpi model

---

## Key Question

Can V. pectenicida communities evolve to be pathogenically active at temperatures
well below their ancestral range (e.g., from 12°C down to 4°C) within ~13 years?

**Short answer: Probably not. There's likely a hard floor around 8-10°C for pathogenic activity.**

---

## 1. V. pectenicida Growth Range

### Original Species Description (Lambert et al. 1998)
- Isolated from scallop (*Pecten maximus*) larvae in a French hatchery
- Growth range: **15-37°C** (original description)
- DOI: 10.1099/00207713-48-2-481

### ABIS Encyclopedia / BacDive Data
- Growth range: **4-30°C**, optimal 20°C
- Facultatively anaerobic
- This wider range (4-30°C) likely reflects SURVIVAL not active pathogenesis
- The discrepancy (15°C min in Lambert vs 4°C in ABIS) is critical: **growth ≠ pathogenicity**

### FHCF-3 Strain (Prentice et al. 2025, Nature E&E)
- Cultured at ~21°C (room temperature) for 5-7 days
- Koch's postulates experiments conducted at ambient seawater temperatures
  at Friday Harbor Labs (Salish Sea summer, likely ~11-14°C)
- **No published data on FHCF-3 growth at temperatures <10°C**
- Genome contains 3 aerolysin-like toxin genes — these are protein-based
  virulence factors whose expression is temperature-dependent

### Hewson 2025 (bioRxiv)
- V. pectenicida FHCF-3 found NOT in non-Pycnopodia echinoderms
- Suggests host-specificity, not universal marine pathogen
- No temperature-specific data reported

---

## 2. V. splendidus Clade — Cold Adaptation

V. pectenicida is in the splendidus clade. Related species give us thermal context:

### V. splendidus Strain JZ6 (Liu et al. 2013; Zhang et al. 2016)
- **MOST virulent at 10°C, virulence LOST at 28°C** — reversed temperature pattern!
- Pathogenic to Yesso scallop at 10°C, isolated from winter mortality
- This is a COLD-ADAPTED pathogenic Vibrio in the same clade
- Growth range in literature: 10-22°C (Armada et al. 2003)
- At 10°C: upregulated adhesion, quorum sensing, metalloprotease virulence factors
- **Key insight: some splendidus clade members are ALREADY cold-adapted pathogens**

### V. tasmaniensis
- Growth range: 8-23°C, optimum 28°C (Thompson et al. 2003)
- VBNC at 20°C after 157 days in nutrient-poor media

### V. salmonicida  
- Growth range: 4-20°C, optimum 15°C
- **True psychrophilic Vibrio** — causes cold-water vibriosis in salmon
- Proves Vibrio CAN be primarily cold-adapted pathogens

---

## 3. VBNC State and Temperature Thresholds

### General Vibrio VBNC Threshold
- Most pathogenic Vibrio enter VBNC below **~10°C** (multiple reviews)
- VBNC state = alive but not culturable, dormant, reduced metabolic activity
- **VBNC cells retain virulence genes** but are NOT actively pathogenic
- Resuscitation occurs upon temperature upshift (typically >15°C)

### Key Distinction: Three Temperature Zones

| Zone | Temperature | Bacterial State | Pathogenic? |
|------|-------------|----------------|-------------|
| **Active** | >15°C (most Vibrio) | Growing, dividing, shedding toxins | YES |
| **Marginal** | 10-15°C | Slow growth, reduced virulence expression | REDUCED |
| **VBNC/Dormant** | <10°C | Alive but non-replicating, dormant | NO |

### Evidence for ~10°C as a Universal Vibrio Threshold
- V. vulnificus: VBNC below 10°C (multiple studies, Oliver 2005 review)
- V. cholerae: VBNC below 10°C (Colwell 2000)  
- V. splendidus: temperature did NOT significantly affect growth in 4-22°C range
  (Armada et al. 2003) — but this is GROWTH, not virulence
- V. parahaemolyticus: can grow at 8.3°C but pathogenic above ~15°C

---

## 4. Can Vibrio Communities Evolve Cold Tolerance?

### Experimental Evolution Evidence

**Cohen et al. 2019 (Evolution)**: V. fischeri evolved at 8°C for **2,000 generations**
- DID adapt to cold stress
- But evolved lines showed fitness trade-offs in symbiosis
- 2,000 generations of Vibrio ≈ ~3-5 years (generation time ~12-24 hours)
- **Conclusion: adaptation is possible but requires thousands of generations
  and comes with fitness costs**

### Constraints on Cold Adaptation

1. **Membrane fluidity**: Cold requires unsaturated fatty acids in membranes.
   This is a fundamental biochemical constraint — there are limits to how fluid
   membranes can be made before they lose integrity.

2. **Enzyme kinetics**: Virulence factors (proteases, toxins) have temperature
   optima. Evolving cold-active enzymes requires substantial protein engineering
   — this takes many more generations than growth adaptation.

3. **Virulence factor expression**: The V. splendidus JZ6 study showed that
   cold virulence requires specific upregulation of adhesion + quorum sensing +
   secretion pathways simultaneously. This is a complex regulatory rewiring,
   not a single mutation.

4. **Trade-off with warm performance**: Cold-adapted enzymes lose activity
   at warm temperatures (psychrophile trade-off). A community that evolves
   cold activity likely loses warm-water virulence.

---

## 5. SSWD and Temperature — Field Evidence

### Kohl et al. 2016 (PLOS ONE) — Pisaster experiments
- SSWD progressed **faster at summer temperatures (~14°C)** than winter (~10°C)
- **Cooler temperature did NOT prevent mortality** — just slowed it
- All stars died regardless of temperature
- Suggests V. pectenicida CAN kill at ~10°C, just more slowly

### Field Observations
- SSWD outbreaks predominantly in **summer/fall** when water is warmest
- Previous outbreaks (1978, 1982-83, 1997) all associated with warm events
- 2013-14 outbreak unusual: **continued through winter** (Hakai Magazine)
  — but this was the initial massive epidemic, not ongoing chronic disease
- Aleutian Islands outbreak reached in **January 2017** — coldest water
  in the range. Suggests pathogen CAN function at cold temps during
  initial invasion, but chronic maintenance may differ.

---

## 6. Implications for Model T_vbnc_min

### What the Literature Supports

| Parameter | Value | Justification |
|-----------|-------|---------------|
| T_vbnc (ancestral) | 12°C | Reasonable: most Vibrio pathogenicity drops sharply below ~12-15°C |
| T_vbnc_min (hard floor) | **8-10°C** | Below 10°C, most Vibrio enter VBNC. Even cold-adapted V. splendidus JZ6 is only pathogenic AT 10°C, not below it. V. pectenicida grows at 4°C but growth ≠ pathogenicity |
| Adaptation rate | **Very slow** | 2000+ generations (3-5 yrs) needed for cold adaptation in V. fischeri. BUT: community-level selection is faster than clonal evolution because it's selecting from standing variation in the environmental pool |
| Adaptation cost | **Yes** | Cold adaptation likely reduces warm-water virulence (trade-off). Community can't be simultaneously optimized for 4°C and 15°C |

### Recommended Model Parameters

**T_vbnc_min = 9°C** (compromise)
- Below 10°C: most Vibrio VBNC (dormant, non-pathogenic)
- V. splendidus JZ6 peak virulence AT 10°C — this is the cold extreme
- Allow adaptation to 9°C but not below — represents the biophysical
  limit of cold-active virulence factor expression

**Alternative: T_vbnc_min = 10°C** (conservative)
- Hard floor at the well-documented VBNC threshold
- Would give Alaska (winter 3-5°C) ~4-5 months of complete disease-free time

### What This Means for the Gradient

With T_vbnc_min = 9-10°C instead of 4°C:
- **AK-PWS** (winter 3-5°C): Pathogen can NEVER be active in winter regardless
  of adaptation → 3-4 month seasonal refuge → recovery window
- **BC-N** (winter 6-8°C): Pathogen mostly inactive in winter → some recovery
- **SS-S** (winter 8-10°C): Pathogen marginal in winter → slow chronic disease
- **CA-S** (winter 12-15°C): Pathogen active year-round → no recovery

The seasonal refuge from a higher T_vbnc_min floor is exactly what the hosts
need to win the evolutionary arms race in northern waters.

---

## References

1. Lambert C et al. (1998) V. pectenicida sp. nov. IJSB 48:481-487
2. Prentice MB et al. (2025) V. pectenicida FHCF-3 is a causative agent of SSWD. Nature E&E
3. Zhang et al. (2016) Comparative transcriptome of V. splendidus JZ6. Appl Environ Microbiol 82:2050-2061
4. Liu et al. (2013) Pathogenic V. splendidus from Yesso scallop. J Invertebr Pathol 114:120-125
5. Armada et al. (2003) V. splendidus growth at different temps. J Appl Microbiol
6. Cohen et al. (2019) Adaptation to temperature stress by V. fischeri. Evolution 73:1885-1897
7. Kohl et al. (2016) Decreased temperature facilitates short-term SSWD survival. PLOS ONE
8. Vattakaven et al. (2006) VBNC state in Vibrio. Int J Food Microbiol
9. Colwell (2000) V. cholerae VBNC state. Int J Food Microbiol
10. Al-Habsi et al. (2022) Temperature upshift review. Front Mar Sci 9:959830
11. Hewson (2025) Autecological insights into V. pectenicida. bioRxiv
