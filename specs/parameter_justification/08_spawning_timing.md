# Parameter Literature Review: Spawning Timing

**Parameters:** p_spontaneous_female, p_spontaneous_male, peak_width_days, female_max_bouts
**Date:** February 22, 2026
**Author:** Anton 🔬

---

## Summary

The spawning timing parameters control the seasonal reproductive dynamics in the SSWD-EvoEpi model. These four parameters determine when and how frequently individuals initiate spawning during the extended breeding season, affecting population recruitment success and disease transmission dynamics.

## Parameter Details

### 1. p_spontaneous_female (Daily Spontaneous Spawning Probability - Females)

**Current Value:** 0.012 d⁻¹
**Range:** 0.005–0.020 d⁻¹

#### First Principles
The daily probability for a reproductively ready female to initiate spawning spontaneously. At p=0.012, expected wait time is ~83 days. Over a 270-day season with Gaussian seasonal modulation, this ensures most females participate in spawning while maintaining temporal clustering for fertilization success. Too low and many females never spawn; too high and spawning becomes completely asynchronous.

#### Literature Evidence

**Direct Evidence:** No species-specific data available for *Pycnopodia helianthoides*. The literature indicates that sea stars exhibit seasonal breeding patterns but lacks quantitative data on daily spawning probabilities.

**Comparative Evidence:** 
- Schiebelhut et al. (2022) found that earlier reproductive seasons (relative to other asteroids) were associated with higher SSWD susceptibility, suggesting that *Pycnopodia* spawning timing may differ from typical patterns
- Animal Diversity Web reports *P. helianthoides* breeds via broadcast fertilization "between March and July" with main peak in "May and June"
- The crown-of-thorns starfish (*Acanthaster planci*) literature emphasizes spawning synchrony as "fundamental for achieving high rates of fertilization" in broadcast spawners (PMC5371309)

**Allee Effect Context:** Lundquist & Botsford (2004) demonstrated that broadcast spawners like *Pycnopodia* experience fertilization success decline at low population densities. This supports the need for moderate spontaneous spawning probabilities that maintain temporal clustering while ensuring participation.

#### Recommendation
**Confidence:** ★★☆☆☆ (moderate uncertainty)
The current value of 0.012 d⁻¹ appears reasonable based on first principles and the need to balance participation with synchrony. However, captive breeding observations at Friday Harbor Laboratories may provide more precise estimates.

### 2. p_spontaneous_male (Daily Spontaneous Spawning Probability - Males)

**Current Value:** 0.0125 d⁻¹  
**Range:** 0.008–0.020 d⁻¹

#### First Principles
Males can spawn multiple times per season (unlike females), so their base rate should be similar to or slightly higher than females to ensure adequate sperm availability. The model allows males 2–3 bouts per season, requiring higher overall spawning probability across the season.

#### Literature Evidence

**Sexual Differences:** Sea stars show no sexual dimorphism (Animal Diversity Web), and both sexes participate in broadcast spawning simultaneously. However, energetic costs differ dramatically—sperm is metabolically cheap compared to eggs.

**Multi-bout Capacity:** Unlike females who release their entire egg mass, males can potentially spawn repeatedly if metabolically supported. The literature does not provide specific data on male spawning frequency in *Pycnopodia*.

#### Recommendation  
**Confidence:** ★★☆☆☆ (moderate uncertainty)
Current value slightly higher than females (0.0125 vs 0.012) reflects potential for multiple male spawning events. This requires field validation from captive breeding programs.

### 3. peak_width_days (Seasonal Peak Standard Deviation)

**Current Value:** 60.0 days
**Range:** 30–90 days  

#### First Principles
Standard deviation of the Gaussian seasonal readiness curve. At 60 days, 95% of spawning occurs within a ~4-month window. Narrower peaks increase fertilization success but raise extinction risk from mistimed seasons; wider peaks provide bet-hedging but reduce fertilization efficiency.

#### Literature Evidence

**Seasonal Timing:** Animal Diversity Web reports *P. helianthoides* spawning "between March and July" (5 months) with "main peak in May and June" (2 months). This suggests a relatively concentrated breeding season.

**Comparative Context:** 
- *Odontaster validus* (Antarctic sea star) reproduces "once a year during the winter season, between the months of April and June, with peak spawning occurring during June" (Animal Diversity Web) — suggesting 2-3 month breeding windows are typical for cold-water asteroids
- Schiebelhut et al. (2022) found phylogenetic signals in reproductive seasons, indicating evolutionary constraint on timing

**Environmental Drivers:** The literature emphasizes that spawning synchrony is crucial for broadcast fertilizers. Environmental cues likely coordinate timing, but *Pycnopodia*-specific triggers are poorly documented.

#### Recommendation
**Confidence:** ★★★☆☆ (moderate-high confidence)
Current value of 60 days is well-supported by the March-July season with May-June peak reported in Animal Diversity Web. This represents a 4-month effective season (±2σ) with 2-month peak window.

### 4. female_max_bouts (Maximum Female Spawning Bouts)

**Current Value:** 2 bouts per season
**Range:** 1–3 bouts per season

#### First Principles
Each spawning event costs substantial energy (gonad development = 10-30% of body mass). Most asteroids are thought to spawn once per season due to energetic constraints, but larger species like *Pycnopodia* may have capacity for multiple smaller releases.

#### Literature Evidence

**Energetic Constraints:** Astropecten literature (Helgoland Marine Research, 2016) notes that "resources stored in pyloric cecum seem to play an important role in the seasonal production of gonads in some asteroid species," suggesting tight energetic trade-offs.

**Size Advantage:** As the largest known sea star species (up to 5 kg, 80 cm diameter), *Pycnopodia* may have greater energetic capacity for multiple spawning events compared to smaller asteroids.

**Lack of Species Data:** No direct data on *Pycnopodia* spawning frequency. Most asteroid literature assumes single annual spawning, but this may reflect smaller species or limited observation periods.

**Model Context:** The original spawning specification suggested female_max_bouts = 1, but current config shows 2. This discrepancy needs resolution.

#### Recommendation
**Confidence:** ★☆☆☆☆ (low confidence)
Conservative estimate of 1-2 bouts per season based on energetic constraints, with the possibility that large individuals can support multiple smaller releases. Priority candidate for empirical validation through captive breeding programs.

---

## Key Knowledge Gaps

1. **No quantitative spawning data** for *Pycnopodia helianthoides* — all parameters based on comparative literature and first principles
2. **Captive breeding data not yet analyzed** — Friday Harbor Laboratory observations (Hodin et al. 2021) may contain relevant timing data
3. **Sex-specific differences unclear** — literature suggests no dimorphism, but energetic costs clearly differ
4. **Environmental triggers unknown** — what cues coordinate spawning synchrony?
5. **Multiple spawning evidence lacking** — do females actually spawn multiple times?

## Research Priorities

1. **Analyze FHL captive breeding data** for quantitative spawning frequencies and timing
2. **Field observations during breeding season** (March-July) to document natural spawning patterns  
3. **Energetic analysis** of gonad development cycles to constrain maximum bout frequencies
4. **Environmental correlation analysis** with temperature, photoperiod, and food availability

## Model Implications

- **Disease transmission:** Post-spawning immunosuppression creates vulnerability windows that depend on spawning frequency and timing
- **Allee effects:** Low populations may fail to achieve spawning synchrony, creating additional extinction risk beyond fertilization limitations
- **Evolutionary dynamics:** Spawning timing affects which individuals survive to reproduce, potentially driving selection on reproductive strategies
- **Conservation:** Captive breeding success depends on understanding and mimicking natural spawning cues

---

**Sources:** 103 papers in local literature database, Animal Diversity Web, comparative asteroid literature via web search. Priority: obtain FHL captive breeding dataset for quantitative calibration.