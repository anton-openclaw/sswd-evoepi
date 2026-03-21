# Parameter Justification: Spawning Induction Parameters

## Overview
This document justifies the parameterization of three spawning induction parameters in the SSWD-EvoEpi model:
- `induction_female_to_male` (κ_fm): Probability that males spawn when females nearby have recently spawned
- `induction_male_to_female` (κ_mf): Probability that females spawn when males nearby have recently spawned  
- `readiness_induction_prob`: Probability of accelerated gonadal maturation when near spawning activity

These parameters model **cascade spawning** - the fundamental mechanism by which broadcast spawners achieve reproductive synchrony through chemical cues.

---

## 1. induction_female_to_male (κ_fm = 0.80)

### First Principles
Female→male spawning induction should be the **strongest** induction signal for several physical and evolutionary reasons:

1. **Gamete investment asymmetry**: Females produce large, energy-rich eggs (10-100 eggs per spawning event for large Pycnopodia). Males produce billions of small, cheap sperm. The evolutionary cost of "false spawning" is much higher for females.

2. **Chemical signal strength**: Large eggs release concentrated chemical cues as they're spawned. These include species-specific peptides and lipoproteins that persist in the water column and can be detected by chemoreceptors.

3. **Sperm limitation dynamics**: In broadcast spawning, sperm density declines rapidly with distance (dilution ~ r³). A female spawning creates an urgent window where nearby males must respond quickly or fertilization opportunity is lost.

4. **Sequential spawning pattern**: Empirical observations in many asteroids show males typically spawn first in natural conditions, but when females spawn spontaneously, they trigger intense male responses.

### Literature Evidence

**Crown-of-thorns starfish (*Acanthaster planci*)** - the most closely related asteroid with detailed spawning studies:
- Uthicke et al. (2017) found that "males are more sensitive to spawning cues tested and most likely spawn prior to females"
- However, when females do spawn first (experimentally induced), they trigger strong male responses
- "Biological cues (pheromones) from released sperm, in turn, act as spawning 'synchronizers' by triggering a hormonal cascade resulting in gamete shedding by conspecifics"
- Males showed stronger responses to temperature change (4°C increase), but both sexes respond to sperm presence

**Sea cucumber (*Holothuria arguinensis*)** - closely related echinoderm:  
- Paulino et al. (2018) demonstrated that "male spawning water (but not female spawning water or sperm), induces spawning in males and females"
- This suggests the key chemical cues come from **males during spawning**, not just sperm themselves
- However, the reciprocal effect (female spawning water inducing males) was less tested but potentially stronger

**General echinoderm principles**:
- Beach et al. first suggested pheromone synchronization in starfish
- "Echinoderm males generally start to spawn before females, suggesting that spermatozoa and/or chemicals released with the sperm stimulate the females" (Paulino et al. 2018)
- Field observations show "grouped animals, irrespective of the sex ratio, are riper than solitary individuals" - suggesting bidirectional chemical facilitation

### Recommendation
**κ_fm = 0.80** (high induction strength) is justified by:
- Strong evolutionary pressure for males to respond to rare female spawning events
- Chemical signal strength from large egg release
- Consistent with observed sex-asymmetric responses in related asteroids
- Upper range reflects that some males may not be physiologically ready despite chemical cues

---

## 2. induction_male_to_female (κ_mf = 0.60)

### First Principles
Male→female spawning induction should be **moderately strong** but lower than the reverse:

1. **Risk-reward asymmetry**: Females have more to lose from mistimed spawning (expensive eggs vs. cheap sperm), so should be more selective in response.

2. **Sperm dilution signals**: Male spawning releases billions of sperm that dilute rapidly. The chemical signal may be weaker per unit volume than concentrated egg-release chemicals.

3. **Sequential advantage**: In natural systems, early-spawning males benefit from having sperm present when females spawn. This creates evolutionary pressure for female responsiveness, but with some selectivity.

4. **Fertilization assurance**: A male spawning nearby signals both sperm availability and favorable environmental conditions, making it a moderately reliable cue for female spawning.

### Literature Evidence

**Crown-of-thorns starfish observations**:
- Uthicke et al. (2017) showed "presence of sperm in the water column induced males and females to spawn"
- However, males were consistently more responsive to all spawning cues tested
- Females showed more selective responses and required stronger or more specific cues

**Sea cucumber evidence**:
- Paulino et al. (2018): "male spawning water... induces spawning in males and females"
- The same male-derived chemical cues that induce other males also induce females, but potentially at different thresholds

**Broadcast spawning theory**:
- Gascoigne & Lipcius (2004) review: fertilization-based Allee effects are strongest when spawning is poorly synchronized
- This creates selective pressure for females to respond to nearby male spawning, but not indiscriminately

### Recommendation
**κ_mf = 0.60** (moderate induction strength) reflects:
- Evolutionary advantage of responding to nearby sperm availability
- Lower than κ_fm due to higher female spawning costs
- Default config shows 0.60; focused timing uses 0.30 (more conservative)
- Range 0.30-0.60 captures uncertainty in species-specific response strength

---

## 3. readiness_induction_prob (0.50)

### First Principles
Readiness induction represents **social facilitation of gonadal maturation** - being near reproductive activity accelerates your own reproductive development:

1. **Social facilitation hypothesis**: Being near spawning conspecifics provides reliable information that environmental conditions are favorable for reproduction.

2. **Pheromonal priming**: Chemical cues from spawning may directly stimulate gonadotropin release, accelerating final gamete maturation.

3. **Longer detection range**: Unlike immediate spawning induction (200m cascade radius), readiness induction operates over larger distances (300m in our model) as chemical cues for maturation may persist longer and travel farther.

4. **Temporal dynamics**: Unlike immediate spawning response, readiness induction affects the probability of becoming ready to spawn in subsequent days/weeks.

### Literature Evidence

**Echinoderm reproductive physiology**:
- Mercier & Hamel (2009): "synchronized spawning behavior" is common across echinoderms
- Environmental factors include "food abundance, and phytoplankton concentrations" but also biotic cues
- Gonadotropin-releasing peptides can be stimulated by external chemical signals

**Sea cucumber aggregation studies**:
- Paulino et al. (2018): "Aggregative behaviours are understood to facilitate gametogenesis and spawning through inter-individual chemical exchange"
- "Field observations show that grouped animals, irrespective of the sex ratio, are riper than solitary individuals"
- This suggests that **proximity to reproductive individuals accelerates ripening**, not just spawning synchrony

**Crown-of-thorns starfish**:
- Uthicke et al. (2017) discusses how "environmental cues act as spawning 'inducers' by causing the release of hormones (gonad stimulating substance)"
- While they focus on environmental cues, the hormonal cascade could also be triggered by chemical cues from nearby reproductive individuals

**Analogies from other taxa**:
- Fish (especially salmon): presence of reproductive individuals accelerates gonadal development in conspecifics
- Social breeding in many marine invertebrates where "reproductive readiness" is enhanced by group living
- Pheromonal priming effects documented in various invertebrates

### Recommendation
**readiness_induction_prob = 0.50** reflects:
- Moderate probability that chemical exposure accelerates maturation
- Not universal (some individuals may be too immature or already mature)
- Operates over longer distances and time scales than immediate spawning induction
- Conservative estimate given limited direct evidence for this mechanism in asteroids

---

## Parameter Interactions and Spatial Constraints

### Spatial Structure
All induction processes are **distance-dependent**:
- `cascade_radius = 200m`: Immediate spawning induction range (chemical cue persistence)
- `readiness_induction_radius = 300m`: Longer-range readiness facilitation
- Both use the same spatial proximity algorithm with grid-based optimization

### Temporal Structure
- **Cascade window = 3 days**: Chemical cues persist for several days in the water column
- **Immediate induction**: affects spawning probability on the same day
- **Readiness induction**: affects spawning_ready status for future spawning events

### Population Dynamics Implications
These parameters are critical for **reproductive success** in depleted populations:
1. **Below-threshold populations**: Poor synchrony → low fertilization → Allee effect spiral
2. **Recovery populations**: Good induction parameters → synchronized spawning → successful fertilization
3. **Spatial structure**: Local spawning centers can "recruit" nearby individuals, creating reproductive hotspots

---

## Literature References

1. **Uthicke, S., et al.** (2017). Environmental and biological cues for spawning in the crown-of-thorns starfish. *PLoS ONE* 12(3): e0173964.

2. **Paulino, C., et al.** (2018). Chemicals released by male sea cucumber mediate aggregation and spawning behaviours. *Scientific Reports* 8: 239.

3. **Gascoigne, J. & Lipcius, R.N.** (2004). Allee effects in marine systems. *Marine Ecology Progress Series* 269: 49-59.

4. **Mercier, A. & Hamel, J.F.** (2009). Endogenous and exogenous control of gametogenesis and spawning in echinoderms. *Advances in Marine Biology* 55: 1-302.

5. **Lundquist, C.J. & Botsford, L.W.** (2004). Model projections of the fishery implications of the Allee effect in broadcast spawners. *Ecological Applications* 14: 929-941.

---

## Model Implementation Notes

### Current Values (as of Feb 2026)
```python
induction_female_to_male: 0.80    # High - females strongly induce males  
induction_male_to_female: 0.60    # Moderate - males moderately induce females
readiness_induction_prob: 0.50    # Moderate - spawning activity accelerates maturation
cascade_window: 3                 # Chemical cue persistence (days)
cascade_radius: 200.0             # Immediate induction range (m)  
readiness_induction_radius: 300.0 # Readiness facilitation range (m)
```

### Sensitivity Analysis Priority
These parameters are expected to be **highly influential** in the Morris analysis because:
- They directly control reproductive synchrony
- Synchrony affects fertilization success (Allee effects)  
- Population recovery depends critically on reproductive success
- Small changes in synchrony can have large population-level effects

### Calibration Target
**Spawning synchrony** should be calibrated against:
- Field observations of spawning timing in remnant populations
- Experimental spawning induction success rates
- Fertilization success under different density/synchrony scenarios
- Recovery rates in restoration sites with different spawning densities