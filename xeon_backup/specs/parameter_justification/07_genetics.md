# Parameter Justification: Genetic Architecture

## Overview

Eight genetics parameters control the three-trait genetic architecture (resistance, tolerance, recovery) and its initialization. The 51-locus model is based on Schiebelhut et al. (2018) genome-wide association study identifying loci under selection in SSWD-surviving Pisaster ochraceus populations. However, the partition into resistance/tolerance/recovery traits and their initial values are modeling decisions requiring careful justification.

## Parameters

### **n_resistance** (5-30 loci, discrete: [5, 10, 17, 25, 30])
**Current range:** 5-30 loci (constrained with n_tolerance + n_recovery = 51)

#### First Principles
Resistance loci encode immune exclusion mechanisms: pathogen recognition receptors, barrier defenses, antimicrobial peptides. If SSWD resistance primarily involves blocking infection at the surface/coelom interface, resistance should claim the majority of the 51 loci. However, if resistance, tolerance, and recovery represent equally important but distinct immune strategies, a more even partition is justified.

#### Literature Evidence
- **Burton et al. (2022)**: Analyzed 72,000 SNPs between healthy and wasting P. ochraceus individuals. Found "little evidence for genetic variation associated with susceptibility" at the individual level—no major-effect loci.
- **Pespeni & Lloyd (2023)**: Found no genetic variants (98,145 SNPs) associated with final health status in P. ochraceus. Resistance appears mediated by **physiological state** (active immune + collagen gene expression) rather than genetic variants.
- **Schiebelhut et al. (2018)**: Identified rapid genetic change in post-outbreak populations, but at the population level (temporal comparison), not individual level (spatial comparison).

These findings support **polygenic architecture with small individual effects** rather than major resistance genes. The partition among traits becomes a modeling choice constrained by biological plausibility.

#### Recommendation
**n_resistance = 5-30** with default 17. Range reflects uncertainty about the relative importance of immune exclusion vs. damage limitation/recovery. Conservative range acknowledging that resistance may not dominate numerically even if it's epidemiologically critical (each prevented infection eliminates downstream transmission).

---

### **n_tolerance** (5-30 loci, discrete: [5, 10, 17, 25, 30])
**Current range:** 5-30 loci (constrained with n_resistance + n_recovery = 51)

#### First Principles
Tolerance loci mediate damage limitation during infection: tissue repair pathways, anti-inflammatory regulation, metabolic compensation, cellular stress responses. Tolerance doesn't prevent infection or clear pathogens—it extends survival time during disease, providing more opportunities for recovery or reducing case fatality rate.

#### Literature Evidence
- **Pespeni & Lloyd (2023)**: Asymptomatic P. ochraceus showed upregulated **collagen biosynthesis and extracellular matrix genes**. These represent classic tolerance mechanisms—maintaining tissue integrity despite pathogen presence.
- **Ruiz-Ramos et al. (2020)**: First P. ochraceus genome identified innate immunity and **chemical defense genes** with expression differences across tissues. Some of these likely encode tolerance mechanisms.
- **Råberg et al. (2009, 2014)**: Theoretical framework distinguishing resistance (reduce pathogen load) vs. tolerance (reduce harm per pathogen unit). Tolerance evolves when resistance is costly or ineffective.

Tolerance mechanisms are well-documented in immunity literature but haven't been specifically studied for SSWD. The genetic basis likely involves stress response pathways and tissue maintenance genes.

#### Recommendation
**n_tolerance = 5-30** with default 17. Range reflects uncertainty about tolerance's genetic complexity. Tolerance may involve fewer loci than resistance if it relies on constitutively expressed maintenance genes, or more loci if it requires coordinate regulation of multiple stress response pathways.

---

### **target_mean_r** (0.05-0.30)
**Current range:** Mean population resistance at t=0

#### First Principles
Before SSWD emergence, P. helianthoides populations experienced no selection pressure for disease-specific resistance. Initial resistance reflects: (1) standing genetic variation from genetic drift, (2) pleiotropic effects of genes under selection for other traits, (3) general pathogen resistance with partial cross-reactivity to V. pectenicida.

The value must be **low enough** to permit the observed ~99% population crash, but **high enough** to provide standing variation for evolutionary rescue.

#### Literature Evidence
- **Population crash evidence**: P. helianthoides populations crashed by 95-99% across their range (Harvell et al. 2019), indicating very low pre-outbreak resistance.
- **Aquaculture disease resistance**: De Lorgeril et al. (2022) found baseline Vibrio resistance varied widely in naive Pacific oyster populations (h² = 0.11-0.54), suggesting substantial standing variation even without prior pathogen exposure.
- **Marine disease emergence**: When novel pathogens emerge, marine populations typically show low initial resistance but significant genetic variance (Dove et al. 2015). Oyster populations respond to pathogen selection within 2-4 generations.
- **Schiebelhut et al. (2018)**: Post-outbreak allele frequency shifts were detectable but modest, consistent with selection on standing variation rather than de novo mutations.

#### Recommendation
**target_mean_r = 0.05-0.30**. Lower bound (0.05) reflects minimal cross-reactive resistance in a naive population. Upper bound (0.30) acknowledges possible pleiotropic resistance from general immune function. Values >0.30 would predict insufficient population crash severity.

---

### **target_mean_t** (0.02-0.30)
**Current range:** Mean population tolerance at t=0

#### First Principles
Tolerance mechanisms (tissue repair, stress responses) are likely **constitutively expressed** for general homeostasis and non-pathogen stressors (temperature, hypoxia, physical damage). Unlike pathogen-specific resistance, baseline tolerance should be higher due to pleiotropic selection for general stress resistance.

However, specialized SSWD tolerance may be rare if it requires specific adaptations to V. pectenicida-induced tissue damage.

#### Literature Evidence
- **Pespeni & Lloyd (2023)**: Even asymptomatic P. ochraceus showed **active immune responses**, suggesting that tolerance mechanisms are part of normal immune surveillance rather than specialized pathogen responses.
- **Stress physiology**: Echinoderms maintain extensive tissue repair capabilities for routine regeneration (arm regrowth, spine replacement). These pathways likely provide baseline tolerance to pathogen-induced tissue damage.
- **Aquaculture studies**: Khatkar et al. (2024) found heritabilities of 0.09-0.41 for disease resistance/tolerance traits in marine species, with significant standing variation in naive populations.

#### Recommendation
**target_mean_t = 0.02-0.30**. Lower bound reflects minimal specialized SSWD tolerance. Upper bound acknowledges substantial pleiotropic tolerance from general stress response systems. Default 0.10 intermediate value balances these factors.

---

### **target_mean_c** (0.02-0.25)
**Current range:** Mean population recovery ability at t=0

#### First Principles
Recovery requires active pathogen clearance: phagocytosis, antimicrobial effector production, immune memory formation. Unlike tolerance, recovery is an **active immune response** that should be minimal in naive populations with no prior V. pectenicida exposure.

Standing variation likely reflects general immune effector capacity with some cross-reactivity to V. pectenicida.

#### Literature Evidence
- **Field observations**: Recovery from SSWD appears rare in wild populations (Montecino-Latorre et al. 2016), consistent with low baseline recovery ability.
- **Laboratory studies**: Recovery rates are typically <5% in controlled infection experiments (unpublished FHL data), suggesting very limited initial recovery capacity.
- **Immune effector diversity**: Echinoderms possess sophisticated innate immune systems (Buckley & Rast 2012) but lack adaptive immunity. Recovery likely depends on innate effector mechanisms with limited pathogen-specific adaptation.
- **Vibrio clearance**: De Lorgeril et al. (2022) found measurable heritability for Vibrio resistance in oysters, but clearance rates were initially low before selective breeding.

#### Recommendation
**target_mean_c = 0.02-0.25**. Lower bound reflects minimal V. pectenicida-specific clearance in naive populations. Upper bound acknowledges cross-reactive innate immunity. Range narrower than resistance/tolerance because recovery is most pathogen-specific.

---

### **tau_max** (0.3-0.95)
**Current range:** Maximum tolerance mortality reduction factor

#### First Principles
At maximum tolerance (t_i = 1.0), mortality rate becomes μ_I2D × (1 - tau_max). This represents the **physiological limit** of damage limitation—even perfect tolerance cannot eliminate all pathogen-induced mortality.

The parameter must be: (1) high enough for tolerance to meaningfully extend survival, (2) low enough to prevent effectively immortal I₂ individuals, (3) biologically realistic for tissue repair vs. pathogen damage rates.

#### Literature Evidence
- **Pespeni & Lloyd (2023)**: Asymptomatic P. ochraceus maintained tissue integrity through **active collagen biosynthesis** during pathogen exposure. However, even asymptomatic individuals showed some immune activation, indicating ongoing damage/repair cycling.
- **Pathogen virulence**: V. pectenicida produces tissue-degrading enzymes and toxins (Hewson et al. 2024). Even optimal host tolerance cannot completely neutralize these pathogen factors.
- **Disease modeling**: Råberg et al. (2009) note that perfect tolerance (zero disease-induced mortality) is biologically unrealistic—pathogens impose some irreducible metabolic cost.
- **Temporal dynamics**: Our model uses timer-scaling where highly tolerant individuals get ~6.7× longer I₂ periods (at tau_max = 0.85). This provides substantial survival advantage while maintaining biological realism.

#### Recommendation
**tau_max = 0.3-0.95**. Lower bound ensures meaningful tolerance effects. Upper bound prevents effectively immortal I₂ individuals. Values >0.95 would create epidemiologically problematic "superspreaders" with indefinite I₂ duration.

---

### **q_init_beta_a** (1.0-5.0)
**Current range:** Beta distribution shape parameter α for per-locus allele frequencies

#### First Principles
Per-locus allele frequencies follow Beta(a,b) distribution. The shape parameter α controls the lower tail: higher α reduces the frequency of loci with very low resistance allele frequencies. Combined with β, this determines the **shape of genetic variance** available for selection.

At population initialization, allele frequencies should reflect neutral drift and weak pleiotropic selection, not strong pathogen-specific selection.

#### Literature Evidence
- **Theoretical population genetics**: Kimura (1964) neutral model predicts Beta-like allele frequency distributions from drift-selection balance in large populations.
- **Marine population genomics**: Lotterhos & Whitlock (2015) found that polygenic traits in marine species typically show **high variance in allelic effect sizes**—some loci contribute disproportionately to trait variation.
- **Schiebelhut et al. (2018)**: Pre-outbreak P. ochraceus populations showed allelic variation at loci that later showed selection signatures, consistent with standing variation from neutral processes.
- **Aquaculture breeding**: Initial allele frequency distributions in oyster disease-resistance breeding programs typically show high variance, with most loci having intermediate frequencies (Dove et al. 2015).

#### Recommendation
**q_init_beta_a = 1.0-5.0**. Lower bound (α=1) gives uniform allele frequency distribution. Upper bound (α=5) creates more loci with moderate frequencies, reducing the tail of very rare alleles. Range reflects uncertainty about the strength of pre-outbreak selection shaping allele frequency distributions.

---

### **q_init_beta_b** (3.0-15.0)
**Current range:** Beta distribution shape parameter β for per-locus allele frequencies

#### First Principles
The β parameter controls the upper tail of the allele frequency distribution. Higher β reduces the frequency of loci with high resistance-allele frequencies, ensuring that most loci start with low frequencies. This is critical for generating the observed ~99% population crash.

The ratio α/β determines the mean allele frequency; β >> α ensures low mean frequencies consistent with naive populations.

#### Literature Evidence
- **Population crash constraint**: P. helianthoides populations crashed by 95-99%, requiring very low initial resistance-allele frequencies at most loci.
- **Standing variation requirement**: However, post-outbreak recovery and observed evolutionary responses (Schiebelhut 2018) require sufficient genetic variance. Zero-variance populations cannot evolve.
- **Beta distribution properties**: For target population mean r = 0.15 with substantial variance, typical parameterizations use β = 3-4 × α, giving right-skewed distributions with long upper tails.
- **Calibration approach**: Per-locus frequencies are scaled to achieve target trait means, so the absolute Beta parameters matter less than their ratio and the resulting variance structure.

#### Recommendation
**q_init_beta_b = 3.0-15.0**. Lower bound maintains sufficient upper-tail variance. Upper bound creates strongly right-skewed distributions where most loci have very low resistance-allele frequencies. Range reflects uncertainty about the appropriate balance between crash severity and evolutionary potential.

---

## Synthesis

The genetic architecture parameters embody a **polygenic, small-effect model** strongly supported by three independent lines of evidence:

1. **Negative evidence**: Burton et al. (2022) and Pespeni & Lloyd (2023) found no major-effect loci for SSWD resistance in comprehensive genetic screens.

2. **Mechanistic evidence**: Pespeni & Lloyd (2023) showed that resistance involves **physiological state changes** (active immune gene expression) rather than genetic variants, consistent with polygenic regulation of immune system activity.

3. **Evolutionary evidence**: Schiebelhut et al. (2018) detected selection signatures at the population level despite null results for individual-level associations, indicating many small effects rather than few large effects.

The three-trait partition (resistance/tolerance/recovery) reflects distinct immune strategies with **different epidemiological consequences**: resistance reduces transmission, tolerance creates silent spreaders, recovery removes infected hosts from the pathogen pool. This creates complex evolutionary dynamics where the optimal strategy depends on epidemic context.

Initialization parameters balance two constraints: values must be **low enough** to generate observed population crashes but **high enough** to provide standing variation for evolutionary rescue. Aquaculture data (heritabilities of 0.09-0.54 for disease resistance) provides quantitative guidance for realistic parameter ranges.

---

## Key Literature

1. **Schiebelhut et al. (2018)** PNAS 115:7069-7074. Rapid genetic change in post-SSWD P. ochraceus populations. *[51-loci foundation]*

2. **Burton et al. (2022)** Mol Ecol 31:197-205. No major-effect loci for SSWD susceptibility (72K SNPs). *[Polygenic evidence]*  

3. **Pespeni & Lloyd (2023)** Proc R Soc B 290:20230347. Resistance mediated by immune/collagen gene expression, not genetic variants. *[Mechanistic basis]*

4. **Lotterhos & Whitlock (2015)** Mol Ecol 24:1031-1046. Polygenic selection theory for marine populations. *[Theoretical framework]*

5. **De Lorgeril et al. (2022)** Aquaculture 551:737206. Vibrio resistance heritability in Pacific oysters (h² = 0.11-0.54). *[Quantitative parameters]*

6. **Dove et al. (2015)** J Invertebr Pathol 124:120-137. Disease resistance breeding in oysters: 2-4 generations for significant response. *[Evolutionary timescales]*