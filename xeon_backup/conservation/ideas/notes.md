# Conservation Genetics — Future Ideas

Parking lot for ideas that are beyond the current scope but worth
exploring in future work.

---

## Inbreeding Depression Modeling

The current model tracks F and Hₑ but doesn't model fitness
consequences of inbreeding. Future work could add:

- **Lethal equivalents**: w̄(F) = exp(-B × F), B ∈ [2, 12] for marine inverts
- **Purging during breeding**: allow selection against deleterious recessives
- **Environment-dependent depression**: inbreeding × disease stress interaction
- **Empirical calibration**: Need captive breeding data on Pycnopodia
  survival × F to estimate B

Status: Eq. 6.7 in the theory report. Straightforward to add once
we have empirical estimates of B for Pycnopodia or close relatives.

---

## Genomic Selection (GBLUP)

Replace the 51-locus additive model with genomic prediction:

- **GBLUP**: Genomic Best Linear Unbiased Prediction using the full
  GRM instead of trait-specific loci
- **Training population**: Use model output as "training data" to
  build prediction equations
- **Accuracy**: Compare GBLUP accuracy to the 51-locus index
- **Integration**: Could connect to actual genotyping data if
  Pycnopodia get a SNP chip or low-coverage WGS panel

Status: Architecturally clean — the GRM code in `inbreeding.py`
already implements VanRaden Method 1. Need to add mixed model
equations (Henderson's). Could use existing packages (e.g., PyBLUP).

---

## Climate Change Projections

Current analyses use satellite SST climatology (static).
Future work could incorporate:

- **SSP scenarios**: CMIP6 projections for NE Pacific SST
- **Temperature × disease interaction**: Warmer water = faster disease
  progression = stronger selection
- **Range shifts**: Pycnopodia habitat may shift northward
- **Phenological mismatch**: Spawning timing × SST changes
- **Impact on breeding targets**: How does the "optimal" resistance
  level change under warming?

Status: The model already has SST-dependent disease rates (Arrhenius
scaling). Plugging in SSP time series is straightforward. The hard
part is the uncertainty — SSP ensemble variance is large.

---

## Multi-Species Interactions

The current model treats the pathogen reservoir as an abstraction
(P_env). Future work could model:

- **Multi-host dynamics**: Other asteroid species as subclinical
  carriers or alternative hosts
- **Community recovery**: Does Pycnopodia recovery depend on
  community-level Vibrio dynamics?
- **Ecological feedbacks**: Pycnopodia as keystone predator →
  urchin populations → kelp forest recovery
- **Co-evolutionary dynamics**: Host-pathogen arms race across
  multiple species

Status: Architecturally explored (see MEMORY.md parked ideas).
Multi-species model is a 2+ year project. Not for this paper,
but a natural PhD chapter for someone.

---

## Epigenetic Effects

Transgenerational epigenetic inheritance could affect:

- **Disease priming**: Parental exposure → offspring resistance
  (seen in some marine inverts)
- **Captive effects**: Epigenetic changes in captivity vs wild
- **Stress memory**: Temperature stress → heritable methylation changes

Status: No empirical data for Pycnopodia. Speculative but worth
noting as a potential confound for breeding program predictions.
Would require adding a non-Mendelian inheritance layer.

---

## Cost Optimization

Add economic constraints to the analysis:

- **Sampling costs**: Travel + collection + genotyping per individual
  (varies by site accessibility)
- **Breeding costs**: Facility, husbandry, per-generation costs
- **Release costs**: Transport, monitoring, per-individual costs
- **Budget-constrained optimization**: Maximize conservation outcome
  per dollar instead of per individual
- **Multi-year budgeting**: Optimal spending trajectory over
  program lifetime

Status: Requires cost data from aquaculture/conservation practitioners.
The optimization framework (screening allocation, OCS) extends
naturally to cost-weighted objectives.

---

## Genetic Rescue from Related Species

Could alleles from closely related Asteroidea provide resistance?

- **Candidate donors**: *Pisaster* spp., *Dermasterias*, other
  non-susceptible or less-susceptible species
- **Hybridization risk**: Reproductive barriers, hybrid fitness
- **Targeted introgression**: CRISPR or marker-assisted backcrossing

Status: Very speculative. Pycnopodia-specific adaptations may not
transfer. Regulatory and ethical barriers are significant.
Not for current modeling work.

---

## Spatial Optimization of Release Sites

More sophisticated release location selection:

- **Network centrality**: Release at nodes with highest downstream
  connectivity (larval dispersal)
- **Source-sink dynamics**: Avoid releasing at sink populations
- **Stepping-stone optimization**: Minimize number of release sites
  for range-wide allele spread
- **Stochastic optimization**: Account for variable oceanographic
  conditions in dispersal

Status: The 11-node stepping-stone network and distance matrix
are already implemented. Could use graph-theoretic optimization
(betweenness centrality, etc.) once reintroduction analysis is
running.

---

*Last updated: 2026-02-23*
*Add new ideas freely — this is a living document.*
