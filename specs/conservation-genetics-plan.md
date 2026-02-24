# Conservation Genetics Analysis Plan

**Status:** Framework laid, awaiting calibrated model  
**Dependencies:** Sobol R4 → ABC-SMC calibration → parameter lock → this  
**Prototype scripts:** `scripts/trait_distributions.py`, `scripts/breeding_analysis.py`  
**Prototype results:** `results/trait_distributions/`, `results/breeding_analysis/`

---

## Goal

Use the calibrated SSWD-EvoEpi model to generate actionable conservation genetics predictions:
1. What does the genetic landscape of surviving Pycnopodia look like **now** (2026)?
2. How does it vary by **site** and **latitude**?
3. How will it change over the **next 5–20 years** with and without intervention?
4. What is the optimal **breeding program design** (founders, crosses, generations)?
5. Where should we **source founders** vs **release captive-bred stock**?

---

## Phase 1: Post-Epidemic Genetic State (Site × Year)

### What we'll do
Run the calibrated model from pre-epidemic baseline (~2013) → 2026 across all 11 stepping-stone nodes. Extract per-site:

- **Trait distributions** (r, t, c) of surviving individuals
- **Allele frequencies** at each of the 51 loci
- **Population sizes** and genetic diversity (heterozygosity, effective Ne)
- **Selection intensity** experienced (cumulative mortality, selective differential)

### Key outputs
- **Latitude × trait heatmap**: Mean resistance (and other traits) by site and year
- **Site-specific screening tables**: "At Site X, expect 1 in N individuals with r ≥ threshold"
- **Genetic diversity gradient**: Which sites retain the most standing variation?
- **Extinction risk by site**: Which populations are below minimum viable Ne?

### Replication
- Multiple seeds (10–50) to characterize stochastic variation
- Confidence intervals on all predictions
- Ensemble mean + worst-case + best-case scenarios

### What changes from prototypes
- Use calibrated parameters (not defaults)
- Run to 2026 specifically (not generic 20-year)
- Use actual satellite SST time series (not climatology)
- Track full genotype data at endpoints (not just trait means)
- All 11 nodes (not just 5 validation nodes)

---

## Phase 2: Screening Effort Analysis (Site-Specific)

### What we'll do
For each site's predicted 2026 population, compute:

- **Trait distributions** (as fraction of local population)
- **Exceedance curves** with 1-in-N ratios (per site)
- **Screening effort** tables: How many to sample for target resistance level
- **Joint trait probabilities**: Finding individuals good at multiple traits
- **Breeding value distributions**: Combined survival-per-exposure

### Key question
> "If I go to Friday Harbor / Sitka / Howe Sound today and collect 100 stars, what resistance levels will I find?"

### Latitude gradient
Southern sites (warmer → stronger disease → stronger selection):
- Fewer survivors but higher mean resistance
- Less genetic diversity (bottleneck)
- Better individual founders, worse population-level diversity

Northern sites (cooler → weaker disease → weaker selection):
- More survivors but lower mean resistance  
- More genetic diversity preserved
- Better source for diversity, worse for pre-selected resistance

**Conservation implication:** Source founders from BOTH ends of the gradient. Northern sites provide diversity; southern survivors provide pre-selected resistance alleles.

---

## Phase 3: Breeding Program Optimization

### What we'll do
Simulate captive breeding programs with realistic constraints:

#### Founder Selection
- **Sample sizes**: 50, 100, 200, 500 wild-caught from various sites
- **Multi-site sourcing**: Optimal allocation across sites (how many from each?)
- **Complementarity screening**: Find founders with non-overlapping protective loci
- **Genetic diversity vs resistance trade-off**: Maximizing both

#### Crossing Strategies
- **Random mating** within selected pool (baseline)
- **Assortative mating** by resistance (pair high × high)
- **Complementary mating** (pair individuals with different protective loci)
- **Minimum kinship mating** (maximize diversity while selecting for resistance)

#### Selection Schemes
- **Truncation selection**: Keep top K% each generation
- **Family selection**: Keep best from each family (maintains diversity)
- **Optimal contribution selection**: Maximize genetic gain while constraining inbreeding (ΔF)
- **Index selection**: Weight resistance, tolerance, recovery differently

#### Generation parameters  
- **Generation time**: Pycnopodia reach maturity in ~2 years (literature estimate)
- **Fecundity**: Can produce millions of larvae, so family sizes aren't limiting
- **Realistic program sizes**: 50–500 breeding adults per generation
- **Inbreeding constraint**: Track Ne, flag when ΔF > threshold

### Key outputs
- **Generations to target**: How many generations to reach r̄ = 0.30 / 0.50 / 0.70?
- **Founder requirement**: Minimum founders for a viable program
- **Inbreeding trajectory**: When does genetic diversity become limiting?
- **Optimal strategy**: Which crossing scheme gives best gain per generation?
- **Cost-benefit**: Marginal gain per additional generation of breeding

### What changes from prototypes
- Track inbreeding (F coefficient) — prototype ignores it
- Realistic population sizes (not 10K founders)
- Multiple crossing strategies compared head-to-head
- Family structure tracked (not just mass selection)
- Generation overlap (Pycnopodia are iteroparous)

---

## Phase 4: Reintroduction Scenario Modeling

### What we'll do
Use the full spatial model to simulate reintroduction of captive-bred stock:

#### Scenarios
1. **No intervention** (baseline): Let natural evolution + connectivity play out
2. **Single release**: X captive-bred individuals at one site, one time
3. **Repeated supplementation**: Periodic releases over multiple years
4. **Multi-site release**: Distribute across sites simultaneously
5. **Stepping-stone release**: Sequential releases following larval connectivity
6. **Genetic rescue**: Release diverse stock into bottlenecked populations

#### Variables to test
- **Release size**: 100, 500, 1,000, 5,000 per event
- **Release timing**: Relative to seasonal disease peak
- **Release location**: Which node(s)?
- **Genetic composition**: Mean resistance of released stock
- **Frequency**: One-shot vs annual vs every 2 years

#### Success metrics
- Population persistence (10, 20, 50 year horizons)
- Mean resistance trajectory (coast-wide and per-site)
- Genetic diversity (He, allelic richness)
- Population recovery to % of historical carrying capacity
- Disease-free periods (consecutive years without major outbreak)

### This is Willem's domain
Willem has a separate conservation scenario plan. This section should be refined in collaboration with him once the model is calibrated.

---

## Phase 5: Paper Figures & Tables

### Conservation genetics figures for publication
1. **Latitude × year heatmap** of predicted resistance (the "money figure")
2. **Site-specific screening effort curves** (practical guidance)
3. **Breeding program trajectories** (generations to resistance targets)
4. **Founder source optimization** (which sites, how many)
5. **Reintroduction scenario comparison** (what works?)
6. **Genetic diversity vs resistance trade-off** (can we have both?)

### Tables
1. Per-site predicted genetic state in 2026
2. Screening recommendations by site
3. Breeding program design parameters
4. Reintroduction scenario outcomes

---

## Technical Requirements

### Model extensions needed
- [ ] **Genotype snapshots**: Save full genotype arrays at specified time points (currently only trait means recorded)
- [ ] **Inbreeding tracking**: F coefficient calculation from pedigree or genotype homozygosity
- [ ] **Captive population module**: Separate breeding population with controlled mating
- [ ] **Release mechanism**: Introduce captive-bred individuals into wild nodes at specified times
- [ ] **Genetic diversity metrics**: He, allelic richness, Ne estimation in recorder

### Computational requirements
- Phase 1 (genetic state): ~50 seeds × 11 nodes × 13 years ≈ moderate (hours on Xeon)
- Phase 3 (breeding optimization): Fast (seconds per run, many scenarios)
- Phase 4 (reintroduction scenarios): ~100+ scenarios × 50 seeds × 20+ years ≈ heavy (days on Xeon)

### Data needs
- Satellite SST time series 2013–2026 (extend current climatology to actual annual data)
- Actual SSWD outbreak timing by region (for model initialization)
- Captive breeding generation time estimates (literature)
- Current population size estimates by region (if available)

---

## Timeline

1. **Now**: Prototype scripts done, this plan written ✅
2. **After Sobol R4** (~Feb 26): Identify which parameters matter most for conservation outcomes
3. **After ABC-SMC calibration** (~March?): Lock parameters, begin Phase 1
4. **Phase 1–2**: ~1 week computation + analysis
5. **Phase 3**: ~1 week (computationally light, analytically heavy)
6. **Phase 4**: Coordinate with Willem's conservation plan
7. **Phase 5**: Integrate into paper

---

## Notes

- The prototype results (8 generations to r̄ = 1.0) are optimistic because they assume perfect selection and no inbreeding depression. Real programs will be slower.
- Generation time is the binding constraint. Pycnopodia mature in ~2 years, so 8 generations = ~16 years of breeding. That's a long program.
- The model currently doesn't include inbreeding depression, which could be significant in small captive populations. This is a known gap.
- Complementary crossing (pairing individuals with different protective loci) is dramatically more effective than random mating among resistant individuals. This should be a key recommendation.
- The latitude gradient creates a natural tension: source from where survivors are most resistant (south) or where diversity is highest (north)? The answer is probably both.
