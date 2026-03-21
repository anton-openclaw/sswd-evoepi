# Reintroduction Experiment: CA Release in 2025

## Overview
Test captive-bred sunflower star reintroduction into California under varying
genetic resistance levels and release densities. All scenarios branch from
W330 baseline at year 13 (2025), run forward to year 38 (2050).

## Experimental Design

### Axis 1: Resistance Genetics (5 levels)

| ID | Label | Mean resistance (r̄) | Description |
|----|-------|---------------------|-------------|
| R1 | CA-2025 survivors | ~0.31 | Breed from wild CA survivors. Highly selected. |
| R2 | Pre-epidemic naive | ~0.15 | Pre-2013 broodstock. No disease selection. |
| R3 | Mid-epidemic 2019 | ~0.20 | Moderate selection (~6 years evolution) |
| R4 | Artificially selected | ~0.50 | Aggressive breeding program |
| R5 | Immune | 1.0 | Theoretical ceiling — best possible outcome |

### Axis 2: Release Density (4 levels)

| ID | Stars per site | Notes |
|----|---------------|-------|
| D1 | 50 | Small pilot release |
| D2 | 200 | Moderate release |
| D3 | 500 | Large release (at K=1000, this is 50% of K) |
| D4 | 1000 | Full saturation (= K) |

### Full Factorial: 5 × 4 = 20 scenarios

Naming: `REINTRO_{R#}_{D#}` (e.g., `REINTRO_R1_D2` = CA survivors at 200/site)

### Release Parameters
- **Release year**: 13 (= 2025)
- **Release regions**: CA-S, CA-C, CA-N (all California sites)
- **Release genetics**:
  - Tolerance (t̄): 0.10 (pre-epidemic default for all)
  - Recovery (c̄): 0.05 (pre-epidemic default for all)
  - Resistance (r̄): varies per scenario
  - Genetic variance: match initial Va for resistance
- **Baseline run**: W330 (best calibrated), K=1000, K_ref=5000

### Key Questions
1. **Is naive reintroduction viable?** Do pre-epidemic genetics survive in
   disease-endemic warm waters?
2. **Resistance threshold**: What r̄ is needed for a self-sustaining population?
3. **Density threshold**: Does release size matter, or do Allee effects dominate
   at low density?
4. **Interaction**: Is there a genetics × density interaction? (e.g., immune
   animals succeed at any density, but naive animals need critical mass)
5. **P_env barrier**: Does the warm-water pathogen reservoir overwhelm even
   resistant animals in CA? (P_env gating may prevent recovery regardless)
6. **Spillover**: Do released CA populations seed recovery in OR/WA via
   larval transport?
7. **Evolutionary dynamics**: Do released genetics persist or get swamped by
   local survivor genetics? Does Va recover?

### Response Variables
- Population trajectory (CA-S, CA-C, CA-N) 2025-2050
- Recovery fraction at 2050 vs pre-epidemic peak
- Time to self-sustaining population (if ever)
- Mean resistance at 2050 (did selection erode released genetics?)
- Spillover: population change in OR, WA-O, JDF
- Disease prevalence in release regions

### Implementation Notes
- Uses release module in `sswd_evoepi/model.py` (`process_release_event`)
- `ReleaseEvent` config: time_step, node_id, n_individuals, genetics_mode,
  trait_targets (dict with 'resistance', 'tolerance', 'recovery')
- 238 CA nodes: CA-S (78), CA-C (97), CA-N (63) — node IDs 224-885
- Release at year 13 = day 13*365 = 4745 (2025)
- Each scenario needs 238 ReleaseEvent entries (one per CA node)
- genetics_mode='trait_targets' is simplest: specify mean r, use defaults for t,c
- At K=1000, D4 (1000/site) doubles the local population instantly
- Immune (r=1.0) means zero force-of-infection gets through — tests whether
  disease-free growth alone is sufficient for recovery
- Release events go in config JSON under "release_events": [...]
- Each event: {"time_step": 4745, "node_id": N, "n_individuals": D,
  "genetics_mode": "trait_targets", "trait_targets": {"resistance": R}}

### Comparison to Literature
- **Arroyo-Esquivel et al. 2026**: 3-class matrix model, no spatial structure,
  no genetics, no environmental pathogen. Our IBM adds all of these.
- **Hodin et al. 2021**: Captive breeding protocol, spawning Nov-Mar. Our
  release module can match this timing.
- Our experiment directly addresses: "Should we breed for resistance, or just
  breed for numbers?"
