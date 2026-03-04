# Wavefront Long-Range Dispersal Kernel — Spec

## Problem

W42 (CDT=1000) successfully bridges BC→AK but the wavefront is 3.5× too slow:
- W42 AK-FS arrival: month 91 (target: month 26)
- Observed: CA-S June 2013 → AK-WG/AL Jan 2017 = ~42 months, ~5000 km ≈ 120 km/month
- Model: ~55 km/month

## Root Cause

The wavefront speed is limited by a sequential cascade:
1. Node activates → local epidemic takes weeks to ramp up
2. Only after ramp-up does the node shed enough pathogen for meaningful dispersal
3. The D matrix (D_P=50km, max_range=500km) decays exponentially with distance
4. Distant unreached nodes accumulate dose very slowly from each newly-activated node

The wave can only advance as fast as local epidemics develop — each node must build up
pathogen levels before it can push the front forward. This is a Fisher-KPP wave speed
limitation: c = 2√(rD), where both r (epidemic growth) and D (dispersal range) matter.

## Solution: Long-Range Wavefront Dispersal Kernel

**Biological justification**: Local pathogen dynamics are driven by nearby shedding
(D_P ~50km). But long-range transport via oceanic currents, tidal flushing, and 
possibly biotic vectors can carry viable pathogen over much larger distances at LOW
concentrations. These low concentrations aren't enough to sustain local epidemics 
directly, but over time they accumulate and eventually trigger an outbreak (cumulative
dose mechanism). This is "stratified dispersal" — well-established in invasion ecology.

**Implementation**: A SECOND dispersal kernel with longer range, used ONLY for 
cumulative dose accumulation at unreached nodes. The standard D matrix continues to 
govern local pathogen dynamics at reached nodes. This preserves all existing disease 
behavior while allowing realistic wavefront speeds.

## New Config Parameters

In `DiseaseSection`:
```python
wavefront_D_P: float = 0.0              # Wavefront dispersal scale (km). 0 = use standard D_P
wavefront_D_P_max_range: float = 0.0    # Wavefront max range (km). 0 = use standard D_P_max_range
```

When `wavefront_D_P > 0`:
- A second sparse D matrix is constructed with the longer-range kernel
- This matrix is used ONLY for computing dispersal input at UNREACHED nodes
- Reached nodes continue using the standard D matrix for all dynamics

When `wavefront_D_P == 0` (default):
- Behavior is identical to current code (backward compatible)

## Code Changes

### 1. config.py — Add parameters to DiseaseSection
- `wavefront_D_P: float = 0.0`
- `wavefront_D_P_max_range: float = 0.0`

### 2. model.py — Build and use wavefront D matrix
- At simulation init (near line 2132): if `wavefront_D_P > 0`, construct second D matrix
  using `construct_pathogen_dispersal()` with the wavefront parameters, then build sparse D^T
- In daily loop (near line 2276): for unreached node dose accumulation, use wavefront 
  dispersal instead of standard dispersal
- Both sparse matrices co-exist; one extra matmul per day (cheap — sparse)

### 3. tests/test_wavefront.py — Add tests
- Backward compat: wavefront_D_P=0 gives identical results
- Correctness: wavefront_D_P > D_P accumulates dose faster at distant unreached nodes
- Config validation tests for new parameters

## Chunk Plan

### Chunk 1: Config + second D matrix construction
- Add `wavefront_D_P` and `wavefront_D_P_max_range` to DiseaseSection
- In model.py init section: construct second D matrix when wavefront_D_P > 0
- Build sparse D^T for wavefront kernel

### Chunk 2: Model integration + tests
- In daily loop: use wavefront dispersal for unreached node dose accumulation
- Write all tests (backward compat, correctness, config validation)
- Run full test suite

### Chunk 3: Kill W45-W52, generate new configs, launch
- Kill running processes on Xeon
- Generate W53+ configs with wavefront_D_P sweeps
- Push code, launch new calibration batch

## Expected Behavior

With `wavefront_D_P=200` and `wavefront_D_P_max_range=2000`:
- Unreached nodes 1000km ahead of the wavefront front receive dose from distant shedding nodes
- Multiple established epidemic nodes contribute simultaneously to distant dose
- Wavefront speed increases because distant nodes "see" the epidemic building up far behind the front
- Local disease dynamics at reached nodes are UNCHANGED (still use D_P=50km kernel)
