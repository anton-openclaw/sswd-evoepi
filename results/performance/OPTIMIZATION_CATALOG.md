# SSWD-EvoEpi Optimization Opportunity Catalog

**Generated:** February 16, 2026  
**Author:** Anton 🔬  
**Context:** Systematic code review for optimization opportunities that don't change model behavior  
**Status:** Complete analysis of all major modules

## Executive Summary

**Current Performance:** 50-agent, 2-year simulation = 24.93 seconds  
**Primary Bottleneck:** Random number generation (80% of runtime)  
**Secondary Bottlenecks:** Python for-loops throughout the codebase  
**Total Identified Opportunities:** 23 optimizations across 8 modules  
**Expected Combined Speedup:** 30-50× improvement possible

## 🔥 CRITICAL PRIORITY - Massive Speedup Opportunities

### Priority 1: Random Number Generation Bottleneck - **reproduction.py**
- **File**: `sswd_evoepi/reproduction.py`, lines 665-675 (in `_generate_larval_cohort()`)
- **Issue**: 13.5M individual random calls (80% of total runtime = 19.92s)
- **Current**: Nested loops calling `rng.integers(0, 2)` for each offspring × each locus
- **Fix**: Batch random generation with single array call
- **Speedup**: ~15-20×
- **Difficulty**: Medium
- **Risk**: None (same RNG sequence with proper reshaping)

**Implementation:**
```python
# BEFORE: 13.5M individual calls
for offspring in range(n_offspring):
    for locus in range(N_LOCI):
        maternal_allele = rng.integers(0, 2)
        paternal_allele = rng.integers(0, 2)

# AFTER: Single batch call
all_choices = rng.integers(0, 2, size=(n_offspring, N_LOCI, 2))
maternal_choices = all_choices[:, :, 0]
paternal_choices = all_choices[:, :, 1]
```

### Priority 2: Mendelian Inheritance Vectorization - **reproduction.py**
- **File**: `sswd_evoepi/reproduction.py`, lines 350-358
- **Issue**: `for l in range(N_LOCI):` loop for allele copying (part of 20% bottleneck)
- **Current**: Per-locus copying in Python loop
- **Fix**: Vectorize with fancy indexing
- **Speedup**: ~5-8×
- **Difficulty**: Medium
- **Risk**: None

## 🔴 HIGH PRIORITY - Major Performance Gains

### Priority 3: Disease Transmission Loops - **disease.py**
- **File**: `sswd_evoepi/disease.py`, lines 602, 618
- **Issue**: Two sequential for-loops processing infected agents
- **Current**: `for idx in susc_indices[infected_mask]:` and `for idx in diseased_indices:`
- **Fix**: Vectorize with boolean indexing and batch operations
- **Speedup**: ~8-12×
- **Difficulty**: Medium
- **Risk**: Low

### Priority 4: Movement Gravity Loop - **movement.py**
- **File**: `sswd_evoepi/movement.py`, lines 215-256
- **Issue**: `for i in gravity_indices_local:` processing spawning aggregation
- **Current**: Per-agent distance calculations and heading updates
- **Fix**: Vectorized distance matrix + broadcasting for center-of-mass
- **Speedup**: ~10-15×
- **Difficulty**: Hard
- **Risk**: Low

### Priority 5: Settler Placement Loop - **reproduction.py**
- **File**: `sswd_evoepi/reproduction.py`, lines 762-790
- **Issue**: `for j in range(n_slots):` individual agent initialization
- **Current**: Per-settler field assignment
- **Fix**: Batch array assignment with slicing
- **Speedup**: ~5-10×
- **Difficulty**: Easy
- **Risk**: None

### Priority 6: Mutation Application Loop - **genetics.py**
- **File**: `sswd_evoepi/genetics.py`, lines 332-338
- **Issue**: `for idx in mut_flat_idx:` applying individual mutations
- **Current**: Per-mutation bit flip
- **Fix**: Vectorized bit flipping with advanced indexing
- **Speedup**: ~3-5×
- **Difficulty**: Medium
- **Risk**: None

## ⚡ MEDIUM PRIORITY - Solid Performance Improvements

### Priority 7: Ne Calculation Loops - **reproduction.py**
- **File**: `sswd_evoepi/reproduction.py`, lines 438-443
- **Issue**: Two separate loops counting offspring per parent
- **Current**: `for i, f in enumerate(females):` and similar for males
- **Fix**: Use `np.bincount()` or `np.unique()` with return counts
- **Speedup**: ~3-5×
- **Difficulty**: Easy
- **Risk**: None

### Priority 8: Parent ID Processing - **reproduction.py**
- **File**: `sswd_evoepi/reproduction.py`, lines 407-410
- **Issue**: `for pid in all_parent_ids:` counting offspring per parent
- **Current**: Manual counting in Python loop
- **Fix**: `np.bincount()` for direct counting
- **Speedup**: ~2-4×
- **Difficulty**: Easy
- **Risk**: None

### Priority 9: Genotype Bank Processing - **genetics.py**
- **File**: `sswd_evoepi/genetics.py`, lines 788-795
- **Issue**: `for l_idx in range(N_LOCI):` processing loci individually
- **Current**: Per-locus genotype bank operations
- **Fix**: Vectorize across all loci simultaneously
- **Speedup**: ~2-3×
- **Difficulty**: Medium
- **Risk**: None

### Priority 10: Movement Exposure Factor Nested Loops - **movement.py**
- **File**: `sswd_evoepi/movement.py`, lines 404-409
- **Issue**: Nested `for di in range(3): for dj in range(3):` loops
- **Current**: 3×3 grid neighborhood calculation
- **Fix**: Pre-compute offsets array, vectorize with broadcasting
- **Speedup**: ~2-3×
- **Difficulty**: Easy
- **Risk**: None

### Priority 11: Genetic Summary Statistics Loops - **genetics.py**
- **File**: `sswd_evoepi/genetics.py`, lines 607, 616
- **Issue**: Two separate `for l_idx in range(N_LOCI):` loops
- **Current**: Per-locus statistics computation
- **Fix**: Vectorized operations across all loci
- **Speedup**: ~2-3×
- **Difficulty**: Medium
- **Risk**: None

## 💡 OPTIMIZATION OPPORTUNITIES - Algorithmic Improvements

### Priority 12: Distance Matrix Nested Loops - **spatial.py**
- **File**: `sswd_evoepi/spatial.py`, lines 100-103
- **Issue**: O(N²) nested loops for Haversine calculations
- **Current**: `for i in range(N): for j in range(i+1, N):`
- **Fix**: Vectorized Haversine using broadcasting
- **Speedup**: ~5-15×
- **Difficulty**: Medium
- **Risk**: None

### Priority 13: Overwater Distance Matching Loops - **spatial.py**
- **File**: `sswd_evoepi/spatial.py`, lines 154-177
- **Issue**: Multiple sequential loops for coordinate matching
- **Current**: Per-node nearest neighbor search
- **Fix**: Vectorized distance calculation with `scipy.spatial.distance`
- **Speedup**: ~3-8×
- **Difficulty**: Medium
- **Risk**: None

### Priority 14: SST Generation Nested Loops - **environment.py**
- **File**: `sswd_evoepi/environment.py`, lines 101-106
- **Issue**: `for yr_idx in range(n_years): for d in range(365):`
- **Current**: Per-day SST calculation
- **Fix**: Vectorized day array generation
- **Speedup**: ~5-10×
- **Difficulty**: Easy
- **Risk**: None

## 🔧 LOW PRIORITY - Minor Optimizations

### Priority 15-23: Various Model.py Loops
- **File**: `sswd_evoepi/model.py`, lines 221, 295, 318, 451, 653, 657
- **Issue**: Multiple small loops in demographic updates
- **Current**: Per-agent processing for mortality, growth, etc.
- **Fix**: Batch array operations
- **Speedup**: ~1.2-2× each
- **Difficulty**: Easy-Medium
- **Risk**: None

## 🚀 ADVANCED OPTIMIZATION OPPORTUNITIES

### GPU Acceleration Candidates
1. **Distance matrix calculations** (Priority 12, 13) - Perfect for CUDA
2. **Genetic operations** (Priority 2, 6, 9, 11) - Embarrassingly parallel
3. **Disease transmission** (Priority 3) - Spatial operations ideal for GPU

### Numba JIT Compilation Candidates
1. **Movement gravity calculations** (Priority 4) - Hot inner loop
2. **Disease progression** (Priority 3) - Complex state transitions
3. **Reproduction SRS lottery** - Mathematical operations

### Memory Optimization Opportunities
1. **Pre-allocated working arrays** - Reduce allocation overhead
2. **In-place operations** - Minimize memory copies
3. **Structured arrays** - Better cache locality

## 📊 Implementation Roadmap

### Phase 1: Critical Bottlenecks (Expected 20-30× speedup, 2-4 hours)
1. **Priority 1**: Random number generation vectorization
2. **Priority 2**: Mendelian inheritance vectorization
3. **Priority 5**: Settler placement vectorization

**Target**: 24.93s → 1-2s (25× improvement)

### Phase 2: Major Loops (Expected 2-5× additional speedup, 4-6 hours)
1. **Priority 3**: Disease transmission vectorization
2. **Priority 4**: Movement gravity vectorization
3. **Priority 6**: Mutation application vectorization

**Target**: 1-2s → 0.3-0.8s (3× improvement)

### Phase 3: Algorithmic Improvements (Expected 1.5-3× additional speedup, 2-4 hours)
1. **Priority 12**: Distance matrix vectorization
2. **Priority 14**: SST generation vectorization
3. **Priority 7-11**: Medium priority vectorizations

**Target**: 0.3-0.8s → 0.1-0.4s (2× improvement)

### Phase 4: Advanced Optimizations (Expected 2-10× additional speedup, 1-2 weeks)
1. **GPU acceleration** for distance/genetic operations
2. **Numba JIT** for hot loops
3. **Memory optimization** patterns

**Target**: 0.1-0.4s → 0.01-0.1s (5× improvement)

## 🎯 Expected Final Performance

**Current**: 50-agent, 2-year simulation = 24.93s  
**After Phase 1**: ~1s (25× improvement)  
**After Phase 2**: ~0.4s (60× improvement)  
**After Phase 3**: ~0.2s (125× improvement)  
**After Phase 4**: ~0.05s (500× improvement)

**Scaling Implications:**
- **500 agents, 20 years**: Currently ~8 hours → ~2 minutes after optimization
- **5,000 agents, 50 years**: Currently impossible → ~30 minutes after optimization
- **150-node coastline network**: Currently impractical → feasible real-time simulation

## ⚠️ Risk Assessment

**No Risk (Priorities 1, 2, 5, 7, 8, 10, 14)**: Pure computational optimizations
**Low Risk (Priorities 3, 4, 6, 9, 11-13)**: Require careful testing but straightforward
**Medium Risk (Advanced optimizations)**: GPU/Numba may introduce platform dependencies

## 🧪 Validation Strategy

1. **Unit tests**: Every optimization must pass existing test suite
2. **Checksum comparison**: Identical outputs with fixed RNG seeds
3. **Edge case testing**: Boundary conditions, empty arrays, large populations
4. **Performance regression tests**: Ensure no future performance degradation
5. **Cross-platform testing**: Verify optimizations work on different systems

---

**Conclusion**: The SSWD-EvoEpi codebase contains numerous high-impact optimization opportunities. Implementing the top 6 priorities alone could achieve 50-100× speedup, making large-scale simulations computationally feasible for the first time.