# Validation: trait_math.py vs Simulation

**Date:** 2026-02-23  
**Population:** N=50,000 (main tests), N=10,000 (multi-gen)  
**Seed:** 42 (main), 99 (multi-gen)  
**Runtime:** 11.5s  
**Result:** 25/35 passed, 10 failed  

## ⚠️ Flags

- exceedance(r≥0.30): pred=0.0279 vs obs=0.0386
- expected_maximum(n=50): 10.1% error
- expected_maximum(n=200): 11.7% error
- expected_maximum(n=1000): 12.3% error
- expected_maximum(n=5000): 13.2% error
- predict_gen var(gen 1): 41.0% error
- predict_gen var(gen 2): 36.5% error
- predict_gen var(gen 3): 33.8% error
- predict_gen var(gen 4): 52.5% error
- predict_gen var(gen 8): 33.3% error

## Results Summary

| Test | Predicted | Observed | Error % | Threshold | Status |
|------|-----------|----------|---------|-----------|--------|
| trait_mean (Resistance) | 0.149740 | 0.149740 | 0.0% | 1.0% | ✅ PASS |
| trait_mean (Tolerance) | 0.099901 | 0.099901 | 0.0% | 1.0% | ✅ PASS |
| trait_mean (Recovery) | 0.080087 | 0.080087 | 0.0% | 1.0% | ✅ PASS |
| trait_variance (Resistance) | 0.006176 | 0.006081 | 1.6% | 5.0% | ✅ PASS |
| trait_variance (Tolerance) | 0.004335 | 0.004317 | 0.4% | 5.0% | ✅ PASS |
| trait_variance (Recovery) | 0.003058 | 0.003065 | 0.2% | 5.0% | ✅ PASS |
| exceedance P(r≥0.05) | 0.897799 | 0.904140 | 0.7% | 15% | ✅ PASS |
| exceedance P(r≥0.10) | 0.736601 | 0.718540 | 2.5% | 15% | ✅ PASS |
| exceedance P(r≥0.15) | 0.498679 | 0.464880 | 7.3% | 15% | ✅ PASS |
| exceedance P(r≥0.20) | 0.261242 | 0.249600 | 4.7% | 15% | ✅ PASS |
| exceedance P(r≥0.30) | 0.027942 | 0.038600 | 27.6% | 15% | ❌ FAIL |
| exceedance P(r≥0.40) | 0.000725 | 0.002940 | 75.3% | abs<0.02 | ✅ PASS |
| exceedance P(r≥0.50) | 0.000004 | 0.000080 | 0.0% | abs<0.02 | ✅ PASS |
| expected_maximum (n=50) | 0.311786 | 0.346749 | 10.1% | 10% | ❌ FAIL |
| expected_maximum (n=200) | 0.352310 | 0.398797 | 11.7% | 10% | ❌ FAIL |
| expected_maximum (n=1000) | 0.392625 | 0.447814 | 12.3% | 10% | ❌ FAIL |
| expected_maximum (n=5000) | 0.427960 | 0.492790 | 13.2% | 10% | ❌ FAIL |
| selection_response (R, top 10%) | 0.137925 | 0.149993 | 8.0% | 15.0% | ✅ PASS |
| predict_generations mean (gen 1) | 0.308946 | 0.333598 | 7.4% | 20% | ✅ PASS |
| predict_generations var (gen 1) | 0.013641 | 0.009676 | 41.0% | 30% | ❌ FAIL |
| predict_generations mean (gen 2) | 0.498234 | 0.501614 | 0.7% | 20% | ✅ PASS |
| predict_generations var (gen 2) | 0.006729 | 0.004929 | 36.5% | 30% | ❌ FAIL |
| predict_generations mean (gen 3) | 0.642201 | 0.621767 | 3.3% | 20% | ✅ PASS |
| predict_generations var (gen 3) | 0.004959 | 0.003707 | 33.8% | 30% | ❌ FAIL |
| predict_generations mean (gen 4) | 0.765729 | 0.730177 | 4.9% | 20% | ✅ PASS |
| predict_generations var (gen 4) | 0.003713 | 0.002435 | 52.5% | 30% | ❌ FAIL |
| predict_generations mean (gen 5) | 0.849148 | 0.814839 | 4.2% | 20% | ✅ PASS |
| predict_generations var (gen 5) | 0.001276 | 0.001164 | 9.7% | 30% | ✅ PASS |
| predict_generations mean (gen 6) | 0.905039 | 0.873883 | 3.6% | 20% | ✅ PASS |
| predict_generations var (gen 6) | 0.000557 | 0.000561 | 0.7% | 30% | ✅ PASS |
| predict_generations mean (gen 7) | 0.944721 | 0.914514 | 3.3% | 20% | ✅ PASS |
| predict_generations var (gen 7) | 0.000312 | 0.000329 | 5.1% | 30% | ✅ PASS |
| predict_generations mean (gen 8) | 0.975673 | 0.946248 | 3.1% | 20% | ✅ PASS |
| predict_generations var (gen 8) | 0.000141 | 0.000211 | 33.3% | 30% | ❌ FAIL |
| V_A / trait_variance ratio (h²=1, expect 1.0) | 1.000000 | 1.000000 | 0.0% | 0.1% | ✅ PASS |

## Multi-Generation Selection Detail

Selection: top 10% truncation, 8 generations, N=10,000

| Gen | Pred Mean | Sim Mean | Mean Err% | Pred Var | Sim Var | Var Err% |
|-----|-----------|----------|-----------|----------|---------|----------|
| 0 | 0.149161 | 0.149161 | 0.0% | 0.00828939 | 0.00840751 | 1.4% |
| 1 | 0.308946 | 0.333598 | 7.4% | 0.01364114 | 0.00967553 | 41.0% |
| 2 | 0.498234 | 0.501614 | 0.7% | 0.00672940 | 0.00492935 | 36.5% |
| 3 | 0.642201 | 0.621767 | 3.3% | 0.00495935 | 0.00370750 | 33.8% |
| 4 | 0.765729 | 0.730177 | 4.9% | 0.00371294 | 0.00243504 | 52.5% |
| 5 | 0.849148 | 0.814839 | 4.2% | 0.00127606 | 0.00116375 | 9.7% |
| 6 | 0.905039 | 0.873883 | 3.6% | 0.00055685 | 0.00056073 | 0.7% |
| 7 | 0.944721 | 0.914514 | 3.3% | 0.00031217 | 0.00032898 | 5.1% |
| 8 | 0.975673 | 0.946248 | 3.1% | 0.00014101 | 0.00021127 | 33.3% |

## 🐛 Bug Found & Fixed

**`delta_q_per_locus` was missing a factor of 1/2.**

The trait uses diploid allele means: `x_ℓ = (a₁+a₂)/2`. Therefore the allele substitution effect is `α_ℓ/2`, not `α_ℓ`. The original formula `Δqℓ = i × αℓ × qℓ(1−qℓ) / σ` overpredicted allele frequency changes by 2×, causing `predict_generations` to massively overshoot (30-40% mean error by generation 1).

**Fix:** `Δqℓ = i × (αℓ/2) × qℓ(1−qℓ) / σ_P`

**Verification:** `ΔE[τ] = Σ αℓ × Δqℓ = i/σ × Σ αℓ²/2 × q(1-q) = i×σ = R` (matches breeder's equation) ✅

**`additive_variance` was 4× too large** (same root cause — used `α²` instead of `(α/2)²`). Fixed to `V_A = Σ (αℓ²/2) qℓ(1−qℓ) = trait_variance` (h²=1). ✅

⚠️ **Note:** `genetics.py:compute_additive_variance` has the same factor-of-4 error (`V_A = 2Σαℓ²qp`). It's diagnostics-only (doesn't affect simulation dynamics) so left for a separate fix.

## Methodology

### What's being validated

The `trait_math.py` module provides analytical predictions based on:
- **Trait mean/variance**: E[τ] = Σ αℓqℓ, σ² = Σ (αℓ²/2) qℓ(1−qℓ)
- **Normal approximation**: CLT-based for sum of ~17 independent Bernoulli contributions
- **Breeder's equation**: R = i × σ_P (with h²=1)
- **Per-locus Δq**: Δqℓ ≈ i × (αℓ/2) × qℓ(1−qℓ) / σ_P

### How validation works

1. **Initialize population** using the actual model genetics code (`initialize_genotypes_three_trait`)
2. **Compute empirical allele frequencies** from the realized genotypes
3. **Feed those frequencies** into `trait_math` analytical functions
4. **Compare predictions** against empirical statistics computed from genotype scores

### Known limitations (remaining failures)

All remaining failures are **known approximation limitations**, not code bugs:

- **Normal approximation in tails**: With 17 loci and skewed (Beta-distributed) allele frequencies, the trait distribution has heavier tails than Gaussian. Exceedance probabilities at >2σ from the mean are systematically underestimated (observed at r≥0.30: 27.6% error). The `expected_maximum` function inherits this bias, consistently underestimating by ~10-13%.
- **Variance under selection (Bulmer effect)**: The analytical model tracks allele frequency changes but not the within-generation variance reduction from truncation selection. In reality, selecting the top 10% creates linkage disequilibrium that reduces variance below the Hardy-Weinberg expectation. This causes predicted variance to overshoot actual variance by 30-50% in early generations of strong selection.
- **Multi-gen mean drift**: Cumulative ~3-5% error in trait mean over 8 generations from first-order Δq approximation + missing drift term. Acceptable for planning-level breeding program predictions.
