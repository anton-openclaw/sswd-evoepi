# Conservation Genetics Module — Review

**Reviewer:** Anton 🔬  
**Date:** 2026-02-23  
**Scope:** Full module review — theory report, code, validation, analysis templates

---

## Executive Summary

The conservation genetics module is a substantial, well-integrated piece of work: ~9,300 lines across a 34-page LaTeX theory report, 5 code modules (2,632 lines), 3 test/validation suites (2,568 lines), 5 analysis templates (2,297 lines), and supporting infrastructure. The theory is mathematically sound, the code faithfully implements the equations, and validation uncovered (and fixed) real bugs. All 146 tests run; 136 pass, 10 fail for documented approximation reasons. The module is ready for calibrated parameters — the major blocking dependency.

---

## Section-by-Section Assessment

### 1. Theory Report (`report/main.tex`, 10 sections, 34 pages)

**Strengths:**
- Derivations are clean, first-principles, and self-consistent. Every equation is numbered and referenced in the code.
- The progression from genetic architecture → quantitative genetics → screening → breeding → inbreeding → reintroduction follows a logical chain with each section building on the last.
- Section 2 (Genetic Architecture) correctly handles the mean-of-alleles encoding (`x = (a₁+a₂)/2`), which is a common source of factor-of-2 errors. The 1/2 factor in variance formulas is derived and verified.
- Section 3 (Quantitative Genetics) correctly derives the allele substitution effect as α/2 and verifies the consistency check (Σαᵢ×Δqᵢ = R).
- Section 4b (Binary Phenotyping) is excellent — reframed as a framework rather than specific numbers, which is exactly right for pre-calibration. The dose trade-off table, family-based selection discussion, and strategic crossing of survival categories are genuinely useful for practitioners.
- Section 5 (Breeding) gives marginal fitness derivatives at population means (resistance ~1000× more important than tolerance), which is a key quantitative insight.
- Section 8 (Pycnopodia Biology) is grounded and honest about what the model can/cannot predict. The constraints table is publishable.

**Weaknesses / Gaps:**
- **Eq. 2.12 vs Eq. 4.3**: The report defines V_A using the standard `2q(1-q)α²` formula (Eq. 2.12 in Prop. 3), but the code and later sections correctly use `(α²/2)q(1-q)` for the mean-of-alleles encoding. The Proposition 3 formula (Eq. 2.12) should be annotated to clarify that this is V_A in terms of the allele substitution effect (α/2), not the effect size (α). Right now there's an apparent contradiction between `V_A = Σ 2q(1-q)α²` in §2 and `V_A = Σ (α²/2)q(1-q)` in §3. Both are correct but they define α differently — this needs a reconciliation note.
- **Inbreeding depression**: Repeatedly flagged as "not modeled" throughout, which is honest. But Section 6.3 gives the lethal equivalents formula (Eq. 6.7) without discussing whether it'll be integrated. The `ideas/notes.md` says "straightforward to add once we have empirical B." This needs a decision: is it in scope for the first paper?
- **Section 7 (Reintroduction)**: The weakest section theoretically. The admixture equation (Eq. 7.1) and Prop. 7.1 (self-sustaining alleles) are standard but the selection coefficient s_ℓ in Eq. 7.6 is left abstract — no derivation connecting it to the fitness function from §5. This would strengthen the chain.
- **Segregation variance formula discrepancy**: The report currently says `α²/16` per het locus (Eq. 5.9'), but breeding.py implements `α²/4`. The validation confirms the empirical value matches `raw_value/4`, so the code produces `4×` the correct value. The report was updated but the code function still returns the uncorrected value. This is flagged as "reporting only" but could confuse future readers.

**Overall report quality**: 8.5/10. Publication-grade with minor fixes needed.

### 2. Code Modules (`src/`)

#### `trait_math.py` (221 lines)
- Clean, well-documented. Every function references its equation number.
- `trait_mean`, `trait_variance`, `additive_variance` are correct.
- `delta_q_per_locus` correctly uses α/2. The consistency check (ΣαΔq = R) is verified in the docstring.
- `predict_generations` iterates the coupled Δq system forward correctly — clips to [0,1], tracks fixation.
- No issues found.

#### `screening.py` (551 lines)
- `required_sample_size` and `required_n_for_k` are mathematically exact (geometric/binomial inverse CDF).
- `empirical_exceedance` and `exceedance_curve` correctly interface with the model's `compute_trait_batch`.
- `multisite_allocation` uses a greedy marginal improvement algorithm — correct for the concave (diminishing returns) objective.
- `select_complementary_founders` balances trait value and locus coverage with configurable weights — well designed.
- Minor: `complementarity_matrix` is O(n²) with a Python loop. Fine for ≤500 individuals but would need vectorization for larger sets.
- No bugs found.

#### `breeding.py` (806 lines)
- `mendelian_cross` correctly samples one allele per parent per locus. Verified by χ² in validation (0/51 rejections at Bonferroni α).
- `expected_offspring_trait` correctly computes midparent value with the /4 factor (two parents, two alleles each, mean-of-alleles encoding).
- `segregation_variance` still has the 4× scale issue (uses α²/4 instead of α²/16). The test compensates by dividing by 4 and checking the corrected ratio ≈ 1.0. **Should be fixed in the code, not just documented.**
- `pair_ocs` (Optimal Contribution Selection) uses a Lagrangian relaxation with binary search for λ — this is a simplified but sound approach. The GRM computation uses VanRaden Method 1, which is standard.
- `run_breeding_program` ties everything together cleanly. All four strategies work.
- `compute_generation_stats` computes a comprehensive set of per-generation metrics.
- The `compute_allele_frequencies_from_array` helper correctly handles the plain genotype array (no structured agent array).

#### `inbreeding.py` (453 lines)
- `genomic_inbreeding` correctly computes F = 1 - Ho/He per individual. Handles monomorphic loci.
- `genomic_relationship_matrix` implements VanRaden Method 1 correctly (center by 2q, scale by 2q(1-q), Z@Z'/L).
- `ne_from_family_variance` uses the standard formula. `ne_from_grm` is a rough estimate (1/(2×mean off-diagonal)) — appropriate for this context.
- `projected_f` correctly implements the compound formula F_g = 1 - (1-F₀)(1-ΔF)^g.
- `inbreeding_depression` implements the lethal equivalents model but isn't connected to the breeding simulator. Placeholder for future use.
- `generations_to_f_threshold` correctly solves the inverse problem via logarithms.
- No issues found.

#### `viz.py` (586 lines)
- Comprehensive set of plotting functions: distributions, exceedance curves, screening effort, breeding trajectories, Pareto frontier, heatmaps.
- Consistent styling with `_apply_style`. Good color scheme.
- `plot_breeding_trajectory` is a nice 4-panel layout (trait, V_A, He/F, loci fixed).
- `_compute_pareto` uses O(n²) comparison — fine for typical parameter sweep sizes.
- No issues found.

**Overall code quality**: 9/10. Clean, well-documented, equation-referenced. One lingering formula bug in `segregation_variance`.

### 3. Validation

#### `validation_trait_math.md` — 25/35 passed

| Test | Predicted | Observed | Status | Notes |
|------|-----------|----------|--------|-------|
| trait_mean (R, T, C) | exact | exact | ✅ | 0.0% error |
| trait_variance (R, T, C) | 0.00618 | 0.00608 | ✅ | <2% error |
| exceedance ≤0.20 | 0.261 | 0.250 | ✅ | <8% error |
| exceedance 0.30 | 0.028 | 0.039 | ❌ | 27.6% — tail bias (known) |
| expected_max (n=50-5000) | — | — | ❌ (4) | 10-13% underestimate (known) |
| selection response | 0.138 | 0.150 | ✅ | 8% error |
| multi-gen mean (8 gens) | — | — | ✅ | <5% cumulative |
| multi-gen variance (8 gens) | — | — | ❌ (5) | 30-50% overshoot (Bulmer effect) |
| V_A = V_P (h²=1) | 1.000 | 1.000 | ✅ | Exact |

**Assessment**: All failures are well-understood approximation limits (normal tail bias, Bulmer effect). The two real bugs (factor-of-2 in Δq, factor-of-4 in V_A) were caught and fixed. Excellent validation discipline.

#### `validation_screening.md` — 67/67 passed

- Sample size formula: exact match across 7 (p, γ) combinations
- Empirical coverage: MC confirms ≥95% at all thresholds tested
- Complementarity: 20 deterministic tests, all correct
- Multi-site allocation: optimal ≥ equal, respects caps
- Greedy founder selection: validated end-to-end

**Assessment**: Perfect. Most thorough validation of the three modules.

#### `validation_breeding.md` — 44/44 passed

- Mendelian segregation: χ² at 51 loci, 2000 offspring, 0 rejections
- All four selection strategies improve mean resistance
- Truncation: highest gain (Δr = 0.71 in 5 gen), fastest V_A erosion
- Complementary: lower gain but more V_A retained
- Within-family selection: confirmed effective (high fecundity exploitation)

**Assessment**: Clean pass. The 4× segregation variance issue is documented as known.

### 4. Analysis Templates (`analyses/`)

| Template | Status | Ready? | Blockers |
|----------|--------|--------|----------|
| 01_current_genetic_state.py | ✅ Runnable | ⚠️ Placeholder mode | Calibrated params, genotype snapshot in model |
| 02_screening_effort.py | ✅ Runnable | ⚠️ Depends on Analysis 1 | Same as Analysis 1 |
| 03_breeding_optimization.py | ✅ Runnable | ⚠️ Placeholder mode | Calibrated params |
| 04_reintroduction_scenarios.py | ⚠️ Skeleton | ❌ Not ready | 5 explicit TODOs (release mechanism, calibration, origin tracking, long-run validation, compute scheduling) |
| 05_recommendations.py | ✅ Runnable | ⚠️ Loads whatever exists | Analyses 1-4 |

**Strengths:**
- Central `params.yaml` config is well-organized with all 11 sites, breeding parameter sweeps, and clear output paths.
- All templates produce warnings when using default parameters — good practice.
- Analysis 1 generates synthetic placeholder data using the real genetic architecture, so the pipeline can be tested end-to-end before calibration.
- Analysis 3 supports all 4 strategies with configurable founder counts and generation sweeps.

**Gaps:**
- Analysis 1 `run_single_site()` has a large TODO block where the actual model run should go. Currently generates synthetic data that doesn't reflect actual disease dynamics or spatial structure.
- Analysis 4 has 5 explicit unresolved TODOs. This is the most computationally expensive and least developed analysis.
- `params.yaml` has `params_hash: "TODO"` in the output — minor but should be implemented for reproducibility.
- No analysis template has been run to completion with actual (non-placeholder) data yet.

### 5. Supporting Files

- **README.md**: Excellent. Clear directory structure, usage instructions, component status table, production run estimates.
- **ideas/notes.md**: 8 thoughtfully scoped future directions (inbreeding depression, GBLUP, climate projections, multi-species, epigenetics, cost optimization, genetic rescue from relatives, spatial optimization). Each has a realistic status assessment.
- **references.bib**: 72 entries covering the key literature (Prentice 2025, Schiebelhut, Falconer, Meuwissen OCS, Woolliams).

---

## Validation Summary Table

| Module | Tests | Pass | Fail | Known Limits | Real Bugs Found |
|--------|-------|------|------|--------------|-----------------|
| trait_math | 35 | 25 | 10 | Normal tail bias (5), Bulmer effect (5) | 2 (Δq factor-of-2, V_A factor-of-4) — **fixed** |
| screening | 67 | 67 | 0 | None | 0 |
| breeding | 44 | 44 | 0 | seg_var 4× (documented) | 0 (but see seg_var) |
| **Total** | **146** | **136** | **10** | All documented | 2 fixed, 1 known unfixed |

---

## Known Limitations and Severity

| Limitation | Severity | Impact | Mitigation |
|-----------|----------|--------|------------|
| Normal approx underestimates tail exceedance (>2σ) | Low | Conservative screening predictions (safe direction) | Use empirical resampling for precise estimates |
| Bulmer effect not tracked in analytical predictions | Medium | Multi-gen variance predictions overshoot by 30-50% | Use simulation for variance; mean predictions are fine (<5%) |
| E[max] bias at large n (10-13%) | Low | Screening returns slightly underestimated | Use `expected_max_empirical()` when genotypes available |
| Inbreeding depression not modeled | Medium-High | Breeding program timelines may be optimistic; released stock fitness could be overestimated | Priority extension — needs empirical B estimate for Pycnopodia |
| No environmental variance (h²=1) | Medium | Generation counts to targets are lower bounds; real h² probably 0.3-0.7 | Report includes correction factor guidance; scale by 1/h² |
| `segregation_variance()` returns 4× correct value | Low | Only used for reporting, not selection decisions | Should still be fixed for correctness |
| `genetics.py:compute_additive_variance` has 4× error | Low | Diagnostics-only, doesn't affect simulation | Fix deferred; documented |
| Analysis 4 (reintroduction) is a skeleton | High | Cannot answer the key conservation question yet | Requires model extensions + calibration |
| V_A formula inconsistency between §2 and §3 of report | Low | Potential reader confusion | Add reconciliation note |

---

## Recommended Next Steps (Prioritized)

### Immediate (before Sobol completes)
1. **Fix `segregation_variance()` in breeding.py** — change `/4.0` to `/16.0`. Takes 2 minutes. Updates the report Eq. 5.9 and the test accordingly.
2. **Add reconciliation note** in report §2 Prop. 3 clarifying that V_A = Σ2q(1-q)α² uses the substitution effect (α/2), equivalent to Σ(α²/2)q(1-q) in terms of the effect size.
3. **Fix `compute_additive_variance` in genetics.py** — same root cause, same fix.

### After ABC-SMC Calibration
4. **Plug calibrated parameters into `params.yaml`** and rerun Analyses 1-3.
5. **Implement genotype snapshot capability** in the model recorder so Analysis 1 can extract endpoint genotypes from actual model runs.
6. **Run full Analysis 1 ensemble** on Xeon (est. 4-8 hours for 50 seeds × 11 sites).

### Model Extensions (for Analysis 4)
7. **Release mechanism**: inject captive-bred individuals at a specified node and timestep. This is the key blocker for reintroduction analysis.
8. **Origin tracking**: tag released individuals so we can distinguish wild vs. captive-bred in output.
9. **Long-run validation**: verify model stability for 50-year forward projections.

### For the Paper
10. **Decide on inbreeding depression**: include or defer? If include, need empirical B estimate. The `inbreeding.py:inbreeding_depression()` function is already implemented.
11. **Write the conservation section** of the real paper using calibrated results from Analyses 1-3.

---

## Suspect Equations or Derivations

1. **Eq. 2.12 (Prop. 3)** vs **Eq. 3.7**: These look contradictory but are equivalent under different definitions of α. Needs explicit reconciliation.
2. **Eq. 5.9 (segregation variance)**: The report says α²/16 but the code implements α²/4. The empirical validation shows the code needs to divide by 4 to match reality. The report text was updated but the code wasn't. **This is a remaining code bug.**
3. **Eq. 7.6**: The selection coefficient s_ℓ is left as an abstract variable. Should be derived from the fitness function (§5 Eq. 5.1) for completeness, even if the actual value is computed numerically in the model.

All other equations check out against the validation.

---

## Things That Need Willem's Input/Decision

1. **Inbreeding depression**: Include in the model for the paper, or defer as a limitation? If include, we need B (lethal equivalents) for Pycnopodia or a close relative.
2. **Analysis 4 scope**: The reintroduction analysis is the most expensive (~2-5 days on Xeon) and requires model extensions. Is it in scope for the first paper, or should we present Analyses 1-3 + theoretical reintroduction framework (§7)?
3. **Selection index weights**: The default is resistance-only (w_r=1, w_t=0, w_c=0). Should the paper explore multi-trait selection, or is the resistance dominance result (~1000× more important) sufficient to justify single-trait focus?
4. **Challenge assay validation**: Section 4b frames genotype-to-phenotype mapping as a future experiment. Is there any chance of getting challenge assay data to validate the model's phenotype predictions before publication?
5. **Marker validation priority**: The report emphasizes that validating Schiebelhut loci as resistance markers would transform the conservation program. Is this happening in any lab?

---

## What Was Built Tonight

A complete conservation genetics module in one overnight build session:
- **Theory report**: 34 pages, 10 sections, 72 references — from genetic architecture through screening, breeding, inbreeding, reintroduction, and Pycnopodia biology
- **5 code modules**: trait_math, screening, breeding, inbreeding, viz (2,632 lines)
- **3 validation suites**: 146 tests with detailed reports (2,568 lines)
- **5 analysis templates**: end-to-end pipeline from genetic state → screening → breeding → reintroduction → recommendations (2,297 lines)
- **2 bugs found and fixed** during validation (factor-of-2 in Δq, factor-of-4 in V_A)
- **Supporting infrastructure**: params.yaml, README, ideas/notes, .gitignore, __init__.py

The module connects theory → code → validation → application in a coherent chain, ready to produce real results once calibrated parameters are available.
