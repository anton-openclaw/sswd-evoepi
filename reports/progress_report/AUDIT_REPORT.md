# Progress Report Audit — Combined Findings
**Date**: 2026-03-09
**Report**: `reports/progress_report/report.tex` (v8, 23 pages)
**Auditors**: 3 Opus subagents (disease, genetics, spatial) + manual numbers audit
**Method**: Every claim cross-referenced against codebase (`sswd_evoepi/`, `experiments/`, `data/`)

---

## 🔴 HIGH SEVERITY (must fix)

### H1. Resistance (r) ≠ infection avoidance probability
**Lines**: 120, 129, 136, 371-382 (breeding table), Fig 2 caption
**Report says**: "r = probability of fighting off infection", "r = 0.25 means 75% infection"
**Code does**: `λ = a × P/(K_half+P) × (1−r) × S_sal × f_size(L)`, then `p_inf = 1 − exp(−λ)`. Resistance is ONE factor in a hazard rate. At r=0.25 with pathogen at K_half: actual p_inf ≈ 24%, not 75%.
**Impact**: Pervades entire report. The breeding program table (§7) lists "Infection probability = 1−r" which is systematically wrong. Would mislead wildlife managers designing screening protocols.
**Flagged by**: All 3 auditors independently

### H2. Recovery trait (c) described as daily clearing probability — actually 20× lower
**Lines**: 133, 136, 483 (appendix)
**Report says**: "c = 0.02 (only a 2% daily chance of clearing infection)"
**Code does**: `p_rec = ρ_rec × c_i = 0.05 × 0.02 = 0.001` (0.1%/day). The trait c is a multiplier on base rate ρ_rec, not the probability itself.
**Impact**: 20× error in the stated recovery rate. Fundamentally changes understanding of recovery's role.
**Flagged by**: All 3 auditors independently

### H3. δ_env decay rate: "2% per month" vs actual 2% per day
**Lines**: 184 (body) vs 447 (appendix table — correct)
**Report says**: "about 2% of the reservoir decays per month"
**Code does**: `pool_decay = cfg.delta_env * P_env_pool` applied DAILY. 2%/day = ~45%/month. Half-life 35 days, not 35 months.
**Impact**: 30× error. Self-contradicts the appendix table which correctly says "per day".
**Flagged by**: All 3 auditors independently

### H4. Virulence gradient claimed as "emergent" — actually formula-prescribed
**Lines**: 216, 219
**Report says**: "We did not prescribe this gradient—it emerges naturally from the evolutionary dynamics"
**Code does**: `v_optimal = v_max_warm × density_ratio × temp_factor`, then `delta_v = v_adapt_rate × (v_optimal − v_local)`. The gradient IS directly prescribed by a temperature × density formula. No strain competition, no R₀ optimization, no trade-off curve drives the evolution. Simple deterministic drift toward a computed target.
**Impact**: The central narrative claim about the model's flagship result is the opposite of what the code does. A reviewer familiar with Anderson & May would catch this immediately.
**Flagged by**: All 3 auditors independently

### H5. Dispersal kernel described as "power-law" — code uses exponential
**Line**: 79
**Report says**: "connectivity decreases with distance according to a power-law function"
**Code does**: `kernel = np.exp(-distances / D_L)` — exponential kernel (spatial.py line 337). Power-law = d^(-α), exponential = e^(-d/D_L). Fundamentally different tail behavior (fat vs thin tails).
**Impact**: Anyone implementing or evaluating this model would use the wrong dispersal function. Changes all conclusions about long-range connectivity.
**Flagged by**: Spatial auditor

### H6. "Monthly timesteps" — code runs daily
**Line**: 65
**Report says**: "The simulation runs in monthly timesteps over 13 years"
**Code does**: `for day in range(DAYS_PER_YEAR)` with DAYS_PER_YEAR=365. All disease, mortality, and growth processes are daily. Monthly coordination exists only for recording.
**Impact**: Daily vs monthly timesteps produce very different disease dynamics for a disease with ~11-day exposure-to-death timescale. The daily resolution is actually a strength the report undersells.
**Flagged by**: Genetics auditor

---

## 🟡 MEDIUM SEVERITY (should fix)

### M1. Warm-site T_vbnc "evolves higher thresholds" — no mechanism for upward evolution
**Line**: 204
**Report says**: "Bacteria at warm-water sites evolve higher thresholds"
**Code does**: T_vbnc can only evolve DOWN (toward T_vbnc_min) or revert to T_vbnc_initial (at rate 0.0 = disabled). No mechanism exists for evolution above 12°C. Warm sites simply stay near 12°C.
**Flagged by**: Genetics + Disease auditors

### M2. K_half = "50% infection chance" — actually ~27% for typical star
**Line**: ~201
**Report says**: "When the reservoir equals K_half, a susceptible star has a 50% chance of becoming infected"
**Code does**: At P=K_half, dose_response=0.5, but actual p_inf = 1 − exp(−0.75 × 0.5 × 0.85 × ...) ≈ 27%.
**Flagged by**: Disease auditor

### M3. T_vbnc_min described as activity cutoff — it's an evolution floor
**Lines**: 165-168, 441
**Report says**: "Below 9°C, bacteria remain in VBNC dormancy—they persist but are metabolically inactive and cannot cause new infections"
**Code does**: T_vbnc_min only limits how far T_vbnc_local can evolve. At T_vbnc_local=9°C and T=9°C, VBNC sigmoid gives 50% activation, not zero. Infections CAN occur below 9°C.
**Flagged by**: Disease auditor + my new thermal adaptation section has this right

### M4. VBNC "threshold" described as binary — actually gradual sigmoid
**Line**: ~150
**Report says**: "When water temperatures drop below a critical threshold... bacteria essentially go dormant"
**Code does**: Sigmoid with k=2.0 at T=12°C. At 11°C = 12% active. At 10°C = 2% active. Not a hard cutoff.
**Flagged by**: Disease auditor

### M5. Spawning season "November to March" — code runs Nov through mid-July
**Line**: 71
**Report says**: "peak spawning occurs between November and March"
**Code does**: season_start_doy=305 (~Nov 1) to season_end_doy=196 (~Jul 15). Peak is Feb but season extends 4+ months further.
**Flagged by**: Spatial auditor

### M6. Self-recruitment presented as binary (fjord/open) — actually continuous
**Lines**: 500-501
**Report says**: α_open=0.05, α_fjord=0.50 (as if two categories)
**Code does**: Continuous interpolation: `alpha_self = alpha_open + (alpha_fjord - alpha_open) × effective_enclosedness`. Most nodes get intermediate values.
**Flagged by**: Spatial auditor

### M7. P_env_max described as "shedding rate" — it's background Vibrio input
**Line**: 445
**Report says**: "Maximum rate at which infected/dying stars shed pathogen"
**Code does**: P_env_max is background environmental Vibrio input (config.py comment: "Background Vibrio input bact/mL/d"). Shedding is σ₁, σ₂, σ_D. Also: P_env_max is unused in W142's dynamic P_env path.
**Flagged by**: Spatial + Disease auditors

### M8. Wavefront D_P (300km) conflated with daily D_P (15km)
**Lines**: 456-457
**Report says**: "Dispersal range D_P = 300 km" under "Wavefront Disease Spread"
**Code does**: Two separate dispersal systems: standard daily D_P=15km (max 50km) and wavefront dose accumulation D_P=300km (max 3000km). The report only shows the wavefront values.
**Flagged by**: Spatial + Disease auditors

### M9. "Especially large bursts at death" — I₂ shedding is 3.3× higher
**Line**: ~190
**Report says**: "especially when dying stars release large bursts"
**Code does**: σ_D=15 (death) vs σ₂=50 (I₂ shedding). I₂ shedding dominates, not death bursts.
**Flagged by**: Disease auditor

### M10. "Additive" genetics omits variable effect sizes
**Lines**: 129, 390
**Report says**: "Each of the 17 resistance loci contributes additively"
**Code does**: Effect sizes drawn from Exponential(λ=1), sorted descending. First locus has largest effect. Additive in dosage but with very unequal weights.
**Flagged by**: Genetics auditor

### M11. I₁ recovery threshold (c>0.5) not mentioned
**Line**: 133
**Report says**: Recovery works at any c value (implied)
**Code does**: I₁ recovery requires c_i > 0.5 (C_EARLY_THRESH). At initial c=0.02, NO individuals can recover from early-stage infection.
**Flagged by**: Genetics + Disease auditors

### M12. Two independent drift models described as "single optimization"
**Line**: ~258 (new thermal adaptation section)
**Report says**: "These two axes of adaptation are not independent; they reflect a single underlying optimization"
**Code does**: T_vbnc and v_local evolve via completely separate functions with separate parameters. Correlated through shared temperature input, but not mechanistically coupled.
**Flagged by**: Disease auditor (my own text!)

### M13. Pathogen thermal adaptation described as "selection" — deterministic drift
**Lines**: ~258-270
**Report says**: "cold-tolerant strains are selectively favored"
**Code does**: `T_vbnc_local -= adapt_rate × pool_factor × temp_gap`. No strains, no competition. Deterministic drift.
**Flagged by**: Disease auditor

### M14. v_local starts at 0.5, ALWAYS declines — not bidirectional evolution
**Lines**: 216, 219
**Report says**: Virulence "ranges from 0.08 to 0.35" (implies evolution in both directions)
**Code does**: All sites start at v_local=0.5, drift toward v_optimal which is always ≤0.7 and typically well below 0.5. The gradient is differential DECAY from a high starting point, not bidirectional evolution.
**Flagged by**: Genetics auditor

---

## 🟢 LOW SEVERITY (consider fixing)

### L1. Tolerance only affects I₂→D stage — report implies all stages
### L2. genetics.py default c=0.08 vs config default c=0.02 (code inconsistency)
### L3. "GWAS" vs selection scan (Schiebelhut et al. used RADseq, not traditional GWAS)
### L4. "Return to health" omits that recovered stars are immediately re-susceptible (S→I→S, no lasting immunity)
### L5. K=5000 is calibration runner default, not inherent model property
### L6. Enclosedness computation described too simply (actually multi-component)
### L7. Legacy spawning code in reproduction.py contradicts active spawning.py
### L8. Post-spawning immunosuppression omitted from biology narrative (only in parameter table)
### L9. 150km connectivity in figure caption not traceable to any model parameter
### L10. n_connectivity described as inter-site mixing — it's flushing rate mapping

---

## 🔵 PARAMETER TABLE DISCREPANCIES (appendix)

| Parameter | Report Table | Config Default | W154 Override | Status |
|-----------|-------------|----------------|---------------|--------|
| α_env | **0.20** | 0.1 | **0.18** | ❌ Wrong (neither default nor W154) |
| n_conn | **0.5** | N/A | **0.3** | ❌ Wrong (W154 = 0.3) |
| α_self_open | **0.05** | 0.10 | **0.02** | ❌ Wrong (neither default nor W154) |
| α_self_fjord | **0.50** | 0.30 | **0.70** | ❌ Wrong (neither default nor W154) |
| k_vbnc | 2.0 | 1.0 | 2.0 | ✅ Matches W154 |
| K_half | 800,000 | 87,000 | 800,000 | ✅ Matches W154 |
| P_env_max | 2,000 | 500 | 2,000 | ✅ Matches W154 (but mislabeled) |
| δ_env | 0.02 | 0.05 | 0.02 | ✅ Matches W154 |
| P_floor | 500 | 50 | 500 | ✅ Matches W154 |
| c₀ | 0.02 | 0.02 | 0.02 | ✅ (but mislabeled as "daily probability") |

**4 parameters have wrong values** that match neither the config defaults nor the W154 configuration.

---

## 📰 RESIDUAL La Niña REFERENCES

The section title was corrected to "Temperature-Driven Population Dynamics" but 4 La Niña references remain in captions/text:
- **Line 149**: Fig caption "La Niña period (2016-2017), temporary recovery as cold water suppresses pathogen"  
- **Line 255**: Fig caption "coinciding with the La Niña event that brought cooler water"
- **Line 266**: Body text mentions La Niña as example of captured climate variability (acceptable)
- **Line 336**: "La Niña response" listed as model achievement

Lines 149 and 255 propagate the old narrative that was explicitly corrected.

---

## 📊 AUDIT STATISTICS

| Severity | Count |
|----------|-------|
| 🔴 HIGH | 6 |
| 🟡 MEDIUM | 14 |
| 🟢 LOW | 10 |
| Parameter table errors | 4 |
| La Niña residuals | 2-3 |

**Total unique findings: ~30** (after deduplication across auditors)

**Most critical cluster**: The resistance-as-probability error (H1) + breeding table (H1) + recovery probability error (H2) form an interconnected web that affects the entire report's quantitative narrative.

**Confirmed accurate**: 896 nodes / 18 regions ✅, SRS implementation ✅, spawning formula ✅, immunosuppression duration ✅
