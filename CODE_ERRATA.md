# CODE_ERRATA.md — Implementation Errata Tracker

**Purpose:** Track any issues found during implementation that affect other build phases.
Each entry is appended by the phase that discovers it. Subsequent phases MUST read this file first.

---

## Format

```
### CE-<N>. <Title>
**Found by:** Phase <X>
**Date:** YYYY-MM-DD
**Severity:** 🔴 HIGH | 🟡 MEDIUM | 🟢 LOW
**Affects:** <list of phases/modules>

**Issue:** <description>
**Resolution:** <fix or workaround>
**Status:** ⚠️ Open | ✅ Resolved
```

---

## Entries

### CE-1. Cost of Resistance Removed (Willem's Decision)
**Found by:** Phase 1 (Infrastructure)
**Date:** 2026-02-13
**Severity:** 🟡 MEDIUM
**Affects:** All phases — genetics, population, config

**Issue:** Willem decided to remove cost of resistance entirely. The `cost_resistance` parameter
from the specs is NOT included in the config. `fecundity_mod` field in AGENT_DTYPE is retained
but always set to 1.0. All code referencing `cost_resistance` or `fecundity_mod < 1.0` must
skip cost-of-resistance logic.

**Resolution:** 
- `genetics.cost_resistance` removed from default.yaml
- `fecundity_mod` initialized to 1.0 for all agents
- Fecundity calculation uses size only: `F(L) = F0 * (L/L_ref)^b`
- If cost is re-enabled later, add `genetics.cost_resistance` back to config

**Status:** ✅ Resolved (by design)

### CE-2. Both Etiological Scenarios as Config Option
**Found by:** Phase 1 (Infrastructure)
**Date:** 2026-02-13
**Severity:** 🟢 LOW
**Affects:** Disease module (Phase 2+)

**Issue:** Willem wants both "ubiquitous" and "invasion" etiological scenarios testable via config.
The `disease.scenario` field accepts either value, which changes how the pathogen is introduced.

**Resolution:** `disease.scenario` field in config accepts "ubiquitous" or "invasion". Default: "ubiquitous".
- **ubiquitous:** Vibrio always present at low background levels; epidemic triggered by SST crossing threshold
- **invasion:** Vibrio absent until explicit introduction event at a configured year/node

**Status:** ✅ Resolved (config supports both)

### CE-3. Exponential Decay for Effect Sizes (Willem's Decision)
**Found by:** Phase 1 (Infrastructure)
**Date:** 2026-02-13
**Severity:** 🟢 LOW
**Affects:** Genetics module (Phase 3)

**Issue:** Willem confirmed exponential distribution for effect sizes (per Lotterhos & Whitlock 2016).
Already specified in genetics-evolution-spec §2.1. No change needed.

**Resolution:** Effect sizes drawn from Exp(λ), normalized to sum to W_add ≈ 0.840. Sorted descending.

**Status:** ✅ Resolved (matches spec)

### CE-4. Beverton-Holt Denominator Corrected
**Found by:** Phase 3 (Reproduction)
**Date:** 2026-02-13
**Severity:** 🔴 HIGH
**Affects:** Population module, settlement, all demographic calibration

**Issue:** The population-dynamics-spec.md §5.2 BH formula uses denominator
`1 + S/(K × s0)`, which gives asymptote `K × s0² ≈ 0.45` — meaning zero
recruits ever survive (rounds to 0 with s0=0.03, K=500). The spec text
claims asymptote is `K × s0 = 15`, which is the correct biological target.

**Resolution:** Corrected BH formula to `R = S × s0 / (1 + S/K)`, which gives:
- Low S: R ≈ S × s0 (supply-limited, correct)
- High S: R → K × s0 = 15 (correct asymptote)
- spec.md formula had wrong denominator; text description was correct

**Status:** ✅ Resolved in reproduction.py

### CE-6. Saprophytic Carcass Shedding Rate Reduced (σ_D: 150 → 15)
**Found by:** Phase 4 (Disease)
**Date:** 2026-02-13
**Severity:** 🔴 HIGH
**Affects:** Disease module, R₀ computation, epidemic calibration

**Issue:** The spec-derived field-effective σ_D = 150 bact/mL/d/carcass was 10× larger than
σ₂_eff = 50 (symptomatic shedding at 20°C) when accounting for the 3-day carcass duration.
Each carcass produced 450 total bact/mL — more pathogen than the entire I₂ stage (289 at 20°C).
This created a dominant positive feedback loop (death → saprophytic burst → more infections)
that made R₀ >> 1 even at 8°C (cold temperatures where epidemics should not sustain).

The R₀ formula in the spec §6.6 excluded carcass contribution, masking this issue:
- Spec R₀ at 8°C/φ=0.1: 0.72 (sub-threshold, correct target)
- Actual R₀ with carcass (σ_D=150): 2.34 (epidemic, wrong)

**Resolution:**
- Reduced σ_D from 150 to 15 bact/mL/d/carcass (★☆☆ confidence parameter)
- R₀ computation now includes carcass term: σ_D × τ_D in shedding integral
- At σ_D=15: carcass total = 45, ~15% of I₂ stage — biologically minor, as intended
- Corrected R₀ at 8°C/φ=0.1: ~0.80 (sub-threshold ✓)
- Corrected R₀ at 16°C/φ=0.02: ~2.1 (epidemic ✓)
- Updated default.yaml and config.py

**Status:** ✅ Resolved

### CE-5. High-Fecundity Broadcast Spawner Allee Effect
**Found by:** Phase 3 (Reproduction)
**Date:** 2026-02-13
**Severity:** 🟡 MEDIUM
**Affects:** Allee effect interpretation, population recovery dynamics

**Issue:** The population-dynamics-spec.md §6.2 states that at post-SSWD
density (ρ=0.001), per-capita growth `r(D) < 0` — a deterministic Allee
threshold. However, with F0=10⁷ eggs/female, the deterministic threshold
is effectively zero: even 0.01% fertilization produces enough recruits.

The practical Allee effect operates through:
1. Dramatically reduced growth RATE at low density (>100× reduction)
2. Demographic stochasticity at small N (integer effects, drift)
3. Environmental stochasticity overwhelming the tiny growth rate

The spec's simplified formula `r = F(D) × f_recruit × s0 − m` omits
the fecundity multiplier, giving qualitatively correct shape but wrong
magnitude. The actual per-capita growth rate includes fecundity:
`r = F(D) × (eggs/2) × larval_surv × s0 − m_adult`

**Resolution:** Implementation uses the full formula. Tests verify Allee
SHAPE (growth rate increases with density, >100× difference between
post-SSWD and healthy density) rather than absolute zero-crossing.
Population decline at low density emerges from stochastic simulation,
not deterministic formula.

**Status:** ✅ Resolved — documented as biological reality for high-fecundity species
