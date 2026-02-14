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
