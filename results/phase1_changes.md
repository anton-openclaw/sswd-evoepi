# Phase 1: R→S Immunity Fix — Change Summary

**Date**: February 21, 2026
**Issue**: Recovered stars were permanently immune (DiseaseState.R). Echinoderms lack adaptive immunity — stars treated for SSWD can be reinfected. Recovered individuals must return to Susceptible.

## Files Changed

### 1. `sswd_evoepi/types.py`
- **DiseaseState.R comment**: `# Recovered (immune)` → `# Recovered (returns to S — no adaptive immunity in echinoderms)`
- **DiseaseState docstring**: Updated compartment flow `I2 → R (recovery)` → `I2 → S (recovery; no lasting immunity)`
- **AGENT_DTYPE**: Updated `pathogen_virulence` comment `(0 when S or R)` → `(0 when S)`
- **DiseaseState.R enum value kept** (R=5) for backward compatibility with saved snapshots and serialization

### 2. `sswd_evoepi/disease.py` (core changes)
- **I₂ recovery (line ~962)**: `ds[rec_idx] = DiseaseState.R` → `ds[rec_idx] = DiseaseState.S`
- **I₁ early recovery (line ~982)**: `ds[rec_idx] = DiseaseState.R` → `ds[rec_idx] = DiseaseState.S`
- **Module docstring**: Updated from "SEIPD+R" to "SEIPD", recovery target `I₁/I₂ → R` → `I₁/I₂ → S`
- **Recovery function docstrings**: `recovery_probability_I2` and `recovery_probability_I1` updated to document R→S behavior
- Added inline comments at both recovery sites explaining the biological rationale

### 3. `sswd_evoepi/reproduction.py`
- **Spawning eligibility comment**: Updated to note `R→S means recovered agents are already S; R check kept for snapshot compat`
- **No behavioral change**: The `| (disease == DiseaseState.R)` clause is now dead code (n_R will always be 0) but kept for backward compatibility with any loaded snapshots

### 4. `sswd_evoepi/model.py`
- **Same comment update** as reproduction.py for the spawning eligibility check

## What Was NOT Changed (and why)

| Item | Reason |
|------|--------|
| `DiseaseState.R` enum value | Backward compat with saved snapshots/serialization |
| `NodeDiseaseState.n_R` field | Will naturally read 0; removing would break interfaces |
| `SPEED_MODIFIER[5]`, `FEEDING_MODIFIER[5]`, `CAN_SPAWN[5]` | These arrays indexed by disease state still need R=5 slot. Values (1.0, 1.0, True) are correct — recovered stars behave normally, and now they're just S anyway |
| `sswd_evoepi/viz/disease.py` | Visualization code handles R-state agents from old snapshots; no change needed |
| `sswd_evoepi/snapshots.py` | Records raw `disease_state` values; backward compatible |
| `cumulative_recoveries` counter | Still tracks recovery events correctly (counts transitions, not R-state population) |
| Compartment counting (bincount) | `n_R` will naturally be 0 since no agent is ever set to R |

## Expected Behavioral Impact

1. **Recovered stars can be reinfected** — this is the primary change
2. **`n_R` will always be 0** in new simulations (no agents accumulate in R)
3. **Higher reinfection rates** → expect higher cumulative crash rates at population level
4. **Recovery still provides a reproductive window** — a star that recovers can spawn before getting reinfected
5. **Genetic selection for recovery trait still works** — recovering gives fitness advantage even without lasting immunity
6. **`cumulative_recoveries`** still accurate (counts transition events, not compartment size)

## Tests NOT Yet Run
Phase 2 will update test assertions (tests expecting `DiseaseState.R` must change to `DiseaseState.S`) and run the full suite.
