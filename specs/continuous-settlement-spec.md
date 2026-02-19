# Continuous Settlement — Specification

**Date:** 2026-02-18
**Authors:** Willem Weertman & Anton 🔬
**Status:** DRAFT
**Priority:** HIGH — fixes artificial sawtooth epidemic dynamics

---

## Problem

The epidemic curve shows sharp, periodic sawtooth spikes: population crashes to near-zero
between recruitment events, then spikes when a new cohort arrives, then crashes again as
disease burns through fresh susceptibles. This pattern is **artificial** — an artifact of
the model's temporal structure, not biology.

**Root cause:** Spawning occurs daily during the season (Nov–Jul), and larvae accumulate
in a list (`accumulated_cohorts`), but ALL settlement happens as a single annual pulse
during "Phase B" at year-end. The disease module (daily timestep) then encounters a
sudden influx of naïve susceptibles on Day 1 of the next year → instant epidemic spike
→ rapid die-off → wait for next annual pulse.

**What it should look like:** Larvae spawned on different days reach competency at
different times (PLD is temperature-dependent: 48–71 days depending on SST). Settlement
should be spread across weeks/months as cohorts sequentially reach competency, producing
a gradual influx of recruits rather than a single pulse.

---

## Biological Basis (Hodin et al. 2021)

### Larval Development Timeline (10–11°C)
| Stage | Timing |
|-------|--------|
| Fertilization | 0 h |
| Feeding bipinnaria | 6 dpf |
| First juvenile skeleton | 27 dpf |
| Competent larva | 49–55 dpf |
| First spontaneous settlement | ~49 dpf (7 weeks) |

### Temperature Dependence
Current model (from specs):
```
PLD(T) = 63 × exp(−0.05 × (T − 10.5))   # days
```

| SST (°C) | PLD (days) | Context |
|----------|------------|---------|
| 8.0 | 71 | Sitka, AK (cold) |
| 10.0 | 65 | Howe Sound, BC |
| 10.5 | 63 | Reference (Hodin lab) |
| 12.0 | 58 | Newport, OR |
| 14.0 | 53 | Monterey, CA |

### Settlement Biology
- Best inducer: conspecific biofilm (adult *P. helianthoides* tank biofilm)
- Settlement is NOT instantaneous — 4 phases: exploration → attachment → settler → settled
- Delayed metamorphosis possible (larvae held >100 dpf still settled)
- **Competency window**: larvae become competent at ~49 dpf but can delay settlement
  for weeks if conditions are poor

### Spawning Season
- Females: November–May (WA); model uses Nov–Jul (270 days)
- Males: sperm available year-round
- Multiple spawning bouts per female per season (2–3×)
- Spawning is asynchronous — not all females spawn on the same day

### Implication
A female spawning on Day 1 of the season (November) at 10°C produces larvae that
reach competency ~65 days later (mid-January). A female spawning in March produces
competent larvae in May. Settlement should span November through September, with
peak settlement trailing peak spawning by ~2 months.

---

## Design

### Core Change

**Each larval cohort tracks its spawn day and PLD.** When the cohort's age
(current_day − spawn_day) reaches its PLD, the larvae become competent and settle
into the population on that day — during the daily loop, not at year-end.

### Data Structure Changes

Add to `LarvalCohort`:
```python
@dataclass
class LarvalCohort:
    source_node: int
    n_competent: int
    genotypes: np.ndarray       # (n_competent, N_LOCI, 2)
    parent_pairs: np.ndarray    # (n_competent, 2)
    pld_days: float
    spawn_day: int              # NEW: absolute simulation day when spawned
    sst_at_spawn: float         # NEW: SST at spawning (determines PLD)
```

### Settlement Logic

**In the daily loop** (not Phase B), after disease update and before spawning step:

```python
# Check if any accumulated cohorts have reached competency
ready_cohorts = []
still_waiting = []
for cohort in accumulated_cohorts:
    age = current_sim_day - cohort.spawn_day
    if age >= cohort.pld_days:
        ready_cohorts.append(cohort)
    else:
        still_waiting.append(cohort)
accumulated_cohorts = still_waiting

# Settle ready cohorts
if ready_cohorts:
    settle_cohorts(ready_cohorts, agents, genotypes, ...)
```

### Settlement Function

```python
def settle_daily_cohorts(
    cohorts: List[LarvalCohort],
    agents: np.ndarray,
    genotypes: np.ndarray,
    carrying_capacity: int,
    pop_cfg: PopulationSection,
    rng: np.random.Generator,
) -> int:
    """Settle competent larval cohorts into the population.
    
    Applies Beverton-Holt density-dependent recruitment per cohort.
    Returns total number of recruits settled.
    """
    total_settled = 0
    for cohort in cohorts:
        current_alive = int(np.sum(agents['alive']))
        available_slots = max(0, carrying_capacity - current_alive)
        if available_slots == 0:
            break
        
        # Beverton-Holt recruitment
        n_recruits = beverton_holt_recruitment(
            cohort.n_competent, carrying_capacity, pop_cfg.settler_survival
        )
        n_recruits = min(n_recruits, available_slots, cohort.n_competent)
        
        # Settlement cue modifier
        cue_mod = settlement_cue_modifier(current_alive)
        n_recruits = max(0, int(n_recruits * cue_mod))
        
        if n_recruits > 0:
            # Place into dead slots (same logic as current annual_reproduction)
            ...
            total_settled += n_recruits
    
    return total_settled
```

### Spatial Sim Changes

The spatial sim (Phase B in `run_spatial_simulation`) currently does:
1. Annual mortality
2. Annual reproduction → produce larvae
3. Larval dispersal via C matrix
4. Settle at destination nodes

With continuous settlement, this becomes:

**Daily loop additions:**
- Each node maintains a list of pending larval cohorts (local + received via dispersal)
- Daily: check each pending cohort's age vs PLD → settle if ready
- Beverton-Holt recruitment applied per cohort, not aggregated

**Annual dispersal stays the same** (C matrix operates on annual timescale for
larval transport). But the settlement of dispersed larvae is now spread across
the following weeks as each cohort's PLD elapses.

**Alternative (simpler):** For the coupled (single-node) sim, settlement happens
locally with no dispersal matrix. Each `spawning_step()` cohort gets tagged with
`spawn_day` and `pld_days`, accumulated, and settled individually when ready.

### Phase B Simplification

Phase B (annual demographic update) loses the recruitment section entirely. It retains:
- B1: Natural mortality (annual)
- B2: Growth and aging (annual)
- B3: (removed — recruitment now daily)
- B4: Genetics recording
- B5: Disease seeding (kept at year boundary for initial introduction)

### Spawning Step Update

`spawning_step()` currently returns `List[LarvalCohort]`. Each cohort needs:
- `spawn_day = current_sim_day` (absolute day counter)
- `sst_at_spawn = current_sst` → used to compute PLD
- `pld_days = pelagic_larval_duration(sst_at_spawn)`

These are set when the cohort is created, not when it settles.

---

## Settlement Spreading — Expected Behavior

### Example: 5-Node, 20-Year Simulation

Spawning season: Nov 1 (Day 305) through Jul 31 (Day 212 next year) = ~270 days.

At Sitka (8°C winter, 11°C summer):
- Nov 1 spawn → PLD=71d → competent Jan 11
- Dec 1 spawn → PLD=70d → competent Feb 9
- Mar 1 spawn → PLD=65d → competent May 5
- Jun 1 spawn → PLD=58d → competent Jul 29

At Monterey (12°C winter, 16°C summer):
- Nov 1 spawn → PLD=58d → competent Dec 29
- Mar 1 spawn → PLD=53d → competent Apr 23
- Jun 1 spawn → PLD=48d → competent Jul 19

**Result:** Settlement occurs over ~8 months (Jan through Aug at Sitka,
Dec through Jul at Monterey) instead of a single pulse. Each week brings
a small number of settlers → disease encounters them gradually → smooth
endemic dynamics instead of sawtooth.

---

## Competency Window (Future Extension)

Hodin 2021 showed larvae can delay metamorphosis for >100 days past competency.
This means competent larvae don't settle immediately — they search for suitable
substrate (conspecific biofilm, coralline algae).

**For now:** Larvae settle on the day they reach competency (PLD elapsed).
**Future:** Add a competency window [PLD, PLD + max_delay] where settlement
probability depends on conspecific density (settlement cue modifier) and
increases with age (desperation settlement). This would further spread
settlement and add density-dependent feedback.

---

## Impact on Existing Tests

### Tests that will break:
- Tests comparing exact population counts at year boundaries (settlement timing changes)
- Tests that check `repro_diag['n_recruits']` from Phase B (no longer exists there)
- Any test assuming zero-population-change during daily loop (settlers now arrive daily)

### Tests that should still pass:
- Disease module tests (independent of settlement)
- Genetics tests (genotype inheritance unchanged)
- Pathogen evolution tests (independent)
- Equilibrium tests (long-run, will settle to same K)

### New tests needed:
- Settlement timing matches PLD for given SST
- No settlement before PLD elapsed
- Multiple cohorts settle on different days
- Beverton-Holt still regulates at K
- Continuous settlement produces smoother population trajectory than annual
- Spawned-in-November settlers arrive before spawned-in-March settlers
- Temperature gradient → PLD gradient → settlement timing gradient

---

## Performance Considerations

- **Daily cohort checking:** O(n_pending_cohorts) per day per node. Typically
  <100 pending cohorts → negligible.
- **Daily settlement:** One Beverton-Holt call + slot allocation per settling cohort.
  Much cheaper than disease update → negligible.
- **Memory:** Each pending cohort stores genotypes. At peak, maybe 50–100 cohorts
  pending per node × genotype arrays. Already happens with accumulation — no change.

**Expected net overhead: < 1%**

---

## Implementation Phases

### Phase 1: Data Structure + Coupled Sim (2 hours)
- Add `spawn_day` and `sst_at_spawn` to LarvalCohort
- Set them in `spawning_step()` output
- Move settlement from Phase B into daily loop
- Write `settle_daily_cohorts()` function
- Remove annual settlement from Phase B (coupled sim only)
- Tests: settlement timing matches PLD, smoother trajectory

### Phase 2: Spatial Sim (2 hours)
- Same daily settlement logic for spatial sim
- Larval dispersal still annual (C matrix), but received cohorts carry spawn_day/PLD
- Dispersed cohorts settle at destination on their individual PLD date
- Tests: spatial settlement spread, node-level timing

### Phase 3: Validation + Visualization (1 hour)
- Re-run 5-node 20yr simulation
- Generate new epidemic curve, vibrio, recruitment figures
- Compare old (sawtooth) vs new (smooth) side-by-side
- Verify biological plausibility of settlement timing

**Estimated total: ~5 hours**

---

## Key References

1. **Hodin et al. 2021** — Complete life-cycle culturing of *Pycnopodia*;
   PLD = 49–55 dpf at 10–11°C; delayed metamorphosis >100 dpf
2. **Population dynamics spec** — Beverton-Holt recruitment, settler_survival
3. **Spawning overhaul spec** — Daily spawning events, immunosuppression
