# Pathogen Co-Evolution — Build Report

**Date:** 2026-02-17  
**Builder:** Anton (cron pipeline)  
**Repository:** https://github.com/anton-openclaw/sswd-evoepi  

---

## Phase-by-Phase Summary

### Phase 1+2: Data Structures + Tradeoff Functions
**Commit:** `61f0f1a` — `feat: pathogen evolution phases 1+2 — data structures + tradeoff functions`

- Added `pathogen_virulence` field to AGENT_DTYPE (float32, default 0.0)
- Added `PathogenEvolutionSection` dataclass to config.py (13 parameters)
- Implemented 4 tradeoff functions: `sigma_1_strain()`, `sigma_2_strain()`, `mu_I2D_strain()`, `mu_I1I2_strain()`
- All functions return base rate when PE disabled (backward compatible)
- Config validation: alpha_kill > 0, sigma_v_mutation ≥ 0, v_min < v_max, v_init in bounds

### Phase 3: Strain Inheritance at Transmission
**Commit:** `1ab01c7` — `feat: pathogen evolution phase 3 — strain inheritance at transmission`

- Source attribution: new infections inherit virulence from shedder or environmental reservoir
- Shedder-weighted sampling: probability proportional to individual shedding contribution
- Mutation at transmission: v' = clip(v + N(0, σ_v), v_min, v_max)
- Modified `daily_disease_update()` Step 2 (S→E transitions)

### Phase 4: Virulence-Dependent Dynamics
**Commit:** `cc2c0ce` — `feat: pathogen evolution phase 4 — virulence-dependent dynamics`

- Per-individual disease rates based on strain virulence
- High-v strains: faster shedding, faster kill, faster progression
- Low-v strains: slower dynamics, hosts survive longer
- TLO (Total Lifetime Output) approximately constant across virulence range
- Modified `daily_disease_update()` Steps 1 and 3

### Phase 5: Output Recording + pe_cfg Wiring
**Commit:** `774c308` — `feat: pathogen evolution phase 5 — output recording + pe_cfg wiring`

- `yearly_mean_virulence`, `yearly_virulence_new_infections`, `yearly_virulence_of_deaths` in CoupledSimResult
- PE config wired into `run_coupled_simulation()` and `run_spatial_simulation()`
- Per-node virulence tracking in spatial sim
- New infections assigned v_init when PE enabled

### Phase 6: Integration, Validation & Benchmark (this phase)
- Backward compatibility test: ✅ PASS (bit-identical results)
- Evolutionary dynamics test: ✅ PASS (virulence decreases from 0.8)
- No-mutation stability test: ✅ PASS (v stays at 0.5)
- Performance benchmark: ✅ PASS (negative overhead — within noise)
- Full test suite: ✅ 401/401 PASS

---

## Test Results

### Pathogen Evolution Tests (6/6 passed)
| Test | Status | Notes |
|------|--------|-------|
| `test_backward_compatibility_pathogen_evo` | ✅ PASS | Exact match: deaths, pop, resistance |
| `test_pe_disabled_deterministic` | ✅ PASS | 3 runs identical |
| `test_virulence_evolves_toward_intermediate` | ✅ PASS | v decreases from 0.8 |
| `test_pe_enabled_records_virulence` | ✅ PASS | Metrics recorded correctly |
| `test_no_mutation_stable_virulence` | ✅ PASS | σ=0 → v stays at v_init |
| `test_pe_overhead_acceptable` | ✅ PASS | -3.0% overhead |

### Full Test Suite (401/401 passed)
| Module | Tests | Status |
|--------|-------|--------|
| test_types.py | 26 | ✅ |
| test_config.py | 26 | ✅ |
| test_disease.py | 115 | ✅ |
| test_genetics.py | 60 | ✅ |
| test_reproduction.py | 73 | ✅ |
| test_rng.py | 13 | ✅ |
| test_movement.py | 33 | ✅ |
| test_spawning.py | 55 | ✅ |
| **Total** | **401** | **✅ ALL PASS** |

---

## Performance Benchmark

```
=== Pathogen Evolution Performance Benchmark ===
Disabled: 1.454s
Enabled:  1.410s
Overhead: -3.0%
```

**Assessment:** Negligible overhead. The extra per-individual `exp()` calls and weighted sampling at transmission are dominated by existing computation. Well within the <10% target.

---

## Git Log (Pathogen Evolution)

```
774c308 feat: pathogen evolution phase 5 — output recording + pe_cfg wiring
cc2c0ce feat: pathogen evolution phase 4 — virulence-dependent dynamics
1ab01c7 feat: pathogen evolution phase 3 — strain inheritance at transmission
61f0f1a feat: pathogen evolution phases 1+2 — data structures + tradeoff functions
6680f3d spec: pathogen co-evolution with mechanistic virulence-transmission tradeoff
```

---

## Known Issues

1. **Environmental virulence (v_env) is static** — the environmental reservoir doesn't evolve. In reality, shed pathogen virulence should influence environmental pool composition. Low priority: environmental transmission is a small fraction of total FOI.

2. **Years with 0 infections record 0.0 for mean virulence** — not NaN. Tests handle this, but downstream analysis code should filter for `yearly_mean_virulence > 0`.

3. **Lineage tracking (track_strain_history) not tested** — optional feature, disabled by default. Will need tests if used for phylogenetic analysis.

4. **Spatial PE metrics** — per-node virulence tracking exists in `run_spatial_simulation()` but not yet tested in spatial integration tests.

---

## Next Steps

1. **Spatial PE integration tests** — verify virulence evolution across multi-node networks
2. **Density-dependent virulence shift** — test that optimal virulence changes pre vs post epidemic (spec §7.2)
3. **Myxomatosis pattern test** — v_init=1.0, confirm convergence toward intermediate (classic analog)
4. **Joint host-pathogen coevolution** — multi-seed runs showing coexistence vs extinction
5. **Clement 2024 comparison** — once Willem provides the paper, validate against DFTD eco-evo predictions
6. **Conservation module** — model captive-bred outplanting with PE-aware disease dynamics
