# Implementation Plan — Willem's Feb 21 Feedback

**Deadline**: 9 PM PST, February 21, 2026
**Three issues**: (1) R→S immunity fix, (2) Satellite SST, (3) Paper corrections

---

## Phase 1: R→S Immunity Fix (1:30 PM)
**Goal**: Recovered stars return to Susceptible, not permanently immune.

### Changes needed:
1. **`sswd_evoepi/types.py`**: Change `R = 5  # Recovered (immune)` comment to clarify R is transient tracking state, OR simply set recovered individuals back to S
2. **`sswd_evoepi/disease.py`**: Two recovery locations (I2 recovery ~line 958, I1 early recovery ~line 977) — change `ds[rec_idx] = DiseaseState.R` to `ds[rec_idx] = DiseaseState.S`
3. **Keep `n_R` tracking**: Use `cumulative_recoveries` counter (already exists) instead of snapshot R count. Set `node_state.n_R` tracking to count only recently-recovered (or just let it always read 0 since they return to S immediately)
4. **Behavioral modifiers**: Check line 72 area — if R state has different feeding/movement modifiers, those won't apply anymore (which is correct — recovered stars behave normally)
5. **`tests/test_disease.py`**: Update any tests asserting `DiseaseState.R` for recovered agents → now expect `DiseaseState.S`
6. **Recording**: Check if `yearly_recoveries`, `SpatialSimResult`, or snapshots depend on R-state counting
7. **Run full test suite**: `pytest tests/ -x -q`
8. **Commit + push**

### Biological rationale:
- Echinoderms lack adaptive immunity (no antibodies, no immunological memory)
- Oregon Coast Aquarium protocol (iodine dip + probiotic) saves stars, but they can get SSWD again
- Genetic traits (resistance, tolerance, recovery) represent heritable defense — this is the correct mechanism

---

## Phase 2: SST Data Extraction (3:30 PM)
**Goal**: Download real satellite SST for all 11 nodes from ERDDAP.

### Data source: NOAA OISST v2.1
- Resolution: 0.25° daily
- Coverage: 1981–present
- Access: ERDDAP REST API at `coastwatch.pfeg.noaa.gov`
- Dataset: `ncdcOisst21Agg_LonPM180` (OISST v2.1 aggregation)

### 11 Node coordinates:
| Node | Lat | Lon |
|------|-----|-----|
| Sitka | 57.06 | -135.34 |
| Ketchikan | 55.34 | -131.64 |
| Haida Gwaii | 53.25 | -132.07 |
| Bella Bella | 52.16 | -128.15 |
| Howe Sound | 49.52 | -123.25 |
| SJI | 48.53 | -123.02 |
| Westport | 46.89 | -124.10 |
| Newport | 44.63 | -124.05 |
| Crescent City | 41.76 | -124.20 |
| Fort Bragg | 39.45 | -123.80 |
| Monterey | 36.62 | -121.90 |

### Steps:
1. Write `scripts/fetch_sst_data.py` — ERDDAP query script
2. For each node: extract nearest-grid-point daily SST, 2002–2025
3. Compute daily climatology (365-day mean across years) per node
4. Also store raw daily time series for potential future use
5. Save to `data/sst/` as CSV (one per node + one climatology summary)
6. Validate: plot climatology curves, check for gaps/NaNs
7. Commit data + script

### Fallback: MUR SST (1km, 2002-present) if OISST coastal coverage is poor for fjord nodes (Howe Sound, Bella Bella)

---

## Phase 3: SST Model Integration (5:00 PM)
**Goal**: Model can use real satellite SST instead of sinusoidal.

### Changes:
1. **`sswd_evoepi/environment.py`**: Add `satellite_sst()` function that loads climatology CSV and returns SST by day_of_year + node
2. **`sswd_evoepi/config.py`**: Add `sst_source: str = "sinusoidal"` to config (options: "sinusoidal", "satellite")
3. **`sswd_evoepi/spatial.py`**: In the daily SST computation, dispatch based on `sst_source` config
4. **Keep sinusoidal as default**: No breaking changes to existing configs/scripts
5. **Warming trend overlay**: For future projections, apply configurable trend on top of climatology
6. **Tests**: New tests for satellite SST loading, fallback behavior, config validation
7. **Run full test suite**
8. **Commit + push**

### Design principle:
- Same interface: `get_sst(day_of_year, node) → float`
- Climatology = 365 daily means (captures real asymmetry, double peaks, etc.)
- Linear warming trend still applied on top (configurable per node)

---

## Phase 4: Validation Run (6:30 PM)
**Goal**: Quick validation with both fixes, compare to baseline.

### Runs:
1. **5K, 5-node, 20yr** with R→S fix + sinusoidal SST (isolate immunity effect)
2. **5K, 5-node, 20yr** with R→S fix + satellite SST (combined effect)
3. Compare to baseline (previous 5K validation: 81.9% crash, 365 recoveries)

### Expected impacts:
- R→S: More reinfection → higher crash rates, fewer net survivors. Recovery still provides reproductive window.
- Satellite SST: Different seasonal patterns may shift epidemic timing, especially at extreme latitudes

### Output:
- `results/validation_rs_fix/` — comparison tables + plots
- Brief analysis markdown

---

## Phase 5: Paper Corrections + Final Compile (8:00 PM)
**Goal**: Fix all three issues in the paper, recompile, push, share.

### Corrections:
1. **Outplanting history** (Introduction, ~line 223):
   - Add: FHL caged outplanting tests 2023 (Seattle Times, OPB)
   - Add: FHL first open release July/Aug 2024, 20 captive-bred stars off SJI dock (KUOW)
   - Keep: Dec 2025 California SSL outplanting in Monterey
   - Willem was a diver in the 2024 release

2. **Immunity description** (Disease module, ~line 431):
   - Remove: "Immune; functionally equivalent to S but not susceptible to reinfection"
   - Replace: R→S transition, explain echinoderms lack adaptive immunity
   - Add: Oregon Coast Aquarium iodine/probiotic protocol as evidence against acquired immunity
   - Update SEIPD+R diagram to show R→S arrow

3. **SST description** (Environmental forcing, ~line 335):
   - Remove: "sinusoidal annual cycle with warming trend"
   - Replace: satellite-derived SST climatology (NOAA OISST v2.1)
   - Describe: per-node extraction, 2002-2025 daily data → 365-day climatology
   - Note: configurable warming trend overlay for projection scenarios

4. **References**: Add KUOW, Mongabay, NOAA OISST citations

5. **Recompile PDF**: `pdflatex` + `bibtex` cycle

6. **Commit + push to GitHub**

7. **Message Willem**: Summary of all changes + link to updated PDF

---

## Timeline Summary

| Phase | Time | Task | Duration |
|-------|------|------|----------|
| 1 | 1:30 PM | R→S immunity fix | ~90 min |
| 2 | 3:30 PM | SST data extraction | ~90 min |
| 3 | 5:00 PM | SST model integration | ~90 min |
| 4 | 6:30 PM | Validation runs | ~90 min |
| 5 | 8:00 PM | Paper corrections + share | ~60 min |

All complete by 9:00 PM.
