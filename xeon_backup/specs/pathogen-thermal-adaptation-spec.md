# Pathogen Thermal Adaptation Spec

## Motivation

The model currently uses a global `T_vbnc` (12°C) for the VBNC sigmoid that gates pathogen activity. This fails to capture a key empirical pattern: the Salish Sea (JDF, SS-S) was hit *harder* than Alaska despite having intermediate SST. JDF/SS summer peaks barely reach 12°C, so the model gives them marginal disease pressure and too much recovery.

Recent evidence suggests V. pectenicida may be one strain of a broader clade, with close relatives potentially adapted to different thermal niches. If pathogen strains can evolve to exploit colder environments (not faster or more virulent — just *infectious* at lower temperatures), this would explain:

- **CA/OR**: Already warm. No adaptation needed. Year-round pressure from arrival.
- **JDF/SS**: Marginal temps, semi-enclosed waterway. Pathogen adapts relatively quickly → year-round pressure → populations stay crushed.
- **BC open coast**: Moderate temps, open water flushes adapted strains. Slower adaptation.
- **Alaska**: Large thermal gap, very cold winters. Adaptation slow/hard → strong seasonal reprieve → recovery.

This makes the model genuinely **co-evolutionary**: host resistance evolves upward while pathogen thermal niche evolves downward. An evolutionary arms race.

## Design

### Per-Node Local VBNC Threshold

Each node tracks a local effective VBNC midpoint `T_vbnc_local` that can decrease over time as cold-adapted strains are selected.

```
T_vbnc_local[i] ∈ [T_vbnc_min, T_vbnc_initial]
```

The VBNC sigmoid at node `i` uses `T_vbnc_local[i]` instead of the global `cfg.T_vbnc`:

```
vbnc_activation = 1 / (1 + exp(-k × (T - T_vbnc_local[i])))
```

### Adaptation Mechanism

Daily, at each disease-reached node with infected hosts:

```python
if T_celsius < T_vbnc_local[i] and (n_I1 + n_I2 + n_E) > 0:
    # Selection pressure: stronger when more hosts available
    # and when temperature is further below threshold
    prevalence = (n_I1 + n_I2) / max(n_total, 1)
    temp_gap = T_vbnc_local[i] - T_celsius  # always positive here
    
    delta_T = adapt_rate × prevalence × temp_gap
    T_vbnc_local[i] = max(T_vbnc_local[i] - delta_T, T_vbnc_min)
```

**Intuition**: When temperature is below the current threshold, pathogen activity is reduced — but the infected hosts still carry viable bacteria. Among the pathogen population within those hosts, variants that can remain active at lower temperatures have a selective advantage. The rate of adaptation depends on:
1. **Prevalence** — more infected hosts = larger pathogen population = more mutation/selection opportunity
2. **Temperature gap** — larger gap = stronger selection pressure (more reproductive advantage for cold-adapted variants)

### Inheritance via Wavefront

When the wavefront activates a new node (cumulative dose exceeds threshold), the arriving pathogen inherits the thermal adaptation of its source:

```python
# At wavefront activation of node j:
# Weighted average of T_vbnc_local from contributing nodes
# weighted by their dispersal contribution to j's cumulative dose
total_dose = 0
weighted_T = 0
for i in reached_nodes:
    dose_from_i = D_wf[i, j] × P[i]  # or use cumulative contribution
    weighted_T += dose_from_i × T_vbnc_local[i]
    total_dose += dose_from_i

if total_dose > 0:
    T_vbnc_local[j] = weighted_T / total_dose
else:
    T_vbnc_local[j] = T_vbnc_initial
```

**Simplification for Chunk 1**: New nodes inherit `T_vbnc_initial`. Inheritance deferred to Chunk 2. The local adaptation will still produce correct qualitative behavior — it just means each node starts adapting from scratch rather than getting a head start.

### Reversion (Optional, Chunk 3)

Without disease pressure, adapted strains may lose fitness advantage and revert:

```python
if (n_I1 + n_I2 + n_E) == 0 and T_vbnc_local[i] < T_vbnc_initial:
    T_vbnc_local[i] = min(T_vbnc_local[i] + revert_rate, T_vbnc_initial)
```

This prevents permanent adaptation at nodes where disease has been cleared. Biologically: without selection pressure, cold-adapted variants (which likely pay a fitness cost in warm conditions) are outcompeted by wild-type.

## Config Parameters

In `DiseaseSection`:

```python
# Pathogen thermal adaptation
pathogen_adaptation: bool = False       # Enable per-node T_vbnc evolution
T_vbnc_initial: float = 12.0           # Starting VBNC threshold (°C)
T_vbnc_min: float = 6.0                # Biophysical floor (°C) — hard limit on cold adaptation
pathogen_adapt_rate: float = 0.001     # Adaptation speed (°C/day per unit selection pressure)
pathogen_revert_rate: float = 0.0005   # Reversion speed when disease absent (°C/day)
```

**Backward compatibility**: `pathogen_adaptation=False` (default) → behavior identical to current. Global `T_vbnc` used everywhere.

## State Changes

### NodeDiseaseState

```python
T_vbnc_local: float = 12.0   # Local VBNC threshold (initialized from T_vbnc_initial)
```

### SpatialSimResult

```python
T_vbnc_local: Optional[np.ndarray] = None  # Final T_vbnc per node (N,)
# Plus optional yearly snapshots for visualization
```

## Code Changes by File

### 1. `config.py` — DiseaseSection
- Add 5 new fields (see Config Parameters above)
- Validation: `T_vbnc_min <= T_vbnc_initial`, `pathogen_adapt_rate >= 0`

### 2. `disease.py` — NodeDiseaseState + adaptation function
- Add `T_vbnc_local` field to `NodeDiseaseState`
- New function `adapt_pathogen_thermal(node_state, T_celsius, n_I1, n_I2, n_E, n_total, cfg)` → updates `T_vbnc_local` in place
- Modify `environmental_vibrio()` to accept optional `T_vbnc_override` parameter (the local value)
- Modify `update_vibrio_concentration()` and `daily_disease_update()` to pass local T_vbnc through

### 3. `model.py` — simulation loop
- Initialize `T_vbnc_local` per node from `T_vbnc_initial` (or inherit at wavefront activation)
- Call `adapt_pathogen_thermal()` daily after disease update
- Pass `T_vbnc_local[i]` to `environmental_vibrio()` calls
- Store final `T_vbnc_local` array in result

### 4. `tests/test_pathogen_adaptation.py`
- Adaptation reduces T_vbnc_local when temp < threshold and infected hosts present
- No adaptation when pathogen_adaptation=False
- No adaptation when temp >= T_vbnc_local
- No adaptation when no infected hosts
- T_vbnc_local never goes below T_vbnc_min
- Reversion when disease absent
- Inheritance at wavefront activation (Chunk 2)

## Implementation Chunks

### Chunk 1: Config + State + Core Adaptation
- Config fields in DiseaseSection
- T_vbnc_local in NodeDiseaseState
- `adapt_pathogen_thermal()` function
- Modify `environmental_vibrio()` for local T_vbnc
- Wire into simulation loop
- Unit tests for adaptation mechanics
- **Backward compatible**: disabled by default

### Chunk 2: Wavefront Inheritance
- Track dose contributions per source node during wavefront accumulation
- At activation, compute weighted-average T_vbnc from sources
- Tests for inheritance

### Chunk 3: Reversion + Recording
- Reversion when disease clears
- Yearly T_vbnc_local snapshots in result
- Monthly recording for visualization
- Integration tests

## Expected Behavior

With adaptation enabled:
- **Year 1-2**: All nodes at T_vbnc_initial=12°C. Disease behavior identical to current.
- **Year 2-4**: Wavefront spreads north. Southern nodes (warm) don't need adaptation. Marginal nodes (JDF/SS, ~11-12°C) start adapting. T_vbnc drops from 12→10°C over ~2 years.
- **Year 4-8**: JDF/SS now have effective year-round disease (T_vbnc≈9-10°C, SST rarely below that). Populations stay suppressed. AK nodes adapt slowly (big thermal gap, lower prevalence), T_vbnc maybe drops to 10-11°C. Winter reprieve persists.
- **Year 8-13**: Strong differential: south=permanent suppression, JDF/SS=adapted pathogen keeps them down, AK=seasonal reprieve allows recovery.

## Calibration Parameters

This adds 3-4 new tunable parameters:
- `T_vbnc_initial` — likely stays at 12°C (Prentice 2025 experiments at 13°C)
- `T_vbnc_min` — biophysical floor, maybe 5-8°C (needs literature check)
- `pathogen_adapt_rate` — THE key new parameter. Controls how fast gradient develops.
- `pathogen_revert_rate` — secondary, controls permanence of adaptation

The adapt_rate controls the *timescale* of the gradient emerging. Too fast → gradient develops in year 1 (unrealistic). Too slow → no differentiation by year 13. Sweet spot: JDF/SS adapted by year 3-4, matching the observed crash timeline.
