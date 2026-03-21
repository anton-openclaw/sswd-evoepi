# Sentinel Agent Design — Multi-Host Reservoir

## Concept
Model non-Pycnopodia sea star species (Pisaster, Dermasterias, Evasterias, etc.)
as immortal, non-reproducing "sentinel" agents that carry and shed SSWD pathogen
at a low rate. This creates a persistent multi-host disease reservoir that
Pycnopodia must contend with regardless of its own population dynamics.

## Biological Justification
- Multiple asteroid species get SSWD but with varying severity
- Some species are more tolerant — they survive infection and become carriers
- These carriers maintain environmental Vibrio even when Pyc populations crash
- Critical for reintroduction: released Pyc face pathogen from OTHER species

## Implementation Plan

### 1. Agent-Level Changes (`sswd_evoepi/types.py`)
Add `is_sentinel` field to AGENT_DTYPE:
```python
('is_sentinel', np.int8),  # 1 = sentinel (non-Pyc species), 0 = normal
```

### 2. Model Changes (`sswd_evoepi/model.py`)

**Mortality**: Sentinels skip natural mortality AND disease mortality
```python
if agents['is_sentinel'][slot]:
    continue  # immortal
```

**Reproduction**: Sentinels excluded from spawning/fertilization
```python
eligible = alive & ~agents['is_sentinel'].astype(bool)
```

**Disease transmission**: Sentinels CAN get infected (same FOI as Pyc)
- They transition S → E → I1 → I2 normally
- But they NEVER die from disease (effectively infinite tolerance timer)
- They stay in I2 indefinitely, shedding at a reduced rate

**Shedding**: Sentinel shedding rate = `sentinel_shedding_fraction` × normal
```python
shedding = sigma * n_infected_normal + sigma * sentinel_shedding_fraction * n_infected_sentinel
```

**Recovery gating**: Sentinels can recover (I → S) via the same P_env-gated
mechanism, but then can get re-infected. This creates an endemic cycling.

### 3. Config Changes (`sswd_evoepi/config.py`)
```python
@dataclass
class SentinelSection:
    enabled: bool = False
    n_per_site: int = 20                    # Sentinel individuals per eligible site
    shedding_fraction: float = 0.1          # Shed at 10% of Pyc rate
    site_filter: str = 'rocky'              # 'rocky', 'all', or 'custom'
    custom_node_ids: List[int] = field(default_factory=list)
```

### 4. Initialization (`sswd_evoepi/model.py` or `calibration_runner.py`)
At simulation start, for each eligible site:
- Create `n_per_site` sentinel agents
- Set `is_sentinel = 1`, `alive = 1`, `sex = 0` (arbitrary)
- Set disease_state = S (susceptible — they'll get infected naturally)
- Set age = 10000 (old, doesn't matter)
- Do NOT count toward initial K or density calculations

### 5. Metrics Exclusion
**Population counts**: Filter `alive & ~is_sentinel`
**Recovery fractions**: Exclude sentinels
**Genetic tracking**: Sentinels have no meaningful genetics (or fixed at 0)
**Monthly NPZ**: Record sentinel counts separately

### Site Eligibility (from habitat data)
- **758/896 sites (84.6%)** are rocky/kelp/reef → get sentinels
- **138 sites** are non-rocky (fjords, estuarine, sand, glacial) → no sentinels
- Based on `all_sites.json` habitat field containing 'rock', 'kelp', or 'reef'

### Parameter Space for Experiments
| Parameter | Values to test |
|-----------|---------------|
| n_per_site | 10, 20, 50 |
| shedding_fraction | 0.05, 0.10, 0.25 |
| Sentinels ON vs OFF | Baseline comparison |

### Key Predictions
1. Sentinels should INCREASE the difficulty of reintroduction (higher ambient P_env)
2. Effect should be stronger in warm water (sentinels stay infected longer in warm)
3. Cold-water sites might still clear P_env despite sentinels (δ_env dominates)
4. At high shedding_fraction, even immune Pyc might struggle (P_env floor too high)

### Files to Modify
1. `sswd_evoepi/types.py` — add is_sentinel to AGENT_DTYPE
2. `sswd_evoepi/model.py` — sentinel logic in mortality, reproduction, disease, shedding
3. `sswd_evoepi/config.py` — SentinelSection dataclass
4. `experiments/calibration_runner.py` — sentinel initialization
5. `sswd_evoepi/results.py` — exclude sentinels from metrics
6. Tests — sentinel-specific unit tests
