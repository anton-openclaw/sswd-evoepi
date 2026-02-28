# Wavefront Disease Spread — Design Specification

## Problem Statement

The current model has two disease scenarios:

1. **Ubiquitous** (default): P_env is active at ALL nodes from day 1. Disease seeds at all nodes simultaneously on `disease_year`. There is no spatial spread — every node gets sick at once.

2. **Invasion**: P_env returns 0 everywhere. Disease seeds at specified `invasion_nodes` with a dose. But pathogen dispersal via D matrix only updates `vibrio_concentration` for tracking — it does NOT feed into infection probability. New infections are driven purely by local shedding (P_shed from I₁/I₂ at the node) + P_env. So disease can only spread between nodes via **infected individuals physically moving there** (CRW movement).

**The real SSWD wavefront**: First detected at Channel Islands / Santa Cruz ~June 2013, spread north through OR, WA, BC, reaching Alaska by 2014-2015. A ~2-year south-to-north wave. This requires both individual movement AND waterborne pathogen dispersal between nodes.

## Current Code Analysis

### Where infection happens (`disease.py`, `daily_disease_update`)

```python
# STEP 1: Update Vibrio concentration
P_k = update_vibrio_concentration(
    P_k, n_I1, n_I2, n_D_fresh,
    T_celsius, salinity, phi_k, dispersal_input,  # dispersal_input DOES affect P_k
    cfg, ...
)

# STEP 2: Force of infection uses P_k
dose_response = P_k / (cfg.K_half + P_k)
lambda_arr = cfg.a_exposure * dose_response * r_mod * sal_mod * s_mod
```

**Key finding**: `dispersal_input` (from D matrix) DOES feed into `P_k` via `update_vibrio_concentration()`. So waterborne pathogen from neighboring nodes DOES increase local `P_k`, which DOES increase infection probability. **This pathway already works.**

### Where P_env enters (`disease.py`, `update_vibrio_concentration`)

```python
env = environmental_vibrio(T_celsius, salinity, cfg)
dP = shed - decay - flush + env + dispersal_input
```

`environmental_vibrio()` returns `P_env_max * vbnc_activation * g_peak * sal_mod` — a constant background input based on local temperature. This is the "VBNC reservoir" — vibrio in sediment/biofilms that resuscitate when warm enough.

### Where disease is seeded (`model.py`, `run_spatial_simulation`)

```python
# Section 6: "Seed disease if it's the epidemic year"
if disease_year is not None and year == disease_year:
    for i, node in enumerate(network.nodes):
        disease_active_flags[i] = True  # ALL nodes activated at once
        # ... seed initial_infected_per_node at each node
```

### Where P_env initialization happens (`model.py`, ~line 1959)

```python
# At population init, each node gets steady-state P_env:
if dis_cfg.scenario == "ubiquitous":
    env = environmental_vibrio(nd.mean_sst, nd.salinity, dis_cfg)
    xi = vibrio_decay_rate(nd.mean_sst)
    phi = nd.flushing_rate
    node.vibrio_concentration = env / (xi + phi) if (xi + phi) > 0 else 0.0
```

This means in ubiquitous mode, every node has ambient vibrio from the start.

## What Needs to Change

### 1. Disease origin seeding — INSTEAD of all nodes, seed at CA-S only

Currently: `for i, node in enumerate(network.nodes): disease_active_flags[i] = True`

Need: Only seed at specified origin node(s). Other nodes have `disease_active_flags[i] = False`.

### 2. Per-node `disease_reached` tracking

New boolean per node. Flips `True` when:
- Explicitly seeded (origin node on disease_year), OR
- First local infection occurs at that node (from any pathway)

### 3. Gate P_env behind `disease_reached`

Currently `environmental_vibrio()` returns P_env based on temperature everywhere. Need to make it return 0 at nodes where disease hasn't reached yet.

**Rationale**: P_env represents VBNC reservoir — vibrio persisting in local sediment/biofilms. There IS no local reservoir until the pathogen arrives and establishes. Before the wavefront reaches a node, local vibrio = 0.

### 4. Activate disease_active_flags when disease reaches each node

Currently all flags flip on `disease_year`. Need: only origin node(s) flip on `disease_year`. Other nodes activate when they first get pathogen (via D matrix dispersal input raising vibrio above threshold, or via an infected individual arriving via movement).

### 5. Config changes

Add to `DiseaseSection`:
- `disease_origin_nodes: Optional[List[int]]` — node IDs where disease starts. If None, all nodes (backward compat).
- `wavefront_enabled: bool = False` — if True, enables gated P_env and per-node activation. Default False for backward compatibility.
- `activation_threshold: float = 1.0` — vibrio concentration (bact/mL) above which a node's disease becomes "reached" (even without explicit infection). Handles waterborne arrival.

### 6. Track disease arrival time per node

`disease_arrival_day: np.ndarray` shape (N,), initialized to -1. Set to sim_day when disease_reached flips True. Useful diagnostic output.

## Backward Compatibility

- Default: `wavefront_enabled = False`, `disease_origin_nodes = None` → identical to current behavior (all nodes seeded at disease_year, P_env everywhere).
- When `wavefront_enabled = True`: requires `disease_origin_nodes` to be set. P_env gated. Disease spreads as wavefront.

## Implementation Chunks

### Chunk 1: Config + NodeDiseaseState changes
- Add `wavefront_enabled`, `disease_origin_nodes`, `activation_threshold` to `DiseaseSection`
- Add `disease_reached` field to `NodeDiseaseState`
- Add validation in `validate_config()`
- **No behavioral changes yet** — purely additive

### Chunk 2: Gate P_env behind disease_reached
- Modify `daily_disease_update()` to accept `disease_reached: bool` parameter
- When `disease_reached=False`: skip `environmental_vibrio()` call, use 0 for env input
- When `disease_reached=True`: existing behavior
- Default `disease_reached=True` for backward compat
- Modify `update_vibrio_concentration()` or add gating in the caller
- Update the steady-state P_env initialization in model.py (init to 0 when wavefront enabled)

### Chunk 3: Wavefront seeding + activation logic in model.py
- Modify `run_spatial_simulation()`:
  - Disease seeding: only origin nodes get initial infections + `disease_active_flags[i] = True`
  - Non-origin nodes: `disease_active_flags[i] = False` until first pathogen arrival
  - Daily check: if `vibrio_concentration > activation_threshold` or first infection occurs → set `disease_reached[i] = True`, `disease_active_flags[i] = True`
  - Track `disease_arrival_day[i]`
  - Pass `disease_reached[i]` to `daily_disease_update()`
- Modify `run_coupled_simulation()` similarly (single-node: trivial, just add param)
- Add `disease_arrival_day` to `SpatialSimResult`

### Chunk 4: Tests
- Test backward compat: `wavefront_enabled=False` → identical behavior
- Test wavefront seeding: only origin nodes get disease initially
- Test P_env gating: unreached nodes have 0 background vibrio
- Test waterborne spread: pathogen dispersal triggers activation at neighbor node
- Test movement spread: infected individual moves to clean node → activates it
- Test disease_arrival_day tracking
- Test activation_threshold: node activates when vibrio exceeds threshold
- Small multi-node validation: 5-node chain, seed disease at one end, verify wavefront propagation

### Chunk 5: Verification run
- 5-node linear chain (or reuse 11-node SA chain), K=1000
- Seed disease at node 0 (southernmost)
- Verify wavefront spreads north over ~months
- Compare with all-nodes-at-once to confirm behavioral difference
- Output: disease_arrival_day per node, population trajectories
