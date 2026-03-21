# Pathogen Community-Level Thermal Adaptation — Redesign Spec

## Motivation

The v1 pathogen thermal adaptation mechanism (W95-W104) failed because it was gated
by infected host prevalence. After the SSWD crash kills most hosts, prevalence → 0
and adaptation stalls. Over 13 years, T_vbnc barely moved (max Δ=0.14°C).

Willem's insight: bacteria are short-lived individually but the **site's community**
persists. Adaptation gained through host infection cycles should be **sticky** — the
community remembers. New bacteria at the site inherit the local adaptation. The host
infection cycle is the *engine* of adaptation; the site community is the *memory*.

## Design

### 1. Host-Driven Adaptation (KEEP, modify driver)

The existing mechanism stays as the primary driver of adaptation. When hosts are
infected, pathogen cycles through them and experiences selection from local temperature.
This is how cold-adaptation *enters* the system.

**Change**: Replace `prevalence` (infected/total_hosts) with a term based on 
**P_env_pool** — the bacterial community size. As long as vibrio are present at the
site, the community is turning over and experiencing selection.

```python
# OLD (v1) — dies with host crash:
delta_T = adapt_rate * prevalence * temp_gap

# NEW (v2) — persists via bacterial community:
# Pool-driven: community turnover adapts to local conditions
pool_factor = min(P_env_pool / P_adapt_half, 1.0)  # Michaelis-Menten saturation
delta_T = adapt_rate * pool_factor * temp_gap
```

`P_adapt_half`: half-saturation constant for adaptation rate (new config param).
When P_env_pool >> P_adapt_half, adaptation runs at full rate. When pool is thin,
adaptation slows. This captures: more bacteria = more generations = faster community shift.

### 2. Kill Reversion (or make it very slow)

**Change**: Set `pathogen_revert_rate` default to 0.0 (was 0.001).

The adapted community doesn't snap back to naive when hosts die. The bacteria are
still there, still cold-adapted. Reversion should only happen if the site's bacterial
community is actually replaced by immigration from warm-adapted sources — which is
captured by the dispersal mixing (section 3).

Keep the config parameter for flexibility but default to zero.

### 3. Dispersal Carries Thermal Information

**NEW**: When bacteria disperse between nodes, receiving nodes' T_vbnc shifts toward
the weighted mean of incoming bacteria.

After the daily pathogen dispersal step (`D_T_sparse @ P`), update T_vbnc:

```python
# For each node i with P_env_pool > 0:
P_local = node_disease_states[i].vibrio_concentration  # before dispersal
T_local = node_disease_states[i].T_vbnc_local

# Dispersal-weighted incoming T_vbnc
T_incoming_weighted = 0.0
P_incoming_total = 0.0
for each source j contributing to node i:
    P_in_j = D[i,j] * P[j]  # amount of bacteria from source j
    T_in_j = node_disease_states[j].T_vbnc_local
    T_incoming_weighted += P_in_j * T_in_j
    P_incoming_total += P_in_j

# Mix: weighted average of local + incoming
if P_local + P_incoming_total > 0:
    node_disease_states[i].T_vbnc_local = (
        P_local * T_local + T_incoming_weighted
    ) / (P_local + P_incoming_total)
```

**Efficient implementation**: Use sparse D matrix. Build vector of T_vbnc values,
compute `D_T_sparse @ (P * T_vbnc)` to get weighted incoming T, divide by
`D_T_sparse @ P` (which we already compute as `dispersal_in`).

```python
# Vectorized:
T_vec = np.array([nds.T_vbnc_local for nds in node_disease_states])
PT_vec = P * T_vec  # concentration-weighted T_vbnc per node
incoming_PT = D_T_sparse @ PT_vec  # weighted T_vbnc arriving at each node
incoming_P = dispersal_in  # already computed

# For each reached node with nonzero total:
P_before = P.copy()  # vibrio BEFORE adding dispersal
total = P_before + incoming_P
mask = total > 0
T_vec[mask] = (P_before[mask] * T_vec[mask] + incoming_PT[mask]) / total[mask]

# Write back:
for i in range(N):
    if total[i] > 0:
        node_disease_states[i].T_vbnc_local = T_vec[i]
```

### 4. Site Bias Only After Population

T_vbnc_local stays at `T_vbnc_initial` until disease reaches the site (wavefront
activation). At activation, inherit T_vbnc from source nodes (existing `_inherit_T_vbnc`
function). This is already implemented and stays unchanged.

### 5. Floor Bacteria Inherit Site Bias

The P_env_floor (SST-modulated community floor) represents the background vibrio
community at the site. Under this model, the floor community IS locally adapted.
No code change needed — the floor adds to P_env_pool concentration, and T_vbnc_local
already applies to all bacteria at the site regardless of source.

## Config Changes

```python
# In DiseaseSection:
P_adapt_half: float = 500.0      # NEW: half-saturation for pool-driven adaptation
pathogen_revert_rate: float = 0.0  # CHANGED default: 0.001 → 0.0
```

No new config fields beyond P_adapt_half. The `pathogen_adapt_rate` and other
existing fields stay.

## Implementation Chunks

### Chunk 1: Modify adapt_pathogen_thermal() 
- Replace prevalence-based driver with P_env_pool-based driver
- Add `P_adapt_half` config field
- Change `pathogen_revert_rate` default to 0.0
- Update function signature to accept P_env_pool instead of host counts
- Update tests

### Chunk 2: Add dispersal T_vbnc mixing
- After pathogen dispersal step in model.py, compute T_vbnc mixing
- Vectorized implementation using sparse D matrix
- New tests for dispersal mixing

## Expected Behavior

- **During active infection**: adaptation driven by both host cycling AND pool turnover
- **After crash**: adaptation continues (pool persists via floor + dynamic P_env)
- **Spatial spread**: cold-adapted bacteria disperse to neighbors, gradually shifting
  receiving sites' T_vbnc
- **South stays warm**: warm sites have T_vbnc ≈ SST, no gap to adapt into
- **North gets colder**: sustained adaptation over 13 years, NOT gated by host survival
- **JDF/SS-S**: moderate pool, moderate gap — key test of whether this resolves the
  structural misfit

## Files Modified

- `sswd_evoepi/config.py` — add P_adapt_half, change revert default
- `sswd_evoepi/disease.py` — modify adapt_pathogen_thermal() signature and logic
- `sswd_evoepi/model.py` — add T_vbnc dispersal mixing after pathogen dispersal step,
  update adapt_pathogen_thermal() call sites
- `tests/test_pathogen_adaptation.py` — update existing tests, add new ones
