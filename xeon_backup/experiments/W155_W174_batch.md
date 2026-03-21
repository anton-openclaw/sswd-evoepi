# Calibration Batch W155-W174: Non-Uniform K + Fine-Tuning

## Rationale
Building on W154 (current best: mean RMSE_log=0.479, n_conn=0.3, α_env=0.18).
Introducing K_cv > 0 for non-uniform carrying capacity (lognormal distribution).

### Why K_cv matters
- With uniform K=5000, all sites are identical in population size
- In reality, Pycnopodia density varies enormously across sites:
  - Productive kelp forests in AK: high density sites
  - Marginal southern habitat: sparse populations
  - Gravem et al. 2021: pre-SSWD densities varied by orders of magnitude
- K_cv creates lognormal variation: E[K] = K_mean, some sites much larger/smaller
- K_cv=0.5 → K roughly 2,000-15,000; K_cv=1.0 → K roughly 500-50,000+
- High-K sites: larger populations → more survivors → better recovery potential
- Low-K sites: Allee effects stronger → harder recovery → amplifies extinction

### W154 baseline config
```
K=5000, K_cv=0.0
disease.K_half=800000
disease.alpha_env=0.18
disease.delta_env=0.02
spatial.n_connectivity=0.3
spatial.alpha_self_open=0.02
spatial.alpha_self_fjord=0.70
disease.pathogen_adaptation=true
disease.virulence_evolution=true
disease.T_vbnc_initial=12.0
disease.T_vbnc_min=9.0
disease.k_vbnc=2.0
disease.P_env_dynamic=true
disease.P_env_floor=500.0
disease.wavefront_enabled=true
+ all other W154 params
```

### W154 regional results (seed 123, RMSE_log=0.479)
- AK-PWS: 8.0% (target 50%) ← main problem
- AK-FN: 12.2% (target 50%) ← main problem
- AK-FS: 10.6% (target 20%) ← low
- BC-N: 5.1% (target 20%) ← low
- SS-S: 8.2% (target 5%) ← slightly high but OK
- JDF: 3.8% (target 2%) ← OK
- OR: 0.8% (target 0.25%) ← high
- CA-N: 0.1% (target 0.1%) ← perfect
- CA-S: 0.0% (target ~0%) ← perfect

## Batch Design

### Tier 1: K_cv sweep (W155-W162) — 8 runs
Fix all W154 params, vary K_cv only.
This isolates the effect of population heterogeneity.

| Run  | K_cv | K_half | α_env | n_conn | Notes |
|------|------|--------|-------|--------|-------|
| W155 | 0.3  | 800K   | 0.18  | 0.3    | Mild variation |
| W156 | 0.5  | 800K   | 0.18  | 0.3    | Moderate variation |
| W157 | 0.8  | 800K   | 0.18  | 0.3    | Strong variation |
| W158 | 1.0  | 800K   | 0.18  | 0.3    | Very strong variation |
| W159 | 0.3  | 1.2M   | 0.18  | 0.3    | Mild + higher K_half |
| W160 | 0.5  | 1.2M   | 0.18  | 0.3    | Moderate + higher K_half |
| W161 | 0.8  | 1.2M   | 0.18  | 0.3    | Strong + higher K_half |
| W162 | 1.0  | 1.2M   | 0.18  | 0.3    | Very strong + higher K_half |

### Tier 2: α_env tuning with best K_cv (W163-W168) — 6 runs
Take the best K_cv from Tier 1, cross with α_env and K_half.

| Run  | K_cv   | K_half | α_env | δ_env | n_conn | Notes |
|------|--------|--------|-------|-------|--------|-------|
| W163 | best   | 800K   | 0.22  | 0.02  | 0.3    | Stronger env accumulation |
| W164 | best   | 800K   | 0.25  | 0.02  | 0.3    | Even stronger |
| W165 | best   | 1.2M   | 0.22  | 0.02  | 0.3    | K_half lift + α |
| W166 | best   | 1.2M   | 0.25  | 0.02  | 0.3    | K_half lift + stronger α |
| W167 | best   | 1.2M   | 0.18  | 0.01  | 0.3    | Slower pathogen decay (t½=69d vs 35d) |
| W168 | best   | 1.5M   | 0.22  | 0.02  | 0.3    | Much higher K_half |

**Literature-justified K_cv range: 0.5-1.0** (marine benthic invertebrate densities follow lognormal; pre-SSWD Pycnopodia varied 1-2 OOM across range).

**Generator script:** `experiments/calibration/generate_W163_W168.sh <best_cv>`
Creates all 6 configs + launcher in one command.

### Tier 3: Exploratory (W169-W174) — 6 runs
Based on Tier 1-2 results. Reserved for:
- delta_env tuning (0.02 vs 0.01 vs 0.03)
- n_conn fine-tuning around 0.3
- Combined best settings
- Any ideas from reviewing Tier 1-2

## Execution Plan

### Timing
- Each run: ~5-7 hours (K=5000, 896 nodes, 13 years)
- 3 seeds each: 42, 123, 999
- Tier 1: 8 configs × 3 seeds = 24 runs → 1 wave at 24 parallel → ~7 hours
- Tier 2: 6 configs × 3 seeds = 18 runs → 1 wave → ~7 hours
- Tier 3: TBD after Tier 1-2

### Parallel capacity
- Xeon: 64 cores, 503GB RAM
- Sweet spot: 24 parallel (~20GB each = 480GB)
- DO NOT exceed 24

### Thread safety
```bash
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1
export NUMBA_NUM_THREADS=1
export NUMBA_THREADING_LAYER=omp
export PYTHONUNBUFFERED=1
```

## Monitoring
- Cron watcher job every 2 hours
- Check process count, memory usage
- Alert on completion or failure
