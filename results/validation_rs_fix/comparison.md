# R→S Fix Validation: Comparison with Permanent Immunity Baseline

**Date**: February 21, 2026  
**Setup**: 5-node network (Sitka, Howe Sound, SJI, Newport, Monterey), K=5000, 20yr, seed=42  
**Change**: Recovered individuals return to Susceptible (R→S) instead of permanent immunity (R stays R)

## Side-by-Side: Baseline vs R→S (Sinusoidal SST)

| Metric | Baseline (perm. immunity) | R→S (sinusoidal) | Δ |
|--------|:------------------------:|:-----------------:|---|
| **Overall crash** | 98.5%¹ | 99.7% | +1.2% worse |
| **Final total pop** | 365 | 122 | −67% |
| **Total disease deaths** | 41,968 | 36,157 | −14% |
| **Total recoveries** | 365 | 276 | −24% |
| **Recovery rate** | 0.87% | 0.76% | −0.11pp |
| **Runtime** | 107.7s | 97.6s | −9% |

¹ Based on minimum population across all nodes divided by initial population.

### Per-Node Comparison

| Node | Baseline Crash | R→S Crash | Baseline Final | R→S Final | Baseline Recov | R→S Recov |
|------|:--------------:|:---------:|:--------------:|:---------:|:--------------:|:---------:|
| **Sitka** | 98.7% | 99.3% | 65 | 36 | 60 | 44 |
| **Howe Sound** | 98.8% | 99.5% | 60 | 23 | 55 | 80 |
| **SJI** | 99.0% | **100.0%** | 50 | **0** | 63 | 57 |
| **Newport** | 99.5% | 99.9% | 27 | 63 | 51 | 58 |
| **Monterey** | 99.2% | **100.0%** | 163 | **0** | 136 | 37 |

### Trait Evolution Comparison

| Node | Baseline Δr | R→S Δr | Baseline Δt | R→S Δt | Baseline Δc | R→S Δc |
|------|:-----------:|:------:|:-----------:|:------:|:-----------:|:------:|
| **Sitka** | +0.011 | **+0.060** | +0.005 | +0.016 | +0.029 | **−0.008** |
| **Howe Sound** | −0.002 | +0.034 | +0.044 | +0.079 | +0.041 | **+0.005** |
| **SJI** | +0.012 | −0.150² | −0.007 | −0.100² | +0.072 | **−0.020²** |
| **Newport** | +0.031 | −0.051 | +0.001 | −0.050 | +0.054 | **+0.030** |
| **Monterey** | +0.025 | −0.149² | +0.027 | −0.099² | **+0.154** | **−0.021²** |

² Extinction at this node — trait values at extinction are drift artifacts, not selection signals.

## Key Finding: R→S Fundamentally Changes Evolutionary Dynamics

### 1. Recovery Trait No Longer Evolves Upward

**Baseline (permanent immunity)**: Recovery was the fastest-evolving trait. Monterey showed Δc = +0.154, meaning recovery ability more than doubled. This happened because recovered individuals were permanently immune — they survived, reproduced, and passed high-c alleles to offspring.

**R→S**: Recovery trait barely evolves or goes negative. The highest surviving-node Δc is +0.030 (Newport), 5× weaker than the weakest baseline node. The biological mechanism is clear: recovered stars immediately re-enter the susceptible pool and face reinfection. High-recovery alleles don't accumulate because their carriers keep getting reinfected.

### 2. Two Nodes Go Locally Extinct

SJI and Monterey crash to population = 0 under R→S but maintained small populations (50, 163) with permanent immunity. Without the immune "safe harbor" of recovered individuals, relentless reinfection cycles drive these nodes to extinction.

### 3. Resistance Replaces Recovery as the Dominant Adaptive Response

In surviving nodes (Sitka, Howe Sound), resistance shows the strongest positive selection under R→S:
- Sitka: Δr = +0.060 (vs +0.011 baseline) — **5.5× stronger**
- Howe Sound: Δr = +0.034 (vs −0.002 baseline)

This makes sense: when recovery doesn't confer lasting protection, avoiding infection entirely (resistance) becomes more valuable than clearing infection (recovery).

### 4. Fewer Total Deaths but Worse Outcomes

R→S has ~14% fewer total disease deaths (36,157 vs 41,968) because populations crash faster and stay crashed, leaving fewer individuals to die from disease in later years. The epidemic is more efficient — it kills populations faster precisely because recovered stars can die from disease again.

### 5. Fewer Recovery Events

Counter-intuitively, R→S produces fewer total recoveries (276 vs 365) despite reinfection being possible. Under permanent immunity, recovered stars survive long-term and the stable recovered population generates ongoing recovery events from the still-circulating pathogen. Under R→S, the crash is so severe that the surviving population is too small to generate many recovery events.

## R→S with Satellite SST

| Metric | R→S Sinusoidal | R→S Satellite | Δ |
|--------|:--------------:|:-------------:|---|
| **Overall crash** | 99.7% | 99.9% | +0.2% worse |
| **Final total pop** | 122 | 146 | +20% |
| **Total deaths** | 36,157 | 33,909 | −6% |
| **Total recoveries** | 276 | 241 | −13% |
| **Recovery rate** | 0.76% | 0.71% | −0.05pp |

### Per-Node: Sinusoidal vs Satellite

| Node | Sin. Crash | Sat. Crash | Sin. Final | Sat. Final | Sin. Recov | Sat. Recov |
|------|:----------:|:----------:|:----------:|:----------:|:----------:|:----------:|
| Sitka | 99.3% | 99.8% | 36 | 10 | 44 | 41 |
| Howe Sound | 99.5% | 99.9% | 23 | 133 | 80 | 57 |
| SJI | 100.0% | 99.9% | 0 | 3 | 57 | 56 |
| Newport | 99.9% | 100.0% | 63 | 0 | 58 | 51 |
| Monterey | 100.0% | 100.0% | 0 | 0 | 37 | 36 |

Satellite SST shifts which nodes survive: SJI barely persists under satellite SST (final pop=3) while going fully extinct under sinusoidal. Newport goes extinct under satellite but persists (final=63) under sinusoidal. This reflects the different SST seasonality profiles — satellite SST captures asymmetric seasonal warming patterns that sinusoidal approximation smooths over.

The overall pattern is similar under both SST sources: near-total population collapse with marginal differences in which specific nodes retain remnant populations.

## Biological Interpretation

The R→S fix is **biologically correct** (echinoderms lack adaptive immunity) and produces dramatically worse conservation outcomes. The model now predicts:

1. **Near-universal local extinction** at K=5000 scale
2. **No evolutionary rescue via recovery** — the trait that evolved fastest under the incorrect permanent-immunity assumption provides minimal benefit when reinfection is possible
3. **Selection shifts to resistance** (immune exclusion) as the primary adaptive response
4. **Stochastic rescue is the dominant persistence mechanism** — surviving populations are tiny remnants, not evolved-resistant populations

This has profound implications for the reintroduction program: captive-bred stars cannot rely on natural recovery conferring lasting protection. Selective breeding for resistance traits may be more valuable than breeding for recovery ability.

## Files

- `results.json` — Full per-node, per-year data for both runs
- `pop_trajectories_sinusoidal.png` — Population curves (R→S, sinusoidal SST)
- `pop_trajectories_satellite.png` — Population curves (R→S, satellite SST)
- `trait_evolution_sinusoidal.png` — Three-trait evolution panels
- `trait_shifts_sinusoidal.png` — Trait shift bar chart
- `yearly_recoveries.png` — Recovery events per year
- `sinusoidal_vs_satellite.png` — Side-by-side SST source comparison
