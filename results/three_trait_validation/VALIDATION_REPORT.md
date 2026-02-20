# Three-Trait Genetic Architecture — Validation Report

**Date:** 2026-02-20
**Phase:** 6/6 (Final Validation)
**Runtime:** 141.8s
**Test suite:** 632 passing

## Configuration

- **Genetic architecture:** 17R / 17T / 17C = 51 loci
- **Initialization:** beta (target_mean_r = 0.15)
- **Nodes:** 5 (Sitka, Howe Sound, SJI, Newport, Monterey)
- **K:** 5000 per node (25000 total)
- **Duration:** 20 years, disease at year 3
- **Seed:** 42

## 5-Node Results

| Node | Pop₀ | Pop_final | Min (yr) | Crash% | Dis.Deaths | Rec. | Δr | Δt | Δc |
|------|------|-----------|----------|--------|------------|------|----|----|----|
| Sitka | 4932 | 310 | 310 (yr19) | 93.8% | 10794 | 215 | -0.0151 | -0.0267 | +0.0616 |
| Howe Sound | 4938 | 180 | 180 (yr19) | 96.4% | 9813 | 143 | +0.0080 | +0.0134 | +0.0488 |
| SJI | 4944 | 201 | 201 (yr19) | 96.0% | 11583 | 195 | -0.0105 | +0.0007 | +0.0434 |
| Newport | 5000 | 333 | 199 (yr16) | 96.0% | 8659 | 178 | +0.0493 | +0.0026 | +0.0479 |
| Monterey | 4999 | 3495 | 1504 (yr8) | 69.9% | 38890 | 1274 | +0.0093 | +0.0311 | +0.0975 |

## Comparison to Single-Trait Baseline (Phase 14)

Previous run: 51 resistance loci, same 5 nodes, K=5000, seed=42, 20yr.

| Node | Crash (old) | Crash (new) | Δ | Δr (old) | Δr (new) | Δt (new) | Δc (new) |
|------|-------------|-------------|---|----------|----------|----------|----------|
| Sitka | 95.4% | 93.8% | -1.6pp | +0.031 | -0.0151 | -0.0267 | +0.0616 |
| Howe Sound | 93.2% | 96.4% | +3.2pp | +0.044 | +0.0080 | +0.0134 | +0.0488 |
| SJI | 99.3% | 96.0% | -3.3pp | +0.061 | -0.0105 | +0.0007 | +0.0434 |
| Newport | 99.1% | 96.0% | -3.1pp | +0.055 | +0.0493 | +0.0026 | +0.0479 |
| Monterey | 98.8% | 69.9% | -28.9pp | +0.051 | +0.0093 | +0.0311 | +0.0975 |

### Key Differences

With three traits (17 loci each) vs one trait (51 loci):
- Each trait has fewer loci → less genetic variance per trait
- Resistance (r_i) alone is weaker → tolerance and recovery provide alternative survival pathways
- Selection acts on all three traits simultaneously

## Emergent Dynamics

### 1. Recovery (c_i) is the Fastest Evolving Trait

- Mean Δr across nodes: +0.0082
- Mean Δt across nodes: +0.0042
- **Mean Δc across nodes: +0.0598** ← 7× faster than resistance, 14× faster than tolerance

c_i is the fastest-evolving trait at **every single node**. This makes biological sense: recovery acts daily on every I₁/I₂ individual (p_rec = ρ_rec × c_i), providing strong per-generation selection. Resistance only matters at the moment of exposure, and tolerance merely delays death without clearing infection.

### 2. Monterey Anomaly: Partial Recovery Instead of Collapse

Monterey dropped from 98.8% crash (single-trait) to **69.9% crash** (three-trait) — a 28.9pp improvement. This is the most dramatic difference from the baseline and suggests:
- Recovery pathway (c_i) enables some individuals to clear infection and survive
- 1,274 recovered individuals at year 20 vs near-total wipeout in single-trait
- Warm SST drives sustained epidemic pressure, which selects strongly for c_i (+0.098)
- Monterey's larger surviving population (3,495) provides more genetic material for selection

### 3. Silent Spreaders — Minimal Signal

Only 5 late-stage infected (I₂) agents remain at year 20 (2 at Newport, 3 at Monterey). Mean timers of 6–7 days and moderate t_i values (0.07–0.12) suggest tolerance is extending I₂ duration slightly but not creating a large reservoir of silent spreaders. This may change with higher initial t_i or pathogen evolution.

### 4. Trait Correlations Are Weak

| Node | r(r,t) | r(r,c) | r(t,c) |
|------|--------|--------|--------|
| Sitka (n=310) | +0.172 | -0.068 | -0.161 |
| Howe Sound (n=180) | +0.097 | -0.058 | -0.125 |
| SJI (n=201) | -0.100 | -0.158 | +0.078 |
| Newport (n=333) | -0.174 | -0.162 | +0.091 |
| Monterey (n=3495) | +0.002 | +0.033 | +0.061 |

Correlations are weak (|r| < 0.18), confirming that the three traits evolve largely independently. The strongest signal: slight negative r(r,c) at bottlenecked northern nodes — possible trade-off under extreme selection (survivors tend to be high-c OR high-r, not both). Monterey's large population shows near-zero correlations, consistent with independent assortment.

### 5. Resistance Signal Is Inconsistent

Unlike the single-trait baseline (where all nodes showed Δr > 0), two northern nodes (Sitka, SJI) show **negative** Δr. With only 17 resistance loci, genetic drift in small post-crash populations can overpower selection — especially when recovery provides an alternative survival route.

### 6. All Nodes Have Recovered Individuals

Every node has R-state individuals at year 20 (143–1,274). This never happened in the single-trait model because recovery was entirely driven by the same resistance score. With a dedicated c_i trait, recovery is now a distinct evolutionary pathway.

## Figures

1. `population_trajectories.png` — Population over time per node
2. `trait_evolution.png` — Three-panel trait mean trajectories
3. `trait_comparison.png` — Δr vs Δt vs Δc bar chart per node

## Concerns & Notes

1. **Sitka and SJI show negative Δr** — drift vs selection at small N. Need convergence testing (multiple seeds) to assess reliability.
2. **Monterey recovery is unexpected** — 69.9% crash vs 98.8% baseline. Could indicate that c_i is too powerful, or that this is genuine emergent behavior from multi-trait selection. Needs sensitivity analysis.
3. **Runtime increased** from ~51s (single-trait) to **141.8s** (three-trait) — largely due to 5× K (5000 vs 1000 default). The three-trait computation itself is minimal overhead.
4. **Tolerance shows weakest selection** — Δt ≈ 0 at most nodes. Timer-scaling mechanism may need recalibration, or tolerance may genuinely be weakly selected compared to clearance.

## Summary

| Metric | Value |
|--------|-------|
| Total initial population | 25,000 |
| Total final population | 4,519 |
| Overall crash | 81.9% |
| Total disease deaths | 79,739 |
| Total recovered (alive) | 2,005 |
| Fastest evolving trait | Recovery (c_i), Δ = +0.060 |
| Runtime | 141.8s |
| Tests passing | 632 |
