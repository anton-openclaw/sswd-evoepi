# Three-Trait Genetic Architecture — Validation Report

**Date:** 2026-02-20
**Phase:** 6/6 (Final Validation)
**Runtime:** 116.0s
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
| Sitka | 4932 | 135 | 61 (yr12) | 98.8% | 8238 | 6 | -0.0637 | -0.0046 | +0.0458 |
| Howe Sound | 4938 | 28 | 28 (yr19) | 99.4% | 7530 | 4 | +0.0583 | -0.0236 | +0.0202 |
| SJI | 4944 | 10 | 10 (yr19) | 99.8% | 9918 | 9 | -0.0244 | -0.0088 | +0.0898 |
| Newport | 5000 | 7 | 7 (yr19) | 99.9% | 9116 | 6 | -0.0253 | +0.0209 | +0.0696 |
| Monterey | 4999 | 586 | 33 (yr9) | 99.3% | 10003 | 22 | +0.0186 | +0.0609 | +0.0864 |

## Comparison to Single-Trait Baseline (Phase 14)

Previous run: 51 resistance loci, same 5 nodes, K=5000, seed=42, 20yr.

| Node | Crash (old) | Crash (new) | Δ | Δr (old) | Δr (new) | Δt (new) | Δc (new) |
|------|-------------|-------------|---|----------|----------|----------|----------|
| Sitka | 95.4% | 98.8% | +3.4pp | +0.031 | -0.0637 | -0.0046 | +0.0458 |
| Howe Sound | 93.2% | 99.4% | +6.2pp | +0.044 | +0.0583 | -0.0236 | +0.0202 |
| SJI | 99.3% | 99.8% | +0.5pp | +0.061 | -0.0244 | -0.0088 | +0.0898 |
| Newport | 99.1% | 99.9% | +0.8pp | +0.055 | -0.0253 | +0.0209 | +0.0696 |
| Monterey | 98.8% | 99.3% | +0.5pp | +0.051 | +0.0186 | +0.0609 | +0.0864 |

### Key Differences

With three traits (17 loci each) vs one trait (51 loci):
- Each trait has fewer loci → less genetic variance per trait
- Resistance (r_i) alone is weaker → tolerance and recovery provide alternative survival pathways
- Selection acts on all three traits simultaneously

## Emergent Dynamics

### Trait Evolution Rates

- Mean Δr across nodes: -0.0073
- Mean Δt across nodes: +0.0089
- Mean Δc across nodes: +0.0624

**Fastest evolving trait:** recovery (c_i)

### Recovery Signal

Nodes with recovered (R-state) individuals at end: ['Sitka', 'Howe Sound', 'SJI', 'Newport', 'Monterey']

### Trait Correlations

See console output for r(r,t), r(r,c), r(t,c) per node.
Independent loci → expect low correlations unless shared selection pressure creates linkage.

## Figures

1. `population_trajectories.png` — Population over time per node
2. `trait_evolution.png` — Three-panel trait mean trajectories
3. `trait_comparison.png` — Δr vs Δt vs Δc bar chart per node

## Summary

Total initial population: 25000
Total final population: 766
Overall crash: 96.9%
