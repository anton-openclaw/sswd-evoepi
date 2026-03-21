# Wavefront Calibration W17-W28 Results

**Date:** 2026-03-01
**Shared parameters:** k_vbnc=2.0, activation_threshold=500, 5 Channel Islands origins, D_P=50km, D_P_max_range=175km

## Key Findings

### 1. Wavefront Stalls at BC-N — Never Reaches Alaska
- In ALL 12 runs, disease never reaches AK-FS, AK-FN, AK-PWS, AK-AL, AK-EG, AK-WG, or AK-OC
- Wave reaches BC-N at ~23-34 months but cannot bridge the gap to Alaska
- Alaska populations remain at 89-96% (untouched by disease)
- **Root cause:** activation_threshold=500 bact/mL too high for the diluted vibrio concentrations reaching BC→AK distances

### 2. Strong Southern Gradient Achieved
- CA-N and OR go completely extinct (0%) — matches targets (0.1% and 0.25%)
- CA-S and CA-C destroyed immediately upon disease arrival
- SS-S ranges 0.1-6.7% across runs (target 5%) — some runs nail this
- JDF ranges 0.1-7.4% (target 2%)

### 3. BC-N Overshoot Problem
- BC-N ranges 0-10.2% (target 20%)
- Disease is too severe once it arrives — wipes out populations before recovery can establish
- Higher K_half (200K) helps: W21=10.2%, W22=7.1%

### 4. Wave Timing
- CA-C arrival: ~12 months (target 6) — about 2× too slow
- OR arrival: ~16-20 months (target 15) — close
- BC-N arrival: ~23-34 months (target ~26) — in range
- Alaska: NEVER (target ~42 months)

### 5. Parameter Effects
- **K_half:** Higher K_half → slightly more recovery in mid-coast (BC-N, SS-S, JDF). No effect on wavefront propagation.
- **s0 (settler survival):** Higher s0 → more recovery (more recruits). W22 (s0=0.003, K_half=200K) has best mid-coast numbers.
- **P_env_max:** 1500 vs 2000 vs 2500 — minimal effect in these runs (wave dynamics dominate over endemic pressure).

## What Needs to Change for Next Round

### Priority 1: Get wavefront to Alaska
- **Lower activation_threshold** from 500 to 50-100 bact/mL
- **Increase D_P_max_range** from 175km to 500-1000km (longer-range waterborne dispersal)
- Consider both simultaneously

### Priority 2: Protect BC-N recovery
- Once wavefront reaches AK, need to ensure BC-N still recovers to ~20%
- K_half=200K + s0≥0.002 seems necessary for any mid-coast recovery

### Priority 3: Speed up southern wave
- CA-C at 12mo vs 6mo target — wave propagation too slow near origin
- May need higher local D_P or lower threshold

## Summary Table

| Round | K_half | s0 | P_env | AK-PWS | AK-FN | AK-FS | BC-N | SS-S | JDF | OR | CA-N | AK Arrival | BC-N Arrival |
|-------|--------|------|-------|--------|-------|-------|------|------|-----|-----|------|------------|--------------|
| W17 | 50K | 0.001 | 2000 | 88.9 | 89.4 | 89.2 | 0.0 | 0.1 | 0.1 | 0.0 | 0.0 | NEVER | 24mo |
| W18 | 100K | 0.001 | 2000 | 89.4 | 89.3 | 89.2 | 0.3 | 0.6 | 0.6 | 0.0 | 0.0 | NEVER | 31mo |
| W19 | 50K | 0.002 | 2000 | 92.8 | 93.0 | 92.6 | 0.0 | 0.1 | 0.1 | 0.0 | 0.0 | NEVER | 24mo |
| W20 | 100K | 0.002 | 2000 | 92.6 | 93.3 | 92.4 | 0.4 | 0.6 | 0.9 | 0.0 | 0.0 | NEVER | 31mo |
| W21 | 200K | 0.002 | 2000 | 92.3 | 92.9 | 92.7 | 10.2 | 5.7 | 6.2 | 0.4 | 0.0 | NEVER | 35mo |
| W22 | 200K | 0.003 | 2000 | 95.5 | 96.0 | 95.7 | 7.1 | 6.7 | 7.4 | 0.7 | 0.1 | NEVER | 34mo |
| W23 | 50K | 0.001 | 1500 | 89.1 | 89.4 | 89.4 | 0.0 | 0.1 | 0.1 | 0.0 | 0.0 | NEVER | 29mo |
| W24 | 100K | 0.001 | 1500 | 89.3 | 89.6 | 89.2 | 4.4 | 0.7 | 0.8 | 0.0 | 0.0 | NEVER | 32mo |
| W25 | 50K | 0.002 | 2500 | 92.5 | 93.1 | 92.7 | 0.0 | 0.1 | 0.1 | 0.0 | 0.0 | NEVER | 23mo |
| W26 | 100K | 0.002 | 2500 | 92.9 | 93.3 | 92.7 | 0.2 | 0.5 | 0.6 | 0.0 | 0.0 | NEVER | 30mo |
| W27 | 50K | 0.003 | 2000 | 95.5 | 96.1 | 95.9 | 0.0 | 0.2 | 0.1 | 0.0 | 0.0 | NEVER | 24mo |
| W28 | 100K | 0.003 | 2000 | 95.6 | 96.5 | 96.1 | 0.6 | 1.0 | 1.4 | 0.0 | 0.0 | NEVER | 31mo |

**Targets:** AK-PWS=50%, AK-FN=50%, AK-FS=20%, BC-N=20%, SS-S=5%, JDF=2%, OR=0.25%, CA-N=0.1%

## Figures
- `fig1_recovery_trajectories_all.png` — Recovery trajectories for all 12 rounds (18 regions each)
- `fig2_best_rounds_detail.png` — Detailed trajectories for selected rounds
- `fig3_wavefront_timing.png` — Wavefront arrival timing vs targets
- `fig4_parameter_sensitivity.png` — Parameter sensitivity & diagnostics
- `fig5_summary_heatmap.png` — Summary heatmap of final recovery vs targets