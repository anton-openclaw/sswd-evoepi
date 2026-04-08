# Report Brief — Reintroduction Experiment 2050

## Title
Captive-Bred Sunflower Star Reintroduction: A Full-Factorial Analysis of Resistance Genetics, Release Density, and Release Location

## Audience
Willem Weertman (marine biology PhC, UW) — direct tone, use jargon freely. This is internal research reporting, not a public paper.

## Data Sources
- **40 reintroduction scenarios**: `/home/starbot/.openclaw/workspace/sswd-evoepi/results/reintro_2050/REINTRO_L{1,2}_R{1-5}_D{1-4}/`
  - Each contains: `result_seed42.json`, `combined_results.json`, `monthly_seed42.npz`, `checkpoint_seed42.json`
- **Baseline (no intervention)**: `/home/starbot/.openclaw/workspace/sswd-evoepi/results/forecast_2050/F01_baseline/`
- **Experiment design**: `/home/starbot/.openclaw/workspace/sswd-evoepi/experiments/reintroduction/README.md`

## Experimental Design
Full factorial: 5 resistance genetics × 4 release densities × 2 release locations = 40 scenarios
- **Location**: L1=Monterey CA, L2=Friday Harbor WA
- **Resistance**: R1=CA-2025 survivors (r̄~0.31), R2=Pre-epidemic naive (r̄~0.15), R3=Mid-epidemic 2019 (r̄~0.20), R4=Artificially selected (r̄~0.50), R5=Immune (r̄=1.0)
- **Density**: D1=50/site, D2=200/site, D3=500/site, D4=1000/site (=K)
- All branch from W330 baseline at year 13 (2025), run to year 38 (2050)
- Release into all CA sites (238 nodes) for L1, or all JDF sites for L2
- K=1000, K_ref=5000, seed 42

## Baseline Reference (F01, no reintroduction)
- Final pop: 8.3% of pre-epidemic
- CA-S: 0%, CA-C: 0%, CA-N: 0.03% (effectively extinct in California)
- JDF: 7.2%, AK-PWS: 27.0%, BC-N: 6.3%

## Key Findings to Highlight

### 1. COMPETITIVE DILUTION (⚠️ Major finding)
Releasing low-resistance animals (R1-R3) at HIGH density HARMS the population:
- R1/D1 → CA-C 2.9%; R1/D4 → CA-C 0.2% (WORSE with more animals!)
- R2/D4 → CA-C 0.0% (total extirpation despite releasing 1000/site)
- Mechanism: naive animals compete for resources but die to disease, diluting resistant survivor genetics

### 2. RESISTANCE IS EVERYTHING
- Only R5 (immune, r̄=1.0) achieves substantial CA recovery: CA-N 21.3%, CA-C 14.6%, CA-S 8.4%
- R4 (artificially selected, r̄=0.50) moderate: CA-N 7.6%, CA-C 5.3%, CA-S 2.8%
- R1-R3 negligible to harmful at scale
- Answers the key question: "breed for resistance, not numbers"

### 3. LOCATION × SPILLOVER
- L1 (Monterey) → boosts CA regions, no JDF effect
- L2 (Friday Harbor) → boosts JDF to 27.7% (R5/D4), ZERO CA spillover
- No larval transport from WA to CA (southward transport absent)

### 4. ALASKA IS SELF-SUSTAINING
- AK-PWS ~26-27% across ALL 40 scenarios — completely insensitive to reintroduction
- Northern populations don't need or benefit from intervention

### 5. BC-N REMAINS A PUZZLE
- ~6% across all scenarios — neither CA nor WA releases help

### 6. DENSITY CEILING
- Diminishing returns above D2 (200/site) for R4-R5
- For R1-R3, increasing density is actively counterproductive

## Report Structure
1. **Executive Summary** (1 page) — Key findings, policy recommendations
2. **Experimental Design** (1-2 pages) — Factorial design, baseline, model description
3. **Results**
   3.1 Overview: Full heatmap of all 40 scenarios
   3.2 Location Effect: L1 vs L2 comparison
   3.3 Resistance Genetics: The central role of r̄
   3.4 Density Response: When more is less
   3.5 Competitive Dilution: ⚠️ The counterintuitive finding
   3.6 Temporal Dynamics: Recovery trajectories post-release
   3.7 Spillover Assessment: Regional propagation
   3.8 Alaska Stability: Self-sustaining northern populations
4. **Discussion** — Conservation implications, breeding program recommendations
5. **Conclusions** — 3 bullet actionable recommendations
6. **Appendix** — Full results table, config details

## Figures Needed (14-16)
1. **Full factorial heatmap** (2 panels: L1, L2) — genetics × density → final CA/JDF recovery %
2. **Baseline comparison bar chart** — key regions: baseline vs best L1 vs best L2
3. **CA time series** — CA-C yearly trajectory for R1/R4/R5 at D2, plus baseline
4. **JDF time series** — JDF yearly trajectory for L2 R1/R4/R5 at D2/D4, plus baseline
5. **Competitive dilution panel** — CA-C recovery vs density for each R level (L1)
6. **Resistance threshold curve** — CA-C final recovery vs r̄ at fixed D2
7. **Dose-response curves** — final recovery vs density for each R level
8. **Spillover heatmap** — region × scenario showing propagation distance from release site
9. **AK-PWS stability panel** — overlay of all 40 AK-PWS trajectories showing robustness
10. **Regional recovery map** — geographic scatter/map for best scenario vs baseline
11. **Genetics × density interaction plot** — 2D contour or interaction plot
12. **OR gradient** — Oregon as intermediate zone between CA release and WA populations
13. **Population composition** — released vs wild-born at CA-C over time (if data available in NPZ)
14. **Full results summary table** — all 40 scenarios + baseline

## Min/Max Pages
20-30 pages

## Delivery
Email to wlweert@gmail.com with subject "Reintroduction Experiment Report — 40 Scenarios (2025→2050)"
