# Publishable Report: SSWD-EvoEpi Forecast Analysis
## Target: Publication-quality, first-principles framing

### Audience
- Primary: Jason Hodin (larval biology, echinoderms)
- Secondary: Conservation biologists, wildlife managers
- Tone: Accessible science, not impenetrable modeling jargon

### Structure

## 1. Introduction (~2 pages)
- SSWD devastation: 90.6% decline, Pycnopodia critically endangered
- V. pectenicida identified (Prentice 2025 Koch's postulates)
- The conservation challenge: captive breeding underway, but where to release?
- Why forecasting matters: need to predict 2025-2050 trajectories
- What this report does: 4 forecast scenarios testing climate × evolution

## 2. Key Findings (~3 pages) 
**Put the interesting results UPFRONT**
- Universal crash: ALL scenarios → 94-95% population loss by 2050
- Regional winners: Oregon emerges as #1 natural recovery region
- Climate paradox: warming HELPS Alaska but HURTS southern recovery
- Evolutionary arms race: host evolution is real but pathogen adapts faster
- Virulence evolution barely matters (surprising null result)
- Implications for reintroduction: southern sites need help, northern sites might self-recover

## 3. Scenario Design (~1 page)
- F01: SSP2-4.5 baseline (moderate warming, all evolution)
- F02: SSP2-4.5 + NO pathogen evolution (thermal or virulence)
- F03: SSP2-4.5 + thermal adaptation ON, virulence evolution OFF
- F04: SSP5-8.5 (high warming, all evolution)
- Common setup: W154 production model, 2012-2050, 3 seeds each
- SST projection method: delta method from CMIP6

## 4. Results in Detail (~6-8 pages)
### 4.1 Population Trajectories
- Coast-wide aggregate: crash timing, recovery attempts, final state
- Regional trajectories: North vs South divergence
- The recruitment pulse phenomenon (seasonal ratchet)
  
### 4.2 Regional Recovery Patterns
- Table: All 18 regions × 4 scenarios (recovery fraction, final pop)
- Heat map: Latitude vs time showing the "recovery gradient"
- Winners and losers by scenario

### 4.3 Climate Effects (F01 vs F04)
- Northern shift under warming: AK and BC benefit
- Southern penalty: CA-C and OR decline
- Net effect near-zero (a wash overall)
- SST projections driving the differences

### 4.4 Evolutionary Dynamics
- Host resistance trajectories to 2050 (extension of calibration period)
- Genetic diversity collapse in southern populations
- Pathogen thermal adaptation (reaches floor in north)
- Virulence gradient emergence and persistence

### 4.5 The Pathogen Evolution Question (F01 vs F02 vs F03)
- Full evo OFF (F02): modest benefit to AK (+2-4%), minimal elsewhere
- Virulence evo OFF only (F03): negligible effect (< 1% difference)
- Interpretation: thermal adaptation matters more than virulence evolution
- Why: pathogen community selection operates on different timescale

## 5. Model Description (~4-5 pages, for those who want mechanics)
- Mine from progress report
- 5.1 Agent-based framework (896 sites, 5000 individuals, daily timesteps)
- 5.2 Disease mechanics (SEIR with VBNC, spatial transmission, wavefront)
- 5.3 Host genetics (51 loci, 3 traits, natural selection)
- 5.4 Pathogen evolution (community selection, virulence-transmission trade-off)
- 5.5 Reproduction (broadcast spawning, SRS, Allee effects, larval dispersal)
- 5.6 Environmental forcing (real SST → CMIP6 projections)
- 5.7 Calibration summary (W154, RMSE=0.504, 8 regional targets)

## 6. Discussion (~2 pages)
- What the model says about recovery prospects without intervention
- The Oregon opportunity
- Alaska's conditional optimism
- Southern California as a conservation priority
- Limitations and uncertainties
- Next steps: reintroduction scenarios (F05-F08), sensitivity analysis

## 7. Methods Appendix
- Parameter table (from progress report)
- SST projection methodology
- Calibration targets

## Figures List (Target: 16-20 figures)
1. Network map (896 sites, 18 regions) — from progress report
2. Population trajectory overview (all regions, F01)
3. Regional recovery comparison (bar chart, 4 scenarios)
4. Coast-wide population over time (aggregate, all 4 scenarios)
5. Regional heatmap (latitude × time × population fraction)
6. Climate comparison detail (F01 vs F04, selected regions)
7. SST projections (selected sites, SSP2-4.5 vs SSP5-8.5)
8. Host resistance evolution to 2050
9. Genetic diversity (additive variance) trajectories
10. Pathogen thermal adaptation trajectories
11. Virulence gradient map (final state, F01)
12. Evolutionary toggle comparison (F01 vs F02 vs F03)
13. Within-site disease dynamics (from progress report)
14. Disease wavefront snapshots (4-panel, from progress report)
15. Recovery vs latitude scatter (all 4 scenarios)
16. SST-recovery correlation
