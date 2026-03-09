# SSWD-EvoEpi Calibration Report: W125-W134 Per-Site Flushing Analysis

## Overview

This report analyzes the first implementation of per-site flushing rates derived from enclosedness metrics in the SSWD-EvoEpi model. The analysis covers calibration runs W125-W132 (W133-W134 still running).

## Key Findings

- **Performance Degradation**: Per-site flushing implementation shows worse performance compared to uniform connectivity baselines
- **Best Run**: W129 (α_env=0.25) achieved RMSE = 0.691
- **Baseline Comparison**: 
  - +14.6% worse than W117 (RMSE = 0.603)
  - +16.9% worse than W91 (RMSE = 0.591)
- **Systematic Underperformance**: All regions show reduced recovery compared to targets

## Files Generated

### Main Report
- `main.pdf` - Complete 10-page report with analysis and discussion
- `main.tex` - LaTeX source file
- `generate_report.py` - Python script used to generate all figures

### Figures (PNG format, 300 DPI)
1. `rmse_comparison.png` - RMSE comparison across all runs vs baselines
2. `recovery_trajectories.png` - 13-year recovery trajectories for best run (W129)
3. `target_vs_actual.png` - Target vs actual recovery scatter plot (log-log)
4. `alpha_sensitivity.png` - Environmental transmission parameter sensitivity
5. `recovery_heatmap.png` - Regional recovery percentages across all runs
6. `flushing_impact.png` - Direct comparison W129 vs W117

### Data
- `summary_data.json` - Key statistics and parameters used in analysis

## Key Question Addressed

**Did per-site flushing help or hurt?**

The analysis clearly shows that per-site flushing **hurt** performance compared to uniform connectivity approaches. Possible explanations:

1. **Conservative flushing rates**: Enclosedness-derived rates may be too low
2. **Spatial bottlenecks**: Per-site approach may impede recolonization
3. **Scale mismatch**: Site-level enclosedness may not capture larval dispersal scale
4. **Calibration bias**: Previous runs may have compensated via high uniform connectivity

## Recommendations

1. Validate enclosedness calculation algorithms
2. Assess spatial scale appropriateness 
3. Consider hybrid uniform/per-site approach
4. Adjust other parameters to compensate for realistic connectivity

## Technical Details

- **Generated**: March 8, 2026
- **Data source**: `/home/starbot/.openclaw/workspace/sswd-evoepi/results/calibration/W125-W132/`
- **Tools**: Python 3, matplotlib, numpy, pdflatex
- **Model version**: SSWD-EvoEpi with pathogen evolution enabled