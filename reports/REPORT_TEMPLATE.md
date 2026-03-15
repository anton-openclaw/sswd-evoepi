# SSWD-EvoEpi Calibration Report Template

**Version:** 1.0 (2026-03-15)
**Minimum pages:** 35-45 (30+ figures, ~5 pages text)
**Text-to-figure ratio:** Each figure gets 2-3 sentence caption. No more than 1 paragraph of discussion per subsection.

## Design Principles

1. **Standardized** — Same structure every sweep. Compare across runs by section.
2. **Figure-heavy** — Figures tell the story. Text is connective tissue.
3. **Comprehensive** — Calibration metrics + temporal dynamics + spatial patterns + evolution.
4. **Inspectable** — Every figure answers a specific question a reviewer would ask.

---

## Section Structure

### Section 1: Executive Summary (1 page)
- Sweep name, date, N configs, seed(s), K, years
- Best config: name, RMSE, key recovery values
- 3-5 bullet key findings
- Comparison to previous best

### Section 2: Calibration Metrics (8-10 figures, ~12 pages)

| Fig | ID | Description | Viz Function | Data | Priority |
|-----|-----|-------------|-------------|------|----------|
| 2.1 | `fig_rmse_ranking` | RMSE ranking bar chart, all configs, historical baselines | `calibration.plot_rmse_ranking` | CalibrationData | P0 |
| 2.2 | `fig_region_heatmap` | Configs × 8 scored regions recovery heatmap | `calibration.plot_region_heatmap` | CalibrationData | P0 |
| 2.3 | `fig_best_vs_targets` | Best config: actual vs target bars, all 8 regions | `calibration.plot_config_vs_targets` | CalibrationData | P0 |
| 2.4 | `fig_parameter_effects` | Multi-panel parameter sensitivity curves | `calibration.plot_parameter_effects` | CalibrationData | P0 |
| 2.5 | `fig_tradeoff` | AK-PWS vs CA-N latitude tradeoff scatter | `calibration.plot_tradeoff` | CalibrationData | P0 |
| 2.6 | `fig_arrival_timing` | Wavefront arrival: model vs observed (grouped bars) | **NEW: `plot_arrival_timing_bars`** | JSON arrival_timing + ARRIVAL_TARGETS | P0 |
| 2.7 | `fig_grading_table` | Colored table: within 2×/5× per region per config | **NEW: `plot_grading_table`** | metrics.score_regions() | P1 |
| 2.8 | `fig_seed_variability` | RMSE spread across seeds (box plot) | `calibration.plot_seed_variability` | CalibrationData | P1 |

### Section 3: Best Config — Temporal Dynamics (12-15 figures, ~18 pages)

*All figures in this section show the BEST configuration (lowest RMSE). For multi-seed runs, show mean ± range.*

#### 3A: Population Dynamics

| Fig | ID | Description | Viz Function | Data | Priority |
|-----|-----|-------------|-------------|------|----------|
| 3.1 | `fig_total_population` | Total population trajectory, all 896 sites summed | adapt `calibration.plot_trajectories` | JSON yearly_totals | P0 |
| 3.2 | `fig_regional_pop_grid` | 18-panel grid: population/K per region, N→S | **NEW: `plot_regional_trajectory_grid`** | JSON yearly_totals | P0 |
| 3.3 | `fig_pop_heatmap_monthly` | Heatmap: 159 months × 18 regions, pop as fraction of K | adapt `population.plot_population_heatmap` | NPZ populations | P0 |
| 3.4 | `fig_crash_severity` | Bar chart: crash % per region, N→S | **NEW: `plot_crash_severity_bars`** | JSON crash_pct | P0 |
| 3.5 | `fig_recruitment_vs_deaths` | Annual recruits vs disease deaths, aggregated | adapt `population.plot_recruitment_timeseries` | JSON yearly_recruits, yearly_disease_deaths | P0 |

#### 3B: Disease Epidemiology

| Fig | ID | Description | Viz Function | Data | Priority |
|-----|-----|-------------|-------------|------|----------|
| 3.6 | `fig_infection_prevalence` | Monthly infection prevalence by latitude band | **NEW: `plot_monthly_infection_prevalence`** | NPZ infected/populations | P0 |
| 3.7 | `fig_infection_heatmap` | Heatmap: 159 months × 18 regions, infection prevalence | adapt `disease.plot_disease_state_heatmap` | NPZ infected/populations | P0 |
| 3.8 | `fig_disease_deaths_stacked` | Stacked area: annual disease deaths by region | **NEW: `plot_disease_deaths_stacked`** | JSON yearly_disease_deaths | P0 |
| 3.9 | `fig_epidemic_timeline` | Gantt-style: disease phase per region (free/active/recovery) | adapt `spatial.plot_spatial_epidemic_timeline` | NPZ infected | P1 |

#### 3C: Evolutionary Genetics

| Fig | ID | Description | Viz Function | Data | Priority |
|-----|-----|-------------|-------------|------|----------|
| 3.10 | `fig_resistance_trajectories` | Mean resistance over 13yr, all 18 regions overlaid | adapt `genetics.plot_resistance_trajectory` | JSON yearly_mean_resistance | P0 |
| 3.11 | `fig_tolerance_recovery_trajectories` | Mean tolerance + recovery traits, 2-panel | **NEW: `plot_trait_evolution_panel`** | JSON yearly_mean_tolerance/recovery | P0 |
| 3.12 | `fig_va_resistance` | Additive variance for resistance, all regions | adapt `genetics.plot_additive_variance_over_time` | JSON yearly_va_resistance | P0 |
| 3.13 | `fig_selection_differential` | Year-over-year Δr̄ (selection intensity) | `genetics.plot_selection_differential` | JSON yearly_mean_resistance | P1 |

#### 3D: Pathogen Co-Evolution

| Fig | ID | Description | Viz Function | Data | Priority |
|-----|-----|-------------|-------------|------|----------|
| 3.14 | `fig_pathogen_tvbnc` | Final T_vbnc by region (pathogen thermal adaptation) | **NEW: `plot_pathogen_adaptation_bars`** | JSON final_mean_T_vbnc | P0 |
| 3.15 | `fig_pathogen_virulence` | Final virulence by region | **NEW: `plot_virulence_adaptation_bars`** | JSON final_mean_v_local | P0 |

### Section 4: Best Config — Spatial Dynamics (10-12 figures, ~15 pages)

#### 4A: Geographic Context

| Fig | ID | Description | Viz Function | Data | Priority |
|-----|-----|-------------|-------------|------|----------|
| 4.1 | `fig_site_map` | 896-node network map, colored by region | adapt `spatial.plot_network_map` | NPZ lat/lon/names | P0 |
| 4.2 | `fig_fjord_depth_by_region` | Box plot of fjord depth norm per region | `salinity.plot_fjord_depth_by_region` | NodeDefinition list | P0 |

#### 4B: Disease Wavefront

| Fig | ID | Description | Viz Function | Data | Priority |
|-----|-----|-------------|-------------|------|----------|
| 4.3 | `fig_spatial_snapshots` | 2×4 map grid: population + infection at 4 timepoints | **NEW: `plot_spatial_snapshots_grid`** | NPZ populations, infected, lat/lon | P0 |
| 4.4 | `fig_wavefront_scatter` | Latitude vs first infection month, all 896 sites | **NEW: `plot_wavefront_latitude_scatter`** | NPZ infected, site_lats | P0 |

#### 4C: North-South Gradient

| Fig | ID | Description | Viz Function | Data | Priority |
|-----|-----|-------------|-------------|------|----------|
| 4.5 | `fig_recovery_vs_latitude` | Recovery fraction vs latitude, 18 regions + targets | adapt `spatial.plot_north_south_gradient` | JSON region_recovery | P0 |
| 4.6 | `fig_all_region_trajectories` | 18-panel small multiples: monthly pop/K per region | **NEW: `plot_all_region_trajectories`** | NPZ populations by region | P0 |

#### 4D: Fjord Refuge & Salinity

| Fig | ID | Description | Viz Function | Data | Priority |
|-----|-----|-------------|-------------|------|----------|
| 4.7 | `fig_fjord_vs_open` | Population trajectories: high-fjord vs open-coast sites | adapt `spatial.plot_fjord_vs_open` | NPZ + NodeDefinition | P0 |
| 4.8 | `fig_suppression_map` | Geographic scatter: June disease suppression (sal_mod) | `salinity.plot_suppression_map` | NodeDefinition + config | P0 |
| 4.9 | `fig_suppression_vs_recovery` | Scatter: regional suppression vs recovery fraction | `salinity.plot_suppression_vs_recovery` | NodeDefinition + JSON | P0 |
| 4.10 | `fig_latitude_asymmetry` | Latitude vs suppression, colored by fjord depth | `salinity.plot_latitude_asymmetry` | NodeDefinition + config | P1 |

### Section 5: Discussion & Next Steps (1 page)
- Key achievements vs previous sweep
- Remaining gaps (AK recovery, latitude gradient)
- Parameter space narrowing
- Recommended next sweep design (2-3 sentences)

---

## Figure Count Summary

| Section | P0 | P1 | Total |
|---------|----|----|-------|
| 2. Calibration | 6 | 2 | 8 |
| 3. Temporal | 12 | 2 | 14 |
| 4. Spatial | 9 | 1 | 10 |
| **Total** | **27** | **5** | **32** |

**P0 figures alone = 27 figures = ~32 pages.** With P1 = 32 figures = ~38 pages.

---

## New Viz Functions Required (14 total)

### Must implement before first use:

1. **`plot_regional_trajectory_grid(result_json, regions, targets)`** — 18-panel grid from yearly_totals
2. **`plot_crash_severity_bars(result_json, regions)`** — Horizontal bars of crash_pct N→S
3. **`plot_monthly_infection_prevalence(npz_path, regions)`** — Monthly prevalence by latitude band
4. **`plot_disease_deaths_stacked(result_json, regions)`** — Stacked area of annual deaths
5. **`plot_trait_evolution_panel(result_json, regions, traits)`** — Multi-trait trajectories
6. **`plot_pathogen_adaptation_bars(result_json, regions)`** — T_vbnc bars by region
7. **`plot_virulence_adaptation_bars(result_json, regions)`** — v_local bars by region
8. **`plot_spatial_snapshots_grid(npz_path, timepoints)`** — 2×N map panels at key months
9. **`plot_wavefront_latitude_scatter(npz_path)`** — First infection month vs latitude
10. **`plot_all_region_trajectories(npz_path)`** — 18-panel monthly from NPZ
11. **`plot_arrival_timing_bars(result_json, targets)`** — Model vs observed arrival
12. **`plot_grading_table(results, targets)`** — Colored 2×/5× table

### Adapt existing (minor changes):

13. `population.plot_population_heatmap` → accept NPZ regional aggregation
14. `disease.plot_disease_state_heatmap` → accept NPZ regional aggregation
15. `genetics.plot_resistance_trajectory` → accept multi-region JSON data
16. `genetics.plot_additive_variance_over_time` → accept multi-region JSON data
17. `spatial.plot_north_south_gradient` → accept JSON region_recovery
18. `spatial.plot_fjord_vs_open` → accept NPZ + fjord_depth_norm

---

## Data Pipeline

### Computed Once (data_summary.json):
- All config RMSE, recovery per region, per-seed values
- Historical baseline references
- Parameter values per config
- Rankings

### Computed Per Figure:
- NPZ monthly aggregation to regions (cached after first computation)
- Trait trajectories from JSON
- Salinity suppression from config + NodeDefinition

### NPZ Region Aggregation Helper:
```python
def aggregate_npz_to_regions(npz_path):
    """Load NPZ, aggregate 896 sites to 18 regions.
    Returns dict: region -> {pop: array(159,), infected: array(159,), months: array(159,)}
    """
```
This function should be in `sswd_evoepi/results.py` and called once per run.

---

## LaTeX Template Skeleton

```latex
\documentclass[11pt,a4paper]{article}
\usepackage[margin=1in]{geometry}
\usepackage{graphicx, booktabs, xcolor, hyperref, float}
\usepackage[font=small,labelfont=bf]{caption}

\definecolor{bestgreen}{HTML}{4CAF50}

\title{SSWD-EvoEpi Calibration Report: [SWEEP_NAME]}
\author{Anton Star \& Willem Weertman}
\date{\today}

\begin{document}
\maketitle

\section{Executive Summary}
% 1 page, standalone, key findings

\section{Calibration Metrics}
% Figs 2.1-2.8

\section{Best Configuration: Temporal Dynamics}
\subsection{Population Dynamics}  % 3.1-3.5
\subsection{Disease Epidemiology}  % 3.6-3.9
\subsection{Evolutionary Genetics}  % 3.10-3.13
\subsection{Pathogen Co-Evolution}  % 3.14-3.15

\section{Best Configuration: Spatial Dynamics}
\subsection{Geographic Context}  % 4.1-4.2
\subsection{Disease Wavefront}  % 4.3-4.4
\subsection{North-South Recovery Gradient}  % 4.5-4.6
\subsection{Fjord Refuge \& Salinity}  % 4.7-4.10

\section{Discussion \& Next Steps}
% 1 page

\end{document}
```

---

## Quality Checklist

- [ ] All P0 figures present and render (27 minimum)
- [ ] All cross-references resolve
- [ ] Every number in text traces to data_summary.json
- [ ] Executive summary standalone (no forward references)
- [ ] Figures inline with relevant text
- [ ] Captions descriptive (figure understandable without text)
- [ ] 35+ pages total
- [ ] Consistent color scheme across all figures
- [ ] Region ordering always N→S
- [ ] Historical baselines included in calibration section
