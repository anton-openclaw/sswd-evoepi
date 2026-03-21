# Report Brief: Parametric vs WOA23 Salinity Comparison

## Title
Parametric vs Observational Salinity: Validation of the SSWD-EvoEpi Seasonal Salinity Model Against WOA23

## Audience
Willem Weertman — PhD candidate, marine biology/ecology. Understands oceanography, salinity dynamics, the SSWD model. Write for a scientist, not a layperson.

## Key Questions
1. How well does the parametric baseline (`S_ocean = 31.32 + 0.054 × (lat - 50)`) match WOA23 surface salinity?
2. Where does the parametric model fail most? (coastal vs open ocean, by region, by latitude)
3. How different are the seasonal cycles? The parametric model has NO seasonal variation in its baseline — all seasonality comes from the melt pulse. WOA23 captures real seasonal cycles.
4. When does minimum salinity actually occur? (The parametric model assumes June 15 everywhere via cosine pulse)
5. What fraction of sites needed nearest-neighbor fill at 0.25° resolution? How far are those NN distances?
6. Should the model use WOA23 monthly climatologies instead of the parametric baseline?
7. How does this affect the disease suppression predictions — especially the AK vs CA asymmetry?

## Data Sources
- **WOA23 site data**: `/home/starbot/.openclaw/workspace/salinity_validation/woa23_site_monthly.npz`
  - `names`: (896,) site names
  - `lats`, `lons`: (896,) coordinates
  - `regions`: (896,) subregion codes
  - `fjord_depth_norm`: (896,) normalized fjord depth [0,1]
  - `woa_monthly`: (896, 12) WOA23 monthly salinity (NN-filled)
  - `woa_direct`: (896, 12) WOA23 direct extraction (NaN for land cells)
  - `param_monthly_0`: (896, 12) parametric baseline (fw_strength=0)
  - `param_monthly_15`: (896, 12) parametric with fw_strength=15
  - `nn_distance`: (896,) distance to nearest WOA23 ocean cell in degrees
- **WOA23 site CSV**: `/home/starbot/.openclaw/workspace/salinity_validation/woa23_site_monthly.csv`
- **DFO lighthouse data**: `/home/starbot/.openclaw/workspace/salinity_validation/dfo_monthly_climatology.csv`
- **Site enclosedness**: `/home/starbot/.openclaw/workspace/sswd-evoepi/data/nodes/site_enclosedness.csv`

## Key Summary Stats (from extraction)
- 896 sites, 235 with direct WOA23 cells, 661 needed NN fill (max 0.177°)
- Parametric baseline vs WOA23: Mean bias = -0.130 psu, RMSE = 2.421 psu, Max error = 9.184 psu
- WOA23 seasonal range: mean 1.676 psu. Parametric baseline seasonal range: 0.000 psu (constant!)
- Month of minimum salinity varies widely: May (171 sites), Jun (152), Feb (116), Aug (102), etc.

## Figures Needed (generate with matplotlib, save as PDF)
1. **Scatter: Parametric vs WOA23 annual mean** — each site a point, colored by region. 1:1 line. Show R², RMSE.
2. **Bias map** — geographic scatter of all 896 sites, colored by (parametric - WOA23) annual mean bias.
3. **Seasonal timing map** — geographic scatter colored by WOA23 month of minimum salinity. Shows spatial pattern of seasonal timing.
4. **Regional seasonal profiles** — 4-panel (AK, BC, Salish, OR/CA), each showing: WOA23 mean±SD (shaded), parametric baseline (dashed), parametric+fw_strength=15 (solid). Shows how the melt pulse adds seasonality the baseline lacks.
5. **Seasonal range comparison** — scatter of WOA23 seasonal range vs parametric range, by region. Shows parametric model has zero baseline seasonality.
6. **NN fill quality** — map showing which sites needed NN fill and at what distance.
7. **DFO validation** — overlay DFO lighthouse monthly climatologies on WOA23 values at same locations. Independent check on WOA23 coastal accuracy.
8. **Suppression comparison** — June disease suppression under parametric-only vs WOA23-informed model. Side-by-side maps showing how using real data changes the spatial pattern.
9. **Latitude gradient** — scatter of latitude vs (a) WOA23 June salinity, (b) parametric June salinity, (c) difference. Key plot showing how the AK-CA asymmetry differs.

## Structure
1. Executive Summary (1 page)
2. Data & Methods — WOA23 description, extraction method, parametric model equations, NN fill approach
3. Baseline Comparison — annual mean, bias, RMSE by region (Figs 1-2)
4. Seasonal Dynamics — timing of minima, seasonal amplitude (Figs 3-5)
5. Regional Profiles — detailed AK vs BC vs SS vs CA comparison (Fig 4)
6. Data Quality — NN fill distances, DFO validation (Figs 6-7)
7. Implications for Disease Model — suppression maps, latitude gradient (Figs 8-9)
8. Recommendation — should we switch to WOA23-informed baseline?
9. References

## Compilation
```
export PATH=$HOME/.TinyTeX/bin/x86_64-linux:$PATH
cd /home/starbot/.openclaw/workspace/sswd-evoepi/reports/salinity-comparison
pdflatex -interaction=nonstopmode report.tex
```

## Delivery
Email via: `python3 /home/starbot/.openclaw/workspace/skills/fastmail-email/scripts/jmap_client.py send-attach wlweert@gmail.com "Parametric vs WOA23 Salinity Comparison" "Report comparing the parametric salinity baseline against WOA23 observational data for all 896 SSWD-EvoEpi model sites." report.pdf`

## Min/Max Pages
8-15 pages
