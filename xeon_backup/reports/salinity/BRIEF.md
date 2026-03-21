# Report Brief — Seasonal Salinity Mechanism

- **Title**: Seasonal Salinity as a Latitude-Dependent Disease Suppression Mechanism in the SSWD-EvoEpi Model
- **Audience**: Willem Weertman (PhD candidate, marine biology / larval ecology). Knows the model well. Wants to understand the data behind the salinity mechanism, the validation, and what it means for calibration.
- **Data sources**:
  - DFO lighthouse climatology: `salinity_validation/dfo_monthly_climatology.csv`
  - Site enclosedness: `data/nodes/site_enclosedness.json`
  - All 896 site definitions: `data/nodes/all_sites.json`
  - Salinity implementation: `sswd_evoepi/salinity.py`
  - Salinity viz: `sswd_evoepi/viz/salinity.py`
  - Implementation sketch: `salinity_implementation_sketch.md`
- **Key questions to answer**:
  1. What's the observational basis for the salinity model? (DFO lighthouse data)
  2. How does fjord_depth_norm create the AK-CA asymmetry?
  3. What's the predicted suppression at different fw_strength values?
  4. How does the seasonal melt pulse create temporal selectivity?
  5. Where are the data gaps / limitations?
- **Figures needed**:
  1. Mechanism components (melt pulse, latitude factor, combined profiles)
  2. Salinity modifier transfer function (sal_mod curve)
  3. Fjord depth distribution by region (box plot)
  4. Regional salinity profiles (seasonal curves)
  5. Latitude asymmetry diagnostic (depression + suppression by lat)
  6. fw_strength sensitivity (dose-response)
  7. Regional suppression bars (N→S gradient)
  8. Suppression monthly panels (spatial snapshots)
  9. DFO validation (model vs observed)
  10. Depression heatmap (nodes × months)
- **Structure**:
  1. Executive Summary
  2. Motivation (latitude gradient problem)
  3. Observational Basis (DFO data, fjord_depth_norm)
  4. Model Design (equations, parameters)
  5. Spatial Predictions (asymmetry, regional patterns)
  6. Sensitivity Analysis (fw_strength, fw_depth_exp)
  7. Validation (DFO comparison, limitations)
  8. Implications for Calibration
- **Min/max pages**: 8-15
- **Delivery**: email to wlweert@gmail.com with subject "Seasonal Salinity Mechanism — Analysis Report"
