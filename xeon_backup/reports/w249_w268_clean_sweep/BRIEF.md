# Report Brief

- **Title**: W249–W268 Clean Sweep: Post-Audit Calibration Results
- **Audience**: Willem Weertman — marine biologist / eco-evo modeler. Wants to understand how the 4 bug fixes changed model dynamics and which parameter combinations now look most promising.
- **Data sources**: 
  - Results: `results/calibration/W249-W268/W{249..268}/result_seed42.json`
  - Configs: `experiments/calibration/W{249..268}_config.json`
  - Repo: `/home/starbot/.openclaw/workspace/sswd-evoepi`
- **Key questions to answer**:
  1. How did the 4 audit bug fixes change the baseline dynamics? (W249 vs historical W201)
  2. Which parameter combination performs best? (W257: fw=25, K_half=1.2M)
  3. Is the AK recovery gap closing? What's the trend?
  4. Which parameters matter most in this sweep?
  5. Where should the next sweep focus?
- **Figures needed**:
  1. RMSE ranking bar chart (all 20 configs)
  2. Regional recovery heatmap (configs × scored regions)
  3. AK-PWS vs CA-N latitude tradeoff scatter
  4. Parameter sensitivity panels (fw_strength, K_half, alpha_env, n_conn, r_total)
  5. Recovery bar chart: W257 vs targets
- **Structure**:
  1. Executive Summary
  2. Sweep Design (configs table, parameters explored)
  3. Results Overview (RMSE ranking, regional heatmap)
  4. Bug Fix Impact (W249 vs historical W201)
  5. Parameter Analysis (fw_strength, K_half, alpha_env, interactions)
  6. Best Configuration Deep Dive (W257)
  7. Discussion & Next Steps
- **Min/max pages**: 8-15
- **Delivery**: email to wlweert@gmail.com with subject "W249-W268 Clean Sweep Report"

## Config Descriptions
- W249: Baseline (fw=0, W201 params)
- W250: fw=15
- W251: fw=20
- W252: fw=25
- W253: fw=30
- W254: fw=25, α_env=0.22
- W255: fw=25, α_env=0.25
- W256: fw=25, K½=400K
- W257: fw=25, K½=1.2M ★ BEST
- W258: fw=25, α=0.25, K½=400K
- W259: fw=25, α=0.25, K½=1.2M
- W260: fw=25, n_conn=0.2
- W261: fw=25, n_conn=0.5
- W262: fw=25, α=0.25, n_conn=0.2
- W263: fw=25, α=0.25, n_conn=0.5
- W264: fw=25, α_fjord=0.5
- W265: fw=25, α_fjord=0.9
- W266: fw=25, r_total=0.002
- W267: fw=25, r_total=0.004
- W268: fw=25, α=0.25, r_total=0.002

## Config Groups
- baseline: W249
- fw_sweep: W250, W251, W252, W253
- alpha_env: W254, W255
- k_half: W256, W257
- alpha_x_k: W258, W259
- n_conn: W260, W261, W262, W263
- alpha_fjord: W264, W265
- r_total: W266, W267, W268

## Calibration Targets (from sswd_evoepi.metrics)
- AK-PWS: 30%, AK-FN: 30%, AK-FS: 20%, BC-N: 20%
- SS-S: 5%, JDF: 2%, OR: 0.25%, CA-N: 0.1%

## Historical Reference
- W201 (pre-audit): RMSE=0.663, AK-PWS=3.1%
- W173 (best pre-salinity): RMSE=0.608, AK-PWS=10.4%
- W217 (best spatial levers): RMSE=0.640, AK-PWS=3.7%
- W233 (best salinity): RMSE=1.011, AK-PWS=3.1%
