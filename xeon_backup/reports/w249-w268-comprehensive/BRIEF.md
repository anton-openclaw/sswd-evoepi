# Report Brief — W249-W268 Comprehensive Calibration Report

## Title
W249-W268 Post-Audit Clean Sweep — Calibration Report

## Audience
Willem Weertman — PhD candidate in marine biology/ecology. Understands ABMs, population genetics, disease ecology. Cares about: which parameters move the needle, AK recovery gap, next steps.

## Sweep Context
First calibration sweep after the March 14 codebase audit that fixed 4 bugs:
- Immunosuppression double-decrement
- Ghost cascade induction
- K overshoot in recruitment
- Growth noise 27× too small

20 configurations × 1 seed (42), K=5000, 13 years (2012-2025), 896 nodes.

## Config Groups & Descriptions
| Config | Group | Description |
|--------|-------|-------------|
| W249 | baseline | Pure W201 baseline (fw=0, post-audit) |
| W250 | fw_sweep | fw_strength=15 |
| W251 | fw_sweep | fw_strength=20 |
| W252 | fw_sweep | fw_strength=25 |
| W253 | fw_sweep | fw_strength=30 |
| W254 | alpha_env | alpha_env=0.22, fw=25 |
| W255 | alpha_env | alpha_env=0.25, fw=25 |
| W256 | K_half | K_half=400K, fw=25 |
| W257 | K_half | K_half=1.2M, fw=25 (★BEST) |
| W258 | cross | K_half=400K, alpha_env=0.25, fw=25 |
| W259 | cross | K_half=1.2M, alpha_env=0.25, fw=25 |
| W260 | n_conn | n_connectivity=0.2, fw=25 |
| W261 | n_conn | n_connectivity=0.5, fw=25 |
| W262 | n_conn_alpha | n_conn=0.2, alpha_env=0.25, fw=25 |
| W263 | n_conn_alpha | n_conn=0.5, alpha_env=0.25, fw=25 |
| W264 | alpha_self | alpha_self_fjord=0.5, fw=25 |
| W265 | alpha_self | alpha_self_fjord=0.9, fw=25 |
| W266 | r_total | r_total=0.002, fw=25 |
| W267 | r_total | r_total=0.004, fw=25 |
| W268 | r_total_alpha | r_total=0.002, alpha_env=0.25, fw=25 |

## Group Colors
- baseline: #607D8B (gray)
- fw_sweep: #2196F3 (blue)
- alpha_env: #FF9800 (orange)
- K_half: #4CAF50 (green)
- cross: #9C27B0 (purple)
- n_conn: #00BCD4 (cyan)
- n_conn_alpha: #009688 (teal)
- alpha_self: #795548 (brown)
- r_total: #F44336 (red)
- r_total_alpha: #E91E63 (pink)

## Data Sources
- Results: reports/w249-w268-comprehensive/data/sweep/W{249..268}/result_seed42.json
- Monthly NPZ (top configs): data/sweep/W{249,253,257,259,261,267}/monthly_seed42.npz
- Configs: reports/w249-w268-comprehensive/data/configs/W{249..268}_config.json

## Key Questions
1. What is the impact of the 4 bug fixes? (W249 vs historical W201 RMSE=0.663)
2. Which parameter is most important? (K_half=1.2M is the big lever)
3. Can we close the AK recovery gap? (12.1% vs 30% target)
4. What is the parameter interaction structure? (K_half × alpha_env, K_half × fw_strength)
5. Are southern regions still overshooting? (CA-N: 1.22% vs 0.1% target)

## Calibration Targets
AK-PWS: 30%, AK-FN: 30%, AK-FS: 20%, BC-N: 20%, SS-S: 5%, JDF: 2%, OR: 0.25%, CA-N: 0.1%

## Historical Comparison
- W201 (pre-audit baseline): RMSE=0.663
- W173 (best pre-salinity): RMSE=0.608, AK-PWS=10.4%
- W233 (best salinity): RMSE=1.011, AK-PWS=3.1%
- W257 (post-audit best): RMSE=0.697, AK-PWS=12.1%

## Report Structure
Follow template spec at skills/report-generator/REPORT_TEMPLATE_SPEC.md sections 1-7 + appendices.
Skip sections 8-10 (disease details, sensitivity, comparison) — not enough data for those.

Target: 20-30 pages, 15-20 figures.

## Delivery
Email to wlweert@gmail.com with subject "W249-W268 Calibration Report — Post-Audit Clean Sweep"
