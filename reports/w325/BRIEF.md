# Report Brief — W325 Single-Run Analysis

- **Title**: W325 Calibration Analysis — Best RMSLE Configuration
- **Audience**: Willem (marine biology PhC) — direct tone, use jargon freely
- **Data sources**:
  - Result JSON: `/home/starbot/.openclaw/workspace/sswd-evoepi/results/k1000_scaled_sweep/W325/result_seed42.json`
  - Combined JSON: `/home/starbot/.openclaw/workspace/sswd-evoepi/results/k1000_scaled_sweep/W325/combined_results.json`
  - Monthly NPZ: `/home/starbot/.openclaw/workspace/sswd-evoepi/results/k1000_scaled_sweep/W325/monthly_seed42.npz`
  - Config JSON: `/home/starbot/.openclaw/workspace/sswd-evoepi/experiments/calibration/W325_config.json`
  - Run log: `/home/starbot/.openclaw/workspace/sswd-evoepi/results/k1000_scaled_sweep/W325/run.log`
- **Key questions to answer**:
  1. How does the latitudinal gradient look? Is the southern suppression working?
  2. What does the wavefront propagation pattern look like?
  3. How do the gradient-breaking mechanisms (δ_env(T) + P_env-gated recovery) affect dynamics?
  4. What's the evolutionary response — is resistance still dominant?
  5. Where are the remaining calibration gaps (BC-N=4.3% vs 20% target, JDF=8.2% vs 2%)?
  6. What parameter adjustments should the next sweep target?
- **Structure**: Single-run dynamics report (not multi-config sweep)
  - Executive Summary (1 page)
  - Configuration & Parameters
  - Spatial-Temporal Dynamics (hero section — heatmaps, wavefront, trajectories)
  - Evolutionary Dynamics (resistance, tolerance, recovery, pathogen evolution)
  - Calibration Performance (recovery vs targets, per-region scoring)
  - Conclusions & Next Steps
- **Genetics available**: YES — yearly_mean_resistance, yearly_mean_tolerance, yearly_mean_recovery, yearly_va_*, final_mean_T_vbnc, final_mean_v_local all present
- **Key W325 parameters**:
  - K_half=800K (reduced from 1.5M), T_vbnc_min=10°C, wavefront_D_P=300
  - delta_env_T_dependent=true, delta_env_cold=0.05, delta_env_warm=0.01
  - recovery_P_env_gated=true, recovery_P_env_half=200, recovery_P_env_hill=2
  - alpha_self_open=0.05, alpha_self_fjord=0.5
  - K=1000 with K_ref=5000 density scaling (density_scale=5.0)
- **RMSLE**: 0.348 (previous best: 0.655 W294)
- **Key results**: AK-PWS=29.9% (target 20%), AK-FN=20.4% (target 20%), BC-N=4.3% (target 20%), SS-S=8.6% (target 5%), JDF=8.2% (target 2%), OR=0.4% (target 0.25%), CA-N=0.1% (target 0.1%)
- **Min/max pages**: 15–25 (single run, not sweep)
- **Delivery**: email to wlweert@gmail.com with subject "W325 Analysis — RMSLE 0.348"

## Calibration Targets (for reference)
- AK-PWS: 20%, AK-FN: 20%, AK-FS: 20%, BC-N: 20%, SS-S: 5%, JDF: 2%, OR: 0.25%, CA-N: 0.1%
- 6/8 within 2×, 8/8 within 5×

## Coastline Order (S→N)
BJ → CA-S → CA-C → CA-N → OR → WA-O → JDF → SS-S → SS-N → BC-C → BC-N → AK-FS → AK-FN → AK-OC → AK-PWS → AK-EG → AK-WG → AK-AL
