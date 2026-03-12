#!/bin/bash
# Generate W163-W168 configs and launcher for Tier 2
# Usage: ./generate_W163_W168.sh <best_K_cv>
# Example: ./generate_W163_W168.sh 0.8
set -e

BEST_CV="${1:?Usage: $0 <best_K_cv>}"
DIR="$(cd "$(dirname "$0")" && pwd)"

echo "Generating Tier 2 configs with K_cv=${BEST_CV}..."

# W163: best_cv + K_half=800K + alpha_env=0.22
cat > "${DIR}/W163_config.json" << EOF
{
  "param_overrides": {
    "disease.K_half": 800000.0,
    "disease.P_env_dynamic": true,
    "disease.P_env_floor": 500.0,
    "disease.P_env_max": 2000.0,
    "disease.T_vbnc": 12.0,
    "disease.T_vbnc_initial": 12.0,
    "disease.T_vbnc_min": 9.0,
    "disease.activation_threshold": 50.0,
    "disease.alpha_env": 0.22,
    "disease.cumulative_dose_threshold": 1000.0,
    "disease.delta_env": 0.02,
    "disease.disease_origin_nodes": [322, 319, 632, 633, 634],
    "disease.k_vbnc": 2.0,
    "disease.pathogen_adapt_rate": 0.001,
    "disease.pathogen_adaptation": true,
    "disease.pathogen_revert_rate": 0.0,
    "disease.v_adapt_rate": 0.001,
    "disease.v_max_warm": 0.7,
    "disease.virulence_evolution": true,
    "disease.wavefront_D_P": 300.0,
    "disease.wavefront_D_P_max_range": 3000.0,
    "disease.wavefront_enabled": true,
    "population.settler_survival": 0.001,
    "spatial.n_connectivity": 0.3,
    "spatial.alpha_self_open": 0.02,
    "spatial.alpha_self_fjord": 0.7
  },
  "K": 5000,
  "K_cv": ${BEST_CV},
  "years": 13,
  "disease_year": 1,
  "seeds": [42, 123, 999],
  "sst_start_year": 2012,
  "notes": "Tier2: K_cv=${BEST_CV}, K_half=800K, alpha_env=0.22 (stronger env accumulation)"
}
EOF

# W164: best_cv + K_half=800K + alpha_env=0.25
cat > "${DIR}/W164_config.json" << EOF
{
  "param_overrides": {
    "disease.K_half": 800000.0,
    "disease.P_env_dynamic": true,
    "disease.P_env_floor": 500.0,
    "disease.P_env_max": 2000.0,
    "disease.T_vbnc": 12.0,
    "disease.T_vbnc_initial": 12.0,
    "disease.T_vbnc_min": 9.0,
    "disease.activation_threshold": 50.0,
    "disease.alpha_env": 0.25,
    "disease.cumulative_dose_threshold": 1000.0,
    "disease.delta_env": 0.02,
    "disease.disease_origin_nodes": [322, 319, 632, 633, 634],
    "disease.k_vbnc": 2.0,
    "disease.pathogen_adapt_rate": 0.001,
    "disease.pathogen_adaptation": true,
    "disease.pathogen_revert_rate": 0.0,
    "disease.v_adapt_rate": 0.001,
    "disease.v_max_warm": 0.7,
    "disease.virulence_evolution": true,
    "disease.wavefront_D_P": 300.0,
    "disease.wavefront_D_P_max_range": 3000.0,
    "disease.wavefront_enabled": true,
    "population.settler_survival": 0.001,
    "spatial.n_connectivity": 0.3,
    "spatial.alpha_self_open": 0.02,
    "spatial.alpha_self_fjord": 0.7
  },
  "K": 5000,
  "K_cv": ${BEST_CV},
  "years": 13,
  "disease_year": 1,
  "seeds": [42, 123, 999],
  "sst_start_year": 2012,
  "notes": "Tier2: K_cv=${BEST_CV}, K_half=800K, alpha_env=0.25 (even stronger)"
}
EOF

# W165: best_cv + K_half=1.2M + alpha_env=0.22
cat > "${DIR}/W165_config.json" << EOF
{
  "param_overrides": {
    "disease.K_half": 1200000.0,
    "disease.P_env_dynamic": true,
    "disease.P_env_floor": 500.0,
    "disease.P_env_max": 2000.0,
    "disease.T_vbnc": 12.0,
    "disease.T_vbnc_initial": 12.0,
    "disease.T_vbnc_min": 9.0,
    "disease.activation_threshold": 50.0,
    "disease.alpha_env": 0.22,
    "disease.cumulative_dose_threshold": 1000.0,
    "disease.delta_env": 0.02,
    "disease.disease_origin_nodes": [322, 319, 632, 633, 634],
    "disease.k_vbnc": 2.0,
    "disease.pathogen_adapt_rate": 0.001,
    "disease.pathogen_adaptation": true,
    "disease.pathogen_revert_rate": 0.0,
    "disease.v_adapt_rate": 0.001,
    "disease.v_max_warm": 0.7,
    "disease.virulence_evolution": true,
    "disease.wavefront_D_P": 300.0,
    "disease.wavefront_D_P_max_range": 3000.0,
    "disease.wavefront_enabled": true,
    "population.settler_survival": 0.001,
    "spatial.n_connectivity": 0.3,
    "spatial.alpha_self_open": 0.02,
    "spatial.alpha_self_fjord": 0.7
  },
  "K": 5000,
  "K_cv": ${BEST_CV},
  "years": 13,
  "disease_year": 1,
  "seeds": [42, 123, 999],
  "sst_start_year": 2012,
  "notes": "Tier2: K_cv=${BEST_CV}, K_half=1.2M, alpha_env=0.22 (K_half lift + alpha)"
}
EOF

# W166: best_cv + K_half=1.2M + alpha_env=0.25
cat > "${DIR}/W166_config.json" << EOF
{
  "param_overrides": {
    "disease.K_half": 1200000.0,
    "disease.P_env_dynamic": true,
    "disease.P_env_floor": 500.0,
    "disease.P_env_max": 2000.0,
    "disease.T_vbnc": 12.0,
    "disease.T_vbnc_initial": 12.0,
    "disease.T_vbnc_min": 9.0,
    "disease.activation_threshold": 50.0,
    "disease.alpha_env": 0.25,
    "disease.cumulative_dose_threshold": 1000.0,
    "disease.delta_env": 0.02,
    "disease.disease_origin_nodes": [322, 319, 632, 633, 634],
    "disease.k_vbnc": 2.0,
    "disease.pathogen_adapt_rate": 0.001,
    "disease.pathogen_adaptation": true,
    "disease.pathogen_revert_rate": 0.0,
    "disease.v_adapt_rate": 0.001,
    "disease.v_max_warm": 0.7,
    "disease.virulence_evolution": true,
    "disease.wavefront_D_P": 300.0,
    "disease.wavefront_D_P_max_range": 3000.0,
    "disease.wavefront_enabled": true,
    "population.settler_survival": 0.001,
    "spatial.n_connectivity": 0.3,
    "spatial.alpha_self_open": 0.02,
    "spatial.alpha_self_fjord": 0.7
  },
  "K": 5000,
  "K_cv": ${BEST_CV},
  "years": 13,
  "disease_year": 1,
  "seeds": [42, 123, 999],
  "sst_start_year": 2012,
  "notes": "Tier2: K_cv=${BEST_CV}, K_half=1.2M, alpha_env=0.25 (K_half lift + stronger alpha)"
}
EOF

# W167: best_cv + K_half=1.2M + alpha_env=0.18 + delta_env=0.01 (slower decay)
cat > "${DIR}/W167_config.json" << EOF
{
  "param_overrides": {
    "disease.K_half": 1200000.0,
    "disease.P_env_dynamic": true,
    "disease.P_env_floor": 500.0,
    "disease.P_env_max": 2000.0,
    "disease.T_vbnc": 12.0,
    "disease.T_vbnc_initial": 12.0,
    "disease.T_vbnc_min": 9.0,
    "disease.activation_threshold": 50.0,
    "disease.alpha_env": 0.18,
    "disease.cumulative_dose_threshold": 1000.0,
    "disease.delta_env": 0.01,
    "disease.disease_origin_nodes": [322, 319, 632, 633, 634],
    "disease.k_vbnc": 2.0,
    "disease.pathogen_adapt_rate": 0.001,
    "disease.pathogen_adaptation": true,
    "disease.pathogen_revert_rate": 0.0,
    "disease.v_adapt_rate": 0.001,
    "disease.v_max_warm": 0.7,
    "disease.virulence_evolution": true,
    "disease.wavefront_D_P": 300.0,
    "disease.wavefront_D_P_max_range": 3000.0,
    "disease.wavefront_enabled": true,
    "population.settler_survival": 0.001,
    "spatial.n_connectivity": 0.3,
    "spatial.alpha_self_open": 0.02,
    "spatial.alpha_self_fjord": 0.7
  },
  "K": 5000,
  "K_cv": ${BEST_CV},
  "years": 13,
  "disease_year": 1,
  "seeds": [42, 123, 999],
  "sst_start_year": 2012,
  "notes": "Tier2: K_cv=${BEST_CV}, K_half=1.2M, alpha_env=0.18, delta_env=0.01 (slower pathogen decay, half-life 69d)"
}
EOF

# W168: best_cv + K_half=1.5M + alpha_env=0.22
cat > "${DIR}/W168_config.json" << EOF
{
  "param_overrides": {
    "disease.K_half": 1500000.0,
    "disease.P_env_dynamic": true,
    "disease.P_env_floor": 500.0,
    "disease.P_env_max": 2000.0,
    "disease.T_vbnc": 12.0,
    "disease.T_vbnc_initial": 12.0,
    "disease.T_vbnc_min": 9.0,
    "disease.activation_threshold": 50.0,
    "disease.alpha_env": 0.22,
    "disease.cumulative_dose_threshold": 1000.0,
    "disease.delta_env": 0.02,
    "disease.disease_origin_nodes": [322, 319, 632, 633, 634],
    "disease.k_vbnc": 2.0,
    "disease.pathogen_adapt_rate": 0.001,
    "disease.pathogen_adaptation": true,
    "disease.pathogen_revert_rate": 0.0,
    "disease.v_adapt_rate": 0.001,
    "disease.v_max_warm": 0.7,
    "disease.virulence_evolution": true,
    "disease.wavefront_D_P": 300.0,
    "disease.wavefront_D_P_max_range": 3000.0,
    "disease.wavefront_enabled": true,
    "population.settler_survival": 0.001,
    "spatial.n_connectivity": 0.3,
    "spatial.alpha_self_open": 0.02,
    "spatial.alpha_self_fjord": 0.7
  },
  "K": 5000,
  "K_cv": ${BEST_CV},
  "years": 13,
  "disease_year": 1,
  "seeds": [42, 123, 999],
  "sst_start_year": 2012,
  "notes": "Tier2: K_cv=${BEST_CV}, K_half=1.5M, alpha_env=0.22 (much higher K_half)"
}
EOF

echo "Generated 6 configs: W163-W168"
echo ""
echo "Parameter summary:"
echo "  W163: K_cv=${BEST_CV}, K_half=800K,  α=0.22, δ=0.02"
echo "  W164: K_cv=${BEST_CV}, K_half=800K,  α=0.25, δ=0.02"
echo "  W165: K_cv=${BEST_CV}, K_half=1.2M,  α=0.22, δ=0.02"
echo "  W166: K_cv=${BEST_CV}, K_half=1.2M,  α=0.25, δ=0.02"
echo "  W167: K_cv=${BEST_CV}, K_half=1.2M,  α=0.18, δ=0.01 (slow decay)"
echo "  W168: K_cv=${BEST_CV}, K_half=1.5M,  α=0.22, δ=0.02"

echo ""
echo "NOTE: launch_W163_W168.sh is a smart rolling launcher (not overwritten by this script)."
echo "It launches W163-W166 first, waits for disease crash, then launches W167-W168."
echo ""
echo "To deploy:"
echo "  1. scp W163-W168 configs + launcher to Xeon"
echo "  2. ssh xeon 'cd ~/projects/sswd-evoepi && bash experiments/calibration/launch_W163_W168.sh'"
