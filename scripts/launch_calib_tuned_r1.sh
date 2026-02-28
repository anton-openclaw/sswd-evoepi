#!/bin/bash
# Calibration Tuning Round 1 — 5 parameter sets × 3 seeds = 15 runs
# Based on tuning agent recommendations from round 0-4 results
# Key insight: K_half=200K for Alaska + T_vbnc shift + P_env for gradient

set -e

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1
export NUMBA_THREADING_LAYER=omp
export NUMBA_NUM_THREADS=8   # 15 processes × 8 = 120 threads on 128 cores
export PYTHONUNBUFFERED=1     # Real-time log output!

cd ~/projects/sswd-evoepi

SEEDS="42 123 999"

# Set 06: K_half=200K + T_vbnc=14 (conservative — T_vbnc shift only)
# Set 07: K_half=200K + P_env=3000 (high P_env at high K_half)
# Set 08: K_half=200K + T_vbnc=14 + P_env=3000 (core combination)
# Set 09: K_half=150K + T_vbnc=14 + P_env=2500 (moderate aggressive)
# Set 10: K_half=120K + T_vbnc=15 + P_env=4000 (aggressive — all levers)

declare -A CONFIGS
CONFIGS[round_06]='{"disease.K_half": 200000, "disease.T_vbnc": 14}'
CONFIGS[round_07]='{"disease.K_half": 200000, "disease.P_env_max": 3000}'
CONFIGS[round_08]='{"disease.K_half": 200000, "disease.T_vbnc": 14, "disease.P_env_max": 3000}'
CONFIGS[round_09]='{"disease.K_half": 150000, "disease.T_vbnc": 14, "disease.P_env_max": 2500}'
CONFIGS[round_10]='{"disease.K_half": 120000, "disease.T_vbnc": 15, "disease.P_env_max": 4000}'

echo "=== Calibration Tuning Round 1 ==="
echo "Launching 5 parameter sets × 3 seeds = 15 runs"
echo "Time: $(date)"
echo ""

for round in round_06 round_07 round_08 round_09 round_10; do
    OUTDIR="results/calibration/${round}"
    mkdir -p "$OUTDIR"
    CONFIG_FILE="$OUTDIR/config.json"
    echo "${CONFIGS[$round]}" > "$CONFIG_FILE"
    echo "Config for ${round}: ${CONFIGS[$round]}"
    
    for seed in $SEEDS; do
        LOG="$OUTDIR/seed_${seed}.log"
        nohup python3 experiments/calibration_runner.py \
            --config "$CONFIG_FILE" \
            --seed $seed \
            --output "$OUTDIR" \
            --K 5000 \
            --years 13 \
            --disease-year 1 \
            > "$LOG" 2>&1 &
        echo "  Launched ${round} seed=${seed} (PID $!)"
    done
done

echo ""
echo "All 15 launched at $(date)"
echo "PIDs:"
pgrep -f calibration_runner | head -20
echo ""
echo "Expected: ~25-30h runtime (15 parallel, memory contention)"
echo "Monitor: tail -f results/calibration/round_08/seed_42.log"
