#!/bin/bash
# Launch LD50 load-dependent disease grid search
# 256 configs × 3 seeds = 768 runs, 8 parallel at K=1000
# Estimated: 768 / 8 = 96 waves × 2.2 hrs = ~211 hrs = ~8.8 days
#
# See specs/LD50_GRID_SEARCH_PLAN.md for full design
# Output: results/ld50_sweep/{scenario_name}/result_seed{N}.json

set -euo pipefail

REPO="/home/starbot/.openclaw/workspace/sswd-evoepi"
RUNNER="$REPO/experiments/calibration_runner.py"
CONFIG_DIR="$REPO/experiments/ld50_sweep/configs"
OUTPUT_BASE="$REPO/results/ld50_sweep"
YEARS=38
K=1000
K_REF=5000
SEEDS="42 137 256"
DISEASE_YEAR=1
MAX_JOBS=8

# Thread control — CRITICAL for parallel launches
# Each job must use exactly 1 thread to avoid oversubscription
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export NUMBA_NUM_THREADS=1
export NUMBA_THREADING_LAYER=omp

mkdir -p "$OUTPUT_BASE/logs"

# Collect all LD50 sweep configs
CONFIGS=("$CONFIG_DIR"/LD50_*.json)
N_CONFIGS=${#CONFIGS[@]}
N_SEEDS=0
for _ in $SEEDS; do N_SEEDS=$((N_SEEDS + 1)); done
TOTAL=$((N_CONFIGS * N_SEEDS))

MASTER_LOG="$OUTPUT_BASE/master.log"

echo "$(date '+%Y-%m-%d %H:%M:%S') LD50 Grid Search" | tee "$MASTER_LOG"
echo "  Configs: $N_CONFIGS, Seeds per config: $N_SEEDS, Total runs: $TOTAL" | tee -a "$MASTER_LOG"
echo "  K=$K, K_ref=$K_REF, disease_year=$DISEASE_YEAR, years=$YEARS" | tee -a "$MASTER_LOG"
echo "  Max parallel: $MAX_JOBS" | tee -a "$MASTER_LOG"
echo "  Output: $OUTPUT_BASE/" | tee -a "$MASTER_LOG"
echo "" | tee -a "$MASTER_LOG"

LAUNCHED=0
SKIPPED=0
RUNNING=0

for config in "${CONFIGS[@]}"; do
    scenario=$(basename "$config" .json)
    outdir="$OUTPUT_BASE/$scenario"

    for seed in $SEEDS; do
        result_file="$outdir/result_seed${seed}.json"

        # Skip already-completed runs
        if [ -f "$result_file" ]; then
            SKIPPED=$((SKIPPED + 1))
            continue
        fi

        # Maintain pool of MAX_JOBS concurrent jobs
        while [ "$RUNNING" -ge "$MAX_JOBS" ]; do
            wait -n 2>/dev/null || true
            RUNNING=$((RUNNING - 1))
        done

        mkdir -p "$outdir"
        logfile="$OUTPUT_BASE/logs/${scenario}_seed${seed}.log"

        echo "$(date '+%H:%M:%S') LAUNCH $scenario seed=$seed ($((LAUNCHED+1))/$TOTAL)"
        setsid python3 "$RUNNER" \
            --config "$config" \
            --seed "$seed" \
            --K "$K" \
            --K-ref "$K_REF" \
            --years "$YEARS" \
            --disease-year "$DISEASE_YEAR" \
            --output "$outdir" \
            > "$logfile" 2>&1 &

        RUNNING=$((RUNNING + 1))
        LAUNCHED=$((LAUNCHED + 1))
        echo "$(date '+%H:%M:%S') LAUNCH $scenario seed=$seed (${LAUNCHED}/${TOTAL}, running=$RUNNING)" >> "$MASTER_LOG"
    done
done

# Wait for all remaining jobs
echo ""
echo "$(date '+%H:%M:%S') All jobs submitted. Waiting for $RUNNING remaining..."
wait

echo "$(date '+%Y-%m-%d %H:%M:%S') DONE — launched=$LAUNCHED, skipped=$SKIPPED, total=$TOTAL" | tee -a "$MASTER_LOG"
echo "  Output: $OUTPUT_BASE/" | tee -a "$MASTER_LOG"
echo "  Logs: $OUTPUT_BASE/logs/" | tee -a "$MASTER_LOG"
echo ""
echo "Monitor with:"
echo "  ps aux | grep calibration_runner | grep -v grep | wc -l"
echo "  tail -f $MASTER_LOG"
echo "  ls $OUTPUT_BASE/*/result_seed*.json | wc -l   # completed count"
