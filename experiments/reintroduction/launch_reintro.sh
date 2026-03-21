#!/bin/bash
# Launch reintroduction experiments (40 scenarios, 8 concurrent)
# All use W330 best-calibrated baseline + release events at year 13
# 38 years total, disease at year 1 (2013)

set -euo pipefail

REPO="/home/starbot/.openclaw/workspace/sswd-evoepi"
RUNNER="$REPO/experiments/calibration_runner.py"
CONFIG_DIR="$REPO/experiments/reintroduction/configs"
OUTPUT_BASE="$REPO/results/reintro_2050"
YEARS=38
K=1000
K_REF=5000
SEED=42
DISEASE_YEAR=1
MAX_JOBS=8

# Thread control — CRITICAL for parallel launches
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export NUMBA_NUM_THREADS=1
export NUMBA_THREADING_LAYER=omp

mkdir -p "$OUTPUT_BASE/logs"

# Collect all scenario configs
CONFIGS=("$CONFIG_DIR"/REINTRO_*.json)
TOTAL=${#CONFIGS[@]}

echo "$(date '+%H:%M:%S') Starting reintroduction experiments ($TOTAL scenarios, max $MAX_JOBS concurrent)"
echo "  K=$K, K_ref=$K_REF, seed=$SEED, disease_year=$DISEASE_YEAR, years=$YEARS"
echo ""

LAUNCHED=0
SKIPPED=0
RUNNING=0

for config in "${CONFIGS[@]}"; do
    scenario=$(basename "$config" .json)
    outdir="$OUTPUT_BASE/$scenario"
    logfile="$OUTPUT_BASE/logs/${scenario}.log"

    # Skip already-completed scenarios
    if [ -f "$outdir/result_seed${SEED}.json" ]; then
        echo "$(date '+%H:%M:%S') SKIP $scenario — already complete"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    # Maintain pool of MAX_JOBS concurrent jobs
    while [ "$RUNNING" -ge "$MAX_JOBS" ]; do
        wait -n 2>/dev/null || true
        RUNNING=$((RUNNING - 1))
    done

    mkdir -p "$outdir"

    echo "$(date '+%H:%M:%S') LAUNCH $scenario"
    setsid python3 "$RUNNER" \
        --config "$config" \
        --seed "$SEED" \
        --K "$K" \
        --K-ref "$K_REF" \
        --years "$YEARS" \
        --disease-year "$DISEASE_YEAR" \
        --output "$outdir" \
        > "$logfile" 2>&1 &

    RUNNING=$((RUNNING + 1))
    LAUNCHED=$((LAUNCHED + 1))
    echo "$(date '+%H:%M:%S') LAUNCH $scenario (${LAUNCHED}/${TOTAL}, running=$RUNNING)" >> "$OUTPUT_BASE/master.log"
done

# Wait for all remaining jobs
echo ""
echo "$(date '+%H:%M:%S') All jobs submitted. Waiting for $RUNNING remaining..."
wait

echo "$(date '+%H:%M:%S') DONE — launched=$LAUNCHED, skipped=$SKIPPED, total=$TOTAL"
echo "  Output: $OUTPUT_BASE/"
echo "  Logs: $OUTPUT_BASE/logs/"
echo ""
echo "Monitor with: ps aux | grep calibration_runner | grep -v grep"
