#!/bin/bash
# Launch W155-W162 calibration batch (Tier 1: K_cv sweep)
# 8 configs × 3 seeds = 24 runs at 24 parallel
set -e

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1
export NUMBA_NUM_THREADS=1
export NUMBA_THREADING_LAYER=omp
export PYTHONUNBUFFERED=1

REPO="/home/weertman/projects/sswd-evoepi"
RUNNER="$REPO/experiments/calibration_runner.py"
CONFIG_DIR="$REPO/experiments/calibration"
RESULT_DIR="$REPO/results/calibration"
LOGDIR="$REPO/logs/W155_W162"

mkdir -p "$LOGDIR"
mkdir -p "$RESULT_DIR"

MAX_PARALLEL=24
RUNNING=0

echo "=== W155-W162 Calibration Batch ==="
echo "Start: $(date)"
echo "Max parallel: $MAX_PARALLEL"
echo ""

for W in W155 W156 W157 W158 W159 W160 W161 W162; do
    CONFIG="$CONFIG_DIR/${W}_config.json"
    
    if [ ! -f "$CONFIG" ]; then
        echo "SKIP: $CONFIG not found"
        continue
    fi
    
    # Read K_cv from config
    K_CV=$(python3 -c "import json; d=json.load(open('$CONFIG')); print(d.get('K_cv', 0.0))")
    
    for SEED in 42 123 999; do
        OUTDIR="$RESULT_DIR/$W"
        mkdir -p "$OUTDIR"
        
        # Skip if result already exists
        if [ -f "$OUTDIR/result_seed${SEED}.json" ]; then
            echo "SKIP: $W seed $SEED already complete"
            continue
        fi
        
        LOG="$LOGDIR/${W}_seed${SEED}.log"
        
        echo "LAUNCH: $W seed=$SEED K_cv=$K_CV"
        
        python3 "$RUNNER" \
            --config "$CONFIG" \
            --seed "$SEED" \
            --K 5000 \
            --K-cv "$K_CV" \
            --years 13 \
            --disease-year 1 \
            --output "$OUTDIR" \
            > "$LOG" 2>&1 &
        
        RUNNING=$((RUNNING + 1))
        
        # Throttle
        if [ "$RUNNING" -ge "$MAX_PARALLEL" ]; then
            echo "Waiting for batch ($RUNNING running)..."
            wait -n 2>/dev/null || true
            RUNNING=$((RUNNING - 1))
        fi
        
        # Small delay to stagger Numba compilation
        sleep 2
    done
done

echo ""
echo "All jobs launched. Waiting for completion..."
wait
echo "=== BATCH COMPLETE: $(date) ==="
