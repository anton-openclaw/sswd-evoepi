#!/bin/bash
# Launch W125-W134 calibration runs on Xeon
# 10 configs × 3 seeds = 30 runs, max 24 parallel

set -e
cd ~/projects/sswd-evoepi

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1
export NUMBA_THREADING_LAYER=omp
export PYTHONUNBUFFERED=1

RESULTS_BASE="results/calibration"
MAX_PARALLEL=24
RUNNING=0

for W in $(seq 125 134); do
    CONFIG="experiments/calibration/W${W}_config.json"
    if [ ! -f "$CONFIG" ]; then
        echo "SKIP: $CONFIG not found"
        continue
    fi
    
    for SEED in 42 123 999; do
        OUTDIR="${RESULTS_BASE}/W${W}"
        LOGFILE="${RESULTS_BASE}/W125-W134/W${W}_seed${SEED}.log"
        
        echo "Launching W${W} seed=${SEED} → ${OUTDIR}"
        
        nohup python3 experiments/calibration_runner.py \
            --config "$CONFIG" \
            --output "$OUTDIR" \
            --seed "$SEED" \
            > "$LOGFILE" 2>&1 &
        
        RUNNING=$((RUNNING + 1))
        
        # Throttle to MAX_PARALLEL
        if [ $RUNNING -ge $MAX_PARALLEL ]; then
            echo "Hit $MAX_PARALLEL parallel, waiting for one to finish..."
            wait -n
            RUNNING=$((RUNNING - 1))
        fi
    done
done

echo ""
echo "All 30 runs launched. PIDs:"
jobs -l
echo ""
echo "Monitor with: ps aux | grep calibration_runner | grep -v grep | wc -l"
echo "Logs in: ${RESULTS_BASE}/W125-W134/"
