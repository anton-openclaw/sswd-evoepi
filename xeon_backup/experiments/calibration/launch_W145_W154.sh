#!/bin/bash
# Launch W145-W154: Drill into W142 sweet spot (n_conn=0.5, alpha=0.20)
set -e
cd ~/projects/sswd-evoepi

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1
export NUMBA_NUM_THREADS=1
export NUMBA_THREADING_LAYER=omp
export PYTHONUNBUFFERED=1

RESULTS_BASE="results/calibration"
mkdir -p "${RESULTS_BASE}/W145-W154"
MAX_PARALLEL=24
RUNNING=0

for W in $(seq 145 154); do
    CONFIG="experiments/calibration/W${W}_config.json"
    if [ ! -f "$CONFIG" ]; then
        echo "SKIP: $CONFIG not found"
        continue
    fi
    
    for SEED in 42 123 999; do
        OUTDIR="${RESULTS_BASE}/W${W}"
        LOGFILE="${RESULTS_BASE}/W145-W154/W${W}_seed${SEED}.log"
        
        echo "Launching W${W} seed=${SEED}"
        
        nohup python3 experiments/calibration_runner.py \
            --config "$CONFIG" \
            --output "$OUTDIR" \
            --seed "$SEED" \
            > "$LOGFILE" 2>&1 &
        
        RUNNING=$((RUNNING + 1))
        
        if [ $RUNNING -ge $MAX_PARALLEL ]; then
            echo "Hit $MAX_PARALLEL parallel, waiting..."
            wait -n
            RUNNING=$((RUNNING - 1))
        fi
    done
done

echo "All 30 runs launched."
