#!/bin/bash
# Launch W125-W134 calibration runs — v2 with NUMBA_NUM_THREADS fix
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
MAX_PARALLEL=24
RUNNING=0

# Clean old output files from failed runs
for W in $(seq 125 134); do
    rm -f "${RESULTS_BASE}/W${W}/stdout_seed"*.log 2>/dev/null
    rm -f "${RESULTS_BASE}/W${W}/checkpoint_seed"*.json 2>/dev/null
    rm -f "${RESULTS_BASE}/W${W}/monthly_seed"*.npz 2>/dev/null
    rm -f "${RESULTS_BASE}/W${W}/combined_results.json" 2>/dev/null
done

for W in $(seq 125 134); do
    CONFIG="experiments/calibration/W${W}_config.json"
    if [ ! -f "$CONFIG" ]; then
        echo "SKIP: $CONFIG not found"
        continue
    fi
    
    for SEED in 42 123 999; do
        OUTDIR="${RESULTS_BASE}/W${W}"
        LOGFILE="${RESULTS_BASE}/W125-W134/W${W}_seed${SEED}.log"
        
        mkdir -p "$OUTDIR"
        
        echo "$(date): Launching W${W} seed=${SEED}"
        
        nohup python3 experiments/calibration_runner.py \
            --config "$CONFIG" \
            --output "$OUTDIR" \
            --seed "$SEED" \
            > "$LOGFILE" 2>&1 &
        
        RUNNING=$((RUNNING + 1))
        
        if [ $RUNNING -ge $MAX_PARALLEL ]; then
            echo "Hit $MAX_PARALLEL parallel, waiting for one to finish..."
            wait -n
            RUNNING=$((RUNNING - 1))
        fi
    done
done

echo ""
echo "$(date): All 30 runs launched."
echo "Monitor: ps aux | grep calibration_runner | grep -v grep | wc -l"

# Wait for all to complete
wait
echo "$(date): All runs complete!"
