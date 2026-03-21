#!/bin/bash
# Launch W185-W198 recalibration sweep (14 configs × 3 seeds = 42 runs)
# With genotyping optimization (e97e7eb)

cd /home/weertman/projects/sswd-evoepi

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export NUMBA_NUM_THREADS=1
export NUMBA_THREADING_LAYER=omp
export PYTHONUNBUFFERED=1

for W in $(seq 185 198); do
    CONFIG="experiments/calibration/W${W}_config.json"
    if [ ! -f "$CONFIG" ]; then
        echo "SKIP: $CONFIG not found"
        continue
    fi
    
    OUTDIR="results/calibration/W${W}"
    mkdir -p "$OUTDIR"
    
    for SEED in 42 123 999; do
        RESULT="$OUTDIR/result_seed${SEED}.json"
        
        # Skip if already complete
        if [ -f "$RESULT" ]; then
            echo "SKIP: W${W} seed ${SEED} already complete"
            continue
        fi
        
        echo "LAUNCH: W${W} seed ${SEED}"
        nohup python3 experiments/calibration_runner.py \
            --config "$CONFIG" \
            --seed "$SEED" \
            --output "$OUTDIR" \
            --K 5000 \
            --years 13 \
            > "$OUTDIR/log_seed${SEED}.txt" 2>&1 &
        
        # Small delay to avoid thundering herd
        sleep 0.5
    done
done

sleep 2
echo ""
echo "Launched. Active calibration processes:"
ps aux | grep "[c]alibration_runner" | wc -l
