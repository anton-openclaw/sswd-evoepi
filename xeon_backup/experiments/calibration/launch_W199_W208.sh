#!/bin/bash
cd /home/weertman/projects/sswd-evoepi

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export NUMBA_NUM_THREADS=1
export NUMBA_THREADING_LAYER=omp
export PYTHONUNBUFFERED=1

for W in $(seq 199 208); do
    CONFIG="experiments/calibration/W${W}_config.json"
    [ ! -f "$CONFIG" ] && continue
    
    OUTDIR="results/calibration/W${W}"
    mkdir -p "$OUTDIR"
    
    for SEED in 42 123 999; do
        RESULT="$OUTDIR/result_seed${SEED}.json"
        [ -f "$RESULT" ] && echo "SKIP: W${W} seed ${SEED}" && continue
        
        echo "LAUNCH: W${W} seed ${SEED}"
        nohup python3 experiments/calibration_runner.py \
            --config "$CONFIG" \
            --seed "$SEED" \
            --output "$OUTDIR" \
            --K 5000 \
            --years 13 \
            > "$OUTDIR/log_seed${SEED}.txt" 2>&1 &
        sleep 0.5
    done
done

sleep 3
echo ""
echo "Active calibration processes: $(ps aux | grep '[c]alibration_runner' | wc -l)"
