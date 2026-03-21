#!/bin/bash
# Launch W229-W248 calibration sweep (20 configs × seed 137)
# Xeon: 64 cores, 503GB RAM — 20 parallel runs is comfortable

set -e
cd ~/projects/sswd-evoepi

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export NUMBA_NUM_THREADS=1
export NUMBA_THREADING_LAYER=omp
export PYTHONUNBUFFERED=1

CONFIG_DIR="experiments/configs/W229-W248"
RESULTS_BASE="results/calibration"
SEED=137
RUNNING=0
MAX_PARALLEL=20

for config_file in $CONFIG_DIR/W*.json; do
    name=$(basename "$config_file" .json)
    result_dir="$RESULTS_BASE/$name"
    
    # Skip if already completed
    if [ -f "$result_dir/result_seed${SEED}.json" ]; then
        echo "SKIP $name — already complete"
        continue
    fi
    
    mkdir -p "$result_dir"
    
    echo "LAUNCH $name (seed $SEED)"
    nohup python3 experiments/calibration_runner.py \
        --config "$config_file" \
        --seeds $SEED \
        --output-dir "$result_dir" \
        > "$result_dir/log_seed${SEED}.txt" 2>&1 &
    
    RUNNING=$((RUNNING + 1))
    
    if [ $RUNNING -ge $MAX_PARALLEL ]; then
        echo "Reached $MAX_PARALLEL parallel — waiting..."
        wait -n
        RUNNING=$((RUNNING - 1))
    fi
    
    # Small delay to stagger startup
    sleep 2
done

echo ""
echo "All $RUNNING runs launched. PIDs:"
jobs -l
echo ""
echo "Monitor: ls results/calibration/W2{29..48}/result_seed137.json 2>/dev/null | wc -l"
