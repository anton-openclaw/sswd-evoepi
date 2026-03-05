#!/bin/bash
# Launch W95-W104: Pathogen thermal adaptation sweep
# First runs with pathogen_adaptation=True
# Base: W90/W91 params (α=0.20, dynamic P_env, K½=800K)
set -e

cd ~/projects/sswd-evoepi
git pull

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1
export NUMBA_THREADING_LAYER=omp
export NUMBA_NUM_THREADS=7
export PYTHONUNBUFFERED=1

CONFIGS_DIR="experiments/calibration"
RESULTS_DIR="results/calibration"

for i in $(seq 95 104); do
    NAME="W${i}"
    CONFIG="${CONFIGS_DIR}/${NAME}_config.json"
    OUTDIR="${RESULTS_DIR}/${NAME}"
    LOG="${RESULTS_DIR}/${NAME}_log.txt"

    if [ ! -f "$CONFIG" ]; then
        echo "SKIP: $CONFIG not found"
        continue
    fi

    mkdir -p "$OUTDIR"

    echo "Launching $NAME..."
    nohup python3 experiments/calibration_runner.py \
        --config "$CONFIG" \
        --output "$OUTDIR" \
        --seeds 42 \
        --K 5000 \
        --years 13 \
        --disease-year 1 \
        > "$LOG" 2>&1 &

    echo "  PID: $! → $LOG"
    sleep 2
done

echo ""
echo "All launched. Monitor with: ps aux | grep calibration_runner | grep -v grep"
