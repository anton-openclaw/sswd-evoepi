#!/bin/bash
# Launch W75-W84 calibration batch on Xeon
# K_half / P_env_max sweep to lift recovery levels (W71 base + dynamic P_env)
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

for i in $(seq 75 84); do
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
        --sst-start-year 2012 \
        > "$LOG" 2>&1 &

    echo "  PID: $! → $LOG"
    sleep 2
done

echo ""
echo "All launched. Monitor with: ps aux | grep calibration_runner | grep -v grep"
