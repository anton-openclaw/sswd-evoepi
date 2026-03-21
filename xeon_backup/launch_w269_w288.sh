#!/bin/bash
# Launch W269-W288 calibration sweep (20 parallel runs)
# Xeon: 64 cores / 503GB RAM — 20 runs @ ~20GB each = ~400GB

set -e

cd ~/projects/sswd-evoepi

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export NUMBA_NUM_THREADS=1
export NUMBA_THREADING_LAYER=omp
export PYTHONUNBUFFERED=1

SWEEP_DIR="results/calibration/W269-W288"
CONFIG_DIR="experiments/calibration"

for W in $(seq 269 288); do
    WNAME="W${W}"
    OUT_DIR="${SWEEP_DIR}/${WNAME}"
    CONFIG="${CONFIG_DIR}/${WNAME}_config.json"
    LOG="${SWEEP_DIR}/${WNAME}.log"
    
    mkdir -p "${OUT_DIR}"
    
    echo "[$(date)] Launching ${WNAME}..."
    python3 experiments/calibration_runner.py --config "${CONFIG}" --output "${OUT_DIR}" > "${LOG}" 2>&1 &
done

echo ""
echo "[$(date)] All 20 runs launched."
echo "PIDs:"
ps aux | grep calibration_runner | grep -v grep | awk '{print $2, $NF}'
echo ""
echo "Monitor with: tail -f results/calibration/W269-W288/W*.log"

wait
echo "[$(date)] All runs complete."
