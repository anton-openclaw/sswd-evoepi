#!/bin/bash
# Run 3 diagnostic configs locally at K=2000
set -e
cd /home/starbot/.openclaw/workspace/sswd-evoepi

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export NUMBA_NUM_THREADS=1
export NUMBA_THREADING_LAYER=omp

OUTBASE="results/local_diagnostic"

# W285 is already running, wait for it
echo "Waiting for W285 baseline..."
while [ ! -f "$OUTBASE/W285_K2000/result_seed137.json" ]; do
    sleep 60
done
echo "W285 done!"

# W310: δ_env(T) only
echo "Starting W310 (Phase 1: δ_env(T) only)..."
python3 experiments/calibration_runner.py \
    --config experiments/calibration/W310_config.json \
    --seeds 137 --K 2000 --years 13 --disease-year 1 \
    --output "$OUTBASE/W310_K2000" \
    2>&1 | tee "$OUTBASE/W310_K2000.log"
echo "W310 done!"

# W316: δ_env(T) + P_env gating
echo "Starting W316 (Phase 1+2: δ_env(T) + P_env gating)..."
python3 experiments/calibration_runner.py \
    --config experiments/calibration/W316_config.json \
    --seeds 137 --K 2000 --years 13 --disease-year 1 \
    --output "$OUTBASE/W316_K2000" \
    2>&1 | tee "$OUTBASE/W316_K2000.log"
echo "W316 done!"

echo "=== ALL 3 DIAGNOSTIC RUNS COMPLETE ==="
