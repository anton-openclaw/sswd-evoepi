#!/bin/bash
set -e
cd ~/projects/sswd-evoepi
git pull origin main

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1
export NUMBA_THREADING_LAYER=omp NUMBA_NUM_THREADS=7 PYTHONUNBUFFERED=1

for ROUND in W53 W54 W55 W56 W57 W58 W59 W60 W61 W62 W63 W64; do
    OUTDIR="results/calibration/${ROUND}"
    mkdir -p "$OUTDIR"
    echo "Launching $ROUND..."
    nohup python3 experiments/calibration_runner.py \
        --config "experiments/calibration/${ROUND}_config.json" \
        --seed 42 --output "$OUTDIR" --K 5000 --years 13 --disease-year 1 \
        > "results/calibration/${ROUND}_log.txt" 2>&1 &
    echo "  PID: $!"
    sleep 2
done
echo "All 12 launched."
ps aux | grep calibration_runner | grep -v grep | wc -l
echo "processes running"
