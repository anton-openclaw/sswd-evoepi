#!/bin/bash
# Launch W209-W228: Asymmetric lever sweep (20 configs × 1 seed = 20 runs)
# Base: W201 (s0=1.0, r_total=0.003), seed 137 only
set -e
cd ~/projects/sswd-evoepi

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export NUMBA_NUM_THREADS=1
export NUMBA_THREADING_LAYER=omp
export PYTHONUNBUFFERED=1

for w in $(seq 209 228); do
    CONFIG="experiments/calibration/W${w}_config.json"
    OUTDIR="results/calibration/W${w}"
    mkdir -p "$OUTDIR"

    # Skip if already has result
    if [ -f "$OUTDIR/result_seed137.json" ]; then
        echo "W${w}: already complete, skipping"
        continue
    fi

    echo "Launching W${w}..."
    nohup python3 experiments/calibration_runner.py \
        --config "$CONFIG" \
        --seed 137 \
        --output "$OUTDIR" \
        --K 5000 \
        --years 13 \
        > "$OUTDIR/log_seed137.txt" 2>&1 &

    sleep 0.5
done

echo ""
echo "All launched. $(ps aux | grep calibration_runner | grep -v grep | wc -l) processes running."
