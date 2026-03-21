#!/bin/bash
# W289-W308 calibration sweep
# 20 configs × 1 seed (137), K=5000, 13 years
# Run 20 parallel (each ~20GB RAM, 20×20=400GB / 503GB)

set -e
cd ~/projects/sswd-evoepi

RESULTS_DIR="results/calibration/W289-W308"
mkdir -p "$RESULTS_DIR"

# Thread control
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export NUMBA_NUM_THREADS=1
export NUMBA_THREADING_LAYER=omp

launch() {
    local name=$1
    local config="experiments/calibration/${name}_config.json"
    local outdir="$RESULTS_DIR/$name"
    
    if [ -f "$outdir/result_seed137.json" ]; then
        echo "SKIP $name (already done)"
        return
    fi
    
    mkdir -p "$outdir"
    echo "LAUNCH $name → $outdir"
    
    nohup python3 experiments/calibration_runner.py \
        --config "$config" \
        --seed 137 \
        --K 5000 \
        --output "$outdir" \
        > "$outdir/run.log" 2>&1 &
    
    echo "  PID=$!"
}

echo "=== W289-W308 Sweep ==="
echo "Start: $(date)"
echo ""

for i in $(seq 289 308); do
    launch "W$i"
    sleep 2  # stagger starts
done

echo ""
echo "All 20 launched. Monitor: ps aux | grep calibration_runner | grep -v grep | wc -l"
echo "Check results: for d in $RESULTS_DIR/W*/; do echo \$(basename \$d): \$(ls \$d/result_seed137.json 2>/dev/null && echo DONE || echo running); done"
