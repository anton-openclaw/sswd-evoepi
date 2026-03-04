#!/bin/bash
# Launch W65-W74 calibration runs (dynamic P_env exploration)
# Run on Xeon with 10 parallel processes (1 seed each for speed)

set -e

REPO=~/projects/sswd-evoepi
RUNNER=$REPO/experiments/calibration_runner.py
CONFIGS=$REPO/experiments/calibration
RESULTS=$REPO/results/calibration
NETWORK=$REPO/results/overwater/distance_matrix_900.npz
SST_DIR=$REPO/data/sst/site_sst
SITES=$REPO/data/sites_900.csv

# Threading limits — critical for parallel runs
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1
export NUMBA_THREADING_LAYER=omp
export NUMBA_NUM_THREADS=7
export PYTHONUNBUFFERED=1

echo "Launching W65-W74 (dynamic P_env) at $(date)"

for W in W65 W66 W67 W68 W69 W70 W71 W72 W73 W74; do
    CFG=$CONFIGS/${W}_config.json
    OUTDIR=$RESULTS/$W
    
    if [ ! -f "$CFG" ]; then
        echo "SKIP $W: config not found"
        continue
    fi
    
    if [ -f "$OUTDIR/combined_results.json" ]; then
        echo "SKIP $W: already complete"
        continue
    fi
    
    mkdir -p "$OUTDIR"
    
    echo "Starting $W..."
    nohup python3 "$RUNNER" \
        --config "$CFG" \
        --output "$OUTDIR" \
        --network "$NETWORK" \
        --sst-dir "$SST_DIR" \
        --sites "$SITES" \
        --seeds 42 \
        --K 5000 \
        --years 13 \
        --disease-year 1 \
        --monthly-snapshots \
        > "$OUTDIR/run.log" 2>&1 &
    
    echo "  PID=$! → $OUTDIR"
done

echo ""
echo "All launched. Monitor with: ps aux | grep calibration_runner | grep -v grep"
echo "Check logs: tail -f $RESULTS/W65/run.log"
