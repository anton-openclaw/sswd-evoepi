#!/usr/bin/env bash
# Launch W249-W268: First clean sweep after audit bug fixes
# 20 configs × 1 seed = 20 runs, ~20GB each → 20 parallel safe on 503GB
#
# Run from project root: bash experiments/launch_W249_W268.sh

set -euo pipefail

RESULTS_DIR="results/calibration/W249-W268"
CONFIG_DIR="experiments/calibration"
RUNNER="experiments/calibration_runner.py"
MAX_PARALLEL=20

# Thread pinning (critical for Xeon)
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1
export NUMBA_NUM_THREADS=1
export NUMBA_THREADING_LAYER=omp
export PYTHONUNBUFFERED=1

mkdir -p "$RESULTS_DIR" logs/W249_W268

echo "=== W249-W268 Clean Sweep ==="
echo "Start time: $(date)"
echo "Max parallel: $MAX_PARALLEL"
echo ""

PIDS=()
NAMES=()

for i in $(seq 249 268); do
    NAME="W${i}"
    CONFIG="$CONFIG_DIR/${NAME}_config.json"
    OUTDIR="$RESULTS_DIR/$NAME"

    if [ ! -f "$CONFIG" ]; then
        echo "SKIP $NAME: config not found"
        continue
    fi

    # Skip if already has results
    if ls "$OUTDIR"/result_seed*.json 2>/dev/null | head -1 > /dev/null 2>&1; then
        echo "SKIP $NAME: results exist"
        continue
    fi

    mkdir -p "$OUTDIR"

    echo "LAUNCH $NAME"
    python3 "$RUNNER" \
        --config "$CONFIG" \
        --output "$OUTDIR" \
        > "logs/W249_W268/${NAME}.log" 2>&1 &

    PIDS+=($!)
    NAMES+=("$NAME")

    # Throttle
    while [ $(jobs -rp | wc -l) -ge $MAX_PARALLEL ]; do
        sleep 30
    done
done

echo ""
echo "All ${#PIDS[@]} jobs launched. Waiting..."
echo ""

# Wait and report
FAILED=0
for idx in "${!PIDS[@]}"; do
    pid=${PIDS[$idx]}
    name=${NAMES[$idx]}
    if wait "$pid"; then
        echo "DONE $name (exit 0)"
    else
        echo "FAIL $name (exit $?)"
        FAILED=$((FAILED + 1))
    fi
done

echo ""
echo "=== Complete ==="
echo "End time: $(date)"
echo "Failed: $FAILED / ${#PIDS[@]}"
