#!/bin/bash
# Launch W135-W144 calibration runs on Xeon
# Respects 24 parallel process limit, accounting for any still-running W125-W134 jobs

set -e
cd ~/projects/sswd-evoepi

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1
export NUMBA_THREADING_LAYER=omp
export PYTHONUNBUFFERED=1

MAX_PARALLEL=24
SEEDS="42 123 999"
LOGDIR="results/calibration/W135-W144"
mkdir -p "$LOGDIR"

# Count currently running calibration processes
running=$(ps aux | grep calibration_runner | grep -v grep | wc -l)
echo "Currently running calibration processes: $running"
available=$((MAX_PARALLEL - running))
echo "Available slots: $available"

if [ $available -le 0 ]; then
    echo "ERROR: No available slots ($running/$MAX_PARALLEL running). Wait for W125-W134 to finish."
    exit 1
fi

launched=0
for w in 135 136 137 138 139 140 141 142 143 144; do
    for seed in $SEEDS; do
        if [ $launched -ge $available ]; then
            echo "Hit $available parallel limit (from $MAX_PARALLEL - $running running), pausing..."
            echo "Waiting for a slot to free up..."
            while true; do
                sleep 60
                cur=$(ps aux | grep calibration_runner | grep -v grep | wc -l)
                if [ $cur -lt $MAX_PARALLEL ]; then
                    running=$cur
                    available=$((MAX_PARALLEL - running))
                    echo "  Slots freed: $available available (${cur} running)"
                    break
                fi
            done
        fi
        
        outdir="results/calibration/W${w}"
        mkdir -p "$outdir"
        
        echo "Launching W${w} seed=${seed} → ${outdir}"
        setsid nohup python3 experiments/calibration_runner.py \
            --config experiments/calibration/W${w}_config.json \
            --output "$outdir" \
            --seed $seed \
            > "${LOGDIR}/W${w}_seed${seed}.log" 2>&1 &
        
        launched=$((launched + 1))
        sleep 1  # stagger launches slightly
    done
done

echo ""
echo "Launched $launched processes total."
echo "Logs: $LOGDIR/"
echo "Monitor: ps aux | grep calibration_runner | grep -v grep | wc -l"
