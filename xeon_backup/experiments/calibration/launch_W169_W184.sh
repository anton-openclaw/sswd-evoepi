#!/bin/bash
# Launch W169-W184: Post-dispersal-fix recalibration sweep
# Primary: s0 sweep (0.01-1.0), 8 configs
# Secondary: s0=0.10 × alpha_self combos, 4 configs
# Tertiary: s0=0.10 × n_connectivity, 2 configs
# Quaternary: s0=0.10 × alpha_env, 2 configs
# Total: 16 configs × 3 seeds = 48 runs
set -e
cd ~/projects/sswd-evoepi

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1
export NUMBA_NUM_THREADS=1
export NUMBA_THREADING_LAYER=omp
export PYTHONUNBUFFERED=1

RESULTS_BASE="results/calibration"
mkdir -p "${RESULTS_BASE}/W169-W184"
MAX_PARALLEL=24
RUNNING=0

for W in $(seq 169 184); do
    CONFIG="experiments/calibration/W${W}_config.json"
    if [ ! -f "$CONFIG" ]; then
        echo "SKIP: $CONFIG not found"
        continue
    fi
    
    for SEED in 42 123 999; do
        OUTDIR="${RESULTS_BASE}/W${W}"
        LOGFILE="${RESULTS_BASE}/W169-W184/W${W}_seed${SEED}.log"
        
        echo "Launching W${W} seed=${SEED}"
        
        nohup python3 experiments/calibration_runner.py \
            --config "$CONFIG" \
            --output "$OUTDIR" \
            --seed "$SEED" \
            > "$LOGFILE" 2>&1 &
        
        RUNNING=$((RUNNING + 1))
        
        if [ $RUNNING -ge $MAX_PARALLEL ]; then
            echo "Hit $MAX_PARALLEL parallel, waiting for one to finish..."
            wait -n
            RUNNING=$((RUNNING - 1))
        fi
    done
done

echo "All 48 runs launched (W169-W184 × 3 seeds)."
echo "Monitor: tail -f ${RESULTS_BASE}/W169-W184/W169_seed42.log"
wait
echo "All runs complete."
