#!/bin/bash
# Launch forecast runs F01-F04 (W154 → 2050)
# Waits for W167-W168 to finish before launching
# 12 total runs: 4 configs × 3 seeds, ~20GB each = ~240GB
set -euo pipefail

cd ~/projects/sswd-evoepi

CONFIGS_DIR="experiments/calibration"
RESULTS_DIR="results/calibration"
SEEDS="42 123 999"
YEARS=38

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1
export NUMBA_NUM_THREADS=1
export NUMBA_THREADING_LAYER=omp
export PYTHONUNBUFFERED=1

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

# --- Phase 0: Wait for W167-W168 to finish ---
log "Checking for running W167-W168 processes..."
while true; do
    running=$(ps aux | grep "[c]alibration_runner.*W16[78]" | wc -l)
    if [ "$running" -eq 0 ]; then
        log "W167-W168 complete. Proceeding to forecast launch."
        break
    fi
    log "W167-W168 still running ($running processes). Waiting 5 minutes..."
    sleep 300
done

# Short cool-down for memory reclaim
sleep 30

# --- Phase 1: Launch all forecast runs ---
log "=== Launching F01-F04 forecast batch (38yr to 2050) ==="

PIDS=()
for run_id in F01 F02 F03 F04; do
    config="${CONFIGS_DIR}/${run_id}_config.json"
    outdir="${RESULTS_DIR}/${run_id}"
    
    if [ ! -f "$config" ]; then
        log "ERROR: Config not found: $config"
        continue
    fi
    
    mkdir -p "$outdir"
    
    for seed in $SEEDS; do
        # Skip if result already exists
        if [ -f "${outdir}/result_seed${seed}.json" ]; then
            log "SKIP: ${run_id} seed=${seed} (result exists)"
            continue
        fi
        
        log "Launching ${run_id} seed=${seed} (${YEARS}yr)..."
        setsid python3 experiments/calibration_runner.py \
            --config "$config" \
            --output "$outdir" \
            --seed "$seed" \
            --K 5000 \
            --years "$YEARS" \
            --disease-year 1 \
            > "${outdir}/run_seed${seed}.log" 2>&1 &
        
        PIDS+=($!)
        log "  PID=$! started"
        sleep 2  # stagger startup
    done
done

log "=== All ${#PIDS[@]} forecast runs launched ==="
log "PIDs: ${PIDS[*]}"
log "Monitor: ps aux | grep calibration_runner"
log "Estimated completion: 20-35 hours from now"

# Don't wait — let them run in background
