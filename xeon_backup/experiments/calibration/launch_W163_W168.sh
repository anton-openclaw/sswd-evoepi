#!/bin/bash
# Rolling Launcher for W163-W168: Temperature-asymmetric disease severity sweep
# Wave 1: W163-W166 (12 runs), Wave 2: W167-W168 (6 runs) after disease crash
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
LOG_DIR="${RESULTS_BASE}/W163-W168"
mkdir -p "$LOG_DIR"

SEEDS=(42 123 999)
WAVE1=(163 164 165 166)
WAVE2=(167 168)

CRASH_YEAR=3
CRASH_POP_FRAC=0.20
MEM_FREE_MIN_GB=60
POLL_INTERVAL=120
MAX_WAIT=14400

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${LOG_DIR}/launcher.log"
}

launch_run() {
    local W=$1 SEED=$2
    local CONFIG="experiments/calibration/W${W}_config.json"
    local OUTDIR="${RESULTS_BASE}/W${W}"
    local LOGFILE="${LOG_DIR}/W${W}_seed${SEED}.log"
    if [ ! -f "$CONFIG" ]; then log "SKIP: $CONFIG not found"; return 1; fi
    mkdir -p "$OUTDIR"
    log "Launching W${W} seed=${SEED}"
    nohup python3 experiments/calibration_runner.py \
        --config "$CONFIG" --output "$OUTDIR" --seed "$SEED" \
        > "$LOGFILE" 2>&1 &
}

get_free_mem_gb() { free -g | awk '/^Mem:/ {print $7}'; }

check_wave1_crashed() {
    local total=$((${#WAVE1[@]} * ${#SEEDS[@]})) crashed=0
    for W in "${WAVE1[@]}"; do
        for SEED in "${SEEDS[@]}"; do
            local r="${RESULTS_BASE}/W${W}/result_seed${SEED}.json"
            local f="${RESULTS_BASE}/W${W}/checkpoint_seed${SEED}.json"
            if [ -f "$r" ]; then crashed=$((crashed+1)); continue; fi
            if [ ! -f "$f" ]; then continue; fi
            local ok=$(python3 -c "
import json; d=json.load(open('$f'))
print('yes' if d['year']>=$CRASH_YEAR or d.get('pop_frac',1.0)<$CRASH_POP_FRAC else 'no')" 2>/dev/null)
            [ "$ok" = "yes" ] && crashed=$((crashed+1))
        done
    done
    log "Crash check: ${crashed}/${total} wave 1 runs past disease crash"
    [ $crashed -eq $total ] && return 0 || return 1
}

log "========================================="
log "W163-W168: Temperature-asymmetric disease severity sweep"
log "Wave 1: W163-W166 (α_env, δ_env sweep)"
log "Wave 2: W167-W168 (v_max_warm sweep)"
log "========================================="

log "Launching Wave 1..."
for W in "${WAVE1[@]}"; do
    for SEED in "${SEEDS[@]}"; do launch_run "$W" "$SEED"; sleep 2; done
done
log "Wave 1: 12 runs launched"

log "Waiting for wave 1 disease crash..."
WAITED=0
while true; do
    sleep $POLL_INTERVAL; WAITED=$((WAITED+POLL_INTERVAL))
    if check_wave1_crashed; then
        FREE=$(get_free_mem_gb)
        log "Wave 1 past crash! Free: ${FREE}GB"
        if [ "$FREE" -ge "$MEM_FREE_MIN_GB" ]; then log "Memory OK. Launching wave 2."; break; fi
        log "Memory tight (${FREE}<${MEM_FREE_MIN_GB}GB). Waiting..."
    fi
    [ $WAITED -ge $MAX_WAIT ] && { log "Max wait reached. Launching wave 2 anyway."; break; }
    [ $((WAITED % 600)) -eq 0 ] && log "Still waiting... (${WAITED}s, $(get_free_mem_gb)GB free)"
done

log "Launching Wave 2..."
for W in "${WAVE2[@]}"; do
    for SEED in "${SEEDS[@]}"; do launch_run "$W" "$SEED"; sleep 2; done
done
log "Wave 2: 6 runs launched. All 18 active."
log "Monitor: ls ${RESULTS_BASE}/W16{3,4,5,6,7,8}/result_seed*.json 2>/dev/null | wc -l"
