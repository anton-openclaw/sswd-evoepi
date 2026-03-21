#!/usr/bin/env bash
# Launch W329-W340 sweep: K=1000, density-scaled, 8 parallel
# Based on W325 (RMSLE=0.348, best so far)
# Targeting: BC-N undershoot (4.3% vs 20%), AK-PWS overshoot (30% vs 20%)
set -euo pipefail

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export NUMBA_NUM_THREADS=1
export NUMBA_THREADING_LAYER=omp

REPO="/home/starbot/.openclaw/workspace/sswd-evoepi"
RUNNER="$REPO/experiments/calibration_runner.py"
CONFIGS="$REPO/experiments/calibration"
OUTBASE="$REPO/results/k1000_sweep_w329"
LOGDIR="$OUTBASE/logs"

mkdir -p "$OUTBASE" "$LOGDIR"

MAX_PARALLEL=8
RUNNING=0

for w in W329 W330 W331 W332 W333 W334 W335 W336 W337 W338 W339 W340; do
    OUTDIR="$OUTBASE/$w"
    CONFIG="$CONFIGS/${w}_config.json"
    LOG="$LOGDIR/${w}.log"
    
    if [ -f "$OUTDIR/result_seed42.json" ]; then
        echo "SKIP $w (already complete)"
        continue
    fi
    
    mkdir -p "$OUTDIR"
    echo "$(date '+%H:%M:%S') LAUNCH $w"
    
    python3 "$RUNNER" \
        --config "$CONFIG" \
        --seed 42 \
        --K 1000 \
        --K-ref 5000 \
        --years 13 \
        --disease-year 1 \
        --output "$OUTDIR" \
        > "$LOG" 2>&1 &
    
    RUNNING=$((RUNNING + 1))
    
    if [ "$RUNNING" -ge "$MAX_PARALLEL" ]; then
        wait -n
        RUNNING=$((RUNNING - 1))
    fi
done

wait
echo "$(date '+%H:%M:%S') ALL DONE"

# Collect results
echo ""
echo "=== RESULTS ==="
for w in W329 W330 W331 W332 W333 W334 W335 W336 W337 W338 W339 W340; do
    RESULT="$OUTBASE/$w/result_seed42.json"
    if [ -f "$RESULT" ]; then
        RMSLE=$(python3 -c "import json; r=json.load(open('$RESULT')); print(f'{r[\"scoring\"][\"rmsle\"]:.4f}')" 2>/dev/null || echo "ERR")
        AK_PWS=$(python3 -c "import json; r=json.load(open('$RESULT')); print(f'{r[\"region_recovery\"][\"AK-PWS\"]*100:.1f}%')" 2>/dev/null || echo "?")
        BC_N=$(python3 -c "import json; r=json.load(open('$RESULT')); print(f'{r[\"region_recovery\"][\"BC-N\"]*100:.1f}%')" 2>/dev/null || echo "?")
        OR=$(python3 -c "import json; r=json.load(open('$RESULT')); print(f'{r[\"region_recovery\"][\"OR\"]*100:.1f}%')" 2>/dev/null || echo "?")
        CA_N=$(python3 -c "import json; r=json.load(open('$RESULT')); print(f'{r[\"region_recovery\"][\"CA-N\"]*100:.1f}%')" 2>/dev/null || echo "?")
        JDF=$(python3 -c "import json; r=json.load(open('$RESULT')); print(f'{r[\"region_recovery\"][\"JDF\"]*100:.1f}%')" 2>/dev/null || echo "?")
        SS_S=$(python3 -c "import json; r=json.load(open('$RESULT')); print(f'{r[\"region_recovery\"][\"SS-S\"]*100:.1f}%')" 2>/dev/null || echo "?")
        echo "$w: RMSLE=$RMSLE  AK-PWS=$AK_PWS  BC-N=$BC_N  SS-S=$SS_S  JDF=$JDF  OR=$OR  CA-N=$CA_N"
    else
        echo "$w: NO RESULT"
    fi
done
echo ""
echo "TARGETS: AK-PWS=20%  AK-FN=20%  AK-FS=20%  BC-N=20%  SS-S=5%  JDF=2%  OR=0.25%  CA-N=0.1%"
