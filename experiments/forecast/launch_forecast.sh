#!/bin/bash
# Launch forward projection experiments (2012-2050)
# All use W330 best-calibrated baseline, K=1000 with density scaling
# 38 years total, disease at year 1 (2013)

set -euo pipefail

REPO="/home/starbot/.openclaw/workspace/sswd-evoepi"
RUNNER="$REPO/experiments/calibration_runner.py"
CONFIG_DIR="$REPO/experiments/forecast"
OUTPUT_BASE="$REPO/results/forecast_2050"
YEARS=38
K=1000
K_REF=5000
SEED=42
DISEASE_YEAR=1

# Thread control
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export NUMBA_NUM_THREADS=1
export NUMBA_THREADING_LAYER=omp

mkdir -p "$OUTPUT_BASE/logs"

SCENARIOS=(
    "F01_baseline"
    "F02_disease_removal_2018"
    "F03_no_pathogen_evolution"
    "F04_double_virulence"
    "F05_no_gradient_mechanisms"
    "F06_disease_removal_2025"
)

echo "$(date '+%H:%M:%S') Starting forecast experiments (${#SCENARIOS[@]} scenarios, $YEARS years each)"
echo "  K=$K, K_ref=$K_REF, seed=$SEED, disease_year=$DISEASE_YEAR"
echo ""

# Launch all scenarios in parallel (6 × ~700MB = ~4.2GB, well within 32GB)
for scenario in "${SCENARIOS[@]}"; do
    config="$CONFIG_DIR/${scenario}.json"
    outdir="$OUTPUT_BASE/$scenario"
    logfile="$OUTPUT_BASE/logs/${scenario}.log"
    
    if [ ! -f "$config" ]; then
        echo "$(date '+%H:%M:%S') SKIP $scenario — config not found"
        continue
    fi
    
    if [ -f "$outdir/result_seed${SEED}.json" ]; then
        echo "$(date '+%H:%M:%S') SKIP $scenario — already complete"
        continue
    fi
    
    mkdir -p "$outdir"
    
    echo "$(date '+%H:%M:%S') LAUNCH $scenario"
    setsid python3 "$RUNNER" \
        --config "$config" \
        --seed "$SEED" \
        --K "$K" \
        --K-ref "$K_REF" \
        --years "$YEARS" \
        --disease-year "$DISEASE_YEAR" \
        --output "$outdir" \
        > "$logfile" 2>&1 &
    
    echo "$(date '+%H:%M:%S') LAUNCH $scenario" >> "$OUTPUT_BASE/master.log"
done

echo ""
echo "$(date '+%H:%M:%S') All ${#SCENARIOS[@]} scenarios launched"
echo "  Output: $OUTPUT_BASE/"
echo "  Logs: $OUTPUT_BASE/logs/"
echo "  Expected runtime: ~2.5h per scenario (38 years at K=1000)"
echo ""
echo "Monitor with: ps aux | grep calibration_runner | grep -v grep"
