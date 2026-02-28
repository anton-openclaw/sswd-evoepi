#!/bin/bash
# Launch wavefront calibration batch on Xeon
# Tests 4 D_P values to find the right wavefront speed
# All use K_half=200K (proven good for Alaska crash severity from R02)
# Wavefront enabled: disease seeds at CA-S only, spreads as wave
#
# Rounds:
#   W01: D_P=50, max_range=175   (conservative, short-range)
#   W02: D_P=100, max_range=350  (moderate)
#   W03: D_P=150, max_range=525  (long-range)
#   W04: D_P=200, max_range=700  (very long-range)
#
# Each with 3 seeds (42, 123, 999)
# 12 runs total, ~6 parallel batches of 2

set -e
cd ~/projects/sswd-evoepi

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1
export NUMBA_THREADING_LAYER=omp
export NUMBA_NUM_THREADS=8
export PYTHONUNBUFFERED=1

SEEDS="42 123 999"

# CA-S origin node indices (78 nodes, indices 306-377 + 631-636)
ORIGIN_NODES="306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,324,325,326,327,328,329,330,331,332,333,334,335,336,337,338,339,340,341,342,343,344,345,346,347,348,349,350,351,352,353,354,355,356,357,358,359,360,361,362,363,364,365,366,367,368,369,370,371,372,373,374,375,376,377,631,632,633,634,635,636"

# Common overrides (from R02 best + wavefront)
# K_half=200K nailed AK-PWS at 47% (target 50%)
COMMON_OVERRIDES="disease.K_half=200000,disease.wavefront_enabled=true,disease.activation_threshold=1.0"

declare -A ROUNDS
ROUNDS[W01]="spatial.D_P=50,spatial.D_P_max_range=175"
ROUNDS[W02]="spatial.D_P=100,spatial.D_P_max_range=350"
ROUNDS[W03]="spatial.D_P=150,spatial.D_P_max_range=525"
ROUNDS[W04]="spatial.D_P=200,spatial.D_P_max_range=700"

for ROUND in W01 W02 W03 W04; do
    OVERRIDES="${COMMON_OVERRIDES},${ROUNDS[$ROUND]}"
    OUTDIR="results/calibration/${ROUND}"
    mkdir -p "$OUTDIR"
    
    # Write config JSON
    python3 -c "
import json, sys
overrides = {}
for item in '${OVERRIDES}'.split(','):
    k, v = item.split('=')
    # Parse value types
    if v.lower() == 'true': v = True
    elif v.lower() == 'false': v = False
    else:
        try: v = int(v)
        except ValueError:
            try: v = float(v)
            except ValueError: pass
    overrides[k] = v
# Add origin nodes
overrides['disease.disease_origin_nodes'] = [${ORIGIN_NODES}]
json.dump(overrides, open('${OUTDIR}/config.json', 'w'), indent=2)
print(f'Config for ${ROUND}: {len(overrides)} overrides')
"
    
    for SEED in $SEEDS; do
        LOG="${OUTDIR}/seed_${SEED}.log"
        echo "Launching ${ROUND} seed=${SEED} → ${LOG}"
        
        nohup python3 experiments/calibration_runner.py \
            --config "${OUTDIR}/config.json" \
            --seed "$SEED" \
            --output "$OUTDIR" \
            --K 5000 \
            --years 13 \
            --disease-year 1 \
            > "$LOG" 2>&1 &
    done
done

echo ""
echo "Launched 12 wavefront calibration runs (4 rounds × 3 seeds)"
echo "Monitor: ps aux | grep calibration_runner | grep -v grep | wc -l"
echo "Check progress: tail -5 results/calibration/W*/seed_42.log"
