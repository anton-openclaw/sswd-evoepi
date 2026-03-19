#!/bin/bash
# K sensitivity / CV sweep — W285 baseline at different K values
# Tests how much stochastic noise changes with K
# Sequential runs (one at a time to avoid OOM)
set -e
cd /home/starbot/.openclaw/workspace/sswd-evoepi

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export NUMBA_NUM_THREADS=1
export NUMBA_THREADING_LAYER=omp

CONFIG="experiments/calibration/W285_config.json"
OUTBASE="results/k_cv_sweep"
mkdir -p "$OUTBASE"

# K values to test (ascending — smallest/fastest first)
# K=500: ~2GB, ~30min  |  K=1000: ~4GB, ~1h  |  K=2000: ~8GB, ~2-3h
# K=5000 likely OOMs locally (20GB)
for K in 500 1000 2000; do
    for SEED in 42 137 256; do
        OUTDIR="${OUTBASE}/K${K}_seed${SEED}"
        RESULT="${OUTDIR}/result_seed${SEED}.json"
        
        # Skip if already done
        if [ -f "$RESULT" ]; then
            echo "SKIP: K=${K} seed=${SEED} (already done)"
            continue
        fi
        
        mkdir -p "$OUTDIR"
        echo "$(date): Starting K=${K} seed=${SEED}..."
        python3 experiments/calibration_runner.py \
            --config "$CONFIG" \
            --seeds "$SEED" \
            --K "$K" \
            --years 13 \
            --disease-year 1 \
            --output "$OUTDIR" \
            2>&1 | tee "${OUTDIR}/run.log"
        echo "$(date): Done K=${K} seed=${SEED}"
    done
done

echo ""
echo "=== K CV SWEEP COMPLETE ==="
echo ""

# Extract results
python3 -c "
import json, os, numpy as np

base = '$OUTBASE'
data = {}
for K in [500, 1000, 2000]:
    rmsles = []
    ak_pws = []
    or_vals = []
    ca_n = []
    for seed in [42, 137, 256]:
        path = os.path.join(base, f'K{K}_seed{seed}', f'result_seed{seed}.json')
        if not os.path.isfile(path):
            continue
        with open(path) as f:
            d = json.load(f)
        s = d.get('scoring', {})
        rr = d.get('region_recovery', {})
        rmsle = s.get('rmsle', float('inf'))
        if rmsle == float('inf'):
            continue
        rmsles.append(rmsle)
        ak_pws.append(rr.get('AK-PWS', 0))
        or_vals.append(rr.get('OR', 0))
        ca_n.append(rr.get('CA-N', 0))
    
    if len(rmsles) >= 2:
        print(f'K={K}: {len(rmsles)} seeds')
        print(f'  RMSLE: mean={np.mean(rmsles):.4f} std={np.std(rmsles):.4f} CV={np.std(rmsles)/np.mean(rmsles)*100:.1f}%')
        print(f'  AK-PWS: mean={np.mean(ak_pws):.1%} std={np.std(ak_pws):.1%}')
        print(f'  OR:     mean={np.mean(or_vals):.2%} std={np.std(or_vals):.2%}')
        print(f'  CA-N:   mean={np.mean(ca_n):.2%} std={np.std(ca_n):.2%}')
        print()
    else:
        print(f'K={K}: only {len(rmsles)} valid seeds (need >=2)')
"
