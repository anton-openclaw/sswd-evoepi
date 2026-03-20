#!/bin/bash
# K=1000 density-scaled calibration sweep
# Runs configs in parallel batches of 3 (each ~714MB, total ~2.1GB)
set -e

cd "$(dirname "$0")/.."
OUTDIR="results/k1000_scaled_sweep"
mkdir -p "$OUTDIR"

CONFIGS=(W321 W322 W323 W324 W325 W326 W327 W328)
PARALLEL=3
K=1000
K_REF=5000
YEARS=13
DISEASE_YEAR=1
SEED=42

run_config() {
    local cfg=$1
    local outpath="$OUTDIR/$cfg"
    
    # Skip if already done
    if [ -f "$outpath/result_seed${SEED}.json" ]; then
        echo "[$(date '+%H:%M:%S')] $cfg: SKIPPED (already done)"
        return 0
    fi
    
    mkdir -p "$outpath"
    echo "[$(date '+%H:%M:%S')] $cfg: STARTING"
    python3 -u experiments/calibration_runner.py \
        --config "experiments/calibration/${cfg}_config.json" \
        --seeds "$SEED" \
        --K "$K" \
        --K-ref "$K_REF" \
        --years "$YEARS" \
        --disease-year "$DISEASE_YEAR" \
        --output "$outpath" \
        > "$outpath/run.log" 2>&1
    
    if [ -f "$outpath/result_seed${SEED}.json" ]; then
        RMSLE=$(python3 -c "import json; d=json.load(open('$outpath/result_seed${SEED}.json')); print(f'{d[\"scoring\"][\"rmsle\"]:.3f}')")
        echo "[$(date '+%H:%M:%S')] $cfg: DONE (RMSLE=$RMSLE)"
    else
        echo "[$(date '+%H:%M:%S')] $cfg: FAILED"
    fi
}

echo "=== K=1000 Density-Scaled Calibration Sweep ==="
echo "Configs: ${CONFIGS[*]}"
echo "Parallel: $PARALLEL, K=$K, K_ref=$K_REF"
echo ""

# Run in batches
for ((i=0; i<${#CONFIGS[@]}; i+=PARALLEL)); do
    batch=("${CONFIGS[@]:i:PARALLEL}")
    echo "[$(date '+%H:%M:%S')] === Batch: ${batch[*]} ==="
    
    pids=()
    for cfg in "${batch[@]}"; do
        run_config "$cfg" &
        pids+=($!)
    done
    
    # Wait for batch
    for pid in "${pids[@]}"; do
        wait $pid
    done
    echo ""
done

echo "[$(date '+%H:%M:%S')] === ALL DONE ==="

# Summary
echo ""
echo "=== RESULTS SUMMARY ==="
python3 -c "
import json, os, glob

outdir = '$OUTDIR'
targets = {'AK-PWS':0.20,'AK-FN':0.20,'AK-FS':0.20,'BC-N':0.20,'SS-S':0.05,'JDF':0.02,'OR':0.0025,'CA-N':0.001}
scored = ['AK-PWS','AK-FN','AK-FS','BC-N','SS-S','JDF','OR','CA-N']

results = []
for cfg_dir in sorted(os.listdir(outdir)):
    rpath = os.path.join(outdir, cfg_dir, 'result_seed${SEED}.json')
    if os.path.exists(rpath):
        with open(rpath) as f:
            d = json.load(f)
        rmsle = d['scoring']['rmsle']
        rr = d['region_recovery']
        results.append((cfg_dir, rmsle, rr, d['wall_time_seconds']))

results.sort(key=lambda x: x[1])

header = f'{\"Config\":>8s} | {\"RMSLE\":>7s}'
for r in scored:
    header += f' | {r:>7s}'
header += f' | {\"Time\":>5s}'
print(header)
print('-' * len(header))

for name, rmsle, rr, wt in results:
    line = f'{name:>8s} | {rmsle:7.3f}'
    for r in scored:
        line += f' | {rr[r]*100:6.1f}%'
    line += f' | {wt/60:4.0f}m'
    print(line)

print()
print('Targets:          ', end='')
for r in scored:
    print(f' | {targets[r]*100:6.2f}%', end='')
print()
"
