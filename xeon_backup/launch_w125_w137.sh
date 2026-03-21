#!/bin/bash
cd ~/projects/sswd-evoepi
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export NUMBA_THREADING_LAYER=omp
export PYTHONUNBUFFERED=1

echo "[$(date)] Launching W125-W137 calibration round"
for w in W125 W126 W127 W128 W129 W130 W131 W132 W133 W134 W135 W136 W137; do
    cfg="experiments/calibration/${w}_config.json"
    [ ! -f "$cfg" ] && continue
    for seed in 42 123 999; do
        outdir="results/calibration/${w}"
        mkdir -p "$outdir"
        nohup python3 experiments/calibration_runner.py \
            --config "$cfg" --seeds "$seed" --output "$outdir/" \
            > "${outdir}/stdout_seed${seed}.log" 2>&1 &
    done
done
echo "[$(date)] All runs launched"
