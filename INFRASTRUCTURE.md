# Infrastructure & Compute Environment

## Current Setup (April 2026)

All experiments run locally on a single workstation:

| Component | Spec |
|-----------|------|
| CPU | AMD Ryzen 7 5700 (8 cores / 16 threads) |
| RAM | 32 GB DDR4 |
| GPU | NVIDIA RTX 5060 (8 GB VRAM) |
| OS | Ubuntu GNOME |

### Previous Environment

A Xeon W-3365 workstation (64 cores, 503 GB RAM, WSL2) was used during
early development. It became permanently unavailable in March 2026. All
work since then targets the Ryzen box above.

## Memory & K Sensitivity

The simulation's particle count `K` (carrying-capacity proxy) controls both
memory usage and statistical fidelity. Benchmarks on the Ryzen workstation:

| K | Peak RSS | Notes |
|------|----------|-------|
| 500 | ~499 MB | Comfortable parallel |
| 1000 | ~714 MB | Comfortable parallel |
| 2000 | fits | Sequential only (OOM if parallelised) |
| 5000 | — | OOM-killed |

### Fidelity vs K (RMSLE against K=10000 reference)

| K | RMSLE |
|------|-------|
| 500 | 1.558 |
| 1000 | 1.441 |
| 2000 | 1.193 |

Higher K yields meaningfully better density estimates. K=2000 is the
practical ceiling on this hardware.

### Recommended Approach

- **Publication experiments:** K=2000, run sequentially (one replicate at a
  time). Maximises fidelity within the 32 GB envelope.
- **Diagnostics / iteration:** K=1000, run in parallel. Fast turnaround,
  acceptable fidelity.

## Thread Pinning (Parallel Runs)

When running multiple replicates in parallel, pin each worker to avoid
Numba/OpenBLAS thread over-subscription:

```bash
export NUMBA_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
```

Without these, parallel workers fight over cores and wall-clock time
actually increases.

## Agent Harness

The AI research-assistant workflow migrated from **OpenClaw** to **Hermes**
in March 2026. All automation scripts and cron-driven experiment pipelines
now use the Hermes agent harness.
