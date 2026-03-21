#!/bin/bash
# Generate formatted results for the tuning agent
# Run from project root

cd "$(dirname "$0")/../.."

echo "=== Formatting results for tuning agent ==="
ROUNDS=$(ls -d results/calibration/round_*/ 2>/dev/null | sort)
python3 experiments/calibration/format_results.py $ROUNDS
