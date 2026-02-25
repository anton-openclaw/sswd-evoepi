#!/bin/bash
# Monitor reintroduction demo progress
# Called by cron to check if the demo finished and run analysis

cd /home/starbot/.openclaw/workspace/sswd-evoepi
LOGFILE=results/reintroduction/demo_log.txt
RESULTS=results/reintroduction/results_demo.json
FIGDIR=results/reintroduction/figures

# Check if demo is still running
PIDS=$(pgrep -f "reintroduction_monterey.py --mode demo" || true)
if [ -n "$PIDS" ]; then
    # Still running — check progress
    STARTED=$(grep -c 'Starting\.\.\.' "$LOGFILE" 2>/dev/null || echo 0)
    DONE=$(grep -c 'Done in' "$LOGFILE" 2>/dev/null || echo 0)
    MEM=$(free -h | awk '/^Mem:/{print $3"/"$2}')
    echo "RUNNING: $DONE/$STARTED scenarios done. Memory: $MEM. PIDs: $PIDS"
    exit 0
fi

# Demo process gone — check if results exist
if [ ! -f "$RESULTS" ]; then
    echo "ERROR: Demo process ended but no results file at $RESULTS"
    # Check for crash in log
    tail -20 "$LOGFILE"
    exit 1
fi

echo "Demo complete! Results at $RESULTS"
echo "Running analysis..."

# Run analysis
python3 experiments/analyze_reintroduction.py \
    --results results/reintroduction/results_demo.json \
    --outdir results/reintroduction/figures \
    2>&1

echo "Analysis done. Figures in $FIGDIR:"
ls -la "$FIGDIR/"

# Git commit
git add results/reintroduction/
git add experiments/analyze_reintroduction.py
git commit -m "Reintroduction demo results + 7 analysis figures

- 907-node network, K=5000, 19 scenarios × 3 reps = 57 runs
- Release at CA-C-043 (Monterey), 2026-2050 CMIP6 SSP2-4.5
- 7 figures: trajectory, heatmap, trait evolution, survival, regional, dose-response
- Summary report at results/reintroduction/demo_summary.md" || echo "Nothing to commit"

git push origin main 2>&1 || echo "Push failed (non-critical)"

echo "COMPLETE: $(date)"
