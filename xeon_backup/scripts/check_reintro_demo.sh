#!/bin/bash
# Check reintroduction demo progress
LOG=/home/starbot/.openclaw/workspace/sswd-evoepi/results/reintroduction/demo_log.txt
RESULTS=/home/starbot/.openclaw/workspace/sswd-evoepi/results/reintroduction/results_demo.json

started=$(grep -c "Starting" "$LOG" 2>/dev/null || echo 0)
done=$(grep -c "Done in" "$LOG" 2>/dev/null || echo 0)
failed=$(grep -c "FAILED" "$LOG" 2>/dev/null || echo 0)
finished=$(grep -c "DONE" "$LOG" 2>/dev/null || echo 0)
workers=$(ps aux | grep "reintroduction_monterey" | grep -v grep | wc -l)

echo "Reintro demo: $done/$started done ($failed failed), $workers active workers"

if [ "$finished" -gt 0 ]; then
    echo "STATUS: COMPLETE"
    grep "EXIT_CODE" "$LOG" 2>/dev/null
    if [ -f "$RESULTS" ]; then
        n_results=$(python3 -c "import json; print(len(json.load(open('$RESULTS'))))" 2>/dev/null)
        echo "Results file: $n_results entries"
    fi
elif [ "$workers" -eq 0 ] && [ "$done" -lt "$started" ]; then
    echo "STATUS: CRASHED (no active workers but not all done)"
else
    echo "STATUS: RUNNING"
    # Show latest completion time
    grep "Done in" "$LOG" | tail -1
fi
