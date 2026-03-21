#!/bin/bash
cd /home/starbot/.openclaw/workspace/sswd-evoepi
python3 scripts/sensitivity/run_morris_r3.py > results/sensitivity_r3/morris_run.log 2>&1
