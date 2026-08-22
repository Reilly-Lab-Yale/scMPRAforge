#!/bin/bash
# Sample pending synfac_worker jobs every 5 min and log Slurm's Reason.
# Reason is only readable while a job is PENDING -- sacct clears it to "None"
# once the job starts -- so it has to be caught live.
OUT="$(dirname "$0")"/pending_reasons.log
while true; do
    p=$(squeue --user=mcn26 --name=synfac_worker -h -t PD -o '%i %.12M %R' 2>/dev/null)
    if [ -n "$p" ]; then
        echo "=== $(date '+%Y-%m-%d %H:%M:%S') ===" >> "$OUT"
        echo "$p" >> "$OUT"
    fi
    sleep 900
done
