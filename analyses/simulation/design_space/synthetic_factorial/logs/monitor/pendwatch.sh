#!/bin/bash
# Sample pending synfac_worker jobs every 15 min and log Slurm's Reason.
# Reason is only readable while a job is PENDING -- sacct clears it to "None"
# once the job starts -- so it has to be caught live.
#
# When our workers ARE starved, also snapshot who else holds the partition.
# Our pending Reason is (Priority), i.e. we are queued behind other users, so
# the useful evidence is what those other jobs are. squeue without -u/--me
# shows everybody. The extra queries only fire during a starvation window,
# which is rare, so this stays gentle on the scheduler.
OUT="$(dirname "$0")"/pending_reasons.log
while true; do
    p=$(squeue --user=mcn26 --name=synfac_worker -h -t PD -o '%i %.12M %R' 2>/dev/null)
    if [ -n "$p" ]; then
        {
            echo "=== $(date '+%Y-%m-%d %H:%M:%S') ==="
            echo "$p"
            echo "--- top users running on priority (user cpus jobs) ---"
            squeue -p priority -t R -h -o '%u %C' 2>/dev/null \
              | awk '{c[$1]+=$2; n[$1]++} END {for (u in c) printf "%s %d %d\n", u, c[u], n[u]}' \
              | sort -k2 -rn | head -8
            echo "--- partition ---"
            sinfo -p priority -h -o '%t %D %C'
        } >> "$OUT"
    fi
    sleep 900
done
