#!/bin/bash
# Status of the multi-arm union sweep: how many of the 5000 samples are cached,
# how much scratch the in-flight ones are holding, and what is queued.
#
#   ./progress.sh

set -uo pipefail

CACHE_DIR="${SYNTHETIC_FACTORIAL_CACHE_ROOT:-/nfs/roberts/project/pi_skr2/shared/tabula_data/simulated/synthetic_factorial_caches}/union"
SIM_DIR="${SYNTHETIC_FACTORIAL_SIM_ROOT:-/nfs/roberts/scratch/pi_skr2/mcn26/synthetic_factorial_sims/2026-08-18_synthetic_factorial_arms}_union"

n_done=$(find "$CACHE_DIR" -maxdepth 1 -name '*.parquet' | wc -l)
echo "cached:      $n_done / 5000  ($((n_done * 100 / 5000))%)"

# Non-zero here is normal mid-write; persistently non-zero with nothing running
# means a driver died between the write and the rename.
n_tmp=$(find "$CACHE_DIR" -maxdepth 1 -name '*.parquet.tmp' | wc -l)
if [ "$n_tmp" -gt 0 ]; then
    echo "part-written: $n_tmp (.tmp)"
fi

if [ -d "$SIM_DIR" ]; then
    echo "in flight:   $(du -sh "$SIM_DIR" | cut -f1) on scratch"
else
    echo "in flight:   nothing on scratch yet"
fi
echo
squeue -u "$USER" -n synfac_slow,synfac_worker -o '%.18i %.12j %.2t %.11M %R'
