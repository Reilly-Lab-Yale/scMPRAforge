#!/bin/bash
#SBATCH --job-name=synfac_slow
#SBATCH --partition=priority
#SBATCH --account=prio_skr2
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --output=logs/slow-%A_%a.out
#SBATCH --error=logs/slow-%A_%a.err

# Slow-burn wrapper for the multi-arm union sweep (5000 LHS points x 5 library
# reps x 5 sims, each tested under all 8 arms in ARMS).
#
#   sbatch --array=0-199%<concurrency> wrap_union_slow.sh
#
# 200 slices of ~25 samples each. Each slice driver stands up ONE dask cluster
# and walks its samples serially, so worker sbatches stay far below the YCRC
# 200/hr limit that the 2026-05-11 attempt tripped.
#
# Pace is the %N suffix, deliberately NOT in this file: to speed up or slow
# down, scancel the array and resubmit with a different %N. Nothing is lost.
# Every sample that already has a cache on project is skipped for the price of
# one stat, so resubmitting is always safe and always resumes -- after a
# timeout, a node failure, or a deliberate cancel.
#
# Progress:  ls <cache_root>/union/*.parquet | wc -l   (out of 5000)

module load miniconda
conda activate tz

# Dated scratch root so this run's transient sims are not confused with the
# 2026-05 single-arm sweep's. Sims here are written, tested, cached to project,
# and deleted within one sample; nothing in scratch outlives the run.
export SYNTHETIC_FACTORIAL_SIM_ROOT=/nfs/roberts/scratch/pi_skr2/mcn26/synthetic_factorial_sims/2026-08-18_synthetic_factorial_arms

python synthetic_factorial.py simulate union "$SLURM_ARRAY_TASK_ID" 200
echo "EXITING SHELL"
