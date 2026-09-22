#!/bin/bash
#SBATCH --job-name=synfac_smoke
#SBATCH --partition=priority
#SBATCH --account=prio_skr2
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --output=logs/smoke-%j.out
#SBATCH --error=logs/smoke-%j.err

# Validation run for the multi-arm sweep, on the real union LHS box but a
# handful of points. Confirms all 8 arms produce results, the cache is balanced
# across arms, and the sim tree is pruned -- and times a sample so the slow
# run's array concurrency can be set from a measurement rather than a guess.
#
#   sbatch wrap_union_smoke.sh

module load miniconda
conda activate tz

export SYNTHETIC_FACTORIAL_SIM_ROOT=/nfs/roberts/scratch/pi_skr2/mcn26/synthetic_factorial_sims/2026-08-18_synthetic_factorial_smoke
export SYNTHETIC_FACTORIAL_CACHE_ROOT=/nfs/roberts/scratch/pi_skr2/mcn26/synthetic_factorial_smoke_caches

python synthetic_factorial.py simulate smoke 0 1
echo "EXITING SHELL"
