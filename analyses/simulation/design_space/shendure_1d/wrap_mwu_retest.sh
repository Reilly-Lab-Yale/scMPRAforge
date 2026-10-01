#!/bin/bash
#SBATCH --job-name=mwu_retest
#SBATCH --partition=priority
#SBATCH --account=prio_skr2
#SBATCH --time=14:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

module load miniconda
conda activate tz

python mwu_retest.py "$@"
rc=$?
echo "EXITING SHELL (rc=$rc)"
exit $rc
