#!/bin/bash
#SBATCH --job-name=coarse_reporter_pilot
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=priority
#SBATCH --account=prio_skr2

module load miniconda
conda activate tz
set -euo pipefail

export PYTHONPATH=/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa:${PYTHONPATH:-}

python coarse_reporter_pilot.py
