#!/bin/bash
# Re-run of the Yin et al. canonical fit purely to benchmark it on the same
# CPU generation as the other two datasets in the compute figure.
#
# The original ran on AMD Turin (a1114u01n01) while Lalanne and Zhao ran on
# Intel Emerald Rapids, so their wall times were not comparable. The priority
# partition is now entirely Turin, so this uses `day`, which still carries 92
# Emerald Rapids nodes, with the default (non-priority) pi_skr2 account.
#
# The fitted parameters are written to scratch and are disposable; the point
# is the SLURM accounting record. Logs stay in this directory, which is on
# project and backed up.
#SBATCH --job-name=seel_er_bench
#SBATCH --partition=day
#SBATCH --constraint=cpugen:emeraldrapids
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
#SBATCH --time=06:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --exclude=a1132u18n02

module load miniconda
eval "$(conda shell.bash hook)"
conda activate tz

python3 fit.py
