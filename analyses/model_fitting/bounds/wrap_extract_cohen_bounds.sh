#!/bin/bash -l
#SBATCH -J cohen_bounds
#SBATCH -p priority
#SBATCH -A prio_skr2
#SBATCH -c 4
#SBATCH --mem=96G
#SBATCH -t 4:00:00
#SBATCH -o slurm-bounds-%j.out
#SBATCH -e slurm-bounds-%j.err

module load miniconda
conda activate tz
export PYTHONPATH=/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa:${PYTHONPATH:-}
cd /nfs/roberts/project/pi_skr2/mcn26/tabula-rasa/analyses/model_fitting/bounds
python extract_cohen_bounds.py
rc=$?
echo "EXITING (rc=$rc)"
exit $rc
