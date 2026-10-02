#!/bin/bash
#SBATCH -p priority_gpu
#SBATCH --gpus=h100:1
#SBATCH --exclude=a1127u42n01
#SBATCH -A prio_skr2
#SBATCH -c 4
#SBATCH --mem=128G
#SBATCH -t 4:00:00
#SBATCH -J cohen_5x5_test
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err

module load miniconda
conda activate tz
export PYTHONPATH=/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa:${PYTHONPATH:-}
module load CUDA cuDNN

cd /nfs/roberts/project/pi_skr2/mcn26/tabula-rasa/analyses/simulation/activity_prc/cohen

python cohen_5x5_activity.py test || exit $?
python cohen_5x5_activity.py metrics || exit $?
