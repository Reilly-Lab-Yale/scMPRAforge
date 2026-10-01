#!/bin/bash -l
#SBATCH -J cohen_modelsel
#SBATCH -p priority
#SBATCH -A prio_skr2
#SBATCH -c 4
#SBATCH --mem=96G
#SBATCH -t 6:00:00
#SBATCH -o slurm-modelsel-%j.out
#SBATCH -e slurm-modelsel-%j.err
module load miniconda
conda activate tz
export PYTHONPATH=/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa:${PYTHONPATH:-}
cd /nfs/roberts/project/pi_skr2/mcn26/tabula-rasa/analyses/model_selection
rc=0
for s in cross_family_agreement.py cohen_obs_vs_obsingle.py lrt_nb_vs_zinb.py nonzero_expression.py; do
  echo "=== $s"; python "$s" || { echo "FAILED: $s"; rc=1; }
done
echo "EXITING (rc=$rc)"; exit $rc
