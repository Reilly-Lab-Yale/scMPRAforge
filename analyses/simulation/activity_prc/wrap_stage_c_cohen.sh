#!/bin/bash
#SBATCH -J cohen_stage_c
#SBATCH -p day
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -t 4:00:00
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err
module load miniconda
eval "$(conda shell.bash hook)"; conda activate tz
R=/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa
export PYTHONPATH=$R:${PYTHONPATH:-}
cd $R/analyses/simulation/activity_prc
python add_counts_tests.py cohen || exit $?
python add_ttest.py cohen || exit $?
python - <<'PY' || exit $?
import pandas as pd
d = pd.read_csv("cohen/output/cohen_5x5_activity_summary.tsv", sep="\t")
counts = d["test"].value_counts().to_dict()
print("summary rows", len(d), counts)
assert set(counts) == {"ks", "mwu", "pseudobulk", "ttest", "wald_auto"}, f"tests: {sorted(counts)}"
assert all(v == 25 for v in counts.values()), f"expected 25 per test: {counts}"
PY
python plot_test_comparison.py || exit $?
python plot_curves.py || exit $?
echo "STAGE C DONE"
