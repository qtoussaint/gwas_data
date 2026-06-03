#!/bin/bash
#SBATCH --job-name=spn_transpose
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/unitigs/transpose.%j.out
#SBATCH --error=/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/unitigs/transpose.%j.err

set -euo pipefail
WORKDIR=/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo
cd "$WORKDIR"

source /hps/software/users/jlees/jacqueline/anaconda3/etc/profile.d/conda.sh
conda activate /hps/software/users/jlees/jacqueline/envs/pyseer

echo "python: $(which python); $(python -c 'import numpy; print("numpy", numpy.__version__)')"
echo "start: $(date)"

python unitigs/transpose_rtab.py \
  unitigs/spn_unitigs.rtab \
  unitigs/spn_unitigs.matrix.tsv \
  unitigs/spn_unitigs.unitig_ids.tsv

echo "done: $(date)"
ls -l unitigs/spn_unitigs.matrix.tsv unitigs/spn_unitigs.unitig_ids.tsv
