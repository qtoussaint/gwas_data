#!/bin/bash
#SBATCH --job-name=spn_unitigs
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/unitigs/unitig_caller.%j.out
#SBATCH --error=/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/unitigs/unitig_caller.%j.err

set -euo pipefail

WORKDIR=/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo
cd "$WORKDIR"

source /hps/software/users/jlees/jacqueline/anaconda3/etc/profile.d/conda.sh
conda activate /hps/software/users/jlees/jacqueline/anaconda3/envs/unitig_caller

echo "unitig-caller: $(which unitig-caller)"
unitig-caller --version
echo "refs: $(wc -l < unitigs/refs.txt) samples"
echo "start: $(date)"

unitig-caller --call \
  --refs    unitigs/refs.txt \
  --kmer    31 \
  --pyseer  --rtab \
  --out     unitigs/spn_unitigs \
  --threads "${SLURM_CPUS_PER_TASK:-8}"

echo "unitig-caller done: $(date)"

# Compress the long pyseer-format file for pyseer --kmers
gzip -f unitigs/spn_unitigs.pyseer
echo "gzip done: $(date)"

ls -l unitigs/spn_unitigs.*
