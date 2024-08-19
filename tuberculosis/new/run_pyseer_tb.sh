#!/usr/bin/env bash

#SBATCH --time=96:00:00   # walltime
#SBATCH --ntasks=1   # number of tasks
#SBATCH --cpus-per-task=16   # number of CPUs Per Task
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem=100GB   # memory per node
#SBATCH -J "pyseer_tb"   # job name
#SBATCH -o "pyseer_tb.out"   # job output file
#SBATCH -e "pyseer_tb.err"   # job error file

path=/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/pyseer_data

# get kinship matrix for pyseer
#similarity_pyseer --pres $path/genotypes.rtab $path/sample_list.txt > $path/genotype_kinship.tsv

# get phylogeny similarity matrix from parsnp tree
python scripts/phylogeny_distance.py --lmm trees/parnsp/parsnp.tree > phylogeny_similarity.tsv

# run pyseer (fast LMM)
pyseer --lmm --phenotypes $path/phenotypes_rfb_mic.txt --pres $path/genotypes.rtab --similarity $path/phylogeny_similarity.tsv --output-patterns $path/snp_patterns.txt --cpu 16 --continuous > $path/rfb_mic_snps.txt