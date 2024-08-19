#!/usr/bin/env bash

#SBATCH --time=96:00:00   # walltime
#SBATCH --ntasks=1   # number of tasks
#SBATCH --cpus-per-task=16   # number of CPUs Per Task
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem=100GB   # memory per node
#SBATCH -J "pyseer_s_100"   # job name
#SBATCH -o "pyseer_s_100.out"   # job output file
#SBATCH -e "pyseer_s_100.err"   # job error file

path=/nfs/research/jlees/jacqueline/gwas_data/klebsiella/pyseer_data

# get kinship matrix for pyseer
similarity_pyseer --pres $path/single_genotypes.rtab $path/single_100_sample_list.txt > $path/single_genotype_kinship.tsv

# run pyseer (fast LMM)
pyseer --lmm --phenotypes $path/single_100_phenotypes.txt --pres $path/single_genotypes.rtab --similarity $path/single_genotype_kinship.tsv --output-patterns $path/single_100_snp_patterns.txt --cpu 16 --continuous > $path/single_100_snps.txt