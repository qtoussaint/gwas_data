#!/usr/bin/env bash

#SBATCH --time=96:00:00   # walltime
#SBATCH --ntasks=1   # number of tasks
#SBATCH --cpus-per-task=16   # number of CPUs Per Task
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem=100GB   # memory per node
#SBATCH -J "pyseer_olig_30"   # job name
#SBATCH -o "pyseer_olig_30.out"   # job output file
#SBATCH -e "pyseer_olig_30.err"   # job error file

path=/nfs/research/jlees/jacqueline/gwas_data/klebsiella/pyseer_data

# get kinship matrix for pyseer
similarity_pyseer --pres $path/oligogenic_genotypes.rtab $path/oligogenic_30_sample_list.txt > $path/oligogenic_genotype_kinship.tsv

# run pyseer (fast LMM)
pyseer --lmm --phenotypes $path/oligogenic_30_phenotypes.txt --pres $path/oligogenic_genotypes.rtab --similarity $path/oligogenic_genotype_kinship.tsv --output-patterns $path/oligogenic_30_snp_patterns.txt --cpu 16 --continuous > $path/oligogenic_30_snps.txt