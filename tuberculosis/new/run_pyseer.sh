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

# run pyseer (fast LMM)
pyseer --lmm --phenotypes $path/phenotypes_rfb_mic.txt --pres $path/genotypes_awk.rtab --similarity $path/genotype_kinship.tsv --output-patterns $path/snp_patterns.txt --cpu 16 --continuous > $path/rfb_mic_snps.txt

# generate mash sketch
#mash sketch -s 10000 -o /nfs/research/jlees/jacqueline/gwas_data/tuberculosis/pyseer_data/mash_sketch /nfs/research/jlees/jacqueline/gwas_data/tuberculosis/genotype/assemblies_for_pyseer/*.fa
#mash dist mash_sketch.msh mash_sketch.msh | square_mash > mash.tsv

# generate scree plot to determine max distances
#scree_plot_pyseer mash.tsv

