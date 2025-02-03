#!/usr/bin/env bash

#SBATCH --time=48:00:00   # walltime
#SBATCH --ntasks=1   # number of tasks
#SBATCH --cpus-per-task=4   # number of CPUs Per Task
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem=10GB   # memory per node
#SBATCH -J "pyseer_spn"   # job name
#SBATCH -o "pyseer_spn.out"   # job output file
#SBATCH -e "pyseer_spn.err"   # job error file

# specify path
path=/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo

############

# run similarity to get kinship matrix
#similarity_pyseer --pres $path/pyseer_data/genotypes.rtab $path/pyseer_data/sample_list.txt > $path/pyseer_data/genotype_kinship.tsv

# run pyseer (fast LMM) using kinship matrix
#pyseer --lmm --phenotypes $path/pyseer_data/phenotypes_pen_mic.txt --pres $path/pyseer_data/genotypes.rtab --similarity $path/pyseer_data/genotype_kinship.tsv --output-patterns $path/pyseer_data/snp_patterns.txt --cpu 16 --continuous > $path/pyseer_data/pen_mic_snps.txt

############

# get phylogeny distances from tree
phylogeny_distance.py --lmm $path/pyseer_data/parsnp/core_genome_aln.tree > $path/pyseer_data/phylogeny_similarity.tsv

# run pyseer (fast LMM) using fasttree
pyseer --lmm --phenotypes $path/pyseer_data/phenotypes_pen_mic.txt --pres $path/pyseer_data/genotypes.rtab --similarity $path/pyseer_data/phylogeny_similarity.tsv --output-patterns $path/pyseer_data/snp_patterns.txt --cpu 16 --continuous > $path/pyseer_data/pen_mic_snps.txt

