#!/usr/bin/env bash

#SBATCH --time=96:00:00   # walltime
#SBATCH --ntasks=1   # number of tasks
#SBATCH --cpus-per-task=16   # number of CPUs Per Task
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem=100GB   # memory per node
#SBATCH -J "fasttree_tb"   # job name
#SBATCH -o "fasttree_tb.out"   # job output file
#SBATCH -e "fasttree_tb.err"   # job error file

# run parsnp to generate core genome alignment
parsnp -r /nfs/research/jlees/jacqueline/gwas_data/tuberculosis/pyseer_data/trees/reference/GCF_000195955.2/GCF_000195955.2_ASM19595v2_genomic.fna -d /nfs/research/jlees/jacqueline/gwas_data/tuberculosis/genotype/assemblies_for_pyseer/*.fa -o /nfs/research/jlees/jacqueline/gwas_data/tuberculosis/pyseer_data/trees/parsnp

#fasttree -nt 

#/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/genotype/assemblies_for_pyseer/*.fa