#!/usr/bin/env bash

#SBATCH --time=48:00:00   # walltime
#SBATCH --ntasks=1   # number of tasks
#SBATCH --cpus-per-task=4   # number of CPUs Per Task
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem=10GB   # memory per node
#SBATCH -J "similarity"   # job name
#SBATCH -o "similarity.out"   # job output file
#SBATCH -e "similarity.err"   # job error file

awk -F'\t' '{
    # Print the first field followed by a space
    printf "%s ", $1;
    
    # Print the remaining fields separated by tabs
    for (i=2; i<=NF; i++) {
        printf "%s", $i;
        if (i < NF) {
            printf "\t";
        }
    }
    # Print a newline at the end of each line
    printf "\n";
}' genotypes.Rtab > genotypes_awk.rtab

similarity_pyseer --pres genotypes_awk.rtab sample_list.txt > genotype_kinship.tsv
