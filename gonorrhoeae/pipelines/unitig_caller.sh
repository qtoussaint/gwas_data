#!/usr/bin/env bash

#BSUB -J gon_unitig_caller # job name
#BSUB -n 10          # number of tasksi/CPUs in job
#BSUB -q standard   # queue
#BSUB -e error.%J   # error file name in which %J is replaced by the job ID
#BSUB -o output.%J  # output file name in which %J is replaced by the job ID

#################################################################################

# Note: use unitig_caller conda environment

unitig-caller --call --reads reads.txt --write-graph --out /nfs/research/jlees/jacqueline/gwas_datasets/gonorrhoeae/genotype/unitigs/ --threads 10
