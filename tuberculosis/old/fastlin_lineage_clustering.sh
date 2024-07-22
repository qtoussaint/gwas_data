#!/usr/bin/env bash

#BSUB -J fastlin_tb # job name
#BSUB -n 10         # number of tasksi/CPUs in job
#BSUB -q standard   # queue
#BSUB -e error.%J   # error file name in which %J is replaced by the job ID
#BSUB -o output.%J  # output file name in which %J is replaced by the job ID

#################################################################################

fastlin -d /nfs/research/jlees/jacqueline/gwas_datasets/tuberculosis/genotype/sequence_reads -b MTBC_barcodes_for_fastlin.tsv -x 80
