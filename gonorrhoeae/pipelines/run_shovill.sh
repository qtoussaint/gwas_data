#!/usr/bin/env bash

#BSUB -J shovill[1-1594]%100 # job name
#BSUB -n 10          # number of tasksi/CPUs in job
#BSUB -q standard   # queue
#BSUB -e error.%J_%I   # error file name in which %J is replaced by the job ID
#BSUB -o output.%J_%I  # output file name in which %J is replaced by the job ID 

################################################################################

#SUBMIT WITH %100 limit (in job name) to limit the number of jobs to run at once to 100!!!

# Use shovill to assemble reads for all samples
shovill --cpus 10 --outdir /nfs/research/jlees/jacqueline/gwas_datasets/gonorrhoeae/genotype/assembled_reads/assembly_$LSB_JOBINDEX --R1 /nfs/research/jlees/jacqueline/gwas_datasets/gonorrhoeae/genotype/sequence_reads/$LSB_JOBINDEX*_1.fastq.gz --R2 /nfs/research/jlees/jacqueline/gwas_datasets/gonorrhoeae/genotype/sequence_reads/$LSB_JOBINDEX*_2.fastq.gz
