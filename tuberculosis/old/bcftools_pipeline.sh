#!/usr/bin/env bash

#BSUB -J merge_tb_vcfs # job name
#BSUB -n 4          # number of tasksi/CPUs in job
#BSUB -q standard   # queue
#BSUB -e error.%J   # error file name in which %J is replaced by the job ID
#BSUB -o output.%J  # output file name in which %J is replaced by the job ID

#################################################################################

# DON"T FORGET to activate the bcftools conda environment beforehand

# Pipeline to concatonate individual VCFs (previously generated using snippy) using the bcftools package

bcftools merge /nfs/research/jlees/jacqueline/gwas_datasets/tuberculosis/genotype/variant_calls/results*/snps.vcf.gz -o merged_files.vcf.gz
