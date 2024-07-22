#!/usr/bin/env bash

#BSUB -J tb_gt        # job name
#BSUB -n 20           # number of tasksi/CPUs in job
#BSUB -q standard     # queue
#BSUB -e error.%J     # error file name in which %J is replaced by the job ID
#BSUB -o output.%J    # output file name in which %J is replaced by the job ID
#BSUB -M 100000       # request 100GB of memory (units are MB) as opposed to 2GB default. Try not to burn down cluster.

#################################################################################

# create filtered VCF with only variants >0.01 MAF
bcftools view -q 0.01:minor /nfs/research/zi/mhunt/Cryptic2/Regeno/Cryptic1/Out/merged.vcf.gz -Oz -o tb_mhunt_filtered.vcf.gz

# after filtering, index
bcftools index tb_mhunt_filtered.vcf.gz

# get genotype matrix
bcftools query -H -f '%POS\_%REF\_%ALT [ %GT]\n' tb_mhunt_filtered.vcf.gz --output tb_mhunt_genotypes_cleaned.txt
