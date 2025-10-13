#!/usr/bin/env bash

#SBATCH --job-name=transpose   # job name
#SBATCH --nodes=1                  # number of nodes
#SBATCH --cpus-per-task=10         # CPUs/task (i.e. if your code is multi-threaded)
#SBATCH --mem=100G                 # memory (add enough for CPPRATE)
#SBATCH --time=06:00:00            # time limit
#SBATCH --error=/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/cryptic_regeno_snpeff/transpose.err     # error file name with job ID
#SBATCH --output=/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/cryptic_regeno_snpeff/transpose.out   # output file name with job ID

#################################################################################

wd=/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/cryptic_regeno_snpeff
dataset=RIF_K8

# activate bcftools conda env

# filter out unwanted samples for given dataset
#bcftools view -S $wd/sample_index_$dataset.txt $wd/cryptic_regeno_MAF_0.0004_effannotated_binary.vcf.gz -Oz -o $wd/subsampled_presabs/filtered_$dataset.vcf.gz
#bcftools index $wd/subsampled_presabs/filtered_$dataset.vcf.gz

# extract pres/abs matrix
#bcftools query -f '%CHROM\_%POS\_%REF\_%ALT[\t%GT]\n' "$wd/subsampled_presabs/filtered_$dataset.vcf.gz" > "$wd/subsampled_presabs/presabs_unformatted_$dataset.tsv"

# ACTIVATE DATAMASH ENV

# transpose
datamash transpose < $wd/subsampled_presabs/presabs_unformatted_$dataset.tsv > $wd/subsampled_presabs/presence_absence_T_$dataset.tsv

# remove the variant index and output to separate file
head -n 1 $wd/subsampled_presabs/presence_absence_T_$dataset.tsv > $wd/subsampled_presabs/variant_index_presence_absence_$dataset.tsv

# remove variants
sed '1d' $wd/subsampled_presabs/presence_absence_T_$dataset.tsv > presence_absence_final_$dataset.tsv