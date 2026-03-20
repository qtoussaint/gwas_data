#!/usr/bin/env bash

#SBATCH --job-name=transpose   # job name
#SBATCH --nodes=1                  # number of nodes
#SBATCH --cpus-per-task=48         # CPUs/task (i.e. if your code is multi-threaded)
#SBATCH --mem=290G                 # memory (add enough for CPPRATE)
#SBATCH --time=06:00:00            # time limit
#SBATCH --error=/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/cryptic_regeno_snpeff/transpose.err     # error file name with job ID
#SBATCH --output=/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/cryptic_regeno_snpeff/transpose.out   # output file name with job ID

#################################################################################

wd=/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/cryptic_regeno_snpeff

Rscript /nfs/research/jlees/jacqueline/gwas_data/code/gwas_data/tuberculosis/presence_absence_from_VCF.R

# activate datamash conda env

datamash transpose --no-strict < $wd/cryptic_regeno_MAF_0.0004_effannotated_presabs.tsv > $wd/cryptic_regeno_MAF_0.0004_effannotated_presabs_transposed.tsv

awk 'NR==FNR {order[NR]=$1; pos[$1]=NR; next}
     FNR==1 {print; next}
     (FNR-1) in pos {lines[pos[FNR-1]]=$0}
     END {for (i=1; i<=length(order); i++) print lines[i]}' $wd/sample_index.csv $wd/cryptic_regeno_MAF_0.0004_effannotated_presabs_transposed.tsv > $wd/subset_K8.tsv