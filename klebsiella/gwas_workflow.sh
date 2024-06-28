#!/usr/bin/env bash

#SBATCH --job-name=gwas_klebsiella # job name
#SBATCH --nodes=1                  # number of nodes
#SBATCH --cpus-per-task=10         # CPUs/task (i.e. if your code is multi-threaded)
#SBATCH --mem=30G                  # memory (whole job)
#SBATCH --time=96:00:00            # time limit
#SBATCH --partition=standard       # partition
#SBATCH --error=/nfs/research/jlees/jacqueline/gwas_data/klebsiella/gwas_runs/error.%j     # error file name with job ID
#SBATCH --output=/nfs/research/jlees/jacqueline/gwas_data/klebsiella/gwas_runs/output.%j   # output file name with job ID

#################################################################################

RSCRIPT_PATH="/nfs/research/jlees/jacqueline/gwas_workflow/code/gwas_workflow/gwas_pipeline.R"

DATA="/nfs/research/jlees/jacqueline/gwas_data/klebsiella/gwas_data_bundles/oligogenic_EF_30_H_0_9_associative_30fba35a0eaa67.json"
STAN_MODEL="/nfs/research/jlees/jacqueline/lasso_ridge_analysis/stan_models/standard/BBGWAS-cont-association.stan"
ANALYSIS_TYPE="inference"
ANALYSIS_NICKNAME="test"
OUTPUT_DIR="/nfs/research/jlees/jacqueline/gwas_data/klebsiella/gwas_runs"

R --no-save $RSCRIPT_PATH $DATA $STAN_MODEL $ANALYSIS_TYPE $ANALYSIS_NICKNAME $OUTPUT_DIR
