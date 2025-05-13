#!/usr/bin/env bash

#SBATCH --time=48:00:00             # walltime
#SBATCH --ntasks=1                  # number of tasks
#SBATCH --cpus-per-task=48          # number of CPUs Per Task
#SBATCH --nodes=1                   # number of nodes
#SBATCH --mem=150GB                  # memory per node
#SBATCH -J "genotyping"             # job name
#SBATCH -o "/nfs/research/jlees/jacqueline/gwas_data/code/gwas_data/spn_pneumo/genotyping_spn.out"         # job output file
#SBATCH -e "/nfs/research/jlees/jacqueline/gwas_data/code/gwas_data/spn_pneumo/genotyping_spn.err"         # job error file

#################################################################################

outdir=/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype

cd $outdir

# chose to use bactsnp over snippy and bwa+bcftools based on results of SKA paper fig 3
/hps/software/users/jlees/jacqueline/manual_installs/bin/bactsnp/bactsnp-1.1.0.linux64/bactsnp --fasta_list sequence_paths.txt --reference ATCC_700669.gbff --thread 47 --out_dir $outdir/bactsnp

