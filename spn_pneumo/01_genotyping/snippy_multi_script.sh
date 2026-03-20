#!/usr/bin/env bash

#SBATCH --time=48:00:00             # walltime
#SBATCH --ntasks=1                  # number of tasks
#SBATCH --cpus-per-task=48          # number of CPUs Per Task
#SBATCH --nodes=1                   # number of nodes
#SBATCH --mem=150GB                  # memory per node
#SBATCH -J "snippy"             # job name
#SBATCH -o "snippy.out"         # job output file
#SBATCH -e "snippy.err"         # job error file

#################################################################################

outdir=/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype

cd $outdir

# Make a list of all assemblies
#find /nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/sequences/*.fa -type f > sequence_paths.txt
#cat sequence_paths.txt | cut -f10 -d"/" > sequence_ids.txt
#cat sequence_ids.txt | cut -f1 -d"." > sequence_ids_new.txt
#mv sequence_ids_new.txt sequence_ids.txt # clean up intermediary

# make tab-separated input file for snippy-multi
#paste sequence_ids.txt sequence_paths.txt > snippy_multi_input.tab

# Use snippy to run samples
snippy-multi $outdir/snippy_multi_input.tab --ref /nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/pyseer_data/pyseer_tutorial/Spn23F.fa --cpus 47 --outdir $outdir/snippy_multi > snippy_multi.sh

# looks like snippy-multi just simplifies the process of running each sample individually, rather than running them jointly.