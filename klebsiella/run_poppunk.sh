#!/usr/bin/env bash

#SBATCH --time=96:00:00   # walltime
#SBATCH --ntasks=1   # number of tasks
#SBATCH --cpus-per-task=41   # number of CPUs Per Task i.e if your code is multi-threaded
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem=250GB   # memory per node
#SBATCH -J "poppunk_kleb"   # job name
#SBATCH -o "poppunk_kleb.out"   # job output file
#SBATCH -e "poppunk_kleb.err"   # job error file

# assign path
path=/nfs/research/jlees/jacqueline/gwas_data/klebsiella/lineages

poppunk_assign --db $path/full_pp_ref_database/Klebsiella_pneumoniae_v3_full --query $path/rfile_poppipe.tsv --output $path/poppunk_clusters --threads 40

# generate microreact visualization files

poppunk_visualise --ref-db $path/full_pp_ref_database/Klebsiella_pneumoniae_v3_full --model-dir $path/poppunk_clusters --output $path/poppunk_clusters --microreact --threads 40
