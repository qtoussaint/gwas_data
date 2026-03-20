#!/bin/bash
#SBATCH --job-name=shovill
#SBATCH --output=/nfs/research/jlees/jacqueline/gwas_data/code/gwas_data/tuberculosis/new/shovill_logs/shovill_%A_%a.out
#SBATCH --error=/nfs/research/jlees/jacqueline/gwas_data/code/gwas_data/tuberculosis/new/shovill_logs/shovill_%A_%a.err
#SBATCH --array=292-12325%10
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=250G
#SBATCH --time=00:45:00

# Directory containing the sequence reads
READS_DIR="/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/genotype/sequence_reads"

# Output directory
OUTDIR_BASE="/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/genotype/assemblies_for_pyseer"

# Get the SAMPLEID from the file
SAMPLEID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sample_ids.txt)

# Run shovill
shovill --cpus 10 --outdir ${OUTDIR_BASE}/${SAMPLEID}.fa --R1 ${READS_DIR}/${SAMPLEID}_1.fastq.gz --R2 ${READS_DIR}/${SAMPLEID}_2.fastq.gz

