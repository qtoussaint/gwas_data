#!/usr/bin/env bash

#BSUB -J gon_snippy # job name
#BSUB -n 24          # number of tasksi/CPUs in job
#BSUB -q standard   # queue
#BSUB -e error.%J   # error file name in which %J is replaced by the job ID
#BSUB -o output.%J  # output file name in which %J is replaced by the job ID

#################################################################################

# Bash script to run snippy variant calling on sequence read fastq files

cd /nfs/research/jlees/jacqueline/gwas_datasets/gonorrhoeae/pipelines

# Make a list of all sequence reads in sequence_reads directory (sequence_reads_ids.txt)
find /nfs/research/jlees/jacqueline/gwas_datasets/gonorrhoeae/genotype/sequence_reads/*_1.fastq.gz -type f > sequence_paths.txt
cat sequence_paths.txt | cut -f10 -d"/" >> sequence_paths_suffixes.txt
cat sequence_paths_suffixes.txt | cut -f1 -d"_" >> sequence_reads_ids.txt

# Clean up intermediary files
rm sequence_paths.txt
rm sequence_paths_suffixes.txt

# Use snippy to run samples

cat /nfs/research/jlees/jacqueline/gwas_datasets/gonorrhoeae/pipelines/sequence_reads_ids.txt | while read line || [[ -n $line ]]; do snippy --cpus 24 --outdir /nfs/research/jlees/jacqueline/gwas_datasets/gonorrhoeae/genotype/variant_calls/results_$line --ref /nfs/research/jlees/jacqueline/gwas_datasets/gonorrhoeae/pipelines/reference_NC_002946.fasta --R1 /nfs/research/jlees/jacqueline/gwas_datasets/gonorrhoeae/genotype/sequence_reads/"$line"_1.fastq.gz --R2 /nfs/research/jlees/jacqueline/gwas_datasets/gonorrhoeae/genotype/sequence_reads/"$line"_2.fastq.gz; done

# Completion notice

echo "Job successfully completed. That's pretty cool!"
