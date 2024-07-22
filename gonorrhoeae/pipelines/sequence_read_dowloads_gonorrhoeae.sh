#!/usr/bin/env bash

#BSUB -J gon_download # job name
#BSUB -n 8          # number of tasksi/CPUs in job
#BSUB -q standard   # queue
#BSUB -e error.%J   # error file name in which %J is replaced by the job ID
#BSUB -o output.%J  # output file name in which %J is replaced by the job ID

#################################################################################

# Bash script to download sequence reads (fastq format)

# Make text file with the names of sample accession IDs that you want to download
# Here it's column 2 of a data file called 41564_2018_309_MOESM3_ESM.csv

awk -F "\"*,\"*" '{print $2}' /nfs/research/jlees/jacqueline/gwas_datasets/gonorrhoeae/pipelines/41564_2018_309_MOESM3_ESM.csv > gon_accession_ids.txt

# This line just removes the first 3 lines since there is a header on the .csv I extracted the IDs from
sed -i '1,3d' gon_accession_ids.txt

# Using ffq, read through the sample IDs and make a file with their ftp URLs

cat gon_accession_ids.txt | while read line || [[ -n $line ]]; do ffq --ftp $line | grep -Eo '"url": "[^"]*"' | grep -o '"[^"]*"$' >> complete_file.txt; done

# Remove quotes from the URLs in the file
sed -i 's/"//g' complete_file.txt
# Remove bam files (since I only want the fastq files)
sed -i '*.bam' complete_file.txt

# Use wget to download the fastq files from the URLs in the previously created file

cat complete_file.txt | while read line || [[ -n $line ]]; do wget -nc $line ; done

# Completion notice

echo "Job successfully completed. That's pretty cool!"
